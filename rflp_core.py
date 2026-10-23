import sys
import bisect
import csv
import gzip
import os
import re
import logging
import math
from typing import List, Dict, Tuple, Optional
from datetime import datetime

from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio.Restriction import AllEnzymes
from assay_scoring import GelModel, MODEL_VERSION, score_genotypes, rank_key

APP_VERSION = "1.3.0"

# Analysis defaults of the GUI. Options added after v1.1.0 are off here, so a
# parameter set written for v1.1.0 reproduces its results.
DEFAULT_PARAMS = dict(
    names2_file='', flank=300, min_frag=80, max_frag=800, delta=25, max_cuts=3,
    gain_loss_only=False, suppliers=[], use_primer3=True, prod_min=250, prod_max=600,
    tm_min=58.0, tm_max=62.0, joint_design=True, primer_pairs=10, top_results=10,
    gel_resolution_bp=10.0, gel_resolution_pct=3.0, gel_stress_factor=1.5,
    skip_norm=False, no_save_norm=True,
    commercial_only=False, group_isoschizomers=False, snp_end_distance=0,
    exclude_primer_sites=False, mask_repeats=False,
    background_vcf='', exclude_background_sites=False, dcaps='off', dcaps_mismatches=1,
    # dCAPS designs only: a longer mismatch primer moves the created cut away from the
    # product end, and a short product shows the ~20 bp difference on a gel.
    dcaps_tm_max=70.0, dcaps_prod_min=120, dcaps_prod_max=300, dcaps_delta=20,
)
# Settings chosen from the validation on published assays (GUI and CLI default).
RECOMMENDED_PARAMS = dict(DEFAULT_PARAMS, commercial_only=True, group_isoschizomers=True,
                          snp_end_distance=100, exclude_primer_sites=True, mask_repeats=True,
                          exclude_background_sites=True, dcaps='fallback')


class AnalysisCancelled(Exception):
    """Raised when a running analysis is cancelled by the caller."""


def validate_params(params: dict) -> list[str]:
    errors = []
    required = ('fasta', 'input_file')
    for key in required:
        if not str(params.get(key, '')).strip():
            errors.append(f"Не задан параметр: {key}")
    for key in ('fasta', 'input_file', 'background_vcf'):
        value = str(params.get(key, '') or '').strip()
        if value and not os.path.isfile(value):
            errors.append(f"Файл не найден ({key}): {value}")
    for lo, hi, label in (('min_frag', 'max_frag', 'размер фрагмента'),
                          ('prod_min', 'prod_max', 'размер продукта'),
                          ('tm_min', 'tm_max', 'Tm'),
                          ('dcaps_prod_min', 'dcaps_prod_max', 'размер продукта dCAPS'),
                          ('tm_min', 'dcaps_tm_max', 'Tm праймеров dCAPS')):
        try:
            # Parameter sets written before the dCAPS settings existed get their defaults.
            lower, upper = (float(params.get(key, DEFAULT_PARAMS[key]) if key.startswith('dcaps_')
                                  else params.get(key)) for key in (lo, hi))
            if not all(math.isfinite(x) and x > 0 for x in (lower, upper)):
                errors.append(f"Некорректный параметр: {lo}/{hi}")
            elif lower > upper:
                errors.append(f"Минимальный {label} больше максимального")
        except (TypeError, ValueError):
            errors.append(f"Некорректный параметр: {lo}/{hi}")
    for key, default, low, high in (
        ('primer_pairs', 10, 1, 100), ('top_results', 10, 0, 10000),
        ('flank', 300, 1, 100000), ('max_cuts', 3, 0, 100),
        ('delta', 25, 1, 100000), ('snp_end_distance', 0, 0, 5000), ('dcaps_mismatches', 1, 1, 2),
        ('dcaps_delta', 20, 1, 1000),
    ):
        value = params.get(key, default)
        if isinstance(value, bool) or not isinstance(value, int) or not low <= value <= high:
            errors.append(f"{key}: требуется целое число от {low} до {high}")
    if params.get('dcaps', 'off') not in ('off', 'fallback', 'always'):
        errors.append("dcaps: допустимо off, fallback или always")
    try:
        gel_model_from_params(params)
    except (ValueError, TypeError) as exc:
        errors.append(str(exc))
    return errors


def gel_model_from_params(params):
    return GelModel(params.get('min_frag', 80), params.get('max_frag', 800),
                    params.get('gel_resolution_bp', 10.0),
                    params.get('gel_resolution_pct', 3.0),
                    params.get('gel_stress_factor', 1.5))


# ====== primer3 ======
_PRIMER3_AVAILABLE = False
_PRIMER3_USE_NEW = False
_PRIMER3_IMPORT_ERROR = ""

try:
    import primer3
    _PRIMER3_AVAILABLE = True
    _PRIMER3_USE_NEW = hasattr(primer3.bindings, "design_primers")
except Exception as e:
    _PRIMER3_IMPORT_ERROR = repr(e)
    _PRIMER3_AVAILABLE = False
    _PRIMER3_USE_NEW = False

logger = logging.getLogger('RFLPPickerLogger')
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
# убираем файловый логгер — только консоль
if not logger.handlers:
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

# ====== базовая директория (для PyInstaller и для .py) ======
def get_base_dir() -> str:
    """
    Если собрано в exe (PyInstaller) — использовать папку рядом с exe.
    Иначе — папку, где лежит .py.
    """
    if getattr(sys, "frozen", False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))

# ====== supplier mapping ======
_SUPPLIER_CANON = {
    "A": ["agilent"],
    "E": ["thermo fisher", "thermo scientific", "fermentas", "thermo"],
    "F": ["promega"],
    "H": ["minotech"],
    "I": ["sibenzyme", "sib enzyme"],
    "J": ["nippon gene"],
    "K": ["takara", "takara bio"],
    "M": ["roche", "roche applied science"],
    "N": ["new england biolabs", "neb"],
    "R": ["molecular biology resources", "chimerx"],
    "S": ["sigma", "sigma chemical"],
    "T": ["toyobo"],
    "V": ["vivantis"],
    "Y": ["eurx"],
    "Z": ["sinaclon", "sinaclon bioscience"],
}


def _name_to_codes(name: str) -> set[str]:
    s = (name or "").strip().lower()
    out: set[str] = set()
    if not s:
        return out
    if len(s) == 1 and s.upper() in _SUPPLIER_CANON:
        return {s.upper()}
    if s.isalpha() and s.upper() == s and len(s) <= 8 and " " not in s:
        return {ch for ch in s if ch.upper() in _SUPPLIER_CANON}
    for code, aliases in _SUPPLIER_CANON.items():
        for a in aliases:
            if a in s:
                out.add(code)
                break
    return out


def parse_suppliers_arg(s: Optional[str]) -> set[str]:
    if not s:
        return set()
    tokens = [t.strip() for t in re.split(r"[,\|;]", s) if t.strip()]
    codes: set[str] = set()
    for t in tokens:
        codes |= _name_to_codes(t)
    return {c.upper() for c in codes if c.upper() in _SUPPLIER_CANON}


def enzyme_supplier_codes(e) -> set[str]:
    items = []
    try:
        items = list(e.supplier_list())
    except Exception:
        pass
    if not items:
        val = getattr(e, "suppliers", None)
        try:
            if callable(val):
                val = val()
        except Exception:
            val = None
        if isinstance(val, dict):
            items = list(val.keys()) + list(val.values())
        elif isinstance(val, (list, tuple, set)):
            items = list(val)
        elif isinstance(val, str) and val.strip():
            items = [val]
    codes: set[str] = set()
    for it in items:
        s = str(it).strip()
        if len(s) == 1 and s.upper() in _SUPPLIER_CANON:
            codes.add(s.upper())
        else:
            codes |= _name_to_codes(s)
    return codes


def load_names2(path: str) -> Dict[str, str]:
    m: Dict[str, str] = {}
    try:
        if not path:
            return m
        if not os.path.isfile(path):
            logger.error(f"Файл names не найден: {path}")
            return m
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = re.split(r"[\t, ]+", line)
                if len(parts) < 2:
                    continue
                key, val = parts[0], parts[1]
                m[key] = val
    except Exception as e:
        logger.error(f"Ошибка при загрузке names: {e}")
    return m


def detect_delim(header: str) -> str:
    if "\t" in header:
        return "\t"
    if "," in header:
        return ","
    return "\t"


ID_COLUMNS = ("id", "name", "snp_id", "rsid")


def read_input_table(path: str) -> Tuple[Dict[str, Optional[int]], List[Tuple[int, List[str]]], str]:
    """Columns, data lines and delimiter of a variants file.

    ``columns`` maps chrom/pos/ref/alt/id to a field index, or None when the
    header has no such column (or there is no header). Data lines are
    ``(line_no, fields)`` with blank lines skipped and, after a header, the
    template lines containing «должно быть:». A missing file gives
    ``({}, [], "")``.
    """
    if not os.path.isfile(path):
        logger.error(f"Файл с вариантами не найден: {path}")
        return {}, [], ""
    with open(path, "r", encoding="utf-8-sig", errors="replace", newline="") as f:
        first = f.readline()
        if not first:
            return {}, [], ""
        delim = detect_delim(first)
        f.seek(0)
        reader = csv.reader(f, delimiter=delim)
        first_fields = next(reader, [])
        hdr = [h.strip().lower() for h in first_fields]

        def idx(*cands):
            for c in cands:
                if c in hdr:
                    return hdr.index(c)
            return None

        columns = dict(chrom=idx("chrom", "chr", "chr1", "chromosome"), pos=idx("pos", "position"),
                       ref=idx("ref", "reference"), alt=idx("alt", "alternate", "alt_allele"),
                       id=idx(*ID_COLUMNS))
        has_header = any(
            x in hdr
            for x in ("chrom", "chr", "pos", "ref", "alt", "position", "reference", "alternate")
        )
        lines: List[Tuple[int, List[str]]] = []
        if not has_header:
            # As before, a header-like token in a data line still sets its index;
            # an ID column needs a real header.
            columns["id"] = None
            lines.append((1, first_fields))
        for line_no, parts in enumerate(reader, start=2):
            if not parts or not any(x.strip() for x in parts):
                continue
            if has_header and "должно быть:" in delim.join(parts).lower():
                continue
            lines.append((line_no, parts))
    return columns, lines, delim


def parse_fields(parts: List[str], columns: Dict[str, Optional[int]]) -> Tuple[str, int, str, str]:
    """(chrom, pos, ref, alt) of one data line; raises IndexError/ValueError."""
    def field(name, default_index):
        index = columns.get(name)
        return parts[index if index is not None else default_index]

    return (
        field("chrom", 0).replace("chr", "").strip(),
        int(field("pos", 1)),
        field("ref", 2).strip().upper(),
        field("alt", 3).strip().upper(),
    )


def parse_input(path: str, report: Optional[dict] = None) -> List[Tuple[str, int, str, str]]:
    rows: List[Tuple[str, int, str, str]] = []
    try:
        columns, lines, delim = read_input_table(path)
        for line_no, parts in lines:
            if len(parts) < 4:
                if report is not None:
                    report.setdefault('rejected', []).append({'line': line_no, 'reason': 'parse_error', 'details': 'Ожидалось минимум 4 поля', 'raw': delim.join(parts)})
                continue
            try:
                rows.append(parse_fields(parts, columns))
            except Exception as e2:
                logger.error(f"Ошибка парсинга строки: {e2}")
                if report is not None:
                    report.setdefault('rejected', []).append({'line': line_no, 'reason': 'parse_error', 'details': str(e2), 'raw': delim.join(parts)})
    except Exception as e:
        logger.error(f"Ошибка при чтении файла с вариантами: {e}")
    return rows

# ===== НОРМАЛИЗАЦИЯ ВАРИАНТОВ ЧЕРЕЗ PYFAIDX =====

DNA = set("ACGT")
RC_MAP = str.maketrans("ACGTN", "TGCAN")


def rc(s: str) -> str:
    return s.translate(RC_MAP)


def revcomp(s: str) -> str:
    """Reverse-complement (5'->3' to 5'->3' on opposite strand)."""
    return s.translate(RC_MAP)[::-1]


def alt_ok(alt: str) -> bool:
    parts = [p for p in alt.split(",") if p != ""]
    if not parts:
        return True
    for p in parts:
        if p == "-":
            continue
        if not set(p) <= DNA:
            return False
    return True


def is_simple_snp(ref: str, alt: str) -> bool:
    parts = [p for p in alt.split(",") if p != ""]
    return (len(ref) == 1) and (len(parts) == 1) and (parts[0] in DNA)


class FastaBlockCache:
    """
    Блочный кэш FASTA на основе pyfaidx:
    держим в памяти один континг и окно [blk_start, blk_end] (1-based, включительно).
    """

    def __init__(self, fasta_path: str, block_bp: int = 100_000):
        self.fa = Fasta(
            fasta_path,
            as_raw=True,
            sequence_always_upper=True,
        )
        self.block_bp = max(10000, int(block_bp))
        self.chrom = None
        self.blk_start = 0
        self.blk_end = -1
        self.seq = ""

    def has_contig(self, chrom: str) -> bool:
        return chrom in self.fa

    def contig_len(self, chrom: str) -> int:
        return len(self.fa[chrom])

    def _load_block(self, chrom: str, center_pos: int):
        clen = self.contig_len(chrom)
        half = self.block_bp // 2
        start = max(1, center_pos - half)
        end = min(clen, start + self.block_bp - 1)
        start = max(1, end - self.block_bp + 1)
        self.seq = str(self.fa[chrom][start - 1 : end])
        self.chrom = chrom
        self.blk_start = start
        self.blk_end = end

    def base(self, chrom: str, pos1: int) -> str:
        if chrom != self.chrom or not (self.blk_start <= pos1 <= self.blk_end):
            self._load_block(chrom, pos1)
        return self.seq[pos1 - self.blk_start]

    def slice(self, chrom: str, start1: int, end1: int) -> str:
        need_center = (start1 + end1) // 2
        if chrom != self.chrom or not (self.blk_start <= start1 and end1 <= self.blk_end):
            self._load_block(chrom, need_center)
            if not (self.blk_start <= start1 and end1 <= self.blk_end):
                return str(self.fa[chrom][start1 - 1 : end1])
        s = start1 - self.blk_start
        e = end1 - self.blk_start + 1
        return self.seq[s:e]


def normalize_variants_for_rflp(
    fasta_path: str,
    input_path: str,
    names2_path: Optional[str],
    block_bp: int = 100_000,
    progress_every: int = 100000,
    status_cb=None,
    skip_norm: bool = False,
    save_files: bool = True,
    report: Optional[dict] = None,
    progress_cb=None,
    cancel_cb=None,
) -> Tuple[List[Tuple[str, int, str, str]], List[str], str, str]:
    """
    Нормализация вариантов:

    - skip_norm=True: нормализацию пропускаем, просто читаем файл и возвращаем как есть.
    - skip_norm=False: выполняем нормализацию по FASTA (+names2), опционально пишем cleaned/rejected.

    Возвращает:
      cleaned_variants: List[(chrom,pos,ref,alt)] — CHROM остаётся в формате входного файла
      summary_lines: строки SUMMARY
      cleaned_path: путь к .cleaned.tsv (или "" если не сохраняли)
      rejected_path: путь к .rejected.tsv (или "" если не сохраняли)
    """
    variants_raw = parse_input(input_path, report=report)
    total = len(variants_raw)

    if report is not None:
        report.setdefault('summary', {})['total'] = total

    if skip_norm:
        summary_lines = [
            "=== NORMALIZE SUMMARY ===",
            "normalization_skipped          : 1",
            f"variants_passed_without_change : {total}",
        ]
        if status_cb:
            status_cb("\n" + "\n".join(summary_lines) + "\n")
        # Even when automatic REF/ALT swapping is disabled, only canonical SNPs
        # with a matching reference base are safe for the single-base RFLP path.
        valid = []
        check_cache = FastaBlockCache(fasta_path, block_bp=block_bp)
        names2_skip = load_names2(names2_path) if names2_path else {}
        for line_no, (chrom, pos, ref, alt) in enumerate(variants_raw, start=1):
            mapped = names2_skip.get(chrom, chrom)
            reason = ''
            details = ''
            if not is_simple_snp(ref, alt):
                reason, details = 'unsupported_variant', 'Поддерживаются только одиночные SNP A/C/G/T'
            elif not check_cache.has_contig(mapped):
                reason, details = 'bad_chrom', f"CHROM '{mapped}' отсутствует в FASTA"
            elif pos < 1 or pos > check_cache.contig_len(mapped):
                reason, details = 'pos_oob', 'Позиция вне FASTA'
            elif check_cache.base(mapped, pos).upper() != ref.upper():
                reason, details = 'ref_not_match', 'REF не совпадает с FASTA'
            if reason:
                if report is not None:
                    report.setdefault('rejected', []).append({'line': line_no, 'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt, 'reason': reason, 'details': details})
            else:
                valid.append((chrom, pos, ref, alt))
            if progress_cb:
                progress_cb(line_no, max(1, total))
        if report is not None:
            report.setdefault('summary', {}).update(accepted=len(valid), rejected=len(report.get('rejected', [])))
        try:
            check_cache.fa.close()
        except Exception:
            pass
        return valid, summary_lines, "", ""

    if status_cb:
        status_cb(f"\n[INFO] NORMALIZE: загружено вариантов: {total}")

    cache = FastaBlockCache(fasta_path, block_bp=block_bp)
    names2 = load_names2(names2_path) if names2_path else {}

    stats = dict(
        total=total,
        kept=0,
        kept_swapped_snp=0,
        rejected_bad_chrom=0,
        rejected_pos_oob=0,
        rejected_bad_ref_chars=0,
        rejected_bad_alt_chars=0,
        rejected_minus_strand_like=0,
        rejected_multi_alt_contains_ref=0,
        rejected_ref_not_match=0,
    )

    cleaned: List[Tuple[str, int, str, str]] = []
    rejected: List[Tuple[str, int, str, str, str, str]] = []

    def ref_matches(mapped_chrom: str, pos1: int, ref: str) -> bool:
        end = pos1 + len(ref) - 1
        if pos1 < 1 or end > cache.contig_len(mapped_chrom):
            return False
        return cache.slice(mapped_chrom, pos1, end) == ref

    for i, (chrom, pos, ref, alt) in enumerate(variants_raw, start=1):
        if cancel_cb and cancel_cb():
            try:
                cache.fa.close()
            except Exception:
                pass
            raise AnalysisCancelled()
        if progress_cb:
            progress_cb(i, max(1, total))
        orig_chrom = chrom
        mapped_chrom = names2.get(orig_chrom, orig_chrom)

        if progress_every and (i % progress_every == 0) and status_cb:
            status_cb(f"[NORMALIZE] {i}/{total} {orig_chrom}:{pos}")

        if not cache.has_contig(mapped_chrom):
            stats["rejected_bad_chrom"] += 1
            rejected.append(
                (orig_chrom, pos, ref, alt, "bad_chrom", f"CHROM '{mapped_chrom}' not in FASTA")
            )
            continue

        clen = cache.contig_len(mapped_chrom)
        if pos < 1 or pos > clen:
            stats["rejected_pos_oob"] += 1
            rejected.append(
                (orig_chrom, pos, ref, alt, "pos_oob", f"POS {pos} out of [1,{clen}]")
            )
            continue

        if (not ref) or (not set(ref) <= DNA):
            stats["rejected_bad_ref_chars"] += 1
            rejected.append(
                (orig_chrom, pos, ref, alt, "bad_ref_chars", "REF must be A/C/G/T only")
            )
            continue

        if not is_simple_snp(ref, alt):
            stats["rejected_bad_alt_chars"] += 1
            rejected.append(
                (orig_chrom, pos, ref, alt, "unsupported_variant", "Поддерживаются только одиночные SNP A/C/G/T")
            )
            continue

        if not alt_ok(alt):
            stats["rejected_bad_alt_chars"] += 1
            rejected.append(
                (orig_chrom, pos, ref, alt, "bad_alt_chars", "ALT must be A/C/G/T or '-'")
            )
            continue

        alt_parts = [p for p in alt.split(",") if p != ""]

        if ref_matches(mapped_chrom, pos, ref):
            cleaned.append((orig_chrom, pos, ref, alt))
            stats["kept"] += 1
            continue

        if is_simple_snp(ref, alt) and len(alt_parts) == 1:
            refbase = cache.base(mapped_chrom, pos)
            if alt_parts[0] == refbase and ref != refbase:
                new_ref, new_alt = alt_parts[0], ref
                if ref_matches(mapped_chrom, pos, new_ref):
                    cleaned.append((orig_chrom, pos, new_ref, new_alt))
                    stats["kept_swapped_snp"] += 1
                    continue

        refbase = cache.base(mapped_chrom, pos)
        if len(ref) == 1 and rc(ref) == refbase:
            stats["rejected_minus_strand_like"] += 1
            rejected.append(
                (orig_chrom, pos, ref, alt, "minus_strand_like", "rc(REF) equals reference base")
            )
            continue

        if any(a == refbase for a in alt_parts):
            stats["rejected_multi_alt_contains_ref"] += 1
            rejected.append(
                (orig_chrom, pos, ref, alt, "multi_alt_contains_ref", "one ALT equals reference base")
            )
            continue

        stats["rejected_ref_not_match"] += 1
        rejected.append(
            (orig_chrom, pos, ref, alt, "ref_not_match", "REF does not match reference and no allowed fix")
        )

    # ----- сохранение cleaned/rejected в results/normalize (если включено) -----
    cleaned_path = ""
    rejected_path = ""
    if save_files:
        base_dir = get_base_dir()
        norm_dir = os.path.join(base_dir, "results", "normalize")
        os.makedirs(norm_dir, exist_ok=True)
        ts = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        outprefix = os.path.join(norm_dir, f"variants_{ts}")
        cleaned_path = f"{outprefix}.cleaned.tsv"
        rejected_path = f"{outprefix}.rejected.tsv"

        with open(cleaned_path, "w", encoding="utf-8", newline="") as fw:
            w = csv.writer(fw, delimiter="\t")
            w.writerow(["CHROM", "POS", "REF", "ALT"])
            for chrom, pos, ref, alt in cleaned:
                w.writerow([chrom, pos, ref, alt])

        with open(rejected_path, "w", encoding="utf-8", newline="") as fw:
            w = csv.writer(fw, delimiter="\t")
            w.writerow(["CHROM", "POS", "REF", "ALT", "reason", "details"])
            for row in rejected:
                w.writerow(list(row))

    # ----- SUMMARY -----
    summary_lines = []
    summary_lines.append("=== NORMALIZE SUMMARY ===")
    for k in (
        "total",
        "kept",
        "kept_swapped_snp",
        "rejected_bad_chrom",
        "rejected_pos_oob",
        "rejected_bad_ref_chars",
        "rejected_bad_alt_chars",
        "rejected_minus_strand_like",
        "rejected_multi_alt_contains_ref",
        "rejected_ref_not_match",
    ):
        summary_lines.append(f"{k:28s}: {stats[k]}")
    if save_files:
        summary_lines.append(f"cleaned_outfile           : {cleaned_path}")
        summary_lines.append(f"rejected_with_reasons     : {rejected_path}")
    else:
        summary_lines.append("files_saved               : 0 (by settings)")

    if status_cb:
        status_cb("\n" + "\n".join(summary_lines) + "\n")

    if report is not None:
        report.setdefault('summary', {}).update(accepted=len(cleaned), rejected=len(rejected))
        report['rejected'] = report.get('rejected', []) + [
            {'chrom': r[0], 'pos': r[1], 'ref': r[2], 'alt': r[3], 'reason': r[4], 'details': r[5]}
            for r in rejected
        ]
    try:
        cache.fa.close()
    except Exception:
        pass

    return cleaned, summary_lines, cleaned_path, rejected_path

# ====== дальнейшая логика RFLP =====

def cut_positions(enzyme, seq: str) -> List[int]:
    # A failed search must not be mistaken for an uncut allele.
    sites = enzyme.search(Seq(seq), linear=True)
    cuts = [0]
    cuts += [p - 1 for p in sites]
    L = len(seq)
    if L not in cuts:
        cuts.append(L)
    cuts = sorted(set(max(0, min(L, c)) for c in cuts))
    return cuts


def frag_lengths(cuts: List[int]) -> List[int]:
    return [cuts[i + 1] - cuts[i] for i in range(len(cuts) - 1)]


def frags_to_str(frags: List[int]) -> str:
    return "+".join(str(x) for x in frags) if frags else ""


def snp_site_enzymes(enzymes, seq: str, altseq: str, off: int) -> list:
    """Enzymes with a recognition site covering the SNP in REF or ALT.

    Bio.Restriction finds sites with fixed-width look-ahead patterns over
    ``' ' + sequence``.  A site that does not cover the SNP is therefore found
    identically in both alleles of every amplicon containing the SNP, so an
    enzyme without such a site always gives equal REF/ALT fragments (Δ = 0)
    and cannot become a candidate.  The searched context keeps the leading
    space: a site starting with N may match it at the start of an amplicon.
    """
    ref_data, alt_data = ' ' + seq, ' ' + altseq
    found = []
    for enz in enzymes:
        # Window positions whose sites cover data index off + 1.
        lo, hi = max(0, off - enz.size + 2), off + enz.size + 1
        if enz.compsite.search(ref_data, lo, hi) or enz.compsite.search(alt_data, lo, hi):
            found.append(enz)
    return found


def primer_site_label(enzyme, amplicon: str, left_len, right_len) -> str:
    """Which primers contain a whole recognition site of the enzyme.

    Such a site is in every product whatever the allele, and cutting that
    close to a DNA end is unreliable. ``amplicon`` is the product top strand,
    so the right primer is searched as its binding site.
    """
    found = []
    if left_len and enzyme.compsite.search(amplicon, 0, left_len):
        found.append('left')
    if right_len and enzyme.compsite.search(amplicon, len(amplicon) - right_len):
        found.append('right')
    return '; '.join(found)


def load_background_variants(path: str, regions: Dict[str, list], names: Optional[dict] = None,
                             cancel_cb=None) -> Tuple[Dict[str, list], int]:
    """Known variants of a VCF (plain or gzip) that lie inside ``regions``.

    ``regions`` maps FASTA contigs to 1-based inclusive ``(start, end)``
    windows. VCF contig names are translated with ``names`` (the chromosome
    names file) when listed there. Symbolic and missing ALT alleles are
    skipped. Returns ``{contig: [(pos, REF, (ALT, ...)), ...]}`` sorted by
    position and the number of VCF records read.
    """
    merged = {}
    for contig, spans in regions.items():
        out = []
        for start, end in sorted(spans):
            if out and start <= out[-1][1] + 1:
                out[-1][1] = max(out[-1][1], end)
            else:
                out.append([start, end])
        merged[contig] = ([a for a, _ in out], [b for _, b in out])
    names = names or {}
    found = {contig: [] for contig in regions}
    with open(path, 'rb') as probe:
        opener = gzip.open if probe.read(2) == b'\x1f\x8b' else open
    records = 0
    with opener(path, 'rt', encoding='utf-8', errors='replace') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            records += 1
            if cancel_cb and records % 200_000 == 0 and cancel_cb():
                raise AnalysisCancelled()
            fields = line.split('\t', 5)
            if len(fields) < 5:
                continue
            contig = names.get(fields[0], fields[0])
            spans = merged.get(contig)
            if spans is None:
                continue
            try:
                pos = int(fields[1])
            except ValueError:
                continue
            i = bisect.bisect_right(spans[0], pos) - 1
            if i < 0 or pos > spans[1][i]:
                continue
            ref = fields[3].upper()
            alts = tuple(a for a in fields[4].upper().split(',') if a and set(a) <= set('ACGTN'))
            if alts and set(ref) <= set('ACGTN'):
                found[contig].append((pos, ref, alts))
    for items in found.values():
        items.sort()
    return found, records


def variant_changes_sites(enzyme, amplicon: str, q: int, ref: str, alt: str, resolution_bp: float) -> bool:
    """Whether replacing ``amplicon[q:q + len(ref)]`` by ``alt`` may change the digest.

    Substitutions: compared exactly through the site matches covering the
    change (the remaining matches are identical). Indels: conservatively true
    if a site touches the change or the length shift reaches the resolution.
    """
    size = enzyme.size
    variant = amplicon[:q] + alt + amplicon[q + len(ref):]
    lo = max(0, q - size + 2)
    ref_hits = [(m.start(), m.lastgroup) for m in enzyme.compsite.finditer(' ' + amplicon, lo, q + len(ref) + size)]
    var_hits = [(m.start(), m.lastgroup) for m in enzyme.compsite.finditer(' ' + variant, lo, q + len(alt) + size)]
    if len(ref) == len(alt):
        return ref_hits != var_hits
    return bool(ref_hits or var_hits) or abs(len(alt) - len(ref)) >= resolution_bp


IUPAC_BASES = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT',
               'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
DCAPS_MAX_DESIGNS = 40   # mismatch primers tried per SNP (fewest mismatches first)
DCAPS_MAX_PRIMER = 36


def _site_classes(enzyme) -> list:
    """Allowed bases of each site position, for both strands of a non-palindromic site."""
    site = enzyme_site_string(enzyme).upper()
    if not site or set(site) - set(IUPAC_BASES):
        return []
    forms = {site, site.translate(_IUPAC_COMPLEMENT)[::-1]}
    return [[IUPAC_BASES[c] for c in form] for form in sorted(forms)]


def dcaps_mismatch_sets(enzymes, seq: str, off: int, ref: str, alt: str, max_mismatches: int) -> dict:
    """Primer mismatches next to the SNP that create a site in exactly one allele.

    The mismatch primer ends next to the SNP: a ``left`` primer covers
    ``seq[:off]`` (3' end at ``off - 1``), a ``right`` primer binds
    ``seq[off + 1:]`` (3' end at ``off + 1``). A created site covers the SNP;
    its other bases come from the primer, where up to ``max_mismatches`` bases
    may change (never the 3'-terminal one), or from the template on the other
    side of the SNP, which cannot change. Returns
    ``{(side, ((pos, base), ...)): [enzyme, ...]}`` with top-strand positions.
    """
    found = {}
    for enz in enzymes:
        for classes in _site_classes(enz):
            size = len(classes)
            for s in range(max(0, off - size + 1), min(off, len(seq) - size) + 1):
                snp_class = classes[off - s]
                if (ref in snp_class) == (alt in snp_class):
                    continue
                for side, terminal in (('left', off - 1), ('right', off + 1)):
                    changes = []
                    for i, allowed in enumerate(classes):
                        j = s + i
                        if j == off or seq[j] in allowed:
                            continue
                        if (j > off if side == 'left' else j < off) or j == terminal or seq[j] not in 'ACGT':
                            break
                        changes.append((j, allowed[0]))
                    else:
                        if 0 < len(changes) <= max_mismatches:
                            listed = found.setdefault((side, tuple(changes)), [])
                            if enz not in listed:
                                listed.append(enz)
    return found


def dcaps_primer(template: str, off: int, side: str, tm_min: float, tm_max: float, blocked=()):
    """Longest primer ending next to the SNP with Tm in range: ``(start0, sequence)`` or None.

    A longer primer moves the created cut further from the product end and so
    enlarges the diagnostic difference. ``blocked`` lists ``(start0, length)``
    regions (known variants) the primer may not overlap.
    """
    best = None
    for length in range(PRIMER_MIN_SIZE, DCAPS_MAX_PRIMER + 1):
        lo, hi = (off - length, off) if side == 'left' else (off + 1, off + 1 + length)
        if lo < 0 or hi > len(template) or set(template[lo:hi].upper()) - set('ACGT') \
                or any(q < hi and lo < q + n for q, n in blocked):
            break
        primer = template[lo:hi].upper() if side == 'left' else revcomp(template[lo:hi].upper())
        if tm_min <= primer3.bindings.calc_tm(primer) <= tm_max:
            best = lo, primer
    return best


def describe_mismatches(seq: str, off: int, side: str, changes) -> str:
    """Mismatches as primer bases counted from its 3' end, e.g. "left 3'-3 T>C"."""
    comp = str.maketrans('ACGT', 'TGCA')
    parts = []
    for j, base in sorted(changes, key=lambda change: abs(change[0] - off)):
        old, new = (seq[j], base) if side == 'left' else (seq[j].translate(comp), base.translate(comp))
        parts.append(f"{side} 3'-{abs(off - j)} {old}>{new}")
    return '; '.join(parts)


def diagnostic_delta(fr1: List[int], fr2: List[int]) -> int:
    """Symmetric nearest-band distance used for the Δ threshold."""
    if not fr1 or not fr2:
        return 0
    return int(max(
        max(min(abs(a - b) for b in fr2) for a in fr1),
        max(min(abs(b - a) for a in fr1) for b in fr2),
    ))


def enzyme_site_string(enzyme) -> str:
    site = getattr(enzyme, "site", None)
    if isinstance(site, str):
        return site
    try:
        return str(site) if site is not None else ""
    except Exception:
        return ""


_IUPAC_COMPLEMENT = str.maketrans('ACGTRYKMSWBDHVN', 'TGCAYRMKSWVHDBN')


def site_key(enzyme) -> str:
    """Recognition site independent of the strand it is written on."""
    site = enzyme_site_string(enzyme).upper()
    return min(site, site.translate(_IUPAC_COMPLEMENT)[::-1])


def group_by_site(hits: list) -> list:
    """Keep one hit per recognition site; the others are listed as isoschizomers.

    ``hits`` carry the enzyme class in ``_enz``. A site is represented by a
    commercially available enzyme when there is one (a neoschizomer cutting
    a base further may score marginally higher), then by the best score, then
    by the number of suppliers. The result is in rank order.
    """
    groups = {}
    for hit in hits:
        groups.setdefault((hit.get('design'), site_key(hit['_enz'])), []).append(hit)
    grouped = []
    for members in groups.values():
        keep = min(members, key=lambda h: (not h['_enz'].suppl, rank_key(h)[:4],
                                           -len(h['_enz'].suppl), rank_key(h)[4:]))
        others = sorted({h['enzyme'] for h in members} - {keep['enzyme']})
        grouped.append(dict(keep, isoschizomers='; '.join(others)))
    return sorted(grouped, key=rank_key)


def enzyme_suppliers_human(enzyme) -> str:
    names = set()
    try:
        for v in enzyme.supplier_list():
            if v:
                names.add(str(v).strip())
    except Exception:
        pass
    if not names:
        codes = enzyme_supplier_codes(enzyme)
        for c in sorted(codes):
            aliases = _SUPPLIER_CANON.get(c, [])
            if aliases:
                nm = aliases[0]
                names.add(" ".join(w.capitalize() for w in nm.split()))
            else:
                names.add(c)
    return ", ".join(sorted(names))


# ====== primer3 обёртка ======
PRIMER_MIN_SIZE = 18  # Primer3 default; stated explicitly when a wider target relies on it


def primer3_pick(
    amplicon: str,
    snp_off: int,
    prod_min: int,
    prod_max: int,
    tm_min: float,
    tm_max: float,
    num_return: int = 1,
    target_flank: int = 0,
    lowercase_masking: bool = False,
    excluded_regions: Optional[list] = None,
    fixed_left: Optional[str] = None,
    fixed_right: Optional[str] = None,
) -> Optional[dict]:
    """Pick PCR primers flanking ``[snp_off - target_flank, snp_off + target_flank]``.

    With ``lowercase_masking`` the soft-masked (lower-case) bases of ``amplicon``
    may not be at a primer's 3' end. Primers never overlap ``excluded_regions``
    (``[start0, length]`` pairs, e.g. known variants). A ``fixed_left`` or
    ``fixed_right`` primer (5'->3', present in ``amplicon``) is kept and only
    its partner is designed.
    """
    primer3_pick.last_error = ""
    if not _PRIMER3_AVAILABLE:
        primer3_pick.last_error = "Модуль Primer3 недоступен"
        return None
    seq_args = {"SEQUENCE_TEMPLATE": amplicon,
                "SEQUENCE_TARGET": [snp_off - target_flank, 2 * target_flank + 1]}
    if excluded_regions:
        seq_args["SEQUENCE_EXCLUDED_REGION"] = [list(region) for region in excluded_regions]
    global_args = {
        "PRIMER_TASK": "pick_pcr_primers",
        "PRIMER_NUM_RETURN": num_return,
        "PRIMER_MIN_TM": tm_min,
        "PRIMER_MAX_TM": tm_max,
        "PRIMER_OPT_TM": (tm_min + tm_max) / 2.0,
        "PRIMER_PRODUCT_SIZE_RANGE": [[prod_min, prod_max]],
        "PRIMER_EXPLAIN_FLAG": 1,
    }
    if target_flank:
        global_args["PRIMER_MIN_SIZE"] = PRIMER_MIN_SIZE
    if lowercase_masking:
        global_args["PRIMER_LOWERCASE_MASKING"] = 1
    fixed = fixed_left or fixed_right
    if fixed:
        seq_args["SEQUENCE_PRIMER" if fixed_left else "SEQUENCE_PRIMER_REVCOMP"] = fixed
        global_args["PRIMER_PICK_ANYWAY"] = 1
        global_args["PRIMER_MAX_SIZE"] = max(27, len(fixed))
    try:
        if _PRIMER3_USE_NEW:
            res = primer3.bindings.design_primers(seq_args, global_args)
        else:
            res = primer3.bindings.designPrimers(seq_args, global_args)
        if res.get("PRIMER_PAIR_NUM_RETURNED", 0) > 0:
            return res
        return None
    except Exception as exc:
        primer3_pick.last_error = str(exc)
        return None


primer3_pick.last_error = ""


BASE_HEADER = [
    'variant', 'mapped_id', 'enzyme', 'site', 'pattern', 'frags_ref', 'frags_alt',
    'diag_delta_bp', 'amplicon_start', 'amplicon_end', 'amplicon_len',
    'snp_offset_in_amplicon', 'primers_left', 'primers_right', 'tm_left', 'tm_right',
    'product_size', 'primer3_size', 'primer_left_start', 'primer_left_len',
    'primer_right_start', 'primer_right_len', 'suppliers', 'status', 'reason',
    'delta_threshold_bp',
]
ASSAY_HEADER = [
    'assay_rank', 'primer_pair_index', 'primer_pair_penalty', 'pairs_evaluated',
    'genotype_quality', 'genotype_margin', 'worst_margin', 'weakest_genotypes',
    'bands_ref_ref', 'bands_ref_alt', 'bands_alt_alt', 'worst_bands_ref_ref',
    'worst_bands_ref_alt', 'worst_bands_alt_alt', 'visible_band_count',
    'hidden_ref_count', 'hidden_alt_count', 'search_mode', 'scoring_model',
]
# Added after v1.1.0; appended so that the v1.1.0 columns keep their order.
V12_HEADER = ['isoschizomers', 'primer_site', 'repeat_pct', 'background_variants',
              'background_site_variants', 'design', 'primer_mismatches']


def _primer_candidates(pres, sequence, snp_off):
    """Validate Primer3 coordinates and retain distinct primer pairs."""
    seen = set()
    for i in range(int(pres.get('PRIMER_PAIR_NUM_RETURNED', 0))):
        left = pres[f'PRIMER_LEFT_{i}_SEQUENCE'].upper()
        right = pres[f'PRIMER_RIGHT_{i}_SEQUENCE'].upper()
        ls, ll = map(int, pres[f'PRIMER_LEFT_{i}'])
        rs, rl = map(int, pres[f'PRIMER_RIGHT_{i}'])
        size = int(pres[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'])
        if not (0 <= ls < ls + ll <= snp_off < rs - rl + 1 <= rs < len(sequence)):
            raise ValueError('Primer3 вернул праймеры вне окна или перекрывающие SNP')
        if size != rs - ls + 1 or len(left) != ll or len(right) != rl:
            raise ValueError('Несогласованные размеры продукта/праймеров Primer3')
        if sequence[ls:ls + ll] != left or str(Seq(sequence[rs - rl + 1:rs + 1]).reverse_complement()) != right:
            raise ValueError('Последовательности праймеров не совпадают с референсом')
        key = (ls, ll, rs, rl, left, right)
        if key in seen:
            continue
        seen.add(key)
        penalty = pres.get(f'PRIMER_PAIR_{i}_PENALTY')
        penalty = float(penalty) if penalty is not None else None
        tl, tr = float(pres[f'PRIMER_LEFT_{i}_TM']), float(pres[f'PRIMER_RIGHT_{i}_TM'])
        if not all(math.isfinite(x) for x in (tl, tr)) or (penalty is not None and not math.isfinite(penalty)):
            raise ValueError('Primer3 вернул некорректную оценку')
        yield dict(start0=ls, end0=rs + 1, primers_left=left, primers_right=right,
                   tm_left=tl, tm_right=tr, primer3_size=size,
                   primer_left_start=ls, primer_left_len=ll,
                   primer_right_start=rs, primer_right_len=rl,
                   primer_pair_index=i + 1, primer_pair_penalty=penalty)


def run_rflp_gui_mode(params, status_cb=None, progress_cb=None, cancel_cb=None, report=None):
    params = dict(params)
    for key, value in dict(joint_design=True, primer_pairs=10, top_results=10,
                           gel_resolution_bp=10.0, gel_resolution_pct=3.0,
                           gel_stress_factor=1.5, commercial_only=False,
                           group_isoschizomers=False, snp_end_distance=0,
                           exclude_primer_sites=False, mask_repeats=False,
                           background_vcf='', exclude_background_sites=False, dcaps='off',
                           dcaps_mismatches=1, dcaps_tm_max=70.0, dcaps_prod_min=120, dcaps_prod_max=300,
                           dcaps_delta=20).items():
        params.setdefault(key, value)
    errors = validate_params(params)
    if errors:
        raise ValueError('; '.join(errors))
    model = gel_model_from_params(params)
    report = report if report is not None else {}
    report['parameters'] = dict(params)
    report['scoring_model'] = MODEL_VERSION
    counters = report.setdefault('summary', {})
    counters.update(pairs_evaluated=0, enzyme_evaluations=0, enzyme_errors=0,
                    ambiguous_templates=0, genotype_rejected=0, candidates=0,
                    variants_with_candidates=0, variants_with_robust_candidates=0,
                    candidate_sites=0, pairs_near_end=0, excluded_primer_site=0,
                    background_variants_loaded=0, background_ref_mismatch=0,
                    excluded_background_site=0, dcaps_snps=0, dcaps_designs=0,
                    dcaps_no_primer=0, dcaps_pairs=0, dcaps_candidates=0)

    def check_cancel():
        if cancel_cb and cancel_cb():
            report['cancelled'] = True
            raise AnalysisCancelled()

    cleaned, _, _, _ = normalize_variants_for_rflp(
        fasta_path=params['fasta'], input_path=params['input_file'],
        names2_path=params.get('names2_file', ''), block_bp=100_000,
        progress_every=100000, status_cb=status_cb,
        skip_norm=params.get('skip_norm', False),
        save_files=not params.get('no_save_norm', True), report=report,
        progress_cb=progress_cb, cancel_cb=cancel_cb,
    )
    header = BASE_HEADER + ASSAY_HEADER + V12_HEADER
    out_rows = []
    if not cleaned:
        counters['rows'] = 0
        return header, out_rows
    use_primer3 = params.get('use_primer3', False)
    joint = bool(params['joint_design'] and use_primer3)
    mode = 'joint' if joint else 'single_pair' if use_primer3 else 'window'
    codes = parse_suppliers_arg(','.join(params.get('suppliers', [])))
    enzymes = sorted((e for e in AllEnzymes if not codes or enzyme_supplier_codes(e) & codes),
                     key=lambda e: e.__name__)
    if params['commercial_only']:
        # REBASE lists at least one commercial supplier.
        enzymes = [e for e in enzymes if e.suppl]
    names = load_names2(params.get('names2_file', ''))
    if status_cb:
        status_cb(f'[INFO] Поиск: {mode}; SNP: {len(cleaned)}; ферментов: {len(enzymes)}')
        status_cb('[INFO] Оценка полос: модель присутствия, полная рестрикция; '
                  'яркость и специфичность ПЦР по всему геному не оцениваются.')
        options = [label for key, label in (
            ('commercial_only', 'только продаваемые ферменты'),
            ('group_isoschizomers', 'изошизомеры объединены (одна строка на сайт)'),
            ('exclude_primer_sites', 'без сайта фермента в праймерах'),
            ('mask_repeats', 'повторы (строчные буквы FASTA) не на 3′-концах праймеров'),
            ('exclude_background_sites', 'без фоновых вариантов, меняющих разрезы'),
        ) if params[key] and (key != 'exclude_background_sites' or params['background_vcf'])]
        if use_primer3 and params['snp_end_distance']:
            options.append(f"SNP не ближе {params['snp_end_distance']} п.н. к концам ампликона")
        if use_primer3 and params['dcaps'] != 'off':
            when = 'если нет robust-схемы' if params['dcaps'] == 'fallback' else 'для каждого SNP'
            options.append(f"dCAPS {when}, несовпадений до {params['dcaps_mismatches']}")
        if options:
            status_cb('[INFO] Отбор: ' + '; '.join(options))

    def append_record(record):
        out_rows.append([record.get(key) for key in header])

    background = None
    if params['background_vcf']:
        regions = {}
        for chrom, pos, _, _ in cleaned:
            regions.setdefault(names.get(chrom, chrom), []).append(
                (pos - params['flank'], pos + params['flank']))
        background, records = load_background_variants(params['background_vcf'], regions, names, cancel_cb)
        counters['background_variants_loaded'] = sum(map(len, background.values()))
        if status_cb:
            status_cb(f"[INFO] Фоновые варианты: {counters['background_variants_loaded']} в окнах SNP "
                      f"из {records} записей VCF")

    dcaps = params['dcaps'] if use_primer3 and _PRIMER3_AVAILABLE else 'off'
    # Without a known cut position the fragment lengths are unknown.
    dcaps_enzymes = [e for e in enzymes if not e.is_unknown()]

    # Always release the FASTA handle, including cancellation and exceptions.
    fa = Fasta(params['fasta'], as_raw=True)
    masked_windows = 0
    try:
        for vi, (chrom, pos, ref, alt) in enumerate(cleaned, 1):
            check_cancel()
            mapped_id = names.get(chrom, chrom)
            start = max(1, pos - params['flank'])
            end = min(len(fa[mapped_id]), pos + params['flank'])
            # Soft-masked FASTA: lower case marks repeats; the analysis uses upper case.
            raw = str(fa[mapped_id][start - 1:end])
            seq = raw.upper()
            masked_windows += raw != seq
            off = pos - start
            base = dict(variant=f'{chrom}:{pos} {ref}>{alt}', mapped_id=mapped_id,
                        enzyme='', site='', pattern='no_enzyme_found',
                        amplicon_start=start, amplicon_end=end, amplicon_len=len(seq),
                        snp_offset_in_amplicon=off, delta_threshold_bp=params['delta'],
                        search_mode=mode, scoring_model=MODEL_VERSION, pairs_evaluated=0)
            if not 0 <= off < len(seq) or seq[off] != ref:
                append_record(dict(base, status='analysis_error', reason='REF не совпадает с FASTA'))
                continue
            altseq = seq[:off] + alt + seq[off + 1:]
            # Known variants of the window (offsets), except a record of the SNP itself.
            known = []
            if background is not None:
                items = background.get(mapped_id, [])
                for p, r, alts in items[bisect.bisect_left(items, (start,)):bisect.bisect_left(items, (end + 1,))]:
                    q = p - start
                    if p + len(r) - 1 > end or p <= pos < p + len(r):
                        continue
                    if seq[q:q + len(r)] != r:
                        counters['background_ref_mismatch'] += 1
                        continue
                    known.append((q, r, alts))
            # The remaining enzymes give equal REF/ALT fragments in every amplicon,
            # i.e. Δ = 0, which is below any threshold (delta >= 1).
            site_enzymes = snp_site_enzymes(enzymes, seq, altseq, off)
            mask = params['mask_repeats']
            excluded = [(q, len(r)) for q, r, _ in known]
            end_distance = params['snp_end_distance'] if use_primer3 else 0
            failure, pairs = None, []
            if end_distance and not end_distance <= off < len(seq) - end_distance:
                failure = dict(base, status='no_primers',
                               reason=f'SNP ближе {end_distance} п.н. к краю окна; увеличьте фланг')
            elif use_primer3:
                check_cancel()
                pres = primer3_pick(raw if mask else seq, off, params['prod_min'], params['prod_max'],
                                    params['tm_min'], params['tm_max'],
                                    num_return=params['primer_pairs'] if joint else 1,
                                    target_flank=max(0, end_distance - PRIMER_MIN_SIZE),
                                    lowercase_masking=mask, excluded_regions=excluded)
                if not pres:
                    error = primer3_pick.last_error
                    failure = dict(base, status='primer3_error' if error else 'no_primers',
                                   reason=error or 'Primer3 не вернул подходящих пар')
                else:
                    try:
                        pairs = list(_primer_candidates(pres, seq, off))
                    except (KeyError, TypeError, ValueError) as exc:
                        failure = dict(base, status='primer3_error', reason=str(exc))
                    else:
                        # Guaranteed by the Primer3 target; kept as an explicit check.
                        near = [p for p in pairs if min(off - p['start0'], p['end0'] - 1 - off) < end_distance]
                        counters['pairs_near_end'] += len(near)
                        pairs = [p for p in pairs if p not in near]
            else:
                pairs = [dict(start0=0, end0=len(seq))]
            if failure is not None and dcaps == 'off':
                append_record(failure)
                continue
            base['pairs_evaluated'] = len(pairs) if use_primer3 else 0
            counters['pairs_evaluated'] += base['pairs_evaluated']
            hits = []
            local = dict(errors=0, ambiguous=0, genotype_rejected=0)

            def amplicon_record(pair, a, b, **extra):
                record = dict(base, **{k: v for k, v in pair.items() if k not in ('start0', 'end0')})
                record.update(amplicon_start=start + a, amplicon_end=start + b - 1,
                              amplicon_len=b - a, product_size=b - a,
                              snp_offset_in_amplicon=off - a,
                              repeat_pct=round(100 * sum(c.islower() for c in raw[a:b]) / (b - a), 1),
                              **extra)
                return record

            def evaluate(enz, refamp, altamp, record, a, b):
                """The assay of one enzyme on one amplicon, or None if it fails a filter."""
                try:
                    cr, ca = cut_positions(enz, refamp), cut_positions(enz, altamp)
                except Exception as exc:
                    local['errors'] += 1
                    if local['errors'] <= 3 and status_cb:
                        status_cb(f'[WARN] {base["variant"]} / {enz.__name__}: {exc}')
                    return None
                if max(len(cr), len(ca)) - 2 > params['max_cuts']:
                    return None
                fr, fal = frag_lengths(cr), frag_lengths(ca)
                if params['gain_loss_only'] and len(cr) == len(ca):
                    return None
                # Legacy baseline requires every fragment in range. Joint
                # search treats the same range as a visibility window.
                if not joint and any(x < params['min_frag'] or x > params['max_frag'] for x in fr + fal):
                    return None
                visible_r = [x for x in fr if model.min_visible <= x <= model.max_visible]
                visible_a = [x for x in fal if model.min_visible <= x <= model.max_visible]
                actual_delta = diagnostic_delta(visible_r, visible_a)
                if actual_delta < params['delta']:
                    return None
                scores = score_genotypes(fr, fal, model)
                if joint and scores['genotype_quality'] == 'ambiguous':
                    local['genotype_rejected'] += 1
                    return None
                reason = {
                    'robust': 'Три генотипа различимы во всём заданном диапазоне разрешения',
                    'nominal_only': 'Три генотипа различимы только при исходном разрешении',
                    'ambiguous': 'Три генотипа неразличимы по присутствию полос',
                }[scores['genotype_quality']]
                if not use_primer3:
                    reason = 'Оценка окна без подобранных праймеров. ' + reason
                dcaps_design = record.get('design') == 'dCAPS'
                if dcaps_design:
                    allele = 'REF' if len(cr) > len(ca) else 'ALT'
                    reason = f'dCAPS: несовпадение в праймере создаёт сайт в аллеле {allele}. ' + reason
                primer_site = primer_site_label(enz, refamp, record.get('primer_left_len'),
                                                record.get('primer_right_len'))
                # In a mismatch primer such a site is a design error: every allele is cut.
                if primer_site and (params['exclude_primer_sites'] or dcaps_design):
                    counters['excluded_primer_site'] += 1
                    return None
                inside = [(q - a, r, alts) for q, r, alts in known if a <= q and q + len(r) <= b]
                changing = sorted({f'{start + a + q} {r}>{x}' for q, r, alts in inside for x in alts
                                   if any(variant_changes_sites(enz, amp, q, r, x, model.resolution_bp)
                                          for amp in (refamp, altamp))},
                                  key=lambda text: int(text.split()[0]))
                if changing and params['exclude_background_sites']:
                    counters['excluded_background_site'] += 1
                    return None
                return dict(record, **scores, enzyme=enz.__name__, _enz=enz,
                            site=enzyme_site_string(enz),
                            pattern='gain/loss' if len(cr) != len(ca) else 'shift',
                            frags_ref=frags_to_str(fr), frags_alt=frags_to_str(fal),
                            diag_delta_bp=actual_delta, suppliers=enzyme_suppliers_human(enz),
                            primer_site=primer_site,
                            background_variants=len(inside) if background is not None else None,
                            background_site_variants='; '.join(changing),
                            status='ok' if use_primer3 else 'window_only', reason=reason)

            for pair in pairs:
                check_cancel()
                a, b = pair['start0'], pair['end0']
                refamp, altamp = seq[a:b], altseq[a:b]
                if set(refamp + altamp) - set('ACGT'):
                    local['ambiguous'] += 1
                    continue
                record = amplicon_record(pair, a, b, design='natural')
                counters['enzyme_evaluations'] += len(enzymes)
                for enz in site_enzymes:
                    check_cancel()
                    hit = evaluate(enz, refamp, altamp, record, a, b)
                    if hit is not None:
                        hits.append(hit)

            if dcaps == 'always' or (dcaps == 'fallback'
                                     and not any(h['genotype_quality'] == 'robust' for h in hits)):
                counters['dcaps_snps'] += 1
                designs = dcaps_mismatch_sets(dcaps_enzymes, seq, off, ref, alt, params['dcaps_mismatches'])
                ordered = sorted(designs.items(), key=lambda item: (len(item[0][1]), -len(item[1]), item[0]))
                for (side, changes), design_enzymes in ordered[:DCAPS_MAX_DESIGNS]:
                    check_cancel()
                    template = list(raw if mask else seq)
                    for j, new in changes:
                        template[j] = new
                    template = ''.join(template)
                    terminal = off - 1 if side == 'left' else off + 1
                    primer = None if mask and raw[terminal].islower() else dcaps_primer(
                        template, off, side, params['tm_min'], params['tm_max'], excluded)
                    if primer is None:
                        counters['dcaps_no_primer'] += 1
                        continue
                    counters['dcaps_designs'] += 1
                    pres = primer3_pick(template, off, params['prod_min'], params['prod_max'],
                                        params['tm_min'], params['tm_max'],
                                        num_return=params['primer_pairs'] if joint else 1,
                                        lowercase_masking=mask, excluded_regions=excluded,
                                        fixed_left=primer[1] if side == 'left' else None,
                                        fixed_right=primer[1] if side == 'right' else None)
                    mutated = template.upper()
                    try:
                        dpairs = list(_primer_candidates(pres, mutated, off)) if pres else []
                    except (KeyError, TypeError, ValueError):
                        dpairs = []
                    counters['dcaps_pairs'] += len(dpairs)
                    mutated_alt = mutated[:off] + alt + mutated[off + 1:]
                    lowered = {j for j, _ in changes}
                    for pair in dpairs:
                        a, b = pair['start0'], pair['end0']
                        refamp, altamp = mutated[a:b], mutated_alt[a:b]
                        if set(refamp + altamp) - set('ACGT'):
                            continue
                        # Mismatched bases in lower case, as usual for dCAPS primers.
                        left = ''.join(c.lower() if a + i in lowered else c
                                       for i, c in enumerate(pair['primers_left']))
                        right = ''.join(c.lower() if b - 1 - i in lowered else c
                                        for i, c in enumerate(pair['primers_right']))
                        record = amplicon_record(pair, a, b, design='dCAPS', primers_left=left,
                                                 primers_right=right,
                                                 primer_mismatches=describe_mismatches(seq, off, side, changes))
                        for enz in design_enzymes:
                            hit = evaluate(enz, refamp, altamp, record, a, b)
                            if hit is not None:
                                hits.append(hit)
                                counters['dcaps_candidates'] += 1

            counters['enzyme_errors'] += local['errors']
            counters['ambiguous_templates'] += local['ambiguous']
            counters['genotype_rejected'] += local['genotype_rejected']
            hits.sort(key=rank_key)
            counters['candidates'] += len(hits)
            counters['candidate_sites'] += len({site_key(h['_enz']) for h in hits})
            if params['group_isoschizomers']:
                hits = group_by_site(hits)
            counters['variants_with_candidates'] += bool(hits)
            counters['variants_with_robust_candidates'] += any(h['genotype_quality'] == 'robust' for h in hits)
            # Limit only joint search. The one-pair/window baselines retain all hits.
            limit = params['top_results'] if joint else 0
            for rank, hit in enumerate(hits[:limit] if limit else hits, 1):
                append_record(dict(hit, assay_rank=rank))
            if not hits and failure is not None:
                append_record(failure)
            elif not hits:
                reason = ('Не найдено сочетание для различения трёх генотипов' if joint
                          else 'Подходящий фермент не найден')
                if local['ambiguous']:
                    reason += f"; окон с неопределёнными основаниями: {local['ambiguous']}"
                if local['errors']:
                    reason += f"; ошибок расчёта ферментов: {local['errors']}"
                append_record(dict(base, status='no_discriminating_assay' if joint else 'no_enzyme_found', reason=reason))
            if status_cb:
                status_cb(f'[INFO] {base["variant"]}: пар {base["pairs_evaluated"]}, '
                          f"кандидатов {len(hits)}, неразличимых {local['genotype_rejected']}")
            if progress_cb:
                progress_cb(vi, len(cleaned))
    finally:
        fa.close()
    if params['mask_repeats'] and use_primer3 and not masked_windows and status_cb:
        status_cb('[WARN] В окнах SNP нет строчных букв: FASTA без soft-mask, '
                  'маскирование повторов не действует.')
    counters['rows'] = len(out_rows)
    return header, out_rows
