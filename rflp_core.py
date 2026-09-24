import sys
import io
import csv
import os
import csv
import re
import traceback
import logging
import math
from typing import List, Dict, Tuple, Optional
from collections import namedtuple
from datetime import datetime

from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio.Restriction import AllEnzymes
from assay_scoring import GelModel, MODEL_VERSION, score_genotypes, rank_key


class AnalysisCancelled(Exception):
    """Raised when a running analysis is cancelled by the caller."""


def validate_params(params: dict) -> list[str]:
    errors = []
    required = ('fasta', 'input_file')
    for key in required:
        if not str(params.get(key, '')).strip():
            errors.append(f"Не задан параметр: {key}")
    for key in ('fasta', 'input_file'):
        value = str(params.get(key, '')).strip()
        if value and not os.path.isfile(value):
            errors.append(f"Файл не найден ({key}): {value}")
    for lo, hi, label in (('min_frag', 'max_frag', 'размер фрагмента'),
                          ('prod_min', 'prod_max', 'размер продукта'),
                          ('tm_min', 'tm_max', 'Tm')):
        try:
            lower, upper = float(params.get(lo)), float(params.get(hi))
            if not all(math.isfinite(x) and x > 0 for x in (lower, upper)):
                errors.append(f"Некорректный параметр: {lo}/{hi}")
            elif lower > upper:
                errors.append(f"Минимальный {label} больше максимального")
        except (TypeError, ValueError):
            errors.append(f"Некорректный параметр: {lo}/{hi}")
    for key, default, low, high in (
        ('primer_pairs', 10, 1, 100), ('top_results', 10, 0, 10000),
        ('flank', 300, 1, 100000), ('max_cuts', 3, 0, 100),
        ('delta', 25, 1, 100000),
    ):
        value = params.get(key, default)
        if isinstance(value, bool) or not isinstance(value, int) or not low <= value <= high:
            errors.append(f"{key}: требуется целое число от {low} до {high}")
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


Variant = namedtuple("Variant", "chrom pos ref alt mapped_id")


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


def parse_input(path: str, report: Optional[dict] = None) -> List[Tuple[str, int, str, str]]:
    rows: List[Tuple[str, int, str, str]] = []
    try:
        if not os.path.isfile(path):
            logger.error(f"Файл с вариантами не найден: {path}")
            return rows
        with open(path, "r", encoding="utf-8-sig", errors="replace", newline="") as f:
            first = f.readline()
            if not first:
                return rows
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

            i_chrom = idx("chrom", "chr", "chr1", "chromosome")
            i_pos = idx("pos", "position")
            i_ref = idx("ref", "reference")
            i_alt = idx("alt", "alternate", "alt_allele")

            def parse_line(parts, line_no):
                if len(parts) < 4:
                    if report is not None:
                        report.setdefault('rejected', []).append({'line': line_no, 'reason': 'parse_error', 'details': 'Ожидалось минимум 4 поля', 'raw': delim.join(parts)})
                    return None
                try:
                    ch = parts[i_chrom] if i_chrom is not None else parts[0]
                    po = int(parts[i_pos] if i_pos is not None else parts[1])
                    rf = parts[i_ref] if i_ref is not None else parts[2]
                    al = parts[i_alt] if i_alt is not None else parts[3]
                    return (
                        ch.replace("chr", "").strip(),
                        po,
                        rf.strip().upper(),
                        al.strip().upper(),
                    )
                except Exception as e2:
                    logger.error(f"Ошибка парсинга строки: {e2}")
                    if report is not None:
                        report.setdefault('rejected', []).append({'line': line_no, 'reason': 'parse_error', 'details': str(e2), 'raw': delim.join(parts)})
                    return None

            data: List[Tuple[str, int, str, str]] = []
            has_header = any(
                x in hdr
                for x in ("chrom", "chr", "pos", "ref", "alt", "position", "reference", "alternate")
            )
            if has_header:
                for line_no, parts in enumerate(reader, start=2):
                    if not parts or not any(x.strip() for x in parts):
                        continue
                    if "должно быть:" in delim.join(parts).lower():
                        continue
                    rec = parse_line(parts, line_no)
                    if rec:
                        data.append(rec)
            else:
                rec = parse_line(first_fields, 1)
                if rec:
                    data.append(rec)
                for line_no, parts in enumerate(reader, start=2):
                    if not parts or not any(x.strip() for x in parts):
                        continue
                    rec = parse_line(parts, line_no)
                    if rec:
                        data.append(rec)
            rows.extend(data)
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

    def __init__(self, fasta_path: str, block_bp: int = 2_000_000):
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
    block_bp: int = 2_000_000,
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


def any_diag_delta(fr1: List[int], fr2: List[int], delta: int) -> bool:
    s1 = sorted(fr1)
    s2 = sorted(fr2)
    if len(s1) != len(s2):
        return True
    return any(abs(a - b) >= delta for a, b in zip(s1, s2))


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
def primer3_pick(
    amplicon: str,
    snp_off: int,
    prod_min: int,
    prod_max: int,
    tm_min: float,
    tm_max: float,
    num_return: int = 1,
) -> Optional[dict]:
    primer3_pick.last_error = ""
    if not _PRIMER3_AVAILABLE:
        primer3_pick.last_error = "Модуль Primer3 недоступен"
        return None
    seq_args = {"SEQUENCE_TEMPLATE": amplicon, "SEQUENCE_TARGET": [snp_off, 1]}
    global_args = {
        "PRIMER_TASK": "pick_pcr_primers",
        "PRIMER_NUM_RETURN": num_return,
        "PRIMER_MIN_TM": tm_min,
        "PRIMER_MAX_TM": tm_max,
        "PRIMER_OPT_TM": (tm_min + tm_max) / 2.0,
        "PRIMER_PRODUCT_SIZE_RANGE": [[prod_min, prod_max]],
        "PRIMER_EXPLAIN_FLAG": 1,
    }
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
                           gel_stress_factor=1.5).items():
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
                    variants_with_candidates=0, variants_with_robust_candidates=0)

    def check_cancel():
        if cancel_cb and cancel_cb():
            report['cancelled'] = True
            raise AnalysisCancelled()

    cleaned, _, _, _ = normalize_variants_for_rflp(
        fasta_path=params['fasta'], input_path=params['input_file'],
        names2_path=params.get('names2_file', ''), block_bp=2_000_000,
        progress_every=100000, status_cb=status_cb,
        skip_norm=params.get('skip_norm', False),
        save_files=not params.get('no_save_norm', True), report=report,
        progress_cb=progress_cb, cancel_cb=cancel_cb,
    )
    header = BASE_HEADER + ASSAY_HEADER
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
    names = load_names2(params.get('names2_file', ''))
    if status_cb:
        status_cb(f'[INFO] Поиск: {mode}; SNP: {len(cleaned)}; ферментов: {len(enzymes)}')
        status_cb('[INFO] Оценка полос: модель присутствия, полная рестрикция; '
                  'яркость и специфичность ПЦР по всему геному не оцениваются.')

    def append_record(record):
        out_rows.append([record.get(key) for key in header])

    # Always release the FASTA handle, including cancellation and exceptions.
    fa = Fasta(params['fasta'], as_raw=True)
    try:
        for vi, (chrom, pos, ref, alt) in enumerate(cleaned, 1):
            check_cancel()
            mapped_id = names.get(chrom, chrom)
            start = max(1, pos - params['flank'])
            end = min(len(fa[mapped_id]), pos + params['flank'])
            seq = str(fa[mapped_id][start - 1:end]).upper()
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
            if use_primer3:
                check_cancel()
                pres = primer3_pick(seq, off, params['prod_min'], params['prod_max'],
                                    params['tm_min'], params['tm_max'],
                                    num_return=params['primer_pairs'] if joint else 1)
                if not pres:
                    error = primer3_pick.last_error
                    append_record(dict(base, status='primer3_error' if error else 'no_primers',
                                       reason=error or 'Primer3 не вернул подходящих пар'))
                    continue
                try:
                    pairs = list(_primer_candidates(pres, seq, off))
                except (KeyError, TypeError, ValueError) as exc:
                    append_record(dict(base, status='primer3_error', reason=str(exc)))
                    continue
            else:
                pairs = [dict(start0=0, end0=len(seq))]
            base['pairs_evaluated'] = len(pairs) if use_primer3 else 0
            counters['pairs_evaluated'] += base['pairs_evaluated']
            hits = []
            local_errors = local_ambiguous = local_genotype_rejected = 0
            for pair in pairs:
                check_cancel()
                a, b = pair['start0'], pair['end0']
                refamp, altamp = seq[a:b], altseq[a:b]
                if set(refamp + altamp) - set('ACGT'):
                    local_ambiguous += 1
                    continue
                record = dict(base, **{k: v for k, v in pair.items() if k not in ('start0', 'end0')})
                record.update(amplicon_start=start + a, amplicon_end=start + b - 1,
                              amplicon_len=b - a, product_size=b - a,
                              snp_offset_in_amplicon=off - a)
                for enz in enzymes:
                    check_cancel()
                    counters['enzyme_evaluations'] += 1
                    try:
                        cr, ca = cut_positions(enz, refamp), cut_positions(enz, altamp)
                    except Exception as exc:
                        local_errors += 1
                        if local_errors <= 3 and status_cb:
                            status_cb(f'[WARN] {base["variant"]} / {enz.__name__}: {exc}')
                        continue
                    if max(len(cr), len(ca)) - 2 > params['max_cuts']:
                        continue
                    fr, fal = frag_lengths(cr), frag_lengths(ca)
                    if params['gain_loss_only'] and len(cr) == len(ca):
                        continue
                    # Legacy baseline requires every fragment in range. Joint
                    # search treats the same range as a visibility window.
                    if not joint and any(x < params['min_frag'] or x > params['max_frag'] for x in fr + fal):
                        continue
                    visible_r = [x for x in fr if model.min_visible <= x <= model.max_visible]
                    visible_a = [x for x in fal if model.min_visible <= x <= model.max_visible]
                    actual_delta = diagnostic_delta(visible_r, visible_a)
                    if actual_delta < params['delta']:
                        continue
                    scores = score_genotypes(fr, fal, model)
                    if joint and scores['genotype_quality'] == 'ambiguous':
                        local_genotype_rejected += 1
                        continue
                    reason = {
                        'robust': 'Три генотипа различимы во всём заданном диапазоне разрешения',
                        'nominal_only': 'Три генотипа различимы только при исходном разрешении',
                        'ambiguous': 'Три генотипа неразличимы по присутствию полос',
                    }[scores['genotype_quality']]
                    if not use_primer3:
                        reason = 'Оценка окна без подобранных праймеров. ' + reason
                    hits.append(dict(record, **scores, enzyme=enz.__name__,
                                     site=enzyme_site_string(enz),
                                     pattern='gain/loss' if len(cr) != len(ca) else 'shift',
                                     frags_ref=frags_to_str(fr), frags_alt=frags_to_str(fal),
                                     diag_delta_bp=actual_delta, suppliers=enzyme_suppliers_human(enz),
                                     status='ok' if use_primer3 else 'window_only', reason=reason))
            counters['enzyme_errors'] += local_errors
            counters['ambiguous_templates'] += local_ambiguous
            counters['genotype_rejected'] += local_genotype_rejected
            hits.sort(key=rank_key)
            counters['candidates'] += len(hits)
            counters['variants_with_candidates'] += bool(hits)
            counters['variants_with_robust_candidates'] += any(h['genotype_quality'] == 'robust' for h in hits)
            # Limit only joint search. The one-pair/window baselines retain all hits.
            limit = params['top_results'] if joint else 0
            for rank, hit in enumerate(hits[:limit] if limit else hits, 1):
                append_record(dict(hit, assay_rank=rank))
            if not hits:
                reason = ('Не найдено сочетание для различения трёх генотипов' if joint
                          else 'Подходящий фермент не найден')
                if local_ambiguous:
                    reason += f'; окон с неопределёнными основаниями: {local_ambiguous}'
                if local_errors:
                    reason += f'; ошибок расчёта ферментов: {local_errors}'
                append_record(dict(base, status='no_discriminating_assay' if joint else 'no_enzyme_found', reason=reason))
            if status_cb:
                status_cb(f'[INFO] {base["variant"]}: пар {base["pairs_evaluated"]}, '
                          f'кандидатов {len(hits)}, неразличимых {local_genotype_rejected}')
            if progress_cb:
                progress_cb(vi, len(cleaned))
    finally:
        fa.close()
    counters['rows'] = len(out_rows)
    return header, out_rows
