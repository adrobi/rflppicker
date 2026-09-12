import sys
import io
import csv
import os
import csv
import re
import traceback
import logging
from typing import List, Dict, Tuple, Optional
from collections import namedtuple
from datetime import datetime

from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio.Restriction import AllEnzymes


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
            if float(params.get(lo)) > float(params.get(hi)):
                errors.append(f"Минимальный {label} больше максимального")
        except (TypeError, ValueError):
            errors.append(f"Некорректный параметр: {lo}/{hi}")
    return errors


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
    try:
        sites = enzyme.search(Seq(seq))
    except Exception:
        sites = []
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
) -> Optional[dict]:
    primer3_pick.last_error = ""
    if not _PRIMER3_AVAILABLE:
        primer3_pick.last_error = "Модуль Primer3 недоступен"
        return None
    seq_args = {"SEQUENCE_TEMPLATE": amplicon, "SEQUENCE_TARGET": [snp_off, 1]}
    global_args = {
        "PRIMER_TASK": "pick_pcr_primers",
        "PRIMER_NUM_RETURN": 1,
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


def run_rflp_gui_mode(params, status_cb=None, progress_cb=None, cancel_cb=None, report=None):
    errors = validate_params(params)
    if errors:
        raise ValueError('; '.join(errors))
    flank = params["flank"]
    min_frag = params["min_frag"]
    max_frag = params["max_frag"]
    delta = params["delta"]
    gain_loss_only = params["gain_loss_only"]
    max_cuts = params["max_cuts"]
    suppliers = params["suppliers"]
    fasta = params["fasta"]
    input_file = params["input_file"]
    names2_file = params["names2_file"]
    use_primer3 = params.get("use_primer3", False)
    prod_min = params.get("prod_min", 250)
    prod_max = params.get("prod_max", 600)
    tm_min = params.get("tm_min", 58.0)
    tm_max = params.get("tm_max", 62.0)
    skip_norm = params.get("skip_norm", False)
    no_save_norm = params.get("no_save_norm", True)

    # ==== 1. НОРМАЛИЗАЦИЯ ВАРИАНТОВ ====
    cleaned_variants, summary_lines, cleaned_path, rejected_path = normalize_variants_for_rflp(
        fasta_path=fasta,
        input_path=input_file,
        names2_path=names2_file,
        block_bp=2_000_000,
        progress_every=100000,
        status_cb=status_cb,
        skip_norm=skip_norm,
        save_files=not no_save_norm,
        report=report,
        progress_cb=progress_cb,
        cancel_cb=cancel_cb,
    )

    if status_cb:
        if not skip_norm:
            if cleaned_path:
                status_cb(f"[INFO] NORMALIZE: cleaned сохранён в: {cleaned_path}")
                status_cb(f"[INFO] NORMALIZE: rejected сохранён в: {rejected_path}")
            else:
                status_cb("[INFO] NORMALIZE: файлы cleaned/rejected не сохранялись")
        status_cb(f"[INFO] NORMALIZE: вариантов после шага нормализации: {len(cleaned_variants)}")

    if not cleaned_variants:
        return [], []

    variants_raw = cleaned_variants

    # ==== 2. ЗАГРУЗКА names2 для RFLP ====
    names2 = load_names2(names2_file)

    # ==== 3. Подготовка результатов RFLP ====
    allowed_codes: set[str] = set()
    if suppliers:
        allowed_codes = parse_suppliers_arg(",".join(suppliers))

    try:
        fa = Fasta(fasta, as_raw=True)
    except Exception as e:
        msg = f"\n[ERROR] Не удалось открыть FASTA файл: {fasta}: {e}"
        if status_cb:
            status_cb(msg)
        logger.error(msg)
        return [], []

    out_rows: List[list] = []

    header = [
        "variant",
        "mapped_id",
        "enzyme",
        "site",
        "pattern",
        "frags_ref",
        "frags_alt",
        "diag_delta_bp",
        "amplicon_start",
        "amplicon_end",
        "amplicon_len",
        "snp_offset_in_amplicon",
        "primers_left",
        "primers_right",
        "tm_left",
        "tm_right",
        "product_size",
        "primer3_size",
        "primer_left_start",
        "primer_left_len",
        "primer_right_start",
        "primer_right_len",
        "suppliers",
        "status",
        "reason",
        "delta_threshold_bp",
    ]

    if status_cb:
        status_cb(f"\n[INFO] RFLP: к анализу передано вариантов: {len(variants_raw)}")

    # Enzyme metadata is static; build the filtered list once, outside the SNP loop.
    if allowed_codes:
        enzymes = [e for e in AllEnzymes if enzyme_supplier_codes(e) & allowed_codes]
    else:
        enzymes = list(AllEnzymes)

    for variant_index, (chrom, pos, ref, alt) in enumerate(variants_raw, start=1):
        if cancel_cb and cancel_cb():
            if report is not None:
                report['cancelled'] = True
            raise AnalysisCancelled()
        if progress_cb:
            progress_cb(variant_index, max(1, len(variants_raw)))
        mapped_id = names2.get(chrom, chrom)
        start = max(1, pos - flank)
        end = pos + flank

        try:
            seq_win = str(fa[mapped_id][start - 1 : end])
        except Exception as e:
            msg = f"\n[WARN] Не удалось взять {mapped_id}:{start}-{end}: {e}"
            if status_cb:
                status_cb(msg)
            logger.warning(msg)
            continue

        snp_off = pos - start
        if snp_off < 0 or snp_off >= len(seq_win):
            msg = f"\n[WARN] SNP offset вне окна: {chrom}:{pos}"
            if status_cb:
                status_cb(msg)
            logger.warning(msg)
            continue

        genome_base = seq_win[snp_off].upper()
        if genome_base not in (ref.upper(), alt.upper()):
            msg = (
                f"\n[WARN] REF FASTA не совпал с входом "
                f"для {chrom}:{pos} ({genome_base} vs {ref})"
            )
            if status_cb:
                status_cb(msg)
            logger.warning(msg)

        seq_ref_win = (seq_win[:snp_off] + ref + seq_win[snp_off + 1 :]).upper()
        seq_alt_win = (seq_win[:snp_off] + alt + seq_win[snp_off + 1 :]).upper()

        pres: Optional[dict] = None
        amp_start0 = 0
        amp_end0_excl = len(seq_ref_win)
        amplicon_len = len(seq_ref_win)
        snp_off_amp = snp_off
        primers_left = ""
        primers_right = ""
        tm_left: Optional[float] = None
        tm_right: Optional[float] = None
        p3_size: Optional[int] = None
        pl_start = pl_len = pr_start = pr_len = None
        primer_status = 'window_only' if not use_primer3 else 'no_primers'
        primer_reason = ''

        if use_primer3:
            if not _PRIMER3_AVAILABLE:
                primer_reason = f"Primer3 недоступен: {_PRIMER3_IMPORT_ERROR}"
                primer_status = 'primer3_error'
                out_rows.append([f"{chrom}:{pos} {ref}>{alt}", mapped_id, "", "", "no_enzyme_found", "", "", "", start, end, len(seq_ref_win), snp_off, "", "", None, None, None, None, None, None, None, None, "", primer_status, primer_reason, int(delta)])
                continue
            pres = primer3_pick(
                seq_ref_win,
                snp_off,
                prod_min,
                prod_max,
                tm_min,
                tm_max,
            )
            if pres:
                primers_left = pres["PRIMER_LEFT_0_SEQUENCE"]
                primers_right = pres["PRIMER_RIGHT_0_SEQUENCE"]
                tm_left = float(pres["PRIMER_LEFT_0_TM"])
                tm_right = float(pres["PRIMER_RIGHT_0_TM"])
                pl_start, pl_len = pres["PRIMER_LEFT_0"]
                pr_start, pr_len = pres["PRIMER_RIGHT_0"]
                p3_size = int(pres["PRIMER_PAIR_0_PRODUCT_SIZE"])

                amp_start0 = int(pl_start)
                amp_end0_excl = amp_start0 + int(pres["PRIMER_PAIR_0_PRODUCT_SIZE"])
                amp_start0 = max(0, min(len(seq_ref_win), amp_start0))
                amp_end0_excl = max(amp_start0, min(len(seq_ref_win), amp_end0_excl))

                amplicon_len = amp_end0_excl - amp_start0
                snp_off_amp = snp_off - amp_start0
                primer_status = 'ok'
            else:
                primer_reason = getattr(primer3_pick, 'last_error', '') or 'Primer3 не вернул подходящую пару праймеров'
                if use_primer3:
                    # Do not report restriction candidates as if a PCR product existed.
                    primer_status = 'primer3_error' if getattr(primer3_pick, 'last_error', '') else 'no_primers'
                    row = [f"{chrom}:{pos} {ref}>{alt}", mapped_id, "", "", "no_enzyme_found", "", "", "", start, end, len(seq_ref_win), snp_off, "", "", None, None, None, None, None, None, None, None, "", primer_status, primer_reason, int(delta)]
                    out_rows.append(row)
                    continue

        seq_ref_amp = seq_ref_win[amp_start0:amp_end0_excl]
        seq_alt_amp = seq_alt_win[amp_start0:amp_end0_excl]

        amplicon_start_genome = start + amp_start0
        amplicon_end_genome = start + amp_end0_excl - 1

        hits_here = 0

        for enz in enzymes:
            cuts_r = cut_positions(enz, seq_ref_amp)
            cuts_a = cut_positions(enz, seq_alt_amp)

            if len(cuts_r) - 2 > max_cuts or len(cuts_a) - 2 > max_cuts:
                continue

            fr_r = frag_lengths(cuts_r)
            fr_a = frag_lengths(cuts_a)

            if (fr_r and (min(fr_r) < min_frag or max(fr_r) > max_frag)) or (
                fr_a and (min(fr_a) < min_frag or max(fr_a) > max_frag)
            ):
                continue

            if gain_loss_only and len(cuts_r) == len(cuts_a):
                continue

            actual_delta = diagnostic_delta(fr_r, fr_a)
            if actual_delta < delta:
                continue

            patt = "gain/loss" if (len(cuts_r) != len(cuts_a)) else "shift"
            site = enzyme_site_string(enz)

            row = [
                f"{chrom}:{pos} {ref}>{alt}",
                mapped_id,
                enz.__name__,
                site,
                patt,
                frags_to_str(fr_r),
                frags_to_str(fr_a),
                int(actual_delta),
                int(amplicon_start_genome),
                int(amplicon_end_genome),
                int(amplicon_len),
                int(snp_off_amp),
                primers_left,
                primers_right,
                (round(tm_left, 3) if tm_left is not None else None),
                (round(tm_right, 3) if tm_right is not None else None),
                int(amplicon_len),
                (int(p3_size) if p3_size is not None else None),
                (int(pl_start) if pl_start is not None else None),
                (int(pl_len) if pl_len is not None else None),
                (int(pr_start) if pr_start is not None else None),
                (int(pr_len) if pr_len is not None else None),
                enzyme_suppliers_human(enz),
                primer_status,
                '',
                int(delta),
            ]
            out_rows.append(row)
            hits_here += 1

        if hits_here == 0:
            final_status = 'no_enzyme_found' if primer_status in ('ok', 'window_only') else primer_status
            final_reason = 'Подходящий фермент не найден' if final_status == 'no_enzyme_found' else primer_reason
            row = [
                f"{chrom}:{pos} {ref}>{alt}",
                mapped_id,
                "",
                "",
                "no_enzyme_found",
                "",
                "",
                int(delta),
                int(amplicon_start_genome),
                int(amplicon_end_genome),
                int(amplicon_len),
                int(snp_off_amp),
                primers_left,
                primers_right,
                (round(tm_left, 3) if tm_left is not None else None),
                (round(tm_right, 3) if tm_right is not None else None),
                int(amplicon_len),
                (int(p3_size) if p3_size is not None else None),
                (int(pl_start) if pl_start is not None else None),
                (int(pl_len) if pl_len is not None else None),
                (int(pr_start) if pr_start is not None else None),
                (int(pr_len) if pr_len is not None else None),
                "",
                final_status,
                final_reason,
                None,
            ]
            out_rows.append(row)

    if report is not None:
        report.setdefault('summary', {}).update(candidates=sum(1 for r in out_rows if r[2]), rows=len(out_rows))
        report.setdefault('parameters', dict(params))
    try:
        fa.close()
    except Exception:
        pass
    return header, out_rows


