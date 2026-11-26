import sys
import bisect
import csv
import gzip
import json
import locale
import os
import re
import logging
import math
from itertools import combinations
from dataclasses import dataclass, field
from typing import Callable, List, Dict, Tuple, Optional
from datetime import datetime

from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio.Restriction import AllEnzymes
from allele_specific import design_assays as design_allele_specific
from assay_scoring import AGAROSE_MODEL_VERSION, GelModel, MODEL_VERSION, diploid_genotypes, score_genotypes, score_locus, rank_key
from gel_model import DEFAULT_AGAROSE_PCT, DEFAULT_DIGEST_PCT, DEFAULT_GEL_MARKER, GEL_MARKERS, digest_probability
import specificity
from tabix_index import index_path as tabix_path, vcf_chunks
from variant_input import VCF_COLUMNS, VCF_REJECTION_TEXT, is_vcf_path, parse_region, vcf_targets

APP_VERSION = "2.0.0"

# Analysis defaults of the window and the command line; the checks were chosen
# from the validation on published assays.
DEFAULT_PARAMS = dict(
    names2_file='', assembly_accession='', vcf_region='', vcf_samples=[],
    flank=300, min_frag=80, max_frag=800, delta=25, max_cuts=3,
    gain_loss_only=False, double_digest=True, allele_specific=True, suppliers=[], use_primer3=True, prod_min=250, prod_max=600,
    tm_min=58.0, tm_max=62.0, primer_pairs=10, top_results=10,
    gel_resolution_bp=10.0, gel_resolution_pct=3.0, gel_stress_factor=1.5,
    gel_model_version=AGAROSE_MODEL_VERSION, gel_agarose_pct=DEFAULT_AGAROSE_PCT,
    gel_marker=DEFAULT_GEL_MARKER, gel_digest_pct=DEFAULT_DIGEST_PCT,
    skip_norm=False, no_save_norm=True,
    enzyme_mode='off', lab_enzymes=[], commercial_only=True, group_isoschizomers=True, snp_end_distance=100,
    exclude_primer_sites=True, mask_repeats=True,
    background_vcf='', exclude_background_sites=True, dcaps='fallback', dcaps_mismatches=1,
    # dCAPS designs only: a longer mismatch primer moves the created cut away from the
    # product end, and a short product shows the ~20 bp difference on a gel.
    dcaps_tm_max=70.0, dcaps_prod_min=120, dcaps_prod_max=300, dcaps_delta=20,
    # Primer pairs that would amplify something else in the genome (specificity.py).
    specificity='exclude', specificity_mismatches=3, specificity_max_product=2000,
)
# Parameters a caller may leave out: run_rflp_gui_mode() takes them from DEFAULT_PARAMS.
OPTIONAL_PARAMS = (
    'double_digest', 'allele_specific', 'assembly_accession', 'vcf_region', 'vcf_samples', 'enzyme_mode', 'lab_enzymes', 'primer_pairs', 'top_results',
    'gel_resolution_bp', 'gel_resolution_pct',
    'gel_model_version', 'gel_agarose_pct', 'gel_marker', 'gel_digest_pct',
    'gel_stress_factor', 'commercial_only', 'group_isoschizomers', 'snp_end_distance',
    'exclude_primer_sites', 'mask_repeats', 'background_vcf', 'exclude_background_sites', 'dcaps',
    'dcaps_mismatches', 'dcaps_tm_max', 'dcaps_prod_min', 'dcaps_prod_max', 'dcaps_delta',
    'specificity', 'specificity_mismatches', 'specificity_max_product',
)


ENZYME_NAMES = {e.__name__.casefold(): e.__name__ for e in AllEnzymes}


def read_enzyme_list(path) -> list[str]:
    """Read a UTF-8 inventory, returning unique canonical enzyme names."""
    with open(path, encoding='utf-8-sig') as stream:
        names = re.split(r'[\s,;]+', stream.read().strip())
    names = [name for name in names if name]
    unknown = sorted({name for name in names if name.casefold() not in ENZYME_NAMES})
    if unknown:
        raise ValueError('Неизвестные ферменты: ' + ', '.join(unknown))
    if not names:
        raise ValueError('Список ферментов пуст')
    return sorted({ENZYME_NAMES[name.casefold()] for name in names})


class AnalysisCancelled(Exception):
    """Raised when a running analysis is cancelled by the caller."""


def validate_params(params: dict, *, check_files: bool = True) -> list[str]:
    errors = []
    mode = params.get('enzyme_mode', 'off')
    if mode not in ('off', 'only', 'prefer'):
        errors.append('enzyme_mode: допустимо off, prefer или only')
    if mode != 'off':
        inventory = params.get('lab_enzymes')
        if not isinstance(inventory, list) or not inventory:
            errors.append('lab_enzymes: список ферментов пуст')
        elif any(not isinstance(name, str) or ENZYME_NAMES.get(name.casefold()) != name for name in inventory):
            errors.append('lab_enzymes: список содержит неизвестные ферменты')
    if check_files:
        for key in ('fasta', 'input_file'):
            if not str(params.get(key, '')).strip():
                errors.append(f"Не задан параметр: {key}")
        for key in ('fasta', 'input_file', 'background_vcf'):
            value = str(params.get(key, '') or '').strip()
            if value and not os.path.isfile(value):
                errors.append(f"Файл не найден ({key}): {value}")
    region, samples = params.get('vcf_region', ''), params.get('vcf_samples', ())
    if region or samples:
        input_file = params.get('input_file', '')
        if (check_files or input_file) and not is_vcf_path(input_file):
            errors.append('input_file: регион и образцы применяются к входному VCF')
        try:
            parse_region(region)
        except ValueError as exc:
            errors.append(f'input_file: {exc}')
    def positive(*values):
        try:
            numbers = [float(value) for value in values]
        except (TypeError, ValueError):
            return None
        return numbers if all(math.isfinite(x) and x > 0 for x in numbers) else None

    def integer_in(key, default, low, high):
        value = params.get(key, default)
        if isinstance(value, bool) or not isinstance(value, int) or not low <= value <= high:
            errors.append(f"{key}: требуется целое число от {low} до {high}")

    for lo, hi, label in (('min_frag', 'max_frag', 'размер фрагмента'),
                          ('prod_min', 'prod_max', 'размер продукта'),
                          ('tm_min', 'tm_max', 'Tm')):
        numbers = positive(params.get(lo), params.get(hi))
        if numbers is None:
            errors.append(f"Некорректный параметр: {lo}/{hi}")
        elif numbers[0] > numbers[1]:
            errors.append(f"Минимальный {label} больше максимального")
    for key, default, low, high in (
        ('primer_pairs', 10, 1, 100), ('top_results', 10, 0, 10000),
        ('flank', 300, 1, 100000), ('max_cuts', 3, 0, 100),
        ('delta', 25, 1, 100000), ('snp_end_distance', 0, 0, 5000), ('dcaps_mismatches', 1, 1, 2),
    ):
        integer_in(key, default, low, high)
    # The dCAPS settings matter only when mismatch primers are designed; parameter sets
    # written before they existed get their defaults.
    if params.get('dcaps', 'off') in ('fallback', 'always') and params.get('use_primer3', True):
        dcaps = {key: params.get(key, DEFAULT_PARAMS[key])
                 for key in ('dcaps_prod_min', 'dcaps_prod_max', 'dcaps_tm_max', 'dcaps_delta')}
        numbers = positive(dcaps['dcaps_prod_min'], dcaps['dcaps_prod_max'])
        if numbers is None:
            errors.append("Некорректный параметр: dcaps_prod_min/dcaps_prod_max")
        elif numbers[0] > numbers[1]:
            errors.append("Минимальный размер продукта dCAPS больше максимального")
        ceiling = positive(dcaps['dcaps_tm_max'])
        if ceiling is None:
            errors.append("Некорректный параметр: dcaps_tm_max")
        elif positive(params.get('tm_min')) and float(params['tm_min']) > ceiling[0]:
            errors.append("Минимальный Tm праймеров dCAPS больше максимального")
        integer_in('dcaps_delta', dcaps['dcaps_delta'], 1, 1000)
    if params.get('dcaps', 'off') not in ('off', 'fallback', 'always'):
        errors.append("dcaps: допустимо off, fallback или always")
    if params.get('specificity', 'off') not in ('off', 'warn', 'exclude'):
        errors.append("specificity: допустимо off, warn или exclude")
    integer_in('specificity_mismatches', DEFAULT_PARAMS['specificity_mismatches'], 0, 5)
    integer_in('specificity_max_product', DEFAULT_PARAMS['specificity_max_product'], 100, 20000)
    try:
        gel_model_from_params(params)
        digest_probability(params.get('gel_digest_pct', DEFAULT_DIGEST_PCT))
    except (ValueError, TypeError) as exc:
        errors.append(str(exc))
    marker = params.get('gel_marker', DEFAULT_GEL_MARKER)
    if not isinstance(marker, str) or marker not in GEL_MARKERS:
        errors.append('gel_marker: неизвестный маркер длин')
    return errors


def gel_model_from_params(params):
    return GelModel(params.get('min_frag', 80), params.get('max_frag', 800),
                    params.get('gel_resolution_bp', 10.0),
                    params.get('gel_resolution_pct', 3.0),
                    params.get('gel_stress_factor', 1.5),
                    version=params.get('gel_model_version', MODEL_VERSION),
                    agarose_pct=params.get('gel_agarose_pct', DEFAULT_AGAROSE_PCT))


# ====== primer3 ======
def _use_c_numeric_locale() -> None:
    """Primer3 parses its numbers with the C library and needs a decimal point.

    Qt sets the locale of the system at startup (on Linux and macOS), and with a
    decimal comma Primer3 rejects "60.0"; its thermodynamic tables,
    read at import, would silently get wrong values. Qt itself formats numbers with QLocale.
    """
    if locale.localeconv()['decimal_point'] != '.':
        locale.setlocale(locale.LC_NUMERIC, 'C')


_PRIMER3_AVAILABLE = False
_PRIMER3_USE_NEW = False
_PRIMER3_IMPORT_ERROR = ""

try:
    _use_c_numeric_locale()
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


def read_input_table(path: str, *, vcf_region='', vcf_samples=()) -> Tuple[Dict[str, Optional[int]], List[Tuple[int, List[str]]], str]:
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
    if is_vcf_path(path):
        rows = [(line, fields) for line, fields, _, _ in vcf_targets(path, vcf_region, vcf_samples)]
        return VCF_COLUMNS, rows, '\t'
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


def clean_chrom(text: str) -> str:
    """A chromosome name of the input file without the leading "chr"."""
    return text.strip().removeprefix('chr').strip()


def resolve_contig(chrom: str, names: Dict[str, str], has_contig) -> str:
    """The FASTA contig of an input chromosome name.

    The names file wins; a genome whose contigs are called chr1, chr2, …
    is found without one. An unknown name is returned as it is.
    """
    mapped = names.get(chrom, chrom)
    if not has_contig(mapped) and has_contig('chr' + chrom):
        return 'chr' + chrom
    return mapped


def resolve_input_contig(chrom, names, has_contig, from_vcf=False):
    if from_vcf:
        if has_contig(chrom):
            return chrom
        if chrom in names and has_contig(names[chrom]):
            return names[chrom]
        chrom = clean_chrom(chrom)
    return resolve_contig(chrom, names, has_contig)


def parse_fields(parts: List[str], columns: Dict[str, Optional[int]], preserve_chrom=False) -> Tuple[str, int, str, str]:
    """(chrom, pos, ref, alt) of one data line; raises IndexError/ValueError."""
    def field(name, default_index):
        index = columns.get(name)
        return parts[index if index is not None else default_index]

    return (
        field("chrom", 0).strip() if preserve_chrom else clean_chrom(field("chrom", 0)),
        int(field("pos", 1)),
        field("ref", 2).strip().upper(),
        ','.join(allele.strip() for allele in field("alt", 3).upper().split(',')),
    )


def parse_input(path: str, report: Optional[dict] = None) -> List[Tuple[str, int, str, str]]:
    return [variant for _, variant in parse_input_lines(path, report)]


def parse_input_lines(path: str, report: Optional[dict] = None, *, vcf_region='', vcf_samples=(),
                      cancel_cb=None) -> List[Tuple[int, Tuple[str, int, str, str]]]:
    """The variants of the input file with the numbers of their lines."""
    rows: List[Tuple[int, Tuple[str, int, str, str]]] = []
    if is_vcf_path(path):
        def check_cancel():
            if cancel_cb and cancel_cb():
                raise AnalysisCancelled()

        for line, fields, reason, omitted in vcf_targets(path, vcf_region, vcf_samples, check_cancel):
            if report is not None and omitted:
                report.setdefault('vcf_skipped_alts', []).append(dict(line=line, alleles=list(omitted)))
            try:
                variant = parse_fields(fields, VCF_COLUMNS, preserve_chrom=True) if reason != 'parse_error' else None
            except (IndexError, ValueError):
                variant = None
            if variant and not reason:
                rows.append((line, variant))
            elif report is not None:
                reason = reason or 'parse_error'
                report.setdefault('rejected', []).append(dict(line=line, reason=reason,
                                                               details=VCF_REJECTION_TEXT[reason], raw='\t'.join(fields)))
        return rows
    try:
        columns, lines, delim = read_input_table(path)
        for line_no, parts in lines:
            if len(parts) < 4:
                if report is not None:
                    report.setdefault('rejected', []).append({'line': line_no, 'reason': 'parse_error', 'details': 'Ожидалось минимум 4 поля', 'raw': delim.join(parts)})
                continue
            try:
                rows.append((line_no, parse_fields(parts, columns)))
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


def is_simple_snp(ref: str, alt: str) -> bool:
    parts = [p for p in alt.split(",") if p != ""]
    return (len(ref) == 1) and (len(parts) == 1) and (parts[0] in DNA)


def is_indel(ref: str, alt: str) -> bool:
    """An explicit single insertion or deletion with a nonempty shared allele."""
    if not ref or not alt or len(ref) == len(alt) or not set(ref + alt) <= DNA:
        return False
    short, long = sorted((ref, alt), key=len)
    prefix = len(os.path.commonprefix((short, long)))
    suffix = len(os.path.commonprefix((short[::-1], long[::-1])))
    return prefix + suffix >= len(short)


def multiallelic_alts(ref: str, alt: str) -> tuple[str, ...]:
    alleles = tuple(alt.split(','))
    if (len(alleles) > 1 and len(set(alleles)) == len(alleles) and ref not in alleles
            and ref and set(ref) <= DNA and all(a and set(a) <= DNA for a in alleles)):
        return alleles
    return ()


UNSUPPORTED_VARIANT = 'Укажите поддерживаемый SNP, индел или набор разных ALT из A/C/G/T без повтора REF'


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
    vcf_region='',
    vcf_samples=(),
    accepted_lines=None,
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
    from_vcf = is_vcf_path(input_path)
    variants_raw = parse_input_lines(input_path, report=report, vcf_region=vcf_region,
                                     vcf_samples=vcf_samples, cancel_cb=cancel_cb)
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
        # Skipping normalization disables swapping, but still validates REF.
        valid = []
        check_cache = FastaBlockCache(fasta_path, block_bp=block_bp)
        names2_skip = load_names2(names2_path) if names2_path else {}
        try:
            for index, (line_no, (chrom, pos, ref, alt)) in enumerate(variants_raw, start=1):
                if cancel_cb and cancel_cb():
                    raise AnalysisCancelled()
                mapped = resolve_input_contig(chrom, names2_skip, check_cache.has_contig, from_vcf)
                reason = ''
                details = ''
                if not (is_simple_snp(ref, alt) or is_indel(ref, alt) or multiallelic_alts(ref, alt)):
                    reason, details = 'unsupported_variant', UNSUPPORTED_VARIANT
                elif not check_cache.has_contig(mapped):
                    reason, details = 'bad_chrom', f"CHROM '{mapped}' отсутствует в FASTA"
                elif pos < 1 or pos + len(ref) - 1 > check_cache.contig_len(mapped):
                    reason, details = 'pos_oob', 'Позиция вне FASTA'
                elif check_cache.slice(mapped, pos, pos + len(ref) - 1) != ref:
                    reason, details = 'ref_not_match', 'REF не совпадает с FASTA'
                if reason:
                    if report is not None:
                        report.setdefault('rejected', []).append({'line': line_no, 'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt, 'reason': reason, 'details': details})
                else:
                    valid.append((chrom, pos, ref, alt))
                    if accepted_lines is not None:
                        accepted_lines.append(line_no)
                if progress_cb:
                    progress_cb(index, max(1, total))
        finally:
            check_cache.fa.close()
        if report is not None:
            report.setdefault('summary', {}).update(accepted=len(valid), rejected=len(report.get('rejected', [])))
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
    rejected: List[Tuple[int, str, int, str, str, str, str]] = []   # the input line first

    def ref_matches(mapped_chrom: str, pos1: int, ref: str) -> bool:
        end = pos1 + len(ref) - 1
        if pos1 < 1 or end > cache.contig_len(mapped_chrom):
            return False
        return cache.slice(mapped_chrom, pos1, end) == ref

    for i, (line_no, (chrom, pos, ref, alt)) in enumerate(variants_raw, start=1):
        if cancel_cb and cancel_cb():
            try:
                cache.fa.close()
            except Exception:
                pass
            raise AnalysisCancelled()
        if progress_cb:
            progress_cb(i, max(1, total))
        orig_chrom = chrom
        mapped_chrom = resolve_input_contig(orig_chrom, names2, cache.has_contig, from_vcf)

        if progress_every and (i % progress_every == 0) and status_cb:
            status_cb(f"[NORMALIZE] {i}/{total} {orig_chrom}:{pos}")

        if not cache.has_contig(mapped_chrom):
            stats["rejected_bad_chrom"] += 1
            rejected.append((line_no, orig_chrom, pos, ref, alt, "bad_chrom", f"CHROM '{mapped_chrom}' not in FASTA"))
            continue

        clen = cache.contig_len(mapped_chrom)
        explicit = is_indel(ref, alt) or bool(multiallelic_alts(ref, alt))
        if pos < 1 or pos > clen or (explicit and pos + len(ref) - 1 > clen):
            stats["rejected_pos_oob"] += 1
            rejected.append((line_no, orig_chrom, pos, ref, alt, "pos_oob", f"POS {pos} out of [1,{clen}]"))
            continue

        if (not ref) or (not set(ref) <= DNA):
            stats["rejected_bad_ref_chars"] += 1
            rejected.append((line_no, orig_chrom, pos, ref, alt, "bad_ref_chars", "REF must be A/C/G/T only"))
            continue

        if not explicit and not is_simple_snp(ref, alt):
            stats["rejected_bad_alt_chars"] += 1
            rejected.append((line_no, orig_chrom, pos, ref, alt, "unsupported_variant", UNSUPPORTED_VARIANT))
            continue

        if ref_matches(mapped_chrom, pos, ref):
            cleaned.append((orig_chrom, pos, ref, alt))
            if accepted_lines is not None:
                accepted_lines.append(line_no)
            stats["kept"] += 1
            continue

        if explicit or from_vcf:
            stats['rejected_ref_not_match'] += 1
            rejected.append((line_no, orig_chrom, pos, ref, alt, 'ref_not_match', 'REF не совпадает с FASTA'))
            continue

        refbase = cache.base(mapped_chrom, pos)
        alt_base = alt.replace(",", "")   # "G," is read as G, as before
        if alt_base == refbase:
            cleaned.append((orig_chrom, pos, alt_base, ref))
            if accepted_lines is not None:
                accepted_lines.append(line_no)
            stats["kept_swapped_snp"] += 1
            continue

        if rc(ref) == refbase:
            stats["rejected_minus_strand_like"] += 1
            rejected.append((line_no, orig_chrom, pos, ref, alt, "minus_strand_like", "rc(REF) equals reference base"))
            continue

        stats["rejected_ref_not_match"] += 1
        rejected.append((line_no, orig_chrom, pos, ref, alt, "ref_not_match", "REF does not match reference and no allowed fix"))

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
                w.writerow(list(row[1:]))

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
            dict(zip(('line', 'chrom', 'pos', 'ref', 'alt', 'reason', 'details'), r))
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
                             cancel_cb=None, *, fasta_contigs=None) -> Tuple[Dict[str, list], int]:
    """Known variants of a VCF that lie inside ``regions``, using a current tabix index when present.

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
    records = 0
    contigs = {}  # VCF contig name -> FASTA contig
    has_contig = (fasta_contigs if fasta_contigs is not None else regions).__contains__

    def map_contig(name):
        if name not in contigs:
            contigs[name] = resolve_input_contig(name, names, has_contig, from_vcf=True)
        return contigs[name]

    def windows(name):
        left, right = merged.get(map_contig(name), ((), ()))
        return zip(left, right)

    def lines():
        indexed = tabix_path(path)
        if indexed:
            yield from vcf_chunks(path, indexed, windows)
        else:
            with open(path, 'rb') as probe:
                opener = gzip.open if probe.read(2) == b'\x1f\x8b' else open
            with opener(path, 'rt', encoding='utf-8', errors='replace') as stream:
                yield from stream

    for line in lines():
        if line.startswith('#'):
            continue
        records += 1
        if cancel_cb and (records == 1 or records % 200_000 == 0) and cancel_cb():
            raise AnalysisCancelled()
        fields = line.split('\t', 5)
        if len(fields) < 5:
            continue
        contig = map_contig(fields[0])
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


def _hit_enzymes(hit):
    if '_enz' not in hit:
        return ()
    return (hit['_enz'], hit['_enz2']) if '_enz2' in hit else (hit['_enz'],)


def _hit_sites(hit):
    return tuple(sorted(site_key(enzyme) for enzyme in _hit_enzymes(hit)))


def _inventory_priority(hit, preferred):
    if preferred:
        return (hit['genotype_quality'] != 'robust', not all(e.__name__ in preferred for e in _hit_enzymes(hit)))
    return ()


def group_by_site(hits: list, preferred=frozenset()) -> list:
    """Keep one hit per recognition site; the others are listed as isoschizomers.

    ``hits`` carry the enzyme class in ``_enz``. A site is represented by a
    listed enzyme within its quality class in prefer mode, then by a
    commercially available enzyme when there is one (a neoschizomer cutting
    a base further may score marginally higher), then by the best score, then
    by the number of suppliers. The result is in rank order.
    """
    grouped = []
    groups = {}
    for hit in hits:
        sites = _hit_sites(hit)
        if sites:
            groups.setdefault((hit.get('design'), sites), []).append(hit)
        else:
            grouped.append(hit)
    for members in groups.values():
        keep = min(members, key=lambda h: (*_inventory_priority(h, preferred), _flagged(h),
                                           not all(e.suppl for e in _hit_enzymes(h)), rank_key(h)[:4],
                                           -sum(len(e.suppl) for e in _hit_enzymes(h)), rank_key(h)[4:]))
        others = sorted({h['enzyme'] for h in members} - {keep['enzyme']})
        grouped.append(dict(keep, isoschizomers='; '.join(others)))
    return sorted(grouped, key=lambda h: _order(h, preferred))


def _flagged(hit: dict) -> bool:
    """The specificity check has marked the primer pair of the assay (it stays only with specificity=warn)."""
    return bool(hit.get('_unspecific'))


def _order(hit: dict, preferred=frozenset()) -> tuple:
    """Inventory preferences precede the existing specificity and band ranking."""
    return (*_inventory_priority(hit, preferred), _flagged(hit), *rank_key(hit))


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
    target_length: int = 1,
    force_left_end: Optional[int] = None,
    force_right_end: Optional[int] = None,
) -> Optional[dict]:
    """Pick PCR primers flanking the reference interval and its target flanks.

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
    _use_c_numeric_locale()
    seq_args = {"SEQUENCE_TEMPLATE": amplicon,
                "SEQUENCE_TARGET": [snp_off - target_flank, 2 * target_flank + target_length]}
    if force_left_end is not None:
        seq_args['SEQUENCE_FORCE_LEFT_END'] = force_left_end
    if force_right_end is not None:
        seq_args['SEQUENCE_FORCE_RIGHT_END'] = force_right_end
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
# Added in 1.2 at the end of the row.
V12_HEADER = ['isoschizomers', 'primer_site', 'repeat_pct', 'background_variants',
              'background_site_variants', 'design', 'primer_mismatches']
# Added in 1.6: where else the primers anneal and what the pair amplifies there.
V16_HEADER = ['primer_sites_left', 'primer_sites_right', 'offtarget_products', 'offtarget_sizes']
DOUBLE_HEADER = ['enzyme_1', 'enzyme_2', 'site_1', 'site_2', 'reaction_compatibility']
LOCUS_HEADER = ['allele_products', 'genotype_bands', 'worst_genotype_bands']
GEL_HEADER = ['gel_agarose_pct', 'gel_marker', 'gel_digest_pct', 'gel_resolution_bp',
              'gel_resolution_pct', 'gel_stress_factor']
ASPCR_HEADER = ['aspcr_reactions', 'aspcr_control', 'aspcr_genotypes', 'aspcr_specific_side']
RESULT_HEADER = (BASE_HEADER + ASSAY_HEADER + V12_HEADER + V16_HEADER + DOUBLE_HEADER + ['product_size_alt']
                 + LOCUS_HEADER + GEL_HEADER + ['gel_min_visible', 'gel_max_visible'] + ASPCR_HEADER)
# SNPs whose assays are held in memory for one pass over the genome (about 0.25 MB each).
SPECIFICITY_BATCH = 2000


def _primer_candidates(pres, sequence, snp_off, target_length=1):
    """Validate Primer3 coordinates and retain distinct primer pairs."""
    seen = set()
    for i in range(int(pres.get('PRIMER_PAIR_NUM_RETURNED', 0))):
        left = pres[f'PRIMER_LEFT_{i}_SEQUENCE'].upper()
        right = pres[f'PRIMER_RIGHT_{i}_SEQUENCE'].upper()
        ls, ll = map(int, pres[f'PRIMER_LEFT_{i}'])
        rs, rl = map(int, pres[f'PRIMER_RIGHT_{i}'])
        size = int(pres[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'])
        if not (0 <= ls < ls + ll <= snp_off < snp_off + target_length <= rs - rl + 1 <= rs < len(sequence)):
            target = 'SNP' if target_length == 1 else 'вариант'
            raise ValueError(f'Primer3 вернул праймеры вне окна или перекрывающие {target}')
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


@dataclass
class _Run:
    """What one analysis shares between its SNPs."""
    params: dict
    model: GelModel
    report: dict
    counters: dict
    enzymes: list
    known_cut_enzymes: list          # with a known cut position: the fragment lengths are defined
    use_primer3: bool
    mode: str                    # 'joint' (primer pairs and the enzyme together) or 'window'
    dcaps: str                   # 'off' without Primer3
    specificity: str = 'off'     # 'off' without Primer3
    status_cb: Optional[Callable] = None
    cancel_cb: Optional[Callable] = None
    background: Optional[dict] = None   # known variants by contig; None without a VCF

    def check_cancel(self) -> None:
        if self.cancel_cb and self.cancel_cb():
            self.report['cancelled'] = True
            raise AnalysisCancelled()

    def log(self, text: str) -> None:
        if self.status_cb:
            self.status_cb(text)


@dataclass
class _Window:
    """The genome around one SNP and what was counted while it was searched."""
    ref: str
    alt: str
    pos: int
    start: int                   # 1-based genome coordinate of the first base of the window
    end: int
    raw: str                     # as in the FASTA: lower case marks soft-masked repeats
    seq: str                     # upper case, REF allele
    off: int                     # 0-based offset of the SNP in the window
    base: dict                   # the result row before an assay is found
    altseq: str = ''
    known: list = field(default_factory=list)    # (offset, REF, ALTs) of the known variants of the window
    errors: int = 0
    ambiguous: int = 0
    genotype_rejected: int = 0
    indel: bool = False
    alternates: tuple[str, ...] = ()

    @property
    def length_change(self) -> int:
        return len(self.alt) - len(self.ref) if self.indel else 0

    @property
    def length_changes(self) -> tuple[int, ...]:
        return (0, *(len(alt) - len(self.ref) for alt in self.alternates)) if self.alternates else (0, self.length_change)

    def amplicons(self, start: int, end: int) -> tuple[str, str]:
        return self.seq[start:end], self.altseq[start:end + self.length_change]

    @property
    def excluded(self) -> list:
        """(offset, length) regions that primers must not overlap."""
        return [(q, len(r)) for q, r, _ in self.known]


def _log_settings(run: _Run, snp_count: int) -> None:
    if not run.status_cb:
        return
    params, use_primer3 = run.params, run.use_primer3
    run.log(f'[INFO] Поиск: {run.mode}; SNP: {snp_count}; ферментов: {len(run.enzymes)}')
    if run.specificity == 'off':
        run.log('[INFO] Оценка полос: модель присутствия, полная рестрикция; '
                'яркость и специфичность ПЦР по всему геному не оцениваются.')
    else:
        run.log('[INFO] Оценка полос: модель присутствия, полная рестрикция; яркость не оценивается.')
    options = [label for key, label in (
        ('commercial_only', 'только продаваемые ферменты'),
        ('group_isoschizomers', 'изошизомеры объединены (одна строка на сайт)'),
        ('exclude_primer_sites', 'без сайта фермента в праймерах'),
        ('mask_repeats', 'повторы (строчные буквы FASTA) не на 3′-концах праймеров'),
        ('exclude_background_sites', 'без фоновых вариантов, меняющих разрезы'),
    ) if params[key] and (key != 'commercial_only' or params.get('enzyme_mode', 'off') == 'off')
        and (key != 'exclude_background_sites' or params['background_vcf'])
        and (use_primer3 or key not in ('exclude_primer_sites', 'mask_repeats'))]   # these two need primers
    if params.get('enzyme_mode') in ('prefer', 'only'):
        options.append(f"enzyme_mode={params['enzyme_mode']} ({', '.join(params['lab_enzymes'])})")
    if use_primer3 and params['snp_end_distance']:
        options.append(f"SNP не ближе {params['snp_end_distance']} п.н. к концам ампликона")
    if use_primer3 and params['dcaps'] != 'off':
        when = 'если нет robust-схемы' if params['dcaps'] == 'fallback' else 'для каждого SNP'
        options.append(f"dCAPS {when}, несовпадений до {params['dcaps_mismatches']}, "
                       f"Tm до {params['dcaps_tm_max']:g} °C, "
                       f"продукт {params['dcaps_prod_min']}–{params['dcaps_prod_max']} п.н., "
                       f"Δ {params['dcaps_delta']} п.н.")
    if run.specificity != 'off':
        action = 'исключаются' if run.specificity == 'exclude' else 'отмечаются'
        options.append(f"пары праймеров с посторонним продуктом ПЦР {action} "
                       f"(до {params['specificity_mismatches']} несовпадений вне 3′-конца, "
                       f"продукт до {params['specificity_max_product']} п.н.)")
    if options:
        run.log('[INFO] Отбор: ' + '; '.join(options))


def _snp_window(run: _Run, fa, chrom: str, mapped_id: str, pos: int, ref: str, alt: str) -> _Window:
    params = run.params
    start = max(1, pos - params['flank'])
    end = min(len(fa[mapped_id]), pos + len(ref) - 1 + params['flank'])
    # Soft-masked FASTA: lower case marks repeats; the analysis uses upper case.
    raw = str(fa[mapped_id][start - 1:end])
    seq = raw.upper()
    off = pos - start
    base = dict(variant=f'{chrom}:{pos} {ref}>{alt}', mapped_id=mapped_id,
                enzyme='', site='', pattern='no_enzyme_found',
                amplicon_start=start, amplicon_end=end, amplicon_len=len(seq),
                snp_offset_in_amplicon=off, delta_threshold_bp=params['delta'],
                search_mode=run.mode, scoring_model=run.model.version, pairs_evaluated=0,
                gel_min_visible=run.model.min_visible, gel_max_visible=run.model.max_visible,
                **{key: params.get(key, DEFAULT_PARAMS[key]) for key in GEL_HEADER})
    return _Window(ref=ref, alt=alt, pos=pos, start=start, end=end, raw=raw, seq=seq, off=off,
                   base=base, indel=is_indel(ref, alt), alternates=multiallelic_alts(ref, alt))


def _known_variants(run: _Run, window: _Window, mapped_id: str) -> list:
    """Known variants of the window as (offset, REF, ALTs), except a record of the SNP itself."""
    known = []
    if run.background is not None:
        items = run.background.get(mapped_id, [])
        first, last = bisect.bisect_left(items, (window.start,)), bisect.bisect_left(items, (window.end + 1,))
        for p, r, alts in items[first:last]:
            q = p - window.start
            if p + len(r) - 1 > window.end or (p < window.pos + len(window.ref) and window.pos < p + len(r)):
                continue
            if window.seq[q:q + len(r)] != r:
                run.counters['background_ref_mismatch'] += 1
                continue
            known.append((q, r, alts))
    return known


def _natural_pairs(run: _Run, window: _Window) -> Tuple[Optional[dict], list]:
    """Primer pairs around the SNP (the whole window without Primer3), or the failure row."""
    params, base, seq, off = run.params, window.base, window.seq, window.off
    mask = params['mask_repeats']
    end_distance = params['snp_end_distance'] if run.use_primer3 else 0
    target_length = len(window.ref)
    target_end = off + target_length
    failure, pairs = None, []
    if end_distance and min(off, len(seq) - target_end) < end_distance:
        label = 'Вариант' if window.indel or window.alternates else 'SNP'
        failure = dict(base, status='no_primers', reason=f'{label} ближе {end_distance} п.н. к краю окна; увеличьте фланг')
    elif run.use_primer3:
        run.check_cancel()
        # Primer3 designs on REF. Restrict its range so every allele fits.
        changes = window.length_changes
        prod_min = params['prod_min'] - min(changes)
        prod_max = params['prod_max'] - max(changes)
        if prod_min > prod_max:
            return dict(base, status='no_primers', reason='Длины продуктов REF и ALT не укладываются в заданный диапазон'), []
        pres = primer3_pick(window.raw if mask else seq, off, prod_min, prod_max,
                            params['tm_min'], params['tm_max'],
                            num_return=params['primer_pairs'],
                            target_flank=max(0, end_distance - PRIMER_MIN_SIZE),
                            lowercase_masking=mask, excluded_regions=window.excluded,
                            target_length=target_length)
        if not pres:
            error = primer3_pick.last_error
            failure = dict(base, status='primer3_error' if error else 'no_primers',
                           reason=error or 'Primer3 не вернул подходящих пар')
        else:
            try:
                pairs = list(_primer_candidates(pres, seq, off, target_length))
            except (KeyError, TypeError, ValueError) as exc:
                failure = dict(base, status='primer3_error', reason=str(exc))
            else:
                # Guaranteed by the Primer3 target; kept as an explicit check.
                near = [p for p in pairs if min(off - p['start0'], p['end0'] - target_end) < end_distance]
                run.counters['pairs_near_end'] += len(near)
                pairs = [p for p in pairs if p not in near]
    else:
        pairs = [dict(start0=0, end0=len(seq))]
    return failure, pairs


def _amplicon_record(window: _Window, pair: dict, a: int, b: int, **extra) -> dict:
    """The result row of the amplicon window.seq[a:b] before an enzyme is chosen."""
    record = dict(window.base, **{k: v for k, v in pair.items() if k not in ('start0', 'end0')})
    record.update(amplicon_start=window.start + a, amplicon_end=window.start + b - 1,
                  amplicon_len=b - a, product_size=b - a,
                  snp_offset_in_amplicon=window.off - a,
                  repeat_pct=round(100 * sum(c.islower() for c in window.raw[a:b]) / (b - a), 1),
                  **extra)
    if window.indel:
        record['product_size_alt'] = b - a + window.length_change
    return record


def _enzyme_cuts(run, window, enzyme, seq):
    run.counters['digest_searches'] += 1
    try:
        return cut_positions(enzyme, seq)
    except Exception as exc:
        window.errors += 1
        if window.errors <= 3:
            run.log(f'[WARN] {window.base["variant"]} / {enzyme.__name__}: {exc}')
        return None


def _evaluate(run, window, enz, refamp, altamp, record, a, b):
    cr = _enzyme_cuts(run, window, enz, refamp)
    if cr is None:
        return None
    ca = _enzyme_cuts(run, window, enz, altamp)
    if ca is None or (window.indel and len(cr) == len(ca) == 2):
        return None
    return _score_digest(run, window, (enz,), cr, ca, refamp, altamp, record, a, b)


def _score_digest(run, window, enzymes, cr, ca, refamp, altamp, record, a, b):
    """Apply the same visibility, genotype and background rules to the final digest."""
    params, model = run.params, run.model
    if max(len(cr), len(ca)) - 2 > params['max_cuts']:
        return None
    fr, fal = frag_lengths(cr), frag_lengths(ca)
    if enzymes and params['gain_loss_only'] and len(cr) == len(ca):
        return None
    # min_frag..max_frag is the visibility window of the gel model:
    # fragments outside it are hidden, the assay itself stays.
    visible_r = [x for x in fr if model.min_visible <= x <= model.max_visible]
    visible_a = [x for x in fal if model.min_visible <= x <= model.max_visible]
    actual_delta = diagnostic_delta(visible_r, visible_a)
    threshold = params['dcaps_delta'] if record.get('design') == 'dCAPS' else params['delta']
    if actual_delta < threshold:
        return None
    scores = score_genotypes(fr, fal, model)
    if scores['genotype_quality'] == 'ambiguous':
        window.genotype_rejected += 1
        return None
    scores.update(frags_ref=frags_to_str(fr), frags_alt=frags_to_str(fal))
    return _digest_record(run, window, enzymes, (refamp, altamp), record, a, b, scores, actual_delta,
                          len(cr) != len(ca), 'REF' if len(cr) > len(ca) else 'ALT')


def _digest_record(run, window, enzymes, amplicons, record, a, b, scores, actual_delta, gain_loss, dcaps_allele=''):
    reason = {
        'robust': 'Три генотипа различимы во всём заданном диапазоне разрешения',
        'nominal_only': 'Три генотипа различимы только при исходном разрешении',
    }[scores['genotype_quality']]
    if window.alternates:
        reason = {
            'robust': 'Все генотипы локуса различимы во всём заданном диапазоне разрешения',
            'nominal_only': 'Все генотипы локуса различимы только при исходном разрешении',
        }[scores['genotype_quality']]
    if not run.use_primer3:
        reason = 'Оценка окна без подобранных праймеров. ' + reason
    dcaps_design = record.get('design') == 'dCAPS'
    if dcaps_design:
        reason = f'dCAPS: несовпадение в праймере создаёт сайт в аллеле {dcaps_allele}. ' + reason
    primer_sites = [primer_site_label(enz, amplicons[0], record.get('primer_left_len'), record.get('primer_right_len'))
                    for enz in enzymes]
    if not enzymes:
        primer_site = ''
    elif len(enzymes) == 1:
        primer_site = primer_sites[0]
    else:
        primer_site = '; '.join(f'{enz.__name__}: {label}' for enz, label in zip(enzymes, primer_sites) if label)
    # In a mismatch primer such a site is a design error: every allele is cut.
    if primer_site and (run.params['exclude_primer_sites'] or dcaps_design):
        run.counters['excluded_primer_site'] += 1
        return None
    inside = [(q - a, r, alts) for q, r, alts in window.known if a <= q and q + len(r) <= b]
    changing = set()
    shifts = window.length_changes
    for q, r, alts in inside:
        downstream = a + q > window.off
        for x in alts:
            if enzymes:
                affects_bands = any(variant_changes_sites(enz, amp, q + (shift if downstream else 0),
                                                         r, x, run.model.resolution_bp)
                                    for enz in enzymes for amp, shift in zip(amplicons, shifts))
            else:
                affects_bands = len(r) != len(x)
            if affects_bands:
                changing.add(f'{window.start + a + q} {r}>{x}')
    if changing and run.params['exclude_background_sites']:
        run.counters['excluded_background_site'] += 1
        return None
    names = [enz.__name__ for enz in enzymes]
    sites = [enzyme_site_string(enz) for enz in enzymes]
    suppliers = enzyme_suppliers_human(enzymes[0]) if enzymes else ''
    extra = {'_enz': enzymes[0]} if enzymes else {}
    if len(enzymes) == 2:
        extra.update(_enz2=enzymes[1], enzyme_1=names[0], enzyme_2=names[1], site_1=sites[0], site_2=sites[1],
                     reaction_compatibility='unconfirmed')
        reason += '. Совместимость условий не подтверждена'
        suppliers = '; '.join(f'{enz.__name__}: {enzyme_suppliers_human(enz)}' for enz in enzymes)
    return dict(record, **scores, **extra, enzyme=' + '.join(names),
                site=' + '.join(sites),
                pattern=('gain/loss' if gain_loss else 'shift') if enzymes else 'amplicon_length',
                diag_delta_bp=actual_delta, suppliers=suppliers,
                primer_site=primer_site,
                background_variants=len(inside) if run.background is not None else None,
                background_site_variants='; '.join(sorted(changing, key=lambda text: int(text.split()[0]))),
                status='ok' if run.use_primer3 else 'window_only', reason=reason)


def _score_multiallelic_digest(run, window, enzymes, cuts, amplicons, record, start, end):
    counts = {len(positions) for positions in cuts}
    if (max(counts) - 2 > run.params['max_cuts'] or
            (enzymes and run.params['gain_loss_only'] and len(counts) == 1)):
        return None
    fragments = [frag_lengths(positions) for positions in cuts]
    scores = score_locus(fragments, run.model, run.check_cancel)
    if scores['genotype_quality'] == 'ambiguous':
        window.genotype_rejected += 1
        return None
    visible = [{n for n in sizes if run.model.min_visible <= n <= run.model.max_visible} for sizes in fragments]
    lanes = [visible[i] | visible[j] for _, i, j in diploid_genotypes(len(fragments))]
    actual_delta = math.inf
    for first, second in combinations(lanes, 2):
        run.check_cancel()
        actual_delta = min(actual_delta, diagnostic_delta(first, second))
        if actual_delta < run.params['delta']:
            return None
    alleles = (window.ref, *window.alternates)
    products = [dict(allele='REF' if i == 0 else f'ALT{i}', sequence=allele,
                     length=len(amplicons[i]), fragments=sizes)
                for i, (allele, sizes) in enumerate(zip(alleles, fragments))]
    scores['allele_products'] = json.dumps(products)
    for key in ('genotype_bands', 'worst_genotype_bands'):
        scores[key] = json.dumps(scores[key], ensure_ascii=False)
    return _digest_record(run, window, enzymes, amplicons, record, start, end, scores, actual_delta, len(counts) > 1)


def _multiallelic_hits(run, window, pairs, design):
    enzymes = run.known_cut_enzymes
    counters = run.counters
    if design == 'double_digest':
        counters['double_snps'] += 1
    hits = []
    for pair in pairs:
        run.check_cancel()
        a, b = pair['start0'], pair['end0']
        prefix, suffix = window.seq[a:window.off], window.seq[window.off + len(window.ref):b]
        amplicons = (window.seq[a:b], *(prefix + alt + suffix for alt in window.alternates))
        if any(set(amp) - DNA for amp in amplicons):
            if design == 'natural':
                window.ambiguous += 1
            continue
        record = _amplicon_record(window, pair, a, b, design=design)

        def append_hit(members, cuts):
            hit = _score_multiallelic_digest(run, window, members, cuts, amplicons, record, a, b)
            if hit is not None:
                hits.append(hit)

        if design == 'amplicon_length':
            append_hit((), [(0, len(amp)) for amp in amplicons])
            continue
        if design == 'natural':
            counters['enzyme_evaluations'] += len(enzymes)
        else:
            counters['double_amplicons'] += 1
            counters['double_pairs_possible'] += len(enzymes) * (len(enzymes) - 1) // 2
        maps = {}
        for enzyme in enzymes:
            run.check_cancel()
            cuts = [_enzyme_cuts(run, window, enzyme, amp) for amp in amplicons]
            if any(c is None or len(c) - 2 > run.params['max_cuts'] for c in cuts):
                continue
            if all(len(c) == 2 for c in cuts):
                continue
            if design == 'natural':
                append_hit((enzyme,), cuts)
            else:
                maps[enzyme] = tuple(frozenset(c) for c in cuts)
        diagnostic = {enzyme for enzyme, cuts in maps.items() if len(set(cuts)) > 1}
        for (first, left), (second, right) in combinations(maps.items(), 2):
            run.check_cancel()
            if first not in diagnostic and second not in diagnostic:
                continue
            combined = tuple(l | r for l, r in zip(left, right))
            if combined in (left, right) or any(len(c) - 2 > run.params['max_cuts'] for c in combined):
                continue
            counters['double_pairs_scored'] += 1
            append_hit((first, second), [sorted(c) for c in combined])
    return hits


def _length_hits(run: _Run, window: _Window, pairs: list) -> list:
    if window.alternates:
        return _multiallelic_hits(run, window, pairs, 'amplicon_length')
    hits = []
    for pair in pairs:
        run.check_cancel()
        a, b = pair['start0'], pair['end0']
        ref, alt = window.amplicons(a, b)
        if set(ref + alt) - DNA:
            continue
        record = _amplicon_record(window, pair, a, b, design='amplicon_length')
        hit = _score_digest(run, window, (), (0, len(ref)), (0, len(alt)), ref, alt, record, a, b)
        if hit is not None:
            hits.append(hit)
    return hits


def _natural_hits(run: _Run, window: _Window, pairs: list) -> list:
    """Assays with a natural site: every primer pair with every enzyme whose site covers the SNP."""
    if window.alternates:
        return _multiallelic_hits(run, window, pairs, 'natural')
    # The remaining enzymes give equal REF/ALT fragments in every amplicon,
    # i.e. Δ = 0, which is below any threshold (delta >= 1).
    site_enzymes = (run.known_cut_enzymes if window.indel else
                    snp_site_enzymes(run.enzymes, window.seq, window.altseq, window.off))
    hits = []
    for pair in pairs:
        run.check_cancel()
        a, b = pair['start0'], pair['end0']
        refamp, altamp = window.amplicons(a, b)
        if set(refamp + altamp) - set('ACGT'):
            window.ambiguous += 1
            continue
        record = _amplicon_record(window, pair, a, b, design='natural')
        run.counters['enzyme_evaluations'] += len(run.enzymes)
        for enz in site_enzymes:
            run.check_cancel()
            hit = _evaluate(run, window, enz, refamp, altamp, record, a, b)
            if hit is not None:
                hits.append(hit)
    return hits


def _double_cut_maps(run, window, ref, alt, potential):
    """Keep usable cuts, resolving allele-dependent maps before common partners."""
    maps, diagnostic = {}, {}
    limit = run.params['max_cuts']
    for enzyme in potential:
        run.check_cancel()
        cr = _enzyme_cuts(run, window, enzyme, ref)
        if cr is None or len(cr) - 2 > limit:
            continue
        ca = _enzyme_cuts(run, window, enzyme, alt)
        if ca is None or len(ca) - 2 > limit:
            continue
        if len(cr) > 2 or len(ca) > 2:
            maps[enzyme] = frozenset(cr), frozenset(ca)
            if cr != ca:
                diagnostic[enzyme] = maps[enzyme]
    if diagnostic:
        for enzyme in run.known_cut_enzymes:
            if enzyme in potential:
                continue
            run.check_cancel()
            cuts = _enzyme_cuts(run, window, enzyme, ref)
            if cuts is not None and 0 < len(cuts) - 2 <= limit:
                common = frozenset(cuts)
                maps[enzyme] = common, common
    return maps, diagnostic


def _double_hits(run, window, pairs):
    """Try SNP-dependent cut maps with every productive partner on the same amplicon."""
    if window.alternates:
        return _multiallelic_hits(run, window, pairs, 'double_digest')
    counters = run.counters
    counters['double_snps'] += 1
    enzymes = run.known_cut_enzymes
    # Preserve catalog order when logging errors, with constant-time membership.
    potential = dict.fromkeys(enzymes if window.indel else
                              snp_site_enzymes(enzymes, window.seq, window.altseq, window.off))
    hits = []
    for pair in pairs:
        run.check_cancel()
        a, b = pair['start0'], pair['end0']
        ref, alt = window.amplicons(a, b)
        if set(ref + alt) - set('ACGT'):
            continue
        counters['double_amplicons'] += 1
        counters['double_pairs_possible'] += len(enzymes) * (len(enzymes) - 1) // 2
        record = _amplicon_record(window, pair, a, b, design='double_digest')
        maps, diagnostic = _double_cut_maps(run, window, ref, alt, potential)
        for first, (r1, a1) in diagnostic.items():
            for second, (r2, a2) in maps.items():
                run.check_cancel()
                if second in diagnostic and first.__name__ >= second.__name__:
                    continue
                if (r1 <= r2 and a1 <= a2) or (r2 <= r1 and a2 <= a1):
                    continue
                cr, ca = r1 | r2, a1 | a2
                if cr == ca or max(len(cr), len(ca)) - 2 > run.params['max_cuts']:
                    continue
                counters['double_pairs_scored'] += 1
                members = tuple(sorted((first, second), key=str))
                hit = _score_digest(run, window, members, sorted(cr), sorted(ca), ref, alt, record, a, b)
                if hit is not None:
                    hits.append(hit)
    return hits


def _dcaps_hits(run: _Run, window: _Window) -> list:
    """dCAPS assays: a mismatch primer next to the SNP creates the site in one allele."""
    params, counters, seq, off = run.params, run.counters, window.seq, window.off
    mask, raw, excluded = params['mask_repeats'], window.raw, window.excluded
    hits = []
    counters['dcaps_snps'] += 1
    designs = dcaps_mismatch_sets(run.known_cut_enzymes, seq, off, window.ref, window.alt, params['dcaps_mismatches'])
    ordered = sorted(designs.items(), key=lambda item: (len(item[0][1]), -len(item[1]), item[0]))
    for (side, changes), design_enzymes in ordered[:DCAPS_MAX_DESIGNS]:
        run.check_cancel()
        template = list(raw if mask else seq)
        for j, new in changes:
            template[j] = new
        template = ''.join(template)
        terminal = off - 1 if side == 'left' else off + 1
        primer = None if mask and raw[terminal].islower() else dcaps_primer(
            template, off, side, params['tm_min'], params['dcaps_tm_max'], excluded)
        if primer is None:
            counters['dcaps_no_primer'] += 1
            continue
        counters['dcaps_designs'] += 1
        pres = primer3_pick(template, off, params['dcaps_prod_min'], params['dcaps_prod_max'],
                            params['tm_min'], params['dcaps_tm_max'],
                            num_return=params['primer_pairs'],
                            lowercase_masking=mask, excluded_regions=excluded,
                            fixed_left=primer[1] if side == 'left' else None,
                            fixed_right=primer[1] if side == 'right' else None)
        mutated = template.upper()
        try:
            dpairs = list(_primer_candidates(pres, mutated, off)) if pres else []
        except (KeyError, TypeError, ValueError):
            dpairs = []
        counters['dcaps_pairs'] += len(dpairs)
        mutated_alt = mutated[:off] + window.alt + mutated[off + 1:]
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
            record = _amplicon_record(window, pair, a, b, design='dCAPS', primers_left=left,
                                      primers_right=right, delta_threshold_bp=params['dcaps_delta'],
                                      primer_mismatches=describe_mismatches(seq, off, side, changes))
            for enz in design_enzymes:
                hit = _evaluate(run, window, enz, refamp, altamp, record, a, b)
                if hit is not None:
                    hits.append(hit)
                    counters['dcaps_candidates'] += 1
    return hits


def _snp_records(run: _Run, window: _Window, hits: list, failure: Optional[dict], unspecific: int = 0) -> list:
    """Rank the assays of one SNP and return its result rows (one failure row without assays).

    ``unspecific`` is the number of assays the specificity check has removed."""
    params, counters = run.params, run.counters
    counters['enzyme_errors'] += window.errors
    counters['ambiguous_templates'] += window.ambiguous
    counters['genotype_rejected'] += window.genotype_rejected
    preferred = set(params['lab_enzymes']) if params.get('enzyme_mode') == 'prefer' else set()
    hits.sort(key=lambda h: _order(h, preferred))
    counters['candidates'] += len(hits)
    counters['candidate_sites'] += len({_hit_sites(h) for h in hits})
    if params['group_isoschizomers']:
        hits = group_by_site(hits, preferred)
    counters['variants_with_candidates'] += bool(hits)
    counters['variants_with_robust_candidates'] += any(h['genotype_quality'] == 'robust' for h in hits)
    limit = params['top_results']
    records = [dict(hit, assay_rank=rank) for rank, hit in enumerate(hits[:limit] if limit else hits, 1)]
    if not hits and unspecific:      # also when the natural design failed and only dCAPS assays were removed
        records.append(dict(window.base, status='no_specific_primers',
                            reason='Все пары праймеров дают посторонние продукты ПЦР или лежат в повторах генома'))
    elif not hits and failure is not None:
        records.append(failure)
    elif not hits:
        reason = ('Не найдено сочетание для различения всех генотипов локуса' if window.alternates else
                  'Не найдено сочетание для различения трёх генотипов')
        if window.ambiguous:
            reason += f"; окон с неопределёнными основаниями: {window.ambiguous}"
        if window.errors:
            reason += f"; ошибок расчёта ферментов: {window.errors}"
        records.append(dict(window.base, status='no_discriminating_assay', reason=reason))
    run.log(f'[INFO] {window.base["variant"]}: пар {window.base["pairs_evaluated"]}, '
            f"кандидатов {len(hits)}, неразличимых {window.genotype_rejected}")
    return records


@dataclass
class _Search:
    """One SNP after its assays were found and before they are ranked."""
    window: _Window
    masked: bool                     # the window has soft-masked bases
    source_line: int = 0
    records: Optional[list] = None   # final rows of an SNP that ended early
    hits: list = field(default_factory=list)
    failure: Optional[dict] = None
    pairs: list = field(default_factory=list)
    natural_tried: bool = False
    double_tried: bool = False
    dcaps_tried: bool = False
    aspcr_tried: bool = False
    unspecific: int = 0              # assays removed by the specificity check


def _search_snp(run: _Run, fa, chrom: str, mapped_id: str, pos: int, ref: str, alt: str) -> _Search:
    """The assays of one SNP, not ranked yet."""
    window = _snp_window(run, fa, chrom, mapped_id, pos, ref, alt)
    search = _Search(window, masked=window.raw != window.seq)
    ref_end = window.off + len(ref)
    if window.off < 0 or window.seq[window.off:ref_end] != ref:
        search.records = [dict(window.base, status='analysis_error', reason='REF не совпадает с FASTA')]
        return search
    if not window.alternates:
        window.altseq = window.seq[:window.off] + alt + window.seq[ref_end:]
    can_dcaps = not (window.indel or window.alternates)
    window.known = _known_variants(run, window, mapped_id)
    search.failure, pairs = _natural_pairs(run, window)
    if search.failure is not None and (not can_dcaps or run.dcaps == 'off'):
        search.records = [search.failure]
        return search
    window.base['pairs_evaluated'] = len(pairs) if run.use_primer3 else 0
    run.counters['pairs_evaluated'] += window.base['pairs_evaluated']
    if window.indel or (window.alternates and len(set(window.length_changes)) == len(window.length_changes)):
        search.pairs = pairs
        search.hits = _length_hits(run, window, pairs)
    if not any(h['genotype_quality'] == 'robust' for h in search.hits):
        search.hits += _natural_hits(run, window, pairs)
        search.natural_tried = True
    if run.params['double_digest'] and pairs and len(run.known_cut_enzymes) >= 2:
        search.pairs = pairs
        if not any(h['genotype_quality'] == 'robust' for h in search.hits):
            search.hits += _double_hits(run, window, pairs)
            search.double_tried = True
    if can_dcaps and (run.dcaps == 'always' or (run.dcaps == 'fallback'
                            and not any(h['genotype_quality'] == 'robust' for h in search.hits))):
        search.hits += _dcaps_hits(run, window)
        search.dcaps_tried = True
    if can_dcaps and _aspcr_enabled(run, pairs):
        search.pairs = pairs
        if not any(h['genotype_quality'] == 'robust' for h in search.hits):
            search.hits += _aspcr_hits(run, window, pairs)
            search.aspcr_tried = True
    return search


def _aspcr_enabled(run, pairs):
    return bool(run.use_primer3 and run.params['allele_specific'] and pairs)


def _aspcr_hits(run, window, pairs):
    params = run.params
    candidates = design_allele_specific(
        window.raw, window.off, window.ref, window.alt, control_pairs=pairs,
        pick_primers=primer3_pick, parse_pairs=_primer_candidates, model=run.model,
        excluded_regions=window.excluded, mask_repeats=params['mask_repeats'], check_cancel=run.check_cancel,
        **{key: params[key] for key in ('prod_min', 'prod_max', 'tm_min', 'tm_max')})
    hits = []
    for candidate in candidates:
        control = candidate['aspcr_control']
        record = _amplicon_record(window, control, control['start0'], control['end0'], design='allele_specific')
        for key in ('aspcr_control', 'aspcr_reactions', 'aspcr_genotypes'):
            candidate[key] = json.dumps(candidate[key])
        hits.append(dict(record, **candidate, status='ok', pattern='allele_specific',
                         reason='Различимость полос рассчитана. Аллельная специфичность ПЦР не подтверждена'))
    if hits:
        run.counters['aspcr_candidates'] = run.counters.get('aspcr_candidates', 0) + len(hits)
    return hits


def _specificity_key(hit):
    return (hit['primers_left'].upper(), hit['primers_right'].upper(),
            hit.get('_intended', (hit['mapped_id'], hit['amplicon_start'] - 1, hit['amplicon_end'])))


def _scan_hits(run, hits, checked):
    ordinary = []
    assays = []
    for hit in hits:
        if hit.get('design') != 'allele_specific':
            ordinary.append(hit)
            continue
        control = json.loads(hit['aspcr_control'])
        reactions = json.loads(hit['aspcr_reactions'])
        side = hit['aspcr_specific_side']
        origin = hit['amplicon_start'] - control['start0']
        checks = [hit]
        for reaction in reactions.values():
            nested = dict(primers_left=reaction['primers_left'], primers_right=reaction['primers_right'],
                          amplicon_start=origin + reaction['start0'], amplicon_end=origin + reaction['end0'] - 1,
                          mapped_id=hit['mapped_id'])
            # The outer primer on this side also coexists with the inner primer.
            checks.extend((nested, dict(nested, primers_left=control[f'primers_{side}'],
                                        primers_right=reaction[f'primers_{side}'], _intended=None)))
        ordinary.extend(checks)
        assays.append((hit, checks))
    _scan_genome(run, ordinary, checked)
    for hit, checks in assays:
        hit['_unspecific'] = any(check.get('_unspecific') for check in checks)
        products = {product for check in checks for product in check.get('_products', ())}
        sizes = sorted(end - start for _, start, end in products)
        hit['offtarget_sizes'] = '; '.join(map(str, sizes[:5])) + ('; …' if len(sizes) > 5 else '')
        hit['offtarget_products'] = None if any(check.get('offtarget_products') is None for check in checks) else len(products)
        for key in ('primer_sites_left', 'primer_sites_right'):
            values = [check.get(key) for check in (hit, *checks[1::2])]
            hit[key] = max((value for value in values if value is not None), default=None)


def _scan_genome(run: _Run, hits: list, checked: dict) -> None:
    """Look up the primers of the assays in the whole genome and write the result into each assay.

    Sets primer_sites_left/right, offtarget_products, offtarget_sizes and the private
    '_unspecific' (a product elsewhere, or a primer inside a repeat)."""
    params = run.params
    pending = []
    for hit in hits:
        run.check_cancel()
        key = _specificity_key(hit)
        if key in checked:
            hit.update(checked[key])
        else:
            pending.append(hit)
    hits = pending
    primers = sorted({hit[key].upper() for hit in hits for key in ('primers_left', 'primers_right')})
    if not primers:
        return
    run.counters['specificity_scans'] += 1
    run.counters['specificity_primers'] += len(primers)
    run.log(f'[INFO] Специфичность: проход по геному для {len(primers)} праймеров')
    shown = [0]

    def progress(done, total):
        tenth = 10 * done // total
        if tenth > shown[0] and tenth < 10:
            shown[0] = tenth
            run.log(f'[INFO] Специфичность: просмотрено {tenth * 10}% генома')

    try:
        sites = specificity.find_sites(params['fasta'], primers, max_mismatches=params['specificity_mismatches'],
                                       progress_cb=progress, cancel_cb=run.cancel_cb)
    except specificity.ScanCancelled:
        run.report['cancelled'] = True
        raise AnalysisCancelled() from None
    run.check_cancel()
    def count(entry):
        # "More than the limit" for a primer inside a repeat: its list of sites is cut.
        return specificity.MAX_SITES + 1 if entry.saturated else len(entry.sites)

    for hit in hits:
        run.check_cancel()
        left, right, intended = _specificity_key(hit)
        answer = checked.get((left, right, intended))
        if answer is None:
            if left not in sites or right not in sites:
                # A primer the scan cannot look up (letters other than A, C, G, T): left unchecked.
                answer = dict(primer_sites_left=None, primer_sites_right=None, offtarget_products=None,
                              offtarget_sizes='', _unspecific=False)
            elif sites[left].saturated or sites[right].saturated:
                # A repeat: not specific, and the products of a cut list would be an arbitrary number.
                answer = dict(primer_sites_left=count(sites[left]), primer_sites_right=count(sites[right]),
                              offtarget_products=None, offtarget_sizes='', _unspecific=True)
            else:
                products = specificity.pair_products(sites, left, right, params['specificity_max_product'], intended)
                sizes = sorted(product.size for product in products)
                answer = dict(primer_sites_left=count(sites[left]), primer_sites_right=count(sites[right]),
                              offtarget_products=len(products),
                              offtarget_sizes='; '.join(map(str, sizes[:5])) + ('; …' if len(sizes) > 5 else ''),
                              _unspecific=bool(products), _products={product[:3] for product in products})
            checked[left, right, intended] = answer
        hit.update(answer)


def _drop_unspecific(run: _Run, search: _Search) -> None:
    if run.specificity != 'exclude':
        return
    kept = [hit for hit in search.hits if not hit.get('_unspecific')]
    removed = len(search.hits) - len(kept)
    search.unspecific += removed
    run.counters['excluded_offtarget'] += removed
    search.hits = kept


def _check_specificity(run: _Run, searches: list) -> None:
    """Scan in batches and reopen each fallback stage if exclusion removed its predecessor."""
    pending = [search for search in searches if search.records is None]
    checked = {}
    _scan_hits(run, [hit for search in pending for hit in search.hits], checked)
    for search in pending:
        _drop_unspecific(run, search)
    retry_natural = [s for s in pending if s.pairs and not s.natural_tried
                     and not any(h['genotype_quality'] == 'robust' for h in s.hits)]
    added = []
    for search in retry_natural:
        extra = _natural_hits(run, search.window, search.pairs)
        search.natural_tried = True
        search.hits.extend(extra)
        added.extend(extra)
    if added:
        _scan_genome(run, added, checked)
        for search in retry_natural:
            _drop_unspecific(run, search)
    # A robust single may have been excluded. An explicitly requested dCAPS
    # hit is independent of the decision to try the natural double digest.
    retry_double = [s for s in pending if run.params['double_digest'] and len(run.known_cut_enzymes) >= 2
                    and s.pairs and not s.double_tried and
                    not any(h.get('design') in ('natural', 'amplicon_length') and h['genotype_quality'] == 'robust'
                            for h in s.hits)]
    added = []
    for search in retry_double:
        extra = _double_hits(run, search.window, search.pairs)
        search.double_tried = True
        search.hits.extend(extra)
        added.extend(extra)
    if added:
        _scan_genome(run, added, checked)
        for search in retry_double:
            _drop_unspecific(run, search)
    retry = [search for search in pending if not (search.window.indel or search.window.alternates)
             and run.dcaps == 'fallback' and not search.dcaps_tried
             and not any(hit['genotype_quality'] == 'robust' for hit in search.hits)]
    added = []
    if retry:
        run.log(f'[INFO] Специфичность: dCAPS для {len(retry)} SNP без специфичной надёжной схемы')
    for search in retry:
        run.check_cancel()
        extra = _dcaps_hits(run, search.window)
        search.dcaps_tried = True
        search.hits += extra
        added += extra
    if added:
        _scan_genome(run, added, checked)
        for search in retry:
            _drop_unspecific(run, search)
    aspcr_retry = [search for search in pending if not (search.window.indel or search.window.alternates)
                   and not search.aspcr_tried and _aspcr_enabled(run, search.pairs)
                   and not any(hit['genotype_quality'] == 'robust' for hit in search.hits)]
    added = []
    for search in aspcr_retry:
        extra = _aspcr_hits(run, search.window, search.pairs)
        search.aspcr_tried = True
        search.hits.extend(extra)
        added.extend(extra)
    if added:
        _scan_hits(run, added, checked)
        for search in aspcr_retry:
            _drop_unspecific(run, search)
    if run.specificity == 'exclude':
        run.log(f"[INFO] Специфичность: исключено схем {sum(search.unspecific for search in pending)}")


def run_analysis(params, status_cb=None, progress_cb=None, cancel_cb=None, report=None):
    """Design assays for every SNP of params['input_file']; return (header, rows).

    ``report`` (a dict) receives the parameters used, the summary counters and the
    rejected input lines. The window and the command line both call this.
    """
    params = dict(params)
    # Settings saved before 1.5.0 may carry the switch of the removed single-pair mode.
    params.pop('joint_design', None)
    for key in OPTIONAL_PARAMS:
        params.setdefault(key, DEFAULT_PARAMS[key])
    errors = validate_params(params)
    if errors:
        raise ValueError('; '.join(errors))
    model = gel_model_from_params(params)
    report = report if report is not None else {}
    report['parameters'] = dict(params)
    report['scoring_model'] = model.version
    counters = report.setdefault('summary', {})
    counters.update(pairs_evaluated=0, enzyme_evaluations=0, enzyme_errors=0, digest_searches=0,
                    double_snps=0, double_amplicons=0, double_pairs_possible=0, double_pairs_scored=0,
                    ambiguous_templates=0, genotype_rejected=0, candidates=0,
                    variants_with_candidates=0, variants_with_robust_candidates=0,
                    candidate_sites=0, pairs_near_end=0, excluded_primer_site=0,
                    background_variants_loaded=0, background_ref_mismatch=0,
                    excluded_background_site=0, dcaps_snps=0, dcaps_designs=0,
                    dcaps_no_primer=0, dcaps_pairs=0, dcaps_candidates=0,
                    specificity_scans=0, specificity_primers=0, excluded_offtarget=0)

    accepted_lines = []
    row_lines = report['result_source_lines'] = []
    cleaned, _, _, _ = normalize_variants_for_rflp(
        fasta_path=params['fasta'], input_path=params['input_file'],
        names2_path=params.get('names2_file', ''), block_bp=100_000,
        progress_every=100000, status_cb=status_cb,
        skip_norm=params.get('skip_norm', False),
        save_files=not params.get('no_save_norm', True), report=report,
        progress_cb=progress_cb, cancel_cb=cancel_cb,
        vcf_region=params['vcf_region'], vcf_samples=params['vcf_samples'], accepted_lines=accepted_lines,
    )
    header = RESULT_HEADER
    out_rows = []
    if not cleaned:
        counters['rows'] = 0
        return header, out_rows
    use_primer3 = params.get('use_primer3', False)
    codes = parse_suppliers_arg(','.join(params.get('suppliers', [])))
    inventory = set(params['lab_enzymes'])
    if params['enzyme_mode'] == 'only':
        enzymes = [e for e in AllEnzymes if e.__name__ in inventory]
    else:
        preferred = inventory if params['enzyme_mode'] == 'prefer' else set()
        enzymes = [e for e in AllEnzymes if e.__name__ in preferred or
                   ((not codes or enzyme_supplier_codes(e) & codes) and (not params['commercial_only'] or e.suppl))]
    enzymes.sort(key=lambda e: e.__name__)
    run = _Run(params=params, model=model, report=report, counters=counters, enzymes=enzymes,
               # Without a known cut position the fragment lengths are unknown.
               known_cut_enzymes=[e for e in enzymes if not e.is_unknown()],
               use_primer3=use_primer3,
               # Primer pairs and the enzyme are chosen together; without Primer3 the whole window is the amplicon.
               mode='joint' if use_primer3 else 'window',
               dcaps=params['dcaps'] if use_primer3 and _PRIMER3_AVAILABLE else 'off',
               specificity=params['specificity'] if use_primer3 else 'off',
               status_cb=status_cb, cancel_cb=cancel_cb)
    names = load_names2(params.get('names2_file', ''))
    _log_settings(run, len(cleaned))

    # Always release the FASTA handle, including cancellation and exceptions.
    fa = Fasta(params['fasta'], as_raw=True)
    masked_windows = 0

    def contig_of(chrom):
        return resolve_input_contig(chrom, names, lambda name: name in fa, is_vcf_path(params['input_file']))

    try:
        if params['background_vcf']:
            regions = {}
            for chrom, pos, ref, _ in cleaned:
                regions.setdefault(contig_of(chrom), []).append(
                    (pos - params['flank'], pos + len(ref) - 1 + params['flank']))
            run.background, records = load_background_variants(
                params['background_vcf'], regions, names, cancel_cb, fasta_contigs=fa)
            counters['background_variants_loaded'] = sum(map(len, run.background.values()))
            run.log(f"[INFO] Фоновые варианты: {counters['background_variants_loaded']} в окнах SNP "
                    f"из {records} записей VCF")

        def emit(search):
            run.check_cancel()
            records = search.records
            if records is None:
                records = _snp_records(run, search.window, search.hits, search.failure, search.unspecific)
            out_rows.extend([record.get(key) for key in header] for record in records)
            row_lines.extend([search.source_line] * len(records))

        searches = []
        for vi, (chrom, pos, ref, alt) in enumerate(cleaned, 1):
            run.check_cancel()
            search = _search_snp(run, fa, chrom, contig_of(chrom), pos, ref, alt)
            search.source_line = accepted_lines[vi - 1]
            masked_windows += search.masked
            if run.specificity == 'off':
                emit(search)             # nothing to wait for
            else:
                searches.append(search)
            if progress_cb:
                progress_cb(vi, len(cleaned))
            if len(searches) >= SPECIFICITY_BATCH or (searches and vi == len(cleaned)):
                # The primers of these SNPs are looked up in one pass over the genome.
                _check_specificity(run, searches)
                for search in searches:
                    emit(search)
                searches = []
    finally:
        fa.close()
    if params['mask_repeats'] and use_primer3 and not masked_windows:
        run.log('[WARN] В окнах SNP нет строчных букв: FASTA без soft-mask, '
                'маскирование повторов не действует.')
    run.check_cancel()
    counters['rows'] = len(out_rows)
    return header, out_rows


# The name of this function up to 1.5.0; tools/regression.py calls it in older source trees too.
run_rflp_gui_mode = run_analysis
