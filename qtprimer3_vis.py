import sys
import os
import csv
import re
import traceback
import logging
from typing import List, Dict, Tuple, Optional
from collections import namedtuple
from datetime import datetime

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QPushButton, QLabel, QFileDialog, QTableWidget, QTableWidgetItem, QLineEdit, QCheckBox,
    QSpinBox, QDoubleSpinBox, QTextEdit, QProgressBar, QMessageBox, QComboBox, QHeaderView,
    QSplitter, QGraphicsView, QGraphicsScene, QSizePolicy, QDialog, QGroupBox, QScrollArea
)
from PyQt6.QtCore import QThread, pyqtSignal, Qt
from PyQt6 import QtCore
from PyQt6 import QtGui  # QColor
from PyQt6.QtGui import QIcon
import resources_rc

from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio.Restriction import AllEnzymes

from openpyxl import Workbook
import rflp_core as core_impl
from reporting import export_csv as safe_export_csv, export_excel as safe_export_excel, collect_versions
import visualization as enhanced_visualization
from result_model import LABELS


class NumericTableItem(QTableWidgetItem):
    """Sort numeric results numerically while retaining human-readable text."""
    def __lt__(self, other):
        left = self.text().strip()
        right = other.text().strip() if other is not None else ""
        raw_left = self.data(Qt.ItemDataRole.UserRole)
        raw_right = other.data(Qt.ItemDataRole.UserRole) if other is not None else None
        if isinstance(raw_left, (int, float)) and isinstance(raw_right, (int, float)):
            return raw_left < raw_right
        try:
            return float(left) < float(right)
        except (ValueError, TypeError):
            return left.casefold() < right.casefold()

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


def parse_input(path: str) -> List[Tuple[str, int, str, str]]:
    rows: List[Tuple[str, int, str, str]] = []
    try:
        if not os.path.isfile(path):
            logger.error(f"Файл с вариантами не найден: {path}")
            return rows
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            first = f.readline()
            if not first:
                return rows
            delim = detect_delim(first)
            hdr = [h.strip().lower() for h in first.strip("\r\n").split(delim)]

            def idx(*cands):
                for c in cands:
                    if c in hdr:
                        return hdr.index(c)
                return None

            i_chrom = idx("chrom", "chr", "chr1", "chromosome")
            i_pos = idx("pos", "position")
            i_ref = idx("ref", "reference")
            i_alt = idx("alt", "alternate", "alt_allele")

            def parse_line(line: str):
                parts = line.strip("\r\n").split(delim)
                if len(parts) < 4:
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
                    return None

            data: List[Tuple[str, int, str, str]] = []
            if any(
                x in hdr
                for x in ("chrom", "chr", "pos", "ref", "alt", "position", "reference", "alternate")
            ):
                for line in f:
                    if not line.strip():
                        continue
                    if "должно быть:" in line.lower():
                        continue
                    rec = parse_line(line)
                    if rec:
                        data.append(rec)
            else:
                rec = parse_line(first)
                if rec:
                    data.append(rec)
                for line in f:
                    if not line.strip():
                        continue
                    rec = parse_line(line)
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
    variants_raw = parse_input(input_path)
    total = len(variants_raw)

    if skip_norm:
        summary_lines = [
            "=== NORMALIZE SUMMARY ===",
            "normalization_skipped          : 1",
            f"variants_passed_without_change : {total}",
        ]
        if status_cb:
            status_cb("\n" + "\n".join(summary_lines) + "\n")
        return variants_raw, summary_lines, "", ""

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
    if not _PRIMER3_AVAILABLE:
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
    except Exception:
        return None


class RFLPCalcThread(QThread):
    finished = pyqtSignal(list, list)
    status = pyqtSignal(str)
    progress = pyqtSignal(int, int)

    def __init__(self, params):
        super().__init__()
        self.params = params
        self.report = {}

    def run(self):
        try:
            header, out_rows = core_impl.run_rflp_gui_mode(
                self.params, status_cb=self.status.emit,
                progress_cb=self.progress.emit,
                cancel_cb=self.isInterruptionRequested,
                report=self.report,
            )
            self.finished.emit(header, out_rows)
        except core_impl.AnalysisCancelled:
            self.status.emit("\n[INFO] Анализ отменён пользователем.")
            self.finished.emit([], [])
        except Exception:
            tb = traceback.format_exc()
            self.status.emit(f"\n[ERROR] Exception in calculation thread:\n{tb}")
            logger.error(tb)
            self.finished.emit([], [])


def run_rflp_gui_mode(params, status_cb=None, **kwargs):
    return core_impl.run_rflp_gui_mode(params, status_cb=status_cb, **kwargs)


def _safe_int(s: object) -> Optional[int]:
    try:
        if s is None:
            return None
        t = str(s).strip()
        if t == "" or t.lower() == "none":
            return None
        return int(float(t))
    except Exception:
        return None



class PrimerGraphicPanel(QWidget):
    """Improved visualization panel.

    Features:
      - Two tracks (REF and ALT) for the same amplicon coordinates.
      - Primer arrows (left ->, right <-), SNP marker.
      - Restriction cut sites for the selected enzyme on each track (REF vs ALT).
      - Sequence view with highlighting, theme-adaptive (light/dark).

    The panel is updated when a result row is selected in the table.
    """

    def __init__(self, parent=None):
        super().__init__(parent)

        self._last_payload: Optional[dict] = None

        layout = QVBoxLayout(self)

        self.title = QLabel("Визуализация RFLP: REF/ALT, праймеры, сайты рестрикции")
        tfont = self.title.font()
        tfont.setPointSize(max(10, tfont.pointSize() + 2))
        tfont.setBold(True)
        self.title.setFont(tfont)

        self.meta = QLabel("")
        self.meta.setWordWrap(True)

        self.scene = QGraphicsScene(self)
        self.view = QGraphicsView(self.scene)
        self.view.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        self.view.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.view.setMinimumHeight(220)

        self.seq_view = QTextEdit()
        self.seq_view.setReadOnly(True)
        mono = QtGui.QFont("Consolas")
        mono.setStyleHint(QtGui.QFont.StyleHint.Monospace)
        mono.setPointSize(10)
        self.seq_view.setFont(mono)
        self.seq_view.setMinimumHeight(180)

        layout.addWidget(self.title)
        layout.addWidget(self.meta)
        layout.addWidget(self.view, stretch=1)
        layout.addWidget(self.seq_view)

        self._apply_theme()
        self.show_empty()

    # ----- theme helpers -----
    def _is_dark(self) -> bool:
        pal = self.palette()
        w = pal.color(QtGui.QPalette.ColorRole.Window)
        b = pal.color(QtGui.QPalette.ColorRole.Base)
        return (w.lightness() + b.lightness()) / 2.0 < 128

    def _colors(self):
        pal = self.palette()
        dark = self._is_dark()
        base = pal.color(QtGui.QPalette.ColorRole.Base)
        window = pal.color(QtGui.QPalette.ColorRole.Window)
        text = pal.color(QtGui.QPalette.ColorRole.Text)

        # Use palette-aware colors with explicit fallbacks.
        axis = QtGui.QColor(text)
        axis.setAlpha(210)

        grid = QtGui.QColor(text)
        grid.setAlpha(120)

        ref_track = QtGui.QColor(text)
        ref_track.setAlpha(220)

        alt_track = QtGui.QColor(text)
        alt_track.setAlpha(200)

        if dark:
            primer_left = QtGui.QColor("#3ddc84")
            primer_right = QtGui.QColor("#6aa9ff")
            snp = QtGui.QColor("#ffb020")
            cut_ref = QtGui.QColor("#ff7aa2")
            cut_alt = QtGui.QColor("#8de2ff")
            hl_left_bg = "#1e5b3b"
            hl_right_bg = "#1b3e6f"
            hl_snp_bg = "#6a4a00"
            hl_diff_bg = "#5a1f33"
        else:
            primer_left = QtGui.QColor("#1f9d55")
            primer_right = QtGui.QColor("#2b6cb0")
            snp = QtGui.QColor("#d38b00")
            cut_ref = QtGui.QColor("#c53030")
            cut_alt = QtGui.QColor("#1a7f9a")
            hl_left_bg = "#d7ffd7"
            hl_right_bg = "#d7e7ff"
            hl_snp_bg = "#fff2a8"
            hl_diff_bg = "#ffd7e1"

        return {
            "dark": dark,
            "base": base,
            "window": window,
            "text": text,
            "axis": axis,
            "grid": grid,
            "ref_track": ref_track,
            "alt_track": alt_track,
            "primer_left": primer_left,
            "primer_right": primer_right,
            "snp": snp,
            "cut_ref": cut_ref,
            "cut_alt": cut_alt,
            "hl_left_bg": hl_left_bg,
            "hl_right_bg": hl_right_bg,
            "hl_snp_bg": hl_snp_bg,
            "hl_diff_bg": hl_diff_bg,
        }

    def _apply_theme(self):
        c = self._colors()
        # Graphics view background should follow Window.
        self.view.setBackgroundBrush(QtGui.QBrush(c["window"]))
        # Sequence view follows Base with native palette, but keep a mild border.
        border = "#3a3a3a" if c["dark"] else "#c0c0c0"
        self.seq_view.setStyleSheet(f"QTextEdit {{ border: 1px solid {border}; }}")

    def changeEvent(self, event):
        # React to palette changes (theme switches).
        try:
            if event.type() == QtCore.QEvent.Type.PaletteChange:
                self._apply_theme()
                if self._last_payload:
                    self._render_impl(**self._last_payload)
        except Exception:
            pass
        super().changeEvent(event)

    def resizeEvent(self, event):
        """Адаптация графики при изменении размера панели/окна."""
        try:
            # fitInView keeps the diagram readable when the window is resized
            if self.scene is not None and self.view is not None:
                self.view.fitInView(self.scene.sceneRect(), Qt.AspectRatioMode.KeepAspectRatio)
        except Exception:
            pass
        super().resizeEvent(event)


    # ----- public API -----
    def show_empty(self):
        self._last_payload = None
        self.meta.setText(
            "Выберите строку в таблице, затем откройте окно визуализации или обновите выбранную строку."
        )
        self.scene.clear()
        self.seq_view.setPlainText("(строка не выбрана)")

    def show_message(self, msg: str):
        self._last_payload = None
        self.meta.setText(msg)
        self.scene.clear()
        self.seq_view.setPlainText(msg)

    def render(
        self,
        mapped_id: str,
        genome_start: int,
        genome_end: int,
        amplicon_seq: str,
        snp_offset: Optional[int],
        primer_left: str,
        primer_right: str,
        tm_left: Optional[float],
        tm_right: Optional[float],
        primer_left_start: Optional[int] = None,
        primer_left_len: Optional[int] = None,
        primer_right_start: Optional[int] = None,
        primer_right_len: Optional[int] = None,
        variant_str: str = "",
        enzyme_name: str = "",
        pattern: str = "",
        frags_ref: str = "",
        frags_alt: str = "",
    ):
        payload = dict(
            mapped_id=mapped_id,
            genome_start=genome_start,
            genome_end=genome_end,
            amplicon_seq=amplicon_seq,
            snp_offset=snp_offset,
            primer_left=primer_left,
            primer_right=primer_right,
            tm_left=tm_left,
            tm_right=tm_right,
            primer_left_start=primer_left_start,
            primer_left_len=primer_left_len,
            primer_right_start=primer_right_start,
            primer_right_len=primer_right_len,
            variant_str=variant_str,
            enzyme_name=enzyme_name,
            pattern=pattern,
            frags_ref=frags_ref,
            frags_alt=frags_alt,
        )
        self._last_payload = payload
        self._render_impl(**payload)

    # ----- internal helpers -----
    def _draw_text(self, x: float, y: float, s: str, size: int = 9, bold: bool = False):
        c = self._colors()
        item = self.scene.addText(s)
        f = item.font()
        f.setPointSize(size)
        f.setBold(bold)
        item.setFont(f)
        item.setDefaultTextColor(c["text"])
        item.setPos(x, y)
        return item

    @staticmethod
    def _html_escape(s: str) -> str:
        return (
            s.replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;")
            .replace('"', "&quot;")
        )

    @staticmethod
    def _parse_variant(variant_str: str) -> Tuple[str, str]:
        # Expected examples:
        #   "1H:12345 A>G"
        #   "chr1:100 C>TT" (we will treat as non-SNP)
        s = (variant_str or "").strip()
        m = re.search(r"\s([ACGTN]+)>([ACGTN\-]+)", s)
        if not m:
            return "", ""
        return m.group(1).upper(), m.group(2).upper()

    @staticmethod
    def _enzyme_by_name(name: str):
        nm = (name or "").strip()
        if not nm:
            return None
        for e in AllEnzymes:
            if getattr(e, "__name__", "") == nm:
                return e
        return None

    def _map_primer_positions(
        self,
        seq: str,
        primer_left: str,
        primer_right: str,
        primer_left_start: Optional[int],
        primer_left_len: Optional[int],
        primer_right_start: Optional[int],
        primer_right_len: Optional[int],
    ) -> Tuple[int, int]:
        L = len(seq)
        pl = (primer_left or "").upper().strip()
        pr = (primer_right or "").upper().strip()
        pr_bind = revcomp(pr) if pr else ""

        left_pos = seq.find(pl) if pl else -1
        right_pos = seq.find(pr_bind) if pr_bind else -1

        if left_pos < 0 and pl and primer_left_start is not None and primer_left_len is not None:
            # In current pipeline amplicon is often aligned so that left primer starts at 0.
            left_pos = 0

        if right_pos < 0 and pr_bind and primer_right_start is not None and primer_left_start is not None:
            # Conservative fallback: put it near the right edge.
            right_pos = max(0, min(L - len(pr_bind), L - len(pr_bind)))

        return left_pos, right_pos

    def _make_alt_seq(self, seq_ref: str, snp_offset: Optional[int], ref: str, alt: str) -> Tuple[Optional[str], str]:
        if snp_offset is None:
            return None, "No SNP offset."
        if not (0 <= snp_offset < len(seq_ref)):
            return None, "SNP offset out of range."

        # Only handle simple SNP visualization (1 bp -> 1 bp).
        if len(ref) != 1 or len(alt) != 1:
            return None, "Non-SNP (indel/multi-bp) visualization is not supported in two-track mode."

        genome_base = seq_ref[snp_offset].upper()
        if ref and genome_base != ref:
            # Still allow, but note mismatch.
            pass

        seq_alt = list(seq_ref)
        seq_alt[snp_offset] = alt
        return "".join(seq_alt), ""

    def _render_impl(
        self,
        mapped_id: str,
        genome_start: int,
        genome_end: int,
        amplicon_seq: str,
        snp_offset: Optional[int],
        primer_left: str,
        primer_right: str,
        tm_left: Optional[float],
        tm_right: Optional[float],
        primer_left_start: Optional[int],
        primer_left_len: Optional[int],
        primer_right_start: Optional[int],
        primer_right_len: Optional[int],
        variant_str: str,
        enzyme_name: str,
        pattern: str,
        frags_ref: str,
        frags_alt: str,
    ):
        self._apply_theme()
        c = self._colors()

        seq_ref = (amplicon_seq or "").upper()
        if not seq_ref:
            self.show_message("Failed to obtain amplicon sequence.")
            return

        L = len(seq_ref)
        ref_allele, alt_allele = self._parse_variant(variant_str)
        seq_alt, alt_err = self._make_alt_seq(seq_ref, snp_offset, ref_allele, alt_allele)

        # Primer mapping on REF sequence (positions should be the same for ALT for simple SNP)
        left_pos, right_pos = self._map_primer_positions(
            seq_ref,
            primer_left,
            primer_right,
            primer_left_start,
            primer_left_len,
            primer_right_start,
            primer_right_len,
        )

        pl = (primer_left or "").upper().strip()
        pr = (primer_right or "").upper().strip()
        pr_bind = revcomp(pr) if pr else ""

        # Prepare restriction cuts (optional)
        enz = self._enzyme_by_name(enzyme_name)
        cuts_ref = []
        cuts_alt = []
        if enz is not None:
            try:
                cuts_ref = cut_positions(enz, seq_ref)
                if seq_alt is not None:
                    cuts_alt = cut_positions(enz, seq_alt)
            except Exception:
                cuts_ref = []
                cuts_alt = []

        # ----- meta -----
        meta_parts = [
            f"Континг: {mapped_id}",
            f"Ампликон: {genome_start}-{genome_end}",
            f"Длина: {L} п.н.",
        ]
        if variant_str:
            meta_parts.append(f"Вариант: {variant_str}")
        if snp_offset is not None:
            meta_parts.append(f"SNP offset: {snp_offset} (0-based, ампликон)")
        if enzyme_name:
            meta_parts.append(f"Фермент: {enzyme_name}")
        if pattern:
            meta_parts.append(f"Паттерн: {pattern}")
        if frags_ref:
            meta_parts.append(f"Фрагменты REF: {frags_ref}")
        if frags_alt:
            meta_parts.append(f"Фрагменты ALT: {frags_alt}")
        if pl:
            meta_parts.append(f"Левый праймер (Tm={tm_left if tm_left is not None else '---'}): {pl}")
        if pr:
            meta_parts.append(f"Правый праймер (Tm={tm_right if tm_right is not None else '---'}): {pr}")

        warn = []
        if pl and left_pos < 0:
            warn.append("Левый праймер не найден внутри ампликона.")
        if pr and right_pos < 0:
            warn.append("Сайт посадки правого праймера (revcomp) не найден внутри ампликона.")
        if alt_err:
            warn.append(alt_err)
        self.meta.setText(
            "  |  ".join(meta_parts)
            + ("\n" + "\n".join(warn) if warn else "")
        )

        # ----- graphics -----
        self.scene.clear()

        # Scene geometry (fixed logical canvas, then scaled by fitInView)
        W, H = 1200.0, 360.0
        self.scene.setSceneRect(0, 0, W, H)

        margin_l = 70.0
        margin_r = 40.0
        usable = max(1.0, W - margin_l - margin_r)

        y_ref = 150.0
        y_alt = 240.0

        def x_from_bp(bp: int) -> float:
            bp = max(0, min(L, int(bp)))
            return margin_l + (bp / max(1, L)) * usable

        # Tracks
        pen_track = QtGui.QPen(c["axis"]) 
        pen_track.setWidth(3)
        self.scene.addLine(margin_l, y_ref, W - margin_r, y_ref, pen_track)
        self.scene.addLine(margin_l, y_alt, W - margin_r, y_alt, pen_track)

        self._draw_text(margin_l, y_ref - 30, "REF", size=10, bold=True)
        self._draw_text(margin_l, y_alt - 30, "ALT", size=10, bold=True)

        # Tick marks
        if L <= 400:
            tick_step = 50
        elif L <= 1200:
            tick_step = 100
        else:
            tick_step = 200

        pen_tick = QtGui.QPen(c["grid"]) 
        pen_tick.setWidth(1)

        for bp in range(0, L + 1, tick_step):
            x = x_from_bp(bp)
            self.scene.addLine(x, y_ref - 7, x, y_ref + 7, pen_tick)
            self.scene.addLine(x, y_alt - 7, x, y_alt + 7, pen_tick)
            if bp == 0 or bp == L or (bp % (tick_step * 2) == 0):
                self._draw_text(x - 12, y_alt + 16, str(bp), size=8)

        # Genome coordinate labels
        self._draw_text(margin_l, 50.0, f"Геном: начало {genome_start}", size=9)
        self._draw_text(W - margin_r - 220, 50.0, f"Геном: конец {genome_end}", size=9)

        # SNP marker + explicit REF->ALT change
        if snp_offset is not None and 0 <= snp_offset <= L:
            x = x_from_bp(snp_offset)
            pen_snp = QtGui.QPen(c["snp"]) 
            pen_snp.setWidth(3)
            self.scene.addLine(x, y_ref - 45, x, y_alt + 45, pen_snp)
            self._draw_text(x + 8, 70.0, "SNP", size=10, bold=True)

            if ref_allele and alt_allele:
                self._draw_text(x + 40, 70.0, f"{ref_allele} → {alt_allele}", size=10, bold=True)

        # Helper: draw nicer arrow with centered label
        def draw_arrow(x1: float, x2: float, y: float, color: QtGui.QColor, label: str, direction: str):
            # main line
            pen = QtGui.QPen(color)
            pen.setWidth(4)
            self.scene.addLine(x1, y, x2, y, pen)

            # arrow head
            head = 10.0
            pen_head = QtGui.QPen(color)
            pen_head.setWidth(3)
            if direction == ">":
                self.scene.addLine(x2, y, x2 - head, y - 7, pen_head)
                self.scene.addLine(x2, y, x2 - head, y + 7, pen_head)
            else:
                self.scene.addLine(x2, y, x2 + head, y - 7, pen_head)
                self.scene.addLine(x2, y, x2 + head, y + 7, pen_head)

            # label near the arrow midpoint
            mx = (x1 + x2) / 2.0
            lx = max(margin_l, min(mx - 80, W - margin_r - 200))
            self._draw_text(lx, y - 22, label, size=9)

        # Primers (shown once, aligned to REF track)
        if pl and left_pos >= 0:
            lp_start = left_pos
            lp_end = left_pos + len(pl)
            draw_arrow(
                x_from_bp(lp_start),
                x_from_bp(lp_end),
                y_ref - 55,
                c["primer_left"],
                f"Левый праймер ({len(pl)} п.н.)",
                ">",
            )

        if pr_bind and right_pos >= 0:
            rp_start = right_pos
            rp_end = right_pos + len(pr_bind)
            draw_arrow(
                x_from_bp(rp_end),
                x_from_bp(rp_start),
                y_ref - 55,
                c["primer_right"],
                f"Правый праймер ({len(pr)} п.н.)",
                "<",
            )

        # Restriction cut sites
        def draw_cuts(cuts: List[int], y: float, color: QtGui.QColor):
            if not cuts:
                return
            pen = QtGui.QPen(color)
            pen.setWidth(2)
            for bp in cuts:
                if bp <= 0 or bp >= L:
                    continue
                x = x_from_bp(bp)
                self.scene.addLine(x, y - 16, x, y + 16, pen)

        if enz is not None:
            draw_cuts(cuts_ref, y_ref, c["cut_ref"])
            if seq_alt is not None:
                draw_cuts(cuts_alt, y_alt, c["cut_alt"])

            # Legend with real color swatches (not hex codes)
            legend_y = y_alt + 55
            box = 12

            # REF swatch
            self.scene.addRect(margin_l, legend_y - box + 2, box, box, QtGui.QPen(QtCore.Qt.PenStyle.NoPen), QtGui.QBrush(c["cut_ref"]))
            self._draw_text(margin_l + box + 6, legend_y - box + 2, "Сайты рестрикции (REF)", size=9)

            # ALT swatch
            self.scene.addRect(margin_l + 220, legend_y - box + 2, box, box, QtGui.QPen(QtCore.Qt.PenStyle.NoPen), QtGui.QBrush(c["cut_alt"]))
            self._draw_text(margin_l + 220 + box + 6, legend_y - box + 2, "Сайты рестрикции (ALT)", size=9)

        # Mini sequence window directly on the diagram (around SNP)
        if snp_offset is not None and seq_alt is not None and 0 <= snp_offset < L:
            win = 20
            a2 = max(0, snp_offset - win)
            b2 = min(L, snp_offset + win + 1)
            ref_win = seq_ref[a2:b2]
            alt_win = seq_alt[a2:b2]
            caret = " " * (snp_offset - a2) + "^"
            # HTML box
            html = (
                "<div style='font-family:Consolas,monospace; font-size:10pt;'>"
                f"<b>Окно вокруг SNP</b> (amp[{a2}..{b2-1}])<br/>"
                f"REF: {ref_win}<br/>"
                f"ALT: {alt_win}<br/>"
                f"     {caret}"
                "</div>"
            )
            t = self.scene.addText("")
            t.setHtml(html)
            # place near SNP but keep inside canvas
            x0 = x_from_bp(snp_offset)
            tx = min(max(margin_l, x0 - 180), W - margin_r - 360)
            ty = 95.0
            t.setPos(tx, ty)

        self.view.fitInView(self.scene.sceneRect(), Qt.AspectRatioMode.KeepAspectRatio)

        # ----- sequence view -----
        # Render sequences in a structured, readable way (no "мешанина").
        # If the amplicon is huge, show a focused window around SNP/primers.

        max_full = 1500  # full render threshold for readability
        if L > max_full:
            center = snp_offset if snp_offset is not None else (left_pos if left_pos >= 0 else (right_pos if right_pos >= 0 else L // 2))
            win = 240
            a = max(0, int(center) - win)
            b = min(L, int(center) + win)
            prefix = f"ПРИМЕЧАНИЕ: ампликон {L} п.н. Показано окно amp[{a}..{b}] вокруг ключевой позиции." 
        else:
            a = 0
            b = L
            prefix = ""

        seq_ref_view = seq_ref[a:b]
        seq_alt_view = seq_alt[a:b] if seq_alt is not None else ""

        left_span = (left_pos, left_pos + len(pl)) if pl and left_pos >= 0 else None
        right_span = (right_pos, right_pos + len(pr_bind)) if pr_bind and right_pos >= 0 else None

        def in_span(i: int, span: Optional[Tuple[int, int]]) -> bool:
            return span is not None and span[0] <= i < span[1]

        line_len = 70
        # Build HTML
        html_lines = []
        html_lines.append("<div style='font-family:Consolas,monospace; font-size:10pt;'>")
        if prefix:
            html_lines.append(f"<div style='opacity:0.85; margin-bottom:6px;'><b>{self._html_escape(prefix)}</b></div>")

        html_lines.append(
            "<div style='margin-bottom:6px; opacity:0.9;'>"
            "Подсветка: <span style='background-color:{};'>левый праймер</span>, "
            "<span style='background-color:{};'>правый праймер</span>, "
            "<span style='background-color:{}; font-weight:700;'>SNP</span>, "
            "<span style='background-color:{};'>различия ALT</span>."
            "</div>".format(c['hl_left_bg'], c['hl_right_bg'], c['hl_snp_bg'], c['hl_diff_bg'])
        )

        for off in range(a, b, line_len):
            chunk_ref = seq_ref[off: off + line_len]
            chunk_alt = seq_alt[off: off + line_len] if seq_alt is not None else ""

            g_from = genome_start + off
            g_to = genome_start + min(L, off + line_len) - 1

            header = f"amp[{off:>6d}..{min(L-1, off+line_len-1):<6d}]  genome[{g_from}..{g_to}]"
            html_lines.append(f"<div style='margin-top:10px; opacity:0.85;'><b>{self._html_escape(header)}</b></div>")

            def render_line(seq: str, label: str, is_alt: bool):
                parts = []
                for j, base in enumerate(seq):
                    i = off + j
                    styles = []
                    if snp_offset is not None and i == snp_offset:
                        styles.append(f"background-color:{c['hl_snp_bg']}; font-weight:700;")
                    elif in_span(i, left_span):
                        styles.append(f"background-color:{c['hl_left_bg']};")
                    elif in_span(i, right_span):
                        styles.append(f"background-color:{c['hl_right_bg']};")

                    if is_alt and seq_alt is not None and i < L:
                        if seq_alt[i] != seq_ref[i]:
                            styles.append(f"background-color:{c['hl_diff_bg']};")

                    st = "".join(styles)
                    if st:
                        parts.append(f"<span style='{st}'>{base}</span>")
                    else:
                        parts.append(base)

                html_lines.append(f"<div><b>{label}:</b> {''.join(parts)}</div>")

            render_line(chunk_ref, "REF", is_alt=False)
            if seq_alt is not None:
                render_line(chunk_alt, "ALT", is_alt=True)

            # markers line
            marks = [" "] * len(chunk_ref)
            if left_span is not None:
                for i2 in range(max(left_span[0], off), min(left_span[1], off + len(chunk_ref))):
                    marks[i2 - off] = ">"
            if right_span is not None:
                for i2 in range(max(right_span[0], off), min(right_span[1], off + len(chunk_ref))):
                    marks[i2 - off] = "<"
            if snp_offset is not None and off <= snp_offset < off + len(chunk_ref):
                marks[snp_offset - off] = "^"

            html_lines.append(f"<div style='opacity:0.75;'><b>    :</b> {self._html_escape(''.join(marks))}</div>")

        html_lines.append("</div>")

        html = "\n".join(html_lines)
        self.seq_view.setHtml(html)


class RFLPVisualizationDialog(QDialog):
    """Отдельное окно визуализации (адаптивное, масштабируемое).

    Открывается по кнопке из главного окна.
    Содержит графику (две дорожки REF/ALT) и вывод последовательностей с подсветкой.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Визуализация RFLP (REF / ALT)")
        self.resize(1200, 850)
        self.setMinimumSize(900, 600)

        layout = QVBoxLayout(self)

        self.panel = PrimerGraphicPanel(self)
        # В диалоге панель должна занимать всё пространство
        self.panel.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(self.panel, stretch=1)

        btn_row = QHBoxLayout()
        btn_row.addStretch(1)
        self.close_btn = QPushButton("Закрыть")
        self.close_btn.clicked.connect(self.close)
        btn_row.addWidget(self.close_btn)
        layout.addLayout(btn_row)

    def show_payload(self, payload: Optional[dict], message_if_empty: str = ""):
        if payload:
            self.panel.render(**payload)
        elif message_if_empty:
            self.panel.show_message(message_if_empty)
        else:
            self.panel.show_empty()


class RFLPPickerGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("RFLP Picker 1.1 — совместный подбор")

        self.setWindowIcon(QIcon(":/icons/app_icon.ico"))

        self.resize(1400, 800)

        central = QWidget()
        main_layout = QVBoxLayout()

        # ===== верх: параметры + лог =====
        top_layout = QHBoxLayout()

        # левая панель
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)

        left_layout.addWidget(QLabel("Файл вариантов (CSV/TSV):"))
        self.input_file = QLineEdit()
        left_layout.addWidget(self.input_file)
        btn_input = QPushButton("Выбрать файл вариантов")
        btn_input.clicked.connect(lambda: self.select_file(self.input_file))
        left_layout.addWidget(btn_input)

        left_layout.addWidget(QLabel("FASTA (геном, папка ./genoms рядом с exe):"))
        self.fasta_combo = QComboBox()
        self.fasta_combo.setEditable(True)
        self.fasta_combo.setToolTip("Можно выбрать локальный файл из genoms или ввести полный путь к FASTA.")
        left_layout.addWidget(self.fasta_combo)

        left_layout.addWidget(QLabel("Файл соответствий (names, папка ./names рядом с exe):"))
        self.names_combo = QComboBox()
        self.names_combo.setEditable(True)
        self.names_combo.setToolTip("Можно выбрать файл из names или ввести полный путь; поле можно оставить пустым.")
        left_layout.addWidget(self.names_combo)

        self.fasta_combo.currentIndexChanged.connect(self.update_names_selection)

        grid = QGridLayout()
        left_layout.addLayout(grid)

        self.flank_spin = QSpinBox()
        self.flank_spin.setRange(1, 5000)
        self.flank_spin.setValue(300)

        self.min_frag_spin = QSpinBox()
        self.min_frag_spin.setRange(10, 5000)
        self.min_frag_spin.setValue(80)

        self.max_frag_spin = QSpinBox()
        self.max_frag_spin.setRange(10, 5000)
        self.max_frag_spin.setValue(800)

        self.delta_spin = QSpinBox()
        self.delta_spin.setRange(1, 1000)
        self.delta_spin.setValue(25)

        self.max_cuts_spin = QSpinBox()
        self.max_cuts_spin.setRange(1, 10)
        self.max_cuts_spin.setValue(3)

        self.gloss_chk = QCheckBox("Только gain/loss")
        self.gloss_chk.setChecked(False)

        self.prod_min_spin = QSpinBox()
        self.prod_min_spin.setRange(50, 2000)
        self.prod_min_spin.setValue(250)

        self.prod_max_spin = QSpinBox()
        self.prod_max_spin.setRange(50, 2000)
        self.prod_max_spin.setValue(600)

        self.tm_min_spin = QDoubleSpinBox()
        self.tm_min_spin.setRange(40.0, 80.0)
        self.tm_min_spin.setDecimals(1)
        self.tm_min_spin.setValue(58.0)

        self.tm_max_spin = QDoubleSpinBox()
        self.tm_max_spin.setRange(40.0, 80.0)
        self.tm_max_spin.setDecimals(1)
        self.tm_max_spin.setValue(62.0)

        self.use_primer3_chk = QCheckBox("Использовать primer3")
        self.use_primer3_chk.setChecked(_PRIMER3_AVAILABLE)

        params_labels = [
            QLabel("Фланк (bp):"),
            QLabel("Мин. фрагмент (bp):"),
            QLabel("Макс. фрагмент (bp):"),
            QLabel("Диагностическая дельта (bp):"),
            QLabel("Макс. разрезов в ампликоне:"),
            QLabel("Опция:"),
            QLabel("Мин. размер продукта (bp):"),
            QLabel("Макс. размер продукта (bp):"),
            QLabel("Мин. Tm (°C):"),
            QLabel("Макс. Tm (°C):"),
            QLabel("Primer3:"),
        ]

        params_widgets = [
            self.flank_spin,
            self.min_frag_spin,
            self.max_frag_spin,
            self.delta_spin,
            self.max_cuts_spin,
            self.gloss_chk,
            self.prod_min_spin,
            self.prod_max_spin,
            self.tm_min_spin,
            self.tm_max_spin,
            self.use_primer3_chk,
        ]

        supplier_widgets: List[QCheckBox] = []
        supplier_codes = sorted(_SUPPLIER_CANON.keys())
        for code in supplier_codes:
            name = _SUPPLIER_CANON[code][0].capitalize()
            cb = QCheckBox(f"{code}: {name}")
            supplier_widgets.append(cb)

        max_len = max(len(params_widgets), len(supplier_widgets))
        for i in range(max_len):
            if i < len(params_widgets):
                grid.addWidget(params_labels[i], i, 0)
                grid.addWidget(params_widgets[i], i, 1)
            if i < len(supplier_widgets):
                col = 2 if i < (len(supplier_widgets) + 1) // 2 else 3
                row_sup = i if col == 2 else i - (len(supplier_widgets) + 1) // 2
                grid.addWidget(supplier_widgets[i], row_sup, col)

        self.supplier_checks = supplier_widgets

        design_box = QGroupBox('Совместный подбор и различимость генотипов')
        design_grid = QGridLayout(design_box)
        self.joint_design_chk = QCheckBox('Подбирать праймеры и фермент совместно')
        self.joint_design_chk.setChecked(True)
        self.joint_design_chk.setToolTip('Отключите для сравнения с поиском по первой паре Primer3.')
        design_grid.addWidget(self.joint_design_chk, 0, 0, 1, 4)
        self.primer_pairs_spin = QSpinBox()
        self.primer_pairs_spin.setRange(1, 100)
        self.primer_pairs_spin.setValue(10)
        self.top_results_spin = QSpinBox()
        self.top_results_spin.setRange(0, 10000)
        self.top_results_spin.setValue(10)
        self.top_results_spin.setSpecialValueText('Все')
        self.gel_bp_spin = QDoubleSpinBox()
        self.gel_bp_spin.setRange(0.1, 1000)
        self.gel_bp_spin.setValue(10)
        self.gel_pct_spin = QDoubleSpinBox()
        self.gel_pct_spin.setRange(0, 100)
        self.gel_pct_spin.setValue(3)
        self.gel_stress_spin = QDoubleSpinBox()
        self.gel_stress_spin.setRange(1, 5)
        self.gel_stress_spin.setSingleStep(0.1)
        self.gel_stress_spin.setValue(1.5)
        controls = [('Пар Primer3 на SNP:', self.primer_pairs_spin),
                    ('Результатов на SNP:', self.top_results_spin),
                    ('Разрешение, п.н.:', self.gel_bp_spin),
                    ('Разрешение, %:', self.gel_pct_spin),
                    ('Ухудшение разрешения, ×:', self.gel_stress_spin)]
        for i, (label, widget) in enumerate(controls):
            row, col = 1 + i // 2, (i % 2) * 2
            design_grid.addWidget(QLabel(label), row, col)
            design_grid.addWidget(widget, row, col + 1)
        note = QLabel('Модель полос без учёта яркости. Различие должно превышать максимум '
                      'из порога в п.н. и процента длины. Параметры требуют проверки на вашем геле.')
        note.setWordWrap(True)
        design_grid.addWidget(note, 4, 0, 1, 4)
        self.min_frag_spin.setToolTip('В совместном поиске — нижняя граница видимых полос. '
                                     'В поиске по одной паре/окну — ограничение на каждый фрагмент.')
        self.max_frag_spin.setToolTip('В совместном поиске — верхняя граница видимых полос. '
                                     'В поиске по одной паре/окну — ограничение на каждый фрагмент.')
        self.use_primer3_chk.toggled.connect(self._update_design_controls)
        self.joint_design_chk.toggled.connect(self._update_design_controls)
        self._update_design_controls()

        btn_row = QHBoxLayout()
        self.run_btn = QPushButton("Запуск анализа")
        self.cancel_btn = QPushButton("Отмена")
        self.cancel_btn.setEnabled(False)
        self.export_csv_btn = QPushButton("Выгрузить в CSV")
        self.export_xlsx_btn = QPushButton("Выгрузить в Excel")
        btn_row.addWidget(self.run_btn)
        btn_row.addWidget(self.cancel_btn)
        btn_row.addWidget(self.export_csv_btn)
        btn_row.addWidget(self.export_xlsx_btn)

        self.progress = QProgressBar()

        # правая часть: лог
        self.status_log = QTextEdit()
        self.status_log.setReadOnly(True)
        self.status_log.setMaximumWidth(380)
        self.status_log.setMinimumWidth(280)
        self.status_log.setMaximumHeight(500)
        self.status_log.setMinimumHeight(120)

        log_norm_layout = QVBoxLayout()

        # ЧЕКБОКСЫ НОРМАЛИЗАЦИИ
        norm_layout = QHBoxLayout()
        self.skip_norm_chk = QCheckBox("Пропускать нормализацию")
        self.skip_norm_chk.setChecked(False)
        self.no_save_norm_chk = QCheckBox("Не сохранять нормализацию")
        self.no_save_norm_chk.setChecked(True)
        norm_layout.addWidget(self.skip_norm_chk)
        norm_layout.addWidget(self.no_save_norm_chk)

        log_norm_layout.addWidget(self.status_log)
        log_norm_layout.addLayout(norm_layout)
        log_norm_layout.addWidget(design_box)

        settings_scroll = QScrollArea()
        settings_scroll.setWidgetResizable(True)
        settings_scroll.setWidget(left_panel)
        left_column = QVBoxLayout()
        left_column.addWidget(settings_scroll)
        left_column.addLayout(btn_row)
        left_column.addWidget(self.progress)
        top_layout.addLayout(left_column, stretch=4)
        top_layout.addLayout(log_norm_layout, stretch=2)

        # ===== низ: фильтры + таблица + визуализация праймеров =====
        bottom_layout = QVBoxLayout()
        self.assay_summary = QLabel('Выберите результат: здесь появятся оценка и полосы трёх генотипов.')
        self.assay_summary.setWordWrap(True)
        self.assay_summary.setTextFormat(Qt.TextFormat.PlainText)
        bottom_layout.addWidget(self.assay_summary)

        filter_layout = QHBoxLayout()
        filter_layout.addWidget(QLabel("Фермент:"))
        self.filter_enzyme_combo = QComboBox()
        self.filter_enzyme_combo.addItem("Все ферменты")
        filter_layout.addWidget(self.filter_enzyme_combo)

        filter_layout.addWidget(QLabel("Паттерн:"))
        self.filter_pattern_combo = QComboBox()
        self.filter_pattern_combo.addItem("Все паттерны")
        filter_layout.addWidget(self.filter_pattern_combo)

        self.open_vis_btn = QPushButton("Открыть визуализацию")
        self.open_vis_btn.setEnabled(False)
        filter_layout.addWidget(self.open_vis_btn)

        filter_layout.addStretch(1)
        bottom_layout.addLayout(filter_layout)

        # Таблица результатов
        self.result_table = QTableWidget()
        self.result_table.setSortingEnabled(True)
        header = self.result_table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        header.setStretchLastSection(True)

        bottom_layout.addWidget(self.result_table)

        main_layout.addLayout(top_layout, stretch=4)
        main_layout.addLayout(bottom_layout, stretch=3)

        central.setLayout(main_layout)
        self.setCentralWidget(central)

        # сигналы
        self.run_btn.clicked.connect(self.run_calc)
        self.cancel_btn.clicked.connect(self.cancel_calc)
        self.export_csv_btn.clicked.connect(self.export_csv)
        self.export_xlsx_btn.clicked.connect(self.export_excel)
        self.filter_enzyme_combo.currentIndexChanged.connect(self.apply_filters)
        self.filter_pattern_combo.currentIndexChanged.connect(self.apply_filters)
        self.result_table.itemSelectionChanged.connect(self.on_table_selection_changed)
        self.open_vis_btn.clicked.connect(self.open_visualization_window)

        self.load_genomes_to_combo()
        self.load_names_to_combo()

        # окно визуализации создаётся по кнопке
        self.vis_dialog = None

        # состояние для визуализации
        self.last_params: dict = {}
        self.run_report: dict = {}
        self.fa_view: Optional[Fasta] = None
        self.col_index: Dict[str, int] = {}
        self._current_vis_payload: Optional[dict] = None
        self._current_vis_error: str = ""
        self.vis_dialog: Optional[RFLPVisualizationDialog] = None
    
    if not _PRIMER3_AVAILABLE:
        logger.warning(f"primer3 недоступен, ошибка импорта: {_PRIMER3_IMPORT_ERROR}")

    def load_genomes_to_combo(self):
        base_dir = get_base_dir()
        folder = os.path.join(base_dir, "genoms")
        self.fasta_combo.clear()
        if not os.path.isdir(folder):
            self.log_msg(f"Папка с геномами не найдена: {folder}", "error")
            return
        files = [f for f in os.listdir(folder) if f.endswith((".fa", ".fasta"))]
        if not files:
            self.log_msg(f"FASTA файлы не найдены в {folder}", "error")
            return
        for f in files:
            self.fasta_combo.addItem(f)
        self.log_msg(f"Найдено {len(files)} FASTA файлов в папке genoms")

    def load_names_to_combo(self):
        base_dir = get_base_dir()
        folder = os.path.join(base_dir, "names")
        self.names_combo.clear()
        if not os.path.isdir(folder):
            self.log_msg(f"Папка с names не найдена: {folder}", "error")
            return
        files = [f for f in os.listdir(folder) if f.endswith(".txt")]
        if not files:
            self.log_msg(f"Файлы names не найдены в {folder}", "error")
        for f in files:
            self.names_combo.addItem(f)
        self.log_msg(f"Загружено {len(files)} файлов names")

    def update_names_selection(self, index):
        fasta_name = self.fasta_combo.currentText()
        base_name = os.path.splitext(fasta_name)[0]
        found = False
        for i in range(self.names_combo.count()):
            names_name = self.names_combo.itemText(i)
            if os.path.splitext(names_name)[0] == base_name:
                self.names_combo.setCurrentIndex(i)
                self.log_msg(
                    f"Выбран файл соответствий: {names_name} для генома {fasta_name}"
                )
                found = True
                break
        if not found and self.names_combo.count() > 0:
            self.log_msg(
                f"Файл соответствий для {fasta_name} не найден среди файлов names",
                "warning",
            )

    def select_file(self, line_edit):
        try:
            path, _ = QFileDialog.getOpenFileName(
                self,
                "Выбрать файл",
                "",
                "Текстовые файлы (*.txt *.tsv *.csv);;Все файлы (*)",
            )
            if path:
                line_edit.setText(path)
                self.log_msg(f"Выбран файл: {path}")
        except Exception:
            tb = traceback.format_exc()
            self.log_msg(f"[ERROR] Ошибка при выборе файла:\n{tb}", "error")

    def _update_design_controls(self):
        self.joint_design_chk.setEnabled(self.use_primer3_chk.isChecked())
        joint = self.use_primer3_chk.isChecked() and self.joint_design_chk.isChecked()
        self.primer_pairs_spin.setEnabled(joint)
        self.top_results_spin.setEnabled(joint)

    def run_calc(self):
        try:
            fasta_selected = self.fasta_combo.currentText()
            names_selected = self.names_combo.currentText()
            if not fasta_selected:
                QMessageBox.warning(self, "Ошибка", "Не выбран FASTA файл")
                return
            base_dir = get_base_dir()
            params = {
                "flank": int(self.flank_spin.value()),
                "min_frag": int(self.min_frag_spin.value()),
                "max_frag": int(self.max_frag_spin.value()),
                "delta": int(self.delta_spin.value()),
                "gain_loss_only": self.gloss_chk.isChecked(),
                "max_cuts": int(self.max_cuts_spin.value()),
                "suppliers": [
                    cb.text().split(":", 1)[1].strip()
                    for cb in self.supplier_checks
                    if cb.isChecked()
                ],
                "fasta": os.path.join(base_dir, "genoms", fasta_selected),
                "input_file": self.input_file.text(),
                "names2_file": (os.path.join(base_dir, "names", names_selected) if names_selected else ""),
                "use_primer3": self.use_primer3_chk.isChecked() and _PRIMER3_AVAILABLE,
                "prod_min": int(self.prod_min_spin.value()),
                "prod_max": int(self.prod_max_spin.value()),
                "tm_min": float(self.tm_min_spin.value()),
                "tm_max": float(self.tm_max_spin.value()),
                "joint_design": self.joint_design_chk.isChecked(),
                "primer_pairs": self.primer_pairs_spin.value(),
                "top_results": self.top_results_spin.value(),
                "gel_resolution_bp": self.gel_bp_spin.value(),
                "gel_resolution_pct": self.gel_pct_spin.value(),
                "gel_stress_factor": self.gel_stress_spin.value(),
                "skip_norm": self.skip_norm_chk.isChecked(),
                "no_save_norm": self.no_save_norm_chk.isChecked(),
            }

            errors = core_impl.validate_params(params)
            if errors:
                QMessageBox.warning(self, "Проверьте параметры", "\n".join(f"• {e}" for e in errors))
                return

            # сохраняем параметры последнего запуска (нужны для визуализации)
            self.last_params = dict(params)
            self.fa_view = None
            self.col_index = {}
            if getattr(self, "vis_dialog", None) is not None:
                self.vis_dialog.panel.show_empty()
            for k in ("flank", "min_frag", "max_frag", "fasta", "input_file"):
                if not params[k]:
                    QMessageBox.warning(self, "Ошибка", f"Не заполнено обязательное поле: {k}")
                    return

            self.status_log.clear()
            self.progress.setRange(0, 0)
            self.run_btn.setDisabled(True)
            self.cancel_btn.setEnabled(True)

            self.worker = RFLPCalcThread(params)
            self.worker.status.connect(self.log_status)
            self.worker.progress.connect(self._on_progress)
            self.worker.finished.connect(self.show_results)
            self.worker.start()
        except Exception:
            tb = traceback.format_exc()
            self.log_msg(f"[ERROR] Exception in run_calc:\n{tb}", "error")

    def cancel_calc(self):
        if getattr(self, "worker", None) is not None and self.worker.isRunning():
            self.worker.requestInterruption()
            self.cancel_btn.setEnabled(False)
            self.log_msg("[INFO] Запрошена отмена анализа…")

    def _on_progress(self, done: int, total: int):
        self.progress.setRange(0, max(1, total))
        self.progress.setValue(done)
        self.progress.setFormat(f"{done} / {total}")

    def log_status(self, msg):
        self.log_msg(msg)

    def log_msg(self, msg, level="info"):
        if level == "info":
            logger.info(msg)
        elif level == "warning":
            logger.warning(msg)
        elif level == "error":
            logger.error(msg)
        elif level == "debug":
            logger.debug(msg)
        self.status_log.append(msg)
        QApplication.processEvents()

    def show_results(self, header, rows):
        try:
            self.assay_summary.setText('Выберите результат: оценка > 1 означает различимость в модели полос.')
            self.progress.setRange(0, 1)
            self.run_btn.setEnabled(True)
            self.cancel_btn.setEnabled(False)
            self.out_header = header
            self.out_rows = rows
            if getattr(self, "worker", None) is not None:
                self.run_report = dict(getattr(self.worker, "report", {}) or {})

            # индекс колонок по имени (нужен для визуализации)
            self.col_index = {str(name): i for i, name in enumerate(header)}

            # подготавливаем FASTA для визуализации (ленивое открытие)
            self.fa_view = None
            fasta_path = (self.last_params or {}).get("fasta")
            if fasta_path:
                try:
                    self.fa_view = Fasta(fasta_path, as_raw=True, sequence_always_upper=True)
                except TypeError:
                    # fallback для старых версий pyfaidx
                    self.fa_view = Fasta(fasta_path, as_raw=True)
                except Exception as e:
                    self.log_msg(f"[WARN] Не удалось открыть FASTA для визуализации: {fasta_path}: {e}", "warning")
                    self.fa_view = None

            if getattr(self, "vis_dialog", None) is not None:
                self.vis_dialog.panel.show_empty()
            self.result_table.clear()
            self.result_table.setRowCount(0)
            self.result_table.setColumnCount(0)
            self.filter_enzyme_combo.clear()
            self.filter_enzyme_combo.addItem("Все ферменты")
            self.filter_pattern_combo.clear()
            self.filter_pattern_combo.addItem("Все паттерны")
            if not rows:
                return

            self.result_table.setColumnCount(len(header))
            self.result_table.setHorizontalHeaderLabels([LABELS.get(name, name) for name in header])
            for c, name in enumerate(header):
                self.result_table.horizontalHeaderItem(c).setToolTip(LABELS.get(name, name))
            table_header = self.result_table.horizontalHeader()
            for position, name in enumerate(('variant', 'enzyme', 'assay_rank',
                                              'genotype_quality', 'worst_margin',
                                              'genotype_margin', 'primer_pair_index')):
                if name in self.col_index:
                    logical = self.col_index[name]
                    table_header.moveSection(table_header.visualIndex(logical), position)
                    self.result_table.setColumnWidth(logical, 150 if position > 2 else 120)
            self.result_table.setRowCount(len(rows))

            self.result_table.setSortingEnabled(False)
            for r, row in enumerate(rows):
                for c, val in enumerate(row):
                    display = f'{val:.4f}'.rstrip('0').rstrip('.') if isinstance(val, float) else str(val)
                    item = NumericTableItem('' if val is None else display)
                    item.setData(Qt.ItemDataRole.UserRole, val)
                    item.setToolTip(f'{LABELS.get(header[c], header[c])}: {val if val is not None else "—"}')
                    self.result_table.setItem(r, c, item)

                pattern_val = str(rows[r][4]) if len(rows[r]) > 4 else ""
                if pattern_val == "no_enzyme_found":
                    for c in range(len(header)):
                        item = self.result_table.item(r, c)
                        if item:
                            item.setForeground(QtGui.QColor("#777777"))

            self.result_table.setSortingEnabled(True)
            self.populate_filters_from_results()
        except Exception:
            tb = traceback.format_exc()
            self.log_msg(f"[ERROR] Exception in show_results:\n{tb}", "error")

    def populate_filters_from_results(self):
        enzymes = set()
        patterns = set()
        for row in getattr(self, "out_rows", []):
            if len(row) > 2 and row[2]:
                enzymes.add(str(row[2]))
            if len(row) > 4 and row[4]:
                patterns.add(str(row[4]))

        current_enzyme = self.filter_enzyme_combo.currentText()
        current_pattern = self.filter_pattern_combo.currentText()

        self.filter_enzyme_combo.blockSignals(True)
        self.filter_pattern_combo.blockSignals(True)

        self.filter_enzyme_combo.clear()
        self.filter_enzyme_combo.addItem("Все ферменты")
        for e in sorted(enzymes):
            self.filter_enzyme_combo.addItem(e)

        self.filter_pattern_combo.clear()
        self.filter_pattern_combo.addItem("Все паттерны")
        for p in sorted(patterns):
            self.filter_pattern_combo.addItem(p)

        idx_e = self.filter_enzyme_combo.findText(current_enzyme)
        if idx_e != -1:
            self.filter_enzyme_combo.setCurrentIndex(idx_e)

        idx_p = self.filter_pattern_combo.findText(current_pattern)
        if idx_p != -1:
            self.filter_pattern_combo.setCurrentIndex(idx_p)

        self.filter_enzyme_combo.blockSignals(False)
        self.filter_pattern_combo.blockSignals(False)

        self.apply_filters()

    def apply_filters(self):
        if not hasattr(self, "out_rows"):
            return
        sel_enzyme = self.filter_enzyme_combo.currentText()
        sel_pattern = self.filter_pattern_combo.currentText()

        # Read values from the visible table order, because sorting changes
        # rows without changing out_rows.
        enzyme_col = self.col_index.get("enzyme", 2)
        pattern_col = self.col_index.get("pattern", 4)
        for r in range(self.result_table.rowCount()):
            enzyme_item = self.result_table.item(r, enzyme_col)
            pattern_item = self.result_table.item(r, pattern_col)
            enzyme = enzyme_item.text() if enzyme_item else ""
            pattern = pattern_item.text() if pattern_item else ""

            if sel_enzyme != "Все ферменты" and enzyme != sel_enzyme:
                self.result_table.setRowHidden(r, True)
                continue
            if sel_pattern != "Все паттерны" and pattern != sel_pattern:
                self.result_table.setRowHidden(r, True)
                continue
            self.result_table.setRowHidden(r, False)

    def on_table_selection_changed(self):
        """Обновляет состояние визуализации при выборе строки (данные готовятся для отдельного окна)."""
        try:
            selected = self.result_table.selectedItems()
            if not selected:
                self.assay_summary.setText('Выберите результат: оценка > 1 означает различимость в модели полос.')
                self._current_vis_payload = None
                self._current_vis_error = ""
                self.open_vis_btn.setEnabled(False)
                if self.vis_dialog is not None and self.vis_dialog.isVisible():
                    self.vis_dialog.panel.show_empty()
                return

            row = selected[0].row()
            self.assay_summary.setText(self._table_cell(row, 'reason'))
            if 'genotype_quality' in self.col_index:
                quality = self._table_cell(row, 'genotype_quality')
                if quality:
                    self.assay_summary.setText(
                        f"{self._table_cell(row, 'reason')}. "
                        f"Оценка: {self._table_cell(row, 'genotype_margin')}; "
                        f"при ухудшении: {self._table_cell(row, 'worst_margin')}.\n"
                        f"Полосы REF/REF: {self._table_cell(row, 'bands_ref_ref') or 'нет'} | "
                        f"REF/ALT: {self._table_cell(row, 'bands_ref_alt') or 'нет'} | "
                        f"ALT/ALT: {self._table_cell(row, 'bands_alt_alt') or 'нет'}")
            status_col = self.col_index.get("status")
            row_status = self.result_table.item(row, status_col).text() if status_col is not None and self.result_table.item(row, status_col) else ""
            if row_status in ("no_primers", "primer3_error", "no_discriminating_assay", "analysis_error"):
                self._current_vis_payload = None
                self._current_vis_error = self._table_cell(row, "reason") or "Визуализация недоступна: праймеры не подобраны."
                self.open_vis_btn.setEnabled(False)
                return
            payload, err = self._build_vis_payload_for_row(row)
            self._current_vis_payload = payload
            self._current_vis_error = err or ""
            self.open_vis_btn.setEnabled(payload is not None)

            # Если окно визуализации уже открыто — обновляем сразу
            if self.vis_dialog is not None and self.vis_dialog.isVisible():
                self.vis_dialog.show_payload(payload, message_if_empty=err or "")
        except Exception:
            tb = traceback.format_exc()
            self.log_msg(f"[ERROR] Ошибка при подготовке визуализации\n{tb}", "error")

    def _table_cell(self, row: int, col_name: str) -> str:
        idx = self.col_index.get(col_name)
        if idx is None:
            return ""
        item = self.result_table.item(row, idx)
        return item.text() if item is not None else ""

    @staticmethod
    def _safe_float(s: object) -> Optional[float]:
        try:
            if s is None:
                return None
            t = str(s).strip()
            if t == "" or t.lower() == "none":
                return None
            return float(t)
        except Exception:
            return None

    def _build_vis_payload_for_row(self, row: int) -> Tuple[Optional[dict], str]:
        """Собирает данные для визуализации (для отдельного окна).

        Возвращает (payload, error_message). Если payload=None — визуализация невозможна.
        """
        fasta_path = (self.last_params or {}).get("fasta")
        if not fasta_path:
            return None, "Не задан FASTA для визуализации. Запустите анализ заново."

        if self.fa_view is None:
            try:
                self.fa_view = Fasta(fasta_path, as_raw=True, sequence_always_upper=True)
            except TypeError:
                self.fa_view = Fasta(fasta_path, as_raw=True)
            except Exception as e:
                return None, f"Не удалось открыть FASTA: {fasta_path}\n{e}"

        mapped_id = self._table_cell(row, "mapped_id").strip()
        amp_start = _safe_int(self._table_cell(row, "amplicon_start"))
        amp_end = _safe_int(self._table_cell(row, "amplicon_end"))
        snp_offset = _safe_int(self._table_cell(row, "snp_offset_in_amplicon"))

        primer_left = self._table_cell(row, "primers_left").strip()
        primer_right = self._table_cell(row, "primers_right").strip()
        tm_left = self._safe_float(self._table_cell(row, "tm_left"))
        tm_right = self._safe_float(self._table_cell(row, "tm_right"))

        p_left_start = _safe_int(self._table_cell(row, "primer_left_start"))
        p_left_len = _safe_int(self._table_cell(row, "primer_left_len"))
        p_right_start = _safe_int(self._table_cell(row, "primer_right_start"))
        p_right_len = _safe_int(self._table_cell(row, "primer_right_len"))

        if not mapped_id or amp_start is None or amp_end is None:
            return None, "В выбранной строке отсутствуют данные для визуализации (mapped_id / amplicon_start / amplicon_end)."

        if amp_start < 1 or amp_end < amp_start:
            return None, f"Некорректные координаты ампликона: {mapped_id}:{amp_start}-{amp_end}"

        try:
            # pyfaidx: end в срезе является end-exclusive, поэтому подаём amp_end как inclusive
            amplicon_seq = str(self.fa_view[mapped_id][amp_start - 1 : amp_end])
        except Exception as e:
            return None, f"Не удалось получить последовательность из FASTA: {mapped_id}:{amp_start}-{amp_end}\n{e}"

        variant_str = self._table_cell(row, "variant").strip()
        enzyme_name = self._table_cell(row, "enzyme").strip()
        pattern = self._table_cell(row, "pattern").strip()
        frags_ref = self._table_cell(row, "frags_ref").strip()
        frags_alt = self._table_cell(row, "frags_alt").strip()

        payload = dict(
            mapped_id=mapped_id,
            genome_start=int(amp_start),
            genome_end=int(amp_end),
            amplicon_seq=amplicon_seq,
            snp_offset=int(snp_offset) if snp_offset is not None else None,
            primer_left=primer_left,
            primer_right=primer_right,
            tm_left=tm_left,
            tm_right=tm_right,
            primer_left_start=p_left_start,
            primer_left_len=p_left_len,
            primer_right_start=p_right_start,
            primer_right_len=p_right_len,
            variant_str=variant_str,
            enzyme_name=enzyme_name,
            pattern=pattern,
            frags_ref=frags_ref,
            frags_alt=frags_alt,
        )

        return payload, ""

    def open_visualization_window(self):
        """Открывает отдельное окно визуализации (адаптивное)."""
        if self.vis_dialog is None:
            self.vis_dialog = RFLPVisualizationDialog(self)

        # Показать окно и поднять наверх
        self.vis_dialog.show()
        self.vis_dialog.raise_()
        self.vis_dialog.activateWindow()

        if self._current_vis_payload is not None:
            self.vis_dialog.show_payload(self._current_vis_payload)
        else:
            msg = self._current_vis_error or "Выберите строку в таблице результатов."
            self.vis_dialog.show_payload(None, message_if_empty=msg)


    def _ensure_results_dir(self) -> str:
        base_dir = get_base_dir()
        results_dir = os.path.join(base_dir, "results")
        if not os.path.isdir(results_dir):
            os.makedirs(results_dir, exist_ok=True)
        return results_dir

    def _make_timestamp_name(self, ext: str) -> str:
        ts = datetime.now().strftime("%Y%m%d%H%M%S")
        return f"rflp_candidates_{ts}.{ext}"

    def _collect_table_data(self) -> Tuple[list, list]:
        rows = self.result_table.rowCount()
        cols = self.result_table.columnCount()
        if rows == 0 or cols == 0:
            return [], []

        header = list(self.out_header)
        data = []
        any_visible = False
        for r in range(rows):
            if self.result_table.isRowHidden(r):
                continue
            any_visible = True
            data.append([
                self.result_table.item(r, c).data(Qt.ItemDataRole.UserRole)
                if self.result_table.item(r, c) is not None else None
                for c in range(cols)
            ])

        if not any_visible:
            return [], []
        return header, data

    def export_csv(self):
        try:
            header, data = self._collect_table_data()
            if not data:
                QMessageBox.information(self, "Информация", "Нет видимых строк для экспорта")
                return
            results_dir = self._ensure_results_dir()
            filename = self._make_timestamp_name("csv")
            path = os.path.join(results_dir, filename)
            safe_export_csv(path, header, data)
            self.log_msg(f"[EXPORT] CSV сохранён: {path}")
            QMessageBox.information(self, "Экспорт CSV", f"Результаты сохранены:\n{path}")
        except Exception:
            tb = traceback.format_exc()
            self.log_msg(f"[ERROR] Exception in export_csv:\n{tb}", "error")

    def export_excel(self):
        try:
            header, data = self._collect_table_data()
            if not data:
                QMessageBox.information(self, "Информация", "Нет видимых строк для экспорта")
                return
            results_dir = self._ensure_results_dir()
            filename = self._make_timestamp_name("xlsx")
            path = os.path.join(results_dir, filename)

            safe_export_excel(path, header, data, metadata={
                "app_version": "1.1.0",
                "scoring_model": core_impl.MODEL_VERSION,
                "run_timestamp": datetime.now().isoformat(timespec="seconds"),
                "genome_path": (self.last_params or {}).get("fasta", ""),
                "parameters": self.last_params,
                "status_summary": (self.run_report or {}).get("summary", {"rows": len(data)}),
            }, rejected=(self.run_report or {}).get("rejected", []))
            self.log_msg(f"[EXPORT] Excel сохранён: {path}")
            QMessageBox.information(self, "Экспорт Excel", f"Результаты сохранены:\n{path}")
        except Exception:
            tb = traceback.format_exc()
            self.log_msg(f"[ERROR] Exception in export_excel:\n{tb}", "error")

# Use the standalone, exportable renderer in the application. The legacy
# implementation above remains available for compatibility with old imports.
RFLPVisualizationDialog = enhanced_visualization.RFLPVisualizationDialog


def configure_system_theme(app: QApplication):
    """Use the Windows/system color scheme and refresh when it changes.

    The previous stylesheet forced white widgets, which made all controls
    unreadable when Windows was using a dark theme.  We set semantic palette
    roles once and generate the stylesheet from the detected palette instead.
    """
    state = {"updating": False}

    def refresh(*_):
        if state["updating"]:
            return
        state["updating"] = True
        try:
            hints = app.styleHints()
            scheme = getattr(hints, "colorScheme", None)
            scheme_value = scheme() if callable(scheme) else None
            scheme_name = str(scheme_value).lower()
            palette = app.palette()
            dark = "dark" in scheme_name or palette.color(QtGui.QPalette.ColorRole.Window).lightness() < 128

            if dark:
                colors = dict(window="#20252b", base="#171b20", alternate="#20262d", text="#f1f5f9",
                               muted="#a9b7c5", border="#4a5968", button="#29323b", hover="#344452",
                               selected="#315c70", selected_text="#ffffff", header="#2a3540",
                               disabled="#71808e", progress="#35b7a6")
            else:
                colors = dict(window="#f4f7fb", base="#ffffff", alternate="#f3f7fa", text="#22354b",
                               muted="#61758b", border="#c6d4df", button="#ffffff", hover="#e7f2f4",
                               selected="#cdebe7", selected_text="#123c40", header="#e8f0f5",
                               disabled="#8998a5", progress="#2c9c91")

            pal = QtGui.QPalette(palette)
            roles = {
                QtGui.QPalette.ColorRole.Window: colors["window"],
                QtGui.QPalette.ColorRole.Base: colors["base"],
                QtGui.QPalette.ColorRole.AlternateBase: colors["alternate"],
                QtGui.QPalette.ColorRole.Text: colors["text"],
                QtGui.QPalette.ColorRole.WindowText: colors["text"],
                QtGui.QPalette.ColorRole.Button: colors["button"],
                QtGui.QPalette.ColorRole.ButtonText: colors["text"],
                QtGui.QPalette.ColorRole.Highlight: colors["selected"],
                QtGui.QPalette.ColorRole.HighlightedText: colors["selected_text"],
                QtGui.QPalette.ColorRole.PlaceholderText: colors["muted"],
                QtGui.QPalette.ColorRole.ToolTipBase: colors["header"],
                QtGui.QPalette.ColorRole.ToolTipText: colors["text"],
            }
            for role, value in roles.items():
                pal.setColor(role, QtGui.QColor(value))
            pal.setColor(QtGui.QPalette.ColorGroup.Disabled, QtGui.QPalette.ColorRole.Text, QtGui.QColor(colors["disabled"]))
            pal.setColor(QtGui.QPalette.ColorGroup.Disabled, QtGui.QPalette.ColorRole.ButtonText, QtGui.QColor(colors["disabled"]))
            app.setPalette(pal)
            app.setStyleSheet(f"""
                QMainWindow, QDialog {{ background: {colors['window']}; color: {colors['text']}; }}
                QWidget {{ color: {colors['text']}; }}
                QGroupBox {{ font-weight: 600; border: 1px solid {colors['border']}; border-radius: 7px; margin-top: 10px; padding: 8px; }}
                QGroupBox::title {{ subcontrol-origin: margin; left: 10px; padding: 0 4px; color: {colors['text']}; }}
                QPushButton {{ padding: 6px 10px; border: 1px solid {colors['border']}; border-radius: 5px; background: {colors['button']}; color: {colors['text']}; }}
                QPushButton:hover {{ background: {colors['hover']}; border-color: {colors['progress']}; }}
                QPushButton:disabled {{ color: {colors['disabled']}; background: {colors['window']}; }}
                QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {{ padding: 5px; border: 1px solid {colors['border']}; border-radius: 4px; background: {colors['base']}; color: {colors['text']}; selection-background-color: {colors['selected']}; selection-color: {colors['selected_text']}; }}
                QTableWidget {{ background: {colors['base']}; alternate-background-color: {colors['alternate']}; color: {colors['text']}; gridline-color: {colors['border']}; selection-background-color: {colors['selected']}; selection-color: {colors['selected_text']}; }}
                QHeaderView::section {{ background: {colors['header']}; color: {colors['text']}; padding: 6px; border: none; border-bottom: 1px solid {colors['border']}; font-weight: 600; }}
                QTextEdit {{ background: {colors['base']}; color: {colors['text']}; border: 1px solid {colors['border']}; border-radius: 5px; }}
                QProgressBar {{ border: none; border-radius: 4px; background: {colors['header']}; color: {colors['text']}; text-align: center; }}
                QProgressBar::chunk {{ background: {colors['progress']}; border-radius: 4px; }}
                QToolTip {{ background: {colors['header']}; color: {colors['text']}; border: 1px solid {colors['border']}; }}
            """)
        finally:
            state["updating"] = False

    app.setStyle("Fusion")
    refresh()
    # Windows emits paletteChanged when the user changes the system theme.
    app.paletteChanged.connect(refresh)


def main():
    app = QApplication(sys.argv)
    configure_system_theme(app)
    win = RFLPPickerGUI()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
