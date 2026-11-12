"""Laboratory sheets: one A4 page per SNP with the chosen assay, for printing (PDF).

The page holds what goes to the bench — enzyme, primers, expected bands, the
virtual gel and the amplicon map — and fields to fill in by hand. Reaction
conditions are not printed: the program does not know them.
"""
import os
from datetime import datetime
from typing import List, Tuple

from PyQt6.QtCore import QMarginsF, QRectF, Qt
from PyQt6.QtGui import QColor, QFont, QPageLayout, QPageSize, QPainter, QPdfWriter, QPen

import palette
import rflp_core as core
from i18n import fmt, t
from reporting import _atomic_destination
from snp_summary import _DESIGN_TITLES, SnpSummary, status_title
from visualization import GelPanel, PrimerGraphicPanel, payload_for_scheme, render_payload

DPI = 150
MARGIN_MM = 14
INK, MUTED, LINE = '#1F1F1F', '#616161', '#B0B0B0'


def sheet_fields(summary: SnpSummary, dcaps_mode: str) -> List[Tuple[str, str]]:
    """The (title, value) lines of the assay of one SNP; empty values are left out."""
    scheme = summary.scheme or {}
    design = scheme.get('design')
    mismatches = scheme.get('primer_mismatches')
    rows = [
        (t('Статус'), status_title(summary.status)),
        (t('Фермент'), scheme.get('enzyme')),
        (t('Сайт'), scheme.get('site')),
        (t('Изошизомеры'), scheme.get('isoschizomers')),
        (t('Поставщики'), scheme.get('suppliers')),
        (t('Тип схемы'), t(_DESIGN_TITLES[design]) if design in _DESIGN_TITLES else design),
        (t('Несовпадения в праймере'), mismatches),
        (t('Продукт, п.н.'), scheme.get('product_size')),
        (t('Полосы REF/REF'), scheme.get('bands_ref_ref')),
        (t('Полосы REF/ALT'), scheme.get('bands_ref_alt')),
        (t('Полосы ALT/ALT'), scheme.get('bands_alt_alt')),
        (t('Примечание'), summary.note(dcaps_mode)),
    ]
    return [(title, str(value)) for title, value in rows if value not in (None, '')]


def primer_rows(summary: SnpSummary) -> List[Tuple[str, str, str, str]]:
    """(name, sequence, length, Tm) of the two primers."""
    scheme = summary.scheme or {}
    rows = []
    for side, primer, tm in (('F', 'primers_left', 'tm_left'), ('R', 'primers_right', 'tm_right')):
        sequence = scheme.get(primer) or ''
        if sequence:
            rows.append((f'{summary.primer_base}_{side}', sequence, str(len(sequence)), fmt(scheme.get(tm))))
    return rows


def handwritten_fields() -> List[str]:
    return [t('Дата'), t('Оператор'), t('№ геля')]


class _Page:
    """Text and pictures laid out from the top of a page; coordinates in device pixels."""

    def __init__(self, painter: QPainter, width: float, height: float):
        self.p, self.width, self.height, self.y = painter, width, height, 0.0

    def font(self, points: float, bold: bool = False, mono: bool = False) -> QFont:
        font = QFont('DejaVu Sans Mono' if mono else 'DejaVu Sans')
        if mono:
            font.setStyleHint(QFont.StyleHint.Monospace)
        else:
            font.setStyleHint(QFont.StyleHint.SansSerif)
        font.setPointSizeF(points)
        font.setBold(bold)
        return font

    def text(self, text: str, points: float = 9.5, bold: bool = False, mono: bool = False, color: str = INK,
             x: float = 0.0, width: float = None, advance: bool = True) -> float:
        """Draw wrapped text at the current height; return its height."""
        self.p.setFont(self.font(points, bold, mono))
        self.p.setPen(QColor(color))
        rect = QRectF(x, self.y, (self.width - x) if width is None else width, self.height)
        flags = Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop | Qt.TextFlag.TextWordWrap
        if mono:                      # a sequence has no spaces to wrap at
            flags |= Qt.TextFlag.TextWrapAnywhere
        flags = int(flags.value) if hasattr(flags, 'value') else int(flags)
        used = self.p.boundingRect(rect, flags, text).height()
        self.p.drawText(rect, flags, text)
        if advance:
            self.y += used
        return used

    def gap(self, pixels: float) -> None:
        self.y += pixels

    def rule(self) -> None:
        self.p.setPen(QPen(QColor(LINE), 1))
        self.p.drawLine(0, round(self.y), round(self.width), round(self.y))
        self.y += 6

    def scene(self, scene, x: float, y: float, width: float, height: float, right: bool = False) -> float:
        """A scene scaled into the box, keeping its proportions, in a thin frame; return the drawn height."""
        source = scene.sceneRect()
        if source.width() <= 0 or source.height() <= 0:
            return 0.0
        scale = min(width / source.width(), height / source.height())
        drawn = source.width() * scale
        target = QRectF(x + (width - drawn if right else 0), y, drawn, source.height() * scale)
        scene.render(self.p, target, source)
        self.p.setPen(QPen(QColor(LINE), 1))
        self.p.setBrush(Qt.BrushStyle.NoBrush)
        self.p.drawRect(target)
        return target.height()


def draw_sheet(page: _Page, summary: SnpSummary, payload, message: str, dcaps_mode: str, footer: str,
               panels) -> None:
    map_panel, gel_panel = panels
    mm = DPI / 25.4
    snp_id = summary.record.snp_id

    page.text(t('RFLP Picker · лабораторный лист'), 8, color=MUTED)
    page.gap(2)
    page.text(f'{snp_id} · {summary.label}' if snp_id else summary.label, 16, bold=True)
    page.gap(6)
    # To fill in by hand.
    column = page.width / len(handwritten_fields())
    for index, title in enumerate(handwritten_fields()):
        page.text(f'{title}: ' + '_' * 18, 9.5, x=index * column, width=column, advance=False)
    page.gap(26)
    page.rule()
    page.gap(4)

    # Left: the assay; right: the virtual gel.
    render_payload(map_panel, gel_panel, payload, message)
    top = page.y
    gel_width = 80 * mm
    left_width = page.width - gel_width - 6 * mm
    title_width = 52 * mm
    for title, value in sheet_fields(summary, dcaps_mode):
        used_title = page.text(title, 9, color=MUTED, width=title_width, advance=False)
        used = page.text(value, 9.5, bold=title == t('Фермент'), x=title_width, width=left_width - title_width,
                         advance=False)
        page.gap(max(used, used_title) + 4)
    left_bottom = page.y
    gel_height = page.scene(gel_panel.scene, page.width - gel_width, top, gel_width, 70 * mm, right=True) \
        if payload else 0
    page.y = max(left_bottom, top + gel_height) + 8
    page.rule()

    page.text(t('Праймеры'), 10, bold=True)
    page.gap(3)
    columns = (0, 34 * mm, page.width - 36 * mm, page.width - 18 * mm)
    for index, title in enumerate((t('Имя'), t('Последовательность 5′→3′'), t('Длина'), 'Tm, °C')):
        page.text(title, 8, color=MUTED, x=columns[index], advance=False)
    page.gap(19)
    for name, sequence, length, tm in primer_rows(summary):
        page.text(name, 9.5, x=columns[0], width=columns[1] - 4, advance=False)
        page.text(length, 9.5, x=columns[2], advance=False)
        page.text(tm, 9.5, x=columns[3], advance=False)
        page.gap(page.text(sequence, 9.5, mono=True, x=columns[1], width=columns[2] - columns[1] - 6,
                           advance=False) + 4)
    page.gap(4)
    page.rule()

    page.text(t('Карта ампликона'), 10, bold=True)
    page.gap(3)
    if payload:
        page.gap(page.scene(map_panel.scene, 0, page.y, page.width, 95 * mm) + 8)
    else:
        page.text(message, 9, color=MUTED)
        page.gap(8)

    # The rest of the page: a frame for the photo of the gel.
    footer_height = 30
    box_bottom = page.height - footer_height - 8
    if box_bottom - page.y > 30 * mm:
        page.text(t('Фото геля'), 10, bold=True)
        page.gap(4)
        page.p.setPen(QPen(QColor(LINE), 1, Qt.PenStyle.DashLine))
        page.p.drawRect(QRectF(0, page.y, page.width, box_bottom - page.y))
    page.y = page.height - footer_height
    page.rule()
    page.text(footer, 7.5, color=MUTED)


def export_lab_sheets(path: str, summaries: List[SnpSummary], fasta_path: str, params: dict, dcaps_mode: str) -> int:
    """Write one page per SNP that has an assay; return the number of pages.

    Without such SNPs no file is written. The destination is replaced only after
    a complete write. A missing genome leaves the pictures out, not the sheet.
    """
    printable = [summary for summary in summaries if summary.scheme is not None]
    if not printable:
        return 0
    fasta = None
    if fasta_path and os.path.isfile(fasta_path):
        from pyfaidx import Fasta
        fasta = Fasta(fasta_path, as_raw=True, sequence_always_upper=True)
    footer = t('RFLP Picker {version} · геном {genome} · продукт {prod_min}–{prod_max} п.н., Tm {tm_min}–{tm_max} °C · '
               'лист создан {date}', version=core.APP_VERSION, genome=os.path.basename(fasta_path or '—'),
               prod_min=params.get('prod_min'), prod_max=params.get('prod_max'), tm_min=fmt(params.get('tm_min')),
               tm_max=fmt(params.get('tm_max')), date=datetime.now().strftime('%Y-%m-%d'))
    # Drawn on paper: the light palette whatever the theme of the window.
    panels = PrimerGraphicPanel(compact=True), GelPanel(compact=True)
    panels[0].apply_colors(palette.colors(False))
    try:
        with _atomic_destination(path) as temporary:
            writer = QPdfWriter(str(temporary))
            writer.setResolution(DPI)
            writer.setTitle(t('RFLP Picker · лабораторный лист'))
            writer.setCreator(f'RFLP Picker {core.APP_VERSION}')
            writer.setPageLayout(QPageLayout(QPageSize(QPageSize.PageSizeId.A4), QPageLayout.Orientation.Portrait,
                                             QMarginsF(MARGIN_MM, MARGIN_MM, MARGIN_MM, MARGIN_MM),
                                             QPageLayout.Unit.Millimeter))
            painter = QPainter(writer)
            try:
                painter.setRenderHints(QPainter.RenderHint.Antialiasing | QPainter.RenderHint.TextAntialiasing)
                area = writer.pageLayout().paintRectPixels(DPI)
                for index, summary in enumerate(printable):
                    if index:
                        writer.newPage()
                    payload, message = payload_for_scheme(summary.scheme, fasta_path, params, fasta)
                    draw_sheet(_Page(painter, area.width(), area.height()), summary, payload, message, dcaps_mode,
                               footer, panels)
            finally:
                painter.end()
    finally:
        if fasta is not None:
            fasta.close()
        for panel in panels:
            panel.deleteLater()
    return len(printable)
