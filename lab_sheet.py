"""Laboratory sheets: one A4 page per SNP with the chosen assay, for printing (PDF).

The page holds what goes to the bench — enzyme, primers, expected bands, the
virtual gel and the amplicon map — and fields to fill in by hand and source-attributed product properties.
"""
import os
from datetime import datetime
from typing import List, Tuple

from PyQt6.QtCore import QMarginsF, QRectF, Qt
from PyQt6.QtGui import QColor, QFont, QPageLayout, QPageSize, QPainter, QPdfWriter, QPen

from allele_specific import reaction_table
from aspcr_display import genotype_bands, genotype_headings, genotype_rows, tube_bands, tube_heading, unresolved_text
from enzyme_catalog import product_sections
import palette
import rflp_core as core
import snp_summary
from i18n import fmt, t
from reporting import _atomic_destination
from snp_summary import _DESIGN_TITLES, SnpSummary, allele_products, locus_bands, status_title
from visualization import GelPanel, PrimerGraphicPanel, gel_view_options, payload_for_scheme, render_payload

DPI = 150
MARGIN_MM = 14
INK, MUTED, LINE = '#1F1F1F', '#616161', '#B0B0B0'
MAX_VALUE = 220   # characters of one value: a longer list would push the primers and the map off the page


def sheet_fields(summary: SnpSummary, dcaps_mode: str) -> List[Tuple[str, str]]:
    """The (title, value) lines of the assay of one SNP; empty values are left out."""
    scheme = summary.scheme or {}
    design = scheme.get('design')
    mismatches = scheme.get('primer_mismatches')
    suppliers = scheme.get('suppliers')
    if scheme.get('enzyme_2') and suppliers:
        # The narrow printed column should leave room for suppliers of both enzymes.
        suppliers = '\n'.join(part if len(part) <= 60 else part[:59].rstrip() + '…'
                              for part in suppliers.split('; '))
    rows = [
        (t('Статус'), status_title(summary.status)),
        (t('Фермент'), scheme.get('enzyme')),
        (t('Сайт'), scheme.get('site')),
        (t('Изошизомеры'), scheme.get('isoschizomers')),
        (t('Поставщики'), suppliers),
        (t('Тип схемы'), t(_DESIGN_TITLES[design]) if design in _DESIGN_TITLES else design),
        (t('Несовпадения в праймере'), mismatches),
        (t('Продукт, п.н.'), scheme.get('product_size')),
        (t('Продукт ALT, п.н.'), scheme.get('product_size_alt')),
        (t('Полосы REF/REF'), scheme.get('bands_ref_ref')),
        (t('Полосы REF/ALT'), scheme.get('bands_ref_alt')),
        (t('Полосы ALT/ALT'), scheme.get('bands_alt_alt')),
        (t('Примечание'), summary.note(dcaps_mode)),
    ]
    texts = [(title, str(value)) for title, value in rows if value not in (None, '')]
    return [(title, value if len(value) <= MAX_VALUE else value[:MAX_VALUE - 1].rstrip() + '…')
            for title, value in texts]


def primer_rows(summary: SnpSummary) -> List[Tuple[str, str, str, str]]:
    """(name, sequence, length, Tm) of the primers."""
    scheme = summary.scheme or {}
    if scheme.get('design') == 'allele_specific':
        return [(name, sequence, str(len(sequence)), fmt(tm))
                for name, sequence, tm in snp_summary.aspcr_primer_rows(summary)]
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
             x: float = 0.0, width: float = None, advance: bool = True, draw: bool = True) -> float:
        """Draw wrapped text at the current height; return its height."""
        self.p.setFont(self.font(points, bold, mono))
        self.p.setPen(QColor(color))
        rect = QRectF(x, self.y, (self.width - x) if width is None else width, self.height)
        flags = Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop | Qt.TextFlag.TextWordWrap
        if mono:                      # a sequence has no spaces to wrap at
            flags |= Qt.TextFlag.TextWrapAnywhere
        flags = int(flags.value) if hasattr(flags, 'value') else int(flags)
        used = self.p.boundingRect(rect, flags, text).height()
        if draw:
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
    label = (f'{summary.record.chrom}:{summary.record.pos}' if summary.scheme.get('genotype_bands')
             else summary.label)
    page.text(f'{snp_id} · {label}' if snp_id else label, 16, bold=True)
    page.gap(6)
    # To fill in by hand.
    column = page.width / len(handwritten_fields())
    for index, title in enumerate(handwritten_fields()):
        page.text(f'{title}: ' + '_' * 18, 9.5, x=index * column, width=column, advance=False)
    page.gap(26)
    page.rule()
    page.gap(4)

    if summary.scheme.get('design') == 'allele_specific':
        _draw_aspcr(page, summary, dcaps_mode, gel_panel, payload or {})
        _finish_page(page, footer)
        return

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

    sections = product_sections(summary.scheme, map_panel.allele_sequences() if payload else None)
    if sections:
        top = bottom = page.y
        column = page.width / len(sections)
        for index, lines in enumerate(sections):
            page.y = top
            for line in lines:
                page.text(line, 8.5, mono=line.startswith('https://'), x=index * column,
                          width=column - (8 if len(sections) == 2 else 0))
            bottom = max(bottom, page.y)
        page.y = bottom
        page.gap(4)
        page.rule()

    _draw_primers(page, summary)

    footer_height = 48 if '\n' in footer else 30
    page.text(t('Карта ампликона'), 10, bold=True)
    page.gap(3)
    if payload:
        # Reserve the footer, photo caption and a writable photo frame on crowded sheets.
        map_height = min(95 * mm, page.height - page.y - 40 * mm - footer_height - 8)
        page.gap(page.scene(map_panel.scene, 0, page.y, page.width, max(1, map_height)) + 8)
    else:
        page.text(message, 9, color=MUTED)
        page.gap(8)

    _finish_page(page, footer)


def _draw_primers(page: _Page, summary: SnpSummary) -> None:
    mm = DPI / 25.4
    primers = primer_rows(summary)
    if primers:                                   # a window-search assay has none
        page.text(t('Праймеры'), 10, bold=True)
        page.gap(3)
        columns = (0, 34 * mm, page.width - 36 * mm, page.width - 18 * mm)
        for index, title in enumerate((t('Имя'), t('Последовательность 5′→3′'), t('Длина'), 'Tm, °C')):
            page.text(title, 8, color=MUTED, x=columns[index], advance=False)
        page.gap(19)
        for name, sequence, length, tm in primers:
            page.text(name, 9.5, x=columns[0], width=columns[1] - 4, advance=False)
            page.text(length, 9.5, x=columns[2], advance=False)
            page.text(tm, 9.5, x=columns[3], advance=False)
            page.gap(page.text(sequence, 9.5, mono=True, x=columns[1], width=columns[2] - columns[1] - 6,
                               advance=False) + 4)
        page.gap(4)
        page.rule()


def _draw_aspcr(page: _Page, summary: SnpSummary, dcaps_mode: str, gel_panel: GelPanel, options: dict) -> None:
    scheme = summary.scheme
    for title, value in sheet_fields(summary, dcaps_mode):
        page.text(f'{title} · {value}', 9.5)
    page.gap(8)
    _draw_primers(page, summary)
    tubes = reaction_table(scheme)
    if not tubes:
        return
    names = {sequence: name for name, sequence, _, _ in primer_rows(summary)}
    for tube in tubes:
        page.text(tube_heading(tube), 10, bold=True)
        page.text(tube_bands(tube), 9.5)
        page.text(' · '.join(names[primer['sequence']] for primer in tube['primers']), 9.5, mono=True)
        page.gap(6)
    column = page.width / 3
    for index, row in enumerate([genotype_headings(), *genotype_rows(tubes)]):
        height = max(page.text(value, 9.5, bold=index == 0, x=i * column, width=column - 6, advance=False)
                     for i, value in enumerate(row))
        page.gap(height + 3)
    page.gap(8)
    page.text(unresolved_text(), 9.5, bold=True)
    page.gap(8)
    gel_panel.render_aspcr(genotype_bands(tubes), scheme.get('variant', ''), **options)
    page.gap(page.scene(gel_panel.scene, 0, page.y, page.width, 60 * DPI / 25.4))
    page.gap(8)


def _finish_page(page: _Page, footer: str) -> None:
    mm = DPI / 25.4
    footer_height = 48 if '\n' in footer else 30
    # The rest of the page: a frame for the photo of the gel.
    box_bottom = page.height - footer_height - 8
    if box_bottom - page.y > 30 * mm:
        page.text(t('Фото геля'), 10, bold=True)
        page.gap(4)
        page.p.setPen(QPen(QColor(LINE), 1, Qt.PenStyle.DashLine))
        page.p.drawRect(QRectF(0, page.y, page.width, box_bottom - page.y))
    page.y = page.height - footer_height
    page.rule()
    page.text(footer, 7.5, color=MUTED)


def footer_text(fasta_path: str, params: dict, metadata: dict | None = None, *, show_digest=True) -> str:
    """The program version, the genome and the search settings that applied."""
    metadata = metadata or {}
    genome = metadata.get('assembly_accession') or metadata.get('genome_name') or os.path.basename(fasta_path or '—')
    parts = [f"RFLP Picker {metadata.get('app_version', core.APP_VERSION)}", t('геном {genome}', genome=genome)]
    if params.get('use_primer3', True):           # the window search has no product and no Tm
        parts.append(t('продукт {prod_min}–{prod_max} п.н., Tm {tm_min}–{tm_max} °C',
                       prod_min=params.get('prod_min'), prod_max=params.get('prod_max'),
                       tm_min=fmt(params.get('tm_min')), tm_max=fmt(params.get('tm_max'))))
    parts.append(t('лист создан {date}', date=datetime.now().strftime('%Y-%m-%d')))
    text = ' · '.join(parts)
    model = metadata.get('scoring_model') or params.get('gel_model_version')
    if model and 'gel_agarose_pct' in params:
        text += '\n' + t('Гель {model} · агароза {agarose}% · маркер {marker}',
                         model=model, agarose=fmt(params['gel_agarose_pct']), marker=params.get('gel_marker', '—'))
        if show_digest:
            text += t(' · вероятность разреза {digest}%', digest=fmt(params.get('gel_digest_pct', 100)))
    return text


def _draw_locus_tables(writer, painter, area, summary):
    scheme = summary.scheme
    worst = locus_bands(scheme, stressed=True)
    products = allele_products(scheme)
    product_rows = []
    for product in products:
        sequence = product['sequence']
        product_rows.append((product['allele'], sequence[:240], str(product['length']),
                             ' + '.join(map(str, product['fragments']))))
        product_rows.extend(('', sequence[start:start + 240], '', '') for start in range(240, len(sequence), 240))
    genotype_rows = [(name, bands, worst[name]) for name, bands in locus_bands(scheme).items()]
    locus = f'{summary.record.chrom}:{summary.record.pos}'
    name = f'{summary.record.snp_id} · {locus}' if summary.record.snp_id else locus
    page = None
    pages = 0
    for headings, rows, fractions, sequence_column in (
            ((t('Аллель'), t('Последовательность аллеля'), t('Продукт, п.н.'), t('Фрагменты, п.н.')),
             product_rows, (0.13, 0.42, 0.17, 0.28), 1),
            ((t('Генотип'), t('Ожидаемые полосы, п.н.'), t('Полосы при ухудшении, п.н.')),
             genotype_rows, (0.19, 0.40, 0.41), None)):
        widths = [area.width() * fraction for fraction in fractions]
        x_positions = [sum(widths[:i]) for i in range(len(widths))]

        def table_row(values, bold=False, draw=True):
            return max(page.text(value, 9.5, bold=bold, mono=(i == sequence_column and not bold),
                                 x=x, width=width - 7, advance=False, draw=draw)
                       for i, (value, x, width) in enumerate(zip(values, x_positions, widths))) + 6

        def start_page():
            if not writer.newPage():
                raise OSError(t('Не удалось создать страницу PDF'))
            new_page = _Page(painter, area.width(), area.height())
            new_page.text(name, 16, bold=True)
            new_page.gap(12)
            return new_page

        if page is None or page.y > area.height() - 150:
            page = start_page()
            pages += 1
        page.gap(table_row(headings, bold=True))
        page.rule()
        for values in rows:
            if page.y + table_row(values, draw=False) > area.height() - 24:
                page = start_page()
                pages += 1
                page.gap(table_row(headings, bold=True))
                page.rule()
            page.gap(table_row(values))
        page.gap(20)
    return pages


def export_lab_sheets(path: str, summaries: List[SnpSummary], fasta_path: str, params: dict, dcaps_mode: str,
                      *, payload_lookup=None, metadata=None) -> int:
    """Write sheets for loci with an assay; return the number of pages.

    Without such SNPs no file is written. The destination is replaced only after
    a complete write. A missing genome leaves the pictures out, not the sheet.
    """
    printable = [summary for summary in summaries if summary.scheme is not None]
    if not printable:
        return 0
    fasta = None
    if (payload_lookup is None and fasta_path and os.path.isfile(fasta_path) and
            any(summary.scheme.get('design') != 'allele_specific' for summary in printable)):
        from pyfaidx import Fasta
        fasta = Fasta(fasta_path, as_raw=True, sequence_always_upper=True)
    footer = footer_text(fasta_path, params, metadata)
    pcr_footer = footer_text(fasta_path, params, metadata, show_digest=False)
    # Drawn on paper: the light palette whatever the theme of the window.
    panels = PrimerGraphicPanel(compact=True), GelPanel(compact=True)
    panels[0].apply_colors(palette.colors(False))
    pages = len(printable)
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
            if not painter.isActive():
                raise OSError(t('Не удалось создать PDF: {path}', path=path))
            try:
                painter.setRenderHints(QPainter.RenderHint.Antialiasing | QPainter.RenderHint.TextAntialiasing)
                area = writer.pageLayout().paintRectPixels(DPI)
                for index, summary in enumerate(printable):
                    if index:
                        writer.newPage()
                    if summary.scheme.get('design') == 'allele_specific':
                        payload, message = gel_view_options(summary.scheme, params), ''
                        payload.pop('digest_pct', None)
                    elif payload_lookup is None:
                        payload, message = payload_for_scheme(summary.scheme, fasta_path, params, fasta)
                    else:
                        payload, message = payload_lookup(summary.scheme)
                    draw_sheet(_Page(painter, area.width(), area.height()), summary, payload, message, dcaps_mode,
                               pcr_footer if summary.scheme.get('design') == 'allele_specific' else footer, panels)
                    if summary.scheme.get('genotype_bands'):
                        pages += _draw_locus_tables(writer, painter, area, summary)
            finally:
                painter.end()
    finally:
        if fasta is not None:
            fasta.close()
        for panel in panels:
            panel.deleteLater()
    return pages
