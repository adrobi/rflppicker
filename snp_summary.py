"""One summary per SNP line of the input file, independent of Qt.

The core returns candidate rows per normalized variant and a list of
rejected lines. The window, the Excel export and the CLI need the opposite
view: every line of the user's file in file order, with its schemes, its
status and a plain-language note. This module builds that view.
"""
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from i18n import N_, core_text, t
from palette import EXCEL_STATUS_FILLS
from reporting import ExtraSheet
from rflp_core import parse_fields, read_input_table


@dataclass(frozen=True)
class InputRecord:
    line: int          # line number in the file, as in report['rejected']['line']
    snp_id: str        # '' without an ID
    chrom: str
    pos: Optional[int]  # None when the line could not be parsed
    ref: str
    alt: str
    raw: str           # fields joined by the delimiter


def _record(line_no: int, parts: List[str], columns: Dict[str, Optional[int]], delim: str) -> InputRecord:
    def text(name: str, default_index: int) -> str:
        index = columns.get(name)
        index = default_index if index is None else index
        return parts[index] if index < len(parts) else ''

    id_index = columns.get('id')
    snp_id = parts[id_index].strip() if id_index is not None and id_index < len(parts) else ''
    raw = delim.join(parts)
    if len(parts) >= 4:
        try:
            chrom, pos, ref, alt = parse_fields(parts, columns)
            return InputRecord(line_no, snp_id, chrom, pos, ref, alt, raw)
        except (IndexError, ValueError):
            pass
    # The core rejects this line; keep what can be read for display.
    return InputRecord(line_no, snp_id, text('chrom', 0).replace('chr', '').strip(), None,
                       text('ref', 2).strip().upper(), text('alt', 3).strip().upper(), raw)


def read_input_records(path: str) -> Tuple[List[InputRecord], bool]:
    """Every data line the core reads, in file order, and whether the file has an ID column."""
    columns, lines, delim = read_input_table(path)
    return [_record(n, parts, columns, delim) for n, parts in lines], columns.get('id') is not None


ROBUST, NOMINAL, NONE, ERROR = 'robust', 'nominal', 'none', 'error'
STATUS_ORDER = (ROBUST, NOMINAL, NONE, ERROR)
STATUS_TITLES = {ROBUST: N_('Надёжно'), NOMINAL: N_('Только на хорошем геле'), NONE: N_('Нет схемы'),
                 ERROR: N_('Ошибка данных')}

_QUALITY_STATUS = {'robust': ROBUST, 'nominal_only': NOMINAL}
_ERROR_TEXTS = {
    'bad_chrom': N_('Хромосома не найдена в геноме'),
    'pos_oob': N_('Позиция за пределами хромосомы'),
    'bad_ref_chars': N_('Не одиночная замена A/C/G/T'),
    'bad_alt_chars': N_('Не одиночная замена A/C/G/T'),
    'minus_strand_like': N_('REF похож на обратную цепь'),
    'multi_alt_contains_ref': N_('ALT совпадает с геномом'),
    'ref_not_match': N_('REF не совпадает с геномом'),
}
_FAILURE_HINTS = {
    'no_primers': N_('Попробуйте расширить диапазон продукта или Tm'),
    'primer3_error': N_('Проверьте параметры Primer3 в окне «Дополнительно»'),
}


def status_title(status: str) -> str:
    return t(STATUS_TITLES[status])


def _rejection_text(rejection: dict) -> str:
    reason = rejection.get('reason') or ''
    details = core_text(rejection.get('details') or '')
    if reason == 'parse_error':
        return t('Строка не прочитана: {details}', details=details)
    return t(_ERROR_TEXTS[reason]) if reason in _ERROR_TEXTS else details or reason


@dataclass
class SnpSummary:
    record: InputRecord
    schemes: List[dict] = field(default_factory=list)   # core rows with an enzyme, '_index' = row index
    failure: Optional[dict] = None                       # core row without an enzyme
    rejection: Optional[dict] = None                     # entry of report['rejected']
    chosen: int = 0
    primer_base: str = ''
    swapped: bool = False   # normalization swapped REF and ALT: bands follow the genome's REF

    @property
    def label(self) -> str:
        r = self.record
        if r.pos is None:
            return r.raw
        return f'{r.chrom}:{r.pos} {r.ref}>{r.alt}'

    @property
    def scheme(self) -> Optional[dict]:
        return self.schemes[self.chosen] if self.schemes else None

    @property
    def manual(self) -> bool:
        return self.chosen != 0

    @property
    def status(self) -> str:
        if self.rejection is not None:
            return ERROR
        scheme = self.scheme
        if scheme is not None:
            return _QUALITY_STATUS.get(scheme.get('genotype_quality'), NONE)
        if self.failure is None or self.failure.get('status') == 'analysis_error':
            return ERROR
        return NONE

    def choose(self, index: int) -> None:
        if not 0 <= index < len(self.schemes):
            raise IndexError(f'No scheme {index} for {self.label}')
        self.chosen = index

    def swap_text(self) -> str:
        return t('REF и ALT переставлены по геному: REF = {ref}', ref=self.record.alt) if self.swapped else ''

    def note(self, dcaps_mode: str) -> str:
        """The «Примечание» text of the spec, §6.6."""
        if self.rejection is not None:
            return _rejection_text(self.rejection)
        scheme = self.scheme
        if scheme is None:
            return '; '.join(text for text in (self.swap_text(), self._failure_note(dcaps_mode)) if text)
        quality = scheme.get('genotype_quality')
        parts = [self.swap_text()] if self.swapped else []
        if quality not in _QUALITY_STATUS:
            parts.append(t('Три генотипа неразличимы по присутствию полос'))
        if self.manual:
            parts.append(t('Выбрано вручную, лучшая — {enzyme}', enzyme=self.schemes[0].get('enzyme')))
        if quality == 'nominal_only':
            parts.append(t('Только при исходном разрешении; труднее всего: {pair}',
                           pair=scheme.get('weakest_genotypes')))
        if scheme.get('design') == 'dCAPS':
            parts.append(f"dCAPS: {scheme.get('primer_mismatches')}")
        if scheme.get('primer_site'):
            parts.append(t('Сайт фермента в праймере'))
        if scheme.get('background_site_variants'):
            parts.append(t('Фоновые варианты меняют разрезы: {variants}',
                           variants=scheme['background_site_variants']))
        elif (scheme.get('background_variants') or 0) > 0:
            parts.append(t('Известных вариантов в ампликоне: {count}', count=scheme['background_variants']))
        if (scheme.get('repeat_pct') or 0) >= 50:
            parts.append(t('Повторы {percent}% ампликона', percent=scheme['repeat_pct']))
        if scheme.get('status') == 'window_only':
            parts.append(t('Поиск по окну, праймеры не подобраны'))
        return '; '.join(parts)

    def _failure_note(self, dcaps_mode: str) -> str:
        if self.failure is None:
            return t('Нет результата для этой строки')
        status = self.failure.get('status')
        reason = core_text(self.failure.get('reason') or '')
        if status == 'analysis_error':
            return reason
        hint = t(_FAILURE_HINTS[status]) if status in _FAILURE_HINTS else ''
        if status in ('no_discriminating_assay', 'no_enzyme_found'):
            hint = (t('Попробуйте другой диапазон продукта') if dcaps_mode == 'always'
                    else t('Попробуйте dCAPS «Для каждого SNP»'))
        return '. '.join(text for text in (reason, hint) if text)


def build_summaries(records: List[InputRecord], header: List[str], rows: List[list],
                    rejected: List[dict]) -> List[SnpSummary]:
    # Only parse errors carry the file line. With skip_norm the core numbers
    # rejected variants 1, 2, ..., so entries with coordinates match by value.
    by_line = {item['line']: item for item in rejected
               if item.get('line') is not None and item.get('chrom') is None}
    by_values = {}
    for item in rejected:
        if item.get('chrom') is not None and item.get('pos') is not None:
            by_values.setdefault((str(item['chrom']), int(item['pos']), item.get('ref'), item.get('alt')), item)
    # The core does not deduplicate variants: a repeated line (or the same SNP
    # with REF/ALT swapped) yields an identical block of rows. Keep one copy.
    groups: Dict[str, List[dict]] = {}
    seen = set()
    for index, row in enumerate(rows):
        key = tuple(row)
        if key in seen:
            continue
        seen.add(key)
        values = dict(zip(header, row), _index=index)
        groups.setdefault(values.get('variant'), []).append(values)

    summaries = []
    for r in records:
        rejection = by_line.get(r.line) or by_values.get((r.chrom, r.pos, r.ref, r.alt))
        found, swapped = [], False
        if rejection is None and r.pos is not None:
            found = groups.get(f'{r.chrom}:{r.pos} {r.ref}>{r.alt}') or []
            if not found:
                found = groups.get(f'{r.chrom}:{r.pos} {r.alt}>{r.ref}') or []
                swapped = bool(found)
        summaries.append(SnpSummary(
            record=r,
            schemes=[values for values in found if values.get('enzyme')],
            failure=next((values for values in found if not values.get('enzyme')), None),
            rejection=rejection,
            swapped=swapped,
        ))
    assign_primer_bases(summaries)
    return summaries


def _clean_name(text: str) -> Tuple[str, int]:
    """Primer-name characters only, and how many characters had to be replaced."""
    cleaned, replaced = re.subn(r'[^A-Za-z0-9_.-]', '_', text)
    return cleaned, replaced


def assign_primer_bases(summaries: List[SnpSummary]) -> None:
    """Unique primer name stems: the SNP ID, or the coordinate when the ID is empty or mostly unusable."""
    seen: Dict[str, int] = {}
    for s in summaries:
        base, replaced = _clean_name(s.record.snp_id)
        if not s.record.snp_id or replaced > len(re.findall(r'[A-Za-z0-9]', base)):
            base = _clean_name(f'{s.record.chrom}_{s.record.pos}')[0]
        seen[base] = seen.get(base, 0) + 1
        s.primer_base = base if seen[base] == 1 else f'{base}_{seen[base]}'


def order_lines(summary: SnpSummary) -> str:
    scheme = summary.scheme
    if scheme is None or not scheme.get('primers_left'):
        return ''
    return (f"{summary.primer_base}_F\t{scheme['primers_left']}\n"
            f"{summary.primer_base}_R\t{scheme['primers_right']}")


def summary_text(summary: SnpSummary, dcaps_mode: str) -> str:
    """The «Копировать всё» text for a lab notebook."""
    snp_id = summary.record.snp_id
    lines = [f'{summary.label} ({snp_id})' if snp_id else summary.label,
             t('Статус: {status}', status=status_title(summary.status))]
    scheme = summary.scheme
    if scheme is not None:
        lines.append(t('Фермент: {enzyme}, сайт {site}', enzyme=scheme.get('enzyme'), site=scheme.get('site')))
        for side, primer, tm in (('F', 'primers_left', 'tm_left'), ('R', 'primers_right', 'tm_right')):
            if scheme.get(primer):
                tm_text = f" (Tm {scheme[tm]:.1f} °C)" if scheme.get(tm) is not None else ''
                lines.append(f'{summary.primer_base}_{side}: {scheme[primer]}{tm_text}')
        lines.append(t('Продукт: {size} п.н.', size=scheme.get('product_size')))
        lines.append(t('Полосы: REF/REF {ref_ref}; REF/ALT {ref_alt}; ALT/ALT {alt_alt}',
                       ref_ref=scheme.get('bands_ref_ref'), ref_alt=scheme.get('bands_ref_alt'),
                       alt_alt=scheme.get('bands_alt_alt')))
    note = summary.note(dcaps_mode)
    if note:
        lines.append(t('Примечание: {note}', note=note))
    return '\n'.join(lines)


def weakest_text(scheme: dict) -> str:
    """The line under the gel: which genotype pair the score rests on."""
    text = t('Труднее всего различить: {pair}', pair=scheme.get('weakest_genotypes'))
    if scheme.get('worst_margin') is not None:
        text += t('; оценка при ухудшении {score}', score=f"{scheme['worst_margin']:.1f}")
    return text


SUMMARY_HEADER = ['ID', 'SNP', N_('Хромосома в геноме'), N_('Статус'), N_('Выбор'), N_('Фермент'), N_('Сайт'),
                  N_('Изошизомеры'), N_('Поставщики'), N_('Тип схемы'), N_('Праймер F 5′→3′'),
                  N_('Праймер R 5′→3′'), 'Tm F', 'Tm R', N_('Продукт, п.н.'), N_('Полосы REF/REF'),
                  N_('Полосы REF/ALT'), N_('Полосы ALT/ALT'), N_('Оценка'), N_('Оценка при ухудшении'),
                  N_('Примечание'), N_('Строка файла')]
PRIMER_HEADER = [N_('Имя'), N_('Последовательность 5′→3′'), N_('Длина'), 'Tm', 'SNP', N_('Фермент')]
_DESIGN_TITLES = {'natural': N_('природный сайт'), 'dCAPS': 'dCAPS'}
_SCHEME_COLUMNS = ('enzyme', 'site', 'isoschizomers', 'suppliers', 'design', 'primers_left', 'primers_right',
                   'tm_left', 'tm_right', 'product_size', 'bands_ref_ref', 'bands_ref_alt', 'bands_alt_alt',
                   'genotype_margin', 'worst_margin')


def summary_table(summaries: List[SnpSummary], dcaps_mode: str) -> Tuple[List[str], List[list]]:
    rows = []
    for s in summaries:
        scheme = s.scheme or {}
        source = s.scheme or s.failure or {}
        values = [scheme.get(name) for name in _SCHEME_COLUMNS]
        design = scheme.get('design')
        values[_SCHEME_COLUMNS.index('design')] = t(_DESIGN_TITLES[design]) if design in _DESIGN_TITLES else design
        choice = (t('вручную') if s.manual else t('авто')) if s.scheme is not None else ''
        rows.append([s.record.snp_id, s.label, source.get('mapped_id'), status_title(s.status), choice,
                     *values, s.note(dcaps_mode), s.record.line])
    return [t(title) for title in SUMMARY_HEADER], rows


def primer_table(summaries: List[SnpSummary]) -> Tuple[List[str], List[list]]:
    rows = []
    for s in summaries:
        scheme = s.scheme
        if scheme is None or not scheme.get('primers_left'):
            continue
        for side, primer, tm in (('F', 'primers_left', 'tm_left'), ('R', 'primers_right', 'tm_right')):
            sequence = scheme.get(primer) or ''
            rows.append([f'{s.primer_base}_{side}', sequence, len(sequence), scheme.get(tm), s.label,
                         scheme.get('enzyme')])
    return [t(title) for title in PRIMER_HEADER], rows


def excel_sheets(summaries: List[SnpSummary], dcaps_mode: str) -> List[ExtraSheet]:
    """The «Сводка» and «Праймеры» worksheets written before Results, in the interface language."""
    header, rows = summary_table(summaries, dcaps_mode)
    primer_header, primer_rows = primer_table(summaries)
    fills = {status_title(status): EXCEL_STATUS_FILLS[status] for status in STATUS_ORDER}
    return [ExtraSheet(t('Сводка'), header, rows, 'RFLPSummary', {t('Статус'): fills}),
            ExtraSheet(t('Праймеры'), primer_header, primer_rows, 'RFLPPrimers')]


def results_table(summaries: List[SnpSummary], header: List[str], rows: List[list]) -> Tuple[List[str], List[list]]:
    """Core rows of these SNPs in core order, with the SNP ID and whether the row is the chosen scheme."""
    owner: Dict[str, SnpSummary] = {}
    selected = set()
    for s in summaries:
        for values in s.schemes + ([s.failure] if s.failure else []):
            owner.setdefault(values.get('variant'), s)
        if s.scheme is not None:
            selected.add(s.scheme['_index'])
    variant = header.index('variant')
    table = [list(row) + [owner[row[variant]].record.snp_id, index in selected]
             for index, row in enumerate(rows) if row[variant] in owner]
    return list(header) + ['snp_id', 'selected'], table
