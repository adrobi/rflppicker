"""One summary per SNP line of the input file, independent of Qt.

The core returns candidate rows per normalized variant and a list of
rejected lines. The window, the Excel export and the CLI need the opposite
view: every line of the user's file in file order, with its schemes, its
status and a plain-language note. This module builds that view.
"""
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

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
STATUS_TITLES = {ROBUST: 'Надёжно', NOMINAL: 'Только на хорошем геле', NONE: 'Нет схемы', ERROR: 'Ошибка данных'}

_QUALITY_STATUS = {'robust': ROBUST, 'nominal_only': NOMINAL}
_ERROR_TEXTS = {
    'bad_chrom': 'Хромосома не найдена в геноме',
    'pos_oob': 'Позиция за пределами хромосомы',
    'bad_ref_chars': 'Не одиночная замена A/C/G/T',
    'bad_alt_chars': 'Не одиночная замена A/C/G/T',
    'minus_strand_like': 'REF похож на обратную цепь',
    'multi_alt_contains_ref': 'ALT совпадает с геномом',
    'ref_not_match': 'REF не совпадает с геномом',
}
_FAILURE_HINTS = {
    'no_primers': 'Попробуйте расширить диапазон продукта или Tm',
    'primer3_error': 'Проверьте параметры Primer3 в окне «Дополнительно»',
}


def _rejection_text(rejection: dict) -> str:
    reason = rejection.get('reason') or ''
    details = rejection.get('details') or ''
    if reason == 'parse_error':
        return f'Строка не прочитана: {details}'
    return _ERROR_TEXTS.get(reason) or details or reason


@dataclass
class SnpSummary:
    record: InputRecord
    schemes: List[dict] = field(default_factory=list)   # core rows with an enzyme, '_index' = row index
    failure: Optional[dict] = None                       # core row without an enzyme
    rejection: Optional[dict] = None                     # entry of report['rejected']
    chosen: int = 0
    primer_base: str = ''

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

    def note(self, dcaps_mode: str) -> str:
        """The «Примечание» text of the spec, §6.6."""
        if self.rejection is not None:
            return _rejection_text(self.rejection)
        scheme = self.scheme
        if scheme is None:
            return self._failure_note(dcaps_mode)
        quality = scheme.get('genotype_quality')
        parts = []
        if quality not in _QUALITY_STATUS:
            parts.append('Три генотипа неразличимы по присутствию полос')
        if self.manual:
            parts.append(f"Выбрано вручную, лучшая — {self.schemes[0].get('enzyme')}")
        if quality == 'nominal_only':
            parts.append(f"Только при исходном разрешении; труднее всего: {scheme.get('weakest_genotypes')}")
        if scheme.get('design') == 'dCAPS':
            parts.append(f"dCAPS: {scheme.get('primer_mismatches')}")
        if scheme.get('primer_site'):
            parts.append('Сайт фермента в праймере')
        if scheme.get('background_site_variants'):
            parts.append(f"Фоновые варианты меняют разрезы: {scheme['background_site_variants']}")
        elif (scheme.get('background_variants') or 0) > 0:
            parts.append(f"Известных вариантов в ампликоне: {scheme['background_variants']}")
        if (scheme.get('repeat_pct') or 0) >= 50:
            parts.append(f"Повторы {scheme['repeat_pct']}% ампликона")
        if scheme.get('status') == 'window_only':
            parts.append('Поиск по окну, праймеры не подобраны')
        return '; '.join(parts)

    def _failure_note(self, dcaps_mode: str) -> str:
        if self.failure is None:
            return 'Нет результата для этой строки'
        status = self.failure.get('status')
        reason = self.failure.get('reason') or ''
        if status == 'analysis_error':
            return reason
        hint = _FAILURE_HINTS.get(status)
        if status in ('no_discriminating_assay', 'no_enzyme_found'):
            hint = ('Попробуйте другой диапазон продукта' if dcaps_mode == 'always'
                    else 'Попробуйте dCAPS «Для каждого SNP»')
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
    groups: Dict[str, List[dict]] = {}
    for index, row in enumerate(rows):
        values = dict(zip(header, row), _index=index)
        groups.setdefault(values.get('variant'), []).append(values)

    summaries = []
    for r in records:
        rejection = by_line.get(r.line) or by_values.get((r.chrom, r.pos, r.ref, r.alt))
        found = []
        if rejection is None and r.pos is not None:
            found = (groups.get(f'{r.chrom}:{r.pos} {r.ref}>{r.alt}')
                     or groups.get(f'{r.chrom}:{r.pos} {r.alt}>{r.ref}') or [])
        summaries.append(SnpSummary(
            record=r,
            schemes=[values for values in found if values.get('enzyme')],
            failure=next((values for values in found if not values.get('enzyme')), None),
            rejection=rejection,
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
             f'Статус: {STATUS_TITLES[summary.status]}']
    scheme = summary.scheme
    if scheme is not None:
        lines.append(f"Фермент: {scheme.get('enzyme')}, сайт {scheme.get('site')}")
        for side, primer, tm in (('F', 'primers_left', 'tm_left'), ('R', 'primers_right', 'tm_right')):
            if scheme.get(primer):
                tm_text = f" (Tm {scheme[tm]:.1f} °C)" if scheme.get(tm) is not None else ''
                lines.append(f'{summary.primer_base}_{side}: {scheme[primer]}{tm_text}')
        lines.append(f"Продукт: {scheme.get('product_size')} п.н.")
        lines.append(f"Полосы: REF/REF {scheme.get('bands_ref_ref')}; REF/ALT {scheme.get('bands_ref_alt')}; "
                     f"ALT/ALT {scheme.get('bands_alt_alt')}")
    note = summary.note(dcaps_mode)
    if note:
        lines.append(f'Примечание: {note}')
    return '\n'.join(lines)


def weakest_text(scheme: dict) -> str:
    """The line under the gel: which genotype pair the score rests on."""
    text = f"Труднее всего различить: {scheme.get('weakest_genotypes')}"
    if scheme.get('worst_margin') is not None:
        text += f"; оценка при ухудшении {scheme['worst_margin']:.1f}"
    return text


SUMMARY_HEADER = ['ID', 'SNP', 'Хромосома в геноме', 'Статус', 'Выбор', 'Фермент', 'Сайт', 'Изошизомеры',
                  'Поставщики', 'Тип схемы', 'Праймер F 5′→3′', 'Праймер R 5′→3′', 'Tm F', 'Tm R',
                  'Продукт, п.н.', 'Полосы REF/REF', 'Полосы REF/ALT', 'Полосы ALT/ALT', 'Оценка',
                  'Оценка при ухудшении', 'Примечание', 'Строка файла']
PRIMER_HEADER = ['Имя', 'Последовательность 5′→3′', 'Длина', 'Tm', 'SNP', 'Фермент']
_DESIGN_TITLES = {'natural': 'природный сайт', 'dCAPS': 'dCAPS'}
_SCHEME_COLUMNS = ('enzyme', 'site', 'isoschizomers', 'suppliers', 'design', 'primers_left', 'primers_right',
                   'tm_left', 'tm_right', 'product_size', 'bands_ref_ref', 'bands_ref_alt', 'bands_alt_alt',
                   'genotype_margin', 'worst_margin')


def summary_table(summaries: List[SnpSummary], dcaps_mode: str) -> Tuple[List[str], List[list]]:
    rows = []
    for s in summaries:
        scheme = s.scheme or {}
        source = s.scheme or s.failure or {}
        values = [scheme.get(name) for name in _SCHEME_COLUMNS]
        values[_SCHEME_COLUMNS.index('design')] = _DESIGN_TITLES.get(scheme.get('design'), scheme.get('design'))
        choice = ('вручную' if s.manual else 'авто') if s.scheme is not None else ''
        rows.append([s.record.snp_id, s.label, source.get('mapped_id'), STATUS_TITLES[s.status], choice,
                     *values, s.note(dcaps_mode), s.record.line])
    return list(SUMMARY_HEADER), rows


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
    return list(PRIMER_HEADER), rows


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
