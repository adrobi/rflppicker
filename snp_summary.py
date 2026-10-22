"""One summary per SNP line of the input file, independent of Qt.

The core returns candidate rows per normalized variant and a list of
rejected lines. The window, the Excel export and the CLI need the opposite
view: every line of the user's file in file order, with its schemes, its
status and a plain-language note. This module builds that view.
"""
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
    return summaries
