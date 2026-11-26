"""VCF targets. Preserve the file line number and the order of ALT alleles."""

import gzip
import re


VCF_COLUMNS = dict(chrom=0, pos=1, ref=2, alt=3, id=4)
VCF_DNA = frozenset('ACGT')
VCF_REJECTION_TEXT = {
    'parse_error': 'Ожидалось минимум 8 полей VCF с целой позицией',
    'vcf_symbolic': 'У ALT нет явной последовательности в VCF',
    'vcf_no_selected_alt': 'В известных генотипах выбранных образцов нет ALT',
    'vcf_no_alt': 'В записи VCF нет ALT',
}


def is_vcf_path(path):
    return str(path).lower().endswith(('.vcf', '.vcf.gz', '.vcf.bgz', '.vcf.bgzf'))


def open_vcf(path):
    with open(path, 'rb') as stream:
        compressed = stream.read(2) == b'\x1f\x8b'
    return (gzip.open if compressed else open)(path, 'rt', encoding='utf-8-sig')


def parse_region(value):
    if not value:
        return None
    name, separator, bounds = value.rpartition(':')
    if not separator:
        return value, 1, None
    match = re.fullmatch(r'([1-9]\d*)(?:-([1-9]\d*))?', bounds)
    if not name or not match:
        raise ValueError('Некорректный регион VCF. Ожидается контиг или контиг:начало-конец')
    left = int(match[1])
    right = int(match[2]) if match[2] else left
    if left > right:
        raise ValueError('Некорректный регион VCF. Ожидается контиг или контиг:начало-конец')
    return name, left, right


def _observed_indices(fields, selected, allele_count):
    if len(fields) < 9 or 'GT' not in fields[8].split(':'):
        return set()
    gt = fields[8].split(':').index('GT')
    found = set()
    for index in selected:
        if index >= len(fields):
            continue
        values = fields[index].split(':')
        if gt < len(values):
            for allele in re.split(r'[/|]', values[gt]):
                if allele.isascii() and allele.isdecimal() and len(allele) <= len(str(allele_count)):
                    number = int(allele)
                    if 1 <= number <= allele_count:
                        found.add(number)
    return found


def vcf_targets(path, region='', samples=(), check_cancel=None):
    """Yield (line, table fields, rejection reason, skipped symbolic ALTs)."""
    scope = parse_region(region)
    with open_vcf(path) as stream:
        selected = None
        for line_no, line in enumerate(stream, 1):
            if check_cancel is not None:
                check_cancel()
            if line.startswith('#CHROM\t'):
                names = line.rstrip('\r\n').split('\t')[9:]
                unknown = sorted(set(samples) - set(names))
                if unknown:
                    raise ValueError('Неизвестные образцы VCF: ' + ', '.join(unknown))
                selected = [9 + names.index(name) for name in samples]
            if line.startswith('#') or not line.strip():
                continue
            if selected is None:
                raise ValueError('В VCF нет строки #CHROM')
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) < 8:
                yield line_no, fields, 'parse_error', ()
                continue
            chrom, pos, name, ref, alt = fields[:5]
            try:
                position = int(pos)
            except ValueError:
                position = 0
            if scope and (chrom != scope[0] or position < scope[1] or
                          (scope[2] is not None and position > scope[2])):
                continue
            row = [chrom, pos, ref, alt, '' if name == '.' else name]
            alleles = alt.upper().split(',')
            included = _observed_indices(fields, selected, len(alleles)) if selected else None
            observed = [a for index, a in enumerate(alleles, 1) if included is None or index in included]
            explicit = [a for a in observed if a and set(a) <= VCF_DNA]
            omitted = tuple(a for a in observed if a != '.' and a not in explicit)
            if explicit:
                row[3] = ','.join(explicit)
                yield line_no, row, '', omitted
            else:
                reason = 'vcf_symbolic' if omitted else 'vcf_no_selected_alt' if selected else 'vcf_no_alt'
                yield line_no, row, reason, omitted
        if selected is None:
            raise ValueError('В VCF нет строки #CHROM')
