"""Translate lookup status and save a verified rsID as an ordinary input row."""

from pathlib import Path

from i18n import N_, core_text, t
from reporting import _atomic_destination


_STATUSES = {
    'accepted': N_('Соответствие подтверждено'),
    'ambiguous': N_('Неоднозначное соответствие'),
    'ref_mismatch': N_('REF не совпадает с геномом'),
    'contig_not_found': N_('Контиг не найден в FASTA'),
    'no_match': N_('Соответствие не найдено'),
    'assembly_mismatch': N_('Сборка не совпадает'),
    'species_mismatch': N_('Вид не совпадает'),
    'unsupported_allele': N_('Нет явных аллелей A/C/G/T'),
    'source_error': N_('Ошибка источника rsID'),
}


def describe_lookup(result):
    lines = [result.rsid + ' · ' + t(_STATUSES.get(result.status, result.status))]
    if result.detail:
        lines.append(core_text(result.detail))
    lines.extend(f'{item.species + " · " if item.species else ""}TaxID {item.taxid} · {item.assembly} · {item.contig}:{item.position} '
                 f'{item.ref}>{",".join(item.alts)}' for item in result.matches)
    return '\n'.join(lines)


def write_rsid_target(path, result, *sources):
    item = result.match
    if item is None:
        raise ValueError('Нет подтверждённого соответствия rsID')
    if any(source and Path(source).resolve() == Path(path).resolve() for source in sources):
        raise ValueError('Выберите выходной файл отдельно от источника rsID и FASTA')
    with _atomic_destination(path) as temporary:
        with temporary.open('w', encoding='utf-8', newline='\n') as stream:
            stream.write(f'##fileformat=VCFv4.3\n##reference={item.assembly}\n'
                         f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                         f'{item.contig}\t{item.position}\t{result.rsid}\t{item.ref}\t{",".join(item.alts)}\t.\t.\t.\n')
