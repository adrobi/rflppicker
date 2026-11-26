"""Labels and band tables shared by the two-reaction PCR card and lab sheet."""

from i18n import N_, t

_PRIMER_LABELS = {
    'allele_specific': N_('Аллель-специфический праймер'),
    'common': N_('Общий праймер'),
    'control_only': N_('Праймер внутреннего контроля'),
}


def primer_label(role):
    return t(_PRIMER_LABELS[role])


def tube_heading(tube):
    return t('Пробирка {reaction} · аллель {allele}', reaction=tube['reaction'], allele=tube['allele'])


def tube_bands(tube):
    return t('Контроль {control} п.н. · диагностическая полоса {product} п.н.',
             control=tube['control_size'], product=tube['product_size'])


def genotype_headings():
    return t('Генотип'), t('Пробирка {reaction}', reaction='REF'), t('Пробирка {reaction}', reaction='ALT')


def genotype_bands(tubes):
    return {genotype: {tube['reaction']: tube['genotype_bands'][genotype] for tube in tubes}
            for genotype in ('REF/REF', 'REF/ALT', 'ALT/ALT')}


def genotype_rows(tubes):
    return [(genotype, ' + '.join(map(str, bands['REF'])), ' + '.join(map(str, bands['ALT'])))
            for genotype, bands in genotype_bands(tubes).items()]


def unresolved_text():
    return t('Если хотя бы одна контрольная полоса отсутствует или отсутствуют обе диагностические полосы, '
             'генотип не определён.')
