"""Published observations for identifiable commercial preparations.

Source papers cited below use CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/).
The entries are abbreviated factual summaries of the cited methods or tables.
Missing properties are deliberately left unset. Do not infer methylation
sensitivity from the recognition sequence or transfer conditions to an isoschizomer.
"""
from dataclasses import dataclass

from Bio.Restriction import AllEnzymes
from i18n import N_, t


@dataclass(frozen=True)
class Product:
    name: str
    supplier: str
    number: str
    citation: str
    doi: str
    temperature: int | None = None
    buffer: str = ''
    cpg: str = ''
    dam: str = ''
    dcm: str = ''
    star: str = ''
    star_source: str = ''

    @property
    def label(self):
        return f'{self.name} · {self.supplier} {self.number}'


# The main citation covers all values except separately cited star activity.
# CutSmart is the buffer named in the 2014 paper, not a later formulation.
RRBS = ('Martinez-Arguelles et al. (2014), Table 1', '10.1186/1756-0500-7-534')
PRODUCTS = {
    'EcoRI': (
        Product('EcoRI', 'NEB', 'R0101S', 'Hoog et al. (2024)', '10.1038/s41467-024-48037-2',
                37, 'EcoRI Buffer'),
        Product('EcoRI-HF', 'NEB', 'R3101', 'Díaz-Viraqué et al. (2023)', '10.1038/s41564-023-01483-y',
                37, 'rCutSmart', star=N_('низкая'),
                star_source='Lefrançois et al. (2016), doi.org/10.1371/journal.pgen.1006347'),
    ),
    'MspI': (Product('MspI', 'NEB', 'R0106S', *RRBS, 37, 'CutSmart'),),
    'MseI': (Product('MseI', 'NEB', 'R0525S', *RRBS, 37, 'CutSmart'),),
    'AciI': (Product('AciI', 'NEB', 'R0551S', 'Feng et al. (2025)', '10.1155/jdr/4342250',
                     cpg=N_('чувствителен')),)
}


def property_lines(scheme, sequences):
    """Identical localized product data for the card and its printed sheet."""
    name = scheme.get('enzyme')
    products = PRODUCTS.get(name, ())
    selected = scheme.get('enzyme_product')
    product = next((p for p in products if p.number == selected), None)
    if len(products) == 1 and selected is None:
        product = products[0]
    missing = t('нет данных')
    if product is None:
        lines = [t('Препарат не выбран') if products else t('Свойства фермента: нет данных')]
    else:
        lines = [
            product.label,
            t('В публикации: инкубация {temp}, буфер {buffer}',
              temp=f'{product.temperature} °C' if product.temperature is not None else missing,
              buffer=product.buffer or missing),
            t('Метилирование: CpG — {cpg}, Dam — {dam}, Dcm — {dcm}',
              cpg=t(product.cpg) or missing, dam=t(product.dam) or missing, dcm=t(product.dcm) or missing),
            t('Star-активность: {value}', value=t(product.star) or missing)
            + (f' · {product.star_source}' if product.star_source else ''),
            f'{product.citation} · doi.org/{product.doi} · CC BY 4.0',
        ]
    if sequences:
        enzyme = next((e for e in AllEnzymes if e.__name__ == name), None)
        if enzyme and any('CG' in seq[max(0, match.start()-1):match.start()+enzyme.size+1]
                          for seq in sequences if seq for match in enzyme.compsite.finditer(seq)):
            if scheme.get('status') == 'window_only':
                lines.append(t('Сайт перекрывает CpG. Для геномной ДНК учитывайте метилирование образца.'))
            else:
                lines.append(t('Сайт перекрывает CpG. Обычный ПЦР-продукт не сохраняет метилирование матрицы.'))
    return lines
