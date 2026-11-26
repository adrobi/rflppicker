"""Sourced properties for identifiable commercial preparations.

Source papers cited below use CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/).
The entries summarize facts from the cited papers and individual product cards.
Missing properties are deliberately left unset. Do not infer methylation
sensitivity from the recognition sequence or transfer conditions to an isoschizomer.
"""
from dataclasses import dataclass

from Bio.Restriction import AllEnzymes
from i18n import N_, t

CATALOG_VERSION = 'product-observations-v2'


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
    url: str = ''

    @property
    def label(self):
        return f'{self.name} · {self.supplier} {self.number}'


# The main citation covers all values except separately cited star activity.
# CutSmart is the buffer named in the 2014 paper, not a later formulation.
RRBS = ('Martinez-Arguelles et al. (2014), Table 1', '10.1186/1756-0500-7-534')
SIB_ECORI = 'https://sibenzyme.com/product/ecor-i/'
TF_ECORI = 'https://documents.thermofisher.com/TFS-Assets/LSG/manuals/MAN0012438_FastDigest_EcoRI_800uL_UG.pdf'
TF_MSPI = 'https://documents.thermofisher.com/TFS-Assets/LSG/manuals/MAN0012468_FastDigest_MspI_UG.pdf'
_NO_EFFECT = N_('нет эффекта')
PRODUCTS = {
    'EcoRI': (
        Product('EcoRI', 'NEB', 'R0101S', 'Hoog et al. (2024)', '10.1038/s41467-024-48037-2',
                37, 'EcoRI Buffer'),
        Product('EcoRI-HF', 'NEB', 'R3101', 'Díaz-Viraqué et al. (2023)', '10.1038/s41564-023-01483-y',
                37, 'rCutSmart', star=N_('низкая'),
                star_source='Lefrançois et al. (2016), doi.org/10.1371/journal.pgen.1006347'),
        Product('EcoR I', 'SibEnzyme', 'E057', 'SibEnzyme', '', 37, 'SE EcoRI',
                star=N_('возможна при избытке фермента или неоптимальном буфере'), url=SIB_ECORI),
        Product('EcoR I Turbo', 'SibEnzyme', 'E057T', 'SibEnzyme', '', 37, 'ROSE', url=SIB_ECORI),
        Product('FastDigest EcoRI', 'Thermo Scientific', 'FD0274', 'Thermo Scientific', '', 37,
                'FastDigest / FastDigest Green', cpg=N_('может снижать активность'), dam=_NO_EFFECT, dcm=_NO_EFFECT,
                star=N_('возможна при длительной инкубации'), url=TF_ECORI),
    ),
    'MspI': (
        Product('MspI', 'NEB', 'R0106S', *RRBS, 37, 'CutSmart'),
        Product('FastDigest MspI', 'Thermo Scientific', 'FD0544', 'Thermo Scientific', '', 37,
                'FastDigest / FastDigest Green', cpg=_NO_EFFECT, dam=_NO_EFFECT, dcm=_NO_EFFECT,
                star=N_('не обнаружена в испытании производителя'), url=TF_MSPI),
    ),
    'HpaII': (Product('Hpa II', 'SibEnzyme', 'E161', 'SibEnzyme', '', 37, 'SE B', cpg=N_('чувствителен'),
                      url='https://sibenzyme.com/product/hpa-ii/'),),
    'SspI': (Product('Ssp I', 'SibEnzyme', 'E041', 'SibEnzyme', '', 37, 'SE K',
                     star=N_('возможна при избытке фермента'), url='https://sibenzyme.com/product/ssp-i/'),),
    'MseI': (Product('MseI', 'NEB', 'R0525S', *RRBS, 37, 'CutSmart'),),
    'AciI': (Product('AciI', 'NEB', 'R0551S', 'Feng et al. (2025)', '10.1155/jdr/4342250',
                     cpg=N_('чувствителен')),)
}


def property_lines(scheme, sequences):
    """Identical localized product data for the card and its printed sheet."""
    return [line for group in product_sections(scheme, sequences) for line in group]


def product_sections(scheme, sequences):
    if scheme.get('design') == 'amplicon_length':
        return []
    names = (scheme['enzyme_1'], scheme['enzyme_2']) if scheme.get('enzyme_2') else (scheme.get('enzyme'),)
    sections = []
    for name, key in zip(names, ('enzyme_product', 'enzyme_product_2')):
        lines = _product_lines(name, scheme.get(key), sequences, scheme.get('status'))
        if len(names) == 2 and not lines[0].startswith(name):
            lines[0] = f'{name} · {lines[0]}'
        sections.append(lines)
    return sections


def _product_lines(name, selected, sequences, status):
    products = PRODUCTS.get(name, ())
    product = next((p for p in products if p.number == selected), None)
    if len(products) == 1 and selected is None:
        product = products[0]
    missing = t('нет данных')
    if product is None:
        lines = [t('Препарат не выбран') if products else t('Свойства фермента: нет данных')]
    else:
        lines = [
            product.label,
            t(N_('У производителя: инкубация {temp}, буфер {buffer}') if product.url else
              N_('В публикации: инкубация {temp}, буфер {buffer}'),
              temp=f'{product.temperature} °C' if product.temperature is not None else missing, buffer=product.buffer or missing),
            t('Метилирование: CpG — {cpg}, Dam — {dam}, Dcm — {dcm}',
              cpg=t(product.cpg) or missing, dam=t(product.dam) or missing, dcm=t(product.dcm) or missing),
            t('Star-активность: {value}', value=t(product.star) or missing)
            + (f' · {product.star_source}' if product.star_source else ''),
            product.url or f'{product.citation} · doi.org/{product.doi} · CC BY 4.0',
        ]
    if sequences:
        enzyme = next((e for e in AllEnzymes if e.__name__ == name), None)
        if enzyme and any('CG' in seq[max(0, match.start()-1):match.start()+enzyme.size+1]
                          for seq in sequences if seq for match in enzyme.compsite.finditer(seq)):
            if status == 'window_only':
                lines.append(t('Сайт перекрывает CpG. Для геномной ДНК учитывайте метилирование образца.'))
            else:
                lines.append(t('Сайт перекрывает CpG. Обычный ПЦР-продукт не сохраняет метилирование матрицы.'))
    return lines
