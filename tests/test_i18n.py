import ast
import re
import string
import unittest
from pathlib import Path
from unittest.mock import patch

import i18n
import translations_en

ROOT = Path(__file__).resolve().parent.parent
CYRILLIC = re.compile('[А-Яа-яЁё]')

# Every Russian text of rflp_core / assay_scoring that reaches a person (i18n spec §4.3).
CORE_SAMPLES = [
    ('Primer3 не вернул подходящих пар', 'Primer3 returned no suitable pairs'),
    ('SNP ближе 100 п.н. к краю окна; увеличьте фланг',
     'SNP is closer than 100 bp to the window edge; increase the flank'),
    ('Не найдено сочетание для различения трёх генотипов',
     'No primer–enzyme combination distinguishes the three genotypes'),
    ('Подходящий фермент не найден', 'No suitable enzyme found'),
    ('REF не совпадает с FASTA', 'REF does not match the FASTA'),
    ('; окон с неопределёнными основаниями: 2', '; windows with ambiguous bases: 2'),
    ('; ошибок расчёта ферментов: 1', '; enzyme calculation errors: 1'),
    ('Три генотипа различимы во всём заданном диапазоне разрешения',
     'The three genotypes are distinguishable over the whole resolution range'),
    ('Три генотипа различимы только при исходном разрешении',
     'The three genotypes are distinguishable at the nominal resolution only'),
    ('Три генотипа неразличимы по присутствию полос',
     'The three genotypes cannot be told apart by band presence'),
    ('Оценка окна без подобранных праймеров. ', 'Window estimate without designed primers. '),
    ('dCAPS: несовпадение в праймере создаёт сайт в аллеле ALT. ',
     'dCAPS: the primer mismatch creates the site in the ALT allele. '),
    ('Модуль Primer3 недоступен', 'The Primer3 module is not available'),
    ('Primer3 вернул праймеры вне окна или перекрывающие SNP',
     'Primer3 returned primers outside the window or over the SNP'),
    ('Несогласованные размеры продукта/праймеров Primer3', 'Inconsistent Primer3 product or primer sizes'),
    ('Последовательности праймеров не совпадают с референсом', 'Primer sequences do not match the reference'),
    ('Primer3 вернул некорректную оценку', 'Primer3 returned an invalid score'),
    ('Ожидалось минимум 4 поля', 'At least 4 fields expected'),
    ('Поддерживаются только одиночные SNP A/C/G/T', 'Only single-base A/C/G/T SNPs are supported'),
    ("CHROM 'NC_1.1' отсутствует в FASTA", "CHROM 'NC_1.1' is not in the FASTA"),
    ('Позиция вне FASTA', 'Position outside the FASTA'),
    ('Не задан параметр: fasta', 'Parameter not set: fasta'),
    (r'Файл не найден (fasta): D:\x.fa', r'File not found (fasta): D:\x.fa'),
    ('Некорректный параметр: prod_min/prod_max', 'Invalid parameter: prod_min/prod_max'),
    ('Минимальный размер продукта больше максимального', 'Minimum product size is greater than maximum'),
    ('Минимальный размер фрагмента больше максимального', 'Minimum fragment size is greater than maximum'),
    ('Минимальный Tm больше максимального', 'Minimum Tm is greater than maximum'),
    ('Минимальный размер продукта dCAPS больше максимального', 'Minimum dCAPS product size is greater than maximum'),
    ('Минимальный Tm праймеров dCAPS больше максимального', 'Minimum dCAPS primer Tm is greater than maximum'),
    ('Некорректный параметр: dcaps_prod_min/dcaps_prod_max', 'Invalid parameter: dcaps_prod_min/dcaps_prod_max'),
    ('dcaps_delta: требуется целое число от 1 до 1000', 'dcaps_delta: an integer from 1 to 1000 is required'),
    ('primer_pairs: требуется целое число от 1 до 100', 'primer_pairs: an integer from 1 to 100 is required'),
    ('dcaps: допустимо off, fallback или always', 'dcaps: allowed values are off, fallback or always'),
    ('Параметры модели полос должны быть конечными числами', 'Band model parameters must be finite numbers'),
    ('Некорректный диапазон видимых фрагментов', 'Invalid range of visible fragments'),
    ('Некорректное разрешение геля', 'Invalid gel resolution'),
    ('Коэффициент ухудшения разрешения должен быть от 1 до 5',
     'The resolution stress factor must be between 1 and 5'),
    ('Длины фрагментов должны быть положительными', 'Fragment lengths must be positive'),
    ('[INFO] Поиск: joint; SNP: 6; ферментов: 623', '[INFO] Search: joint; SNPs: 6; enzymes: 623'),
    ('[INFO] Оценка полос: модель присутствия, полная рестрикция; '
     'яркость и специфичность ПЦР по всему геному не оцениваются.',
     '[INFO] Band model: presence of bands, complete digestion; '
     'band brightness and genome-wide PCR specificity are not assessed.'),
    ('[INFO] Отбор: только продаваемые ферменты; изошизомеры объединены (одна строка на сайт); '
     'без сайта фермента в праймерах; повторы (строчные буквы FASTA) не на 3′-концах праймеров; '
     'без фоновых вариантов, меняющих разрезы; SNP не ближе 100 п.н. к концам ампликона; '
     'dCAPS если нет robust-схемы, несовпадений до 1, Tm до 70 °C, продукт 120–300 п.н., Δ 20 п.н.',
     '[INFO] Filters: commercial enzymes only; isoschizomers grouped (one row per site); '
     'no enzyme site in primers; repeats (lower-case FASTA) not at primer 3′ ends; '
     'no background variants that change cuts; SNP at least 100 bp from amplicon ends; '
     'dCAPS if no robust assay, mismatches up to 1, Tm up to 70 °C, product 120–300 bp, Δ 20 bp'),
    ('dCAPS для каждого SNP, несовпадений до 2', 'dCAPS for every SNP, mismatches up to 2'),
    ('[INFO] Фоновые варианты: 12 в окнах SNP из 15000 записей VCF',
     '[INFO] Background variants: 12 in SNP windows out of 15000 VCF records'),
    ('[INFO] 11:37525128 C>G: пар 10, кандидатов 1, неразличимых 0',
     '[INFO] 11:37525128 C>G: pairs 10, candidates 1, indistinguishable 0'),
    ('[WARN] В окнах SNP нет строчных букв: FASTA без soft-mask, маскирование повторов не действует.',
     '[WARN] No lower-case letters in SNP windows: the FASTA is not soft-masked, repeat masking has no effect.'),
    ('\n[INFO] NORMALIZE: загружено вариантов: 6', '\n[INFO] NORMALIZE: variants loaded: 6'),
    # Logger lines (shown by the command line).
    (r'Файл names не найден: D:\names.txt', r'Names file not found: D:\names.txt'),
    ('Ошибка при загрузке names: bad line', 'Error loading the names file: bad line'),
    ('Файл с вариантами не найден: snps.tsv', 'Variants file not found: snps.tsv'),
    ('Ошибка парсинга строки: bad value', 'Line parse error: bad value'),
    ('Ошибка при чтении файла с вариантами: denied', 'Error reading the variants file: denied'),
]


def translatable_literals() -> set:
    """First literal arguments of t(...) and N_(...) calls in the program's modules."""
    found = set()
    for path in ROOT.glob('*.py'):
        for node in ast.walk(ast.parse(path.read_text(encoding='utf-8'))):
            if not (isinstance(node, ast.Call) and node.args):
                continue
            func = node.func
            name = func.id if isinstance(func, ast.Name) else func.attr if isinstance(func, ast.Attribute) else ''
            first = node.args[0]
            if name in ('t', 'N_') and isinstance(first, ast.Constant) and isinstance(first.value, str):
                found.add(first.value)
    return found


def placeholders(text: str) -> set:
    return {field for _, field, _, _ in string.Formatter().parse(text) if field is not None}


class I18nTests(unittest.TestCase):
    def setUp(self):
        self.addCleanup(i18n.set_language, i18n.language())

    def test_t_returns_source_in_russian_and_translation_in_english(self):
        with patch.dict(translations_en.EN, {'Показано {shown} из {total}': 'Shown {shown} of {total}'}):
            self.assertEqual(i18n.t('Показано {shown} из {total}', shown=2, total=4), 'Показано 2 из 4')
            i18n.set_language('en')
            self.assertEqual(i18n.t('Показано {shown} из {total}', shown=2, total=4), 'Shown 2 of 4')
            self.assertEqual(i18n.t('Нет такой строки'), 'Нет такой строки')

    def test_text_without_values_is_not_formatted(self):
        self.assertEqual(i18n.t('Скобки {не подстановка}'), 'Скобки {не подстановка}')

    def test_marker_returns_the_source(self):
        i18n.set_language('en')
        self.assertEqual(i18n.N_('Надёжно'), 'Надёжно')

    def test_unknown_language_is_rejected(self):
        with self.assertRaises(ValueError):
            i18n.set_language('fr')

    def test_fmt(self):                                   # Review Focus 5
        self.assertEqual(i18n.fmt(16.13), '16,1')
        self.assertEqual((i18n.fmt(475), i18n.fmt(None), i18n.fmt('')), ('475', '', ''))
        self.assertEqual(i18n.fmt(59.04, 2), '59,04')
        i18n.set_language('en')
        self.assertEqual((i18n.fmt(16.13), i18n.fmt(59.04, 1)), ('16.1', '59.0'))

    def test_language_detection(self):
        self.assertEqual(i18n.system_language(['ru-RU', 'en-US']), 'ru')
        self.assertEqual(i18n.system_language(['en-US', 'ru-RU']), 'en')
        self.assertEqual(i18n.system_language(['de-DE']), 'en')
        self.assertEqual(i18n.system_language([]), 'en')
        self.assertEqual(i18n.locale_language('Russian_Russia'), 'ru')
        self.assertEqual(i18n.locale_language('ru_RU'), 'ru')
        self.assertEqual(i18n.locale_language('en_US'), 'en')
        self.assertEqual(i18n.locale_language(None), 'en')
        self.assertEqual(i18n.resolve_language('en', 'ru'), 'en')
        self.assertEqual(i18n.resolve_language('ru', 'en'), 'ru')
        self.assertEqual(i18n.resolve_language('system', 'ru'), 'ru')
        self.assertEqual(i18n.resolve_language('fr', 'en'), 'en')
        self.assertEqual(i18n.resolve_language(None, 'ru'), 'ru')
        self.assertEqual(i18n.resolve_language(3, 'en'), 'en')

    def test_core_texts_are_translated(self):
        i18n.set_language('en')
        for russian, english in CORE_SAMPLES:
            self.assertEqual(i18n.core_text(russian), english, russian)

    def test_compound_reason_is_translated_whole(self):   # Review Focus 1
        i18n.set_language('en')
        text = i18n.core_text('Оценка окна без подобранных праймеров. dCAPS: несовпадение в праймере создаёт '
                              'сайт в аллеле REF. Три генотипа различимы только при исходном разрешении')
        self.assertIsNone(CYRILLIC.search(text), text)
        text = i18n.core_text('Не найдено сочетание для различения трёх генотипов; '
                              'окон с неопределёнными основаниями: 2; ошибок расчёта ферментов: 1')
        self.assertEqual(text, 'No primer–enzyme combination distinguishes the three genotypes; '
                               'windows with ambiguous bases: 2; enzyme calculation errors: 1')

    def test_core_text_keeps_russian_and_unknown_text(self):
        self.assertEqual(i18n.core_text('REF не совпадает с FASTA'), 'REF не совпадает с FASTA')
        i18n.set_language('en')
        self.assertEqual(i18n.core_text('Совсем новое сообщение'), 'Совсем новое сообщение')
        self.assertEqual(i18n.core_text(''), '')

    def test_english_plural_follows_count(self):
        i18n.set_language('en')
        self.assertEqual(i18n.t('Готово: {count} SNP · {parts}', count=1, parts='robust 1'), 'Done: 1 SNP · robust 1')
        self.assertEqual(i18n.t('Готово: {count} SNP · {parts}', count=4, parts='robust 4'), 'Done: 4 SNPs · robust 4')
        self.assertEqual(i18n.t('{done} из {count} SNP', done=1, count=1), '1 of 1 SNP')
        self.assertEqual(i18n.t('{done} из {count} SNP', done=1, count=3), '1 of 3 SNPs')
        self.assertEqual(i18n.t('Идёт расчёт… {done} из {count} SNP', done=0, count=1), 'Running… 0 of 1 SNP')
        self.assertEqual(i18n.t('Ручной выбор у {count} SNP будет сброшен. Запустить?', count=1),
                         'The choice made by hand for 1 SNP will be reset. Run?')

    def test_english_uses_english_quotation_marks(self):
        texts = [x for value in translations_en.EN.values() for x in (value if isinstance(value, tuple) else [value])]
        texts += [replacement for _, replacement in translations_en.CORE]
        self.assertEqual([x for x in texts if '«' in x or '»' in x], [])

    def test_catalog_is_complete(self):                   # Review Focus 4
        used = translatable_literals()
        self.assertEqual(sorted(used - set(translations_en.EN)), [])
        self.assertEqual(sorted(set(translations_en.EN) - used), [])
        for source, target in translations_en.EN.items():
            for form in target if isinstance(target, tuple) else (target,):
                self.assertEqual(placeholders(source), placeholders(form), source)
                self.assertIsNone(CYRILLIC.search(form), source)
            if isinstance(target, tuple):        # (one, other), chosen by the count value
                self.assertEqual(len(target), 2, source)
                self.assertIn('count', placeholders(source), source)


if __name__ == '__main__':
    unittest.main()
