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

    def test_catalog_is_complete(self):                   # Review Focus 4
        used = translatable_literals()
        self.assertEqual(sorted(used - set(translations_en.EN)), [])
        self.assertEqual(sorted(set(translations_en.EN) - used), [])
        for source, target in translations_en.EN.items():
            self.assertEqual(placeholders(source), placeholders(target), source)
            self.assertIsNone(CYRILLIC.search(target), source)


if __name__ == '__main__':
    unittest.main()
