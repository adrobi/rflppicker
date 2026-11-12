import importlib.util
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location(
    'release_notes', Path(__file__).resolve().parent.parent / 'tools' / 'release_notes.py')
release_notes = importlib.util.module_from_spec(spec)
spec.loader.exec_module(release_notes)

CHANGELOG = '''# Изменения

## Не выпущено

- ещё не готово

## 1.5.0

**Убрано**

- режим одной пары

## 1.4.0

- старое
'''


class ReleaseNotesTests(unittest.TestCase):
    def test_section_of_the_version(self):
        self.assertEqual(release_notes.section(CHANGELOG, '1.5.0'), '**Убрано**\n\n- режим одной пары')
        self.assertEqual(release_notes.section(CHANGELOG, '1.4.0'), '- старое')

    def test_missing_version_is_an_error(self):
        for version in ('1.5', '2.0.0', 'Не выпущено'):
            with self.assertRaises(SystemExit):
                release_notes.section(CHANGELOG, version)

    def test_notes_list_the_archives_with_checksums(self):
        text = release_notes.notes(CHANGELOG, '1.5.0', {'a-windows-x64.zip': 'ab12', 'a-linux-x86_64.tar.gz': 'cd34'})
        self.assertTrue(text.startswith('**Убрано**'))
        self.assertIn('`a-windows-x64.zip` — SHA-256 `ab12`', text)
        self.assertIn('`a-linux-x86_64.tar.gz` — SHA-256 `cd34`', text)
        self.assertIn('https://github.com/adrobi/rflppicker/blob/v1.5.0/README.md', text)

    def test_tag_must_match_the_program_version(self):
        import rflp_core
        self.assertEqual(release_notes.version_of_tag('v' + rflp_core.APP_VERSION), rflp_core.APP_VERSION)
        with self.assertRaises(SystemExit):
            release_notes.version_of_tag('v0.0.1')
