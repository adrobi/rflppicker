import os
import tempfile
import unittest

from PyQt6.QtCore import QSettings

import rflp_core
from app_settings import FILE_KEYS, AppSettings

RECOMMENDED = {k: v for k, v in rflp_core.RECOMMENDED_PARAMS.items() if k not in FILE_KEYS}


class AppSettingsTests(unittest.TestCase):
    def setUp(self):
        self._dir = tempfile.TemporaryDirectory()
        self.path = os.path.join(self._dir.name, 's.ini')
        self.qs = QSettings(self.path, QSettings.Format.IniFormat)
        self.settings = AppSettings(self.qs)

    def tearDown(self):
        self.qs.sync()
        self._dir.cleanup()

    def reopen(self):
        self.qs.sync()
        return QSettings(self.path, QSettings.Format.IniFormat)

    def test_defaults_are_recommended(self):
        self.assertEqual(self.settings.params(), RECOMMENDED)
        self.assertNotIn('background_vcf', self.settings.params())

    def test_params_roundtrip_and_bad_values(self):
        self.settings.set_params(dict(self.settings.params(), prod_min=150, dcaps='always', suppliers=['E']))
        params = AppSettings(self.reopen()).params()
        self.assertEqual((params['prod_min'], params['dcaps'], params['suppliers']), (150, 'always', ['E']))
        self.qs.setValue('params', '{"prod_min": "много", "unknown": 1, "tm_max": 61.5, "tm_min": 57,'
                                   ' "mask_repeats": 1, "group_isoschizomers": false}')
        p = AppSettings(self.reopen()).params()
        self.assertEqual((p['prod_min'], p['tm_max'], p['tm_min']), (250, 61.5, 57.0))
        self.assertIsInstance(p['tm_min'], float)
        self.assertEqual((p['mask_repeats'], p['group_isoschizomers']), (True, False))
        self.assertNotIn('unknown', p)
        self.qs.setValue('params', 'not json')
        self.assertEqual(AppSettings(self.reopen()).params(), RECOMMENDED)

    def test_set_params_ignores_file_keys(self):
        self.settings.set_params(dict(RECOMMENDED, background_vcf='x.vcf', fasta='g.fa'))
        self.assertEqual(self.settings.params(), RECOMMENDED)

    def test_recent_genomes(self):
        for path in ['a.fa', 'b.fa', 'A.FA', 'c.fa', 'd.fa', 'e.fa', 'f.fa']:
            self.settings.add_recent_genome(path)
        self.assertEqual(self.settings.recent_genomes(), ['f.fa', 'e.fa', 'd.fa', 'c.fa', 'A.FA'])
        self.assertEqual(AppSettings(self.reopen()).recent_genomes()[0], 'f.fa')

    def test_names_per_genome(self):
        self.assertIsNone(self.settings.names_for_genome(r'D:\g.fa'))
        self.settings.set_names_for_genome(r'D:\g.fa', '')
        self.assertEqual(self.settings.names_for_genome(r'D:\g.fa'), '')
        self.settings.set_names_for_genome(r'D:\g.fa', r'D:\names.txt')
        self.assertEqual(AppSettings(self.reopen()).names_for_genome(r'D:\g.fa'), r'D:\names.txt')
        self.settings.clear_names_for_genome(r'D:\g.fa')
        self.assertIsNone(self.settings.names_for_genome(r'D:\g.fa'))

    @unittest.skipUnless(os.name == 'nt', 'Windows paths are case-insensitive; elsewhere they are not')
    def test_names_per_genome_ignore_path_case_on_windows(self):
        self.settings.set_names_for_genome(r'D:\g.fa', r'D:\names.txt')
        self.assertEqual(self.settings.names_for_genome(r'd:\G.FA'), r'D:\names.txt')

    def test_values_keep_their_types(self):
        self.settings.set_value('launch_collapsed', True)
        self.settings.set_value('splitter_sizes', [300, 700, 500])
        self.settings.set_value('theme_mode', 'dark')
        reopened = AppSettings(self.reopen())
        self.assertIs(reopened.value('launch_collapsed'), True)
        self.assertEqual(reopened.value('splitter_sizes'), [300, 700, 500])
        self.assertEqual(reopened.value('theme_mode'), 'dark')
        self.assertEqual(reopened.value('missing', 'fallback'), 'fallback')


if __name__ == '__main__':
    unittest.main()
