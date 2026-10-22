import os
import tempfile
import unittest
from unittest.mock import patch

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
from PyQt6.QtCore import QSettings
from PyQt6.QtWidgets import QApplication, QComboBox, QLabel, QPushButton

import rflp_core
from app_settings import AppSettings
from launch_panel import ADVANCED_PARAMS, KEY_PARAMS, LaunchPanel, SupplierPicker, field_errors
from tests.helpers import write_text


class LaunchPanelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self):
        self._dir = tempfile.TemporaryDirectory()
        self.tmp = self._dir.name
        os.makedirs(os.path.join(self.tmp, 'names'))
        self.names_name = 'sheep.txt'
        write_text(os.path.join(self.tmp, 'names'), self.names_name, '1\tNC_1\n')
        self.genome = write_text(self.tmp, 'sheep.fa', '>NC_1\nACGT\n')
        self.ini = os.path.join(self.tmp, 's.ini')
        self.qs = QSettings(self.ini, QSettings.Format.IniFormat)
        self.settings = AppSettings(self.qs)
        patcher = patch('launch_panel.get_base_dir', return_value=self.tmp)
        patcher.start()
        self.addCleanup(patcher.stop)
        self.panel = LaunchPanel(self.settings)

    def tearDown(self):
        self.panel.close()
        self.qs.sync()
        self._dir.cleanup()

    def reopen(self):
        self.qs.sync()
        return QSettings(self.ini, QSettings.Format.IniFormat)

    def test_parameter_groups_cover_the_core(self):
        files = {'input_file', 'fasta', 'names2_file', 'background_vcf'}
        self.assertEqual(set(KEY_PARAMS) | set(ADVANCED_PARAMS) | files,
                         set(rflp_core.RECOMMENDED_PARAMS) | {'input_file', 'fasta'})

    def test_defaults_are_recommended(self):
        p = self.panel.collect_params()
        for key, value in rflp_core.RECOMMENDED_PARAMS.items():
            if key not in ('names2_file', 'background_vcf'):
                self.assertEqual(p[key], value, key)
        self.assertEqual(self.panel.changed_count(), 0)
        self.assertTrue(self.panel.findChild(QLabel, 'changed_badge').isHidden())

    def test_supplier_picker_values_and_text(self):
        picker = self.panel.findChild(SupplierPicker, 'supplier_picker')
        self.assertEqual((picker.value(), picker.text()), ((True, []), 'Все продаваемые'))
        picker.set_value(True, ['E', 'I'])
        self.assertEqual(picker.value(), (True, ['E', 'I']))
        self.assertEqual(picker.text(), 'Thermo Fisher, SibEnzyme')
        picker.set_value(True, ['A', 'E', 'I', 'N'])
        self.assertEqual(picker.text(), 'Agilent, Thermo Fisher и ещё 2')
        picker.set_value(False, [])
        self.assertEqual(picker.text(), 'Все, включая непродаваемые')
        self.assertEqual(self.panel.collect_params()['commercial_only'], False)

    def test_changed_badge_and_presets(self):
        self.panel.set_params(dict(self.panel.collect_params(), delta=15, snp_end_distance=0))
        self.assertEqual(self.panel.changed_count(), 2)
        badge = self.panel.findChild(QLabel, 'changed_badge')
        self.assertEqual((badge.text(), badge.isHidden()), ('изменено 2', False))
        self.panel.apply_preset('v1.1')
        self.assertEqual(self.panel.collect_params()['dcaps'], 'off')
        self.panel.apply_preset('recommended')
        self.assertEqual(self.panel.changed_count(), 0)
        self.assertTrue(badge.isHidden())

    def test_presets_keep_the_files(self):
        self.panel.set_files(self.genome, self.genome, self.genome)
        self.panel.apply_preset('v1.1')
        p = self.panel.collect_params()
        self.assertEqual((p['input_file'], p['fasta'], p['background_vcf']), (self.genome,) * 3)

    def test_inline_errors(self):
        self.assertFalse(self.panel.validate())
        error = self.panel.findChild(QLabel, 'error_input_file')
        self.assertEqual((error.isHidden(), error.text()), (False, 'Выберите файл SNP'))
        self.assertEqual(self.panel.findChild(QLabel, 'error_fasta').text(), 'Выберите геном')
        self.panel.set_params(dict(self.panel.collect_params(), delta=0))
        self.panel.validate()
        self.assertEqual(self.panel.findChild(QLabel, 'changed_badge').text(), 'ошибка')
        self.panel.set_files(self.genome, self.genome)
        self.panel.set_params(dict(self.panel.collect_params(), delta=25))
        self.assertTrue(self.panel.validate())
        self.assertTrue(error.isHidden())

    def test_field_errors_mapping(self):
        params = dict(rflp_core.RECOMMENDED_PARAMS, input_file='', fasta=r'C:\missing.fa', prod_min=700,
                      prod_max=600, tm_min=63.0, primer_pairs=0)
        errors = field_errors(params)
        self.assertEqual(errors['input_file'], 'Выберите файл SNP')
        self.assertEqual(errors['fasta'], 'Файл не найден')
        self.assertEqual(errors['prod_min'], 'Минимальный размер продукта больше максимального')
        self.assertEqual(errors['tm_min'], 'Минимальный Tm больше максимального')
        self.assertTrue(errors['primer_pairs'].startswith('primer_pairs: требуется целое число'))

    def test_genome_list_and_names_label(self):
        self.settings.add_recent_genome(self.genome)
        panel = LaunchPanel(self.settings)
        self.addCleanup(panel.close)
        combo = panel.findChild(QComboBox, 'genome_combo')
        self.assertEqual(combo.itemText(0), 'sheep.fa')
        self.assertEqual(combo.itemData(0), self.genome)
        self.assertEqual(combo.itemText(combo.count() - 1), 'Выбрать файл…')
        self.assertEqual(panel.findChild(QLabel, 'names_label').text(), f'Хромосомы: {self.names_name} (авто)')
        self.assertEqual(panel.collect_params()['names2_file'], os.path.join(self.tmp, 'names', self.names_name))
        self.settings.set_names_for_genome(self.genome, '')
        panel.refresh_names()
        self.assertEqual(panel.findChild(QLabel, 'names_label').text(), 'Хромосомы: как в файле SNP')
        self.assertEqual(panel.collect_params()['names2_file'], '')

    def test_missing_remembered_file_is_marked(self):
        self.settings.set_value('last_snp_file', os.path.join(self.tmp, 'gone.tsv'))
        panel = LaunchPanel(self.settings)
        self.addCleanup(panel.close)
        self.assertEqual(panel.findChild(object, 'snp_edit').toolTip(), 'Файл не найден')

    def test_settings_roundtrip(self):
        self.panel.set_files(self.genome, self.genome)
        self.panel.set_params(dict(self.panel.collect_params(), prod_min=150, delta=15))
        self.panel.save_settings()
        panel = LaunchPanel(AppSettings(self.reopen()))
        self.addCleanup(panel.close)
        p = panel.collect_params()
        self.assertEqual((p['prod_min'], p['delta'], p['input_file'], p['fasta']),
                         (150, 15, self.genome, self.genome))

    def test_running_and_collapse(self):
        seen = []
        self.panel.stop_requested.connect(lambda: seen.append('stop'))
        self.panel.set_running(True)
        button = self.panel.findChild(QPushButton, 'run_button')
        self.assertEqual(button.text(), 'Остановить')
        self.assertFalse(self.panel.findChild(QComboBox, 'genome_combo').isEnabled())
        button.click()
        self.assertEqual(seen, ['stop'])
        self.panel.set_progress(37, 100)
        self.assertEqual(self.panel.findChild(QLabel, 'progress_label').text(), '37 из 100 SNP')
        self.panel.set_running(False)
        self.assertEqual(button.text(), 'Запустить')
        self.panel.set_collapsed(True)
        self.assertLessEqual(self.panel.maximumWidth(), 36)
        self.panel.set_collapsed(False)
        self.assertGreater(self.panel.maximumWidth(), 36)


if __name__ == '__main__':
    unittest.main()
