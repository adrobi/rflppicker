import json
import os
import random
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
import openpyxl
from PyQt6.QtCore import QSettings
from PyQt6.QtWidgets import QApplication, QMessageBox

import lab_sheet
import qtprimer3_vis
import rflp_core
from app_settings import AppSettings
from tests.helpers import CHECKS_OFF, dispose, write_text
from tests.test_lab_sheet import page_count
from theme import ThemeController


class ProjectFileTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.folder = Path(directory.name)
        self.window = self.make_window('origin')
        self.project = self.folder / 'transfer.rflp'

    def make_window(self, name):
        settings = AppSettings(QSettings(str(self.folder / f'{name}.ini'), QSettings.Format.IniFormat))
        window = qtprimer3_vis.RFLPPickerGUI(settings, ThemeController(self.app, 'light'))
        self.addCleanup(dispose, window)
        return window

    def analyze(self):
        rng = random.Random(19)
        sequence = ''.join(rng.choice('ACGT') for _ in range(400))
        sequence = sequence[:197] + 'GAATTC' + sequence[203:]
        fasta = write_text(self.folder, 'GCF_012345678.1_test.fa', '>1\n' + sequence + '\n')
        variants = write_text(self.folder, 'targets.tsv', 'chrom\tpos\tref\talt\tid\n1\t200\tA\tG\tEcoRI_test\n'
                              '1\t100\tN\tG\tbad\n')
        params = dict(rflp_core.DEFAULT_PARAMS, **CHECKS_OFF, fasta=fasta, input_file=variants, flank=500,
                      prod_min=200, prod_max=500, tm_min=50, tm_max=70, primer_pairs=4, max_cuts=1,
                      gain_loss_only=True)
        report = {}
        header, rows = rflp_core.run_rflp_gui_mode(params, report=report)
        self.window.on_succeeded(header, rows, report, params)
        self.window.list.table.selectRow(0)
        return fasta, variants

    def test_project_restores_choices_gel_and_exports_without_sources(self):
        fasta, variants = self.analyze()
        summary = self.window.summaries[0]
        eco = next(i for i, scheme in enumerate(summary.schemes) if i and scheme['enzyme'] == 'EcoRI')
        self.window.details.choose_scheme(eco)
        picker = self.window.details.product_picker
        picker.setCurrentIndex(picker.findData('R3101'))
        expected_sequences = self.window.details.map.allele_sequences()
        self.assertTrue(all(expected_sequences))
        self.assertEqual(self.window.save_project(str(self.project)), str(self.project))
        self.window.details.close_fasta()
        Path(fasta).unlink()
        Path(variants).unlink()

        restored = self.make_window('destination')
        with patch('pyfaidx.Fasta', side_effect=AssertionError('project unexpectedly accessed FASTA')), \
                patch('project_file.Fasta', side_effect=AssertionError('project unexpectedly accessed FASTA')):
            self.assertEqual(restored.open_project(str(self.project)), str(self.project))
            restored.list.table.selectRow(0)
            self.assertEqual(restored.summaries[0].chosen, eco)
            self.assertEqual(restored.details.product_picker.currentData(), 'R3101')
            self.assertEqual(restored.details.map.allele_sequences(), expected_sequences)
            self.assertEqual(restored.details.gel_message(), '')
            restored.details.view_scheme(0)
            self.assertEqual(restored.details.gel_message(), '')
            self.assertEqual(restored.launch.collect_params()['prod_min'], 200)
            self.assertEqual(restored.launch.collect_params()['fasta'], '')
            self.assertEqual(restored.launch.collect_params()['input_file'], '')
            pdf = str(self.folder / 'reopened.pdf')
            with patch('lab_sheet.draw_sheet', wraps=lab_sheet.draw_sheet) as draw:
                self.assertEqual(restored.export('pdf', pdf), pdf)
                args = draw.call_args.args
                self.assertEqual(args[1].scheme['enzyme_product'], 'R3101')
                self.assertEqual(args[2]['amplicon_seq'], expected_sequences[0])
            self.assertEqual(page_count(pdf), 1)
            copy = str(self.folder / 'copy.rflp')
            self.assertEqual(restored.save_project(copy), copy)

        workbook = openpyxl.load_workbook(restored.export('all', str(self.folder / 'reopened.xlsx')))
        self.addCleanup(workbook.close)
        metadata = {row[0].value: row[1].value for row in workbook['Run'] if len(row) >= 2}
        self.assertEqual(metadata['assembly_accession'], 'GCF_012345678.1')
        self.assertEqual(metadata['scoring_model'], self.window.report['scoring_model'])
        self.assertTrue(metadata['enzyme_catalog_version'])
        self.assertEqual(restored.summaries[1].status, 'error')

    def test_invalid_project_does_not_replace_current_results(self):
        self.analyze()
        self.project.write_text(json.dumps({'format': 'rflp-picker-project', 'version': 999}), encoding='utf-8')
        with patch.object(QMessageBox, 'warning') as warning:
            self.assertIsNone(self.window.open_project(str(self.project)))
        warning.assert_called_once()
        self.assertEqual(len(self.window.summaries), 2)
        self.assertEqual(self.window.summaries[0].record.snp_id, 'EcoRI_test')

    def test_incomplete_render_payload_is_rejected_before_results_change(self):
        self.analyze()
        self.assertEqual(self.window.save_project(str(self.project)), str(self.project))
        document = json.loads(self.project.read_text(encoding='utf-8'))
        document['payloads'][str(self.window.summaries[0].scheme['_index'])].pop('amplicon_seq')
        document['records'][0]['snp_id'] = 'corrupt'
        self.project.write_text(json.dumps(document), encoding='utf-8')
        with patch.object(QMessageBox, 'warning') as warning:
            self.assertIsNone(self.window.open_project(str(self.project)))
        warning.assert_called_once()
        self.assertEqual(self.window.summaries[0].record.snp_id, 'EcoRI_test')

    def test_invalid_cached_alleles_and_nan_are_rejected(self):
        self.analyze()
        self.assertEqual(self.window.save_project(str(self.project)), str(self.project))
        original = self.project.read_text(encoding='utf-8')
        key = str(self.window.summaries[0].scheme['_index'])
        for field, value in (('locus_products', [{'allele': 'REF', 'sequence': 'A'}]),
                             ('snp_offset', float('nan'))):
            with self.subTest(field=field), patch.object(QMessageBox, 'warning') as warning:
                document = json.loads(original)
                document['payloads'][key][field] = value
                self.project.write_text(json.dumps(document), encoding='utf-8')
                self.assertIsNone(self.window.open_project(str(self.project)))
                warning.assert_called_once()

    def test_bad_launch_parameters_leave_current_run_and_inputs_intact(self):
        fasta, _ = self.analyze()
        launch = self.window.launch
        targets = write_text(self.folder, 'pending.vcf', '##fileformat=VCFv4.2\n')
        background = write_text(self.folder, 'background.vcf', '##fileformat=VCFv4.2\n')
        launch.set_files(targets, fasta, background)
        launch.set_params({'assembly_accession': 'current', 'vcf_region': '1:10-30',
                           'vcf_samples': ['sampleA'], 'prod_min': 220})
        previous_params = launch.collect_params()
        previous_summaries = self.window.summaries
        document = {'format': 'rflp-picker-project', 'version': 1, 'metadata': {}, 'has_ids': False,
                    'results': {'header': [], 'rows': [], 'report': {}},
                    'records': [], 'choices': [], 'products': {}, 'payloads': {}}
        for bad_parameter in ({'vcf_samples': 42}, {'prod_min': 2**80}, {'specificity': 'unknown'},
                              {'gel_digest_pct': -10}, {'vcf_region': '1:0-2'}):
            with self.subTest(bad_parameter=bad_parameter):
                document['parameters'] = {'assembly_accession': 'overwritten', **bad_parameter}
                self.project.write_text(json.dumps(document), encoding='utf-8')
                with patch.object(QMessageBox, 'warning') as warning:
                    self.assertIsNone(self.window.open_project(str(self.project)))
                warning.assert_called_once()
                self.assertEqual(launch.collect_params(), previous_params)
                self.assertIs(self.window.summaries, previous_summaries)
