import contextlib
import csv
import io
import random
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

from PyQt6.QtCore import QSettings, Qt
from PyQt6.QtWidgets import QApplication
from Bio.Restriction import BspLU11I, PciI

import i18n
import rflp_cli
import rflp_core as core
from app_settings import AppSettings
from launch_panel import InventoryDialog, LaunchPanel
from qtprimer3_vis import RFLPPickerGUI
from tests.helpers import dispose
from theme import ThemeController


class LabEnzymeTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.genome = self.root / 'sample.fa'
        self.genome.write_text('>1\n' + 'C' * 147 + 'ACATGT' + 'C' * 147 + '\n')
        self.variants = self.root / 'variants.tsv'
        self.variants.write_text('chrom\tpos\tref\talt\n1\t150\tA\tG\n')
        self.inventory = self.root / 'enzymes.txt'
        self.inventory.write_text('bsplu11i\nBspLU11I\n', encoding='utf-8-sig')
        self.addCleanup(i18n.set_language, i18n.language())

    def test_cli_modes_before_grouping_and_limit(self):
        for mode in ('only', 'prefer'):
            out = self.root / 'out.csv'
            args = ['--lang', 'en', '--fasta', str(self.genome), '--variants', str(self.variants),
                    '--out', str(out), '--enzyme-list', str(self.inventory), '--enzyme-mode', mode,
                    '--no-primer3', '--min-frag', '20', '--top', '1', '-q']
            self.assertEqual(rflp_cli.main(args), 0)
            with out.open(encoding='utf-8-sig') as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual([row['enzyme'] for row in rows], ['BspLU11I'])

    def test_quality_precedes_availability_within_group(self):
        def hit(enzyme, quality, score):
            return dict(_enz=enzyme, enzyme=str(enzyme), genotype_quality=quality,
                        genotype_margin=score, worst_margin=score, visible_band_count=3)
        rows = [hit(PciI, 'robust', 3), hit(BspLU11I, 'nominal_only', 0.8)]
        self.assertEqual(core.group_by_site(rows, {'BspLU11I'})[0]['enzyme'], 'PciI')

    def test_cli_and_gui_rank_quality_before_availability(self):
        # This window has both robust designs and a nominal-only SetI design.
        rng = random.Random(2)
        def flank():
            return ''.join(rng.choice('ACGT') for _ in range(147))
        self.genome.write_text('>1\n' + flank() + 'ACATGT' + flank() + '\n')
        self.inventory.write_text('SetI')
        out = self.root / 'result.csv'
        args = ['--fasta', str(self.genome), '--variants', str(self.variants), '--out', str(out),
                '--enzyme-list', str(self.inventory), '--no-primer3', '--min-frag', '10', '--delta', '1', '-q']
        self.assertEqual(rflp_cli.main([*args, '--top', '0']), 0)
        with out.open(encoding='utf-8-sig') as stream:
            records = list(csv.DictReader(stream))
        self.assertIn(('SetI', 'nominal_only'), [(r['enzyme'], r['genotype_quality']) for r in records])
        self.assertEqual(records[0]['genotype_quality'], 'robust')
        app = QApplication.instance() or QApplication([])
        settings = AppSettings(QSettings(str(self.root / 'window.ini'), QSettings.Format.IniFormat))
        window = RFLPPickerGUI(settings, ThemeController(app, 'light'))
        self.addCleanup(dispose, window)
        window.launch.set_files(str(self.variants), str(self.genome))
        for mode, names, expected in [('prefer', ['SetI'], records[0]['enzyme']),
                                      ('only', ['BspLU11I'], 'BspLU11I'),
                                      ('prefer', ['BspLU11I'], 'BspLU11I')]:
            window.launch.set_params(dict(enzyme_mode=mode, lab_enzymes=names, use_primer3=False,
                                          min_frag=10, delta=1, top_results=1))
            window.start_run()
            self.assertTrue(window.worker.wait(10000))
            app.processEvents()
            self.assertEqual([r[window.header.index('enzyme')] for r in window.rows], [expected])

    def test_bad_file_stops_cli_before_analysis(self):
        self.inventory.write_text('Typo1\nTypo2')
        errors = io.StringIO()
        with patch.object(core, 'run_rflp_gui_mode') as run, contextlib.redirect_stderr(errors):
            code = rflp_cli.main(['--lang', 'en', '--fasta', str(self.genome), '--variants', str(self.variants),
                                  '--out', str(self.root / 'out.csv'), '--enzyme-list', str(self.inventory)])
        run.assert_not_called()
        self.assertEqual(code, 2)
        self.assertIn('Unknown enzymes: Typo1, Typo2', errors.getvalue())

    def test_invalid_import_reports_all_unknown_names(self):
        self.inventory.write_text('EcoRI\nTypo1\nTypo2')
        with self.assertRaisesRegex(ValueError, 'Typo1.*Typo2'):
            core.read_enzyme_list(self.inventory)

    def test_gui_import_and_manual_selection_in_both_languages(self):
        app = QApplication.instance() or QApplication([])
        settings = AppSettings(QSettings(str(self.root / "settings.ini"), QSettings.Format.IniFormat))
        for lang in ('ru', 'en'):
            i18n.set_language(lang)
            panel = LaunchPanel(settings)
            dialog = InventoryDialog([], panel)
            with patch('launch_panel.QFileDialog.getOpenFileName', return_value=(str(self.inventory), '')):
                dialog.load_file()
            self.assertEqual(dialog.selected_names(), ['BspLU11I'])
            for i in range(dialog.list_widget.count()):
                item = dialog.list_widget.item(i)
                if item.text() == 'EcoRI':
                    item.setCheckState(Qt.CheckState.Checked)
                    break
            panel.set_params(dict(enzyme_mode='only', lab_enzymes=dialog.selected_names()))
            self.assertEqual(panel.collect_params()['lab_enzymes'], ['BspLU11I', 'EcoRI'])
            self.assertEqual(panel.collect_params()['enzyme_mode'], 'only')
            self.assertFalse(panel.supplier_picker.isEnabled())
            self.assertEqual(panel.inventory_mode.currentText(), 'Только из списка' if lang == 'ru' else 'Listed only')
            panel.reset_to_defaults()
            self.assertTrue(panel.supplier_picker.isEnabled())
            self.assertEqual(panel.collect_params()['lab_enzymes'], [])
            dispose(panel)
            app.processEvents()
