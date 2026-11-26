import json
import os
import random
import tempfile
import unittest
from pathlib import Path

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
from PyQt6.QtCore import QSettings
from PyQt6.QtWidgets import QApplication

import i18n
import lab_sheet
import qtprimer3_vis
import rflp_core
from app_settings import AppSettings
from details_panel import DetailsPanel
from snp_summary import InputRecord, SnpSummary
from tests.helpers import CHECKS_OFF, dispose, write_text
from tests.test_lab_sheet import page_count
from theme import ThemeController


def _summary():
    control = dict(product_size=330, primers_left='TCTTAGGG', tm_left=59.0,
                   primers_right='GACACCTC', tm_right=59.1)
    results = dict(
        aspcr_specific_side='left',
        aspcr_control=control,
        aspcr_reactions={
            'REF': dict(allele='C', product_size=129, primers_left='GACCTAAC', tm_left=58.0),
            'ALT': dict(allele='G', product_size=130, primers_left='GGGTTCAG', tm_left=58.1),
        },
        aspcr_genotypes={
            'REF/REF': {'REF': [129, 330], 'ALT': [330]},
            'REF/ALT': {'REF': [129, 330], 'ALT': [130, 330]},
            'ALT/ALT': {'REF': [330], 'ALT': [130, 330]},
        },
    )
    for key in ('aspcr_control', 'aspcr_reactions', 'aspcr_genotypes'):
        results[key] = json.dumps(results[key])
    scheme = dict(results, design='allele_specific', genotype_quality='robust')
    return SnpSummary(InputRecord(2, 'target', '1', 42, 'C', 'G', ''), [scheme], primer_base='target')


class AlleleSpecificUiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self):
        i18n.set_language('ru')
        self.addCleanup(i18n.set_language, 'ru')

    def test_details_show_two_tubes_and_failed_control_rule(self):
        panel = DetailsPanel()
        self.addCleanup(dispose, panel)
        panel.show_summary(_summary(), 'off')
        self.assertEqual(panel.tabs.tabText(0), 'Реакции')
        self.assertEqual(panel.tabs.tabText(1), 'Гель')
        self.assertEqual(panel.gel_message(), '')
        text = panel.reaction_details.text()
        for part in ('Пробирка REF · аллель C', 'Пробирка ALT · аллель G', 'Контроль 330 п.н.',
                     'REF/REF', 'REF/ALT', 'ALT/ALT', '129 + 330', '130 + 330', 'генотип не определён'):
            self.assertIn(part, text)
        for sequence in ('GACCTAAC', 'GGGTTCAG', 'GACACCTC', 'TCTTAGGG'):
            self.assertIn(sequence, panel.card_details.text())

    def test_mixed_pdf_contains_one_page_per_assay_with_four_aspcr_oligos(self):
        aspcr = _summary()
        restriction = SnpSummary(InputRecord(3, 'second', '1', 202, 'C', 'T', ''),
                                 [dict(enzyme='EcoRI', design='natural', genotype_quality='robust')])
        self.assertEqual([(name, seq) for name, seq, _, _ in lab_sheet.primer_rows(aspcr)],
                         [('target_REF', 'GACCTAAC'), ('target_ALT', 'GGGTTCAG'),
                          ('target_C', 'GACACCTC'), ('target_IC', 'TCTTAGGG')])
        with tempfile.TemporaryDirectory() as folder:
            path = str(Path(folder) / 'assays.pdf')
            self.assertEqual(lab_sheet.export_lab_sheets(path, [restriction, aspcr], '', {}, 'off'), 2)
            self.assertEqual(page_count(path), 2)

    @unittest.skipUnless(rflp_core._PRIMER3_AVAILABLE, 'Primer3 unavailable')
    def test_project_keeps_restriction_and_aspcr_gels_without_fasta(self):
        folder = tempfile.TemporaryDirectory()
        self.addCleanup(folder.cleanup)
        origin, restored = self._project_windows(Path(folder.name))
        rng = random.Random(19)
        sequence = ''.join(rng.choice('ACGT') for _ in range(500))
        sequence = sequence[:197] + 'GAATTC' + sequence[203:]
        fasta = write_text(folder.name, 'mixed.fa', '>1\n' + sequence + '\n')
        variants = write_text(folder.name, 'targets.tsv',
                              'chrom\tpos\tref\talt\tid\n1\t200\tA\tG\tEcoRI_test\n1\t225\tT\tA\tAS_test\n')
        params = dict(rflp_core.DEFAULT_PARAMS, **CHECKS_OFF)
        params.update(fasta=fasta, input_file=variants, flank=500, prod_min=180, prod_max=450,
                      tm_min=50, tm_max=70, primer_pairs=4, max_cuts=1, gain_loss_only=True,
                      allele_specific=True, enzyme_mode='only', lab_enzymes=['EcoRI'])
        report = {}
        header, rows = rflp_core.run_rflp_gui_mode(params, report=report)
        origin.on_succeeded(header, rows, report, params)
        self.assertEqual([summary.scheme['design'] for summary in origin.summaries], ['natural', 'allele_specific'])
        origin.summaries[1].choose(1)
        project = str(Path(folder.name) / 'mixed.rflp')
        self.assertEqual(origin.save_project(project), project)
        origin.details.close_fasta()
        Path(fasta).unlink()
        Path(variants).unlink()

        self.assertEqual(restored.open_project(project), project)
        self.assertEqual(restored.summaries[1].chosen, 1)
        restored.details.show_summary(restored.summaries[0], 'off')
        self.assertTrue(all(restored.details.map.allele_sequences()))
        restored.details.show_summary(restored.summaries[1], 'off')
        self.assertEqual(restored.details.gel_message(), '')
        self.assertIn('Пробирка ALT · аллель A', restored.details.reaction_details.text())
        pdf = str(Path(folder.name) / 'mixed.pdf')
        self.assertEqual(restored.export('pdf', pdf), pdf)
        self.assertEqual(page_count(pdf), 2)

    def _project_windows(self, folder):
        windows = []
        for name in ('origin', 'restored'):
            qsettings = QSettings(str(folder / f'{name}.ini'), QSettings.Format.IniFormat)
            window = qtprimer3_vis.RFLPPickerGUI(AppSettings(qsettings), ThemeController(self.app, 'light'))
            self.addCleanup(dispose, window)
            windows.append(window)
        return windows
