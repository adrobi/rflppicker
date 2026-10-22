import os
import tempfile
import unittest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
from PyQt6.QtWidgets import QApplication, QLabel, QListWidget, QPushButton, QTabWidget

import snp_summary
from details_panel import DetailsPanel
from snp_summary import InputRecord, order_lines
from tests.helpers import HEADER, core_row, write_text

SEQUENCE = 'AAAAAAAAAA' + 'GAATTC' + 'A' * 84 + 'GAATTC' + 'C' * 194
HOX = dict(variant='1:13 A>G', mapped_id='1', enzyme='EcoRI', site='GAATTC', isoschizomers='FunII; RI',
           genotype_quality='robust', genotype_margin=15.77, worst_margin=10.52, weakest_genotypes='REF/REF–REF/ALT',
           bands_ref_ref='90; 199', bands_ref_alt='90; 101; 199', bands_alt_alt='101; 199', product_size=300,
           amplicon_start=1, amplicon_end=len(SEQUENCE), snp_offset_in_amplicon=12, primers_left='AAAAAAAAAAGAAT',
           primers_right='GGGGGGGGGGGGGG', tm_left=58.2, tm_right=61.0, pattern='gain/loss', assay_rank=1,
           design='natural')


class DetailsPanelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self):
        records = [InputRecord(2, 'HOXB13', '1', 13, 'A', 'G', '1\t13\tA\tG\tHOXB13'),
                   InputRecord(3, 'OLFML3', '1', 91, 'C', 'T', ''),
                   InputRecord(4, 'bad', '1', None, '', '', '1\tx\tA\tG\tbad')]
        rows = [core_row(**HOX),
                core_row(**dict(HOX, enzyme='MluCI', isoschizomers='', design='dCAPS', genotype_margin=2.3,
                                primer_mismatches="right 3'-2 G>C", assay_rank=2)),
                core_row(variant='1:91 C>T', mapped_id='1', enzyme='', status='no_discriminating_assay',
                         reason='Не найдено сочетание для различения трёх генотипов')]
        rejected = [{'line': 4, 'reason': 'parse_error', 'details': "invalid literal for int() with base 10: 'x'"}]
        self.hoxb13, self.olfml3, self.bad = snp_summary.build_summaries(records, HEADER, rows, rejected)
        self.panel = DetailsPanel()

    def tearDown(self):
        self.panel.close_fasta()
        self.panel.close()

    def test_empty_state(self):
        self.panel.show_summary(None, 'fallback')
        self.assertEqual(self.panel.findChild(QLabel, 'message_label').text(), 'Выберите SNP в списке')
        self.assertTrue(self.panel.findChild(QTabWidget, 'tabs').isHidden())

    def test_header_and_scheme_rows(self):
        self.panel.show_summary(self.hoxb13, 'fallback')
        self.assertEqual(self.panel.findChild(QLabel, 'title_label').text(), 'HOXB13 · 1:13 A>G')
        self.assertEqual(self.panel.findChild(QLabel, 'subtitle_label').text(), '1 · строка 2 файла')
        self.assertEqual(self.panel.findChild(QLabel, 'status_badge').text(), 'Надёжно')
        rows = self.panel.scheme_texts()
        self.assertIn('EcoRI', rows[0])
        self.assertIn('+2 изошизомера', rows[0])
        self.assertIn('лучшая', rows[0])
        self.assertIn('dCAPS', rows[1])
        self.assertEqual(self.panel.findChild(QTabWidget, 'tabs').currentIndex(), 0)
        self.assertEqual(self.panel.findChild(QListWidget, 'schemes_list').count(), 2)

    def test_star_changes_choice(self):
        seen = []
        self.panel.choice_changed.connect(seen.append)
        self.panel.show_summary(self.hoxb13, 'fallback')
        self.panel.choose_scheme(1)
        self.assertEqual((self.hoxb13.chosen, seen), (1, [self.hoxb13]))
        self.assertEqual(self.panel.viewed_index(), 1)
        self.assertEqual(self.panel.findChild(QLabel, 'status_badge').text(), 'Надёжно')

    def test_viewing_does_not_choose(self):
        self.panel.show_summary(self.hoxb13, 'fallback')
        self.panel.view_scheme(1)
        self.assertEqual((self.panel.viewed_index(), self.hoxb13.chosen), (1, 0))

    def test_copy_for_order(self):
        self.panel.show_summary(self.hoxb13, 'fallback')
        self.panel.findChild(QPushButton, 'copy_order_button').click()
        self.assertEqual(QApplication.clipboard().text(), order_lines(self.hoxb13))
        self.panel.findChild(QPushButton, 'copy_all_button').click()
        self.assertTrue(QApplication.clipboard().text().startswith('1:13 A>G (HOXB13)'))

    def test_summary_without_scheme_shows_reason_and_hint(self):
        self.panel.show_summary(self.olfml3, 'fallback')
        self.assertEqual(self.panel.findChild(QLabel, 'message_label').text(), self.olfml3.note('fallback'))
        self.assertTrue(self.panel.findChild(QTabWidget, 'tabs').isHidden())
        self.assertEqual(self.panel.findChild(QLabel, 'status_badge').text(), 'Нет схемы')

    def test_parse_error_shows_the_raw_line(self):
        self.panel.show_summary(self.bad, 'fallback')
        self.assertIn('Строка не прочитана', self.panel.findChild(QLabel, 'message_label').text())
        raw = self.panel.findChild(QLabel, 'raw_label')
        self.assertEqual((raw.isHidden(), raw.text()), (False, '1\tx\tA\tG\tbad'))

    def test_gel_and_weakest_line_for_a_real_genome(self):
        with tempfile.TemporaryDirectory() as directory:
            fasta = write_text(directory, 'g.fa', '>1\n' + SEQUENCE + '\n')
            self.panel.set_context(fasta, dict(min_frag=80, max_frag=800))
            self.panel.show_summary(self.hoxb13, 'fallback')
            self.assertEqual(self.panel.gel_message(), '')
            self.assertEqual(self.panel.findChild(QLabel, 'weakest_label').text(),
                             'Труднее всего различить: REF/REF–REF/ALT; оценка при ухудшении 10.5')
            self.panel.close_fasta()

    def test_moved_genome_shows_error_in_tabs(self):        # Review Focus 5
        self.panel.set_context(r'C:\moved\genome.fa', {})
        self.panel.show_summary(self.hoxb13, 'fallback')
        self.assertIn('Не удалось открыть FASTA', self.panel.gel_message())
        self.assertFalse(self.panel.findChild(QTabWidget, 'tabs').isHidden())

    def test_large_window_follows_the_viewed_scheme(self):
        self.panel.show_summary(self.hoxb13, 'fallback')
        self.panel.open_large()
        dialog = self.panel.large_dialog()
        self.addCleanup(dialog.close)
        self.assertTrue(dialog.isVisible())
        self.panel.show_summary(self.olfml3, 'fallback')
        self.assertEqual(dialog.gel.message(), self.olfml3.note('fallback'))


if __name__ == '__main__':
    unittest.main()
