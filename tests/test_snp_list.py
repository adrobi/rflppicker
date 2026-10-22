import os
import tempfile
import unittest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
from PyQt6.QtCore import QSettings, Qt
from PyQt6.QtWidgets import QApplication

import snp_summary
from app_settings import AppSettings
from snp_list import DEFAULT_COLUMNS, SORT_ROLE, SnpListView
from snp_summary import InputRecord
from tests.helpers import HEADER, core_row, dispose

HOX = dict(variant='11:37525128 C>G', mapped_id='NC_056064.1', enzyme='Bse1I', isoschizomers='BsrI',
           genotype_quality='robust', genotype_margin=15.77, bands_ref_ref='448', bands_alt_alt='212; 236',
           product_size=475, assay_rank=1)


def fixture_summaries():
    records = [InputRecord(2, 'HOXB13', '11', 37525128, 'C', 'G', ''),
               InputRecord(3, 'LEPR', '1', 41250567, 'A', 'C', ''),
               InputRecord(4, 'OLFML3', '1', 91308651, 'C', 'T', ''),
               InputRecord(5, 'rs7781', '12', 880121, 'A', 'G', '')]
    rows = [core_row(**HOX),
            core_row(**dict(HOX, enzyme='BsrDI', isoschizomers='', genotype_margin=12.1, assay_rank=2)),
            core_row(variant='1:41250567 A>C', enzyme='AciI', genotype_quality='nominal_only', genotype_margin=1.4,
                     bands_ref_ref='292; 140', bands_alt_alt='432', product_size=432, assay_rank=1,
                     weakest_genotypes='REF/REF–REF/ALT'),
            core_row(variant='1:91308651 C>T', enzyme='', status='no_discriminating_assay',
                     reason='Не найдено сочетание для различения трёх генотипов')]
    rejected = [{'chrom': '12', 'pos': 880121, 'ref': 'A', 'alt': 'G', 'reason': 'ref_not_match',
                 'details': 'REF does not match'}]
    return snp_summary.build_summaries(records, HEADER, rows, rejected)


class SnpListTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self):
        self._dir = tempfile.TemporaryDirectory()
        self.qs = QSettings(os.path.join(self._dir.name, 's.ini'), QSettings.Format.IniFormat)
        self.view = SnpListView(AppSettings(self.qs))
        self.summaries = fixture_summaries()
        self.view.set_summaries(self.summaries, 'fallback', has_ids=True)

    def tearDown(self):
        dispose(self.view)
        self.qs.sync()
        self._dir.cleanup()

    def ids_in_view(self):
        column = self.view.column_index('id')
        return [self.view.proxy.index(r, column).data() for r in range(self.view.proxy.rowCount())]

    def cell(self, row, key):
        return self.view.proxy.index(row, self.view.column_index(key)).data()

    def test_cells(self):
        self.assertEqual(self.cell(0, 'snp'), '11:37525128 C>G')
        self.assertEqual(self.cell(0, 'bands'), '448 | 212; 236')
        self.assertEqual(self.cell(0, 'score'), '15,8')
        self.assertEqual(self.cell(0, 'schemes'), '2')
        self.assertEqual(self.cell(2, 'score'), '')
        self.assertEqual(self.cell(3, 'note'), 'REF не совпадает с геномом')
        score = self.view.proxy.index(0, self.view.column_index('score'))
        self.assertEqual(score.data(SORT_ROLE), 15.77)
        self.assertTrue(score.data(Qt.ItemDataRole.TextAlignmentRole) & Qt.AlignmentFlag.AlignRight)
        self.assertEqual(self.view.table.verticalHeader().defaultSectionSize(), 26)

    def test_counts_and_chip_filter(self):
        self.assertEqual(self.view.model.counts(), {'robust': 1, 'nominal': 1, 'none': 1, 'error': 1})
        self.view.chips.toggled.emit({'robust', 'nominal'})
        self.assertEqual([s.record.snp_id for s in self.view.shown_summaries()], ['HOXB13', 'LEPR'])
        self.assertIn('Только показанные (2)', [a.text() for a in self.view.export_menu.actions()])
        self.assertEqual(self.view.footer.text(), 'Показано 2 из 4 · порядок: как в файле')

    def test_chip_buttons_toggle_groups(self):
        self.view.chips.chip('robust').click()
        self.assertEqual([s.record.snp_id for s in self.view.shown_summaries()], ['HOXB13'])
        self.view.chips.chip('error').click()
        self.assertEqual([s.record.snp_id for s in self.view.shown_summaries()], ['HOXB13', 'rs7781'])
        self.view.chips.chip('all').click()
        self.assertEqual(len(self.view.shown_summaries()), 4)
        self.assertEqual(self.view.chips.chip('robust').text(), 'Надёжно 1')

    def test_error_chip_hidden_without_errors(self):
        self.view.set_summaries(self.summaries[:3], 'fallback', has_ids=True)
        self.assertTrue(self.view.chips.chip('error').isHidden())

    def test_search_matches_isoschizomers(self):
        self.view.search.setText('bsri')
        self.assertEqual([s.record.snp_id for s in self.view.shown_summaries()], ['HOXB13'])
        self.view.search.setText('41250567')
        self.assertEqual([s.record.snp_id for s in self.view.shown_summaries()], ['LEPR'])

    def test_empty_filter_shows_reset(self):
        self.view.search.setText('nothing like this')
        self.assertFalse(self.view.placeholder.isHidden())
        self.assertEqual(self.view.placeholder.text(), 'Нет SNP с такими условиями')
        self.view.reset_filter_button.click()
        self.assertEqual(len(self.view.shown_summaries()), 4)
        self.assertTrue(self.view.placeholder.isHidden())

    def test_score_sort_keeps_empty_last_and_third_click_restores_file_order(self):
        header = self.view.table.horizontalHeader()
        column = self.view.column_index('score')
        header.sectionClicked.emit(column)
        self.assertEqual(self.ids_in_view(), ['LEPR', 'HOXB13', 'OLFML3', 'rs7781'])
        self.assertEqual(self.view.footer.text(), 'Показано 4 из 4 · порядок: по столбцу «Оценка»')
        header.sectionClicked.emit(column)
        self.assertEqual(self.ids_in_view(), ['HOXB13', 'LEPR', 'OLFML3', 'rs7781'])
        header.sectionClicked.emit(self.view.column_index('snp'))
        self.assertEqual(self.ids_in_view(), ['LEPR', 'OLFML3', 'HOXB13', 'rs7781'])   # chromosome 1 < 11 < 12
        header.sectionClicked.emit(column)
        header.sectionClicked.emit(column)
        header.sectionClicked.emit(column)
        self.assertEqual(self.ids_in_view(), ['HOXB13', 'LEPR', 'OLFML3', 'rs7781'])
        self.assertIn('как в файле', self.view.footer.text())

    def test_manual_choice_marks_enzyme(self):
        self.summaries[0].choose(1)
        self.view.refresh_current()
        self.assertEqual(self.cell(0, 'enzyme'), 'BsrDI ✎')

    def test_selection_emits_summary(self):
        seen = []
        self.view.selection_changed.connect(seen.append)
        self.view.table.selectRow(1)
        self.assertIs(seen[-1], self.summaries[1])

    def test_id_column_hidden_without_ids(self):
        self.view.set_summaries(self.summaries, 'fallback', has_ids=False)
        self.assertTrue(self.view.table.isColumnHidden(self.view.column_index('id')))

    def test_visible_columns_are_saved(self):
        self.assertTrue(self.view.table.isColumnHidden(self.view.column_index('mapped_id')))
        self.view.set_column_visible('mapped_id', True)
        view = SnpListView(AppSettings(self.qs))
        self.addCleanup(dispose, view)
        view.set_summaries(self.summaries, 'fallback', has_ids=True)
        self.assertFalse(view.table.isColumnHidden(view.column_index('mapped_id')))
        self.assertEqual(DEFAULT_COLUMNS[0], 'dot')

    def test_export_menu_emits_kind(self):
        seen = []
        self.view.export_requested.connect(seen.append)
        for action in self.view.export_menu.actions():
            action.trigger()
        self.assertEqual(seen, ['all', 'shown', 'csv'])

    def test_placeholder(self):
        self.view.set_placeholder('Идёт расчёт… 3 из 10 SNP')
        self.assertFalse(self.view.placeholder.isHidden())
        self.assertTrue(self.view.table.isHidden())


if __name__ == '__main__':
    unittest.main()
