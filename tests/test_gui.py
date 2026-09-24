import os
import unittest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QApplication

import qtprimer3_vis


class GuiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def test_numeric_sort_filter_and_empty_reset(self):
        window = qtprimer3_vis.RFLPPickerGUI()
        header = ['variant', 'mapped_id', 'enzyme', 'site', 'pattern', 'diag_delta_bp', 'status', 'reason', 'delta_threshold_bp']
        rows = [
            ['1:1 A>G', '1', 'Zeta', 'AC', 'shift', 120, 'window_only', '', 25],
            ['1:2 C>T', '1', 'Alpha', 'GT', 'gain/loss', 50, 'no_primers', 'No pair', 25],
        ]
        window.show_results(header, rows)
        window.result_table.sortItems(5, Qt.SortOrder.AscendingOrder)
        self.assertEqual([window.result_table.item(i, 5).text() for i in range(2)], ['50', '120'])
        window.filter_enzyme_combo.setCurrentText('Alpha')
        self.assertEqual([window.result_table.item(i, 2).text() for i in range(2) if not window.result_table.isRowHidden(i)], ['Alpha'])
        window.show_results([], [])
        self.assertEqual((window.result_table.rowCount(), window.result_table.columnCount()), (0, 0))
        window.close()

    def test_joint_controls_summary_and_full_precision_export(self):
        window = qtprimer3_vis.RFLPPickerGUI()
        self.addCleanup(window.close)
        window.use_primer3_chk.setChecked(True)
        self.assertTrue(window.primer_pairs_spin.isEnabled())
        window.joint_design_chk.setChecked(False)
        self.assertFalse(window.primer_pairs_spin.isEnabled())
        window.joint_design_chk.setChecked(True)
        window.use_primer3_chk.setChecked(False)
        self.assertFalse(window.joint_design_chk.isEnabled())
        header = ['variant', 'mapped_id', 'enzyme', 'site', 'pattern',
                  'status', 'reason', 'genotype_quality', 'genotype_margin',
                  'worst_margin', 'bands_ref_ref', 'bands_ref_alt', 'bands_alt_alt']
        row = ['1:2 A>G', '1', 'EcoRI', 'GAATTC', 'gain/loss', 'ok', 'Model only',
               'robust', 8.123456789, 5.123456789, '100; 200', '100; 200; 300', '300']
        window.show_results(header, [row])
        window.result_table.selectRow(0)
        self.assertIn('REF/ALT: 100; 200; 300', window.assay_summary.text())
        self.assertEqual(window._collect_table_data(), (header, [row]))
        window.show_results([], [])
        self.assertNotIn('100; 200', window.assay_summary.text())


if __name__ == '__main__':
    unittest.main()
