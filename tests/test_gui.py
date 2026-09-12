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


if __name__ == '__main__':
    unittest.main()
