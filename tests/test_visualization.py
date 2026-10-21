import os
import tempfile
import unittest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
from PyQt6.QtWidgets import QApplication

from visualization import PrimerGraphicPanel, RFLPVisualizationDialog, gel_lanes


class VisualizationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def test_render_and_export(self):
        panel = PrimerGraphicPanel()
        panel.render(mapped_id='1', genome_start=1, genome_end=40, amplicon_seq='ACGT' * 10,
                     snp_offset=10, primer_left='ACGTAC', primer_right='GTACGT',
                     tm_left=60.0, tm_right=61.0, primer_left_start=0, primer_left_len=6,
                     primer_right_start=34, primer_right_len=6, variant_str='1:11 G>T',
                     enzyme_name='', pattern='shift', frags_ref='20+20', frags_alt='10+30')
        with tempfile.TemporaryDirectory() as directory:
            png = os.path.join(directory, 'diagram.png')
            svg = os.path.join(directory, 'diagram.svg')
            panel.export_png(png)
            panel.export_svg(svg)
            self.assertGreater(os.path.getsize(png), 1000)
            self.assertGreater(os.path.getsize(svg), 1000)
            self.assertEqual(len(panel._seq_ref), 40)
            self.assertEqual(len(panel._seq_alt), 40)

    def test_gel_lanes_amounts(self):
        lanes = dict(gel_lanes([100, 200, 100], [400], 400))
        self.assertEqual(lanes['Без фермента'], [(400, 1.0)])
        self.assertEqual(lanes['REF/REF'], [(200, 0.5), (100, 0.5)])
        self.assertEqual(lanes['REF/ALT'], [(400, 0.5), (200, 0.25), (100, 0.25)])
        self.assertEqual(lanes['ALT/ALT'], [(400, 1.0)])

    def test_dialog_draws_the_gel_and_exports(self):
        dialog = RFLPVisualizationDialog()
        self.addCleanup(dialog.close)
        sequence = 'AAAAAAAAAA' + 'GAATTC' + 'A' * 84 + 'GAATTC' + 'C' * 194
        dialog.show_payload(dict(mapped_id='1', genome_start=1, genome_end=len(sequence), amplicon_seq=sequence,
                                 snp_offset=12, primer_left='', primer_right='', tm_left=None, tm_right=None,
                                 variant_str='1:13 A>G', enzyme_name='EcoRI', pattern='gain/loss',
                                 gel=dict(min_visible=80, max_visible=800)))
        lanes = dict(dialog.gel._lanes)
        self.assertEqual([length for length, _ in lanes['REF/REF']], [199, 90, 11])
        self.assertEqual([length for length, _ in lanes['ALT/ALT']], [199, 101])
        with tempfile.TemporaryDirectory() as directory:
            for name in ('gel.png', 'gel.svg'):
                path = os.path.join(directory, name)
                (dialog.gel.export_png if name.endswith('png') else dialog.gel.export_svg)(path)
                self.assertGreater(os.path.getsize(path), 1000)
        dialog.show_payload(None, 'нет данных')
        self.assertEqual(dialog.gel._lanes, [])


if __name__ == '__main__':
    unittest.main()
