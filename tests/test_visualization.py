import os
import tempfile
import unittest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
from PyQt6.QtWidgets import QApplication

from visualization import PrimerGraphicPanel


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


if __name__ == '__main__':
    unittest.main()
