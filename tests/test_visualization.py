import os
import tempfile
import unittest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
from PyQt6.QtWidgets import QApplication

import i18n
import palette
from tests.helpers import cyrillic_texts, dispose, scene_texts, write_text
from visualization import (GelPanel, PrimerGraphicPanel, RFLPVisualizationDialog, gel_lanes, overlay_primers,
                           payload_for_scheme, render_payload)

SAMPLE_PAYLOAD = dict(mapped_id='1', genome_start=1, genome_end=40, amplicon_seq='ACGT' * 10,
                      snp_offset=10, primer_left='ACGTAC', primer_right='GTACGT',
                      tm_left=60.0, tm_right=61.0, primer_left_start=0, primer_left_len=6,
                      primer_right_start=34, primer_right_len=6, variant_str='1:11 G>T',
                      enzyme_name='', pattern='shift', frags_ref='20+20', frags_alt='10+30')


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
        self.addCleanup(dispose, dialog)
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

    def test_compact_map_has_no_sequence_view(self):
        panel = PrimerGraphicPanel(compact=True)
        self.addCleanup(dispose, panel)
        self.assertTrue(panel.seq_view.isHidden())
        self.assertTrue(panel.export_png_btn.isHidden())
        full = PrimerGraphicPanel()
        self.addCleanup(dispose, full)
        self.assertFalse(full.seq_view.isHidden())
        gel = GelPanel(compact=True)
        self.addCleanup(dispose, gel)
        self.assertTrue(gel.export_png_btn.isHidden())

    def test_primers_are_bold_without_background(self):
        panel = PrimerGraphicPanel()
        self.addCleanup(dispose, panel)
        colors = palette.colors(False)
        panel.apply_colors(colors)
        panel.render(**SAMPLE_PAYLOAD)
        html = panel.seq_view.toHtml().replace(' ', '')
        self.assertIn('font-weight:700', html)
        self.assertIn('text-decoration:underline', html)
        self.assertNotIn('#d8f3e8', html.lower())
        self.assertIn(colors['snp_bg'].lower(), html.lower())
        self.assertEqual(panel.view.backgroundBrush().color().name().upper(), colors['panel'])
        panel.apply_colors(palette.colors(True))
        self.assertEqual(panel.view.backgroundBrush().color().name().upper(), palette.DARK['panel'])

    def test_gel_is_dark_in_both_themes(self):
        gel = GelPanel()
        self.addCleanup(dispose, gel)
        gel.render([100, 300], [400], 400, '1:1 A>G', 'AluI')
        self.assertEqual(gel.view.backgroundBrush().color().name().upper(), palette.GEL['background'])

    def test_payload_for_missing_fasta(self):
        payload, error = payload_for_scheme(dict(mapped_id='1', amplicon_start=1, amplicon_end=40),
                                            r'C:\missing\genome.fa', {})
        self.assertIsNone(payload)
        self.assertTrue(error.startswith('Не удалось открыть FASTA'))

    def test_payload_for_scheme_reads_the_amplicon(self):
        with tempfile.TemporaryDirectory() as directory:
            fasta = write_text(directory, 'g.fa', '>1\n' + 'ACGT' * 10 + '\n')
            scheme = dict(mapped_id='1', amplicon_start=5, amplicon_end=24, snp_offset_in_amplicon=3,
                          primers_left='ACGTA', primers_right='ACGTA', tm_left=60.0, tm_right=61.0,
                          variant='1:8 T>C', enzyme='AluI', pattern='shift', frags_ref='20', frags_alt='20')
            payload, error = payload_for_scheme(scheme, fasta, dict(min_frag=50, max_frag=700))
        self.assertEqual(error, '')
        self.assertEqual((payload['genome_start'], payload['genome_end'], len(payload['amplicon_seq'])),
                         (5, 24, 20))
        self.assertEqual(payload['gel'], dict(min_visible=50, max_visible=700))
        self.assertEqual(payload['variant_str'], '1:8 T>C')

    def test_render_payload_without_payload_shows_message(self):
        panel, gel = PrimerGraphicPanel(compact=True), GelPanel(compact=True)
        self.addCleanup(dispose, panel)
        self.addCleanup(dispose, gel)
        render_payload(panel, gel, None, 'Не удалось открыть FASTA: x')
        self.assertEqual(gel.message(), 'Не удалось открыть FASTA: x')

    def test_english_gel_map_and_messages(self):
        i18n.set_language('en')
        self.addCleanup(i18n.set_language, 'ru')
        gel = GelPanel()
        self.addCleanup(dispose, gel)
        gel.render([100, 300], [400], 400, '1:1 A>G', 'AluI')
        self.assertIn('No enzyme', scene_texts(gel.scene))
        self.assertEqual(cyrillic_texts(gel), [])
        dialog = RFLPVisualizationDialog()
        self.addCleanup(dispose, dialog)
        self.assertEqual(cyrillic_texts(dialog), [])
        sequence = 'AAAAAAAAAA' + 'GAATTC' + 'A' * 84 + 'GAATTC' + 'C' * 194
        dialog.show_payload(dict(SAMPLE_PAYLOAD, amplicon_seq=sequence, genome_end=len(sequence), snp_offset=12,
                                 variant_str='1:13 A>G', enzyme_name='EcoRI', pattern='gain/loss',
                                 gel=dict(min_visible=80, max_visible=800)))
        self.assertEqual(cyrillic_texts(dialog), [])
        dialog.show_payload(dict(SAMPLE_PAYLOAD, variant_str='1:11 C>T'))   # REF does not match: no ALT lane
        self.assertEqual(cyrillic_texts(dialog), [])
        payload, error = payload_for_scheme(dict(mapped_id='1', amplicon_start=1, amplicon_end=40),
                                            r'C:\missing\genome.fa', {})
        self.assertTrue(error.startswith('Could not open the FASTA'), error)
        dialog.show_payload(None, error)
        self.assertEqual(cyrillic_texts(dialog), [])

    def test_overlay_primers_moved(self):
        self.assertEqual(overlay_primers('AAAACCCCGGGG', 'AAtA', 'cCCC'), 'AATACCCCGGGG')
        self.assertEqual(overlay_primers('acgtACGT', '', ''), 'ACGTACGT')


if __name__ == '__main__':
    unittest.main()
