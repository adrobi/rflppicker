import os
import tempfile
import unittest
from pathlib import Path

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
from PyQt6.QtWidgets import QApplication

import assay_scoring
from assay_scoring import GelModel, score_genotypes, score_locus
from tests.helpers import dispose, write_text
from visualization import GelPanel, gel_lanes, payload_for_scheme


class GelResolutionTests(unittest.TestCase):
    def test_reference_concentration_preserves_legacy_scores(self):
        legacy = GelModel()
        current = GelModel(version=assay_scoring.AGAROSE_MODEL_VERSION, agarose_pct=2.0)
        short, long = [150, 250], [162, 238]
        self.assertEqual(score_genotypes(short, long, current), score_genotypes(short, long, legacy))
        self.assertEqual(score_locus([short, long, [100, 300]], current),
                         score_locus([short, long, [100, 300]], legacy))

    def test_concentration_changes_which_fragment_sizes_are_distinguishable(self):
        version = assay_scoring.AGAROSE_MODEL_VERSION
        high = GelModel(resolution_pct=0, version=version, agarose_pct=4.0)
        low = GelModel(resolution_pct=0, version=version, agarose_pct=1.0)
        self.assertEqual(score_genotypes([150, 250], [162, 238], high)['genotype_quality'], 'robust')
        self.assertEqual(score_genotypes([150, 250], [162, 238], low)['genotype_quality'], 'ambiguous')
        self.assertEqual(score_genotypes([500, 700], [512, 688], low)['genotype_quality'], 'robust')
        self.assertEqual(score_genotypes([500, 700], [512, 688], high)['genotype_quality'], 'ambiguous')

    def test_legacy_is_concentration_independent_and_brightness_never_rescues_a_call(self):
        legacy = GelModel(agarose_pct=4)
        self.assertEqual(score_genotypes([100, 200], [300], legacy),
                         score_genotypes([100, 200], [300], GelModel()))
        current = GelModel(version=assay_scoring.AGAROSE_MODEL_VERSION, agarose_pct=4)
        self.assertEqual(score_genotypes([100, 200], [100, 100, 100], current)['genotype_quality'], 'ambiguous')

    def test_model_rejects_unknown_version_and_nonfinite_concentration(self):
        with self.assertRaises(ValueError):
            GelModel(version='unknown')
        with self.assertRaises(ValueError):
            GelModel(agarose_pct=float('nan'))


class DigestMassTests(unittest.TestCase):
    def test_partial_digest_only_joins_adjacent_fragments_and_adds_equal_length_mass(self):
        bands = dict(dict(gel_lanes([100, 200, 300], [600], 600, digest_pct=50))['REF/REF'])
        self.assertEqual(set(bands), {100, 200, 300, 500, 600})
        for size, share in {100: 1/12, 200: 1/12, 300: 3/8, 500: 5/24, 600: 1/4}.items():
            self.assertAlmostEqual(bands[size], share)
        self.assertAlmostEqual(sum(bands.values()), 1)

    def test_heterozygote_uses_one_copy_of_each_product_even_at_different_lengths(self):
        bands = dict(dict(gel_lanes([100, 200], [400], 300, 400, digest_pct=50))['REF/ALT'])
        self.assertEqual(set(bands), {100, 200, 300, 400})
        for size, share in {100: 1/14, 200: 1/7, 300: 3/14, 400: 4/7}.items():
            self.assertAlmostEqual(bands[size], share)
        self.assertAlmostEqual(sum(bands.values()), 1)

    def test_zero_completion_leaves_products_intact_and_uncut_view_ignores_completion(self):
        lanes = dict(gel_lanes([100, 200], [400], 300, 400, digest_pct=0))
        self.assertEqual(lanes['REF/REF'], [(300, 1.0)])
        self.assertEqual(lanes['REF/ALT'], [(400, 4/7), (300, 3/7)])
        self.assertEqual(lanes['ALT/ALT'], [(400, 1.0)])
        self.assertEqual(gel_lanes([300], [400], 300, 400, uncut=True, digest_pct=0),
                         gel_lanes([300], [400], 300, 400, uncut=True))


class GelDrawingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def test_render_and_export_use_selected_ladder_and_partial_digest(self):
        gel = GelPanel()
        self.addCleanup(dispose, gel)
        gel.render([100, 200], [300], 300, marker='1-kb', digest_pct=50)
        tips = {item.toolTip() for item in gel.scene.items() if item.toolTip()}
        self.assertIn('Маркер: 10000 п.н.', tips)
        self.assertNotIn('Маркер: 50 п.н.', tips)
        self.assertIn('REF/REF: 300 п.н., доля ДНК 50%', tips)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'gel.svg'
            gel.export_svg(path)
            svg = path.read_text(encoding='utf-8')
        self.assertIn('1 kb', svg)
        self.assertIn('Вероятность разреза на сайт', svg)

    def test_nominal_visible_clusters_follow_the_selected_concentration(self):
        gel = GelPanel()
        self.addCleanup(dispose, gel)
        version = assay_scoring.AGAROSE_MODEL_VERSION
        gel.render([150, 250], [162, 238], 400, model_version=version, agarose_pct=1, resolution_pct=0)
        tips = {item.toolTip() for item in gel.scene.items() if item.toolTip()}
        self.assertIn('REF/ALT: 150–162 п.н., доля ДНК 39%', tips)
        gel.render([150, 250], [162, 238], 400, model_version=version, agarose_pct=4, resolution_pct=0)
        tips = {item.toolTip() for item in gel.scene.items() if item.toolTip()}
        self.assertIn('REF/ALT: 150 п.н., доля ДНК 19%', tips)
        self.assertIn('REF/ALT: 162 п.н., доля ДНК 20%', tips)

    def test_multiallelic_lane_includes_partial_bands_from_each_selected_allele(self):
        gel = GelPanel()
        self.addCleanup(dispose, gel)
        gel.render([100, 200], [300], 300, digest_pct=50, locus_products=[
            {'length': 300, 'fragments': [100, 200]},
            {'length': 300, 'fragments': [300]},
            {'length': 340, 'fragments': [140, 200]},
        ])
        lanes = dict(gel._lanes)
        self.assertEqual(len(lanes), 6)
        bands = dict(lanes['ALT1/ALT2'])
        self.assertEqual(set(bands), {140, 200, 300, 340})
        self.assertAlmostEqual(bands[300], 300/640)
        self.assertAlmostEqual(bands[340], 170/640)

    def test_allele_specific_gel_uses_six_reaction_lanes_with_shared_control(self):
        gel = GelPanel()
        self.addCleanup(dispose, gel)
        gel.render_aspcr({
            'REF/REF': {'REF': [140, 300], 'ALT': [300]},
            'REF/ALT': {'REF': [140, 300], 'ALT': [160, 300]},
            'ALT/ALT': {'REF': [300], 'ALT': [160, 300]},
        })
        lanes = dict(gel._lanes)
        self.assertEqual(len(lanes), 6)
        self.assertTrue(all((300, None) in bands for bands in lanes.values()))
        self.assertEqual(lanes['REF/REF · ALT'], [(300, None)])
        self.assertEqual(lanes['REF/ALT · ALT'], [(300, None), (160, None)])
        tips = {item.toolTip() for item in gel.scene.items() if item.toolTip()}
        self.assertIn('ALT/ALT · ALT: 160 п.н.', tips)
        self.assertFalse(any('доля ДНК' in tip for tip in tips))

    def test_payload_preserves_analysis_and_view_settings(self):
        params = dict(gel_model_version=assay_scoring.AGAROSE_MODEL_VERSION, gel_agarose_pct=3,
                      gel_resolution_bp=13, gel_resolution_pct=5, gel_stress_factor=2,
                      gel_marker='100-bp', gel_digest_pct=70)
        with tempfile.TemporaryDirectory() as directory:
            fasta = write_text(directory, 'seq.fa', '>1\nACGT\n')
            payload, error = payload_for_scheme(dict(mapped_id='1', amplicon_start=1, amplicon_end=4), fasta, params)
        self.assertEqual(error, '')
        self.assertEqual(payload['gel'], dict(min_visible=80, max_visible=800, resolution_bp=13,
                                             resolution_pct=5, stress_factor=2, agarose_pct=3,
                                             model_version=assay_scoring.AGAROSE_MODEL_VERSION,
                                             marker='100-bp', digest_pct=70))


    def test_payload_recovers_row_settings_and_respects_supplied_view_parameters(self):
        scheme = dict(mapped_id='1', amplicon_start=1, amplicon_end=4,
                      scoring_model=assay_scoring.AGAROSE_MODEL_VERSION,
                      gel_agarose_pct='4.0', gel_marker='1-kb', gel_digest_pct='75',
                      gel_min_visible='90', gel_max_visible='750', gel_resolution_bp='12.5')
        with tempfile.TemporaryDirectory() as directory:
            fasta = write_text(directory, 'seq.fa', '>1\nACGT\n')
            payload, error = payload_for_scheme(scheme, fasta, {'gel_marker': '100-bp'})
        self.assertEqual(error, '')
        self.assertEqual(payload['gel'], dict(model_version=assay_scoring.AGAROSE_MODEL_VERSION,
                                             agarose_pct=4, marker='100-bp', digest_pct=75,
                                             min_visible=90, max_visible=750, resolution_bp=12.5))


if __name__ == '__main__':
    unittest.main()
