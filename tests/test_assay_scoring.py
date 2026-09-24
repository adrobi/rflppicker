import itertools
import random
import unittest

from assay_scoring import GelModel, score_genotypes, rank_key


class ScoringTests(unittest.TestCase):
    def setUp(self):
        self.model = GelModel(50, 1000, 10, 0, 1.5)

    def test_heterozygote_subset_is_ambiguous_despite_distinct_homozygotes(self):
        result = score_genotypes([100, 200], [100, 100, 100], self.model)
        self.assertEqual(result['genotype_quality'], 'ambiguous')
        self.assertEqual(result['bands_ref_ref'], result['bands_ref_alt'])
        self.assertEqual(result['genotype_margin'], 0)

    def test_known_gain_loss_and_worst_case(self):
        result = score_genotypes([100, 200], [300], self.model)
        self.assertEqual(result['genotype_quality'], 'robust')
        self.assertEqual(result['genotype_margin'], 10)
        self.assertAlmostEqual(result['worst_margin'], 100 / 15)
        self.assertEqual(result['bands_ref_alt'], '100; 200; 300')

    def test_equal_total_length_shift_loses_resolution_under_stress(self):
        result = score_genotypes([100, 200], [114, 186], self.model)
        self.assertEqual(result['genotype_quality'], 'nominal_only')
        self.assertAlmostEqual(result['genotype_margin'], 1.4)
        self.assertEqual(result['worst_margin'], 0)

    def test_boundary_is_unresolved_and_relative_resolution_matters(self):
        self.assertEqual(score_genotypes([100, 200], [110, 190], self.model)['genotype_margin'], 0)
        result = score_genotypes([1000, 2000], [1015, 1985], GelModel(1, 3000, 10, 3, 1))
        self.assertEqual(result['genotype_quality'], 'ambiguous')

    def test_invisible_allele_cannot_be_a_positive_call(self):
        result = score_genotypes([200], [40] * 5, self.model)
        self.assertEqual(result['hidden_alt_count'], 5)
        self.assertEqual(result['genotype_margin'], 0)

    def test_duplicates_and_permutation_do_not_change_score(self):
        expected = score_genotypes([100, 200], [300], self.model)
        actual = score_genotypes([200, 100, 100], [300, 300], self.model)
        self.assertEqual(actual, expected)

    def test_shared_bins_merge_transitively(self):
        result = score_genotypes([100, 118], [109, 109], self.model)
        self.assertEqual(result['bands_ref_ref'], '100–118')
        self.assertEqual(result['bands_ref_alt'], result['bands_alt_alt'])
        self.assertEqual(result['genotype_margin'], 0)

    def test_score_matches_exhaustive_presence_classification(self):
        # Independent oracle: graph connectivity, then compare three bit sets.
        rng = random.Random(17)
        for _ in range(100):
            ref = [rng.randrange(1, 50) * 10 for _ in range(4)]
            alt = [rng.randrange(1, 50) * 10 for _ in range(4)]
            model = GelModel(50, 1000, 20, 4, 1.5)
            visible_r = set(x for x in ref if x >= 50)
            visible_a = set(x for x in alt if x >= 50)
            remaining = visible_r | visible_a
            components = []
            while remaining:
                component = {next(iter(remaining))}
                while True:
                    neighbors = {b for b in remaining for a in component
                                 if abs(a - b) <= max(20, .04 * max(a, b))}
                    if neighbors <= component:
                        break
                    component |= neighbors
                components.append(component)
                remaining -= component
            r = frozenset(i for i, c in enumerate(components) if c & visible_r)
            a = frozenset(i for i, c in enumerate(components) if c & visible_a)
            expected = bool(r and a and len({r, a, r | a}) == 3)
            result = score_genotypes(ref, alt, model)
            self.assertEqual(result['genotype_margin'] > 1, expected)
            # Verify upper endpoint really bounds sampled intermediate scores.
            for scale in (1, 1.1, 1.3, 1.5):
                sample = score_genotypes(ref, alt, GelModel(50, 1000, 20, 4, scale))
                self.assertGreaterEqual(sample['worst_margin'], result['worst_margin'])

    def test_swap_symmetry(self):
        for ref, alt in itertools.product(([100, 200], [300], [100, 100, 100]), repeat=2):
            a, b = score_genotypes(ref, alt, self.model), score_genotypes(alt, ref, self.model)
            self.assertEqual(a['genotype_margin'], b['genotype_margin'])
            self.assertEqual(a['worst_margin'], b['worst_margin'])

    def test_rank_uses_worst_margin_before_primer_penalty(self):
        common = dict(genotype_margin=5, visible_band_count=3, enzyme='EcoRI', primer_pair_index=1)
        best = dict(common, worst_margin=3, primer_pair_penalty=8)
        poor = dict(common, worst_margin=0, primer_pair_penalty=1)
        self.assertEqual(sorted([poor, best], key=rank_key)[0], best)

    def test_invalid_parameters(self):
        for kw in ({'resolution_bp': 0}, {'resolution_pct': float('nan')},
                   {'stress_factor': .5}, {'min_visible': 900}):
            with self.assertRaises(ValueError):
                GelModel(**kw)


if __name__ == '__main__':
    unittest.main()
