"""The specificity check inside the analysis (specificity.py is tested on its own)."""
import random
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import rflp_core as core
import specificity
from tests.helpers import CHECKS_OFF

COMPLEMENT = str.maketrans('ACGT', 'TGCA')


class SpecificityRunTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        rng = random.Random(19)                       # the natural EcoRI scenario of JointDesignTests
        sequence = ''.join(rng.choice('ACGT') for _ in range(400))
        self.sequence = sequence[:197] + 'GAATTC' + sequence[203:]
        self.spacer = ''.join(rng.choice('ACGT') for _ in range(900))
        (self.root / 'variants.tsv').write_text('chrom\tpos\tref\talt\n1\t200\tA\tG\n')

    def run_core(self, contigs: dict, report=None, **overrides):
        fasta = self.root / f'g{len(list(self.root.glob("*.fa")))}.fa'
        fasta.write_text(''.join(f'>{name}\n{sequence}\n' for name, sequence in contigs.items()))
        params = dict(core.DEFAULT_PARAMS, **CHECKS_OFF, fasta=str(fasta), input_file=str(self.root / 'variants.tsv'),
                      flank=500, prod_min=200, prod_max=500, tm_min=50, tm_max=70, primer_pairs=3, top_results=0)
        params.update(overrides)
        header, rows = core.run_analysis(params, report=report)
        return [dict(zip(header, row)) for row in rows]

    def pairs(self, records):
        return sorted({(r['primers_left'], r['primers_right']) for r in records if r['enzyme']})

    def test_columns_are_empty_when_the_check_is_off(self):
        with patch.object(specificity, 'find_sites', side_effect=AssertionError('no scan expected')):
            records = self.run_core({'1': self.sequence}, specificity='off')
        self.assertTrue(records)
        for name in core.V16_HEADER:
            self.assertEqual({r[name] for r in records}, {None}, name)

    def test_unique_primers_pass(self):
        report = {}
        records = self.run_core({'1': self.sequence, '2': self.spacer}, report=report, specificity='exclude')
        self.assertEqual(self.pairs(records), self.pairs(self.run_core({'1': self.sequence}, specificity='off')))
        self.assertEqual({(r['primer_sites_left'], r['primer_sites_right'], r['offtarget_products'],
                           r['offtarget_sizes']) for r in records}, {(1, 1, 0, '')})
        self.assertEqual(report['summary']['excluded_offtarget'], 0)
        self.assertGreater(report['summary']['specificity_primers'], 1)

    def planted(self, left, right, size=700):
        """A second contig where the pair gives a product of ``size`` bp."""
        reverse = right.translate(COMPLEMENT)[::-1]
        return self.spacer[:100] + left + self.spacer[100 + len(left):100 + size - len(reverse)] + reverse + self.spacer[:60]

    def test_pair_with_a_product_elsewhere_is_flagged_or_excluded(self):
        every = self.pairs(self.run_core({'1': self.sequence}, specificity='off'))
        self.assertGreater(len(every), 1)
        bad = every[0]
        genome = {'1': self.sequence, '2': self.planted(*bad)}
        warned = self.run_core(genome, specificity='warn')
        self.assertEqual(self.pairs(warned), every)
        for record in warned:
            if (record['primers_left'], record['primers_right']) == bad:
                self.assertEqual((record['offtarget_products'], record['offtarget_sizes']), (1, '700'))
                self.assertEqual((record['primer_sites_left'], record['primer_sites_right']), (2, 2))
        report = {}
        kept = self.run_core(genome, report=report, specificity='exclude')
        self.assertNotIn(bad, self.pairs(kept))
        self.assertTrue(self.pairs(kept))
        self.assertEqual({r['offtarget_products'] for r in kept}, {0})
        self.assertGreater(report['summary']['excluded_offtarget'], 0)
        self.assertEqual([r['assay_rank'] for r in kept], list(range(1, len(kept) + 1)))

    def test_product_longer_than_the_limit_does_not_count(self):
        every = self.pairs(self.run_core({'1': self.sequence}, specificity='off'))
        genome = {'1': self.sequence, '2': self.planted(*every[0], size=700)}
        records = self.run_core(genome, specificity='exclude', specificity_max_product=600)
        self.assertEqual(self.pairs(records), every)
        self.assertEqual({r['offtarget_products'] for r in records}, {0})

    def test_duplicated_locus_leaves_no_specific_primers(self):
        report = {}
        (record,) = self.run_core({'1': self.sequence, '2': self.sequence}, report=report, specificity='exclude')
        self.assertEqual(record['status'], 'no_specific_primers')
        self.assertIn('посторонние продукты', record['reason'])
        self.assertEqual(record['enzyme'], '')
        self.assertEqual(report['summary']['variants_with_candidates'], 0)
        warned = self.run_core({'1': self.sequence, '2': self.sequence}, specificity='warn')
        self.assertEqual({r['offtarget_products'] for r in warned}, {1})

    def test_window_search_is_not_checked(self):
        with patch.object(specificity, 'find_sites', side_effect=AssertionError('no scan expected')):
            records = self.run_core({'1': self.sequence}, specificity='exclude', use_primer3=False, flank=150, min_frag=20)
        self.assertTrue(records)

    def test_dcaps_is_tried_when_the_robust_natural_assays_are_not_specific(self):
        report, scans = {}, []
        real = specificity.find_sites

        def counting(*args, **kwargs):
            scans.append(len(list(args[1])))
            return real(*args, **kwargs)

        with patch.object(specificity, 'find_sites', counting):
            self.run_core({'1': self.sequence, '2': self.sequence}, report=report, specificity='exclude',
                          dcaps='fallback', dcaps_tm_max=70.0)
        self.assertEqual(report['summary']['dcaps_snps'], 1)       # not tried at first: a robust natural assay existed
        self.assertGreaterEqual(len(scans), 1)

    def test_scan_can_be_cancelled(self):
        calls = []

        def cancel():
            calls.append(1)
            return specificity_started[0]

        specificity_started = [False]
        real = specificity.find_sites

        def started(*args, **kwargs):
            specificity_started[0] = True
            return real(*args, **kwargs)

        fasta = self.root / 'c.fa'
        fasta.write_text(f'>1\n{self.sequence}\n')
        params = dict(core.DEFAULT_PARAMS, **CHECKS_OFF, fasta=str(fasta), input_file=str(self.root / 'variants.tsv'),
                      flank=500, prod_min=200, prod_max=500, tm_min=50, tm_max=70)
        params['specificity'] = 'exclude'
        report = {}
        with patch.object(specificity, 'find_sites', started), self.assertRaises(core.AnalysisCancelled):
            core.run_analysis(params, cancel_cb=cancel, report=report)
        self.assertTrue(report['cancelled'])

    def test_parameters_are_validated(self):
        base = dict(core.DEFAULT_PARAMS, fasta=__file__, input_file=__file__)
        self.assertEqual(core.validate_params(base), [])
        self.assertEqual((base['specificity'], base['specificity_mismatches'], base['specificity_max_product']),
                         ('exclude', 3, 2000))
        self.assertEqual(core.validate_params(dict(base, specificity='yes')),
                         ['specificity: допустимо off, warn или exclude'])
        self.assertEqual(core.validate_params(dict(base, specificity_mismatches=9)),
                         ['specificity_mismatches: требуется целое число от 0 до 5'])
        self.assertEqual(core.validate_params(dict(base, specificity_max_product=10)),
                         ['specificity_max_product: требуется целое число от 100 до 20000'])


if __name__ == '__main__':
    unittest.main()
