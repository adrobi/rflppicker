import json
import random
import unittest
from collections import Counter
from types import SimpleNamespace
from unittest.mock import patch

from Bio.Seq import Seq

import allele_specific
import rflp_core as core
import specificity
from assay_scoring import GelModel


class TwoTubeSpecificityTests(unittest.TestCase):
    def setUp(self):
        self.outer, self.common, self.ref, self.alt = 'A' * 20, 'C' * 20, 'G' * 20, 'T' * 20
        self.run = SimpleNamespace(params=dict(core.DEFAULT_PARAMS, fasta='unused'), counters=Counter(),
                                   cancel_cb=None, check_cancel=lambda: None, log=lambda _: None)

    def scan(self, site_data):
        control = dict(primers_left=self.outer, primers_right=self.common, start0=0, end0=300)
        reactions = {name: dict(control, primers_left=primer, start0=100)
                     for name, primer in (('REF', self.ref), ('ALT', self.alt))}
        hit = dict(control, design='allele_specific', mapped_id='1', amplicon_start=1, amplicon_end=300,
                   aspcr_control=json.dumps(control), aspcr_reactions=json.dumps(reactions), aspcr_specific_side='left')
        with patch.object(specificity, 'find_sites', return_value=site_data):
            core._scan_hits(self.run, [hit], {})
        return hit

    def test_saturated_inner_primer_is_reported_in_warn_mode(self):
        sites = {primer: specificity.PrimerSites() for primer in (self.outer, self.common, self.ref, self.alt)}
        sites[self.ref].saturated = True
        row = self.scan(sites)
        self.assertTrue(row['_unspecific'])
        self.assertGreater(row['primer_sites_left'], specificity.MAX_SITES)
        self.assertIsNone(row['offtarget_products'])

    def test_shared_offtarget_is_counted_once(self):
        sites = {primer: specificity.PrimerSites() for primer in (self.outer, self.common, self.ref, self.alt)}
        sites[self.common].sites = [specificity.Site('2', 100, '+', 0), specificity.Site('2', 200, '-', 0)]
        row = self.scan(sites)
        self.assertEqual((row['offtarget_products'], row['offtarget_sizes']), (1, '120'))


class AlleleSpecificInterpretationTests(unittest.TestCase):
    def test_published_three_primer_reaction_pattern(self):
        # Hoek et al. (1997), Fig. 2, DOI 10.1016/S0022-1759(97)00185-3.
        # REF corresponds to wt and ALT to gld. Only the band pattern is reproduced.
        observed = {
            'REF/REF': {'REF': [536, 815], 'ALT': [815]},
            'REF/ALT': {'REF': [536, 815], 'ALT': [536, 815]},
            'ALT/ALT': {'REF': [815], 'ALT': [536, 815]},
        }
        self.assertEqual(allele_specific.expected_bands(815, 536, 536), observed)
        for label, lanes in observed.items():
            with self.subTest(label=label):
                self.assertEqual(allele_specific.interpret_genotype(
                    ref_control=815 in lanes['REF'], alt_control=815 in lanes['ALT'],
                    ref_band=536 in lanes['REF'], alt_band=536 in lanes['ALT']), label)

    def test_any_missing_control_or_no_diagnostic_bands_leaves_genotype_undetermined(self):
        for ref_control, alt_control in ((False, True), (True, False), (False, False)):
            with self.subTest(controls=(ref_control, alt_control)):
                self.assertIsNone(allele_specific.interpret_genotype(
                    ref_control=ref_control, alt_control=alt_control, ref_band=True, alt_band=True))
        self.assertIsNone(allele_specific.interpret_genotype(
            ref_control=True, alt_control=True, ref_band=False, alt_band=False))


@unittest.skipUnless(core._PRIMER3_AVAILABLE, 'Primer3 unavailable')
class AlleleSpecificDesignTests(unittest.TestCase):
    def setUp(self):
        self.sequence = ''.join(random.Random(234).choices('ACGT', k=500))
        self.off = 250
        self.ref = self.sequence[self.off]
        self.alt = 'A'
        result = core.primer3_pick(self.sequence, self.off, 250, 500, 50, 70, num_return=4)
        self.controls = list(core._primer_candidates(result, self.sequence, self.off))

    def design(self, sequence=None, **options):
        parameters = dict(control_pairs=self.controls, pick_primers=core.primer3_pick,
                          parse_pairs=core._primer_candidates,
                          prod_min=250, prod_max=500, tm_min=50, tm_max=70, model=GelModel())
        parameters.update(options)
        return allele_specific.design_assays(sequence or self.sequence, self.off,
                                             self.ref, self.alt, **parameters)

    def test_two_reactions_share_the_control_and_terminate_at_snp(self):
        designs = self.design()
        self.assertTrue(designs)
        self.assertEqual({d['aspcr_specific_side'] for d in designs}, {'left', 'right'})
        for design in designs:
            self.assertEqual(design['genotype_quality'], 'robust')
            control = design['aspcr_control']
            reactions = design['aspcr_reactions']
            for label, allele in [('REF', self.ref), ('ALT', self.alt)]:
                reaction = reactions[label]
                side = design['aspcr_specific_side']
                common = 'right' if side == 'left' else 'left'
                self.assertEqual(reaction['allele'], allele)
                self.assertEqual(reaction[f'primers_{common}'], control[f'primers_{common}'])
                primer = reaction[f'primers_{side}']
                self.assertEqual(primer[-1], allele if side == 'left' else str(Seq(allele).complement()))
                if side == 'left':
                    self.assertEqual(reaction['start0'] + len(primer) - 1, self.off)
                else:
                    self.assertEqual(reaction['end0'] - len(primer), self.off)
                self.assertEqual(reaction['product_size'], reaction['end0'] - reaction['start0'])
                self.assertLess(reaction['product_size'], control['product_size'])
                self.assertEqual(design['aspcr_genotypes']['REF/ALT'][label],
                                 [reaction['product_size'], control['product_size']])

    def test_unresolved_invisible_and_excluded_designs_are_rejected(self):
        self.assertEqual(self.design(model=GelModel(resolution_bp=500)), [])
        self.assertEqual(self.design(model=GelModel(max_visible=200)), [])
        self.assertEqual(self.design(excluded_regions=[(self.off - 1, 3)]), [])
        masked = self.sequence[:self.off] + self.ref.lower() + self.sequence[self.off + 1:]
        self.assertEqual(self.design(masked, mask_repeats=True), [])

    def test_reaction_table_lists_three_primers_per_tube(self):
        design = self.design()[0]
        ref_tube, alt_tube = allele_specific.reaction_table(design)
        self.assertEqual([ref_tube['reaction'], alt_tube['reaction']], ['REF', 'ALT'])
        self.assertEqual([ref_tube['allele'], alt_tube['allele']], [self.ref, self.alt])
        control = design['aspcr_control']
        for tube in (ref_tube, alt_tube):
            side = design['aspcr_specific_side']
            diagnostic = design['aspcr_reactions'][tube['reaction']]
            roles = {p['role']: p for p in tube['primers']}
            self.assertEqual(len(tube['primers']), 3)
            self.assertEqual(roles['allele_specific']['sequence'], diagnostic[f'primers_{side}'])
            self.assertEqual(roles['allele_specific']['tm'], diagnostic[f'tm_{side}'])
            self.assertEqual(roles['control_only']['sequence'], control[f'primers_{side}'])
            common_side = roles['common']['side']
            self.assertEqual(roles['common']['sequence'], control[f'primers_{common_side}'])
            self.assertEqual(tube['genotype_bands']['REF/ALT'], [tube['product_size'], tube['control_size']])
        self.assertEqual(allele_specific.reaction_table({}), [])
        persisted = dict(design)
        for field in ('aspcr_control', 'aspcr_reactions', 'aspcr_genotypes'):
            persisted[field] = json.dumps(persisted[field])
        self.assertEqual(allele_specific.reaction_table(persisted), [ref_tube, alt_tube])
