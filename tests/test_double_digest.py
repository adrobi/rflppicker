"""Two natural enzymes on the same linear PCR product."""
from collections import Counter
import itertools
import random
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from Bio.Restriction import ApoI, BsaI, EcoRI, HpaII, MspI, MseI, SmaI, XmaI
from Bio.Seq import Seq

import rflp_core as core
from assay_scoring import GelModel, score_genotypes


def template(partner_cut=50, shared_cut=100):
    # A constitutive EcoRI site and an allele-specific one restored by C>G at position 200.
    rng = random.Random(1)
    while True:
        seq = list(''.join(rng.choices('ACGT', k=300)))
        seq[shared_cut - 1:shared_cut + 5] = 'GAATTC'
        seq[199:205] = 'CAATTC'
        seq[partner_cut - 1:partner_cut + 3] = 'CCGG'
        seq = ''.join(seq)
        if seq.count('GAATTC') == seq.count('CCGG') == 1:
            return seq


class DoubleDigestTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.variants = self.root / 'snps.tsv'
        self.variants.write_text('chrom\tpos\tref\talt\n1\t200\tC\tG\n')

    def search(self, partner_cut=50, shared_cut=100, starts=(0,), duplicate_start=None, **options):
        genome = self.root / f'{partner_cut}.fa'
        seq = template(partner_cut, shared_cut)
        contigs = f'>1\n{seq}\n'
        if duplicate_start is not None:
            contigs += f'>2\n{seq[duplicate_start:]}\n'
        genome.write_text(contigs)
        params = dict(core.DEFAULT_PARAMS, fasta=str(genome), input_file=str(self.variants),
                      enzyme_mode='only', lab_enzymes=['EcoRI', 'MspI'], dcaps='off',
                      specificity='off', top_results=0, group_isoschizomers=False)
        params.update(options)
        # Fix the amplicon, exercise the normal validation of Primer3 coordinates.
        result = {'PRIMER_PAIR_NUM_RETURNED': len(starts)}
        for index, start in enumerate(starts):
            result.update({f'PRIMER_LEFT_{index}': (start, 20), f'PRIMER_RIGHT_{index}': (299, 20),
                           f'PRIMER_LEFT_{index}_SEQUENCE': seq[start:start + 20],
                           f'PRIMER_RIGHT_{index}_SEQUENCE': str(Seq(seq[-20:]).reverse_complement()),
                           f'PRIMER_LEFT_{index}_TM': 60, f'PRIMER_RIGHT_{index}_TM': 60,
                           f'PRIMER_PAIR_{index}_PRODUCT_SIZE': 300 - start})
        self.report = {}
        with patch.object(core, 'primer3_pick', return_value=result):
            header, rows = core.run_analysis(params, report=self.report)
        return [dict(zip(header, row)) for row in rows]

    def test_common_cut_rescues_both_ambiguous_singles(self):
        for enzyme in ['EcoRI', 'MspI']:
            self.assertEqual(self.search(lab_enzymes=[enzyme])[0]['status'], 'no_discriminating_assay')
        row, = self.search()
        self.assertEqual(row['design'], 'double_digest')
        self.assertEqual((row['enzyme_1'], row['enzyme_2']), ('EcoRI', 'MspI'))
        self.assertEqual((row['frags_ref'], row['frags_alt']), ('50+50+200', '50+50+100+100'))
        self.assertEqual((row['bands_ref_ref'], row['bands_ref_alt'], row['bands_alt_alt']),
                         ('200', '100; 200', '100'))
        self.assertEqual(row['genotype_quality'], 'robust')
        self.assertEqual(row['reaction_compatibility'], 'unconfirmed')
        self.assertEqual(self.search(double_digest=False)[0]['status'], 'no_discriminating_assay')

    def test_common_cut_leaving_rr_and_ra_equal_is_rejected(self):
        row, = self.search(150)
        self.assertEqual(row['status'], 'no_discriminating_assay')
        self.assertEqual(self.report['summary']['double_pairs_scored'], 1)

    def test_limits_and_background_on_partner(self):
        self.assertEqual(self.search(max_cuts=2)[0]['status'], 'no_discriminating_assay')
        vcf = self.root / 'background.vcf'
        vcf.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                       '1\t50\t.\tC\tA\t.\tPASS\t.\n')
        self.assertEqual(self.search(background_vcf=str(vcf))[0]['status'], 'no_discriminating_assay')
        row, = self.search(background_vcf=str(vcf), exclude_background_sites=False)
        self.assertIn('50 C>A', row['background_site_variants'])

    def test_card_sheet_and_map_use_both_enzymes(self):
        import os
        os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
        from PyQt6.QtWidgets import QApplication
        from lab_sheet import sheet_fields
        from enzyme_catalog import property_lines
        from snp_summary import SnpSummary, read_input_records
        from visualization import PrimerGraphicPanel, GelPanel, payload_for_scheme, render_payload
        from tests.helpers import dispose
        import i18n
        app = QApplication.instance() or QApplication([])
        self.addCleanup(i18n.set_language, i18n.language())
        row, = self.search()
        row['enzyme_product'] = 'R3101'
        row['enzyme_product_2'] = 'R0106S'
        summary = SnpSummary(read_input_records(str(self.variants))[0][0], schemes=[row])
        payload, message = payload_for_scheme(row, str(self.root / '50.fa'), {})
        panel, gel = PrimerGraphicPanel(), GelPanel()
        self.addCleanup(dispose, panel)
        self.addCleanup(dispose, gel)
        render_payload(panel, gel, payload, message)
        self.assertEqual(panel._cuts_ref, [0, 50, 100, 300])
        self.assertEqual(panel._cuts_alt, [0, 50, 100, 200, 300])
        for lang, warning in [('ru', 'Совместимость условий не подтверждена'),
                              ('en', 'Reaction compatibility unconfirmed')]:
            i18n.set_language(lang)
            self.assertIn(warning, summary.note('off'))
            supplier_row = dict(sheet_fields(summary, 'off'))[i18n.t('Поставщики')]
            self.assertIn('EcoRI:', supplier_row)
            self.assertIn('MspI:', supplier_row)
            text = '\n'.join(property_lines(row, panel.allele_sequences()))
            self.assertIn('R3101', text)
            self.assertIn('R0106S', text)
        app.processEvents()

    def test_stage_order_preserves_earlier_nominal_and_explicit_dcaps(self):
        rows = self.search(shared_cut=85, specificity='exclude', dcaps='fallback')
        self.assertEqual([(r['design'], r['genotype_quality']) for r in rows],
                         [('double_digest', 'robust'), ('natural', 'nominal_only')])
        self.assertEqual(self.report['summary']['dcaps_snps'], 0)
        self.assertEqual(self.report['summary']['specificity_scans'], 1)
        self.search(shared_cut=85, dcaps='always')
        self.assertEqual(self.report['summary']['dcaps_snps'], 1)
        rows = self.search(shared_cut=85, lab_enzymes=['EcoRI'], dcaps='fallback')
        self.assertEqual([r['genotype_quality'] for r in rows], ['nominal_only'])
        self.assertEqual(self.report['summary']['dcaps_snps'], 1)

    def test_specificity_reopens_search_on_another_amplicon(self):
        settings = dict(starts=(20, 0), duplicate_start=20, dcaps='fallback')
        rows = self.search(**settings, specificity='warn')
        self.assertEqual({r['design'] for r in rows}, {'natural'})
        self.assertEqual(self.report['summary']['double_snps'], 0)
        rows = self.search(**settings, specificity='exclude')
        self.assertEqual([(r['design'], r['amplicon_start'], r['offtarget_products']) for r in rows],
                         [('double_digest', 1, 0)])
        self.assertEqual(self.report['summary']['dcaps_snps'], 0)
        self.assertGreater(self.report['summary']['excluded_offtarget'], 0)

    def test_pair_membership_before_grouping_and_top_one(self):
        with patch.object(core, 'AllEnzymes', [EcoRI, HpaII, MspI]):
            for grouped in (True, False):
                row, = self.search(enzyme_mode='prefer', top_results=1, group_isoschizomers=grouped)
                self.assertEqual(row['enzyme'], 'EcoRI + MspI')
            row, = self.search(enzyme_mode='prefer', lab_enzymes=['MspI'], top_results=1)
            self.assertEqual(row['enzyme'], 'EcoRI + HpaII')
            rows = self.search(enzyme_mode='only')
            self.assertEqual([r['enzyme'] for r in rows], ['EcoRI + MspI'])

    def test_failed_partner_search_is_not_an_uncut_allele(self):
        original = core.cut_positions

        def cuts(enzyme, seq):
            if enzyme is MspI:
                raise ValueError('unavailable coordinates')
            return original(enzyme, seq)

        with patch.object(core, 'cut_positions', cuts):
            row, = self.search()
        self.assertEqual(row['status'], 'no_discriminating_assay')
        self.assertEqual(self.report['summary']['enzyme_errors'], 1)

    def test_unresolved_enzyme_cannot_be_a_partner(self):
        from Bio.Restriction import AllEnzymes
        unresolved = next(e for e in AllEnzymes if e.is_unknown())
        rows = self.search(lab_enzymes=['EcoRI', 'MspI', str(unresolved)])
        self.assertEqual([r['enzyme'] for r in rows], ['EcoRI + MspI'])
        self.assertEqual(self.report['summary']['double_pairs_possible'], 1)

    def test_partner_site_in_primer(self):
        row, = self.search(partner_cut=10)
        self.assertEqual(row['status'], 'no_discriminating_assay')
        # A 90 bp band is too close to 100 under stress; expose the filter at finer resolution.
        rows = self.search(partner_cut=10, exclude_primer_sites=False, gel_resolution_bp=3, gel_resolution_pct=0)
        self.assertEqual(rows[0]['primer_site'], 'MspI: left')
        row, = self.search(partner_cut=10, gel_resolution_bp=3, gel_resolution_pct=0)
        self.assertEqual(row['status'], 'no_discriminating_assay')
        self.assertEqual(self.report['summary']['excluded_primer_site'], 1)

    def test_removed_robust_double_digest_reopens_dcaps(self):
        # A second copy of the complete amplicon gives each primer pair a second product.
        row, = self.search(specificity='exclude', dcaps='fallback', duplicate_start=0)
        self.assertEqual(row['status'], 'no_specific_primers')
        self.assertEqual(self.report['summary']['double_snps'], 1)
        self.assertEqual(self.report['summary']['dcaps_snps'], 1)

    def test_exhaustive_pairs_agree_before_grouping(self):
        # Independent oracle: search both complete alleles for EVERY pair, with
        # no SNP/anchor filter or cached allele result. Compare coordinates via
        # the ordered fragments, then the full genotype band pattern and class.
        enzymes = [EcoRI, ApoI, MspI, HpaII, MseI, SmaI, XmaI, BsaI]
        params = dict(core.DEFAULT_PARAMS, delta=1, min_frag=10, max_cuts=10)
        model = GelModel(min_visible=10)
        run = core._Run(params, model, {}, Counter(), enzymes, enzymes, False, 'window', 'off')
        seq = list(template())
        seq[0:6] = 'GGTCTC'
        seq[30:36] = 'CCCGGG'
        seq = ''.join(seq)
        cases = [(199, 'G'), (0, 'A'), (30, 'A'), (101, 'C'), (299, 'A' if seq[-1] != 'A' else 'C')]
        coverage = set()
        for off, new_base in cases:
            alt = seq[:off] + new_base + seq[off + 1:]
            window = core._Window(seq[off], new_base, off + 1, 1, 300, seq, seq, off, {}, altseq=alt)
            for lo, hi in [(0, 300), (max(0, min(off, 35)), max(250, off + 1))]:
                observed = core._double_hits(run, window, [dict(start0=lo, end0=hi)])
                keys = ('frags_ref', 'frags_alt', 'bands_ref_ref', 'bands_ref_alt', 'bands_alt_alt',
                        'worst_bands_ref_ref', 'worst_bands_ref_alt', 'worst_bands_alt_alt', 'genotype_quality')
                actual = {(r['enzyme_1'], r['enzyme_2']): tuple(r[k] for k in keys) for r in observed}
                self.assertEqual(len(observed), len(actual))
                expected = {}
                for e1, e2 in itertools.combinations(enzymes, 2):
                    maps = [tuple(set(core.cut_positions(e, allele[lo:hi])) for allele in (seq, alt)) for e in (e1, e2)]
                    combined = tuple(a | b for a, b in zip(*maps))
                    if combined in maps or max(map(len, combined)) - 2 > params['max_cuts']:
                        continue
                    fr, fa = [core.frag_lengths(sorted(cuts)) for cuts in combined]
                    scores = score_genotypes(fr, fa, model)
                    if scores['genotype_quality'] == 'ambiguous':
                        continue
                    # With integer lengths and threshold 1, Q > 1 already implies Δ >= 1.
                    pair = tuple(sorted((str(e1), str(e2))))
                    expected[pair] = (core.frags_to_str(fr), core.frags_to_str(fa), *(scores[k] for k in keys[2:]))
                    if all(cr != ca for cr, ca in maps):
                        coverage.add('two differential enzymes')
                    coverage.update(pair)
                self.assertEqual(actual, expected, (off, lo, hi))
        self.assertTrue({'BsaI', 'two differential enzymes'} <= coverage, coverage)
