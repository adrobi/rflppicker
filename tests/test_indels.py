import random
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from Bio.Seq import Seq

import rflp_core as core
from tests.test_double_digest import template


class IndelTests(unittest.TestCase):
    def search(self, ref, alt, pos=101, seq=None, starts=(0,), duplicate_start=None, **options):
        if seq is None:
            rng = random.Random(12)
            seq = ''.join(rng.choices('ACGT', k=300))
            seq = seq[:100] + 'GCAATTC' + seq[107:]
        primers = {'PRIMER_PAIR_NUM_RETURNED': len(starts)}
        for i, start in enumerate(starts):
            primers.update({f'PRIMER_LEFT_{i}': (start, 20), f'PRIMER_RIGHT_{i}': (299, 20),
                            f'PRIMER_LEFT_{i}_SEQUENCE': seq[start:start + 20],
                            f'PRIMER_RIGHT_{i}_SEQUENCE': str(Seq(seq[-20:]).reverse_complement()),
                            f'PRIMER_LEFT_{i}_TM': 60, f'PRIMER_RIGHT_{i}_TM': 60,
                            f'PRIMER_PAIR_{i}_PRODUCT_SIZE': 300 - start})
        work = tempfile.TemporaryDirectory()
        self.addCleanup(work.cleanup)
        self.fasta, self.variants = Path(work.name) / 'ref.fa', Path(work.name) / 'variants.tsv'
        contigs = f'>1\n{seq}\n'
        if duplicate_start is not None:
            contigs += f'>2\n{seq[duplicate_start:]}\n'
        self.fasta.write_text(contigs)
        self.variants.write_text(f'chrom\tpos\tref\talt\n1\t{pos}\t{ref}\t{alt}\n')
        params = dict(core.DEFAULT_PARAMS, fasta=str(self.fasta), input_file=str(self.variants),
                      enzyme_mode='only', lab_enzymes=['EcoRI'], specificity='off')
        params.update(options)
        report = {}
        with patch.object(core, 'primer3_pick', return_value=primers):
            header, rows = core.run_analysis(params, report=report)
        return [dict(zip(header, row)) for row in rows], report

    def test_long_insertion_distinguishes_three_genotypes_by_product_length(self):
        rows, _ = self.search('G', 'G' + 'T' * 40, dcaps='always')
        self.assertEqual([row['design'] for row in rows], ['amplicon_length'])
        row, = rows
        self.assertEqual((row['product_size'], row['product_size_alt']), (300, 340))
        self.assertEqual((row['bands_ref_ref'], row['bands_ref_alt'], row['bands_alt_alt']),
                         ('300', '300; 340', '340'))

    def test_short_deletion_creates_a_site_and_falls_back_to_restriction(self):
        rows, _ = self.search('GC', 'G', skip_norm=True)
        self.assertEqual([(r['design'], r['enzyme']) for r in rows], [('natural', 'EcoRI')])
        row, = rows
        self.assertEqual((row['product_size'], row['product_size_alt']), (300, 299))
        self.assertEqual((row['frags_ref'], row['frags_alt']), ('300', '101+198'))

    def test_long_deletion_distinguishes_three_genotypes_by_product_length(self):
        rng = random.Random(12)
        seq = ''.join(rng.choices('ACGT', k=300))
        rows, _ = self.search(seq[100:141], seq[100], seq=seq)
        row, = rows
        self.assertEqual(row['design'], 'amplicon_length')
        self.assertEqual((row['bands_ref_ref'], row['bands_ref_alt'], row['bands_alt_alt']),
                         ('300', '260; 300', '260'))

    def test_specificity_exclusion_reopens_digest_search(self):
        rng = random.Random(12)
        seq = ''.join(rng.choices('ACGT', k=300))
        seq = seq[:100] + 'G' + 'C' * 20 + 'AATTC' + seq[126:]
        options = dict(seq=seq, starts=(0, 50), duplicate_start=50, dcaps='always',
                       gel_resolution_pct=5, delta=1, prod_min=200, snp_end_distance=40, top_results=0)
        rows, _ = self.search('G' + 'C' * 20, 'G', specificity='warn', **options)
        self.assertEqual({(r['design'], r['genotype_quality']) for r in rows},
                         {('amplicon_length', 'robust'), ('amplicon_length', 'nominal_only')})
        rows, _ = self.search('G' + 'C' * 20, 'G', specificity='exclude', **options)
        self.assertEqual([(r['design'], r['genotype_quality'], r['offtarget_products']) for r in rows],
                         [('natural', 'robust', 0), ('amplicon_length', 'nominal_only', 0)])

    def test_reference_mismatch_after_anchor_is_rejected(self):
        for skip_norm in (False, True):
            with self.subTest(skip_norm=skip_norm):
                rows, report = self.search('GA', 'G', skip_norm=skip_norm)
                self.assertEqual(rows, [])
                self.assertEqual([r['reason'] for r in report['rejected']], ['ref_not_match'])

    def test_short_deletion_can_require_two_enzymes(self):
        seq = template()
        seq = seq[:198] + 'G' + seq[199:]
        options = dict(pos=199, seq=seq, lab_enzymes=['EcoRI', 'MspI'])
        rows, _ = self.search('GC', 'G', double_digest=False, **options)
        self.assertEqual([r['status'] for r in rows], ['no_discriminating_assay'])
        rows, _ = self.search('GC', 'G', **options)
        row, = rows
        self.assertEqual((row['design'], row['enzyme']), ('double_digest', 'EcoRI + MspI'))
        self.assertEqual((row['frags_ref'], row['frags_alt']), ('50+50+200', '50+50+99+100'))

    def test_length_scheme_survives_summary_and_export(self):
        from snp_summary import build_summaries, read_input_records, summary_table

        rows, report = self.search('G', 'G' + 'T' * 40)
        header = list(rows[0])
        records, _ = read_input_records(str(self.variants))
        summaries = build_summaries(records, header, [list(rows[0].values())], report.get('rejected', []))
        self.assertEqual(summaries[0].status, 'robust')
        titles, exported = summary_table(summaries, 'off')
        from i18n import t
        self.assertEqual(exported[0][titles.index(t('Продукт ALT, п.н.'))], 340)

    def test_map_renders_all_bases_of_the_shorter_allele(self):
        from PyQt6.QtWidgets import QApplication
        from tests.helpers import dispose
        from visualization import GelPanel, PrimerGraphicPanel, payload_for_scheme, render_payload

        app = QApplication.instance() or QApplication([])
        rows, _ = self.search('GC', 'G')
        panel, gel = PrimerGraphicPanel(), GelPanel()
        self.addCleanup(dispose, panel)
        self.addCleanup(dispose, gel)
        payload, message = payload_for_scheme(rows[0], str(self.fasta), {})
        render_payload(panel, gel, payload, message)
        self.assertIsNotNone(panel.allele_sequences()[1])
        self.assertEqual(tuple(map(len, panel.allele_sequences())), (300, 299))
        self.assertEqual(panel._cuts_alt, [0, 101, 299])
        app.processEvents()
