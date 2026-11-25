import json
import tempfile
import unittest
from pathlib import Path

import assay_scoring as scoring
import rflp_core as core


class MultiallelicTests(unittest.TestCase):
    def search(self, pos, ref, alt, **options):
        work = tempfile.TemporaryDirectory()
        self.addCleanup(work.cleanup)
        self.fasta, self.variants = Path(work.name) / 'toy.fa', Path(work.name) / 'alleles.tsv'
        self.fasta.write_text('>1\n' + 'A' * 300 + '\n')
        self.variants.write_text(f'chrom\tpos\tref\talt\n1\t{pos}\t{ref}\t{alt}\n')
        params = dict(core.DEFAULT_PARAMS, fasta=str(self.fasta), input_file=str(self.variants),
                      use_primer3=False, top_results=1, enzyme_mode='only', lab_enzymes=['EcoRI', 'MspI'])
        params.update(options)
        report = {}
        header, rows = core.run_analysis(params, report=report)
        return header, rows, report

    def test_alt_alt_collision_rejects_locus_even_when_each_alt_works_with_ref(self):
        fragments = [[600], [200, 400], [200, 200, 200]]
        model = scoring.GelModel()
        for alt in fragments[1:]:
            self.assertEqual(scoring.score_genotypes(fragments[0], alt, model)['genotype_quality'], 'robust')
        result = scoring.score_locus(fragments, model)
        self.assertEqual(result['genotype_quality'], 'ambiguous')
        self.assertEqual(result['genotype_bands']['ALT1/ALT1'], result['genotype_bands']['ALT1/ALT2'])

    def test_window_search_returns_one_scheme_for_six_genotypes(self):
        header, rows, _ = self.search(150, 'A', f'A{"G" * 60}, A{"C" * 130}')
        self.assertEqual(len(rows), 1)
        row = dict(zip(header, rows[0]))
        self.assertEqual((row['design'], row['genotype_quality']), ('amplicon_length', 'robust'))
        self.assertEqual(json.loads(row['genotype_bands']), {
            'REF/REF': '300', 'REF/ALT1': '300; 360', 'REF/ALT2': '300; 430',
            'ALT1/ALT1': '360', 'ALT1/ALT2': '360; 430', 'ALT2/ALT2': '430',
        })
        self.assertEqual([x['length'] for x in json.loads(row['allele_products'])], [300, 360, 430])

    def test_two_allele_specific_cuts_together_identify_the_locus(self):
        alt1, alt2 = 'GAATTC' + 'A' * 38, 'A' * 40 + 'CCGG'
        for enabled, expected in ((False, 'no_discriminating_assay'), (True, 'window_only')):
            header, rows, _ = self.search(100, 'A' * 44, f'{alt1},{alt2}', double_digest=enabled)
            row = dict(zip(header, rows[0]))
            self.assertEqual(row['status'], expected)
        self.assertEqual((row['design'], row['enzyme']), ('double_digest', 'EcoRI + MspI'))
        self.assertEqual(json.loads(row['genotype_bands'])['ALT1/ALT2'], '100; 140; 160; 200')

    def test_one_enzyme_can_identify_all_genotypes(self):
        alt1, alt2 = 'GAATTC' + 'A' * 40, 'A' * 40 + 'GAATTC'
        header, rows, _ = self.search(100, 'A' * 46, f'{alt1},{alt2}')
        row = dict(zip(header, rows[0]))
        self.assertEqual((row['design'], row['enzyme']), ('natural', 'EcoRI'))
        self.assertEqual(json.loads(row['genotype_bands'])['REF/ALT2'], '140; 160; 300')

    def test_duplicate_and_symbolic_alternates_are_rejected(self):
        for alt in ('A,C', 'C,C', 'C,<DEL>'):
            with self.subTest(alt=alt):
                _, rows, report = self.search(150, 'A', alt)
                self.assertEqual(rows, [])
                self.assertEqual([r['reason'] for r in report['rejected']], ['unsupported_variant'])

    def test_allele_numbering_supports_four_alleles(self):
        result = scoring.score_locus([[180], [260], [340], [430]], scoring.GelModel())
        self.assertEqual(result['genotype_quality'], 'robust')
        self.assertEqual(len(result['genotype_bands']), 10)
        self.assertEqual(result['genotype_bands']['ALT2/ALT3'], '340; 430')

    def test_length_scheme_shows_every_allele_and_genotype(self):
        from PyQt6.QtWidgets import QApplication
        from tests.helpers import dispose
        from visualization import GelPanel, PrimerGraphicPanel, payload_for_scheme, render_payload

        app = QApplication.instance() or QApplication([])
        header, rows, _ = self.search(150, 'A', f'A{"G" * 60}, A{"C" * 130}')
        panel, gel = PrimerGraphicPanel(), GelPanel()
        self.addCleanup(dispose, panel)
        self.addCleanup(dispose, gel)
        payload, message = payload_for_scheme(dict(zip(header, rows[0])), str(self.fasta), {})
        render_payload(panel, gel, payload, message)
        self.assertEqual(tuple(map(len, panel.allele_sequences())), (300, 360, 430))
        self.assertEqual([label for label, _ in gel._lanes],
                         ['REF/REF', 'REF/ALT1', 'REF/ALT2', 'ALT1/ALT1', 'ALT1/ALT2', 'ALT2/ALT2'])
        app.processEvents()
