import random
import re
import tempfile
import unittest
from pathlib import Path

from Bio.Restriction import AllEnzymes

import rflp_core
from tests.helpers import write_text
from rflp_core import (
    AnalysisCancelled,
    cut_positions,
    diagnostic_delta,
    normalize_variants_for_rflp,
    parse_input,
    run_rflp_gui_mode,
    snp_site_enzymes,
)


class CoreTests(unittest.TestCase):
    def test_csv_quotes_and_rejected_multibase(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            fasta = root / 'toy.fa'
            fasta.write_text('>1\nACGTACGTACGT\n', encoding='utf-8')
            variants = root / 'variants.csv'
            variants.write_text('chrom,pos,ref,alt\n1,2,"C","T"\n1,3,GT,AA\n', encoding='utf-8')
            report = {}
            parsed = parse_input(str(variants), report=report)
            self.assertEqual(parsed[0], ('1', 2, 'C', 'T'))
            cleaned, _, _, _ = normalize_variants_for_rflp(str(fasta), str(variants), None, save_files=False, report=report)
            self.assertEqual(cleaned, [('1', 2, 'C', 'T')])
            self.assertTrue(any(x['reason'] == 'unsupported_variant' for x in report['rejected']))

    def test_skip_norm_still_checks_ref(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            fasta = root / 'toy.fa'; fasta.write_text('>1\nACGTACGT\n', encoding='utf-8')
            variants = root / 'v.tsv'; variants.write_text('chrom\tpos\tref\talt\n1\t2\tA\tG\n', encoding='utf-8')
            report = {}
            cleaned, _, _, _ = normalize_variants_for_rflp(str(fasta), str(variants), None, skip_norm=True, save_files=False, report=report)
            self.assertEqual(cleaned, [])
            self.assertEqual(report['rejected'][0]['reason'], 'ref_not_match')

    def test_threshold_is_applied_for_gain_loss_and_reported(self):
        self.assertEqual(diagnostic_delta([500], [100, 400]), 400)
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            fasta = root / 'toy.fa'; fasta.write_text('>1\nACGTACGTACGTACGTACGT\n', encoding='utf-8')
            variants = root / 'v.tsv'; variants.write_text('chrom\tpos\tref\talt\n1\t2\tC\tT\n', encoding='utf-8')
            params = dict(fasta=str(fasta), input_file=str(variants), names2_file='', flank=4,
                          min_frag=1, max_frag=1000, delta=1, gain_loss_only=False, max_cuts=3,
                          suppliers=[], use_primer3=False, prod_min=50, prod_max=100,
                          tm_min=50, tm_max=70, skip_norm=False, no_save_norm=True)
            header, rows = run_rflp_gui_mode(params)
            self.assertIn('status', header); self.assertIn('delta_threshold_bp', header)
            self.assertTrue(rows)
            self.assertEqual(len(header), len(rows[0]))

    def test_cancellation_is_honored(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            fasta = root / 'toy.fa'; fasta.write_text('>1\nACGTACGTACGT\n', encoding='utf-8')
            variants = root / 'v.tsv'; variants.write_text('chrom\tpos\tref\talt\n1\t2\tC\tT\n', encoding='utf-8')
            params = dict(fasta=str(fasta), input_file=str(variants), names2_file='', flank=4,
                          min_frag=1, max_frag=1000, delta=1, gain_loss_only=False, max_cuts=3,
                          suppliers=[], use_primer3=False, prod_min=50, prod_max=100,
                          tm_min=50, tm_max=70, skip_norm=False, no_save_norm=True)
            with self.assertRaises(AnalysisCancelled):
                run_rflp_gui_mode(params, cancel_cb=lambda: True)

    def test_site_patterns_have_fixed_width(self):
        # snp_site_enzymes relies on every site being a look-ahead of enz.size bases.
        lookahead = re.compile(r'\(\?=\(\?P<\w+>((?:\[[A-Z]+\]|[A-Z.])+)\)\)')
        for enz in AllEnzymes:
            parts = enz.compsite.pattern.split('|')
            self.assertIn(len(parts), (1, 2), enz)
            for part in parts:
                match = lookahead.fullmatch(part)
                self.assertIsNotNone(match, (enz, part))
                self.assertEqual(len(re.findall(r'\[[A-Z]+\]|[A-Z.]', match.group(1))), enz.size, enz)

    def test_enzymes_without_snp_site_cut_both_alleles_equally(self):
        rng = random.Random(7)
        enzymes = sorted(AllEnzymes, key=lambda e: e.__name__)
        for _ in range(6):
            seq = ''.join(rng.choice('ACGT') for _ in range(90))
            off = rng.choice([0, 1, 5, 45, 84, 89])
            altseq = seq[:off] + rng.choice([b for b in 'ACGT' if b != seq[off]]) + seq[off + 1:]
            found = set(snp_site_enzymes(enzymes, seq, altseq, off))
            self.assertTrue(found)
            amplicons = {(0, len(seq)), (off, off + 1), (max(0, off - 1), min(len(seq), off + 2)),
                         (rng.randint(0, off), rng.randint(off + 1, len(seq)))}
            for enz in enzymes:
                if enz in found:
                    continue
                for a, b in amplicons:
                    self.assertEqual(cut_positions(enz, seq[a:b]), cut_positions(enz, altseq[a:b]),
                                     (enz, seq, off, a, b))


class InputTableTests(unittest.TestCase):
    def setUp(self):
        self._dir = tempfile.TemporaryDirectory()
        self.tmp = self._dir.name

    def tearDown(self):
        self._dir.cleanup()

    def test_read_input_table_detects_columns_and_id(self):
        path = write_text(self.tmp, 'v.tsv', 'chrom\tpos\tref\talt\tName\n1\t10\tA\tG\tDGAT1\n\n1\t11\tC\n')
        columns, lines, delim = rflp_core.read_input_table(path)
        self.assertEqual(delim, '\t')
        self.assertEqual(columns, {'chrom': 0, 'pos': 1, 'ref': 2, 'alt': 3, 'id': 4})
        self.assertEqual(lines, [(2, ['1', '10', 'A', 'G', 'DGAT1']), (4, ['1', '11', 'C'])])

    def test_read_input_table_headerless_and_template_rows(self):
        path = write_text(self.tmp, 'v.csv', '1,10,A,G\n2,20,C,T\n')
        columns, lines, _ = rflp_core.read_input_table(path)
        self.assertEqual(set(columns.values()), {None})
        self.assertEqual([n for n, _ in lines], [1, 2])
        path = write_text(self.tmp, 'h.csv', 'chrom,pos,ref,alt\nдолжно быть: 1,10,A,G\n1,10,A,G\n')
        self.assertEqual([n for n, _ in rflp_core.read_input_table(path)[1]], [3])

    def test_read_input_table_missing_file(self):
        self.assertEqual(rflp_core.read_input_table(str(Path(self.tmp) / 'missing.tsv')), ({}, [], ''))

    def test_parse_fields(self):
        columns = {'chrom': None, 'pos': None, 'ref': None, 'alt': None, 'id': None}
        self.assertEqual(rflp_core.parse_fields(['chr1 ', '10', 'a', 'g '], columns), ('1', 10, 'A', 'G'))
        with self.assertRaises(ValueError):
            rflp_core.parse_fields(['1', 'x', 'A', 'G'], columns)


if __name__ == '__main__':
    unittest.main()
