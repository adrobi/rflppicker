import csv
import tempfile
import unittest
from pathlib import Path

from rflp_core import (
    AnalysisCancelled,
    diagnostic_delta,
    normalize_variants_for_rflp,
    parse_input,
    run_rflp_gui_mode,
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


if __name__ == '__main__':
    unittest.main()
