import random
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from Bio.Restriction import EcoRI
from Bio.Seq import Seq
import rflp_core as core
from reporting import export_csv, export_excel


class JointDesignTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        rng = random.Random(19)
        sequence = ''.join(rng.choice('ACGT') for _ in range(400))
        self.seq = sequence[:197] + 'GAATTC' + sequence[203:]
        self.assertEqual(self.seq.count('GAATTC'), 1)
        self.fasta = self.root / 'toy.fa'
        self.fasta.write_text('>1\n' + self.seq + '\n')
        variants = self.root / 'variants.tsv'
        variants.write_text('chrom\tpos\tref\talt\n1\t200\tA\tG\n')
        self.params = dict(fasta=str(self.fasta), input_file=str(variants), names2_file='',
                           flank=500, min_frag=80, max_frag=800, delta=25,
                           gain_loss_only=False, max_cuts=3, suppliers=[], use_primer3=True,
                           prod_min=200, prod_max=500, tm_min=50, tm_max=70,
                           no_save_norm=True, joint_design=True, primer_pairs=2, top_results=0)

    def primer_result(self, starts=(150, 0)):
        data = {'PRIMER_PAIR_NUM_RETURNED': len(starts)}
        for i, start in enumerate(starts):
            data.update({f'PRIMER_LEFT_{i}': [start, 20], f'PRIMER_RIGHT_{i}': [399, 20],
                         f'PRIMER_LEFT_{i}_SEQUENCE': self.seq[start:start + 20],
                         f'PRIMER_RIGHT_{i}_SEQUENCE': str(Seq(self.seq[380:]).reverse_complement()),
                         f'PRIMER_LEFT_{i}_TM': 60, f'PRIMER_RIGHT_{i}_TM': 60,
                         f'PRIMER_PAIR_{i}_PRODUCT_SIZE': 400 - start,
                         f'PRIMER_PAIR_{i}_PENALTY': i + 1})
        return data

    def run_design(self, pres=None, **overrides):
        report = {}
        with patch.object(core, 'AllEnzymes', [EcoRI]), patch.object(core, 'primer3_pick', return_value=pres or self.primer_result()):
            header, rows = core.run_rflp_gui_mode(dict(self.params, **overrides), report=report)
        self.assertTrue(all(len(row) == len(header) for row in rows))
        return header, rows, report

    def test_second_pair_wins_with_real_restriction_coordinates_and_export(self):
        header, rows, report = self.run_design()
        first, second = [dict(zip(header, row)) for row in rows]
        self.assertEqual(first['primer_pair_index'], 2)
        self.assertEqual(first['frags_ref'], '198+202')
        self.assertEqual(first['frags_alt'], '400')
        self.assertEqual((first['amplicon_start'], first['amplicon_end']), (1, 400))
        self.assertEqual(second['hidden_ref_count'], 1)
        self.assertGreater(first['worst_margin'], second['worst_margin'])
        self.assertEqual(report['summary']['pairs_evaluated'], 2)
        export_csv(self.root / 'results.csv', header, rows)
        export_excel(self.root / 'results.xlsx', header, rows, metadata=report)

    def test_top_limit_and_deduplication(self):
        header, rows, report = self.run_design(top_results=1)
        self.assertEqual(len(rows), 1)
        self.assertEqual(report['summary']['candidates'], 2)
        header, rows, report = self.run_design(self.primer_result((0, 0)))
        self.assertEqual(report['summary']['pairs_evaluated'], 1)
        self.assertEqual(len(rows), 1)

    def test_primer_overlap_rejected(self):
        header, rows, _ = self.run_design(self.primer_result((190,)))
        self.assertEqual(dict(zip(header, rows[0]))['status'], 'primer3_error')

    def test_unavailable_primer3_does_not_return_window_hits(self):
        with patch.object(core, '_PRIMER3_AVAILABLE', False):
            header, rows = core.run_rflp_gui_mode(self.params)
        row = dict(zip(header, rows[0]))
        self.assertEqual(row['status'], 'primer3_error')
        self.assertFalse(row['enzyme'])

    def test_no_primers_status_and_empty_input_schema(self):
        with patch.object(core, 'primer3_pick', return_value=None) as mocked:
            mocked.last_error = ''
            header, rows = core.run_rflp_gui_mode(self.params)
        self.assertEqual(dict(zip(header, rows[0]))['status'], 'no_primers')
        Path(self.params['input_file']).write_text('chrom\tpos\tref\talt\n')
        header, rows = core.run_rflp_gui_mode(self.params)
        self.assertIn('worst_margin', header)
        self.assertEqual(rows, [])

    def test_enzyme_error_is_not_treated_as_uncut(self):
        with patch.object(core, 'cut_positions', side_effect=RuntimeError('broken search')):
            header, rows, report = self.run_design()
        self.assertFalse(dict(zip(header, rows[0]))['enzyme'])
        self.assertEqual(report['summary']['enzyme_errors'], 2)

    def test_cancellation_inside_pair_search_closes_fasta(self):
        original_cut = core.cut_positions
        cut_calls = []
        def cut(*args):
            cut_calls.append(True)
            return original_cut(*args)
        with patch.object(core, 'AllEnzymes', [EcoRI]), patch.object(core, 'primer3_pick', return_value=self.primer_result()), patch.object(core, 'cut_positions', side_effect=cut):
            with self.assertRaises(core.AnalysisCancelled):
                core.run_rflp_gui_mode(self.params, cancel_cb=lambda: bool(cut_calls))
        # Windows refuses to replace a file while an ordinary reader holds it.
        self.fasta.replace(self.root / 'closed.fa')

    def test_primer3_receives_requested_number_of_pairs(self):
        with patch.object(core, '_PRIMER3_AVAILABLE', True), patch.object(core, '_PRIMER3_USE_NEW', True), patch.object(core.primer3.bindings, 'design_primers', return_value=self.primer_result()) as design:
            core.primer3_pick(self.seq, 199, 200, 500, 58, 62, num_return=12)
        self.assertEqual(design.call_args.args[1]['PRIMER_NUM_RETURN'], 12)
        self.assertEqual(design.call_args.args[0]['SEQUENCE_TARGET'], [199, 1])


if __name__ == '__main__':
    unittest.main()
