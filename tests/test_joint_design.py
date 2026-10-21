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
        self.assertNotIn('PRIMER_MIN_SIZE', design.call_args.args[1])
        with patch.object(core, '_PRIMER3_AVAILABLE', True), patch.object(core, '_PRIMER3_USE_NEW', True), patch.object(core.primer3.bindings, 'design_primers', return_value=self.primer_result()) as design:
            core.primer3_pick(self.seq, 199, 200, 500, 58, 62, target_flank=62)
        self.assertEqual(design.call_args.args[0]['SEQUENCE_TARGET'], [137, 125])
        self.assertEqual(design.call_args.args[1]['PRIMER_MIN_SIZE'], core.PRIMER_MIN_SIZE)

    def test_site_inside_primer_is_flagged_and_can_be_excluded(self):
        # A second EcoRI site inside the left primer of pair 1 (bases 150-169).
        self.seq = self.seq[:152] + 'GAATTC' + self.seq[158:]
        self.fasta.write_text('>1\n' + self.seq + '\n')
        header, rows, report = self.run_design()
        by_pair = {r['primer_pair_index']: r for r in (dict(zip(header, row)) for row in rows)}
        self.assertEqual((by_pair[1]['primer_site'], by_pair[2]['primer_site']), ('left', ''))
        header, rows, report = self.run_design(exclude_primer_sites=True)
        self.assertEqual([dict(zip(header, row))['primer_pair_index'] for row in rows], [2])
        self.assertEqual(report['summary']['excluded_primer_site'], 1)

    def test_primer_site_label_searches_the_right_binding_site(self):
        amplicon = 'GAATTCAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCGGGAATTCTT'
        self.assertEqual(core.primer_site_label(EcoRI, amplicon, 20, 10), 'left; right')
        self.assertEqual(core.primer_site_label(EcoRI, amplicon, 5, 7), '')
        self.assertEqual(core.primer_site_label(EcoRI, amplicon, None, None), '')

    def test_soft_masked_repeats_are_not_at_primer_3_ends(self):
        rng = random.Random(23)
        seq = ''.join(rng.choice('ACGT') for _ in range(1001))
        raw = seq[:380] + seq[380:480].lower() + seq[480:520] + seq[520:620].lower() + seq[620:]
        self.fasta.write_text('>1\n' + raw + '\n')
        alt = 'A' if seq[500] != 'A' else 'C'
        Path(self.params['input_file']).write_text(f'chrom\tpos\tref\talt\n1\t501\t{seq[500]}\t{alt}\n')
        params = dict(self.params, prod_min=250, prod_max=600, tm_min=55, tm_max=65, primer_pairs=10,
                      flank=500, top_results=0, mask_repeats=True)
        log = []
        header, rows = core.run_rflp_gui_mode(params, status_cb=log.append)
        records = [dict(zip(header, row)) for row in rows if dict(zip(header, row))['primers_left']]
        self.assertTrue(records)
        for r in records:
            a = r['amplicon_start'] - 1
            left_end = a + r['primer_left_len'] - 1
            right_end = r['amplicon_end'] - r['primer_right_len']  # 3' end of the right primer (0-based)
            self.assertTrue(raw[left_end].isupper() and raw[right_end].isupper(), r['variant'])
            span = raw[a:r['amplicon_end']]
            self.assertEqual(r['repeat_pct'], round(100 * sum(c.islower() for c in span) / len(span), 1))
        self.assertFalse(any('soft-mask' in line for line in log))
        self.fasta.write_text('>1\n' + seq + '\n')
        log = []
        core.run_rflp_gui_mode(params, status_cb=log.append)
        self.assertTrue(any('без soft-mask' in line for line in log))

    def test_snp_end_distance_with_real_primer3(self):
        rng = random.Random(11)
        seq = ''.join(rng.choice('ACGT') for _ in range(1001))
        self.fasta.write_text('>1\n' + seq + '\n')
        Path(self.params['input_file']).write_text(f'chrom\tpos\tref\talt\n1\t501\t{seq[500]}\tN\n'.replace(
            'N', 'A' if seq[500] != 'A' else 'C'))
        params = dict(self.params, prod_min=250, prod_max=600, tm_min=55, tm_max=65, primer_pairs=10,
                      snp_end_distance=100, flank=500)
        report = {}
        header, rows = core.run_rflp_gui_mode(params, report=report)
        pairs = {(r['amplicon_len'], r['snp_offset_in_amplicon']) for r in (dict(zip(header, x)) for x in rows)
                 if r['primers_left']}
        self.assertTrue(pairs)
        for length, offset in pairs:
            self.assertGreaterEqual(min(offset, length - 1 - offset), 100)
        self.assertEqual(report['summary']['pairs_near_end'], 0)
        header, rows = core.run_rflp_gui_mode(dict(params, flank=80))
        row = dict(zip(header, rows[0]))
        self.assertEqual(row['status'], 'no_primers')
        self.assertIn('увеличьте фланг', row['reason'])



class SiteGroupingTests(unittest.TestCase):
    """MspI, HpaII, HapII, BsiSI and the unsold neoschizomer Sth302II share CCGG."""

    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        root = Path(temp.name)
        rng = random.Random(5)
        while True:
            seq = ''.join(rng.choice('ACGT') for _ in range(400))
            seq = seq[:197] + 'CCGG' + seq[201:]
            if seq.count('CCGG') == 1:
                break
        (root / 'toy.fa').write_text('>1\n' + seq + '\n')
        (root / 'v.tsv').write_text('chrom\tpos\tref\talt\n1\t199\tC\tT\n')
        self.params = dict(core.DEFAULT_PARAMS, fasta=str(root / 'toy.fa'), input_file=str(root / 'v.tsv'),
                           use_primer3=False, flank=500)
        self.enzymes = [e for e in core.AllEnzymes if e.site == 'CCGG']

    def run_design(self, **overrides):
        with patch.object(core, 'AllEnzymes', self.enzymes):
            header, rows = core.run_rflp_gui_mode(dict(self.params, **overrides))
        return [dict(zip(header, row)) for row in rows]

    def test_one_row_per_site_with_commercial_representative(self):
        self.assertEqual(len(self.enzymes), 5)
        self.assertEqual(len(self.run_design()), 5)
        (row,) = self.run_design(group_isoschizomers=True)
        self.assertEqual(row['enzyme'], 'MspI')
        self.assertEqual(row['isoschizomers'], 'BsiSI; HapII; HpaII; Sth302II')
        self.assertEqual(row['frags_ref'], '198+202')

    def test_commercial_filter(self):
        rows = self.run_design(commercial_only=True)
        self.assertEqual(sorted(r['enzyme'] for r in rows), ['BsiSI', 'HapII', 'HpaII', 'MspI'])
        (row,) = self.run_design(commercial_only=True, group_isoschizomers=True)
        self.assertEqual(row['isoschizomers'], 'BsiSI; HapII; HpaII')


if __name__ == '__main__':
    unittest.main()
