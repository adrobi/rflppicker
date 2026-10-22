import tempfile
import unittest

import snp_summary
from snp_summary import InputRecord, SnpSummary
from tests.helpers import HEADER, core_row, write_text


def rec(line, chrom, pos, ref, alt, snp_id=''):
    return InputRecord(line=line, snp_id=snp_id, chrom=chrom, pos=pos, ref=ref, alt=alt, raw='')


class InputRecordTests(unittest.TestCase):
    def setUp(self):
        self._dir = tempfile.TemporaryDirectory()
        self.tmp = self._dir.name

    def tearDown(self):
        self._dir.cleanup()

    def test_records_with_ids_and_parse_errors(self):
        path = write_text(self.tmp, 'v.tsv',
                          'chrom\tpos\tref\talt\tid\nchr1\t10\tA\tG\tDGAT1\n1\tx\tA\tG\tbad\n2\t5\tC\tT\t\n')
        records, has_ids = snp_summary.read_input_records(path)
        self.assertTrue(has_ids)
        self.assertEqual([(r.line, r.snp_id, r.chrom, r.pos) for r in records],
                         [(2, 'DGAT1', '1', 10), (3, 'bad', '1', None), (4, '', '2', 5)])
        self.assertEqual(records[1].raw, '1\tx\tA\tG\tbad')

    def test_records_without_id_column(self):
        path = write_text(self.tmp, 'v.csv', 'chrom,pos,ref,alt\n1,10,A,G\n')
        records, has_ids = snp_summary.read_input_records(path)
        self.assertFalse(has_ids)
        self.assertEqual(records[0].snp_id, '')

    def test_short_line_is_kept_as_a_record(self):
        path = write_text(self.tmp, 'v.csv', 'chrom,pos,ref,alt\n1,10\n')
        records, _ = snp_summary.read_input_records(path)
        self.assertEqual([(r.line, r.pos) for r in records], [(2, None)])


class JoinTests(unittest.TestCase):
    def test_join_by_line_values_variant_and_swap(self):
        records = [rec(2, '1', 10, 'A', 'G'), rec(3, '1', 11, 'C', 'T'), rec(4, '1', 12, 'G', 'A'),
                   rec(5, '1', 13, 'T', 'C'), rec(6, '1', 10, 'A', 'G')]
        rows = [core_row(variant='1:10 A>G', enzyme='BamHI', assay_rank=1),
                core_row(variant='1:10 A>G', enzyme='AciI', assay_rank=2),
                core_row(variant='1:12 A>G', enzyme='', status='no_primers',
                         reason='Primer3 не вернул подходящих пар')]
        rejected = [{'line': 3, 'reason': 'parse_error', 'details': 'x'},
                    {'chrom': '1', 'pos': 13, 'ref': 'T', 'alt': 'C', 'reason': 'bad_chrom', 'details': 'y'}]
        s = snp_summary.build_summaries(records, HEADER, rows, rejected)
        self.assertEqual([x['enzyme'] for x in s[0].schemes], ['BamHI', 'AciI'])
        self.assertEqual(s[0].schemes[1]['_index'], 1)
        self.assertEqual(s[1].rejection['reason'], 'parse_error')
        self.assertEqual(s[2].failure['status'], 'no_primers')      # REF/ALT swapped by normalization
        self.assertEqual(s[3].rejection['reason'], 'bad_chrom')
        self.assertEqual([x['enzyme'] for x in s[4].schemes], ['BamHI', 'AciI'])   # repeated line
        self.assertEqual(s[0].label, '1:10 A>G')
        self.assertIsInstance(s[0], SnpSummary)

    def test_skip_norm_rejections_match_by_values_not_by_variant_index(self):
        # With skip_norm the core numbers rejected variants 1, 2, ... — not file lines.
        records = [rec(2, '1', 10, 'A', 'G'), rec(3, '1', 13, 'T', 'C')]
        rows = [core_row(variant='1:10 A>G', enzyme='BamHI', assay_rank=1)]
        rejected = [{'line': 2, 'chrom': '1', 'pos': 13, 'ref': 'T', 'alt': 'C', 'reason': 'ref_not_match'}]
        s = snp_summary.build_summaries(records, HEADER, rows, rejected)
        self.assertIsNone(s[0].rejection)
        self.assertEqual([x['enzyme'] for x in s[0].schemes], ['BamHI'])
        self.assertEqual(s[1].rejection['reason'], 'ref_not_match')


if __name__ == '__main__':
    unittest.main()
