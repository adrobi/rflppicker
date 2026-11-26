import gzip
import tempfile
import unittest
from pathlib import Path

import rflp_core as core


class VcfInputTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        self.fasta = self.root / 'genome.fa'
        self.fasta.write_text('>chr1\n' + 'A' * 300 + '\n')

    def run_vcf(self, text, **options):
        path = self.root / 'targets.vcf.gz'
        with gzip.open(path, 'wt') as stream:
            stream.write('##fileformat=VCFv4.3\n' + text)
        params = dict(core.DEFAULT_PARAMS, fasta=str(self.fasta), input_file=str(path), use_primer3=False,
                      enzyme_mode='only', lab_enzymes=['EcoRI'], top_results=1)
        params.update(options)
        report = {}
        header, rows = core.run_analysis(params, report=report)
        return [dict(zip(header, row)) for row in rows], report

    def test_without_samples_keeps_all_alts_and_filters_by_region(self):
        rows, report = self.run_vcf('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                                    'chr1\t150\trs1\tA\tC,G\t.\tPASS\t.\n'
                                    'chr1\t151\trs2\tA\tT\t.\tPASS\t.\n', vcf_region='chr1:150-150')
        self.assertEqual([r['variant'] for r in rows], ['chr1:150 A>C,G'])
        self.assertEqual(report['summary']['accepted'], 1)

    def test_selected_samples_keep_only_observed_alts_in_record_order(self):
        rows, report = self.run_vcf('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tx\ty\tz\n'
                                    'chr1\t150\trs1\tA\tC,G,T\t.\tPASS\t.\tDP:GT\t7:3|.\t4:./.\t5:0/2\n',
                                    vcf_samples=['z', 'x', 'y'])
        self.assertEqual([r['variant'] for r in rows], ['chr1:150 A>G,T'])
        self.assertEqual(report['summary']['accepted'], 1)

    def test_missing_symbolic_and_reference_mismatch_have_separate_reasons(self):
        rows, report = self.run_vcf('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tx\n'
                                    'chr1\t150\t.\tA\tC\t.\tPASS\t.\tGT\t./.\n'
                                    'chr1\t151\t.\tA\t<DEL>\t.\tPASS\t.\tGT\t0/1\n'
                                    'chr1\t152\t.\tC\tA\t.\tPASS\t.\tGT\t0/1\n', vcf_samples=['x'])
        self.assertEqual(rows, [])
        self.assertEqual([r['reason'] for r in report['rejected']],
                         ['vcf_no_selected_alt', 'vcf_symbolic', 'ref_not_match'])

    def test_input_vcf_and_background_vcf_are_independent(self):
        background = self.root / 'background.vcf'
        background.write_text('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                              'chr1\t155\t.\tA\tG\t.\tPASS\t.\n')
        rows, report = self.run_vcf('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                                    'chr1\t150\trs1\tA\tC,G\t.\tPASS\t.\n', background_vcf=str(background))
        self.assertEqual([r['variant'] for r in rows], ['chr1:150 A>C,G'])
        self.assertEqual(report['summary']['background_variants_loaded'], 1)

    def test_exact_vcf_contig_wins_when_prefixed_and_plain_names_both_exist(self):
        self.fasta.write_text('>chr1\n' + 'A' * 300 + '\n>1\n' + 'G' * 300 + '\n')
        rows, report = self.run_vcf('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                                    'chr1\t150\trs1\tA\tG\t.\tPASS\t.\n')
        self.assertEqual(rows[0]['mapped_id'], 'chr1')
        self.assertEqual(rows[0]['variant'], 'chr1:150 A>G')
        self.assertFalse(report.get('rejected'))
