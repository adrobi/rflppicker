import gzip
import os
import random
import struct
import tempfile
import unittest
from pathlib import Path

from Bio.bgzf import BgzfWriter
from pyfaidx import Fasta

import rflp_core as core
from tabix_index import index_path


class IndexedDataTests(unittest.TestCase):
    def setUp(self):
        work = tempfile.TemporaryDirectory()
        self.addCleanup(work.cleanup)
        self.path = Path(work.name)

    def test_bgzf_fasta_with_indices_matches_plain_fasta_and_specificity(self):
        rng = random.Random(784)
        seq = ''.join(rng.choices('ACGT', k=5000))
        plain, compressed = self.path / 'genome.fa', self.path / 'genome.fa.gz'
        variants = self.path / 'loci.tsv'
        fasta = '>chr1\n' + '\n'.join(seq[i:i + 60] for i in range(0, len(seq), 60)) + '\n'
        plain.write_text(fasta)
        with BgzfWriter(compressed, 'wb') as writer:
            writer.write(fasta.encode())
        Fasta(compressed).close()
        self.assertTrue(Path(str(compressed) + '.fai').is_file())
        self.assertTrue(Path(str(compressed) + '.gzi').is_file())
        variants.write_text(f'chrom\tpos\tref\talt\n1\t2500\t{seq[2499]}\t{seq[2499]}{"A" * 50}\n')
        params = dict(core.DEFAULT_PARAMS, input_file=str(variants), top_results=1, specificity='exclude')
        reports = [{}, {}]
        first = core.run_analysis(dict(params, fasta=str(plain)), report=reports[0])
        second = core.run_analysis(dict(params, fasta=str(compressed)), report=reports[1])
        self.assertEqual(first, second)
        self.assertEqual(reports[0]['summary'], reports[1]['summary'])
        self.assertEqual(reports[1]['summary']['specificity_scans'], 1)

    def test_tabix_queries_merge_overlaps_and_keep_inclusive_boundaries(self):
        plain, compressed = self.path / 'background.vcf', self.path / 'background.vcf.gz'
        head = b'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        positions = [9, 10, 15, 20, 21, 40000]
        records = [f'1\t{pos}\t.\tA\tC\t.\tPASS\t.\n'.encode() for pos in positions]
        plain.write_bytes(head + b''.join(records))
        with BgzfWriter(compressed, 'wb') as writer:
            writer.write(head)
            first = writer.tell()
            writer.write(b''.join(records[:-1]))
            writer.flush()
            second = writer.tell()
            writer.write(records[-1])
            end = writer.tell()
        with gzip.open(str(compressed) + '.tbi', 'wb') as index:
            index.write(b'TBI\x01' + struct.pack('<8I', 1, 2, 1, 2, 0, ord('#'), 0, 2) + b'1\x00')
            index.write(struct.pack('<IIIQQIIQQIQQQ', 2, 4681, 1, first, second, 4683, 1, second, end,
                                    3, first, 0, second))

        regions = {'NC_000001.1': [(10, 15), (14, 20)], 'chr2': [(1, 10)]}
        names = {'1': 'NC_000001.1'}
        ordinary, ordinary_read = core.load_background_variants(plain, regions, names)
        indexed, indexed_read = core.load_background_variants(compressed, regions, names)
        self.assertEqual(indexed, ordinary)
        self.assertEqual([r[0] for r in indexed['NC_000001.1']], [10, 15, 20])
        self.assertLess(indexed_read, ordinary_read)

        tbi = Path(str(compressed) + '.tbi')
        csi = Path(str(compressed) + '.csi')
        csi.write_bytes(b'old index')
        os.utime(csi, (1, 1))
        self.assertEqual(index_path(compressed), str(tbi))
        os.utime(tbi, (1, 1))
        with BgzfWriter(compressed, 'ab') as writer:
            writer.write(b'1\t50000\t.\tA\tT\t.\tPASS\t.\n')
        found, _ = core.load_background_variants(compressed, {'NC_000001.1': [(49999, 50001)]}, names)
        self.assertEqual(found, {'NC_000001.1': [(50000, 'A', ('T',))]})

    def test_bg_vcf_uses_exact_fasta_name_before_alias(self):
        vcf = self.path / 'background.vcf'
        vcf.write_text('chr1\t15\t.\tA\tT\t.\tPASS\t.\n')
        regions = {'NC_000001.1': [(10, 20)]}
        names = {'1': 'NC_000001.1'}
        found, _ = core.load_background_variants(vcf, regions, names)
        self.assertEqual(found, {'NC_000001.1': [(15, 'A', ('T',))]})
        found, _ = core.load_background_variants(vcf, regions, names, fasta_contigs={'chr1', 'NC_000001.1'})
        self.assertEqual(found, {'NC_000001.1': []})
