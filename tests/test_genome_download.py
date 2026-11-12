import gzip
import hashlib
import os
import random
import tempfile
import unittest

import genome_download as gd

FOLDER = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/000/123/'
STEM = 'GCF_000000123.4_Toy_asm_v1.0'
_RNG = random.Random(7)       # random bases, so that the archive is several kilobytes; lower case is soft-mask
FASTA = (b'>NC_000001.1 toy chromosome 1\n' + ''.join(_RNG.choice('ACGTacgt') for _ in range(40000)).encode()
         + b'\n>NC_000002.1 toy MT\nACGT\n')
ARCHIVE = gzip.compress(FASTA)
REPORT = '''# Assembly name:  Toy_asm_v1.0
# Sequence-Name	Sequence-Role	Assigned-Molecule	Assigned-Molecule-Location/Type	GenBank-Accn	Relationship	RefSeq-Accn	Assembly-Unit	Sequence-Length	UCSC-style-name
1	assembled-molecule	1	Chromosome	CM000001.1	=	NC_000001.1	Primary Assembly	40000	chr1
X	assembled-molecule	X	Chromosome	CM000002.1	=	NC_000003.1	Primary Assembly	500	na
MT	assembled-molecule	MT	Mitochondrion	CM000009.1	=	NC_000002.1	non-nuclear	4	chrM
scaffold7	unplaced-scaffold	na	na	JA000007.1	=	NW_000007.1	Primary Assembly	90	na
scaffold8	unlocalized-scaffold	1	Chromosome	JA000008.1	=	NW_000008.1	Primary Assembly	90	na
onlygb	assembled-molecule	2	Chromosome	CM000004.1	<>	na	Primary Assembly	90	na
'''


class FakeTransport:
    """The files of one assembly folder, served from memory; counts the requested ranges."""

    def __init__(self, archive=ARCHIVE, md5=None, fail_after=None):
        self.archive = archive
        self.md5 = md5 or hashlib.md5(archive).hexdigest()
        self.fail_after = fail_after          # bytes sent before the connection breaks
        self.starts = []

    def text(self, url):
        if url == FOLDER:
            return f'<a href="/genomes/all/GCF/000/000/">Parent</a> <a href="{STEM}/">{STEM}/</a>'
        if url == f'{FOLDER}{STEM}/md5checksums.txt':
            return (f'0000  ./{STEM}_cds_from_genomic.fna.gz\n{self.md5}  ./{STEM}_genomic.fna.gz\n'
                    f'1111  ./{STEM}_rna_from_genomic.fna.gz\n')
        if url == f'{FOLDER}{STEM}/{STEM}_assembly_report.txt':
            return REPORT
        raise gd.GenomeDownloadError(f'404 {url}')

    def size(self, url):
        assert url == f'{FOLDER}{STEM}/{STEM}_genomic.fna.gz', url
        return len(self.archive)

    def stream(self, url, start):
        self.starts.append(start)
        data = self.archive[start:]
        for offset in range(0, len(data), 1000):
            if self.fail_after is not None and offset >= self.fail_after:
                raise gd.GenomeDownloadError('соединение прервано')
            yield data[offset:offset + 1000]


class AddressTests(unittest.TestCase):
    def test_folder_of_an_accession(self):
        self.assertEqual(gd.accession_folder('GCF_016772045.1'),
                         'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/772/045/')
        self.assertEqual(gd.accession_folder(' gca_000003195.3 '),
                         'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/003/195/')

    def test_malformed_accessions_are_rejected(self):
        for text in ('', 'sheep', 'GCF_123', 'GCF_016772045', 'NC_056054.1', 'GCF_016772045.1/../x'):
            with self.subTest(text=text), self.assertRaises(gd.GenomeDownloadError):
                gd.accession_folder(text)

    def test_assembly_is_found_in_the_folder_listing(self):
        assembly = gd.find_assembly('GCF_000000123.4', FakeTransport())
        self.assertEqual((assembly.accession, assembly.name, assembly.stem), ('GCF_000000123.4', 'Toy_asm_v1.0', STEM))
        self.assertEqual(assembly.url, f'{FOLDER}{STEM}/')
        with self.assertRaises(gd.GenomeDownloadError):     # another version of the same assembly
            gd.find_assembly('GCF_000000123.5', FakeTransport())

    def test_presets_are_well_formed_and_unique(self):
        accessions = [preset[2] for preset in gd.PRESETS]
        self.assertEqual(len(accessions), len(set(accessions)))
        self.assertGreaterEqual(len(accessions), 10)        # nine species, two sheep assemblies
        for accession in accessions:
            gd.accession_folder(accession)


class NamesTests(unittest.TestCase):
    def test_chromosomes_of_a_refseq_assembly(self):
        self.assertEqual(gd.chromosome_names(REPORT, 'GCF_000000123.4'),
                         [('1', 'NC_000001.1'), ('X', 'NC_000003.1'), ('MT', 'NC_000002.1')])

    def test_chromosomes_of_a_genbank_assembly(self):
        self.assertEqual(gd.chromosome_names(REPORT, 'GCA_000000123.4'),
                         [('1', 'CM000001.1'), ('X', 'CM000002.1'), ('MT', 'CM000009.1'), ('2', 'CM000004.1')])


class DownloadTests(unittest.TestCase):
    def setUp(self):
        self._dir = tempfile.TemporaryDirectory()
        self.addCleanup(self._dir.cleanup)
        self.base = self._dir.name
        self.fasta = os.path.join(self.base, 'genoms', 'Toy_asm_v1.0.fa')
        self.names = os.path.join(self.base, 'names', 'Toy_asm_v1.0.txt')
        self.part = self.fasta + '.gz.part'

    def test_download_unpacks_indexes_and_writes_the_names(self):
        progress, stages = [], []
        result = gd.download_genome('GCF_000000123.4', self.base, transport=FakeTransport(),
                                    progress_cb=lambda done, total: progress.append((done, total)),
                                    status_cb=stages.append)
        self.assertEqual(result, (self.fasta, self.names))
        with open(self.fasta, 'rb') as stream:
            self.assertEqual(stream.read(), FASTA)                    # soft-mask (lower case) is kept
        self.assertTrue(os.path.isfile(self.fasta + '.fai'))
        with open(self.names, encoding='utf-8') as stream:
            self.assertEqual(stream.read(), '1\tNC_000001.1\nX\tNC_000003.1\nMT\tNC_000002.1\n')
        self.assertEqual(sorted(os.listdir(os.path.dirname(self.fasta))), ['Toy_asm_v1.0.fa', 'Toy_asm_v1.0.fa.fai'])
        self.assertEqual(progress[-1], (len(ARCHIVE), len(ARCHIVE)))
        self.assertEqual(len(stages), 4)                              # download, checksum, unpack, index

    def test_names_file_works_with_the_core(self):
        import rflp_core
        gd.download_genome('GCF_000000123.4', self.base, transport=FakeTransport())
        self.assertEqual(rflp_core.load_names2(self.names)['MT'], 'NC_000002.1')

    def test_broken_download_is_resumed(self):
        with self.assertRaises(gd.GenomeDownloadError):
            gd.download_genome('GCF_000000123.4', self.base, transport=FakeTransport(fail_after=3000))
        self.assertEqual(os.path.getsize(self.part), 3000)
        self.assertFalse(os.path.exists(self.fasta))
        transport = FakeTransport()
        gd.download_genome('GCF_000000123.4', self.base, transport=transport)
        self.assertEqual(transport.starts, [3000])
        with open(self.fasta, 'rb') as stream:
            self.assertEqual(stream.read(), FASTA)

    def test_cancel_keeps_the_part_for_a_later_resume(self):
        calls = []

        def cancel():
            calls.append(1)
            return len(calls) > 3

        with self.assertRaises(gd.DownloadCancelled):
            gd.download_genome('GCF_000000123.4', self.base, transport=FakeTransport(), cancel_cb=cancel)
        self.assertTrue(0 < os.path.getsize(self.part) < len(ARCHIVE))
        self.assertFalse(os.path.exists(self.fasta))

    def test_wrong_checksum_discards_the_download(self):
        with self.assertRaises(gd.GenomeDownloadError) as raised:
            gd.download_genome('GCF_000000123.4', self.base, transport=FakeTransport(md5='0' * 32))
        self.assertIn('Контрольная сумма', str(raised.exception))
        self.assertEqual(os.listdir(os.path.dirname(self.fasta)), [])

    def test_existing_names_file_is_not_overwritten(self):
        os.makedirs(os.path.dirname(self.names))
        with open(self.names, 'w', encoding='utf-8') as stream:
            stream.write('27\tNC_000003.1\n')
        gd.download_genome('GCF_000000123.4', self.base, transport=FakeTransport())
        with open(self.names, encoding='utf-8') as stream:
            self.assertEqual(stream.read(), '27\tNC_000003.1\n')

    def test_existing_genome_is_not_downloaded_again(self):
        gd.download_genome('GCF_000000123.4', self.base, transport=FakeTransport())
        transport = FakeTransport()
        self.assertEqual(gd.download_genome('GCF_000000123.4', self.base, transport=transport),
                         (self.fasta, self.names))
        self.assertEqual(transport.starts, [])

    def test_not_enough_disk_space_is_reported_before_the_download(self):
        transport = FakeTransport()
        with self.assertRaises(gd.GenomeDownloadError) as raised:
            gd.download_genome('GCF_000000123.4', self.base, transport=transport, free_space=lambda _path: 10)
        self.assertIn('Недостаточно места', str(raised.exception))
        self.assertEqual(transport.starts, [])


if __name__ == '__main__':
    unittest.main()
