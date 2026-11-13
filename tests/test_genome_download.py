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


class ReviewFixTests(unittest.TestCase):
    def setUp(self):
        self._dir = tempfile.TemporaryDirectory()
        self.addCleanup(self._dir.cleanup)
        self.base = self._dir.name
        self.fasta = os.path.join(self.base, 'genoms', 'Toy_asm_v1.0.fa')
        self.names = os.path.join(self.base, 'names', 'Toy_asm_v1.0.txt')
        self.part = self.fasta + '.gz.part'

    def test_connection_closed_early_keeps_the_part(self):
        class Truncating(FakeTransport):
            def stream(self, url, start):          # the server closes the connection without an error
                self.starts.append(start)
                yield self.archive[start:start + 2500]

        with self.assertRaises(gd.GenomeDownloadError) as raised:
            gd.download_genome('GCF_000000123.4', self.base, transport=Truncating())
        self.assertIn('прервано', str(raised.exception))
        self.assertEqual(os.path.getsize(self.part), 2500)
        transport = FakeTransport()
        gd.download_genome('GCF_000000123.4', self.base, transport=transport)
        self.assertEqual(transport.starts, [2500])

    def test_genbank_assembly_has_its_own_files(self):
        class GenBank(FakeTransport):
            def text(self, url):
                return super().text(url.replace('GCA', 'GCF')).replace('GCF', 'GCA') if 'GCA' in url else super().text(url)

            def size(self, url):
                return len(self.archive)

        refseq = gd.download_genome('GCF_000000123.4', self.base, transport=FakeTransport())
        genbank = gd.download_genome('GCA_000000123.4', self.base, transport=GenBank())
        self.assertEqual(refseq, (self.fasta, self.names))
        self.assertEqual([os.path.basename(path) for path in genbank], ['Toy_asm_v1.0_GCA.fa', 'Toy_asm_v1.0_GCA.txt'])
        with open(genbank[1], encoding='utf-8') as stream:
            self.assertTrue(stream.read().startswith('1\tCM000001.1\n'))

    def test_existing_genome_gets_its_missing_names_and_index(self):
        os.makedirs(os.path.dirname(self.fasta))
        with open(self.fasta, 'wb') as stream:        # an earlier run stopped right after unpacking
            stream.write(FASTA)
        with open(self.part, 'wb') as stream:
            stream.write(b'left over')
        transport = FakeTransport()
        self.assertEqual(gd.download_genome('GCF_000000123.4', self.base, transport=transport),
                         (self.fasta, self.names))
        self.assertEqual(transport.starts, [])                         # nothing is downloaded again
        self.assertTrue(os.path.isfile(self.names))
        self.assertTrue(os.path.isfile(self.fasta + '.fai'))
        self.assertFalse(os.path.exists(self.part))

    def test_unknown_size_is_downloaded_without_the_space_check(self):
        class NoLength(FakeTransport):
            def size(self, url):
                return 0

        gd.download_genome('GCF_000000123.4', self.base, transport=NoLength(), free_space=lambda _path: 10)
        with open(self.fasta, 'rb') as stream:
            self.assertEqual(stream.read(), FASTA)

    def test_assembly_name_must_be_a_plain_file_name(self):
        class Hostile(FakeTransport):
            def text(self, url):
                return '<a href="GCF_000000123.4_..\\evil:name/">x</a>' if url == FOLDER else super().text(url)

        with self.assertRaises(gd.GenomeDownloadError):
            gd.find_assembly('GCF_000000123.4', Hostile())

    def test_missing_assembly_folder_is_not_found(self):
        class Missing(FakeTransport):
            def text(self, url):
                raise gd.GenomeDownloadError('NCBI ответил ошибкой 404', code=404)

        with self.assertRaises(gd.GenomeDownloadError) as raised:
            gd.find_assembly('GCF_000000999.1', Missing())
        self.assertEqual(str(raised.exception), 'Сборка GCF_000000999.1 не найдена на NCBI')


class HttpTransportTests(unittest.TestCase):
    """The real transport against a local server."""

    @classmethod
    def setUpClass(cls):
        import http.server
        import threading
        data = bytes(range(256)) * 40

        class Handler(http.server.BaseHTTPRequestHandler):
            def log_message(self, *_args):
                pass

            def do_HEAD(self):
                self.send_response(200)
                self.send_header('Content-Length', str(len(data)))
                self.end_headers()

            def do_GET(self):
                if self.path == '/missing':
                    self.send_error(404)
                    return
                start = int(self.headers['Range'].split('=')[1].rstrip('-')) if self.headers.get('Range') else 0
                ignore = self.path == '/norange'
                self.send_response(206 if start and not ignore else 200)
                self.send_header('Content-Length', str(len(data) - (0 if ignore else start)))
                self.end_headers()
                self.wfile.write(data if ignore else data[start:])

        cls.data = data
        cls.server = http.server.ThreadingHTTPServer(('127.0.0.1', 0), Handler)
        cls.url = f'http://127.0.0.1:{cls.server.server_address[1]}'
        threading.Thread(target=cls.server.serve_forever, daemon=True).start()

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()

    def test_size_and_ranges(self):
        transport = gd.HttpTransport()
        self.assertEqual(transport.size(self.url + '/file'), len(self.data))
        self.assertEqual(b''.join(transport.stream(self.url + '/file', 0)), self.data)
        self.assertEqual(b''.join(transport.stream(self.url + '/file', 1000)), self.data[1000:])

    def test_server_ignoring_the_range_is_an_error(self):
        with self.assertRaises(gd.GenomeDownloadError):
            b''.join(gd.HttpTransport().stream(self.url + '/norange', 1000))

    def test_http_error_carries_its_code(self):
        with self.assertRaises(gd.GenomeDownloadError) as raised:
            gd.HttpTransport().text(self.url + '/missing')
        self.assertEqual(raised.exception.code, 404)

if __name__ == '__main__':
    unittest.main()
