import contextlib
import csv
import io
import random
import tempfile
import unittest
from pathlib import Path

from openpyxl import load_workbook

import rflp_cli
import rflp_core as core


class CliTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        rng = random.Random(3)
        sequence = ''.join(rng.choice('ACGT') for _ in range(400))
        sequence = sequence[:197] + 'GAATTC' + sequence[203:]
        self.fasta = self.root / 'toy.fa'
        self.fasta.write_text('>1\n' + sequence + '\n')
        self.variants = self.root / 'variants.tsv'
        self.variants.write_text('chrom\tpos\tref\talt\n1\t200\tA\tG\n1\t2\tA\tA\n')

    def run_cli(self, *args):
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            code = rflp_cli.main(['--fasta', str(self.fasta), '--variants', str(self.variants), *args])
        return code, stderr.getvalue()

    def test_csv_matches_core_run(self):
        out = self.root / 'out.csv'
        rejected = self.root / 'rejected.csv'
        code, log = self.run_cli('--out', str(out), '--rejected', str(rejected), '--no-primer3',
                                 '--flank', '150', '--min-frag', '20', '--preset', 'v1.1')
        self.assertEqual(code, 0, log)
        params = dict(core.DEFAULT_PARAMS, fasta=str(self.fasta.resolve()),
                      input_file=str(self.variants.resolve()), use_primer3=False, flank=150, min_frag=20)
        header, rows = core.run_rflp_gui_mode(params)
        with out.open(encoding='utf-8-sig', newline='') as stream:
            written = list(csv.reader(stream))
        self.assertEqual(written[0], header)
        self.assertEqual(written[1:], [['' if v is None else str(v) for v in row] for row in rows])
        self.assertTrue(any(r['enzyme'] == 'EcoRI' for r in (dict(zip(header, row)) for row in rows)))
        with rejected.open(encoding='utf-8-sig', newline='') as stream:
            self.assertEqual(len(list(csv.DictReader(stream))), 1)
        self.assertIn('Сохранено', log)

    def test_excel_has_run_and_rejected_sheets(self):
        out = self.root / 'out.xlsx'
        code, log = self.run_cli('--out', str(out), '--no-primer3', '--min-frag', '20', '-q')
        self.assertEqual((code, log), (0, ''))
        workbook = load_workbook(out, read_only=True)
        self.assertEqual(workbook.sheetnames, ['Results', 'Run', 'Rejected'])
        run = {row[0]: row[1] for row in workbook['Run'].iter_rows(min_row=2, values_only=True)}
        self.assertEqual(run['app_version'], core.APP_VERSION)
        workbook.close()

    def test_invalid_parameters_and_extension(self):
        code, log = self.run_cli('--out', str(self.root / 'x.csv'), '--prod-min', '700', '--prod-max', '600')
        self.assertEqual(code, 2)
        self.assertIn('размер продукта', log)
        with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
            self.run_cli('--out', str(self.root / 'x.txt'))


if __name__ == '__main__':
    unittest.main()
