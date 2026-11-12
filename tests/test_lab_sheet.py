import os
import random
import re
import tempfile
import unittest
from unittest.mock import patch

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
from PyQt6.QtWidgets import QApplication

import i18n
import lab_sheet
import rflp_core
import snp_summary
from tests.helpers import CHECKS_OFF, write_text


def page_count(path: str) -> int:
    with open(path, 'rb') as stream:
        return len(re.findall(rb'/Type\s*/Page\b(?!s)', stream.read()))


class LabSheetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])
        cls._dir = tempfile.TemporaryDirectory()
        rng = random.Random(19)                       # the natural EcoRI scenario of JointDesignTests
        sequence = ''.join(rng.choice('ACGT') for _ in range(400))
        sequence = sequence[:197] + 'GAATTC' + sequence[203:]
        cls.fasta = write_text(cls._dir.name, 'toy.fa', '>1\n' + sequence + '\n')
        variants = write_text(cls._dir.name, 'v.tsv', 'chrom\tpos\tref\talt\tid\n1\t200\tA\tG\tDGAT1\n'
                                                      '1\t100\tN\tG\tbad\n1\t201\tA\tG\tnone\n')
        cls.params = dict(rflp_core.DEFAULT_PARAMS, **CHECKS_OFF, fasta=cls.fasta, input_file=variants, flank=500,
                          prod_min=200, prod_max=500, tm_min=50, tm_max=70, primer_pairs=2, max_cuts=1,
                          gain_loss_only=True)
        report = {}
        header, rows = rflp_core.run_rflp_gui_mode(cls.params, report=report)
        records, _ = snp_summary.read_input_records(variants)
        cls.summaries = snp_summary.build_summaries(records, header, rows, report.get('rejected', []))
        snp_summary.assign_primer_bases(cls.summaries)

    @classmethod
    def tearDownClass(cls):
        cls._dir.cleanup()

    def setUp(self):
        self.addCleanup(i18n.set_language, 'ru')
        self.out = os.path.join(self._dir.name, f'{self.id().rsplit(".", 1)[-1]}.pdf')

    def test_fixture_has_one_snp_with_an_assay(self):
        self.assertEqual([s.scheme is not None for s in self.summaries], [True, False, False])

    def test_fields_of_a_sheet(self):
        fields = dict(lab_sheet.sheet_fields(self.summaries[0], 'off'))
        scheme = self.summaries[0].scheme
        self.assertEqual(fields['Фермент'], scheme['enzyme'])
        self.assertEqual(fields['Сайт'], scheme['site'])
        self.assertEqual(fields['Продукт, п.н.'], str(scheme['product_size']))
        self.assertEqual(fields['Полосы REF/REF'], scheme['bands_ref_ref'])
        self.assertEqual(fields['Полосы ALT/ALT'], scheme['bands_alt_alt'])
        self.assertEqual(fields['Статус'], 'Надёжно')
        self.assertEqual(fields['Тип схемы'], 'природный сайт')
        self.assertNotIn('Примечание', fields)             # nothing to warn about

    def test_primer_rows(self):
        rows = lab_sheet.primer_rows(self.summaries[0])
        scheme = self.summaries[0].scheme
        self.assertEqual([row[0] for row in rows], ['DGAT1_F', 'DGAT1_R'])
        self.assertEqual(rows[0][1:3], (scheme['primers_left'], str(len(scheme['primers_left']))))
        self.assertRegex(rows[1][3], r'^\d\d,\d$')            # Tm with the decimal comma of the language

    def test_english_fields(self):
        i18n.set_language('en')
        fields = dict(lab_sheet.sheet_fields(self.summaries[0], 'off'))
        self.assertEqual(fields['Status'], 'Robust')
        self.assertIsNone(re.search('[А-Яа-яЁё]', ' '.join(f'{k} {v}' for k, v in fields.items())))
        self.assertIsNone(re.search('[А-Яа-яЁё]', ' '.join(lab_sheet.handwritten_fields())))

    def test_one_page_per_snp_with_an_assay(self):
        count = lab_sheet.export_lab_sheets(self.out, self.summaries, self.fasta, self.params, 'off')
        self.assertEqual(count, 1)
        with open(self.out, 'rb') as stream:
            self.assertEqual(stream.read(5), b'%PDF-')
        self.assertEqual(page_count(self.out), 1)
        self.assertGreater(os.path.getsize(self.out), 5000)

    def test_several_pages(self):
        count = lab_sheet.export_lab_sheets(self.out, [self.summaries[0]] * 3, self.fasta, self.params, 'off')
        self.assertEqual((count, page_count(self.out)), (3, 3))

    def test_nothing_to_print_writes_no_file(self):
        self.assertEqual(lab_sheet.export_lab_sheets(self.out, self.summaries[1:], self.fasta, self.params, 'off'), 0)
        self.assertFalse(os.path.exists(self.out))

    def test_missing_genome_still_gives_the_sheet(self):
        count = lab_sheet.export_lab_sheets(self.out, self.summaries, os.path.join(self._dir.name, 'gone.fa'),
                                            self.params, 'off')
        self.assertEqual((count, page_count(self.out)), (1, 1))

    def test_failed_export_keeps_the_previous_file(self):
        with open(self.out, 'wb') as stream:
            stream.write(b'previous')
        with patch.object(lab_sheet, 'draw_sheet', side_effect=RuntimeError('boom')), \
                self.assertRaises(RuntimeError):
            lab_sheet.export_lab_sheets(self.out, self.summaries, self.fasta, self.params, 'off')
        with open(self.out, 'rb') as stream:
            self.assertEqual(stream.read(), b'previous')
        self.assertEqual([name for name in os.listdir(self._dir.name) if name.startswith('.')], [])   # no temp file


if __name__ == '__main__':
    unittest.main()
