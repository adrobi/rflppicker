import random
import tempfile
import unittest

import i18n
import rflp_core
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

    def test_swapped_alleles_are_flagged(self):
        rows = [core_row(variant='1:12 A>G', enzyme='AluI', genotype_quality='robust', assay_rank=1)]
        swapped = snp_summary.build_summaries([rec(2, '1', 12, 'G', 'A')], HEADER, rows, [])[0]
        self.assertTrue(swapped.swapped)
        self.assertEqual(swapped.note('off'), 'REF и ALT переставлены по геному: REF = A')
        straight = snp_summary.build_summaries([rec(2, '1', 12, 'A', 'G')], HEADER, rows, [])[0]
        self.assertFalse(straight.swapped)
        self.assertEqual(straight.note('off'), '')
        failure = [core_row(variant='1:12 A>G', enzyme='', status='no_primers', reason='Primer3 не вернул подходящих пар')]
        self.assertTrue(snp_summary.build_summaries([rec(2, '1', 12, 'G', 'A')], HEADER, failure, [])[0]
                        .note('off').startswith('REF и ALT переставлены по геному: REF = A; Primer3'))

    def test_skip_norm_rejections_match_by_values_not_by_variant_index(self):
        # With skip_norm the core numbers rejected variants 1, 2, ... — not file lines.
        records = [rec(2, '1', 10, 'A', 'G'), rec(3, '1', 13, 'T', 'C')]
        rows = [core_row(variant='1:10 A>G', enzyme='BamHI', assay_rank=1)]
        rejected = [{'line': 2, 'chrom': '1', 'pos': 13, 'ref': 'T', 'alt': 'C', 'reason': 'ref_not_match'}]
        s = snp_summary.build_summaries(records, HEADER, rows, rejected)
        self.assertIsNone(s[0].rejection)
        self.assertEqual([x['enzyme'] for x in s[0].schemes], ['BamHI'])
        self.assertEqual(s[1].rejection['reason'], 'ref_not_match')


class RealCoreJoinTests(unittest.TestCase):
    def test_repeated_and_swapped_lines_do_not_multiply_schemes(self):
        rng = random.Random(3)
        sequence = ''.join(rng.choice('ACGT') for _ in range(400))
        sequence = sequence[:197] + 'GAATTC' + sequence[203:]
        ref = sequence[199]
        alt = next(base for base in 'ACGT' if base != ref)
        with tempfile.TemporaryDirectory() as directory:
            fasta = write_text(directory, 'toy.fa', '>1\n' + sequence + '\n')
            variants = write_text(directory, 'v.tsv', f'chrom\tpos\tref\talt\n1\t200\t{ref}\t{alt}\n'
                                                      f'1\t200\t{ref}\t{alt}\n1\t200\t{alt}\t{ref}\n')
            params = dict(rflp_core.DEFAULT_PARAMS, fasta=fasta, input_file=variants, use_primer3=False,
                          flank=150, min_frag=20)
            report = {}
            header, rows = rflp_core.run_rflp_gui_mode(params, report=report)
            records, _ = snp_summary.read_input_records(variants)
            summaries = snp_summary.build_summaries(records, header, rows, report.get('rejected', []))
        self.assertEqual(len(rows) % 3, 0)
        single = len(rows) // 3
        self.assertGreater(single, 0)
        for s in summaries:
            self.assertEqual(len(s.schemes), single)
            ranks = [x['assay_rank'] for x in s.schemes]
            self.assertEqual(len(ranks), len(set(ranks)))
        self.assertEqual([s.swapped for s in summaries], [False, False, True])


def summary(schemes=(), failure=None, rejection=None, snp_id=''):
    return SnpSummary(record=rec(2, '1', 10, 'A', 'G', snp_id), schemes=[dict(s) for s in schemes],
                      failure=failure, rejection=rejection)


class StatusTests(unittest.TestCase):
    def test_status_follows_chosen_scheme(self):
        s = summary(schemes=[dict(enzyme='Bse1I', genotype_quality='robust'),
                             dict(enzyme='HinfI', genotype_quality='nominal_only',
                                  weakest_genotypes='REF/REF–REF/ALT')])
        self.assertEqual(s.status, snp_summary.ROBUST)
        self.assertFalse(s.manual)
        self.assertEqual(s.scheme['enzyme'], 'Bse1I')
        s.choose(1)
        self.assertTrue(s.manual)
        self.assertEqual(s.status, snp_summary.NOMINAL)
        self.assertEqual(s.note('fallback'), 'Выбрано вручную, лучшая — Bse1I; '
                                             'Только при исходном разрешении; труднее всего: REF/REF–REF/ALT')
        with self.assertRaises(IndexError):
            s.choose(2)
        with self.assertRaises(IndexError):
            s.choose(-1)

    def test_statuses(self):
        self.assertEqual(summary(failure=dict(status='no_primers', reason='x')).status, snp_summary.NONE)
        self.assertEqual(summary(failure=dict(status='primer3_error', reason='x')).status, snp_summary.NONE)
        self.assertEqual(summary(failure=dict(status='no_enzyme_found', reason='x')).status, snp_summary.NONE)
        self.assertEqual(summary(rejection=dict(reason='bad_chrom')).status, snp_summary.ERROR)
        self.assertEqual(summary(failure=dict(status='analysis_error', reason='REF не совпадает с FASTA')).status,
                         snp_summary.ERROR)
        self.assertIsNone(summary().scheme)
        self.assertEqual(snp_summary.STATUS_TITLES[snp_summary.NOMINAL], 'Только на хорошем геле')
        self.assertEqual(snp_summary.STATUS_ORDER, ('robust', 'nominal', 'none', 'error'))

    def test_line_without_any_result_is_an_error(self):            # Review Focus 1
        s = summary()
        self.assertEqual((s.status, s.note('off')), (snp_summary.ERROR, 'Нет результата для этой строки'))

    def test_ambiguous_scheme_has_no_assay_status(self):            # Review Focus 2
        s = summary(schemes=[dict(enzyme='AluI', genotype_quality='ambiguous')])
        self.assertEqual(s.status, snp_summary.NONE)
        self.assertTrue(s.note('off').startswith('Три генотипа неразличимы по присутствию полос'))


class NoteTests(unittest.TestCase):
    def test_note_order_for_scheme_warnings(self):
        s = summary(schemes=[dict(enzyme='HinfI', genotype_quality='robust', design='dCAPS',
                                  primer_mismatches="right 3'-2 G>C", primer_site='HinfI', background_variants=2,
                                  background_site_variants='37525140 A>G', repeat_pct=52.5)])
        self.assertEqual(s.note('always'), "dCAPS: right 3'-2 G>C; Сайт фермента в праймере; "
                         "Фоновые варианты меняют разрезы: 37525140 A>G; Повторы 52.5% ампликона")

    def test_other_scheme_notes(self):
        robust = dict(enzyme='BamHI', genotype_quality='robust', design='natural')
        self.assertEqual(summary(schemes=[robust]).note('fallback'), '')
        self.assertEqual(summary(schemes=[dict(robust, background_variants=3, background_site_variants='')])
                         .note('fallback'), 'Известных вариантов в ампликоне: 3')
        self.assertEqual(summary(schemes=[dict(robust, repeat_pct=49.9)]).note('fallback'), '')
        self.assertEqual(summary(schemes=[dict(robust, status='window_only')]).note('fallback'),
                         'Поиск по окну, праймеры не подобраны')

    def test_failure_hints(self):
        reason = 'Не найдено сочетание для различения трёх генотипов'
        s = summary(failure=dict(status='no_discriminating_assay', reason=reason))
        self.assertEqual(s.status, snp_summary.NONE)
        self.assertEqual(s.note('fallback'), f'{reason}. Попробуйте dCAPS «Для каждого SNP»')
        self.assertEqual(s.note('always'), f'{reason}. Попробуйте другой диапазон продукта')
        self.assertEqual(summary(failure=dict(status='no_enzyme_found', reason='Подходящий фермент не найден'))
                         .note('off'), 'Подходящий фермент не найден. Попробуйте dCAPS «Для каждого SNP»')
        self.assertEqual(summary(failure=dict(status='no_primers', reason='Primer3 не вернул подходящих пар'))
                         .note('off'), 'Primer3 не вернул подходящих пар. '
                                       'Попробуйте расширить диапазон продукта или Tm')
        self.assertEqual(summary(failure=dict(status='primer3_error', reason='bad input')).note('off'),
                         'bad input. Проверьте параметры Primer3 в окне «Дополнительно»')

    def test_errors(self):
        def note(**rejection):
            return summary(rejection=rejection).note('off')
        self.assertEqual(note(reason='bad_chrom', details='x'), 'Хромосома не найдена в геноме')
        self.assertEqual(note(reason='parse_error', details='Ожидалось минимум 4 поля'),
                         'Строка не прочитана: Ожидалось минимум 4 поля')
        self.assertEqual(note(reason='pos_oob', details='x'), 'Позиция за пределами хромосомы')
        self.assertEqual(note(reason='bad_ref_chars', details='x'), 'Не одиночная замена A/C/G/T')
        self.assertEqual(note(reason='bad_alt_chars', details='x'), 'Не одиночная замена A/C/G/T')
        self.assertEqual(note(reason='minus_strand_like', details='x'), 'REF похож на обратную цепь')
        self.assertEqual(note(reason='multi_alt_contains_ref', details='x'), 'ALT совпадает с геномом')
        self.assertEqual(note(reason='ref_not_match', details='x'), 'REF не совпадает с геномом')
        self.assertEqual(note(reason='unsupported_variant', details='Поддерживаются только одиночные SNP A/C/G/T'),
                         'Поддерживаются только одиночные SNP A/C/G/T')
        self.assertEqual(note(reason='something_new'), 'something_new')
        self.assertEqual(summary(failure=dict(status='analysis_error', reason='REF не совпадает с FASTA'))
                         .note('off'), 'REF не совпадает с FASTA')


HOX_SCHEME = dict(variant='11:37525128 C>G', mapped_id='NC_056064.1', enzyme='Bse1I', site='ACTGG',
                  isoschizomers='BsrI', suppliers='Thermo', design='natural', genotype_quality='robust',
                  primers_left='TATGAGCTTCTCTCCGCCAG', primers_right='CGCCTCCAAAGTAGCCGTAA',
                  tm_left=59.04, tm_right=60.1, product_size=475, bands_ref_ref='448',
                  bands_ref_alt='212; 236; 448', bands_alt_alt='212; 236', genotype_margin=15.77,
                  worst_margin=10.52, weakest_genotypes='REF/REF–REF/ALT', assay_rank=1)


class ExportTests(unittest.TestCase):
    def setUp(self):
        records = [rec(2, '11', 37525128, 'C', 'G', 'HOXB13'), rec(3, '1', 91308651, 'C', 'T', 'OLFML3')]
        self.core_rows = [core_row(**HOX_SCHEME),
                          core_row(**dict(HOX_SCHEME, enzyme='BsrDI', assay_rank=2, isoschizomers='')),
                          core_row(variant='1:91308651 C>T', mapped_id='NC_056054.1', enzyme='',
                                   status='no_discriminating_assay',
                                   reason='Не найдено сочетание для различения трёх генотипов')]
        self.summaries = snp_summary.build_summaries(records, HEADER, self.core_rows, [])

    def test_primer_bases(self):                       # Review Focus 3
        s = [summary(snp_id='DGAT1'), summary(snp_id='DGAT1'), summary(snp_id='ДГАТ1'),
             summary(snp_id='rs 12/3'), summary(snp_id='')]
        snp_summary.assign_primer_bases(s)
        self.assertEqual([x.primer_base for x in s], ['DGAT1', 'DGAT1_2', '1_10', 'rs_12_3', '1_10_2'])

    def test_build_assigns_primer_bases(self):
        self.assertEqual([s.primer_base for s in self.summaries], ['HOXB13', 'OLFML3'])

    def test_order_lines_and_text(self):
        hox = self.summaries[0]
        self.assertEqual(snp_summary.order_lines(hox),
                         'HOXB13_F\tTATGAGCTTCTCTCCGCCAG\nHOXB13_R\tCGCCTCCAAAGTAGCCGTAA')
        self.assertEqual(snp_summary.order_lines(self.summaries[1]), '')
        self.assertEqual(snp_summary.summary_text(hox, 'fallback').splitlines(), [
            '11:37525128 C>G (HOXB13)',
            'Статус: Надёжно',
            'Фермент: Bse1I, сайт ACTGG',
            'HOXB13_F: TATGAGCTTCTCTCCGCCAG (Tm 59.0 °C)',
            'HOXB13_R: CGCCTCCAAAGTAGCCGTAA (Tm 60.1 °C)',
            'Продукт: 475 п.н.',
            'Полосы: REF/REF 448; REF/ALT 212; 236; 448; ALT/ALT 212; 236',
        ])
        self.assertEqual(snp_summary.summary_text(self.summaries[1], 'fallback').splitlines(), [
            '1:91308651 C>T (OLFML3)',
            'Статус: Нет схемы',
            'Примечание: Не найдено сочетание для различения трёх генотипов. Попробуйте dCAPS «Для каждого SNP»',
        ])
        self.assertEqual(snp_summary.weakest_text(hox.scheme),
                         'Труднее всего различить: REF/REF–REF/ALT; оценка при ухудшении 10.5')

    def test_tables(self):
        header, rows = snp_summary.summary_table(self.summaries, 'fallback')
        self.assertEqual(len(header), 22)
        self.assertEqual((header[0], header[-1]), ('ID', 'Строка файла'))
        self.assertEqual(rows[0][:6], ['HOXB13', '11:37525128 C>G', 'NC_056064.1', 'Надёжно', 'авто', 'Bse1I'])
        self.assertEqual(rows[1][:6], ['OLFML3', '1:91308651 C>T', 'NC_056054.1', 'Нет схемы', '', None])
        self.assertEqual(rows[0][-1], 2)
        header, rows = snp_summary.primer_table(self.summaries)
        self.assertEqual(header, snp_summary.PRIMER_HEADER)
        self.assertEqual(rows, [['HOXB13_F', 'TATGAGCTTCTCTCCGCCAG', 20, 59.04, '11:37525128 C>G', 'Bse1I'],
                                ['HOXB13_R', 'CGCCTCCAAAGTAGCCGTAA', 20, 60.1, '11:37525128 C>G', 'Bse1I']])
        header, rows = snp_summary.results_table(self.summaries, HEADER, self.core_rows)
        self.assertEqual(header[-2:], ['snp_id', 'selected'])
        self.assertEqual([r[-1] for r in rows], [True, False, False])
        self.assertEqual([r[-2] for r in rows], ['HOXB13', 'HOXB13', 'OLFML3'])

    def test_manual_choice_in_tables(self):
        self.summaries[0].choose(1)
        _, rows = snp_summary.summary_table(self.summaries, 'fallback')
        self.assertEqual(rows[0][4:6], ['вручную', 'BsrDI'])
        _, rows = snp_summary.results_table(self.summaries, HEADER, self.core_rows)
        self.assertEqual([r[-1] for r in rows], [False, True, False])

    def test_english_note_summary_and_tables(self):
        i18n.set_language('en')
        self.addCleanup(i18n.set_language, 'ru')
        s = summary(schemes=[dict(enzyme='Bse1I', genotype_quality='robust'),
                             dict(enzyme='HinfI', genotype_quality='nominal_only',
                                  weakest_genotypes='REF/REF–REF/ALT')])
        s.choose(1)
        self.assertEqual(s.note('fallback'), 'Chosen by hand, best is Bse1I; '
                                             'Nominal resolution only; hardest pair: REF/REF–REF/ALT')
        header, rows = snp_summary.summary_table(self.summaries, 'fallback')
        self.assertEqual(header[:5], ['ID', 'SNP', 'Genome chromosome', 'Status', 'Choice'])
        self.assertEqual(rows[0][3:5], ['Robust', 'auto'])
        self.assertEqual(snp_summary.primer_table(self.summaries)[0][:2], ['Name', 'Sequence 5′→3′'])
        sheets = snp_summary.excel_sheets(self.summaries, 'fallback')
        self.assertEqual([x.title for x in sheets], ['Summary', 'Primers'])
        self.assertEqual(sheets[0].fills['Status']['Robust'], 'EAF3DE')       # Review Focus 3
        self.assertEqual(sheets[0].fills['Status']['Data error'], 'FCEBEB')

    def test_english_texts_for_the_lab_notebook(self):
        i18n.set_language('en')
        self.addCleanup(i18n.set_language, 'ru')
        hox = self.summaries[0]
        self.assertEqual(snp_summary.summary_text(hox, 'fallback').splitlines(), [
            '11:37525128 C>G (HOXB13)',
            'Status: Robust',
            'Enzyme: Bse1I, site ACTGG',
            'HOXB13_F: TATGAGCTTCTCTCCGCCAG (Tm 59.0 °C)',
            'HOXB13_R: CGCCTCCAAAGTAGCCGTAA (Tm 60.1 °C)',
            'Product: 475 bp',
            'Bands: REF/REF 448; REF/ALT 212; 236; 448; ALT/ALT 212; 236',
        ])
        self.assertEqual(snp_summary.summary_text(self.summaries[1], 'fallback').splitlines()[-1],
                         'Note: No primer pair and enzyme tell the three genotypes apart. Try dCAPS «For every SNP»')
        self.assertEqual(snp_summary.weakest_text(hox.scheme),
                         'Hardest to tell apart: REF/REF–REF/ALT; stressed score 10.5')

    def test_english_failure_and_error_notes(self):
        i18n.set_language('en')
        self.addCleanup(i18n.set_language, 'ru')
        self.assertEqual(summary(failure=dict(status='no_primers', reason='Primer3 не вернул подходящих пар'))
                         .note('off'), 'Primer3 returned no suitable pairs. Try a wider product size range or Tm')
        self.assertEqual(summary(rejection=dict(reason='parse_error', details='Ожидалось минимум 4 поля'))
                         .note('off'), 'Line not read: At least 4 fields expected')
        self.assertEqual(summary(rejection=dict(reason='bad_chrom', details='x')).note('off'),
                         'Chromosome not found in the genome')
        self.assertEqual(summary(rejection=dict(reason='unsupported_variant',
                                                details='Поддерживаются только одиночные SNP A/C/G/T')).note('off'),
                         'Only single-base A/C/G/T SNPs are supported')
        self.assertEqual(summary().note('off'), 'No result for this line')

    def test_results_table_for_a_subset(self):
        _, rows = snp_summary.results_table(self.summaries[1:], HEADER, self.core_rows)
        self.assertEqual([r[HEADER.index('variant')] for r in rows], ['1:91308651 C>T'])


if __name__ == '__main__':
    unittest.main()
