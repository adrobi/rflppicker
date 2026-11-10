import locale
import random
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from Bio.Restriction import AllEnzymes

import rflp_core
from tests.helpers import write_text
from rflp_core import (
    AnalysisCancelled,
    cut_positions,
    diagnostic_delta,
    normalize_variants_for_rflp,
    parse_input,
    run_rflp_gui_mode,
    snp_site_enzymes,
)


class NumericLocaleTests(unittest.TestCase):
    def test_primer3_works_with_decimal_comma_locale(self):
        # Qt switches LC_NUMERIC to the system locale; Primer3 then reads "60.0" as invalid.
        saved = locale.setlocale(locale.LC_NUMERIC)
        self.addCleanup(locale.setlocale, locale.LC_NUMERIC, saved)
        for name in ('ru_RU.UTF-8', 'de_DE.UTF-8', 'Russian_Russia'):
            try:
                locale.setlocale(locale.LC_NUMERIC, name)
                break
            except locale.Error:
                continue
        else:
            self.skipTest('no decimal-comma locale installed')
        self.assertEqual(locale.localeconv()['decimal_point'], ',')
        rng = random.Random(1)
        sequence = ''.join(rng.choice('ACGT') for _ in range(601))
        result = rflp_core.primer3_pick(sequence, 300, 120, 300, 55.0, 65.0, num_return=2)
        self.assertIsNotNone(result, rflp_core.primer3_pick.last_error)


class CoreTests(unittest.TestCase):
    def test_csv_quotes_and_rejected_multibase(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            fasta = root / 'toy.fa'
            fasta.write_text('>1\nACGTACGTACGT\n', encoding='utf-8')
            variants = root / 'variants.csv'
            variants.write_text('chrom,pos,ref,alt\n1,2,"C","T"\n1,3,GT,AA\n', encoding='utf-8')
            report = {}
            parsed = parse_input(str(variants), report=report)
            self.assertEqual(parsed[0], ('1', 2, 'C', 'T'))
            cleaned, _, _, _ = normalize_variants_for_rflp(str(fasta), str(variants), None, save_files=False, report=report)
            self.assertEqual(cleaned, [('1', 2, 'C', 'T')])
            self.assertTrue(any(x['reason'] == 'unsupported_variant' for x in report['rejected']))

    def test_skip_norm_still_checks_ref(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            fasta = root / 'toy.fa'; fasta.write_text('>1\nACGTACGT\n', encoding='utf-8')
            variants = root / 'v.tsv'; variants.write_text('chrom\tpos\tref\talt\n1\t2\tA\tG\n', encoding='utf-8')
            report = {}
            cleaned, _, _, _ = normalize_variants_for_rflp(str(fasta), str(variants), None, skip_norm=True, save_files=False, report=report)
            self.assertEqual(cleaned, [])
            self.assertEqual(report['rejected'][0]['reason'], 'ref_not_match')

    def test_rejected_variants_carry_the_line_of_the_input_file(self):
        with tempfile.TemporaryDirectory() as d:
            fasta = write_text(d, 'toy.fa', '>1\nACGTACGTACGT\n')
            variants = write_text(d, 'variants.tsv', 'chrom\tpos\tref\talt\n\nbad\tline\n1\t2\tN\tG\n1\t3\tG\tA\n')
            for skip_norm, reason in ((True, 'ref_not_match'), (False, 'bad_ref_chars')):
                with self.subTest(skip_norm=skip_norm):
                    report = {}
                    cleaned, _, _, _ = normalize_variants_for_rflp(fasta, variants, None, skip_norm=skip_norm,
                                                                   save_files=False, report=report)
                    self.assertEqual(cleaned, [('1', 3, 'G', 'A')])
                    self.assertEqual([(item['line'], item['reason']) for item in report['rejected']],
                                     [(3, 'parse_error'), (4, reason)])
                    self.assertEqual(list(report['rejected'][1])[0], 'line')

    def test_optional_params_take_their_defaults_from_default_params(self):
        with tempfile.TemporaryDirectory() as d:
            fasta = write_text(d, 'toy.fa', '>1\nACGTACGTACGT\n')
            variants = write_text(d, 'variants.tsv', 'chrom\tpos\tref\talt\n1\t2\tC\tT\n')
            given = {key: value for key, value in rflp_core.DEFAULT_PARAMS.items()
                     if key not in rflp_core.OPTIONAL_PARAMS}
            given.update(fasta=fasta, input_file=variants, use_primer3=False, flank=5, min_frag=1)
            report = {}
            run_rflp_gui_mode(given, report=report)
        self.assertEqual(len(rflp_core.OPTIONAL_PARAMS), 19)
        self.assertEqual(set(report['parameters']), set(given) | set(rflp_core.OPTIONAL_PARAMS))
        for key in rflp_core.OPTIONAL_PARAMS:
            self.assertEqual(report['parameters'][key], rflp_core.DEFAULT_PARAMS[key], key)

    def test_skip_norm_can_be_cancelled(self):
        with tempfile.TemporaryDirectory() as d:
            fasta = write_text(d, 'toy.fa', '>1\nACGTACGTACGT\n')
            variants = write_text(d, 'variants.tsv', 'chrom\tpos\tref\talt\n1\t2\tC\tT\n')
            with self.assertRaises(AnalysisCancelled):
                normalize_variants_for_rflp(fasta, variants, None, skip_norm=True, save_files=False,
                                            cancel_cb=lambda: True)

    def test_threshold_is_applied_for_gain_loss_and_reported(self):
        self.assertEqual(diagnostic_delta([500], [100, 400]), 400)
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            fasta = root / 'toy.fa'; fasta.write_text('>1\nACGTACGTACGTACGTACGT\n', encoding='utf-8')
            variants = root / 'v.tsv'; variants.write_text('chrom\tpos\tref\talt\n1\t2\tC\tT\n', encoding='utf-8')
            params = dict(fasta=str(fasta), input_file=str(variants), names2_file='', flank=4,
                          min_frag=1, max_frag=1000, delta=1, gain_loss_only=False, max_cuts=3,
                          suppliers=[], use_primer3=False, prod_min=50, prod_max=100,
                          tm_min=50, tm_max=70, skip_norm=False, no_save_norm=True)
            header, rows = run_rflp_gui_mode(params)
            self.assertIn('status', header); self.assertIn('delta_threshold_bp', header)
            self.assertTrue(rows)
            self.assertEqual(len(header), len(rows[0]))

    def test_cancellation_is_honored(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            fasta = root / 'toy.fa'; fasta.write_text('>1\nACGTACGTACGT\n', encoding='utf-8')
            variants = root / 'v.tsv'; variants.write_text('chrom\tpos\tref\talt\n1\t2\tC\tT\n', encoding='utf-8')
            params = dict(fasta=str(fasta), input_file=str(variants), names2_file='', flank=4,
                          min_frag=1, max_frag=1000, delta=1, gain_loss_only=False, max_cuts=3,
                          suppliers=[], use_primer3=False, prod_min=50, prod_max=100,
                          tm_min=50, tm_max=70, skip_norm=False, no_save_norm=True)
            with self.assertRaises(AnalysisCancelled):
                run_rflp_gui_mode(params, cancel_cb=lambda: True)

    def test_site_patterns_have_fixed_width(self):
        # snp_site_enzymes relies on every site being a look-ahead of enz.size bases.
        lookahead = re.compile(r'\(\?=\(\?P<\w+>((?:\[[A-Z]+\]|[A-Z.])+)\)\)')
        for enz in AllEnzymes:
            parts = enz.compsite.pattern.split('|')
            self.assertIn(len(parts), (1, 2), enz)
            for part in parts:
                match = lookahead.fullmatch(part)
                self.assertIsNotNone(match, (enz, part))
                self.assertEqual(len(re.findall(r'\[[A-Z]+\]|[A-Z.]', match.group(1))), enz.size, enz)

    def test_enzymes_without_snp_site_cut_both_alleles_equally(self):
        rng = random.Random(7)
        enzymes = sorted(AllEnzymes, key=lambda e: e.__name__)
        for _ in range(6):
            seq = ''.join(rng.choice('ACGT') for _ in range(90))
            off = rng.choice([0, 1, 5, 45, 84, 89])
            altseq = seq[:off] + rng.choice([b for b in 'ACGT' if b != seq[off]]) + seq[off + 1:]
            found = set(snp_site_enzymes(enzymes, seq, altseq, off))
            self.assertTrue(found)
            amplicons = {(0, len(seq)), (off, off + 1), (max(0, off - 1), min(len(seq), off + 2)),
                         (rng.randint(0, off), rng.randint(off + 1, len(seq)))}
            for enz in enzymes:
                if enz in found:
                    continue
                for a, b in amplicons:
                    self.assertEqual(cut_positions(enz, seq[a:b]), cut_positions(enz, altseq[a:b]),
                                     (enz, seq, off, a, b))


class InputTableTests(unittest.TestCase):
    def setUp(self):
        self._dir = tempfile.TemporaryDirectory()
        self.tmp = self._dir.name

    def tearDown(self):
        self._dir.cleanup()

    def test_read_input_table_detects_columns_and_id(self):
        path = write_text(self.tmp, 'v.tsv', 'chrom\tpos\tref\talt\tName\n1\t10\tA\tG\tDGAT1\n\n1\t11\tC\n')
        columns, lines, delim = rflp_core.read_input_table(path)
        self.assertEqual(delim, '\t')
        self.assertEqual(columns, {'chrom': 0, 'pos': 1, 'ref': 2, 'alt': 3, 'id': 4})
        self.assertEqual(lines, [(2, ['1', '10', 'A', 'G', 'DGAT1']), (4, ['1', '11', 'C'])])

    def test_read_input_table_headerless_and_template_rows(self):
        path = write_text(self.tmp, 'v.csv', '1,10,A,G\n2,20,C,T\n')
        columns, lines, _ = rflp_core.read_input_table(path)
        self.assertEqual(set(columns.values()), {None})
        self.assertEqual([n for n, _ in lines], [1, 2])
        path = write_text(self.tmp, 'h.csv', 'chrom,pos,ref,alt\nдолжно быть: 1,10,A,G\n1,10,A,G\n')
        self.assertEqual([n for n, _ in rflp_core.read_input_table(path)[1]], [3])

    def test_read_input_table_missing_file(self):
        self.assertEqual(rflp_core.read_input_table(str(Path(self.tmp) / 'missing.tsv')), ({}, [], ''))

    def test_parse_fields(self):
        columns = {'chrom': None, 'pos': None, 'ref': None, 'alt': None, 'id': None}
        self.assertEqual(rflp_core.parse_fields(['chr1 ', '10', 'a', 'g '], columns), ('1', 10, 'A', 'G'))
        with self.assertRaises(ValueError):
            rflp_core.parse_fields(['1', 'x', 'A', 'G'], columns)


class ContigNameTests(unittest.TestCase):
    """Genomes whose contigs are named chr1, chr2, … work without a names file."""

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        rng = random.Random(3)
        sequence = ''.join(rng.choice('ACGT') for _ in range(400))
        self.sequence = sequence[:197] + 'GAATTC' + sequence[203:]

    def run_core(self, fasta_contig, input_chrom, report=None, **extra):
        fasta = write_text(self.temp.name, 'toy.fa', f'>{fasta_contig}\n{self.sequence}\n')
        variants = write_text(self.temp.name, 'variants.tsv', f'chrom\tpos\tref\talt\n{input_chrom}\t200\tA\tG\n')
        params = dict(rflp_core.DEFAULT_PARAMS, fasta=fasta, input_file=variants, use_primer3=False,
                      flank=150, min_frag=20, **extra)
        header, rows = run_rflp_gui_mode(params, report=report)
        return [dict(zip(header, row)) for row in rows]

    def test_chr_prefixed_fasta_without_names_file(self):
        for chrom in ('chr1', '1'):
            with self.subTest(input_chrom=chrom):
                records = self.run_core('chr1', chrom)
                self.assertTrue(records)
                self.assertEqual({r['variant'] for r in records}, {'1:200 A>G'})
                self.assertEqual({r['mapped_id'] for r in records}, {'chr1'})
                self.assertIn('EcoRI', {r['enzyme'] for r in records})

    def test_bare_fasta_with_chr_input_is_unchanged(self):
        records = self.run_core('1', 'chr1')
        self.assertEqual({r['variant'] for r in records}, {'1:200 A>G'})
        self.assertEqual({r['mapped_id'] for r in records}, {'1'})

    def test_chr_inside_a_contig_name_is_kept(self):
        self.assertEqual(rflp_core.clean_chrom(' scaffold_chr5 '), 'scaffold_chr5')
        self.assertEqual(rflp_core.clean_chrom('chr5'), '5')
        records = self.run_core('scaffold_chr5', 'scaffold_chr5')
        self.assertEqual({r['mapped_id'] for r in records}, {'scaffold_chr5'})

    def test_unknown_contig_is_still_rejected_with_its_name(self):
        report = {}
        self.assertEqual(self.run_core('chr1', '7', report=report), [])
        self.assertEqual(report['rejected'][0]['reason'], 'bad_chrom')
        self.assertIn("'7'", report['rejected'][0]['details'])

    def test_background_vcf_contig_names_with_and_without_chr(self):
        for fasta_contig, vcf_contig in (('chr1', 'chr1'), ('chr1', '1'), ('1', 'chr1')):
            with self.subTest(fasta=fasta_contig, vcf=vcf_contig):
                vcf = write_text(self.temp.name, 'known.vcf', '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                                 f'{vcf_contig}\t150\t.\t{self.sequence[149]}\tN\t.\t.\t.\n')
                report = {}
                self.run_core(fasta_contig, '1', report=report, background_vcf=vcf)
                self.assertEqual(report['summary']['background_variants_loaded'], 1)


class ReviewFixTests(unittest.TestCase):
    def test_swap_with_a_stray_comma_in_alt_as_in_1_3_0(self):
        with tempfile.TemporaryDirectory() as d:
            fasta = write_text(d, 'toy.fa', '>1\nACGTACGTACGT\n')
            variants = write_text(d, 'variants.tsv', 'chrom\tpos\tref\talt\n1\t3\tT\tG,\n1\t2\tA\t,C\n')
            cleaned, _, _, _ = normalize_variants_for_rflp(fasta, variants, None, save_files=False, report={})
        self.assertEqual(cleaned, [('1', 3, 'G', 'T'), ('1', 2, 'C', 'A')])

    def test_primer3_tables_are_read_with_a_decimal_point(self):
        # Primer3 parses its thermodynamic tables when it is imported: with a decimal
        # comma active at that moment the hairpin and dimer values were silently wrong.
        code = ("import locale, sys\n"
                "if sys.argv[1] != 'C':\n"
                "    locale.setlocale(locale.LC_ALL, sys.argv[1])\n"
                "import rflp_core, primer3\n"
                "print(repr(primer3.bindings.calc_hairpin('CCCCCATCCGATCAGGGGG').dg))\n")
        root = str(Path(rflp_core.__file__).resolve().parent)

        def hairpin(name):
            done = subprocess.run([sys.executable, '-c', code, name], cwd=root, capture_output=True, text=True)
            return done.returncode, done.stdout.strip()

        for name in ('ru_RU.UTF-8', 'de_DE.UTF-8', 'Russian_Russia'):
            code_comma, value = hairpin(name)
            if code_comma == 0:
                break
        else:
            self.skipTest('no decimal-comma locale installed')
        self.assertEqual(value, hairpin('C')[1])


if __name__ == '__main__':
    unittest.main()
