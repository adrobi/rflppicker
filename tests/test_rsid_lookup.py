import gzip
import tempfile
import unittest
from pathlib import Path

from rsid_lookup import lookup_rsid


VCF_COLUMNS = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'


def eva_hit(build='GCA_001433935.1', taxid=4530, contig='AP014957.1', pos=4, ref='T', alt='C'):
    return {'data': dict(referenceSequenceAccession=build, taxonomyAccession=taxid, contig=contig,
                         start=pos, referenceAllele=ref, alternateAllele=alt)}


class RsidLookupTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.fasta = self.root / 'rice.fa'
        self.fasta.write_text('>NC_029256.1\nACGTACGT\n>chr2\nACGTACGT\n', encoding='utf-8')
        self.vcf = self.root / 'variants.vcf.gz'
        self.report = self.root / 'report.txt'
        self.report.write_text('# Assembly name: IRGSP-1.0\n# Organism name: Oryza sativa (rice)\n# Taxid: 4530\n'
                               '# GenBank assembly accession: GCA_001433935.1\n'
                               '# RefSeq assembly accession: GCF_001433935.1\n'
                               '1\tassembled-molecule\t1\tChromosome\tAP014957.1\t=\tNC_029256.1\tPrimary Assembly\t8\tna\n'
                               '2\tassembled-molecule\t2\tChromosome\tAP014958.1\t<>\tNC_029257.1\tPrimary Assembly\t8\tna\n',
                               encoding='utf-8')

    def write_vcf(self, *rows, metadata=''):
        with gzip.open(self.vcf, 'wt', encoding='utf-8') as stream:
            stream.write('##fileformat=VCFv4.3\n' + metadata + VCF_COLUMNS)
            for record in rows:
                stream.write(record + '\t.\tPASS\t.\n')

    def lookup(self, **overrides):
        params = dict(fasta=self.fasta, taxid=4530, assembly='GCF_001433935.1',
                      local_vcf=self.vcf, names={'1': 'NC_029256.1'})
        params.update(overrides)
        return lookup_rsid('RS42', **params)

    def test_local_semicolon_ids_deduplicate_and_merge_alt(self):
        self.write_vcf('1\t4\tss3;rs42\tT\tG,C', 'NC_029256.1\t4\trs42\tT\tC',
                       '2\t4\trs420\tT\tA')
        result = self.lookup()
        self.assertEqual((result.rsid, result.status), ('rs42', 'accepted'))
        self.assertEqual((result.match.taxid, result.match.assembly, result.match.contig,
                          result.match.position, result.match.ref, result.match.alts),
                         (4530, 'GCF_001433935.1', 'NC_029256.1', 4, 'T', ('G', 'C')))

    def test_local_rejects_mismatched_metadata(self):
        self.write_vcf('1\t4\trs42\tT\tC', metadata='##reference=file:///GCF_001433935.2_genomic.fna\n')
        self.assertEqual(self.lookup().status, 'assembly_mismatch')
        self.write_vcf('1\t4\trs42\tT\tC', metadata='##reference=file:///GCF_001433935_genomic.fna\n')
        self.assertEqual(self.lookup().status, 'assembly_mismatch')
        self.write_vcf('1\t4\trs42\tT\tC', metadata='##reference=IRGSP-0.9\n')
        self.assertEqual(self.lookup(assembly_report=self.report).status, 'assembly_mismatch')
        self.write_vcf('1\t4\trs42\tT\tC', metadata='##reference=IRGSP-1.0\n')
        self.assertEqual(self.lookup(assembly_report=self.report).status, 'accepted')
        self.write_vcf('1\t4\trs42\tT\tC', metadata='##contig=<ID=1,assembly=GCF_001433935.1,species=NCBITaxon:9940>\n')
        self.assertEqual(self.lookup().status, 'species_mismatch')
        self.write_vcf('1\t4\trs42\tT\tC', metadata='##contig=<ID=1,assembly=IRGSP-1.0,species="Oryza sativa">\n')
        self.assertEqual(self.lookup(assembly_report=self.report).status, 'accepted')

    def test_report_taxid_mismatch_rejects_local_vcf_without_taxid(self):
        self.write_vcf('1\t4\trs42\tT\tC', metadata='##reference=GCF_001433935.1\n')
        self.assertEqual(self.lookup(taxid=9940, assembly_report=self.report).status, 'species_mismatch')

    def test_rsid_export_preserves_exact_contig_and_ref_in_analysis(self):
        from rsid_export import write_rsid_target
        import rflp_core as core

        self.fasta.write_text('>chr1\nATCG\n>1\nGTCA\n')
        self.write_vcf('chr1\t1\trs42\tA\tG')
        result = self.lookup(names={})
        resolved = self.root / 'resolved.vcf'
        write_rsid_target(resolved, result, self.fasta, self.vcf)
        params = dict(core.DEFAULT_PARAMS, fasta=str(self.fasta), input_file=str(resolved), use_primer3=False)
        header, rows = core.run_analysis(params)
        row = dict(zip(header, rows[0]))
        self.assertEqual((row['variant'], row['mapped_id']), ('chr1:1 A>G', 'chr1'))

    def test_ambiguity_precedes_ref_validation(self):
        self.write_vcf('1\t4\trs42\tT\tC', 'chr2\t5\trs42\tT\tG')
        result = self.lookup()
        self.assertEqual(result.status, 'ambiguous')
        self.assertIsNone(result.match)
        self.assertEqual(len(result.matches), 2)

    def test_unavailable_contig_mismatch_symbolic_and_absent_id(self):
        for row, expected in [('3\t4\trs42\tT\tC', 'contig_not_found'), ('1\t4\trs42\tA\tC', 'ref_mismatch'),
                              ('1\t4\trs42\tT\t<DEL>', 'unsupported_allele'), ('1\t4\trs420\tT\tC', 'no_match')]:
            with self.subTest(expected=expected):
                self.write_vcf(row)
                self.assertEqual(self.lookup().status, expected)

    def test_eva_filters_exact_species_and_build_and_uses_verified_pair(self):
        urls = []

        def fetch(url):
            urls.append(url)
            return [eva_hit(), eva_hit(), eva_hit(alt='G'), eva_hit(taxid=39947, pos=8),
                    eva_hit(build='GCA_000005425.2', pos=5)]

        result = self.lookup(local_vcf=None, fetch_json=fetch, assembly_report=self.report)
        self.assertTrue(urls[0].endswith('/clustered-variants/42/submitted'))
        self.assertEqual(result.status, 'accepted')
        self.assertEqual((result.match.contig, result.match.source_contig, result.match.alts),
                         ('NC_029256.1', 'AP014957.1', ('C', 'G')))
        mismatch = self.lookup(local_vcf=None, fetch_json=lambda _: [eva_hit(contig='AP014958.1')],
                               assembly_report=self.report)
        self.assertEqual(mismatch.status, 'assembly_mismatch')
        self.assertEqual(self.lookup(local_vcf=None, fetch_json=lambda _: [eva_hit()]).status, 'assembly_mismatch')

    def test_invalid_queries_and_incomplete_online_payload(self):
        for value in ('42', 'rs0', 'rs42/anything'):
            with self.subTest(value=value), self.assertRaises(ValueError):
                lookup_rsid(value, fasta=self.fasta, taxid=4530, assembly='GCF_001433935.1')
        with self.assertRaises(ValueError):
            self.lookup(assembly='GCF_001433935')
        self.assertEqual(self.lookup(local_vcf=None, fetch_json=lambda _: []).status, 'no_match')
        self.assertEqual(self.lookup(local_vcf=None, fetch_json=lambda _: {'data': []}).status, 'source_error')
