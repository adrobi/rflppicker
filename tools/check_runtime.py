"""Functional smoke test for the packaged application."""

import random
from pathlib import Path
from tempfile import TemporaryDirectory

from Bio.bgzf import BgzfWriter
from PyQt6.QtWidgets import QApplication

import lab_sheet
import project_file
import reporting
import rflp_core as core
from snp_summary import build_summaries, read_input_records, results_table


def check_runtime() -> int:
    if not core._PRIMER3_AVAILABLE:
        raise RuntimeError(f'Primer3 unavailable: {core._PRIMER3_IMPORT_ERROR}')
    app = QApplication.instance() or QApplication(['RFLP Picker'])
    with TemporaryDirectory() as work:
        root = Path(work)
        fasta, vcf, saved = root / 'genome.fa.gz', root / 'targets.vcf', root / 'run.rflp'
        rng = random.Random(19)
        sequence = ''.join(rng.choice('ACGT') for _ in range(400))
        sequence = sequence[:197] + 'GAATTC' + sequence[203:]
        with BgzfWriter(fasta, 'wb') as stream:
            stream.write(f'>chr1\n{sequence}\n'.encode('ascii'))
        vcf.write_text('##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
                       'chr1\t200\ttest\tA\tG\t.\tPASS\t.\n', encoding='ascii')
        params = dict(core.DEFAULT_PARAMS, fasta=str(fasta), input_file=str(vcf),
                      prod_min=200, prod_max=500, tm_min=50, tm_max=70, primer_pairs=4,
                      enzyme_mode='only', lab_enzymes=['EcoRI'], max_cuts=1,
                      snp_end_distance=0, exclude_primer_sites=False, mask_repeats=False,
                      double_digest=False, allele_specific=False, dcaps='off', specificity='off')
        report = {}
        header, rows = core.run_analysis(params, report=report)
        records, has_ids = read_input_records(str(vcf))
        summaries = build_summaries(records, header, rows, report.get('rejected', []), report=report)
        if not (summaries and summaries[0].scheme and summaries[0].scheme['primers_left']):
            raise RuntimeError('The test variant produced no assay with Primer3 primers')
        metadata = project_file.run_metadata(params, report)
        payloads = project_file.capture_payloads(summaries, params)
        project_file.save_project(saved, header=header, rows=rows, report=report, params=params,
                                  metadata=metadata, has_ids=has_ids, summaries=summaries, payloads=payloads)
        fasta.unlink()
        reopened = project_file.open_project(saved)
        table = results_table(reopened.summaries, reopened.header, reopened.rows,
                              row_lines=reopened.report['result_source_lines'])
        reporting.export_excel(root / 'run.xlsx', *table, metadata=reopened.metadata)
        pdf = root / 'run.pdf'
        if not lab_sheet.export_lab_sheets(
            str(pdf), reopened.summaries, '', reopened.params, 'off',
            payload_lookup=lambda scheme: project_file.cached_payload(reopened.payloads, scheme),
            metadata=reopened.metadata,
        ) or not pdf.read_bytes().startswith(b'%PDF'):
            raise RuntimeError('Could not export a PDF from the saved project')
    print(f'Runtime OK (Qt {app.platformName()}): Primer3, BGZF FASTA, VCF, project, XLSX and PDF')
    return 0
