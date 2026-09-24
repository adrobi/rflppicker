"""Reproducible comparison on user-supplied SNPs (no scientific data invented)."""
import argparse
import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path

from reporting import collect_versions, export_csv
from rflp_core import run_rflp_gui_mode, MODEL_VERSION


def sha256(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def benchmark(params, output):
    output = Path(output)
    # Do not overwrite a previous experiment, even with the same parameters.
    output.mkdir(parents=True, exist_ok=False)
    manifest = dict(timestamp=datetime.now(timezone.utc).isoformat(),
                    model=MODEL_VERSION, versions=collect_versions(),
                    input_sha256=sha256(params['input_file']),
                    fasta_sha256=sha256(params['fasta']),
                    names_sha256=sha256(params['names2_file']) if params.get('names2_file') else None,
                    source_sha256={name: sha256(Path(__file__).with_name(name))
                                   for name in ('rflp_core.py', 'assay_scoring.py', 'benchmark_design.py')},
                    note='Timing is a single run per mode; rerun in varied order for publication.',
                    runs={})
    summary = []
    for mode, joint, pairs in (('legacy_single', False, 1),
                               ('joint_one', True, 1),
                               ('joint_many', True, params.get('primer_pairs', 10))):
        actual = dict(params, use_primer3=True, joint_design=joint,
                      primer_pairs=pairs, top_results=0, no_save_norm=True)
        report = {}
        started = time.perf_counter()
        header, rows = run_rflp_gui_mode(actual, report=report)
        elapsed = time.perf_counter() - started
        records = [dict(zip(header, row)) for row in rows]
        hits = [r for r in records if r['enzyme']]
        nominal = {r['variant'] for r in hits if r['genotype_margin'] > 1}
        robust = {r['variant'] for r in hits if r['worst_margin'] > 1}
        summary.append([mode, report['summary'].get('accepted', 0), len(hits),
                        len(nominal), len(robust), report['summary']['pairs_evaluated'],
                        elapsed, report['summary']['enzyme_errors']])
        export_csv(output / f'{mode}.csv', header, rows)
        manifest['runs'][mode] = dict(report, elapsed_seconds=elapsed)
        print(f'{mode}: nominal={len(nominal)}, robust={len(robust)}, {elapsed:.2f}s', flush=True)
    export_csv(output / 'summary.csv',
               ['mode', 'accepted_snps', 'candidate_assays', 'snps_nominal', 'snps_robust',
                'primer_pairs_evaluated', 'elapsed_seconds', 'enzyme_errors'], summary)
    (output / 'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--variants', required=True)
    parser.add_argument('--names', default='')
    parser.add_argument('--output', required=True, help='New directory for this experiment')
    parser.add_argument('--pairs', type=int, default=10)
    parser.add_argument('--flank', type=int, default=300)
    parser.add_argument('--min-frag', type=int, default=80)
    parser.add_argument('--max-frag', type=int, default=800)
    parser.add_argument('--delta', type=int, default=25)
    parser.add_argument('--max-cuts', type=int, default=3)
    parser.add_argument('--prod-min', type=int, default=250)
    parser.add_argument('--prod-max', type=int, default=600)
    parser.add_argument('--tm-min', type=float, default=58)
    parser.add_argument('--tm-max', type=float, default=62)
    parser.add_argument('--resolution-bp', type=float, default=10)
    parser.add_argument('--resolution-pct', type=float, default=3)
    parser.add_argument('--stress-factor', type=float, default=1.5)
    args = parser.parse_args()
    benchmark(dict(fasta=str(Path(args.fasta).resolve()), input_file=str(Path(args.variants).resolve()),
                   names2_file=str(Path(args.names).resolve()) if args.names else '',
                   primer_pairs=args.pairs, flank=args.flank, min_frag=args.min_frag,
                   max_frag=args.max_frag, delta=args.delta, max_cuts=args.max_cuts,
                   prod_min=args.prod_min, prod_max=args.prod_max, tm_min=args.tm_min,
                   tm_max=args.tm_max, gel_resolution_bp=args.resolution_bp,
                   gel_resolution_pct=args.resolution_pct, gel_stress_factor=args.stress_factor,
                   gain_loss_only=False, suppliers=[]), args.output)


if __name__ == '__main__':
    main()
