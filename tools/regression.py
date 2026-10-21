"""Check that two RFLP Picker source trees give identical analysis results.

Typical use: compare the working copy with the v1.1.0 release before a change
that must not alter results (refactoring, speed-ups)::

    git worktree add ../rflp-v1.1.0 v1.1.0
    python tools/regression.py --baseline ../rflp-v1.1.0 --work ../regression \
        --dataset sheep=genome.fa,variants.tsv --random 100 --jobs 8

Every tree is imported in its own process.  For each dataset and parameter
profile the canonical JSON of the header, rows and report (summary counters,
rejected variants, parameters) is compared: ``same`` means byte-identical,
``extended`` that the candidate only appends columns and report fields while
every baseline value is unchanged; anything else fails.  The log and progress
calls are compared separately.  Inputs are copied into ``--work`` so that no
index or result file is written next to the baseline tree.
"""
import argparse
import hashlib
import json
import random
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

BASE_PARAMS = dict(names2_file='', flank=300, min_frag=80, max_frag=800, delta=25, max_cuts=3,
                   prod_min=250, prod_max=600, tm_min=58.0, tm_max=62.0,
                   gel_resolution_bp=10.0, gel_resolution_pct=3.0, gel_stress_factor=1.5,
                   gain_loss_only=False, suppliers=[], use_primer3=True, joint_design=True,
                   top_results=0, primer_pairs=10, skip_norm=False, no_save_norm=True)
# Parameter sets of the article validation plus the other search modes and filters.
PROFILES = {
    'joint': {},
    'joint_wide': dict(prod_min=150, prod_max=1000),
    'joint_top10': dict(top_results=10),
    'single_pair': dict(joint_design=False),
    'window': dict(use_primer3=False),
    'filters': dict(gain_loss_only=True, max_cuts=1, delta=40, flank=500, primer_pairs=5,
                    min_frag=50, max_frag=1200, gel_resolution_bp=5.0, gel_resolution_pct=2.0,
                    gel_stress_factor=2.0, suppliers=['New england biolabs', 'Thermo fisher']),
    'skip_norm_window': dict(skip_norm=True, use_primer3=False, flank=150),
}


def sha256_text(text):
    return hashlib.sha256(text.encode('utf-8')).hexdigest()


def _canon(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True)


def _extends(base, cand):
    """cand equals base, except that dictionaries may have additional keys."""
    if isinstance(base, dict):
        return isinstance(cand, dict) and all(k in cand and _extends(v, cand[k]) for k, v in base.items())
    if isinstance(base, list):
        return (isinstance(cand, list) and len(base) == len(cand)
                and all(_extends(b, c) for b, c in zip(base, cand)))
    return _canon(base) == _canon(cand)


def compatible(base, cand):
    """Only appended columns and new report fields: every baseline value is unchanged."""
    n = len(base['header'])
    return (cand['header'][:n] == base['header'] and len(cand['rows']) == len(base['rows'])
            and all(_canon(c[:n]) == _canon(b) for b, c in zip(base['rows'], cand['rows']))
            and base['error'] == cand['error'] and _extends(base['report'], cand['report']))


# ---------------------------------------------------------------- worker side

def run_case(code_dir, case_file, out_file):
    """Run one case with the rflp_core found in code_dir (called in a subprocess)."""
    sys.dont_write_bytecode = True  # keep the compared trees free of __pycache__
    code_dir = Path(code_dir).resolve()
    sys.path.insert(0, str(code_dir))
    import rflp_core
    if Path(rflp_core.__file__).resolve().parent != code_dir:
        raise RuntimeError(f'rflp_core imported from {rflp_core.__file__}, expected {code_dir}')
    case_text = Path(case_file).read_text(encoding='utf-8')
    case = json.loads(case_text)
    log, progress, report = [], [], {}
    started = time.perf_counter()
    try:
        header, rows = rflp_core.run_rflp_gui_mode(
            case['params'], status_cb=log.append,
            progress_cb=lambda done, total: progress.append([done, total]), report=report)
        error = None
    except Exception as exc:  # the error itself is part of the compared result
        header, rows, error = [], [], f'{type(exc).__name__}: {exc}'
    elapsed = time.perf_counter() - started
    result = json.dumps(dict(header=header, rows=rows, report=report, error=error),
                        ensure_ascii=False, sort_keys=True)
    trace = json.dumps(dict(log=log, progress=progress), ensure_ascii=False)
    out = Path(out_file)
    out.with_suffix('.result.json').write_text(result, encoding='utf-8')
    out.with_suffix('.log.json').write_text(trace, encoding='utf-8')
    out.write_text(json.dumps(dict(result_sha256=sha256_text(result), log_sha256=sha256_text(trace),
                                   rows=len(rows), error=error, elapsed_seconds=elapsed,
                                   case_sha256=sha256_text(case_text))),
                   encoding='utf-8')


# ---------------------------------------------------------------- case setup

def random_variants(fasta, path, count, seed):
    """Random SNPs plus edge cases (swapped/foreign REF, contig ends, bad rows)."""
    from pyfaidx import Fasta
    rng = random.Random(seed)
    fa = Fasta(str(fasta), as_raw=True, sequence_always_upper=True)
    contigs = [(name, len(fa[name])) for name in fa.keys()]
    total = sum(length for _, length in contigs)
    rows = []

    def pick_contig():
        x = rng.randrange(total)
        for name, length in contigs:
            if x < length:
                return name, length
            x -= length
        return contigs[-1]

    def snp(name, pos):
        ref = fa[name][pos - 1]
        alt = rng.choice([b for b in 'ACGT' if b != ref])
        return [name, pos, ref, alt]

    for _ in range(count):
        rows.append(snp(*((lambda n, l: (n, rng.randint(1, l)))(*pick_contig()))))
    name, length = max(contigs, key=lambda c: c[1])
    short, short_len = min(contigs, key=lambda c: c[1])
    edge = [snp(name, 1), snp(name, 7), snp(name, 150), snp(name, length - 4), snp(name, length),
            snp(short, short_len // 2)]
    swapped = snp(name, rng.randint(1000, length - 1000))
    swapped[2], swapped[3] = swapped[3], swapped[2]
    comp = snp(name, rng.randint(1000, length - 1000))
    comp[2] = comp[2].translate(str.maketrans('ACGT', 'TGCA'))
    rows += edge + [swapped, comp,
                    [name, length + 10, 'A', 'G'], ['no_such_contig', 100, 'A', 'G'],
                    [name, 5000, 'AC', 'A'], [name, 5001, 'A', 'N'], [name, 'x', 'A', 'G']]
    fa.close()
    with open(path, 'w', encoding='utf-8', newline='') as fh:
        fh.write('chrom\tpos\tref\talt\n')
        for row in rows:
            fh.write('\t'.join(map(str, row)) + '\n')


def build_cases(args, work):
    inputs = work / 'inputs'
    inputs.mkdir(parents=True, exist_ok=True)
    datasets = []
    for spec in args.dataset:
        name, _, files = spec.partition('=')
        fasta, variants, *names = files.split(',')
        names_file = names[0] if names else ''
        fasta = Path(fasta).resolve()
        if args.copy_fasta:
            local = inputs / fasta.name
            if not local.exists():
                shutil.copyfile(fasta, local)
            fasta = local
        # Index once here: parallel workers would otherwise race to write the .fai.
        from pyfaidx import Fasta
        Fasta(str(fasta)).close()
        local_variants = inputs / f'{name}.tsv'
        shutil.copyfile(variants, local_variants)
        datasets.append((name, fasta, local_variants, names_file))
        if args.random:
            rnd = inputs / f'{name}_random.tsv'
            if not rnd.exists():
                random_variants(fasta, rnd, args.random, seed=f'{name}:{args.random}')
            datasets.append((f'{name}_random', fasta, rnd, names_file))
    cases = []
    profiles = args.profile or list(PROFILES)
    for name, fasta, variants, names_file in datasets:
        for profile in profiles:
            params = dict(BASE_PARAMS, **PROFILES[profile], fasta=str(fasta),
                          input_file=str(variants), names2_file=str(names_file))
            cases.append(dict(name=f'{name}.{profile}', params=params))
    return cases


# ---------------------------------------------------------------- driver side

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--baseline', required=True, help='Source tree of the reference version')
    parser.add_argument('--candidate', default=str(ROOT), help='Source tree to check (default: this one)')
    parser.add_argument('--work', required=True, help='Directory for inputs and results')
    parser.add_argument('--dataset', action='append', required=True,
                        help='NAME=FASTA,VARIANTS[,NAMES]; repeat for several genomes')
    parser.add_argument('--random', type=int, default=0, help='Also test N random SNPs per genome')
    parser.add_argument('--profile', action='append', choices=list(PROFILES),
                        help='Parameter profiles to run (default: all)')
    parser.add_argument('--copy-fasta', action='store_true',
                        help='Copy small FASTA files into --work (their .fai is written there)')
    parser.add_argument('--jobs', type=int, default=1, help='Parallel processes (timings then overlap)')
    parser.add_argument('--reuse-baseline', action='store_true',
                        help='Keep baseline results of unchanged cases from a previous run')
    parser.add_argument('--worker', nargs=3, metavar=('CODE', 'CASE', 'OUT'), help=argparse.SUPPRESS)
    if '--worker' in sys.argv:
        i = sys.argv.index('--worker')
        return run_case(*sys.argv[i + 1:i + 4])
    args = parser.parse_args()
    work = Path(args.work).resolve()
    trees = {'baseline': Path(args.baseline).resolve(), 'candidate': Path(args.candidate).resolve()}
    cases = build_cases(args, work)
    for label in trees:
        (work / label).mkdir(parents=True, exist_ok=True)
    (work / 'cases').mkdir(exist_ok=True)
    for case in cases:
        (work / 'cases' / f"{case['name']}.json").write_text(
            json.dumps(case, ensure_ascii=False, indent=1), encoding='utf-8')

    def launch(job):
        label, case = job
        out = work / label / f"{case['name']}.json"
        case_file = work / 'cases' / f"{case['name']}.json"
        if label == 'baseline' and args.reuse_baseline and out.exists():
            previous = json.loads(out.read_text(encoding='utf-8'))
            if previous.get('case_sha256') == sha256_text(case_file.read_text(encoding='utf-8')):
                return label, case['name'], previous
        proc = subprocess.run([sys.executable, str(Path(__file__).resolve()), '--worker',
                               str(trees[label]), str(case_file), str(out)],
                              capture_output=True, text=True, encoding='utf-8', errors='replace')
        if proc.returncode:
            raise RuntimeError(f"{label} {case['name']} failed:\n{proc.stderr}")
        return label, case['name'], json.loads(out.read_text(encoding='utf-8'))

    jobs = [(label, case) for case in cases for label in trees]
    results = {}
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        for label, name, summary in pool.map(launch, jobs):
            results[label, name] = summary
            print(f"{name:32s} {label:9s} rows={summary['rows']:6d} {summary['elapsed_seconds']:9.2f}s",
                  flush=True)

    report, failed = [], 0
    for case in cases:
        b, c = results['baseline', case['name']], results['candidate', case['name']]
        same, same_log = b['result_sha256'] == c['result_sha256'], b['log_sha256'] == c['log_sha256']
        if same:
            verdict = 'same'
        else:
            load = lambda label: json.loads((work / label / f"{case['name']}.result.json").read_text(encoding='utf-8'))  # noqa: E731
            verdict = 'extended' if compatible(load('baseline'), load('candidate')) else 'DIFFERS'
        failed += verdict == 'DIFFERS'
        report.append(dict(case=case['name'], rows=b['rows'], result=verdict, identical_log=same_log,
                           baseline_seconds=b['elapsed_seconds'], candidate_seconds=c['elapsed_seconds'],
                           result_sha256=b['result_sha256'] if same else None))
    (work / 'comparison.json').write_text(json.dumps(dict(trees={k: str(v) for k, v in trees.items()},
                                                          cases=report), indent=1), encoding='utf-8')
    print(f"\n{'case':32s} {'rows':>6s} {'result':>9s} {'log':>9s} {'baseline':>9s} {'candidate':>9s} {'speed-up':>8s}")
    for r in report:
        speed = r['baseline_seconds'] / max(r['candidate_seconds'], 1e-9)
        print(f"{r['case']:32s} {r['rows']:6d} {r['result']:>9s} "
              f"{'same' if r['identical_log'] else 'differs':>9s} {r['baseline_seconds']:8.1f}s "
              f"{r['candidate_seconds']:8.1f}s {speed:7.1f}x")
    print(f"\n{len(report) - failed}/{len(report)} cases match the baseline "
          "(same: byte-identical; extended: only new columns/report fields)")
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
