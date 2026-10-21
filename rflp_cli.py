"""RFLP Picker из командной строки: тот же анализ, что в окне программы.

Пример:
    python rflp_cli.py --fasta genome.fa --variants snps.tsv --out results.xlsx

Результат сохраняется в CSV или Excel (по расширению --out). В Excel, как и
при экспорте из окна, добавляются листы Run (параметры, версии, сводка) и
Rejected (отклонённые варианты). Журнал выводится в stderr.
"""
import argparse
import sys
from datetime import datetime
from pathlib import Path

import rflp_core as core
from reporting import export_csv, export_excel

PRESETS = {'recommended': core.RECOMMENDED_PARAMS, 'v1.1': core.DEFAULT_PARAMS}


def build_parser():
    parser = argparse.ArgumentParser(
        prog='rflp_cli', description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version=f'RFLP Picker {core.APP_VERSION}')
    data = parser.add_argument_group('данные')
    data.add_argument('--fasta', required=True, help='Референсный геном FASTA (рядом создаётся .fai)')
    data.add_argument('--variants', required=True, help='CSV/TSV/TXT с колонками chrom, pos, ref, alt')
    data.add_argument('--names', dest='names2_file', metavar='FILE', help='Файл соответствий названий хромосом')
    data.add_argument('--out', required=True, help='Файл результата: .csv или .xlsx')
    data.add_argument('--rejected', help='CSV с отклонёнными вариантами (в .xlsx они всегда на листе Rejected)')
    data.add_argument('--preset', choices=list(PRESETS), default='recommended',
                      help='Исходный набор параметров: recommended — настройки по итогам '
                           'проверки (как в окне программы), v1.1 — воспроизведение v1.1.0 '
                           '(по умолчанию: %(default)s)')

    search = parser.add_argument_group('поиск (по умолчанию — значения пресета)')
    search.add_argument('--flank', type=int, help='Фланг вокруг SNP, п.н.')
    search.add_argument('--min-frag', dest='min_frag', type=int, metavar='BP', help='Мин. видимый фрагмент, п.н.')
    search.add_argument('--max-frag', dest='max_frag', type=int, metavar='BP', help='Макс. видимый фрагмент, п.н.')
    search.add_argument('--delta', type=int, help='Диагностическая дельта, п.н.')
    search.add_argument('--max-cuts', dest='max_cuts', type=int, metavar='N', help='Макс. число разрезов в ампликоне')
    search.add_argument('--gain-loss-only', dest='gain_loss_only', action=argparse.BooleanOptionalAction,
                        help='Только появление/исчезновение сайта')
    search.add_argument('--suppliers', type=lambda s: [t.strip() for t in s.split(',') if t.strip()],
                        help='Поставщики через запятую: коды REBASE (N,E) или названия')

    choice = parser.add_argument_group('отбор и проверки')
    choice.add_argument('--commercial-only', dest='commercial_only', action=argparse.BooleanOptionalAction,
                        help='Только ферменты, у которых в REBASE есть поставщик')
    choice.add_argument('--group-isoschizomers', dest='group_isoschizomers',
                        action=argparse.BooleanOptionalAction,
                        help='Одна строка на сайт узнавания; остальные ферменты — в столбце isoschizomers')
    choice.add_argument('--snp-end-distance', dest='snp_end_distance', type=int, metavar='BP',
                        help='SNP не ближе BP п.н. к каждому концу ампликона (0 — не ограничивать)')

    primers = parser.add_argument_group('праймеры и совместный подбор')
    primers.add_argument('--primer3', dest='use_primer3', action=argparse.BooleanOptionalAction,
                         help='Подбирать праймеры Primer3 (--no-primer3 — поиск по окну)')
    primers.add_argument('--joint', dest='joint_design', action=argparse.BooleanOptionalAction,
                         help='Совместный подбор праймеров и фермента')
    primers.add_argument('--pairs', dest='primer_pairs', type=int, metavar='K', help='Пар Primer3 на SNP')
    primers.add_argument('--top', dest='top_results', type=int, metavar='N', help='Результатов на SNP (0 — все)')
    primers.add_argument('--prod-min', dest='prod_min', type=int, metavar='BP', help='Мин. размер продукта, п.н.')
    primers.add_argument('--prod-max', dest='prod_max', type=int, metavar='BP', help='Макс. размер продукта, п.н.')
    primers.add_argument('--tm-min', dest='tm_min', type=float, metavar='C', help='Мин. Tm, °C')
    primers.add_argument('--tm-max', dest='tm_max', type=float, metavar='C', help='Макс. Tm, °C')

    gel = parser.add_argument_group('модель геля')
    gel.add_argument('--resolution-bp', dest='gel_resolution_bp', type=float, metavar='BP', help='Разрешение, п.н.')
    gel.add_argument('--resolution-pct', dest='gel_resolution_pct', type=float, metavar='PCT', help='Разрешение, %%')
    gel.add_argument('--stress-factor', dest='gel_stress_factor', type=float, metavar='X', help='Ухудшение разрешения, ×')

    norm = parser.add_argument_group('нормализация')
    norm.add_argument('--skip-norm', dest='skip_norm', action=argparse.BooleanOptionalAction,
                      help='Не исправлять REF/ALT (проверка REF по FASTA остаётся)')
    norm.add_argument('--save-norm', dest='save_norm', action=argparse.BooleanOptionalAction,
                      help='Сохранять cleaned/rejected в results/normalize')
    parser.add_argument('-q', '--quiet', action='store_true', help='Не выводить журнал')
    return parser


def params_from_args(args):
    params = dict(PRESETS[args.preset])
    for key, value in vars(args).items():
        if key in params and value is not None:
            params[key] = value
    if args.save_norm is not None:
        params['no_save_norm'] = not args.save_norm
    params.update(fasta=str(Path(args.fasta).resolve()), input_file=str(Path(args.variants).resolve()))
    if params['names2_file']:
        params['names2_file'] = str(Path(params['names2_file']).resolve())
    return params


def main(argv=None):
    # A redirected Windows console uses a legacy code page without "×", "′" etc.
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, 'reconfigure'):
            stream.reconfigure(errors='replace')
    parser = build_parser()
    args = parser.parse_args(argv)
    out = Path(args.out)
    if out.suffix.lower() not in ('.csv', '.xlsx'):
        parser.error('--out: требуется расширение .csv или .xlsx')
    params = params_from_args(args)
    errors = core.validate_params(params)
    if errors:
        for error in errors:
            print(f'Ошибка параметров: {error}', file=sys.stderr)
        return 2

    def log(message):
        if not args.quiet:
            print(message, file=sys.stderr, flush=True)

    report = {}
    started = datetime.now()
    try:
        header, rows = core.run_rflp_gui_mode(params, status_cb=log, report=report)
    except KeyboardInterrupt:
        print('Анализ прерван.', file=sys.stderr)
        return 130
    summary = report.get('summary', {})
    if out.suffix.lower() == '.xlsx':
        export_excel(out, header, rows, metadata={
            'app_version': core.APP_VERSION, 'scoring_model': core.MODEL_VERSION,
            'run_timestamp': started.isoformat(timespec='seconds'), 'genome_path': params['fasta'],
            'parameters': report.get('parameters', params), 'status_summary': summary,
        }, rejected=report.get('rejected', []))
    else:
        export_csv(out, header, rows)
    if args.rejected:
        rejected = report.get('rejected', [])
        names = []
        for item in rejected:
            names += [key for key in item if key not in names]
        names = names or ['line', 'chrom', 'pos', 'ref', 'alt', 'reason', 'details']
        export_csv(args.rejected, names, [[item.get(name) for name in names] for item in rejected])
    log(f"[INFO] Сохранено: {out} · строк {len(rows)} · SNP принято {summary.get('accepted', 0)}, "
        f"с кандидатами {summary.get('variants_with_candidates', 0)}, "
        f"с robust-кандидатами {summary.get('variants_with_robust_candidates', 0)}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
