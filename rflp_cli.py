"""RFLP Picker from the command line: the same analysis as in the window.

The help, the messages and the Summary/Primers sheets follow --lang or, without
it, the language of the system (i18n spec §3.2, §4.4).
"""
import argparse
import locale
import logging
import os
import sys
from datetime import datetime
from pathlib import Path

from pyfaidx import FastaIndexingError

import rflp_core as core
from i18n import LANGS, N_, core_text, locale_language, set_language, t
from reporting import export_csv, export_excel
from snp_summary import build_summaries, excel_sheets, read_input_records, results_table

DESCRIPTION = N_('''RFLP Picker из командной строки: тот же анализ, что в окне программы.

Пример:
    python rflp_cli.py --fasta genome.fa --variants snps.tsv --out results.xlsx

Результат сохраняется в CSV или Excel (по расширению --out). В Excel, как и
при экспорте из окна, добавляются листы Run (параметры, версии, сводка) и
Rejected (отклонённые варианты). Журнал выводится в stderr.

Скачать референсный геном с NCBI в папку genoms рядом с программой:
    python rflp_cli.py --download-genome GCF_016772045.1''')


def build_parser():
    parser = argparse.ArgumentParser(
        prog='rflp_cli', description=t(DESCRIPTION),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version=f'RFLP Picker {core.APP_VERSION}')
    parser.add_argument('--lang', choices=LANGS,
                        help=t('Язык справки, сообщений и листов Сводка/Праймеры (по умолчанию — язык системы)'))
    data = parser.add_argument_group(t('данные'))
    data.add_argument('--fasta', required=True, help=t('Референсный геном FASTA (рядом создаётся .fai)'))
    data.add_argument('--variants', required=True, help=t('CSV/TSV/TXT с колонками chrom, pos, ref, alt'))
    data.add_argument('--names', dest='names2_file', metavar='FILE',
                      help=t('Файл соответствий названий хромосом'))
    data.add_argument('--out', required=True, help=t('Файл результата: .csv или .xlsx'))
    data.add_argument('--rejected', help=t('CSV с отклонёнными вариантами (в .xlsx они всегда на листе Rejected)'))

    search = parser.add_argument_group(t('поиск (по умолчанию — как в окне программы)'))
    search.add_argument('--flank', type=int, help=t('Фланг вокруг SNP, п.н.'))
    search.add_argument('--min-frag', dest='min_frag', type=int, metavar='BP',
                        help=t('Мин. видимый фрагмент, п.н.'))
    search.add_argument('--max-frag', dest='max_frag', type=int, metavar='BP',
                        help=t('Макс. видимый фрагмент, п.н.'))
    search.add_argument('--delta', type=int, help=t('Диагностическая дельта, п.н.'))
    search.add_argument('--max-cuts', dest='max_cuts', type=int, metavar='N',
                        help=t('Макс. число разрезов в ампликоне'))
    search.add_argument('--gain-loss-only', dest='gain_loss_only', action=argparse.BooleanOptionalAction,
                        help=t('Только появление/исчезновение сайта'))
    search.add_argument('--suppliers', type=lambda s: [part.strip() for part in s.split(',') if part.strip()],
                        help=t('Поставщики через запятую: коды REBASE (N,E) или названия'))

    choice = parser.add_argument_group(t('отбор и проверки'))
    choice.add_argument('--commercial-only', dest='commercial_only', action=argparse.BooleanOptionalAction,
                        help=t('Только ферменты, у которых в REBASE есть поставщик'))
    choice.add_argument('--group-isoschizomers', dest='group_isoschizomers',
                        action=argparse.BooleanOptionalAction,
                        help=t('Одна строка на сайт узнавания; остальные ферменты — в столбце isoschizomers'))
    choice.add_argument('--snp-end-distance', dest='snp_end_distance', type=int, metavar='BP',
                        help=t('SNP не ближе BP п.н. к каждому концу ампликона (0 — не ограничивать)'))
    choice.add_argument('--exclude-primer-sites', dest='exclude_primer_sites',
                        action=argparse.BooleanOptionalAction,
                        help=t('Исключать схемы, где праймер сам содержит сайт фермента'))
    choice.add_argument('--mask-repeats', dest='mask_repeats', action=argparse.BooleanOptionalAction,
                        help=t('Повторы (строчные буквы soft-mask FASTA) не на 3′-концах праймеров'))
    choice.add_argument('--background-vcf', dest='background_vcf', metavar='FILE',
                        help=t('VCF/VCF.gz известных вариантов: праймеры их избегают, влияние на разрезы '
                               'отмечается (названия хромосом — через --names)'))
    choice.add_argument('--exclude-background-sites', dest='exclude_background_sites',
                        action=argparse.BooleanOptionalAction,
                        help=t('Исключать схемы, где фоновые варианты меняют разрезы или длины полос'))
    choice.add_argument('--dcaps', choices=('off', 'fallback', 'always'),
                        help=t('Праймеры с несовпадением (dCAPS): fallback — только если нет robust-схемы'))
    choice.add_argument('--dcaps-mismatches', dest='dcaps_mismatches', type=int, choices=(1, 2),
                        help=t('Наибольшее число несовпадений в праймере dCAPS'))
    choice.add_argument('--dcaps-tm-max', dest='dcaps_tm_max', type=float, metavar='C',
                        help=t('dCAPS: верхняя граница Tm обоих праймеров, °C (нижняя — --tm-min)'))
    choice.add_argument('--dcaps-prod-min', dest='dcaps_prod_min', type=int, metavar='BP',
                        help=t('dCAPS: мин. размер продукта, п.н.'))
    choice.add_argument('--dcaps-prod-max', dest='dcaps_prod_max', type=int, metavar='BP',
                        help=t('dCAPS: макс. размер продукта, п.н.'))
    choice.add_argument('--dcaps-delta', dest='dcaps_delta', type=int, metavar='BP',
                        help=t('dCAPS: порог Δ, п.н. (вместо --delta для схем dCAPS)'))

    primers = parser.add_argument_group(t('праймеры'))
    primers.add_argument('--primer3', dest='use_primer3', action=argparse.BooleanOptionalAction,
                         help=t('Подбирать праймеры Primer3 (--no-primer3 — поиск по окну)'))
    primers.add_argument('--pairs', dest='primer_pairs', type=int, metavar='K', help=t('Пар Primer3 на SNP'))
    primers.add_argument('--top', dest='top_results', type=int, metavar='N', help=t('Результатов на SNP (0 — все)'))
    primers.add_argument('--prod-min', dest='prod_min', type=int, metavar='BP',
                         help=t('Мин. размер продукта, п.н.'))
    primers.add_argument('--prod-max', dest='prod_max', type=int, metavar='BP',
                         help=t('Макс. размер продукта, п.н.'))
    primers.add_argument('--tm-min', dest='tm_min', type=float, metavar='C', help=t('Мин. Tm, °C'))
    primers.add_argument('--tm-max', dest='tm_max', type=float, metavar='C', help=t('Макс. Tm, °C'))

    gel = parser.add_argument_group(t('модель геля'))
    gel.add_argument('--resolution-bp', dest='gel_resolution_bp', type=float, metavar='BP',
                     help=t('Разрешение, п.н.'))
    gel.add_argument('--resolution-pct', dest='gel_resolution_pct', type=float, metavar='PCT',
                     help=t('Разрешение, %%'))
    gel.add_argument('--stress-factor', dest='gel_stress_factor', type=float, metavar='X',
                     help=t('Ухудшение разрешения, ×'))

    norm = parser.add_argument_group(t('нормализация'))
    norm.add_argument('--skip-norm', dest='skip_norm', action=argparse.BooleanOptionalAction,
                      help=t('Не исправлять REF/ALT (проверка REF по FASTA остаётся)'))
    norm.add_argument('--save-norm', dest='save_norm', action=argparse.BooleanOptionalAction,
                      help=t('Сохранять cleaned/rejected в results/normalize'))
    parser.add_argument('-q', '--quiet', action='store_true', help=t('Не выводить журнал'))
    return parser


def params_from_args(args):
    params = dict(core.DEFAULT_PARAMS)
    for key, value in vars(args).items():
        if key in params and value is not None:
            params[key] = value
    if args.save_norm is not None:
        params['no_save_norm'] = not args.save_norm
    params.update(fasta=str(Path(args.fasta).resolve()), input_file=str(Path(args.variants).resolve()))
    for key in ('names2_file', 'background_vcf'):
        if params[key]:
            params[key] = str(Path(params[key]).resolve())
    return params


def choose_language(argv) -> None:
    """--lang if given, otherwise the language of the system; before the parser is built."""
    early = argparse.ArgumentParser(add_help=False)
    early.add_argument('--lang', nargs='?')   # a missing value is reported by the full parser
    known, _ = early.parse_known_args(argv)
    if known.lang in LANGS:
        set_language(known.lang)
        return
    try:
        name = locale.getlocale()[0]
    except ValueError:   # a locale name Python does not know
        name = None
    set_language(locale_language(name))


class CoreLogTranslation(logging.Filter):
    """The core logs in Russian; on the command line its messages follow the chosen language."""

    def filter(self, record):
        record.msg, record.args = core_text(record.getMessage()), None
        return True


def translate_core_log() -> None:
    if not any(isinstance(item, CoreLogTranslation) for item in core.logger.filters):
        core.logger.addFilter(CoreLogTranslation())


def download_command(argv) -> int:
    """--download-genome ACCESSION: the genome and its chromosome names from NCBI; the path goes to stdout."""
    import genome_download
    parser = argparse.ArgumentParser(prog='rflp_cli', add_help=False)
    parser.add_argument('--download-genome', required=True, metavar='ACCESSION')
    parser.add_argument('--lang', choices=LANGS)
    parser.add_argument('-q', '--quiet', action='store_true')
    args = parser.parse_args(argv)
    shown = [-1]

    def progress(done, total):
        tenth = 10 * done // total if total else 0
        if not args.quiet and tenth != shown[0]:
            shown[0] = tenth
            print(t('Загрузка с NCBI: {done} из {total} МБ', done=round(done / 1e6), total=round(total / 1e6)),
                  file=sys.stderr, flush=True)

    def status(text):
        if not args.quiet:
            print(text, file=sys.stderr, flush=True)

    try:
        fasta, names = genome_download.download_genome(args.download_genome, core.get_base_dir(),
                                                       progress_cb=progress, status_cb=status)
    except KeyboardInterrupt:
        print(t('Загрузка остановлена; при повторе она продолжится с этого места.'), file=sys.stderr)
        return 130
    except (genome_download.GenomeDownloadError, OSError) as exc:
        print(t('Ошибка: {error}', error=str(exc) or type(exc).__name__), file=sys.stderr)
        return 1
    status(t('Геном: {fasta}\nХромосомы: {names}', fasta=fasta, names=names))
    print(fasta)
    return 0


def main(argv=None):
    # A redirected Windows console uses a legacy code page without "×", "′" etc.
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, 'reconfigure'):
            stream.reconfigure(errors='replace')
    choose_language(argv)
    translate_core_log()
    if any(str(arg).split('=')[0] == '--download-genome' for arg in (sys.argv[1:] if argv is None else argv)):
        return download_command(argv)
    parser = build_parser()
    args = parser.parse_args(argv)
    out = Path(args.out)
    if out.suffix.lower() not in ('.csv', '.xlsx'):
        parser.error(t('--out: требуется расширение .csv или .xlsx'))
    # Before the analysis: a mistyped folder must not cost a finished run.
    for key, target in (('--out', args.out), ('--rejected', args.rejected)):
        if not target:
            continue
        target = Path(target).expanduser().absolute()
        if not target.parent.is_dir():
            parser.error(t('{key}: папка не существует: {folder}', key=key, folder=target.parent))
        if target.is_dir():
            parser.error(t('{key}: это папка, а не файл: {path}', key=key, path=target))
        if not os.access(target.parent, os.W_OK):
            parser.error(t('{key}: нет права записи в папку: {folder}', key=key, folder=target.parent))
    params = params_from_args(args)
    errors = core.validate_params(params)
    if errors:
        for error in errors:
            print(t('Ошибка параметров: {error}', error=core_text(error)), file=sys.stderr)
        return 2

    def log(message):
        if not args.quiet:
            print(core_text(message), file=sys.stderr, flush=True)

    report = {}
    started = datetime.now()
    try:
        header, rows = core.run_rflp_gui_mode(params, status_cb=log, report=report)
    except KeyboardInterrupt:
        print(t('Анализ прерван.'), file=sys.stderr)
        return 130
    except (OSError, ValueError, FastaIndexingError) as exc:   # an unreadable FASTA or VCF, a lost file, …
        print(t('Ошибка: {error}', error=core_text(str(exc) or type(exc).__name__)), file=sys.stderr)
        return 1
    try:
        save_results(args, out, params, report, header, rows, started)
    except (OSError, ValueError) as exc:   # the destination is open in Excel, over the Excel limits, …
        print(t('Ошибка: {error}', error=str(exc) or type(exc).__name__), file=sys.stderr)
        return 1
    summary = report.get('summary', {})
    if not summary.get('accepted') and not report.get('rejected'):
        print(t('[WARN] В файле нет SNP.'), file=sys.stderr)
    elif not summary.get('accepted'):
        # Shown even with --quiet: the files are written, but they hold no assay.
        print(t('[WARN] Ни один SNP не принят: проверьте геном и файл соответствий названий хромосом.'),
              file=sys.stderr)
    log(t('[INFO] Сохранено: {out} · строк {rows} · SNP принято {accepted}, с кандидатами {candidates}, '
          'с robust-кандидатами {robust}', out=out, rows=len(rows), accepted=summary.get('accepted', 0),
          candidates=summary.get('variants_with_candidates', 0),
          robust=summary.get('variants_with_robust_candidates', 0)))
    return 0


def save_results(args, out, params, report, header, rows, started) -> None:
    summary = report.get('summary', {})
    # One summary per line of the variants file; the best scheme is the chosen one.
    records, _ = read_input_records(params['input_file'])
    summaries = build_summaries(records, header, rows, report.get('rejected', []))
    table_header, table_rows = results_table(summaries, header, rows)
    if out.suffix.lower() == '.xlsx':
        export_excel(out, table_header, table_rows, metadata={
            'app_version': core.APP_VERSION, 'scoring_model': core.MODEL_VERSION,
            'run_timestamp': started.isoformat(timespec='seconds'), 'genome_path': params['fasta'],
            'parameters': report.get('parameters', params), 'status_summary': summary,
        }, rejected=report.get('rejected', []), extra_sheets=excel_sheets(summaries, params['dcaps']))
    else:
        export_csv(out, table_header, table_rows)
    if args.rejected:
        rejected = report.get('rejected', [])
        names = []
        for item in rejected:
            names += [key for key in item if key not in names]
        names = names or ['line', 'chrom', 'pos', 'ref', 'alt', 'reason', 'details']
        export_csv(args.rejected, names, [[item.get(name) for name in names] for item in rejected])


if __name__ == '__main__':
    sys.exit(main())
