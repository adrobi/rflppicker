"""Screenshots of the RFLP Picker window for the manual check and the README.

Runs one analysis in the chosen interface language, selects the first SNP
with a robust scheme and saves DIR/<lang>_<theme>_<width>x<height>.png for both
themes and three window sizes (logical pixels), plus
DIR/<lang>_start_light_1536x816.png before the run. Uses its own settings file,
so the user's settings stay untouched.

    python tools/gui_screenshots.py --lang en --fasta genome.fa --variants snps.tsv --out shots
"""
import argparse
import os
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from PyQt6.QtCore import QSettings, QTimer  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

import i18n  # noqa: E402
from app_settings import AppSettings  # noqa: E402
from qtprimer3_vis import RFLPPickerGUI, apply_language  # noqa: E402
from snp_summary import ROBUST  # noqa: E402
from theme import ThemeController  # noqa: E402

SIZES = ((1920, 1040), (1536, 816), (1280, 680))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--variants', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--lang', choices=i18n.LANGS, default='ru', help='Interface language (default: ru)')
    args = parser.parse_args(argv)
    apply_language(args.lang)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    app = QApplication(sys.argv)
    workdir = tempfile.TemporaryDirectory()
    settings = AppSettings(QSettings(os.path.join(workdir.name, 'shots.ini'), QSettings.Format.IniFormat))
    window = RFLPPickerGUI(settings, ThemeController(app, 'light'))
    window.launch.set_files(str(Path(args.variants).resolve()), str(Path(args.fasta).resolve()))
    window.resize(*SIZES[1])
    window.show()
    saved = []
    result = {'code': 0}

    def shot(name):
        app.processEvents()
        path = out / name
        window.grab().save(str(path))
        saved.append(path)

    def after_run():
        if not window.summaries:
            print(window.status_message.text(), file=sys.stderr)
            result['code'] = 1
            app.quit()
            return
        statuses = [s.status for s in window.list.shown_summaries()]
        window.list.table.selectRow(statuses.index(ROBUST) if ROBUST in statuses else 0)
        for mode in ('light', 'dark'):
            window.set_theme_mode(mode)
            for width, height in SIZES:
                window.resize(width, height)
                shot(f'{args.lang}_{mode}_{width}x{height}.png')
        app.quit()

    def start():
        shot(f'{args.lang}_start_light_1536x816.png')
        window.start_run()
        if window.worker is None:
            print('The analysis did not start: check the files.', file=sys.stderr)
            result['code'] = 2
            app.quit()
            return
        window.worker.finished.connect(lambda: QTimer.singleShot(500, after_run))

    QTimer.singleShot(500, start)
    app.exec()
    window.close()
    workdir.cleanup()
    for path in saved:
        print(path)
    return result['code']


if __name__ == '__main__':
    sys.exit(main())
