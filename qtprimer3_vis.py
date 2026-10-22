"""RFLP Picker main window: launch panel, SNP list and details in three columns (spec §4)."""
import os
import sys
import traceback
from datetime import datetime
from typing import List, Optional

from PyQt6.QtCore import QByteArray, QLocale, Qt, QThread, pyqtSignal
from PyQt6.QtGui import QIcon
from PyQt6.QtWidgets import (
    QApplication, QDialog, QHBoxLayout, QLabel, QMainWindow, QMenu, QMessageBox, QPlainTextEdit, QPushButton,
    QSplitter, QStatusBar, QToolButton, QVBoxLayout, QWidget,
)

import i18n
import resources_rc  # noqa: F401  (window icon)
import reporting
import rflp_core as core
from app_settings import AppSettings
from details_panel import DetailsPanel
from i18n import core_text, t
from launch_panel import LaunchPanel
from snp_list import SnpListView
from snp_summary import (
    ERROR, STATUS_ORDER, build_summaries, excel_sheets, read_input_records, results_table, status_title,
)
from theme import MODE_TITLES, MODES, ThemeController

LAUNCH_WIDTH = 300
LANGUAGE_MODES = ('system', 'ru', 'en')
# The names of the languages are always written in their own language.
LANGUAGE_NAMES = {'ru': 'Русский', 'en': 'English'}


def detected_language() -> str:
    return i18n.system_language(QLocale.system().uiLanguages())


def startup_language(settings: AppSettings) -> str:
    """The saved language, or the system one for 'system' and unknown values (i18n spec §3.2)."""
    return i18n.resolve_language(settings.value('language', 'system'), detected_language())


class CalcThread(QThread):
    status = pyqtSignal(str)
    progress = pyqtSignal(int, int)
    succeeded = pyqtSignal(list, list, dict)
    failed = pyqtSignal(str)
    cancelled = pyqtSignal()

    def __init__(self, params: dict):
        super().__init__()
        self.params = params

    def run(self):
        report: dict = {}
        try:
            header, rows = core.run_rflp_gui_mode(self.params, status_cb=self.status.emit,
                                                  progress_cb=self.progress.emit,
                                                  cancel_cb=self.isInterruptionRequested, report=report)
        except core.AnalysisCancelled:
            self.cancelled.emit()
        except Exception as exc:
            self.status.emit(traceback.format_exc())
            core.logger.error(traceback.format_exc())
            self.failed.emit(str(exc) or type(exc).__name__)
        else:
            self.succeeded.emit(header, rows, report)


class LogWindow(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(t('Журнал расчёта'))
        self.resize(760, 480)
        layout = QVBoxLayout(self)
        self.text = QPlainTextEdit()
        self.text.setReadOnly(True)
        layout.addWidget(self.text, 1)
        row = QHBoxLayout()
        row.addStretch()
        copy = QPushButton(t('Копировать'))
        copy.clicked.connect(lambda: QApplication.clipboard().setText(self.text.toPlainText()))
        row.addWidget(copy)
        layout.addLayout(row)

    def append(self, text: str) -> None:
        self.text.appendPlainText(text)


class RFLPPickerGUI(QMainWindow):
    def __init__(self, settings: Optional[AppSettings] = None, theme: Optional[ThemeController] = None):
        super().__init__()
        self.settings = settings or AppSettings()
        self.theme = theme or ThemeController(QApplication.instance(), self._stored_theme_mode())
        self.worker: Optional[CalcThread] = None
        self.header: List[str] = []
        self.rows: List[list] = []
        self.report: dict = {}
        self.params: dict = {}
        self.summaries = []
        self._has_ids = False
        self.setWindowTitle(f'RFLP Picker {core.APP_VERSION}')
        self.setWindowIcon(QIcon(':/icons/app_icon.ico'))
        self.setMinimumSize(1000, 640)
        self.resize(1600, 900)

        self.launch = LaunchPanel(self.settings)
        self.list = SnpListView(self.settings)
        self.details = DetailsPanel()
        self.splitter = QSplitter(Qt.Orientation.Horizontal)
        self.splitter.setChildrenCollapsible(False)
        for widget, margins in ((self.launch, None), (self.list, (8, 12, 8, 8)), (self.details, (8, 12, 12, 8))):
            if margins:
                box = QWidget()
                layout = QVBoxLayout(box)
                layout.setContentsMargins(*margins)
                layout.addWidget(widget)
                self.splitter.addWidget(box)
            else:
                self.splitter.addWidget(widget)
        self.details.parentWidget().setMinimumWidth(360)
        self.splitter.setStretchFactor(1, 55)
        self.splitter.setStretchFactor(2, 45)
        self.setCentralWidget(self.splitter)

        self.log = LogWindow(self)
        status = QStatusBar()
        self.status_message = QLabel()
        status.addWidget(self.status_message, 1)
        log_button = QPushButton(t('Журнал'))
        log_button.clicked.connect(self.log.show)
        status.addPermanentWidget(log_button)
        self.theme_button = QToolButton()
        self.theme_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        theme_menu = QMenu(self.theme_button)
        for mode in MODES:
            theme_menu.addAction(t(MODE_TITLES[mode]), lambda m=mode: self.set_theme_mode(m))
        self.theme_button.setMenu(theme_menu)
        status.addPermanentWidget(self.theme_button)
        self.language_button = QToolButton()
        self.language_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        language_menu = QMenu(self.language_button)
        for mode in LANGUAGE_MODES:
            title = t('Системный') if mode == 'system' else LANGUAGE_NAMES[mode]
            language_menu.addAction(title, lambda m=mode: self.set_language_mode(m))
        self.language_button.setMenu(language_menu)
        status.addPermanentWidget(self.language_button)
        status.addPermanentWidget(QLabel(f'RFLP Picker {core.APP_VERSION}'))
        self.setStatusBar(status)

        self.launch.run_requested.connect(self.start_run)
        self.launch.stop_requested.connect(self.stop_run)
        self.launch.collapse_toggled.connect(self._on_collapse_toggled)
        self.list.selection_changed.connect(self._on_selection)
        self.list.export_requested.connect(self.export)
        self.details.choice_changed.connect(lambda _summary: self.list.refresh_current())
        self.theme.changed.connect(self._apply_colors)
        self._apply_colors(self.theme.colors)
        self._update_theme_button()
        self._update_language_button()
        self._restore_layout()

    # ----- theme -----
    def _stored_theme_mode(self) -> str:
        mode = self.settings.value('theme_mode', 'system')
        return mode if mode in MODES else 'system'

    def set_theme_mode(self, mode: str) -> None:
        self.theme.set_mode(mode)
        self._update_theme_button()

    def _update_theme_button(self):
        self.theme_button.setText(t('Тема: {mode} ▾', mode=t(MODE_TITLES[self.theme.mode]).lower()))

    # ----- language -----
    def _language_mode(self) -> str:
        mode = self.settings.value('language', 'system')
        return mode if mode in LANGUAGE_MODES else 'system'

    def set_language_mode(self, mode: str) -> None:
        """Save the language for the next launch; the window keeps its language until then."""
        self.settings.set_value('language', mode)
        self._update_language_button()
        if startup_language(self.settings) != i18n.language():
            QMessageBox.information(self, t('Язык'), t('Язык изменится после перезапуска программы.'))

    def _update_language_button(self):
        mode = self._language_mode()
        name = t('системный') if mode == 'system' else LANGUAGE_NAMES[mode]
        self.language_button.setText(t('Язык: {name} ▾', name=name))

    def _apply_colors(self, colors: dict):
        for widget in (self.launch, self.list, self.details):
            widget.apply_colors(colors)

    # ----- layout -----
    def _restore_layout(self):
        geometry = self.settings.value('geometry')
        if isinstance(geometry, str):
            self.restoreGeometry(QByteArray.fromBase64(geometry.encode()))
        # A new session has no results yet, so «Запуск» always starts open.
        self.launch.set_collapsed(False)
        sizes = self.settings.value('splitter_sizes')
        if isinstance(sizes, list) and len(sizes) == 3 and all(isinstance(x, int) for x in sizes):
            rest = max(0, sum(sizes) - max(sizes[0], LAUNCH_WIDTH))
            ratio = sizes[1] / max(1, sizes[1] + sizes[2])
            self.splitter.setSizes([max(sizes[0], LAUNCH_WIDTH), round(rest * ratio), rest - round(rest * ratio)])
        else:
            rest = max(600, self.width() - LAUNCH_WIDTH)
            self.splitter.setSizes([LAUNCH_WIDTH, rest * 55 // 100, rest * 45 // 100])

    def _on_collapse_toggled(self, collapsed: bool):
        sizes = self.splitter.sizes()
        launch = 36 if collapsed else LAUNCH_WIDTH
        rest = max(0, sum(sizes) - launch)
        self.splitter.setSizes([launch, rest * 55 // 100, rest - rest * 55 // 100])

    def closeEvent(self, event):
        self.launch.save_settings()
        self.settings.set_value('theme_mode', self.theme.mode)
        self.settings.set_value('splitter_sizes', self.splitter.sizes())
        self.settings.set_value('geometry', bytes(self.saveGeometry().toBase64()).decode())
        if self.worker is not None:
            self.worker.requestInterruption()
            self.worker.wait(5000)
        self.log.close()
        super().closeEvent(event)

    # ----- running -----
    def manual_choice_count(self) -> int:
        return sum(s.manual for s in self.list.model.summaries())

    def confirm_rerun(self, manual: int) -> bool:
        answer = QMessageBox.question(self, t('Новый расчёт'),
                                      t('Ручной выбор у {count} SNP будет сброшен. Запустить?', count=manual))
        return answer == QMessageBox.StandardButton.Yes

    def start_run(self) -> None:
        if self.worker is not None or not self.launch.validate():
            return
        manual = self.manual_choice_count()
        if manual and not self.confirm_rerun(manual):
            return
        params = self.launch.collect_params()
        self.launch.save_settings()
        self.log.text.clear()
        self.status_message.setText(t('Идёт расчёт…'))
        self.launch.set_running(True)
        worker = CalcThread(params)
        worker.status.connect(self._on_status)
        worker.progress.connect(self._on_progress)
        worker.succeeded.connect(lambda header, rows, report: self.on_succeeded(header, rows, report, params))
        worker.failed.connect(self._on_failed)
        worker.cancelled.connect(self._on_cancelled)
        worker.finished.connect(self._on_finished)
        self.worker = worker
        worker.start()

    def stop_run(self) -> None:
        if self.worker is not None:
            self.worker.requestInterruption()
            self.status_message.setText(t('Останавливаю расчёт…'))

    def _on_status(self, message: str):
        message = core_text(message)
        self.log.append(message)
        last = next((line.strip() for line in reversed(message.splitlines()) if line.strip()), '')
        if last.startswith('[INFO]'):
            self.status_message.setText(last[len('[INFO]'):].strip())

    def _on_progress(self, done: int, total: int):
        self.launch.set_progress(done, total)
        if self.list.model.rowCount() == 0:
            self.list.set_placeholder(t('Идёт расчёт… {done} из {total} SNP', done=done, total=total))

    def _on_finished(self):
        self.worker = None
        self.launch.set_running(False)

    def _on_failed(self, message: str):
        self.status_message.setText(t('Ошибка расчёта: {message}', message=core_text(message)))
        self._restore_list()

    def _on_cancelled(self):
        self.status_message.setText(t('Расчёт остановлен'))
        self._restore_list()

    def _restore_list(self):
        if self.list.model.rowCount():
            self.list.set_summaries(self.list.model.summaries(), self.params.get('dcaps', 'off'), self._has_ids)
        else:
            self.list.set_placeholder(t('Выберите файл SNP и геном слева и нажмите «Запустить»'))

    def on_succeeded(self, header: list, rows: list, report: dict, params: Optional[dict] = None) -> None:
        self.header, self.rows, self.report = list(header), list(rows), dict(report)
        self.params = dict(params or self.launch.collect_params())
        records, self._has_ids = read_input_records(self.params['input_file'])
        self.summaries = build_summaries(records, self.header, self.rows, self.report.get('rejected', []))
        self.details.set_context(self.params['fasta'], self.params)
        self.details.show_summary(None, self.params['dcaps'])
        self.list.set_summaries(self.summaries, self.params['dcaps'], self._has_ids)
        if not self.summaries:
            self.list.set_placeholder(t('В файле нет SNP'))
            self.status_message.setText(t('В файле нет SNP'))
            return
        counts = self.list.model.counts()
        parts = [f'{status_title(s).lower()} {counts[s]}' for s in STATUS_ORDER if counts[s]]
        self.status_message.setText(t('Готово: {count} SNP · {parts}', count=len(self.summaries),
                                      parts=', '.join(parts)))
        if counts[ERROR] < len(self.summaries):
            # When every SNP is a data error, the fix (file, genome, names) is in «Запуск»: keep it open.
            self.launch.set_collapsed(True)
            self._on_collapse_toggled(True)

    # ----- export -----
    def _default_path(self, extension: str) -> str:
        folder = os.path.join(core.get_base_dir(), 'results')
        os.makedirs(folder, exist_ok=True)
        return os.path.join(folder, f"rflp_candidates_{datetime.now().strftime('%Y%m%d%H%M%S')}.{extension}")

    def export(self, kind: str, path: Optional[str] = None) -> Optional[str]:
        """Write 'all' or 'shown' SNPs to Excel, or 'csv' (all schemes). Return the path, or None."""
        if not self.summaries:
            QMessageBox.information(self, t('Выгрузка'), t('Сначала выполните расчёт.'))
            return None
        interactive = path is None
        folder = os.path.dirname(path) if path else os.path.join(core.get_base_dir(), 'results')
        try:
            path = path or self._default_path('csv' if kind == 'csv' else 'xlsx')
            summaries = self.list.shown_summaries() if kind == 'shown' else self.list.model.summaries()
            table_header, table_rows = results_table(summaries, self.header, self.rows)
            if kind == 'csv':
                reporting.export_csv(path, table_header, table_rows)
            else:
                reporting.export_excel(path, table_header, table_rows, metadata={
                    'app_version': core.APP_VERSION, 'scoring_model': core.MODEL_VERSION,
                    'run_timestamp': datetime.now().isoformat(timespec='seconds'),
                    'genome_path': self.params.get('fasta', ''), 'parameters': self.params,
                    'status_summary': self.report.get('summary', {}),
                }, rejected=self.report.get('rejected', []),
                    extra_sheets=excel_sheets(summaries, self.params.get('dcaps', 'off')))
        except PermissionError:
            if path and os.path.exists(path):
                message = t('Файл занят другой программой. Закройте его и повторите выгрузку.')
            else:
                message = t('Нет доступа к папке {folder}. Запустите программу из папки, куда разрешена запись.',
                            folder=folder)
            QMessageBox.warning(self, t('Выгрузка'), message)
            return None
        except Exception as exc:
            # Never let an export error close the window with the user's choices.
            self.log.append(traceback.format_exc())
            core.logger.error(traceback.format_exc())
            QMessageBox.warning(self, t('Выгрузка'), t('Не удалось сохранить: {error}', error=exc))
            return None
        self.status_message.setText(t('Сохранено: {path}', path=path))
        if interactive:
            QMessageBox.information(self, t('Выгрузка'), t('Результаты сохранены:\n{path}', path=path))
        return path

    # ----- selection -----
    def _on_selection(self, summary):
        self.details.show_summary(summary, self.params.get('dcaps', 'off'))


def install_exception_hook(window: RFLPPickerGUI) -> None:
    """Log unexpected errors instead of letting PyQt abort the process (and lose the results)."""
    def hook(exc_type, exc, tb):
        text = ''.join(traceback.format_exception(exc_type, exc, tb))
        core.logger.error(text)
        window.log.append(text)
        window.status_message.setText(t('Ошибка: {error}. Подробности — в журнале.', error=exc))
        QMessageBox.critical(window, t('Ошибка'), t('{error}\n\nПодробности — в окне «Журнал».', error=exc))

    sys.excepthook = hook


def main() -> None:
    if not core._PRIMER3_AVAILABLE:
        core.logger.warning(f'primer3 недоступен, ошибка импорта: {core._PRIMER3_IMPORT_ERROR}')
    app = QApplication(sys.argv)
    settings = AppSettings()
    i18n.set_language(startup_language(settings))   # before any window is built
    window = RFLPPickerGUI(settings)
    install_exception_hook(window)
    window.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
