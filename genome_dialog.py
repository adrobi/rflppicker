"""The «Скачать геном…» dialog: a reference genome from NCBI into genoms/ (see genome_download.py)."""
from typing import Optional

from PyQt6.QtCore import QThread, pyqtSignal
from PyQt6.QtWidgets import (
    QComboBox, QDialog, QHBoxLayout, QLabel, QLineEdit, QProgressBar, QPushButton, QVBoxLayout,
)

import genome_download
from i18n import fmt, t


class DownloadThread(QThread):
    progress = pyqtSignal(float, float)   # bytes: more than a 32-bit int holds
    status = pyqtSignal(str)
    succeeded = pyqtSignal(str, str)
    failed = pyqtSignal(str)
    cancelled = pyqtSignal()

    def __init__(self, download, accession: str, base_dir: str):
        super().__init__()
        self._download, self._accession, self._base_dir = download, accession, base_dir

    def run(self):
        try:
            fasta, names = self._download(self._accession, self._base_dir,
                                          progress_cb=lambda done, total: self.progress.emit(done, total),
                                          status_cb=self.status.emit, cancel_cb=self.isInterruptionRequested)
        except genome_download.DownloadCancelled:
            self.cancelled.emit()
        except Exception as error:   # shown in the dialog; the program goes on
            self.failed.emit(str(error) or type(error).__name__)
        else:
            self.succeeded.emit(fasta, names)


class GenomeDownloadDialog(QDialog):
    def __init__(self, base_dir: str, parent=None, download=genome_download.download_genome):
        super().__init__(parent)
        self.setWindowTitle(t('Скачать геном'))
        self.setMinimumWidth(520)
        self._base_dir, self._download = base_dir, download
        self._thread: Optional[DownloadThread] = None
        self.fasta_path = ''
        layout = QVBoxLayout(self)

        layout.addWidget(QLabel(t('Сборка генома на NCBI')))
        self.assembly_combo = QComboBox()
        for title, _species, accession, name, size in genome_download.PRESETS:
            self.assembly_combo.addItem(t('{species} — {name} ({accession}, {size} ГБ)', species=t(title), name=name,
                                          accession=accession, size=fmt(size)), accession)
        self.assembly_combo.addItem(t('Другая сборка…'), None)
        self.assembly_combo.currentIndexChanged.connect(self._on_assembly_changed)
        layout.addWidget(self.assembly_combo)
        self.accession_edit = QLineEdit()
        self.accession_edit.setPlaceholderText(t('Номер сборки, например GCF_016772045.1'))
        self.accession_edit.setEnabled(False)
        layout.addWidget(self.accession_edit)
        note = QLabel(t('Координаты в файле SNP должны относиться к той же сборке. Геном сохраняется в папку '
                        'genoms рядом с программой; файл соответствий названий хромосом создаётся автоматически.'))
        note.setWordWrap(True)
        note.setProperty('role', 'muted')
        layout.addWidget(note)

        self.progress = QProgressBar()
        self.progress.setRange(0, 1000)
        self.progress.setTextVisible(False)
        layout.addWidget(self.progress)
        self.status_label = QLabel()
        layout.addWidget(self.status_label)
        self.error_label = QLabel()
        self.error_label.setWordWrap(True)
        self.error_label.setProperty('role', 'error')
        self.error_label.hide()
        layout.addWidget(self.error_label)

        buttons = QHBoxLayout()
        buttons.addStretch()
        self.cancel_button = QPushButton(t('Закрыть'))
        self.cancel_button.clicked.connect(self.reject)
        buttons.addWidget(self.cancel_button)
        self.download_button = QPushButton(t('Скачать'))
        self.download_button.setProperty('role', 'primary')
        self.download_button.clicked.connect(self.start)
        buttons.addWidget(self.download_button)
        layout.addLayout(buttons)

    def _on_assembly_changed(self, index: int):
        self.accession_edit.setEnabled(self.assembly_combo.itemData(index) is None)

    def accession(self) -> str:
        return self.assembly_combo.currentData() or self.accession_edit.text().strip()

    def _set_running(self, running: bool):
        self.download_button.setEnabled(not running)
        self.assembly_combo.setEnabled(not running)
        self.accession_edit.setEnabled(not running and self.assembly_combo.currentData() is None)
        self.cancel_button.setText(t('Отмена') if running else t('Закрыть'))

    def start(self) -> None:
        if self._thread is not None:
            return
        self.error_label.hide()
        self.progress.setValue(0)
        self._set_running(True)
        thread = DownloadThread(self._download, self.accession(), self._base_dir)
        thread.progress.connect(self._on_progress)
        thread.status.connect(self.status_label.setText)
        thread.succeeded.connect(self._on_succeeded)
        thread.failed.connect(self._on_failed)
        thread.cancelled.connect(self._on_cancelled)
        thread.finished.connect(self._on_finished)
        self._thread = thread
        thread.start()

    def cancel(self) -> None:
        if self._thread is not None:
            self._thread.requestInterruption()
            self.status_label.setText(t('Останавливаю загрузку…'))

    def wait_for_thread(self) -> None:
        if self._thread is not None:
            self._thread.wait()

    def _on_progress(self, done: float, total: float):
        self.progress.setValue(round(1000 * done / total) if total else 0)
        self.status_label.setText(t('Загрузка с NCBI: {done} из {total} МБ', done=round(done / 1e6),
                                    total=round(total / 1e6)))

    def _on_succeeded(self, fasta: str, _names: str):
        self.fasta_path = fasta
        self.accept()

    def _on_failed(self, message: str):
        self.status_label.setText('')
        self.error_label.setText(message)
        self.error_label.show()

    def _on_cancelled(self):
        self.status_label.setText(t('Загрузка остановлена; при повторе она продолжится с этого места.'))

    def _on_finished(self):
        self._thread = None
        self._set_running(False)

    def reject(self):
        """«Отмена» during a download stops it; the dialog closes when nothing is running."""
        if self._thread is not None:
            self.cancel()
            return
        super().reject()
