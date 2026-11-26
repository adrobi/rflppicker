"""Choose an rsID source and inspect verified matches before saving a target."""

from pathlib import Path

from PyQt6.QtCore import QThread, pyqtSignal
from PyQt6.QtWidgets import (QDialog, QFileDialog, QFormLayout, QHBoxLayout, QLabel, QLineEdit,
                             QPlainTextEdit, QPushButton, QVBoxLayout)

from i18n import core_text, t
from rflp_core import load_names2
from rsid_export import describe_lookup, write_rsid_target


class _LookupThread(QThread):
    done = pyqtSignal(object)
    failed = pyqtSignal(str)

    def __init__(self, query, options, parent):
        super().__init__(parent)
        self.query, self.options = query, options

    def run(self):
        from rsid_lookup import lookup_rsid

        def check_cancel():
            if self.isInterruptionRequested():
                raise InterruptedError()

        try:
            result = lookup_rsid(self.query, check_cancel=check_cancel, **self.options)
            check_cancel()
        except InterruptedError:
            return
        except (OSError, ValueError) as exc:
            self.failed.emit(core_text(str(exc)))
        else:
            self.done.emit(result)


class RsidDialog(QDialog):
    def __init__(self, fasta, assembly, names, parent=None):
        super().__init__(parent)
        self.setWindowTitle(t('Поиск по rsID…'))
        self.resize(690, 440)
        self.fasta, self.names = fasta, names
        self.target_path = ''
        self.result = None
        self.worker = None
        layout = QVBoxLayout(self)
        form = QFormLayout()
        layout.addLayout(form)
        self.query, self.taxid = QLineEdit(), QLineEdit()
        self.query.setPlaceholderText('rs422668280')
        self.taxid.setPlaceholderText('9940')
        self.assembly = QLineEdit(assembly)
        form.addRow('rsID', self.query)
        form.addRow(t('Вид (NCBI TaxID)'), self.taxid)
        form.addRow(t('Сборка'), self.assembly)
        self.vcf = QLineEdit()
        self.vcf.setPlaceholderText(t('Оставьте пустым для поиска в EVA'))
        self.vcf_button = QPushButton(t('Выбрать файл…'))
        self.vcf_button.clicked.connect(self._browse)
        source = QHBoxLayout()
        source.addWidget(self.vcf)
        source.addWidget(self.vcf_button)
        form.addRow(t('Локальный VCF'), source)
        self.result_text = QPlainTextEdit()
        self.result_text.setReadOnly(True)
        layout.addWidget(self.result_text, 1)
        self.status = QLabel()
        self.status.setWordWrap(True)
        layout.addWidget(self.status)
        buttons = QHBoxLayout()
        self.search = QPushButton(t('Найти rsID'))
        self.search.clicked.connect(self._search)
        buttons.addWidget(self.search)
        buttons.addStretch()
        self.use = QPushButton(t('Сохранить как входной VCF…'))
        self.use.setEnabled(False)
        self.use.clicked.connect(self._save)
        buttons.addWidget(self.use)
        close = QPushButton(t('Закрыть'))
        close.clicked.connect(self.reject)
        buttons.addWidget(close)
        layout.addLayout(buttons)
        for field in (self.query, self.taxid, self.assembly, self.vcf):
            field.textChanged.connect(self._invalidate)

    def _browse(self):
        path, _ = QFileDialog.getOpenFileName(self, t('Локальный VCF'), '',
                                             t('VCF (*.vcf *.vcf.gz *.gz);;Все файлы (*)'))
        if path:
            self.vcf.setText(path)

    def _invalidate(self):
        self.result = None
        self.use.setEnabled(False)
        self.result_text.clear()

    def _search(self):
        self._invalidate()
        self.status.setText(t('Поиск rsID…'))
        fields = (self.query, self.taxid, self.assembly, self.vcf, self.vcf_button, self.search)
        for field in fields:
            field.setEnabled(False)
        options = dict(fasta=self.fasta, taxid=self.taxid.text().strip(), assembly=self.assembly.text().strip(),
                       local_vcf=self.vcf.text().strip() or None, names=load_names2(self.names))
        self.worker = _LookupThread(self.query.text().strip(), options, self)
        self.worker.done.connect(self._receive)
        self.worker.failed.connect(self.status.setText)

        def finished():
            worker, self.worker = self.worker, None
            worker.deleteLater()
            for field in fields:
                field.setEnabled(True)
            self.use.setEnabled(self.result is not None and self.result.match is not None)

        self.worker.finished.connect(finished)
        self.worker.start()

    def _receive(self, result):
        self.result = result
        self.result_text.setPlainText(describe_lookup(result))
        self.status.clear()

    def _save(self):
        path, _ = QFileDialog.getSaveFileName(self, t('Сохранить как входной VCF…'), self.result.rsid + '.vcf',
                                             t('Файл VCF (*.vcf)'))
        if not path:
            return
        if Path(path).suffix.lower() != '.vcf':
            path += '.vcf'
        try:
            write_rsid_target(path, self.result, self.fasta, self.vcf.text().strip())
        except (OSError, ValueError) as exc:
            self.status.setText(core_text(str(exc)))
            return
        self.target_path = path
        self.accept()

    def reject(self):
        if self.worker and self.worker.isRunning():
            self.worker.requestInterruption()
            self.worker.finished.connect(super().reject)
            return
        super().reject()
