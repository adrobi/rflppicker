"""Shared helpers for the test suite."""
import os
import re

import rflp_core

HEADER = rflp_core.BASE_HEADER + rflp_core.ASSAY_HEADER + rflp_core.V12_HEADER


def core_row(**values) -> list:
    """A core result row in HEADER order; unspecified columns are None."""
    unknown = set(values) - set(HEADER)
    if unknown:
        raise KeyError(f'Unknown columns: {sorted(unknown)}')
    return [values.get(name) for name in HEADER]


def dispose(widget) -> None:
    """Close and really destroy a test widget now, in the main thread.

    Closed but undeleted windows piled up across the suite and left Qt's style
    sheet repolish touching freed memory (segfault on Linux, Python 3.12).
    """
    import gc

    from PyQt6 import sip
    from PyQt6.QtCore import QCoreApplication, QEvent

    if sip.isdeleted(widget):   # already destroyed with its parent
        return
    widget.close()
    widget.deleteLater()
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)
    gc.collect()


CYRILLIC = re.compile('[А-Яа-яЁё]')

# A fixed analysis of four SNPs for the window tests: a robust SNP with two assays, a «good gel only» SNP,
# an SNP without an assay and a data error.
FOUR_SNPS = ('chrom\tpos\tref\talt\tid\n11\t37525128\tC\tG\tHOXB13\n1\t41250567\tA\tC\tLEPR\n'
             '1\t91308651\tC\tT\tOLFML3\n12\t880121\tA\tG\trs7781\n')
_HOX = dict(variant='11:37525128 C>G', mapped_id='NC_056064.1', enzyme='Bse1I', isoschizomers='BsrI',
            genotype_quality='robust', genotype_margin=15.77, bands_ref_ref='448', bands_alt_alt='212; 236',
            product_size=475, assay_rank=1, primers_left='TATGAGCTTCTCTCCGCCAG', primers_right='CGCCTCCAAAGTAGCCGTAA',
            tm_left=59.0, tm_right=60.1)
FIXED_OUTPUT = (
    HEADER,
    [core_row(**_HOX),
     core_row(**dict(_HOX, enzyme='BsrDI', isoschizomers='', genotype_margin=12.1, assay_rank=2)),
     core_row(variant='1:41250567 A>C', enzyme='AciI', genotype_quality='nominal_only', genotype_margin=1.4,
              bands_ref_ref='292; 140', bands_alt_alt='432', product_size=432, assay_rank=1,
              weakest_genotypes='REF/REF–REF/ALT'),
     core_row(variant='1:91308651 C>T', enzyme='', status='no_discriminating_assay',
              reason='Не найдено сочетание для различения трёх генотипов')],
    {'summary': {'rows': 4}, 'rejected': [{'chrom': '12', 'pos': 880121, 'ref': 'A', 'alt': 'G',
                                            'reason': 'ref_not_match', 'details': 'x'}]},
)


def scene_texts(scene) -> list:
    """Texts and tooltips of the items of a QGraphicsScene."""
    texts = []
    for item in scene.items():
        if hasattr(item, 'toPlainText'):
            texts.append(item.toPlainText())
        texts.append(item.toolTip())
    return texts


def cyrillic_texts(root) -> list:
    """Texts under a widget that contain Cyrillic letters: titles, tooltips, placeholders, menus, tabs, scenes."""
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import (
        QAbstractButton, QAbstractSpinBox, QComboBox, QGraphicsView, QLabel, QLineEdit, QPlainTextEdit, QTableView,
        QTabWidget, QTextEdit, QToolButton, QWidget,
    )

    texts = [root.windowTitle(), root.toolTip()]
    for widget in [root] + root.findChildren(QWidget):
        texts.append(widget.toolTip())
        if isinstance(widget, (QLabel, QAbstractButton)):
            texts.append(widget.text())
        if isinstance(widget, QAbstractSpinBox):
            texts.append(widget.specialValueText())
        if isinstance(widget, QComboBox):
            texts += [widget.itemText(i) for i in range(widget.count())]
        if isinstance(widget, (QLineEdit, QTextEdit)):
            texts.append(widget.placeholderText())
        if isinstance(widget, QTextEdit):
            texts.append(widget.toPlainText())
        if isinstance(widget, QToolButton) and widget.menu():
            texts += [action.text() for action in widget.menu().actions()]
        if isinstance(widget, QTabWidget):
            texts += [widget.tabText(i) for i in range(widget.count())]
        if isinstance(widget, QGraphicsView) and widget.scene() is not None:
            texts += scene_texts(widget.scene())
        if isinstance(widget, QPlainTextEdit):
            texts += [widget.toPlainText(), widget.placeholderText()]
        if isinstance(widget, QTableView) and widget.model() is not None:
            model = widget.model()
            for column in range(model.columnCount()):
                texts.append(model.headerData(column, Qt.Orientation.Horizontal) or '')
                for row in range(model.rowCount()):
                    for role in (Qt.ItemDataRole.DisplayRole, Qt.ItemDataRole.ToolTipRole):
                        value = model.index(row, column).data(role)
                        texts.append(value if isinstance(value, str) else '')
    return [text for text in texts if CYRILLIC.search(text or '')]


def run_fake_analysis(window, folder: str, fasta: str, output, variants: str = FOUR_SNPS, log=()) -> None:
    """Run the window's analysis with the core replaced by a fixed output (or exception) and log lines."""
    from unittest.mock import patch

    from PyQt6.QtWidgets import QApplication

    path = write_text(folder, 'variants.tsv', variants)
    window.launch.set_files(path, fasta)

    def fake(params, status_cb=None, progress_cb=None, cancel_cb=None, report=None):
        for line in log:
            status_cb(line)
        if isinstance(output, Exception):
            raise output
        header, rows, fixed_report = output
        if report is not None:
            report.update(fixed_report)
        if progress_cb:
            progress_cb(1, 4)
        return header, rows

    with patch('rflp_core.run_rflp_gui_mode', side_effect=fake):
        window.start_run()
        worker = window.worker
        if worker is None:
            raise AssertionError('The analysis did not start')
        worker.wait(10000)
        for _ in range(3):
            QApplication.processEvents()
    if window.worker is not None:
        raise AssertionError('The analysis did not finish')


def write_text(directory: str, name: str, text: str) -> str:
    """Write UTF-8 text to directory/name and return the path."""
    path = os.path.join(directory, name)
    with open(path, 'w', encoding='utf-8', newline='') as handle:
        handle.write(text)
    return path
