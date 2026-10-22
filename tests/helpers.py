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
    from PyQt6.QtWidgets import (
        QAbstractButton, QAbstractSpinBox, QComboBox, QGraphicsView, QLabel, QLineEdit, QTabWidget, QTextEdit,
        QToolButton, QWidget,
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
    return [text for text in texts if CYRILLIC.search(text or '')]


def write_text(directory: str, name: str, text: str) -> str:
    """Write UTF-8 text to directory/name and return the path."""
    path = os.path.join(directory, name)
    with open(path, 'w', encoding='utf-8', newline='') as handle:
        handle.write(text)
    return path
