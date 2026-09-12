"""Interactive, exportable RFLP diagrams; coordinates are relative to the amplicon."""

from __future__ import annotations

import html
import math
import re
from pathlib import Path
from typing import Optional, Tuple

from Bio.Restriction import AllEnzymes
from PyQt6 import QtCore, QtGui
from PyQt6.QtCore import Qt
from PyQt6.QtSvg import QSvgGenerator
from PyQt6.QtWidgets import (
    QApplication, QDialog, QFileDialog, QGraphicsScene, QGraphicsView,
    QHBoxLayout, QLabel, QMessageBox, QPushButton, QSizePolicy, QSplitter,
    QTextEdit, QVBoxLayout, QWidget,
)

from rflp_core import cut_positions, revcomp


def _safe_int(value: object) -> Optional[int]:
    try:
        number = float(str(value).strip())
        return int(number) if math.isfinite(number) and number.is_integer() else None
    except (ValueError, TypeError, OverflowError):
        return None


class DiagramView(QGraphicsView):
    """Fit on resize until the user explicitly zooms the diagram."""

    zoom_changed = QtCore.pyqtSignal(float)

    def __init__(self, scene, parent=None):
        super().__init__(scene, parent)
        self.auto_fit = True
        self._fit_scale = 1.0
        self.setRenderHints(QtGui.QPainter.RenderHint.Antialiasing |
                            QtGui.QPainter.RenderHint.TextAntialiasing)
        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setDragMode(QGraphicsView.DragMode.ScrollHandDrag)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.setMinimumHeight(300)

    def fit_diagram(self):
        self.auto_fit = True
        if self.scene() and not self.scene().sceneRect().isEmpty():
            self.resetTransform()
            self.fitInView(self.scene().sceneRect(), Qt.AspectRatioMode.KeepAspectRatio)
            self._fit_scale = self.transform().m11() or 1.0
        self.zoom_changed.emit(1.0)

    def zoom_by(self, factor: float):
        if not self.scene() or not self.scene().items():
            return
        relative = self.transform().m11() / self._fit_scale
        target = max(0.25, min(8.0, relative * factor))
        self.auto_fit = False
        self.scale(target / relative, target / relative)
        self.zoom_changed.emit(target)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if self.auto_fit:
            self.fit_diagram()

    def wheelEvent(self, event):
        if event.modifiers() & Qt.KeyboardModifier.ControlModifier:
            self.zoom_by(1.2 if event.angleDelta().y() > 0 else 1 / 1.2)
            event.accept()
        else:
            super().wheelEvent(event)


class PrimerGraphicPanel(QWidget):
    """REF/ALT tracks, verified primer binding positions and full sequence copying.

    ``genome_start`` / ``genome_end`` are one-based inclusive. ``snp_offset`` is
    zero-based. Primer3 starts refer to its original template; the amplicon
    begins at ``primer_left_start``. Restriction cuts are interbase boundaries.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self._last_payload = None
        self._seq_ref = ""
        self._seq_alt = None
        self._primers = ("", "")
        self._cuts_ref = []
        self._cuts_alt = []
        self._primer_positions = (-1, -1)
        self._warnings = []
        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 12, 16, 12)
        layout.setSpacing(10)
        self.title = QLabel("Карта ампликона")
        font = self.title.font()
        font.setPointSize(17)
        font.setBold(True)
        self.title.setFont(font)
        self.meta = QLabel()
        self.meta.setTextFormat(Qt.TextFormat.PlainText)
        self.meta.setWordWrap(True)
        layout.addWidget(self.title)
        layout.addWidget(self.meta)
        self.scene = QGraphicsScene(self)
        self.view = DiagramView(self.scene, self)
        toolbar = QHBoxLayout()
        self.zoom_out_btn = self._button("−", "Уменьшить схему", lambda: self.view.zoom_by(1 / 1.25))
        self.zoom_in_btn = self._button("+", "Увеличить схему; Ctrl + колесо мыши", lambda: self.view.zoom_by(1.25))
        self.fit_btn = self._button("Вписать", "Сбросить масштаб и показать всю схему", self.view.fit_diagram)
        self.zoom_label = QLabel("100%")
        self.zoom_label.setMinimumWidth(48)
        self.view.zoom_changed.connect(lambda value: self.zoom_label.setText(f"{value:.0%}"))
        for widget in (self.zoom_out_btn, self.zoom_label, self.zoom_in_btn, self.fit_btn):
            toolbar.addWidget(widget)
        toolbar.addStretch()
        self.export_png_btn = self._button("Сохранить PNG", "Вся схема в высоком разрешении", lambda: self._export_dialog("png"))
        self.export_svg_btn = self._button("Сохранить SVG", "Векторная схема для печати и редактирования", lambda: self._export_dialog("svg"))
        toolbar.addWidget(self.export_png_btn)
        toolbar.addWidget(self.export_svg_btn)
        layout.addLayout(toolbar)
        sequence_panel = QWidget()
        sequence_layout = QVBoxLayout(sequence_panel)
        sequence_layout.setContentsMargins(0, 4, 0, 0)
        copy_row = QHBoxLayout()
        copy_row.addWidget(QLabel("Последовательности 5′ → 3′"))
        copy_row.addStretch()
        self.copy_ref_btn = self._button("Копировать REF", "Полная последовательность REF без разметки", lambda: self.copy_sequence("REF"))
        self.copy_alt_btn = self._button("Копировать ALT", "Полная последовательность ALT без разметки", lambda: self.copy_sequence("ALT"))
        self.copy_primers_btn = self._button("Копировать праймеры", "Оба праймера в направлении 5′ → 3′", self.copy_primers)
        for widget in (self.copy_ref_btn, self.copy_alt_btn, self.copy_primers_btn):
            copy_row.addWidget(widget)
        sequence_layout.addLayout(copy_row)
        self.seq_view = QTextEdit()
        self.seq_view.setReadOnly(True)
        self.seq_view.setFont(QtGui.QFont("Consolas", 10))
        self.seq_view.setMinimumHeight(125)
        sequence_layout.addWidget(self.seq_view)
        self.status_label = QLabel("Ctrl + колесо — масштаб; перетаскивание — перемещение схемы.")
        sequence_layout.addWidget(self.status_label)
        splitter = QSplitter(Qt.Orientation.Vertical)
        splitter.setChildrenCollapsible(False)
        splitter.addWidget(self.view)
        splitter.addWidget(sequence_panel)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)
        layout.addWidget(splitter, 1)
        self._apply_theme()
        self.show_empty()

    @staticmethod
    def _button(text, tooltip, callback):
        button = QPushButton(text)
        button.setToolTip(tooltip)
        button.setMinimumHeight(32)
        button.clicked.connect(callback)
        return button

    def _is_dark(self):
        return self.palette().color(QtGui.QPalette.ColorRole.Base).lightness() < 128

    def _colors(self):
        dark = self._is_dark()
        return {
            "dark": dark,
            "background": "#18212d" if dark else "#ffffff",
            "text": "#e6edf5" if dark else "#233247",
            "muted": "#adbbcd" if dark else "#5b6c81",
            "line": "#73859a" if dark else "#8c9bae",
            "grid": "#3c4b60" if dark else "#e1e8f0",
            "card": "#222f40" if dark else "#f5f8fc",
            "left": "#62d6b4" if dark else "#13846b",
            "right": "#87b8ff" if dark else "#2867bf",
            "snp": "#f5c367" if dark else "#aa6300",
            "ref": "#ff9aa8" if dark else "#c7435a",
            "alt": "#72cfdf" if dark else "#167d92",
            "hl_left": "#215344" if dark else "#d8f3e8",
            "hl_right": "#24466b" if dark else "#dfebff",
            "hl_snp": "#694d20" if dark else "#ffebba",
        }

    def _apply_theme(self):
        colors = self._colors()
        self.view.setBackgroundBrush(QtGui.QColor(colors["background"]))
        self.scene.setBackgroundBrush(QtGui.QColor(colors["background"]))
        border = colors["grid"]
        self.view.setStyleSheet(f"QGraphicsView {{ border: 1px solid {border}; border-radius: 8px; }}")
        self.seq_view.setStyleSheet(f"QTextEdit {{ border: 1px solid {border}; border-radius: 8px; padding: 8px; }}")
        self.status_label.setStyleSheet(f"color: {colors['muted']};")

    def changeEvent(self, event):
        super().changeEvent(event)
        if event.type() in (QtCore.QEvent.Type.PaletteChange, QtCore.QEvent.Type.ApplicationPaletteChange) and hasattr(self, "view"):
            self._apply_theme()
            if self._last_payload:
                self._render_impl(**self._last_payload)

    def _set_enabled(self, available):
        for button in (self.zoom_out_btn, self.zoom_in_btn, self.fit_btn, self.export_png_btn, self.export_svg_btn, self.copy_ref_btn):
            button.setEnabled(available)
        self.copy_alt_btn.setEnabled(available and self._seq_alt is not None)
        self.copy_primers_btn.setEnabled(available and any(self._primers))

    def show_empty(self):
        self.show_message("Выберите результат в таблице, чтобы увидеть ампликон и сравнить REF с ALT.")

    def show_message(self, msg: str):
        self._last_payload = None
        self._seq_ref, self._seq_alt = "", None
        self._primers = ("", "")
        self._cuts_ref, self._cuts_alt = [], []
        self._primer_positions = (-1, -1)
        self._warnings = []
        self.meta.setText(msg)
        self.scene.clear()
        self.scene.setSceneRect(0, 0, 1200, 590)
        self._draw_text(50, 220, "Карта появится после выбора результата", 17, True)
        self.seq_view.setPlainText(msg)
        self.status_label.setText("Ctrl + колесо — масштаб; перетаскивание — перемещение схемы.")
        self._set_enabled(False)
        self.view.fit_diagram()

    def render(self, mapped_id: str, genome_start: int, genome_end: int,
               amplicon_seq: str, snp_offset: Optional[int], primer_left: str,
               primer_right: str, tm_left: Optional[float], tm_right: Optional[float],
               primer_left_start: Optional[int] = None, primer_left_len: Optional[int] = None,
               primer_right_start: Optional[int] = None, primer_right_len: Optional[int] = None,
               variant_str: str = "", enzyme_name: str = "", pattern: str = "",
               frags_ref: str = "", frags_alt: str = ""):
        payload = dict(locals())
        payload.pop("self")
        self._last_payload = payload
        self._render_impl(**payload)

    def _draw_text(self, x, y, text, size=10, bold=False, color=None, width=None):
        item = self.scene.addText(str(text), QtGui.QFont("Segoe UI", size))
        font = item.font()
        font.setBold(bold)
        item.setFont(font)
        item.setDefaultTextColor(QtGui.QColor(color or self._colors()["text"]))
        if width is not None:
            item.setTextWidth(width)
        item.setPos(x, y)
        return item

    @staticmethod
    def _html_escape(text):
        return html.escape(str(text))

    @staticmethod
    def _parse_variant(variant_str: str) -> Tuple[str, str]:
        match = re.search(r"(?:^|\s)([ACGTN]+)\s*>\s*([ACGTN-]+)(?:$|\s)", (variant_str or "").strip(), re.IGNORECASE)
        return (match[1].upper(), match[2].upper()) if match else ("", "")

    @staticmethod
    def _enzyme_by_name(name):
        name = (name or "").strip()
        return next((enzyme for enzyme in AllEnzymes if enzyme.__name__ == name), None)

    def _map_primer_positions(self, seq, primer_left, primer_right,
                             primer_left_start, primer_left_len,
                             primer_right_start, primer_right_len):
        """Prefer verified Primer3 offsets; never invent a binding site."""
        left = (primer_left or "").upper().strip()
        right_binding = revcomp((primer_right or "").upper().strip())
        left_pos = seq.find(left) if left else -1
        right_pos = seq.rfind(right_binding) if right_binding else -1
        left_start, right_start = _safe_int(primer_left_start), _safe_int(primer_right_start)
        if left_start is not None and left and seq.startswith(left):
            left_pos = 0
        if left_start is not None and right_start is not None and right_binding:
            expected = right_start - left_start - len(right_binding) + 1
            if 0 <= expected <= len(seq) - len(right_binding) and seq[expected:expected + len(right_binding)] == right_binding:
                right_pos = expected
        return left_pos, right_pos

    @staticmethod
    def _make_alt_seq(seq_ref, snp_offset, ref, alt):
        if snp_offset is None or not 0 <= snp_offset < len(seq_ref):
            return None, "ALT недоступна: позиция SNP отсутствует или находится за пределами ампликона."
        if len(ref) != 1 or len(alt) != 1 or ref not in "ACGT" or alt not in "ACGT":
            return None, "ALT недоступна: для двух дорожек требуется одиночная замена A/C/G/T."
        if seq_ref[snp_offset] != ref:
            return None, f"ALT недоступна: REF {ref} не совпадает с основанием {seq_ref[snp_offset]} в FASTA."
        return seq_ref[:snp_offset] + alt + seq_ref[snp_offset + 1:], ""

    def _render_impl(self, **payload):
        self._apply_theme()
        sequence = re.sub(r"\s+", "", payload["amplicon_seq"] or "").upper()
        if not sequence:
            self.show_message("Не удалось получить последовательность ампликона.")
            return
        colors = self._colors()
        length = len(sequence)
        offset = _safe_int(payload["snp_offset"])
        ref, alt = self._parse_variant(payload["variant_str"])
        alt_sequence, alt_error = self._make_alt_seq(sequence, offset, ref, alt)
        left, right = ((payload[key] or "").upper().strip() for key in ("primer_left", "primer_right"))
        left_pos, right_pos = self._map_primer_positions(sequence, left, right, payload["primer_left_start"], payload["primer_left_len"], payload["primer_right_start"], payload["primer_right_len"])
        enzyme = self._enzyme_by_name(payload["enzyme_name"])
        self._seq_ref, self._seq_alt = sequence, alt_sequence
        self._primers = left, right
        self._primer_positions = left_pos, right_pos
        self._cuts_ref = cut_positions(enzyme, sequence) if enzyme else [0, length]
        self._cuts_alt = cut_positions(enzyme, alt_sequence) if enzyme and alt_sequence is not None else ([0, length] if alt_sequence else [])
        self._warnings = [alt_error] if alt_error else []
        if left and left_pos < 0:
            self._warnings.append("Левый праймер не найден в ампликоне.")
        if right and right_pos < 0:
            self._warnings.append("Сайт посадки правого праймера не найден в ампликоне.")
        if payload["enzyme_name"] and enzyme is None:
            self._warnings.append("Фермент не найден в базе Biopython; разрезы не показаны.")
        if payload["genome_end"] - payload["genome_start"] + 1 != length:
            self._warnings.append("Длина последовательности не совпадает с геномными координатами.")
        self.meta.setText("\n".join(self._warnings) if self._warnings else "Сравните разрезы и длины фрагментов; наведите курсор на метку для точных координат.")
        self.status_label.setText("Ctrl + колесо — масштаб; перетаскивание — перемещение схемы.")
        self.scene.clear()
        self.scene.setSceneRect(0, 0, 1200, 590)
        c = colors
        x0, x1 = 100.0, 1150.0

        def x_at(boundary):
            return x0 + float(boundary) / length * (x1 - x0)

        no_pen = QtGui.QPen(Qt.PenStyle.NoPen)
        self._draw_text(30, 15, payload["variant_str"] or "Карта ампликона", 17, True, width=790)
        self._draw_text(865, 19, payload["enzyme_name"] or "Фермент не выбран", 13, True, width=305)
        self._draw_text(30, 52, f"{payload['mapped_id']}:{payload['genome_start']}–{payload['genome_end']}  ·  {length} п.н.  ·  геномные координаты 1-based", 10, color=c["muted"], width=800)
        pattern = {"gain/loss": "Появление / исчезновение сайта", "shift": "Изменение фрагментов"}.get(payload["pattern"], payload["pattern"])
        self._draw_text(865, 51, pattern, 9, color=c["muted"], width=305)

        def tm_text(value):
            try:
                number = float(value)
                return f" · Tm {number:.1f} °C" if math.isfinite(number) else ""
            except (ValueError, TypeError):
                return ""

        def primer_arrow(position, primer, y, direction, color, label, tm):
            self._draw_text(100, y - 34, f"{label}  5′ → 3′  ·  {len(primer)} п.н.{tm_text(tm)}" if primer else f"{label}: не подобран", 10, True, color=color, width=1020)
            if not primer or position < 0:
                return
            start, end = x_at(position), x_at(position + len(primer))
            pen = QtGui.QPen(QtGui.QColor(color), 5)
            pen.setCapStyle(Qt.PenCapStyle.RoundCap)
            line = self.scene.addLine(start, y, end, y, pen)
            line.setToolTip(f"{label}: {primer}\nПосадка: {position + 1}–{position + len(primer)} в ампликоне (1-based)")
            tip = end if direction == 1 else start
            head_size = min(10, max(4, abs(end - start) / 2))
            polygon = QtGui.QPolygonF([QtCore.QPointF(tip, y), QtCore.QPointF(tip - direction * head_size, y - 7), QtCore.QPointF(tip - direction * head_size, y + 7)])
            self.scene.addPolygon(polygon, no_pen, QtGui.QBrush(QtGui.QColor(color)))

        primer_arrow(left_pos, left, 124, 1, c["left"], "Левый праймер →", payload["tm_left"])
        primer_arrow(right_pos, right, 174, -1, c["right"], "← Правый праймер", payload["tm_right"])

        def track(label, y, cuts, color, available):
            self.scene.addRect(25, y - 39, 1150, 106, no_pen, QtGui.QBrush(QtGui.QColor(c["card"])))
            self._draw_text(37, y - 17, label, 14, True, color=color)
            pen = QtGui.QPen(QtGui.QColor(c["line"]), 3)
            if not available:
                pen.setStyle(Qt.PenStyle.DashLine)
            self.scene.addLine(x0, y, x1, y, pen)
            if not available:
                self._draw_text(x0, y + 15, "Дорожка недоступна — см. сообщение над схемой", 10, color=c["muted"])
                return
            positions = sorted(set(position for position in cuts if 0 <= position <= length))
            fragments = [b - a for a, b in zip(positions, positions[1:])]
            last_label_right = -1000.0
            for position in positions:
                x = x_at(position)
                marker = self.scene.addLine(x, y - 13, x, y + 13, QtGui.QPen(QtGui.QColor(color), 2))
                marker.setToolTip(f"{label}: после {position} п.н. от начала ампликона" if 0 < position < length else f"Граница ампликона: {position} п.н.")
                if 0 < position < length and x - 19 >= last_label_right:
                    text = self._draw_text(x - 19, y - 39, str(position), 9, color=color)
                    text.setToolTip(marker.toolTip())
                    last_label_right = text.sceneBoundingRect().right() + 8
            for start, end in zip(positions, positions[1:]):
                left_x, right_x = x_at(start), x_at(end)
                value = str(end - start)
                metrics = QtGui.QFontMetricsF(QtGui.QFont("Segoe UI", 10))
                text_width = metrics.horizontalAdvance(value) + 8
                if right_x - left_x > text_width + 12:
                    self._draw_text((left_x + right_x - text_width) / 2, y + 7, value, 10, True, color=color)
            detail = " + ".join(str(value) for value in fragments)
            summary = f"Фрагменты: {detail} п.н." if enzyme else "Разрезы не рассчитаны: фермент не выбран или неизвестен."
            summary_item = self._draw_text(x0, y + 36, summary, 9, color=c["muted"], width=1040)
            summary_item.setToolTip(summary)
            if summary_item.boundingRect().height() > 28:
                summary_item.setPlainText(f"Фрагментов: {len(fragments)}; полные длины — в последовательностях и подсказке.")

        track("REF", 245, self._cuts_ref, c["ref"], True)
        track("ALT", 360, self._cuts_alt, c["alt"], alt_sequence is not None)
        if offset is not None and 0 <= offset < length:
            x = x_at(offset + 0.5)
            marker = self.scene.addLine(x, 203, x, 382, QtGui.QPen(QtGui.QColor(c["snp"]), 2, Qt.PenStyle.DashLine))
            marker.setToolTip(f"SNP: основание {offset + 1} в ампликоне; геномная позиция {payload['genome_start'] + offset}")
            for y in (245, 360):
                self.scene.addEllipse(x - 4, y - 4, 8, 8, no_pen, QtGui.QBrush(QtGui.QColor(c["snp"])))
        axis_y = 445
        self.scene.addLine(x0, axis_y, x1, axis_y, QtGui.QPen(QtGui.QColor(c["grid"]), 1))
        raw_step = max(1, length / 10)
        power = 10 ** math.floor(math.log10(raw_step))
        step = next(multiplier * power for multiplier in (1, 2, 5, 10) if multiplier * power >= raw_step)
        ticks = list(range(0, length + 1, int(step)))
        if ticks[-1] != length:
            if length - ticks[-1] < step * 0.6 and len(ticks) > 1:
                ticks.pop()
            ticks.append(length)
        for position in ticks:
            x = x_at(position)
            self.scene.addLine(x, axis_y - 4, x, axis_y + 4, QtGui.QPen(QtGui.QColor(c["line"]), 1))
            self._draw_text(x - 13, axis_y + 6, str(position), 9, color=c["muted"])
        self._draw_text(100, 477, "Шкала и метки разрезов: п.н. от начала ампликона · длины фрагментов указаны под дорожками", 9, color=c["muted"])
        if alt_sequence is not None:
            start, stop = max(0, offset - 12), min(length, offset + 13)
            self._draw_text(30, 518, f"SNP  {ref} → {alt}  ·  основание {offset + 1}", 11, True, color=c["snp"], width=390)
            context = f"REF  {sequence[start:stop]}\nALT  {alt_sequence[start:stop]}"
            item = self._draw_text(475, 510, context, 11)
            item.setFont(QtGui.QFont("Consolas", 11))
            item.setToolTip(f"Основания {start + 1}–{stop} ампликона (1-based)")
        else:
            self._draw_text(30, 520, "REF: референсная последовательность · ALT: последовательность после замены", 10, color=c["muted"])
        self._render_sequences(payload, offset, left_pos, right_pos, colors)
        self._set_enabled(True)
        self.view.fit_diagram()

    def _render_sequences(self, payload, offset, left_pos, right_pos, colors):
        sequence, alt = self._seq_ref, self._seq_alt
        length = len(sequence)
        left, right = self._primers
        begin, end = 0, length
        if length > 1500:
            center = offset if offset is not None and 0 <= offset < length else length // 2
            begin, end = max(0, center - 240), min(length, center + 241)
        lines = [f"<div style='color:{colors['text']};font-family:Consolas,monospace;font-size:10pt'>"]
        if begin or end != length:
            lines.append(f"<p>Показаны основания {begin + 1}–{end} из {length}. Кнопки копируют полную последовательность.</p>")
        lines.append(f"<p><span style='background-color:{colors['hl_left']}'>Левый праймер</span> · <span style='background-color:{colors['hl_right']}'>Правый праймер</span> · <span style='background-color:{colors['hl_snp']}'>SNP</span></p>")
        if left or right:
            lines.append(f"<pre>Левый  5′→3′: {html.escape(left) or '—'}\nПравый 5′→3′: {html.escape(right) or '—'}</pre>")
        if payload["enzyme_name"] and self._enzyme_by_name(payload["enzyme_name"]):
            for label, cuts in (("REF", self._cuts_ref), ("ALT", self._cuts_alt)):
                if cuts:
                    fragments = " + ".join(str(b - a) for a, b in zip(cuts, cuts[1:]))
                    lines.append(f"<p>{label}, фрагменты: {fragments} п.н.</p>")
        for start in range(begin, end, 70):
            stop = min(start + 70, end)
            lines.append(f"<p style='color:{colors['muted']}'>Ампликон {start + 1}–{stop} · геном {payload['genome_start'] + start}–{payload['genome_start'] + stop - 1}</p><pre>")
            for label, bases in (("REF", sequence), ("ALT", alt)):
                if bases is None:
                    continue
                tokens = []
                for position in range(start, stop):
                    background = None
                    if position == offset:
                        background = colors["hl_snp"]
                    elif left_pos >= 0 and left_pos <= position < left_pos + len(left):
                        background = colors["hl_left"]
                    elif right_pos >= 0 and right_pos <= position < right_pos + len(right):
                        background = colors["hl_right"]
                    base = html.escape(bases[position])
                    tokens.append(f"<span style='background-color:{background};font-weight:bold'>{base}</span>" if background else base)
                lines.append(f"{label}  {''.join(tokens)}\n")
            lines.append("</pre>")
        lines.append("</div>")
        self.seq_view.setHtml("".join(lines))

    def copy_sequence(self, allele: str):
        if allele not in ("REF", "ALT"):
            raise ValueError("Allele must be REF or ALT")
        sequence = self._seq_ref if allele == "REF" else self._seq_alt
        if not sequence:
            return False
        QApplication.clipboard().setText(sequence)
        self.status_label.setText(f"{allele}: скопировано {len(sequence)} п.н.")
        return True

    def copy_primers(self):
        left, right = self._primers
        if not left and not right:
            return False
        text = "\n".join(f">{label} (5'->3')\n{sequence}" for label, sequence in (("LEFT", left), ("RIGHT", right)) if sequence)
        QApplication.clipboard().setText(text)
        self.status_label.setText("Праймеры скопированы в формате FASTA (5′ → 3′).")
        return True

    def export_png(self, path, scale=2.0):
        """Save the whole scene, independently of zoom/pan. Return the saved path."""
        if not self._seq_ref:
            raise ValueError("Сначала выберите результат для визуализации.")
        if not math.isfinite(float(scale)) or not 0.5 <= scale <= 8:
            raise ValueError("Масштаб экспорта должен быть от 0.5 до 8.")
        rect = self.scene.sceneRect()
        image = QtGui.QImage(math.ceil(rect.width() * scale), math.ceil(rect.height() * scale), QtGui.QImage.Format.Format_ARGB32)
        image.fill(QtGui.QColor(self._colors()["background"]))
        painter = QtGui.QPainter(image)
        painter.setRenderHints(self.view.renderHints())
        try:
            self.scene.render(painter, QtCore.QRectF(0, 0, image.width(), image.height()), rect)
        finally:
            painter.end()
        if not image.save(str(path), "PNG"):
            raise OSError(f"Не удалось сохранить PNG: {path}")
        return str(path)

    def export_svg(self, path):
        if not self._seq_ref:
            raise ValueError("Сначала выберите результат для визуализации.")
        output = QtCore.QFile(str(path))
        if not output.open(QtCore.QIODevice.OpenModeFlag.WriteOnly):
            raise OSError(f"Не удалось сохранить SVG: {output.errorString()}")
        generator = QSvgGenerator()
        generator.setOutputDevice(output)
        rect = self.scene.sceneRect()
        generator.setSize(rect.size().toSize())
        generator.setViewBox(rect)
        generator.setTitle("RFLP · REF / ALT")
        generator.setDescription("Ампликон, праймеры 5′→3′, SNP и разрезы рестриктазы. Координаты генома 1-based; разрезы — п.н. от начала ампликона.")
        painter = QtGui.QPainter()
        try:
            if not painter.begin(generator):
                raise OSError(f"Не удалось начать экспорт SVG: {path}")
            painter.setRenderHints(self.view.renderHints())
            self.scene.render(painter, rect, rect)
        finally:
            if painter.isActive():
                painter.end()
            output.close()
        return str(path)

    def _export_dialog(self, extension):
        payload = self._last_payload or {}
        name = re.sub(r"[^\w.-]+", "_", f"{payload.get('variant_str', 'amplicon')}_{payload.get('enzyme_name', '')}").strip("_")
        path, _ = QFileDialog.getSaveFileName(self, "Сохранить карту ампликона", f"{name}.{extension}", f"{extension.upper()} (*.{extension})")
        if not path:
            return
        if not Path(path).suffix:
            path += f".{extension}"
        try:
            (self.export_png if extension == "png" else self.export_svg)(path)
        except (OSError, ValueError) as error:
            QMessageBox.warning(self, "Экспорт схемы", str(error))
        else:
            self.status_label.setText(f"Сохранено: {path}")


class RFLPVisualizationDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Карта ампликона · RFLP REF / ALT")
        self.resize(1240, 930)
        self.setMinimumSize(900, 700)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 12)
        self.panel = PrimerGraphicPanel(self)
        layout.addWidget(self.panel, 1)
        row = QHBoxLayout()
        row.addStretch()
        self.close_btn = QPushButton("Закрыть")
        self.close_btn.setMinimumHeight(32)
        self.close_btn.clicked.connect(self.close)
        row.addWidget(self.close_btn)
        row.addSpacing(16)
        layout.addLayout(row)

    def show_payload(self, payload: Optional[dict], message_if_empty: str = ""):
        if payload:
            self.panel.render(**payload)
        elif message_if_empty:
            self.panel.show_message(message_if_empty)
        else:
            self.panel.show_empty()
