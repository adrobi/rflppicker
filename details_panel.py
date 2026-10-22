"""The «Детали» column: schemes of the selected SNP, the gel, the map and the primer card (spec §7)."""
import dataclasses
from typing import List, Optional

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (
    QApplication, QGridLayout, QHBoxLayout, QLabel, QListWidget, QListWidgetItem, QPushButton, QSizePolicy,
    QTabWidget, QToolButton, QVBoxLayout, QWidget,
)

import icons
import palette
from i18n import N_, fmt, language, t
from snp_summary import SnpSummary, order_lines, status_title, summary_text, weakest_text
from theme import mono_font
from visualization import GelPanel, PrimerGraphicPanel, RFLPVisualizationDialog, payload_for_scheme, render_payload

EMPTY_TEXT = N_('Выберите SNP в списке')


def _isoschizomer_text(count: int) -> str:
    if language() == 'en':
        return t('+{count} изошизомер', count=count) if count == 1 else t('+{count} изошизомеров', count=count)
    if count % 10 == 1 and count % 100 != 11:
        word = 'изошизомер'
    elif 2 <= count % 10 <= 4 and not 12 <= count % 100 <= 14:
        word = 'изошизомера'
    else:
        word = 'изошизомеров'
    return f'+{count} {word}'


def _score(value) -> str:
    return fmt(float(value)) if isinstance(value, (int, float)) else ''


class DetailsPanel(QWidget):
    choice_changed = pyqtSignal(object)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._summary: Optional[SnpSummary] = None
        self._dcaps_mode = 'off'
        self._viewed = 0
        self._fasta_path = ''
        self._params: dict = {}
        self._fasta = None
        self._payload: Optional[dict] = None
        self._message = ''
        self._dialog: Optional[RFLPVisualizationDialog] = None
        self._colors = palette.colors(False)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.title_label = self._label('title_label', 'title')
        self.subtitle_label = self._label('subtitle_label', 'muted')
        header = QHBoxLayout()
        header.addWidget(self.title_label, 1)
        self.status_badge = self._label('status_badge')
        self.status_badge.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Fixed)
        header.addWidget(self.status_badge, 0, Qt.AlignmentFlag.AlignTop)
        layout.addLayout(header)
        layout.addWidget(self.subtitle_label)

        self.message_label = self._label('message_label')
        self.message_label.setWordWrap(True)
        layout.addWidget(self.message_label)
        self.raw_label = self._label('raw_label', 'muted')
        self.raw_label.setFont(mono_font())
        layout.addWidget(self.raw_label)

        self.schemes_caption = self._label('schemes_caption', 'caption')
        layout.addWidget(self.schemes_caption)
        self.schemes_list = QListWidget()
        self.schemes_list.setObjectName('schemes_list')
        self.schemes_list.setMaximumHeight(170)
        self.schemes_list.currentRowChanged.connect(self._on_row_changed)
        layout.addWidget(self.schemes_list)

        self.tabs = QTabWidget()
        self.tabs.setObjectName('tabs')
        gel_tab = QWidget()
        gel_layout = QVBoxLayout(gel_tab)
        gel_layout.setContentsMargins(0, 6, 0, 0)
        self.gel = GelPanel(compact=True)
        gel_layout.addWidget(self.gel, 1)
        self.weakest_label = self._label('weakest_label', 'muted')
        self.weakest_label.setWordWrap(True)
        gel_layout.addWidget(self.weakest_label)
        self.tabs.addTab(gel_tab, t('Гель'))
        self.map = PrimerGraphicPanel(compact=True)
        self.tabs.addTab(self.map, t('Карта'))
        self.tabs.addTab(self._build_primer_card(), t('Праймеры'))
        self.large_button = QToolButton()
        self.large_button.setObjectName('large_button')
        self.large_button.setText(t('⤢ Крупно'))
        self.large_button.setProperty('role', 'link')
        self.large_button.setToolTip(t('Карта с последовательностями и гель в отдельном окне, с сохранением PNG и SVG'))
        self.large_button.clicked.connect(self.open_large)
        self.tabs.setCornerWidget(self.large_button)
        layout.addWidget(self.tabs, 1)
        # Takes the free height when there are no tabs (SNP without a scheme or with an error).
        self._filler = QWidget()
        self._filler.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Expanding)
        layout.addWidget(self._filler, 1)

        self.apply_colors(self._colors)
        self.show_summary(None, 'off')

    @staticmethod
    def _label(name: str, role: str = '') -> QLabel:
        label = QLabel()
        label.setObjectName(name)
        if role:
            label.setProperty('role', role)
        label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        return label

    def _build_primer_card(self) -> QWidget:
        card = QWidget()
        card.setObjectName('primer_card')
        layout = QVBoxLayout(card)
        grid = QGridLayout()
        self._card_rows = []
        for row in range(2):
            name, sequence = QLabel(), QLabel()
            sequence.setFont(mono_font(10))
            sequence.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
            info = QLabel()
            info.setProperty('role', 'muted')
            grid.addWidget(name, row * 2, 0)
            grid.addWidget(sequence, row * 2, 1)
            grid.addWidget(info, row * 2 + 1, 1)
            self._card_rows.append((name, sequence, info))
        grid.setColumnStretch(1, 1)
        layout.addLayout(grid)
        self.card_details = QLabel()
        self.card_details.setWordWrap(True)
        layout.addWidget(self.card_details)
        buttons = QHBoxLayout()
        self.copy_order_button = QPushButton(t('Копировать для заказа'))
        self.copy_order_button.setObjectName('copy_order_button')
        self.copy_order_button.setToolTip(t('Две строки «имя ⇥ последовательность» для бланка заказа синтеза'))
        self.copy_order_button.clicked.connect(self._copy_order)
        self.copy_all_button = QPushButton(t('Копировать всё'))
        self.copy_all_button.setObjectName('copy_all_button')
        self.copy_all_button.setToolTip(t('Текстовая сводка схемы для лабораторного журнала'))
        self.copy_all_button.clicked.connect(self._copy_all)
        buttons.addWidget(self.copy_order_button)
        buttons.addWidget(self.copy_all_button)
        buttons.addStretch()
        layout.addLayout(buttons)
        layout.addStretch()
        return card

    # ----- genome -----
    def set_context(self, fasta_path: str, params: dict) -> None:
        self.close_fasta()
        self._fasta_path = fasta_path or ''
        self._params = dict(params)

    def close_fasta(self) -> None:
        if self._fasta is not None:
            self._fasta.close()
            self._fasta = None

    def _open_fasta(self):
        if self._fasta is None and self._fasta_path:
            try:
                from pyfaidx import Fasta
                self._fasta = Fasta(self._fasta_path, as_raw=True, sequence_always_upper=True)
            except Exception:
                self._fasta = None   # payload_for_scheme reports the error
        return self._fasta

    # ----- summary -----
    def _viewed_summary(self) -> Optional[SnpSummary]:
        """The summary as if the viewed scheme were chosen (for the card and copy texts)."""
        if self._summary is None or not self._summary.schemes:
            return self._summary
        return dataclasses.replace(self._summary, chosen=self._viewed)

    def show_summary(self, summary: Optional[SnpSummary], dcaps_mode: str) -> None:
        self._summary, self._dcaps_mode = summary, dcaps_mode
        has_schemes = bool(summary and summary.schemes)
        for widget in (self.schemes_caption, self.schemes_list, self.tabs):
            widget.setVisible(has_schemes)
        self._filler.setVisible(not has_schemes)
        self.raw_label.hide()
        if summary is None:
            self.title_label.setText('')
            self.subtitle_label.setText('')
            self.status_badge.hide()
            self.message_label.setText(t(EMPTY_TEXT))
            self.message_label.show()
            self._show_payload(None, t(EMPTY_TEXT))
            return
        r = summary.record
        self.title_label.setText(f'{r.snp_id} · {summary.label}' if r.snp_id else summary.label)
        source = summary.scheme or summary.failure or {}
        mapped = source.get('mapped_id')
        parts = [mapped, t('строка {line} файла', line=r.line), summary.swap_text()]
        self.subtitle_label.setText(' · '.join(part for part in parts if part))
        self._update_badge()
        if not has_schemes:
            self.message_label.setText(summary.note(dcaps_mode))
            self.message_label.show()
            if summary.rejection and summary.rejection.get('reason') == 'parse_error':
                self.raw_label.setText(r.raw)
                self.raw_label.show()
            self._show_payload(None, summary.note(dcaps_mode))
            return
        self.message_label.hide()
        self._fill_schemes()
        self.view_scheme(summary.chosen)

    def _update_badge(self):
        status = self._summary.status
        c = self._colors
        self.status_badge.setText(status_title(status))
        self.status_badge.setStyleSheet(f"color: {c[f'{status}_fg']}; background: {c[f'{status}_bg']}; "
                                        f"border-radius: 4px; padding: 2px 8px;")
        self.status_badge.show()

    def _fill_schemes(self):
        s = self._summary
        self.schemes_caption.setText(t('Схемы ({count}) · ★ идёт в Excel', count=len(s.schemes)))
        self.schemes_list.blockSignals(True)
        self.schemes_list.clear()
        self._row_texts: List[str] = []
        c = self._colors
        for index, scheme in enumerate(s.schemes):
            row = QWidget()
            layout = QHBoxLayout(row)
            layout.setContentsMargins(4, 2, 4, 2)
            star = QToolButton()
            star.setAutoRaise(True)
            chosen = index == s.chosen
            star.setIcon(icons.icon('star', c['accent'] if chosen else c['text2'], filled=chosen))
            star.setToolTip(t('Выбрана для Excel') if chosen else t('Сделать выбранной'))
            star.clicked.connect(lambda _=False, i=index: self.choose_scheme(i))
            layout.addWidget(star)
            parts = [scheme.get('enzyme') or '']
            others = [n for n in (scheme.get('isoschizomers') or '').split('; ') if n]
            if others:
                parts.append(_isoschizomer_text(len(others)))
            if scheme.get('design') == 'dCAPS':
                parts.append('dCAPS')
            if index == 0:
                parts.append(t('лучшая'))
            first = ' · '.join(parts)
            second = t('{ref_ref} | {alt_alt} · продукт {size}', ref_ref=scheme.get('bands_ref_ref') or '',
                       alt_alt=scheme.get('bands_alt_alt') or '', size=scheme.get('product_size'))
            text = QLabel(f'{first}\n{second}')
            if others:
                text.setToolTip(t('Изошизомеры: {names}', names=', '.join(others)))
            layout.addWidget(text, 1)
            note = dataclasses.replace(s, chosen=index).note(self._dcaps_mode)
            warning = scheme.get('genotype_quality') == 'nominal_only' or any(
                scheme.get(k) for k in ('primer_site', 'background_site_variants'))
            if warning:
                alert = QLabel()
                alert.setPixmap(icons.icon('alert-triangle', c['nominal_fg']).pixmap(16, 16))
                alert.setToolTip(note)
                layout.addWidget(alert)
            score = QLabel(_score(scheme.get('genotype_margin')))
            layout.addWidget(score)
            item = QListWidgetItem()
            item.setSizeHint(row.sizeHint())
            self.schemes_list.addItem(item)
            self.schemes_list.setItemWidget(item, row)
            self._row_texts.append(f'{first}\n{second}')
        self.schemes_list.blockSignals(False)

    def scheme_texts(self) -> List[str]:
        return list(getattr(self, '_row_texts', []))

    def viewed_index(self) -> int:
        return self._viewed

    def _on_row_changed(self, row: int):
        if row >= 0:
            self.view_scheme(row)

    def view_scheme(self, index: int) -> None:
        s = self._summary
        if s is None or not s.schemes:
            return
        self._viewed = index
        self.schemes_list.blockSignals(True)
        self.schemes_list.setCurrentRow(index)
        self.schemes_list.blockSignals(False)
        scheme = s.schemes[index]
        self.weakest_label.setText(weakest_text(scheme))
        self._fill_card()
        payload, error = payload_for_scheme(scheme, self._fasta_path, self._params, fasta=self._open_fasta())
        self._show_payload(payload, error)

    def choose_scheme(self, index: int) -> None:
        s = self._summary
        if s is None:
            return
        s.choose(index)
        self._update_badge()
        self._fill_schemes()
        self.view_scheme(index)
        self.choice_changed.emit(s)

    def _fill_card(self):
        viewed = self._viewed_summary()
        scheme = viewed.scheme
        for (name, sequence, info), side, primer, tm in zip(
                self._card_rows, ('F', 'R'), ('primers_left', 'primers_right'), ('tm_left', 'tm_right')):
            text = scheme.get(primer) or ''
            name.setText(f'{viewed.primer_base}_{side}' if text else '')
            sequence.setText(text or '—')
            tm_value = scheme.get(tm)
            info.setText(t('{length} нт · Tm {tm} °C', length=len(text), tm=f'{tm_value:.1f}')
                         if text and tm_value is not None else '')
        details = [t('Продукт {size} п.н. · {enzyme} ({site})', size=scheme.get('product_size'),
                     enzyme=scheme.get('enzyme'), site=scheme.get('site'))]
        if scheme.get('suppliers'):
            details.append(t('Поставщики: {suppliers}', suppliers=scheme['suppliers']))
        if scheme.get('design') == 'dCAPS':
            details.append(t('Строчные буквы — намеренные несовпадения.'))
        chosen = self._viewed == self._summary.chosen
        if not chosen:
            details.append(t('Эта схема не выбрана: в Excel идёт ★ {enzyme}. Отметьте её звёздочкой, чтобы копировать.',
                             enzyme=self._summary.scheme.get('enzyme')))
        self.card_details.setText('\n'.join(details))
        # Copy only the ★ scheme: its primer names are the ones Excel uses.
        self.copy_order_button.setEnabled(chosen and bool(scheme.get('primers_left')))
        self.copy_all_button.setEnabled(chosen)

    def _copy_order(self):
        QApplication.clipboard().setText(order_lines(self._viewed_summary()))

    def _copy_all(self):
        QApplication.clipboard().setText(summary_text(self._viewed_summary(), self._dcaps_mode))

    # ----- drawing -----
    def _show_payload(self, payload: Optional[dict], message: str):
        self._payload, self._message = payload, message
        render_payload(self.map, self.gel, payload, message)
        if self._dialog is not None and self._dialog.isVisible():
            self._dialog.show_payload(payload, message)

    def gel_message(self) -> str:
        return self.gel.message()

    def open_large(self) -> None:
        if self._dialog is None:
            self._dialog = RFLPVisualizationDialog(self)
            self._dialog.apply_colors(self._colors)
        self._dialog.show()
        self._dialog.raise_()
        self._dialog.activateWindow()
        self._dialog.show_payload(self._payload, self._message)

    def large_dialog(self) -> Optional[RFLPVisualizationDialog]:
        return self._dialog

    def apply_colors(self, colors: dict) -> None:
        self._colors = dict(colors)
        self.map.apply_colors(colors)
        if self._dialog is not None:
            self._dialog.apply_colors(colors)
        if self._summary is not None:
            self._update_badge()
            if self._summary.schemes:
                self._fill_schemes()
                self.schemes_list.setCurrentRow(self._viewed)
