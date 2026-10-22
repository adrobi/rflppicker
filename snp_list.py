"""The «Список SNP» column: one row per input SNP with status filters and search (spec §6)."""
from typing import Dict, List, Optional, Set

from PyQt6.QtCore import QAbstractTableModel, QModelIndex, QSortFilterProxyModel, Qt, pyqtSignal
from PyQt6.QtGui import QColor, QIcon, QPainter, QPixmap
from PyQt6.QtWidgets import (
    QHBoxLayout, QHeaderView, QLabel, QLineEdit, QMenu, QPushButton, QTableView, QToolButton, QVBoxLayout, QWidget,
)

import icons
import palette
from app_settings import AppSettings
from i18n import N_, fmt, t
from result_model import LABELS
from snp_summary import ERROR, NONE, STATUS_ORDER, SnpSummary, status_title

SORT_ROLE = Qt.ItemDataRole.UserRole
_HIDDEN_KEYS = ('mapped_id', 'site', 'isoschizomers', 'suppliers', 'primers_left', 'primers_right', 'tm_left',
                'tm_right', 'bands_ref_alt', 'worst_margin', 'diag_delta_bp', 'design', 'repeat_pct',
                'amplicon_start', 'amplicon_end')
# Titles are Russian sources (LABELS too), translated with t() when shown.
COLUMNS = [('dot', ''), ('id', 'ID'), ('snp', 'SNP'), ('enzyme', N_('Фермент')),
           ('bands', N_('Полосы REF | ALT')), ('product', N_('Продукт')), ('score', N_('Оценка')),
           ('schemes', N_('Схем')), ('note', N_('Примечание'))] + [
    (key, LABELS.get(key, key)) for key in _HIDDEN_KEYS]
DEFAULT_COLUMNS = ('dot', 'id', 'snp', 'enzyme', 'bands', 'product', 'score', 'schemes', 'note')
_KEYS = [key for key, _ in COLUMNS]
_NUMERIC = {'product', 'score', 'schemes', 'tm_left', 'tm_right', 'worst_margin', 'diag_delta_bp', 'repeat_pct',
            'amplicon_start', 'amplicon_end'}
_WIDTHS = {'dot': 24, 'id': 90, 'snp': 140, 'enzyme': 90, 'bands': 130, 'product': 70, 'score': 64,
           'schemes': 52}


def _chrom_key(chrom: str):
    return (0, int(chrom), '') if chrom.isdigit() else (1, 0, chrom)


class SnpListModel(QAbstractTableModel):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._summaries: List[SnpSummary] = []
        self._dcaps_mode = 'off'
        self._colors = palette.colors(False)
        self._dots: Dict[str, QIcon] = {}

    def set_summaries(self, summaries: List[SnpSummary], dcaps_mode: str) -> None:
        self.beginResetModel()
        self._summaries = list(summaries)
        self._dcaps_mode = dcaps_mode
        self.endResetModel()

    def set_colors(self, colors: dict) -> None:
        self._colors = dict(colors)
        self._dots = {}
        if self._summaries:
            self.dataChanged.emit(self.index(0, 0), self.index(self.rowCount() - 1, self.columnCount() - 1))

    def summary(self, row: int) -> SnpSummary:
        return self._summaries[row]

    def summaries(self) -> List[SnpSummary]:
        return list(self._summaries)

    def refresh_row(self, row: int) -> None:
        self.dataChanged.emit(self.index(row, 0), self.index(row, self.columnCount() - 1))

    def counts(self) -> Dict[str, int]:
        counts = dict.fromkeys(STATUS_ORDER, 0)
        for s in self._summaries:
            counts[s.status] += 1
        return counts

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self._summaries)

    def columnCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(COLUMNS)

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):
        if orientation == Qt.Orientation.Horizontal and role in (Qt.ItemDataRole.DisplayRole,
                                                                  Qt.ItemDataRole.ToolTipRole):
            return t(COLUMNS[section][1])
        return None

    def _dot(self, status: str) -> QIcon:
        if status not in self._dots:
            pixmap = QPixmap(20, 20)
            pixmap.fill(Qt.GlobalColor.transparent)
            painter = QPainter(pixmap)
            painter.setRenderHint(QPainter.RenderHint.Antialiasing)
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QColor(self._colors[f'{status}_fg']))
            painter.drawEllipse(5, 5, 10, 10)
            painter.end()
            self._dots[status] = QIcon(pixmap)
        return self._dots[status]

    def _value(self, s: SnpSummary, key: str):
        """(display text, sort key) of one cell."""
        scheme = s.scheme or {}
        if key == 'dot':
            return '', STATUS_ORDER.index(s.status)
        if key == 'id':
            return s.record.snp_id, s.record.snp_id.lower() or None
        if key == 'snp':
            r = s.record
            return s.label, (_chrom_key(r.chrom), r.pos if r.pos is not None else -1)
        if key == 'enzyme':
            enzyme = scheme.get('enzyme') or ''
            return (f'{enzyme} ✎' if s.manual else enzyme), enzyme.lower() or None
        if key == 'bands':
            if not scheme:
                return '', None
            text = f"{scheme.get('bands_ref_ref') or ''} | {scheme.get('bands_alt_alt') or ''}"
            return text, text
        if key == 'product':
            value = scheme.get('product_size')
            return fmt(value), value
        if key == 'score':
            value = scheme.get('genotype_margin')
            return fmt(value), value
        if key == 'schemes':
            return (str(len(s.schemes)), len(s.schemes)) if s.schemes else ('', None)
        if key == 'note':
            note = s.note(self._dcaps_mode)
            return note, note.lower() or None
        value = scheme.get(key)
        if value is None or value == '':
            return '', None
        return (fmt(value), value) if isinstance(value, (int, float)) else (str(value), str(value).lower())

    def data(self, index, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None
        s = self._summaries[index.row()]
        key = _KEYS[index.column()]
        if role == Qt.ItemDataRole.DisplayRole:
            return self._value(s, key)[0]
        if role == SORT_ROLE:
            return self._value(s, key)[1]
        if role == Qt.ItemDataRole.DecorationRole and key == 'dot':
            return self._dot(s.status)
        if role == Qt.ItemDataRole.ToolTipRole:
            if key == 'dot':
                return status_title(s.status)
            if key == 'note':
                return s.note(self._dcaps_mode) or None
        if role == Qt.ItemDataRole.ForegroundRole:
            if key == 'note' and s.status == ERROR:
                return QColor(self._colors['error_fg'])
            if s.status == NONE:
                return QColor(self._colors['text2'])
        if role == Qt.ItemDataRole.TextAlignmentRole and key in _NUMERIC:
            return Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        return None


class SnpFilterProxy(QSortFilterProxyModel):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._statuses: Set[str] = set()
        self._search = ''
        self.setSortRole(SORT_ROLE)

    def set_statuses(self, statuses: Set[str]) -> None:
        self._statuses = set(statuses)
        self.invalidateFilter()

    def set_search(self, text: str) -> None:
        self._search = text.strip().lower()
        self.invalidateFilter()

    def filterAcceptsRow(self, source_row, source_parent):
        s = self.sourceModel().summary(source_row)
        if self._statuses and s.status not in self._statuses:
            return False
        if not self._search:
            return True
        scheme = s.scheme or {}
        haystack = ' '.join([s.record.snp_id, s.label, scheme.get('enzyme') or '',
                             scheme.get('isoschizomers') or '']).lower()
        return self._search in haystack

    def lessThan(self, left, right):
        a, b = left.data(SORT_ROLE), right.data(SORT_ROLE)
        # Empty cells stay at the end in both directions.
        if a is None or b is None:
            if a is None and b is None:
                return False
            return (a is None) == (self.sortOrder() == Qt.SortOrder.DescendingOrder)
        return a < b


class StatusChips(QWidget):
    """Counters that filter the list: «Все N» and one per status."""
    toggled = pyqtSignal(set)

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)
        self._buttons: Dict[str, QPushButton] = {}
        self._colors = palette.colors(False)
        for key in ('all',) + STATUS_ORDER:
            button = QPushButton()
            button.setCheckable(True)
            button.clicked.connect(lambda checked, k=key: self._on_clicked(k, checked))
            layout.addWidget(button)
            self._buttons[key] = button
        layout.addStretch()
        self._buttons['all'].setChecked(True)
        self.set_counts(dict.fromkeys(STATUS_ORDER, 0), 0)

    def chip(self, status: str) -> QPushButton:
        return self._buttons[status]

    def selected(self) -> Set[str]:
        return {k for k in STATUS_ORDER if self._buttons[k].isChecked()}

    def _on_clicked(self, key: str, checked: bool):
        if key == 'all' or not self.selected():
            self.reset()
        else:
            self._buttons['all'].setChecked(False)
        self.toggled.emit(self.selected())

    def reset(self) -> None:
        for key in STATUS_ORDER:
            self._buttons[key].setChecked(False)
        self._buttons['all'].setChecked(True)

    def set_counts(self, counts: Dict[str, int], total: int) -> None:
        self._buttons['all'].setText(t('Все {count}', count=total))
        for key in STATUS_ORDER:
            self._buttons[key].setText(f'{status_title(key)} {counts.get(key, 0)}')
        self._buttons[ERROR].setVisible(counts.get(ERROR, 0) > 0)

    def apply_colors(self, colors: dict) -> None:
        self._colors = dict(colors)
        c = colors
        for key, button in self._buttons.items():
            if key == 'all':
                style = (f"QPushButton {{ border: 1px solid {c['line']}; background: {c['panel']}; }}"
                         f"QPushButton:checked {{ border-color: {c['accent']}; color: {c['accent']}; }}")
            else:
                style = (f"QPushButton {{ color: {c[f'{key}_fg']}; background: {c[f'{key}_bg']}; border: 1px solid "
                         f"{c[f'{key}_bg']}; }} QPushButton:checked {{ border-color: {c[f'{key}_fg']}; }}")
            button.setStyleSheet(style)


class SnpListView(QWidget):
    selection_changed = pyqtSignal(object)
    export_requested = pyqtSignal(str)

    def __init__(self, settings: AppSettings, parent=None):
        super().__init__(parent)
        self._settings = settings
        self._has_ids = False
        self._sort_column: Optional[int] = None
        self._sort_order = Qt.SortOrder.AscendingOrder
        stored = settings.value('visible_columns')
        self._visible = set(stored) if isinstance(stored, list) else set(DEFAULT_COLUMNS)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        top = QHBoxLayout()
        self.chips = StatusChips()
        top.addWidget(self.chips, 1)
        self.export_button = QToolButton()
        self.export_button.setText(t('Выгрузить в Excel'))
        self.export_button.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self.export_button.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.export_menu = QMenu(self.export_button)
        self._export_all = self.export_menu.addAction(t('Все SNP'), lambda: self.export_requested.emit('all'))
        self._export_shown = self.export_menu.addAction(t('Только показанные ({count})', count=0),
                                                        lambda: self.export_requested.emit('shown'))
        self.export_menu.addAction(t('CSV: все схемы'), lambda: self.export_requested.emit('csv'))
        self.export_button.setMenu(self.export_menu)
        top.addWidget(self.export_button)
        layout.addLayout(top)

        self.search = QLineEdit()
        self.search.setPlaceholderText(t('Поиск по ID, координате или ферменту'))
        self.search.setClearButtonEnabled(True)
        layout.addWidget(self.search)

        self.model = SnpListModel(self)
        self.proxy = SnpFilterProxy(self)
        self.proxy.setSourceModel(self.model)
        self.table = QTableView()
        self.table.setModel(self.proxy)
        self.table.setSelectionBehavior(QTableView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QTableView.SelectionMode.SingleSelection)
        self.table.setAlternatingRowColors(False)
        self.table.setShowGrid(False)
        self.table.setWordWrap(False)
        self.table.verticalHeader().hide()
        self.table.verticalHeader().setDefaultSectionSize(26)
        header = self.table.horizontalHeader()
        header.setMinimumHeight(28)
        header.setSectionsClickable(True)
        header.setSortIndicatorShown(False)
        header.setStretchLastSection(False)
        header.sectionClicked.connect(self._on_header_clicked)
        header.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        header.customContextMenuRequested.connect(self._columns_menu)
        for index, (key, _) in enumerate(COLUMNS):
            if key in _WIDTHS:
                self.table.setColumnWidth(index, _WIDTHS[key])
        header.setSectionResizeMode(self.column_index('note'), QHeaderView.ResizeMode.Stretch)
        self.table.selectionModel().currentRowChanged.connect(self._on_current_changed)
        layout.addWidget(self.table, 1)

        placeholder_box = QWidget()
        placeholder_layout = QVBoxLayout(placeholder_box)
        placeholder_layout.addStretch()
        self.placeholder = QLabel(t('Выберите файл SNP и геном слева и нажмите «Запустить»'))
        self.placeholder.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.placeholder.setWordWrap(True)
        self.placeholder.setProperty('role', 'muted')
        placeholder_layout.addWidget(self.placeholder)
        self.reset_filter_button = QPushButton(t('Сбросить фильтр'))
        self.reset_filter_button.setProperty('role', 'link')
        self.reset_filter_button.clicked.connect(self.reset_filter)
        placeholder_layout.addWidget(self.reset_filter_button, 0, Qt.AlignmentFlag.AlignHCenter)
        placeholder_layout.addStretch()
        self._placeholder_box = placeholder_box
        layout.addWidget(placeholder_box, 1)

        self.footer = QLabel()
        self.footer.setProperty('role', 'muted')
        layout.addWidget(self.footer)

        self.chips.toggled.connect(self._on_statuses)
        self.search.textChanged.connect(self._on_search)
        self.apply_colors(palette.colors(False))
        self._apply_column_visibility()
        self.set_placeholder(self.placeholder.text())

    # ----- columns -----
    def column_index(self, key: str) -> int:
        return _KEYS.index(key)

    def _apply_column_visibility(self):
        for index, key in enumerate(_KEYS):
            hidden = key not in self._visible and key != 'dot'
            if key == 'id' and not self._has_ids:
                hidden = True
            self.table.setColumnHidden(index, hidden)

    def set_column_visible(self, key: str, visible: bool) -> None:
        (self._visible.add if visible else self._visible.discard)(key)
        self._settings.set_value('visible_columns', [k for k in _KEYS if k in self._visible])
        self._apply_column_visibility()

    def _columns_menu(self, position):
        menu = QMenu(self)
        menu.addSection(t('Столбцы…'))
        for key, title in COLUMNS:
            if key == 'dot' or (key == 'id' and not self._has_ids):
                continue
            action = menu.addAction(t(title))
            action.setCheckable(True)
            action.setChecked(key in self._visible)
            action.toggled.connect(lambda checked, k=key: self.set_column_visible(k, checked))
        menu.exec(self.table.horizontalHeader().mapToGlobal(position))

    # ----- data -----
    def set_summaries(self, summaries: List[SnpSummary], dcaps_mode: str, has_ids: bool) -> None:
        self._has_ids = has_ids
        self.model.set_summaries(summaries, dcaps_mode)
        self._apply_column_visibility()
        self._update_counts()
        self._update_view()

    def shown_summaries(self) -> List[SnpSummary]:
        return [self.model.summary(self.proxy.mapToSource(self.proxy.index(r, 0)).row())
                for r in range(self.proxy.rowCount())]

    def current_summary(self) -> Optional[SnpSummary]:
        index = self.table.currentIndex()
        if not index.isValid():
            return None
        return self.model.summary(self.proxy.mapToSource(index).row())

    def refresh_current(self) -> None:
        index = self.table.currentIndex()
        rows = ([self.proxy.mapToSource(index).row()] if index.isValid()
                else range(self.model.rowCount()))
        for row in rows:
            self.model.refresh_row(row)
        self.proxy.invalidateFilter()
        self._update_counts()
        self._update_view()

    def set_placeholder(self, text: str) -> None:
        self.placeholder.setText(text)
        self.reset_filter_button.hide()
        self._show_placeholder(True)

    def _show_placeholder(self, shown: bool):
        self.table.setVisible(not shown)
        self._placeholder_box.setVisible(shown)
        self.placeholder.setVisible(shown)

    def reset_filter(self) -> None:
        self.chips.reset()
        self.search.clear()
        self.proxy.set_statuses(set())
        self._update_view()

    # ----- reactions -----
    def _on_statuses(self, statuses: set):
        self.proxy.set_statuses(statuses)
        self._update_view()

    def _on_search(self, text: str):
        self.proxy.set_search(text)
        self._update_view()

    def _on_header_clicked(self, column: int):
        if column != self._sort_column:
            self._sort_column, self._sort_order = column, Qt.SortOrder.AscendingOrder
        elif self._sort_order == Qt.SortOrder.AscendingOrder:
            self._sort_order = Qt.SortOrder.DescendingOrder
        else:
            self._sort_column = None
        header = self.table.horizontalHeader()
        if self._sort_column is None:
            self.proxy.sort(-1)
            header.setSortIndicatorShown(False)
        else:
            self.proxy.sort(self._sort_column, self._sort_order)
            header.setSortIndicatorShown(True)
            header.setSortIndicator(self._sort_column, self._sort_order)
        self._update_view()

    def _on_current_changed(self, current, _previous):
        self.selection_changed.emit(self.current_summary() if current.isValid() else None)

    def _update_counts(self):
        self.chips.set_counts(self.model.counts(), self.model.rowCount())

    def _update_view(self):
        total, shown = self.model.rowCount(), self.proxy.rowCount()
        self._export_shown.setText(t('Только показанные ({count})', count=shown))
        order = (t('как в файле') if self._sort_column is None
                 else t('по столбцу «{column}»', column=t(COLUMNS[self._sort_column][1])))
        self.footer.setText(t('Показано {shown} из {total} · порядок: {order}', shown=shown, total=total, order=order)
                            if total else '')
        if total and not shown:
            self.placeholder.setText(t('Нет SNP с такими условиями'))
            self.reset_filter_button.show()
            self._show_placeholder(True)
        elif total:
            self._show_placeholder(False)
        if not self.table.currentIndex().isValid():
            self.selection_changed.emit(None)

    def apply_colors(self, colors: dict) -> None:
        self.model.set_colors(colors)
        self.chips.apply_colors(colors)
        self.export_button.setIcon(icons.icon('download', colors['text']))
