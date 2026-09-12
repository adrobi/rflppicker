"""Typed, read-only results with sorting and filtering in a single proxy model."""
from PyQt6.QtCore import QAbstractTableModel, QModelIndex, Qt, QSortFilterProxyModel
from PyQt6.QtGui import QColor


LABELS = {
    'variant': 'Вариант', 'mapped_id': 'Контиг', 'enzyme': 'Фермент',
    'site': 'Сайт', 'pattern': 'Паттерн', 'frags_ref': 'Фрагменты REF, п.н.',
    'frags_alt': 'Фрагменты ALT, п.н.', 'diag_delta_bp': 'Различие, п.н.',
    'amplicon_start': 'Начало (1-based)', 'amplicon_end': 'Конец (1-based)',
    'amplicon_len': 'Длина, п.н.', 'snp_offset_in_amplicon': 'SNP (0-based)',
    'primers_left': 'Левый праймер 5′→3′', 'primers_right': 'Правый праймер 5′→3′',
    'tm_left': 'Tm лев., °C', 'tm_right': 'Tm прав., °C', 'product_size': 'Продукт, п.н.',
    'primer3_size': 'Продукт Primer3', 'primer_left_start': 'Начало L',
    'primer_left_len': 'Длина L', 'primer_right_start': 'Начало R',
    'primer_right_len': 'Длина R', 'suppliers': 'Поставщики', 'status': 'Статус',
    'reason': 'Пояснение', 'delta_threshold_bp': 'Порог Δ, п.н.',
    'line': 'Строка', 'chrom': 'Хромосома', 'pos': 'Позиция', 'ref': 'REF',
    'alt': 'ALT', 'details': 'Подробности', 'raw': 'Исходная строка',
}
STATUS_LABELS = {
    'ok': 'Праймеры подобраны', 'window_only': 'Поиск по окну',
    'no_primers': 'Праймеры не найдены', 'primer3_error': 'Ошибка Primer3',
    'no_enzyme_found': 'Фермент не найден', 'analysis_error': 'Ошибка анализа',
}
MAIN_COLUMNS = {'variant', 'enzyme', 'pattern', 'frags_ref', 'frags_alt',
                'diag_delta_bp', 'product_size', 'tm_left', 'tm_right', 'status'}


class ResultModel(QAbstractTableModel):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.header = []
        self.rows = []

    def replace(self, header, rows):
        self.beginResetModel()
        self.header = list(header)
        self.rows = [list(row) for row in rows]
        self.endResetModel()

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self.rows)

    def columnCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self.header)

    def record(self, row):
        return dict(zip(self.header, self.rows[row]))

    def data(self, index, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid() or index.row() >= len(self.rows):
            return None
        row = self.rows[index.row()]
        value = row[index.column()] if index.column() < len(row) else None
        name = self.header[index.column()]
        if role == Qt.ItemDataRole.UserRole:
            return value
        if role == Qt.ItemDataRole.DisplayRole:
            if value is None:
                return ''
            if name == 'status':
                return STATUS_LABELS.get(value, str(value))
            if isinstance(value, float):
                return f'{value:.3f}'.rstrip('0').rstrip('.')
            return str(value)
        if role == Qt.ItemDataRole.ToolTipRole:
            return f'{LABELS.get(name, name)} ({name})\n{value if value is not None else "—"}'
        if role == Qt.ItemDataRole.TextAlignmentRole and isinstance(value, (int, float)):
            return Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        if role == Qt.ItemDataRole.ForegroundRole and name == 'status':
            if value == 'ok':
                return QColor('#087d6a')
            if value in ('primer3_error', 'analysis_error'):
                return QColor('#bc3545')
            return QColor('#916116')
        return None

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):
        if role == Qt.ItemDataRole.DisplayRole:
            if orientation == Qt.Orientation.Horizontal:
                return LABELS.get(self.header[section], self.header[section])
            return str(section + 1)
        return None


class ResultFilter(QSortFilterProxyModel):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.enzyme = self.pattern = self.status = self.query = ''
        self.setSortRole(Qt.ItemDataRole.UserRole)
        self.setDynamicSortFilter(True)

    def set_filters(self, enzyme='', pattern='', status='', query=''):
        self.enzyme, self.pattern, self.status = enzyme, pattern, status
        self.query = query.casefold().strip()
        self.invalidateFilter()

    def filterAcceptsRow(self, row, parent):
        record = self.sourceModel().record(row)
        for key in ('enzyme', 'pattern', 'status'):
            target = getattr(self, key)
            if target and record.get(key) != target:
                return False
        return not self.query or self.query in ' '.join(str(v) for v in record.values()).casefold()

    def lessThan(self, left, right):
        a, b = left.data(Qt.ItemDataRole.UserRole), right.data(Qt.ItemDataRole.UserRole)
        if a is None or b is None:
            return a is None and b is not None
        if isinstance(a, (int, float)) and isinstance(b, (int, float)):
            return a < b
        return str(a).casefold() < str(b).casefold()
