"""The «Запуск» column: data files, the key parameters and the run button (spec §5)."""
import os
from typing import Dict, List, Optional, Tuple

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDoubleSpinBox, QFileDialog, QGridLayout, QHBoxLayout, QLabel, QLineEdit,
    QMenu, QPushButton, QScrollArea, QSizePolicy, QSpinBox, QToolButton, QVBoxLayout, QWidget,
)

import icons
import palette
from app_settings import AppSettings
from i18n import N_, core_text, t
from rflp_core import (
    DEFAULT_PARAMS, RECOMMENDED_PARAMS, _PRIMER3_AVAILABLE, _SUPPLIER_CANON, get_base_dir, validate_params,
)

KEY_PARAMS = ('prod_min', 'prod_max', 'tm_min', 'tm_max', 'commercial_only', 'suppliers', 'dcaps')
ADVANCED_PARAMS = ('min_frag', 'max_frag', 'delta', 'gel_resolution_bp', 'gel_resolution_pct', 'gel_stress_factor',
                   'flank', 'primer_pairs', 'top_results', 'max_cuts', 'group_isoschizomers', 'snp_end_distance',
                   'exclude_primer_sites', 'mask_repeats', 'exclude_background_sites', 'dcaps_mismatches',
                   'joint_design', 'use_primer3', 'gain_loss_only', 'skip_norm', 'no_save_norm')
PANEL_ERROR_KEYS = ('input_file', 'fasta', 'background_vcf', 'prod_min', 'tm_min')
SUPPLIER_TITLES = {'A': 'Agilent', 'E': 'Thermo Fisher', 'F': 'Promega', 'H': 'Minotech', 'I': 'SibEnzyme',
                   'J': 'Nippon Gene', 'K': 'Takara', 'M': 'Roche', 'N': 'NEB',
                   'R': 'Molecular Biology Resources', 'S': 'Sigma', 'T': 'Toyobo', 'V': 'Vivantis', 'Y': 'EURx',
                   'Z': 'SinaClon'}
DCAPS_CHOICES = ((N_('Нет'), 'off'), (N_('Если нет надёжной схемы'), 'fallback'), (N_('Для каждого SNP'), 'always'))
CHOOSE_FILE = N_('Выбрать файл…')
COLLAPSED_WIDTH = 36
_MISSING = {'input_file': N_('Выберите файл SNP'), 'fasta': N_('Выберите геном')}


def field_errors(params: dict) -> Dict[str, str]:
    """Field key → message in the interface language, from the messages of rflp_core.validate_params."""
    errors: Dict[str, str] = {}
    for message in validate_params(params):
        if message.startswith('Не задан параметр: '):
            key = message.split(': ', 1)[1]
            text = t(_MISSING[key]) if key in _MISSING else core_text(message)
        elif message.startswith('Файл не найден ('):
            key, text = message[len('Файл не найден ('):].split(')', 1)[0], t('Файл не найден')
        elif 'min_frag/max_frag' in message or 'размер фрагмента' in message:
            key, text = 'min_frag', core_text(message)
        elif 'prod_min/prod_max' in message or 'размер продукта' in message:
            key, text = 'prod_min', core_text(message)
        elif 'tm_min/tm_max' in message or ' Tm' in message:
            key, text = 'tm_min', core_text(message)
        else:
            key = next((k for k in params if message.startswith(f'{k}:')), 'gel_resolution_bp')
            text = core_text(message)
        errors.setdefault(key, text)
    return errors


def _label(text: str, role: str = '', name: str = '') -> QLabel:
    label = QLabel(text)
    if role:
        label.setProperty('role', role)
    if name:
        label.setObjectName(name)
    return label


def _error_label(key: str) -> QLabel:
    label = _label('', 'error', f'error_{key}')
    label.setWordWrap(True)
    label.hide()
    return label


class SupplierPicker(QToolButton):
    """One button with a checklist of suppliers instead of fifteen check boxes."""
    changed = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setObjectName('supplier_picker')
        self.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        menu = QMenu(self)
        self._actions = {}
        for code in sorted(_SUPPLIER_CANON):
            action = menu.addAction(SUPPLIER_TITLES.get(code, _SUPPLIER_CANON[code][0].capitalize()))
            action.setCheckable(True)
            action.toggled.connect(lambda checked, c=code: self._on_supplier(c, checked))
            self._actions[code] = action
        menu.addSeparator()
        self._all = menu.addAction(t('Включая непродаваемые'))
        self._all.setCheckable(True)
        self._all.toggled.connect(self._on_all)
        self.setMenu(menu)
        self._update_text()

    def _on_supplier(self, code: str, checked: bool):
        if checked and self._all.isChecked():
            self._all.blockSignals(True)
            self._all.setChecked(False)
            self._all.blockSignals(False)
        self._update_text()
        self.changed.emit()

    def _on_all(self, checked: bool):
        if checked:
            for action in self._actions.values():
                action.blockSignals(True)
                action.setChecked(False)
                action.blockSignals(False)
        self._update_text()
        self.changed.emit()

    def value(self) -> Tuple[bool, List[str]]:
        return not self._all.isChecked(), [c for c, a in self._actions.items() if a.isChecked()]

    def set_value(self, commercial_only: bool, suppliers: List[str]) -> None:
        wanted = set() if not commercial_only else {s.upper() for s in suppliers}
        for code, action in self._actions.items():
            action.blockSignals(True)
            action.setChecked(code in wanted)
            action.blockSignals(False)
        self._all.blockSignals(True)
        self._all.setChecked(not commercial_only)
        self._all.blockSignals(False)
        self._update_text()

    def _update_text(self):
        commercial_only, suppliers = self.value()
        if not commercial_only:
            text = t('Все, включая непродаваемые')
        elif not suppliers:
            text = t('Все продаваемые')
        else:
            names = [self._actions[c].text() for c in suppliers]
            text = ', '.join(names[:2])
            if len(names) > 2:
                text = t('{names} и ещё {count}', names=text, count=len(names) - 2)
        self.setText(text)
        self.setToolTip(text)


# key -> (label, kind, minimum, maximum, step/decimals, tooltip)
_ADVANCED_FIELDS = {
    'min_frag': (N_('Видимые фрагменты от, п.н.'), 'int', 10, 5000, 1,
                 N_('В совместном поиске — нижняя граница видимых полос. В поиске по одной паре/окну — '
                    'ограничение на каждый фрагмент.')),
    'max_frag': (N_('до, п.н.'), 'int', 10, 5000, 1,
                 N_('В совместном поиске — верхняя граница видимых полос. В поиске по одной паре/окну — '
                    'ограничение на каждый фрагмент.')),
    'delta': (N_('Порог Δ, п.н.'), 'int', 1, 1000, 1, N_('Наименьшее различие видимых полос REF и ALT.')),
    'gel_resolution_bp': (N_('Разрешение, п.н.'), 'float', 0.1, 1000, 2, ''),
    'gel_resolution_pct': (N_('Разрешение, %'), 'float', 0, 100, 2, ''),
    'gel_stress_factor': (N_('Ухудшение разрешения, ×'), 'float', 1, 5, 2,
                          N_('Модель полос без учёта яркости. Различие должно превышать максимум из порога в п.н. '
                             'и процента длины. Параметры требуют проверки на вашем геле.')),
    'flank': (N_('Фланк, п.н.'), 'int', 1, 5000, 1, N_('Число оснований по обе стороны SNP.')),
    'primer_pairs': (N_('Пар Primer3 на SNP'), 'int', 1, 100, 1, ''),
    'top_results': (N_('Схем на SNP'), 'int', 0, 10000, 1, N_('0 — все схемы.')),
    'max_cuts': (N_('Макс. разрезов'), 'int', 1, 10, 1, ''),
    'group_isoschizomers': (N_('Объединять изошизомеры'), 'bool', 0, 0, 0,
                            N_('Одна строка на сайт узнавания для SNP: лучшая схема, остальные ферменты с тем же '
                               'сайтом — в столбце «Изошизомеры».')),
    'snp_end_distance': (N_('SNP от края ампликона, п.н.'), 'int', 0, 1000, 1,
                         N_('Primer3 ставит SNP не ближе заданного числа п.н. к каждому концу ампликона, чтобы '
                            'диагностический фрагмент был видим. Значение должно быть не меньше мин. фрагмента. '
                            'Не действует в поиске по окну.')),
    'exclude_primer_sites': (N_('Исключать сайт фермента в праймере'), 'bool', 0, 0, 0,
                             N_('Сайт внутри праймера есть в продукте при любом аллеле, а разрез у самого конца ДНК '
                                'ненадёжен.')),
    'mask_repeats': (N_('Учитывать повторы (soft-mask)'), 'bool', 0, 0, 0,
                     N_('Строчные буквы FASTA (повторы в геномах NCBI/Ensembl) не допускаются на 3′-конце праймера.')),
    'exclude_background_sites': (N_('Исключать схемы с фоновыми вариантами'), 'bool', 0, 0, 0,
                                 N_('Схемы, где известные варианты из VCF меняют разрезы или длины полос.')),
    'dcaps_mismatches': (N_('Несовпадений dCAPS, до'), 'int', 1, 2, 1,
                         N_('Наибольшее число несовпадений в праймере dCAPS.')),
    'joint_design': (N_('Совместный подбор праймеров и фермента'), 'bool', 0, 0, 0,
                     N_('Отключите для сравнения с поиском по первой паре Primer3.')),
    'use_primer3': ('Primer3', 'bool', 0, 0, 0, N_('Без Primer3 анализируется окно вокруг SNP.')),
    'gain_loss_only': (N_('Только gain/loss'), 'bool', 0, 0, 0, ''),
    'skip_norm': (N_('Пропускать нормализацию'), 'bool', 0, 0, 0, N_('Отключает только перестановку REF/ALT.')),
    'no_save_norm': (N_('Не сохранять нормализацию'), 'bool', 0, 0, 0, ''),
}
_ADVANCED_GROUPS = (
    (N_('Гель'), ('min_frag', 'max_frag', 'delta', 'gel_resolution_bp', 'gel_resolution_pct', 'gel_stress_factor')),
    (N_('Поиск'), ('flank', 'primer_pairs', 'top_results', 'max_cuts')),
    (N_('Проверки'), ('group_isoschizomers', 'snp_end_distance', 'exclude_primer_sites', 'mask_repeats',
                  'exclude_background_sites', 'dcaps_mismatches')),
    (N_('Для сравнения и отладки'), ('joint_design', 'use_primer3', 'gain_loss_only', 'skip_norm', 'no_save_norm')),
)


class AdvancedDialog(QDialog):
    """Rarely changed parameters, grouped (spec §5.2)."""
    preset_requested = pyqtSignal(str)

    def __init__(self, params: dict, parent=None):
        super().__init__(parent)
        self.setWindowTitle(t('Дополнительные параметры'))
        self._widgets: Dict[str, QWidget] = {}
        self._errors: Dict[str, QLabel] = {}
        layout = QVBoxLayout(self)
        columns = QHBoxLayout()
        layout.addLayout(columns)
        for index, (title, keys) in enumerate(_ADVANCED_GROUPS):
            if index % 2 == 0:
                column = QVBoxLayout()
                columns.addLayout(column)
            column.addWidget(_label(t(title), 'caption'))
            grid = QGridLayout()
            column.addLayout(grid)
            for row, key in enumerate(keys):
                text, kind, low, high, step, tip = _ADVANCED_FIELDS[key]
                if kind == 'bool':
                    widget = QCheckBox(t(text))
                    grid.addWidget(widget, row * 2, 0, 1, 2)
                else:
                    widget = QSpinBox() if kind == 'int' else QDoubleSpinBox()
                    widget.setRange(low, high)
                    if kind == 'float':
                        widget.setDecimals(step)
                        widget.setSingleStep(0.1)
                    grid.addWidget(QLabel(t(text)), row * 2, 0)
                    grid.addWidget(widget, row * 2, 1)
                widget.setObjectName(f'adv_{key}')
                widget.setToolTip(t(tip) if tip else '')
                error = _label('', 'error', f'dialog_error_{key}')
                error.setWordWrap(True)
                error.hide()
                grid.addWidget(error, row * 2 + 1, 0, 1, 2)
                self._widgets[key] = widget
                self._errors[key] = error
            column.addStretch()
        self._widgets['top_results'].setSpecialValueText(t('Все'))
        self._widgets['snp_end_distance'].setSpecialValueText(t('Не ограничивать'))
        buttons = QHBoxLayout()
        for text, preset in ((t('Сбросить к рекомендуемым'), 'recommended'), (t('Как в v1.1.0'), 'v1.1')):
            button = QPushButton(text)
            button.clicked.connect(lambda _=False, p=preset: self.preset_requested.emit(p))
            buttons.addWidget(button)
        buttons.addStretch()
        done = QPushButton(t('Готово'))
        done.setProperty('role', 'primary')
        done.clicked.connect(self.accept)
        buttons.addWidget(done)
        layout.addLayout(buttons)
        self.set_params(params)

    def params(self) -> dict:
        values = {}
        for key, widget in self._widgets.items():
            values[key] = widget.isChecked() if isinstance(widget, QCheckBox) else widget.value()
        return values

    def set_params(self, params: dict) -> None:
        for key, widget in self._widgets.items():
            if key in params:
                if isinstance(widget, QCheckBox):
                    widget.setChecked(bool(params[key]))
                else:
                    widget.setValue(params[key])

    def show_errors(self, errors: Dict[str, str]) -> None:
        for key, label in self._errors.items():
            label.setText(errors.get(key, ''))
            label.setVisible(key in errors)


class LaunchPanel(QWidget):
    run_requested = pyqtSignal()
    stop_requested = pyqtSignal()
    collapse_toggled = pyqtSignal(bool)

    def __init__(self, settings: AppSettings, parent=None):
        super().__init__(parent)
        self._settings = settings
        self._colors = palette.colors(False)
        self._advanced = {k: RECOMMENDED_PARAMS[k] for k in ADVANCED_PARAMS}
        self._advanced_errors: Dict[str, str] = {}
        self._running = False
        self._collapsed = False
        self._names_path = ''

        outer = QVBoxLayout(self)
        outer.setContentsMargins(12, 12, 12, 12)
        # Full panel
        self._full = QWidget()
        full = QVBoxLayout(self._full)
        full.setContentsMargins(0, 0, 0, 0)
        header = QHBoxLayout()
        header.addWidget(_label(t('Запуск'), 'title'))
        header.addStretch()
        self._collapse_button = QToolButton()
        self._collapse_button.setObjectName('collapse_button')
        self._collapse_button.setToolTip(t('Свернуть панель'))
        self._collapse_button.clicked.connect(lambda: self._toggle_collapsed(True))
        header.addWidget(self._collapse_button)
        full.addLayout(header)

        self._form = QWidget()
        form = QGridLayout(self._form)
        form.setContentsMargins(0, 0, 4, 0)
        form.setColumnStretch(1, 1)
        row = 0

        def caption(text):
            nonlocal row
            form.addWidget(_label(text, 'caption'), row, 0, 1, 3)
            row += 1

        def field(text, widget, button=None, error_key=None):
            nonlocal row
            form.addWidget(QLabel(text), row, 0)
            form.addWidget(widget, row, 1, 1, 1 if button else 2)
            if button:
                form.addWidget(button, row, 2)
            row += 1
            if error_key:
                form.addWidget(self._errors[error_key], row, 1, 1, 2)
                row += 1

        self._errors = {key: _error_label(key) for key in PANEL_ERROR_KEYS}
        caption(t('Данные'))
        self.snp_edit = QLineEdit()
        self.snp_edit.setObjectName('snp_edit')
        self.snp_edit.setPlaceholderText(t('CSV или TSV: chrom, pos, ref, alt'))
        self.snp_edit.textChanged.connect(lambda: self._mark_missing(self.snp_edit))
        field('SNP', self.snp_edit, self._browse_button(self._browse_snp), 'input_file')
        self.genome_combo = QComboBox()
        self.genome_combo.setObjectName('genome_combo')
        self.genome_combo.activated.connect(self._on_genome_activated)
        field(t('Геном'), self.genome_combo, None, 'fasta')
        names_row = QHBoxLayout()
        self.names_label = _label('', 'muted', 'names_label')
        self.names_label.setWordWrap(True)
        names_row.addWidget(self.names_label, 1)
        names_change = QToolButton()
        names_change.setText(t('изменить'))
        names_change.setProperty('role', 'link')
        names_change.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        names_menu = QMenu(names_change)
        names_menu.addAction(t(CHOOSE_FILE), self._choose_names_file)
        names_menu.addAction(t('Без соответствий'), lambda: self._set_names(''))
        names_menu.addAction(t('Автоматически'), lambda: self._set_names(None))
        names_change.setMenu(names_menu)
        names_row.addWidget(names_change)
        form.addLayout(names_row, row, 1, 1, 2)
        row += 1
        self.vcf_edit = QLineEdit()
        self.vcf_edit.setObjectName('vcf_edit')
        self.vcf_edit.setPlaceholderText(t('не обязательно'))
        self.vcf_edit.setToolTip(t('VCF или VCF.gz известных вариантов породы или популяции.'))
        self.vcf_edit.textChanged.connect(lambda: self._mark_missing(self.vcf_edit))
        field('VCF', self.vcf_edit, self._browse_button(self._browse_vcf), 'background_vcf')

        caption(t('Праймеры'))
        self.prod_min_spin, self.prod_max_spin = self._range_spins('prod', QSpinBox, 50, 2000)
        field(t('Продукт, п.н.'), self._pair(self.prod_min_spin, self.prod_max_spin), None, 'prod_min')
        self.tm_min_spin, self.tm_max_spin = self._range_spins('tm', QDoubleSpinBox, 40.0, 80.0)
        field('Tm, °C', self._pair(self.tm_min_spin, self.tm_max_spin), None, 'tm_min')

        caption(t('Ферменты'))
        self.supplier_picker = SupplierPicker()
        field(t('Поставщики'), self.supplier_picker)
        self.dcaps_combo = QComboBox()
        self.dcaps_combo.setObjectName('dcaps_combo')
        for text, value in DCAPS_CHOICES:
            self.dcaps_combo.addItem(t(text), value)
        self.dcaps_combo.setToolTip(t('dCAPS: праймер с 1–2 несовпадениями у 3′-конца создаёт сайт фермента только '
                                      'в одном аллеле. Различие полос ≈ длине праймера.'))
        field('dCAPS', self.dcaps_combo)

        advanced_row = QHBoxLayout()
        self.advanced_link = QPushButton(t('Дополнительно…'))
        self.advanced_link.setObjectName('advanced_link')
        self.advanced_link.setProperty('role', 'link')
        self.advanced_link.clicked.connect(self.open_advanced)
        advanced_row.addWidget(self.advanced_link)
        self.changed_badge = _label('', '', 'changed_badge')
        self.changed_badge.hide()
        advanced_row.addWidget(self.changed_badge)
        advanced_row.addStretch()
        form.addLayout(advanced_row, row, 0, 1, 3)
        row += 1
        form.setRowStretch(row, 1)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setWidget(self._form)
        full.addWidget(scroll, 1)
        self.run_button = QPushButton(t('Запустить'))
        self.run_button.setObjectName('run_button')
        self.run_button.setProperty('role', 'primary')
        self.run_button.clicked.connect(self._on_run_clicked)
        full.addWidget(self.run_button)
        self.progress_label = _label('', 'muted', 'progress_label')
        self.progress_label.hide()
        full.addWidget(self.progress_label)
        outer.addWidget(self._full)

        # Collapsed strip
        self._strip = QWidget()
        strip = QVBoxLayout(self._strip)
        strip.setContentsMargins(0, 0, 0, 0)
        self._expand_button = QToolButton()
        self._expand_button.setToolTip(t('Развернуть панель «Запуск»'))
        self._expand_button.clicked.connect(lambda: self._toggle_collapsed(False))
        strip.addWidget(self._expand_button)
        vertical = QLabel('\n'.join(t('Запуск')))
        vertical.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        vertical.setProperty('role', 'muted')
        strip.addWidget(vertical)
        strip.addStretch()
        self._strip.hide()
        outer.addWidget(self._strip)

        self.supplier_picker.changed.connect(self._update_badge)
        self._load_settings()
        self.apply_colors(self._colors)

    # ----- building helpers -----
    def _browse_button(self, slot) -> QToolButton:
        button = QToolButton()
        button.setText('…')
        button.clicked.connect(slot)
        return button

    def _range_spins(self, name, cls, low, high):
        spins = []
        for side in ('min', 'max'):
            spin = cls()
            spin.setObjectName(f'{name}_{side}_spin')
            spin.setRange(low, high)
            if cls is QDoubleSpinBox:
                spin.setDecimals(1)
            spins.append(spin)
        return spins

    @staticmethod
    def _pair(first, second) -> QWidget:
        box = QWidget()
        layout = QHBoxLayout(box)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(first)
        layout.addWidget(QLabel('–'))
        layout.addWidget(second)
        return box

    # ----- genome and names -----
    def _genome_folder_files(self) -> List[str]:
        folder = os.path.join(get_base_dir(), 'genoms')
        if not os.path.isdir(folder):
            return []
        return [os.path.join(folder, f) for f in sorted(os.listdir(folder))
                if f.lower().endswith(('.fa', '.fasta', '.fna'))]

    def _fill_genomes(self, current: str = ''):
        self.genome_combo.blockSignals(True)
        self.genome_combo.clear()
        seen = set()
        for path in self._settings.recent_genomes() + self._genome_folder_files():
            key = os.path.normcase(path)
            if key in seen:
                continue
            seen.add(key)
            self.genome_combo.addItem(os.path.basename(path), path)
            self.genome_combo.setItemData(self.genome_combo.count() - 1, path, Qt.ItemDataRole.ToolTipRole)
        self.genome_combo.addItem(t(CHOOSE_FILE), None)
        index = self.genome_combo.findData(current) if current else (0 if self.genome_combo.count() > 1 else -1)
        self.genome_combo.setCurrentIndex(index if index >= 0 else -1)
        self.genome_combo.blockSignals(False)
        self._mark_missing(self.genome_combo)
        self.refresh_names()

    def _on_genome_activated(self, index: int):
        if self.genome_combo.itemData(index) is None:
            path, _ = QFileDialog.getOpenFileName(self, t('Геном FASTA'), '',
                                                  t('FASTA (*.fa *.fasta *.fna);;Все файлы (*)'))
            if path:
                self._settings.add_recent_genome(path)
            self._fill_genomes(path or self._previous_genome)
            return
        self._previous_genome = self.genome_path()
        self._mark_missing(self.genome_combo)
        self.refresh_names()

    def genome_path(self) -> str:
        return self.genome_combo.currentData() or ''

    def refresh_names(self) -> None:
        genome = self.genome_path()
        manual = self._settings.names_for_genome(genome) if genome else None
        if manual is not None:
            self._names_path = manual
            text = (t('Хромосомы: {name}', name=os.path.basename(manual)) if manual
                    else t('Хромосомы: как в файле SNP'))
        else:
            stem = os.path.splitext(os.path.basename(genome))[0]
            auto = os.path.join(get_base_dir(), 'names', f'{stem}.txt') if genome else ''
            self._names_path = auto if auto and os.path.isfile(auto) else ''
            text = (t('Хромосомы: {name} (авто)', name=os.path.basename(auto)) if self._names_path
                    else t('Хромосомы: как в файле SNP'))
        self.names_label.setText(text)
        self.names_label.setToolTip(self._names_path)

    def _set_names(self, names: Optional[str]):
        genome = self.genome_path()
        if not genome:
            return
        if names is None:
            self._settings.clear_names_for_genome(genome)
        else:
            self._settings.set_names_for_genome(genome, names)
        self.refresh_names()

    def _choose_names_file(self):
        path, _ = QFileDialog.getOpenFileName(self, t('Файл соответствий хромосом'),
                                              os.path.join(get_base_dir(), 'names'), t('Текст (*.txt);;Все файлы (*)'))
        if path:
            self._set_names(path)

    # ----- files -----
    def _browse_snp(self):
        start = self._settings.value('last_snp_dir', '') or ''
        path, _ = QFileDialog.getOpenFileName(self, t('Файл SNP'), start,
                                              t('Таблицы (*.tsv *.csv *.txt);;Все файлы (*)'))
        if path:
            self.snp_edit.setText(path)
            self._settings.set_value('last_snp_dir', os.path.dirname(path))

    def _browse_vcf(self):
        path, _ = QFileDialog.getOpenFileName(self, t('VCF известных вариантов'), '',
                                              t('VCF (*.vcf *.vcf.gz *.gz);;Все файлы (*)'))
        if path:
            self.vcf_edit.setText(path)

    def set_files(self, snp: str, genome: str, vcf: str = '') -> None:
        self.snp_edit.setText(snp)
        if genome:
            self._settings.add_recent_genome(genome)
        self._fill_genomes(genome)
        self.vcf_edit.setText(vcf)

    def _mark_missing(self, widget):
        path = widget.currentData() if isinstance(widget, QComboBox) else widget.text().strip()
        missing = bool(path) and not os.path.isfile(path)
        widget.setStyleSheet(f"color: {self._colors['error_fg']};" if missing else '')
        widget.setToolTip(t('Файл не найден') if missing else (path or ''))

    # ----- parameters -----
    def collect_params(self) -> dict:
        commercial_only, suppliers = self.supplier_picker.value()
        params = dict(self._advanced)
        params.update(
            input_file=self.snp_edit.text().strip(), fasta=self.genome_path(), names2_file=self._names_path,
            background_vcf=self.vcf_edit.text().strip(),
            prod_min=self.prod_min_spin.value(), prod_max=self.prod_max_spin.value(),
            tm_min=self.tm_min_spin.value(), tm_max=self.tm_max_spin.value(),
            commercial_only=commercial_only, suppliers=suppliers, dcaps=self.dcaps_combo.currentData(),
        )
        params['use_primer3'] = bool(params['use_primer3']) and _PRIMER3_AVAILABLE
        return params

    def set_params(self, params: dict) -> None:
        for key in ADVANCED_PARAMS:
            if key in params:
                self._advanced[key] = params[key]
        for key, spin in (('prod_min', self.prod_min_spin), ('prod_max', self.prod_max_spin),
                          ('tm_min', self.tm_min_spin), ('tm_max', self.tm_max_spin)):
            if key in params:
                spin.setValue(params[key])
        if 'commercial_only' in params or 'suppliers' in params:
            current = self.supplier_picker.value()
            self.supplier_picker.set_value(params.get('commercial_only', current[0]),
                                           params.get('suppliers', current[1]))
        if 'dcaps' in params:
            index = self.dcaps_combo.findData(params['dcaps'])
            if index >= 0:
                self.dcaps_combo.setCurrentIndex(index)
        self._advanced_errors = {}
        self._update_badge()

    def apply_preset(self, name: str) -> None:
        preset = {'recommended': RECOMMENDED_PARAMS, 'v1.1': DEFAULT_PARAMS}[name]
        self.set_params({k: v for k, v in preset.items() if k in KEY_PARAMS + ADVANCED_PARAMS})

    def changed_count(self) -> int:
        return sum(self._advanced[k] != RECOMMENDED_PARAMS[k] for k in ADVANCED_PARAMS)

    def _update_badge(self):
        c = self._colors
        if self._advanced_errors:
            text, fg, bg = t('ошибка'), c['error_fg'], c['error_bg']
        elif self.changed_count():
            text, fg, bg = t('изменено {count}', count=self.changed_count()), c['nominal_fg'], c['nominal_bg']
        else:
            self.changed_badge.hide()
            return
        self.changed_badge.setText(text)
        self.changed_badge.setStyleSheet(f'color: {fg}; background: {bg}; border-radius: 4px; padding: 1px 6px;')
        self.changed_badge.show()

    def validate(self) -> bool:
        errors = field_errors(self.collect_params())
        for key, label in self._errors.items():
            label.setText(errors.get(key, ''))
            label.setVisible(key in errors)
        self._advanced_errors = {k: v for k, v in errors.items() if k not in self._errors}
        self._update_badge()
        return not errors

    def open_advanced(self) -> None:
        dialog = AdvancedDialog(self._advanced, self)
        dialog.show_errors(self._advanced_errors)

        def preset(name):
            self.apply_preset(name)
            dialog.set_params(self._advanced)
            dialog.show_errors({})

        dialog.preset_requested.connect(preset)
        if dialog.exec():
            self.set_params(dialog.params())

    # ----- running and collapsing -----
    def _on_run_clicked(self):
        (self.stop_requested if self._running else self.run_requested).emit()

    def set_running(self, running: bool) -> None:
        self._running = running
        self._form.setEnabled(not running)
        self.run_button.setText(t('Остановить') if running else t('Запустить'))
        self.run_button.setIcon(icons.icon('player-stop' if running else 'player-play', self._colors['on_accent']))
        self.progress_label.setVisible(running)
        if running:
            self.progress_label.setText(t('Подготовка…'))

    def set_progress(self, done: int, total: int) -> None:
        self.progress_label.setText(t('{done} из {total} SNP', done=done, total=total))

    def _toggle_collapsed(self, collapsed: bool):
        self.set_collapsed(collapsed)
        self.collapse_toggled.emit(collapsed)

    def set_collapsed(self, collapsed: bool) -> None:
        self._collapsed = collapsed
        self._full.setVisible(not collapsed)
        self._strip.setVisible(collapsed)
        self.layout().setContentsMargins(*(4, 12, 4, 12) if collapsed else (12, 12, 12, 12))
        self.setMaximumWidth(COLLAPSED_WIDTH if collapsed else 16777215)
        self.setMinimumWidth(COLLAPSED_WIDTH if collapsed else 260)

    def is_collapsed(self) -> bool:
        return self._collapsed

    # ----- settings and colors -----
    def _load_settings(self):
        params = self._settings.params()
        self.set_params(params)
        self.snp_edit.setText(self._settings.value('last_snp_file', '') or '')
        self.vcf_edit.setText(self._settings.value('vcf', '') or '')
        self._fill_genomes()
        self._previous_genome = self.genome_path()

    def save_settings(self) -> None:
        params = self.collect_params()
        self._settings.set_params(params)
        self._settings.set_value('last_snp_file', params['input_file'])
        self._settings.set_value('vcf', params['background_vcf'])
        if params['fasta']:
            self._settings.add_recent_genome(params['fasta'])

    def apply_colors(self, colors: dict) -> None:
        self._colors = dict(colors)
        self._collapse_button.setIcon(icons.icon('chevrons-left', colors['text2']))
        self._expand_button.setIcon(icons.icon('chevrons-right', colors['text2']))
        self.run_button.setIcon(icons.icon('player-stop' if self._running else 'player-play', colors['on_accent']))
        for widget in (self.snp_edit, self.vcf_edit, self.genome_combo):
            self._mark_missing(widget)
        self._update_badge()
