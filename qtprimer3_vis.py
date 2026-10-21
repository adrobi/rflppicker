import sys
import os
import traceback
from typing import List, Dict, Tuple, Optional
from datetime import datetime

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QPushButton, QLabel, QFileDialog, QTableWidget, QTableWidgetItem, QLineEdit, QCheckBox,
    QSpinBox, QDoubleSpinBox, QTextEdit, QProgressBar, QMessageBox, QComboBox, QHeaderView,
    QGroupBox, QScrollArea
)
from PyQt6.QtCore import QThread, pyqtSignal, Qt
from PyQt6 import QtGui  # QColor
from PyQt6.QtGui import QIcon
import resources_rc

from pyfaidx import Fasta

import rflp_core as core_impl
from rflp_core import _PRIMER3_AVAILABLE, _PRIMER3_IMPORT_ERROR, _SUPPLIER_CANON, get_base_dir, logger
from reporting import export_csv as safe_export_csv, export_excel as safe_export_excel
from result_model import LABELS
from visualization import RFLPVisualizationDialog


class NumericTableItem(QTableWidgetItem):
    """Sort numeric results numerically while retaining human-readable text."""
    def __lt__(self, other):
        left = self.text().strip()
        right = other.text().strip() if other is not None else ""
        raw_left = self.data(Qt.ItemDataRole.UserRole)
        raw_right = other.data(Qt.ItemDataRole.UserRole) if other is not None else None
        if isinstance(raw_left, (int, float)) and isinstance(raw_right, (int, float)):
            return raw_left < raw_right
        try:
            return float(left) < float(right)
        except (ValueError, TypeError):
            return left.casefold() < right.casefold()


class RFLPCalcThread(QThread):
    finished = pyqtSignal(list, list)
    status = pyqtSignal(str)
    progress = pyqtSignal(int, int)

    def __init__(self, params):
        super().__init__()
        self.params = params
        self.report = {}

    def run(self):
        try:
            header, out_rows = core_impl.run_rflp_gui_mode(
                self.params, status_cb=self.status.emit,
                progress_cb=self.progress.emit,
                cancel_cb=self.isInterruptionRequested,
                report=self.report,
            )
            self.finished.emit(header, out_rows)
        except core_impl.AnalysisCancelled:
            self.status.emit("\n[INFO] Анализ отменён пользователем.")
            self.finished.emit([], [])
        except Exception:
            tb = traceback.format_exc()
            self.status.emit(f"\n[ERROR] Exception in calculation thread:\n{tb}")
            logger.error(tb)
            self.finished.emit([], [])


def _safe_int(s: object) -> Optional[int]:
    try:
        if s is None:
            return None
        t = str(s).strip()
        if t == "" or t.lower() == "none":
            return None
        return int(float(t))
    except Exception:
        return None


class RFLPPickerGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("RFLP Picker 1.1 — совместный подбор")

        self.setWindowIcon(QIcon(":/icons/app_icon.ico"))

        self.resize(1400, 800)

        central = QWidget()
        main_layout = QVBoxLayout()

        # ===== верх: параметры + лог =====
        top_layout = QHBoxLayout()

        # левая панель
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)

        left_layout.addWidget(QLabel("Файл вариантов (CSV/TSV):"))
        self.input_file = QLineEdit()
        left_layout.addWidget(self.input_file)
        btn_input = QPushButton("Выбрать файл вариантов")
        btn_input.clicked.connect(lambda: self.select_file(self.input_file))
        left_layout.addWidget(btn_input)

        left_layout.addWidget(QLabel("FASTA (геном, папка ./genoms рядом с exe):"))
        self.fasta_combo = QComboBox()
        self.fasta_combo.setEditable(True)
        self.fasta_combo.setToolTip("Можно выбрать локальный файл из genoms или ввести полный путь к FASTA.")
        left_layout.addWidget(self.fasta_combo)

        left_layout.addWidget(QLabel("Файл соответствий (names, папка ./names рядом с exe):"))
        self.names_combo = QComboBox()
        self.names_combo.setEditable(True)
        self.names_combo.setToolTip("Можно выбрать файл из names или ввести полный путь; поле можно оставить пустым.")
        left_layout.addWidget(self.names_combo)

        self.fasta_combo.currentIndexChanged.connect(self.update_names_selection)

        grid = QGridLayout()
        left_layout.addLayout(grid)

        self.flank_spin = QSpinBox()
        self.flank_spin.setRange(1, 5000)
        self.flank_spin.setValue(300)

        self.min_frag_spin = QSpinBox()
        self.min_frag_spin.setRange(10, 5000)
        self.min_frag_spin.setValue(80)

        self.max_frag_spin = QSpinBox()
        self.max_frag_spin.setRange(10, 5000)
        self.max_frag_spin.setValue(800)

        self.delta_spin = QSpinBox()
        self.delta_spin.setRange(1, 1000)
        self.delta_spin.setValue(25)

        self.max_cuts_spin = QSpinBox()
        self.max_cuts_spin.setRange(1, 10)
        self.max_cuts_spin.setValue(3)

        self.gloss_chk = QCheckBox("Только gain/loss")
        self.gloss_chk.setChecked(False)

        self.prod_min_spin = QSpinBox()
        self.prod_min_spin.setRange(50, 2000)
        self.prod_min_spin.setValue(250)

        self.prod_max_spin = QSpinBox()
        self.prod_max_spin.setRange(50, 2000)
        self.prod_max_spin.setValue(600)

        self.tm_min_spin = QDoubleSpinBox()
        self.tm_min_spin.setRange(40.0, 80.0)
        self.tm_min_spin.setDecimals(1)
        self.tm_min_spin.setValue(58.0)

        self.tm_max_spin = QDoubleSpinBox()
        self.tm_max_spin.setRange(40.0, 80.0)
        self.tm_max_spin.setDecimals(1)
        self.tm_max_spin.setValue(62.0)

        self.use_primer3_chk = QCheckBox("Использовать primer3")
        self.use_primer3_chk.setChecked(_PRIMER3_AVAILABLE)

        params_labels = [
            QLabel("Фланк (bp):"),
            QLabel("Мин. фрагмент (bp):"),
            QLabel("Макс. фрагмент (bp):"),
            QLabel("Диагностическая дельта (bp):"),
            QLabel("Макс. разрезов в ампликоне:"),
            QLabel("Опция:"),
            QLabel("Мин. размер продукта (bp):"),
            QLabel("Макс. размер продукта (bp):"),
            QLabel("Мин. Tm (°C):"),
            QLabel("Макс. Tm (°C):"),
            QLabel("Primer3:"),
        ]

        params_widgets = [
            self.flank_spin,
            self.min_frag_spin,
            self.max_frag_spin,
            self.delta_spin,
            self.max_cuts_spin,
            self.gloss_chk,
            self.prod_min_spin,
            self.prod_max_spin,
            self.tm_min_spin,
            self.tm_max_spin,
            self.use_primer3_chk,
        ]

        supplier_widgets: List[QCheckBox] = []
        supplier_codes = sorted(_SUPPLIER_CANON.keys())
        for code in supplier_codes:
            name = _SUPPLIER_CANON[code][0].capitalize()
            cb = QCheckBox(f"{code}: {name}")
            supplier_widgets.append(cb)

        max_len = max(len(params_widgets), len(supplier_widgets))
        for i in range(max_len):
            if i < len(params_widgets):
                grid.addWidget(params_labels[i], i, 0)
                grid.addWidget(params_widgets[i], i, 1)
            if i < len(supplier_widgets):
                col = 2 if i < (len(supplier_widgets) + 1) // 2 else 3
                row_sup = i if col == 2 else i - (len(supplier_widgets) + 1) // 2
                grid.addWidget(supplier_widgets[i], row_sup, col)

        self.supplier_checks = supplier_widgets

        recommended = core_impl.RECOMMENDED_PARAMS
        checks_box = QGroupBox('Отбор и проверки')
        self.checks_grid = QGridLayout(checks_box)
        self.commercial_chk = QCheckBox('Только продаваемые ферменты')
        self.commercial_chk.setChecked(recommended['commercial_only'])
        self.commercial_chk.setToolTip('Ферменты, для которых в REBASE указан хотя бы один поставщик.')
        self.group_iso_chk = QCheckBox('Объединять изошизомеры')
        self.group_iso_chk.setChecked(recommended['group_isoschizomers'])
        self.group_iso_chk.setToolTip('Одна строка на сайт узнавания для SNP: лучшая схема, '
                                      'остальные ферменты с тем же сайтом — в столбце «Изошизомеры».')
        self.checks_grid.addWidget(self.commercial_chk, 0, 0)
        self.checks_grid.addWidget(self.group_iso_chk, 0, 1)
        self.end_distance_spin = QSpinBox()
        self.end_distance_spin.setRange(0, 1000)
        self.end_distance_spin.setValue(recommended['snp_end_distance'])
        self.end_distance_spin.setSpecialValueText('Не ограничивать')
        self.end_distance_spin.setToolTip('Primer3 ставит SNP не ближе заданного числа п.н. к каждому концу '
                                          'ампликона, чтобы диагностический фрагмент был видим. Значение '
                                          'должно быть не меньше мин. фрагмента. Не действует в поиске по окну.')
        self.checks_grid.addWidget(QLabel('SNP не ближе к краю ампликона, п.н.:'), 1, 0)
        self.checks_grid.addWidget(self.end_distance_spin, 1, 1)
        self.primer_site_chk = QCheckBox('Исключать сайт фермента в праймере')
        self.primer_site_chk.setChecked(recommended['exclude_primer_sites'])
        self.primer_site_chk.setToolTip('Сайт внутри праймера есть в продукте при любом аллеле, а разрез '
                                        'у самого конца ДНК ненадёжен. Найденные сайты показаны в столбце '
                                        '«Сайт в праймере».')
        self.checks_grid.addWidget(self.primer_site_chk, 2, 0, 1, 2)
        self.mask_chk = QCheckBox('Учитывать повторы (soft-mask FASTA)')
        self.mask_chk.setChecked(recommended['mask_repeats'])
        self.mask_chk.setToolTip('Строчные буквы FASTA (повторы в геномах NCBI/Ensembl) не допускаются '
                                 'на 3′-конце праймера. Доля повторов в ампликоне — столбец «Повторы, %».')
        self.checks_grid.addWidget(self.mask_chk, 3, 0, 1, 2)
        left_layout.addWidget(checks_box)

        design_box = QGroupBox('Совместный подбор и различимость генотипов')
        design_grid = QGridLayout(design_box)
        self.joint_design_chk = QCheckBox('Подбирать праймеры и фермент совместно')
        self.joint_design_chk.setChecked(True)
        self.joint_design_chk.setToolTip('Отключите для сравнения с поиском по первой паре Primer3.')
        design_grid.addWidget(self.joint_design_chk, 0, 0, 1, 4)
        self.primer_pairs_spin = QSpinBox()
        self.primer_pairs_spin.setRange(1, 100)
        self.primer_pairs_spin.setValue(10)
        self.top_results_spin = QSpinBox()
        self.top_results_spin.setRange(0, 10000)
        self.top_results_spin.setValue(10)
        self.top_results_spin.setSpecialValueText('Все')
        self.gel_bp_spin = QDoubleSpinBox()
        self.gel_bp_spin.setRange(0.1, 1000)
        self.gel_bp_spin.setValue(10)
        self.gel_pct_spin = QDoubleSpinBox()
        self.gel_pct_spin.setRange(0, 100)
        self.gel_pct_spin.setValue(3)
        self.gel_stress_spin = QDoubleSpinBox()
        self.gel_stress_spin.setRange(1, 5)
        self.gel_stress_spin.setSingleStep(0.1)
        self.gel_stress_spin.setValue(1.5)
        controls = [('Пар Primer3 на SNP:', self.primer_pairs_spin),
                    ('Результатов на SNP:', self.top_results_spin),
                    ('Разрешение, п.н.:', self.gel_bp_spin),
                    ('Разрешение, %:', self.gel_pct_spin),
                    ('Ухудшение разрешения, ×:', self.gel_stress_spin)]
        for i, (label, widget) in enumerate(controls):
            row, col = 1 + i // 2, (i % 2) * 2
            design_grid.addWidget(QLabel(label), row, col)
            design_grid.addWidget(widget, row, col + 1)
        note = QLabel('Модель полос без учёта яркости. Различие должно превышать максимум '
                      'из порога в п.н. и процента длины. Параметры требуют проверки на вашем геле.')
        note.setWordWrap(True)
        design_grid.addWidget(note, 4, 0, 1, 4)
        self.min_frag_spin.setToolTip('В совместном поиске — нижняя граница видимых полос. '
                                     'В поиске по одной паре/окну — ограничение на каждый фрагмент.')
        self.max_frag_spin.setToolTip('В совместном поиске — верхняя граница видимых полос. '
                                     'В поиске по одной паре/окну — ограничение на каждый фрагмент.')
        self.use_primer3_chk.toggled.connect(self._update_design_controls)
        self.joint_design_chk.toggled.connect(self._update_design_controls)
        self._update_design_controls()

        btn_row = QHBoxLayout()
        self.run_btn = QPushButton("Запуск анализа")
        self.cancel_btn = QPushButton("Отмена")
        self.cancel_btn.setEnabled(False)
        self.export_csv_btn = QPushButton("Выгрузить в CSV")
        self.export_xlsx_btn = QPushButton("Выгрузить в Excel")
        btn_row.addWidget(self.run_btn)
        btn_row.addWidget(self.cancel_btn)
        btn_row.addWidget(self.export_csv_btn)
        btn_row.addWidget(self.export_xlsx_btn)

        self.progress = QProgressBar()

        # правая часть: лог
        self.status_log = QTextEdit()
        self.status_log.setReadOnly(True)
        self.status_log.setMaximumWidth(380)
        self.status_log.setMinimumWidth(280)
        self.status_log.setMaximumHeight(500)
        self.status_log.setMinimumHeight(120)

        log_norm_layout = QVBoxLayout()

        # ЧЕКБОКСЫ НОРМАЛИЗАЦИИ
        norm_layout = QHBoxLayout()
        self.skip_norm_chk = QCheckBox("Пропускать нормализацию")
        self.skip_norm_chk.setChecked(False)
        self.no_save_norm_chk = QCheckBox("Не сохранять нормализацию")
        self.no_save_norm_chk.setChecked(True)
        norm_layout.addWidget(self.skip_norm_chk)
        norm_layout.addWidget(self.no_save_norm_chk)

        log_norm_layout.addWidget(self.status_log)
        log_norm_layout.addLayout(norm_layout)
        log_norm_layout.addWidget(design_box)

        settings_scroll = QScrollArea()
        settings_scroll.setWidgetResizable(True)
        settings_scroll.setWidget(left_panel)
        left_column = QVBoxLayout()
        left_column.addWidget(settings_scroll)
        left_column.addLayout(btn_row)
        left_column.addWidget(self.progress)
        top_layout.addLayout(left_column, stretch=4)
        top_layout.addLayout(log_norm_layout, stretch=2)

        # ===== низ: фильтры + таблица + визуализация праймеров =====
        bottom_layout = QVBoxLayout()
        self.assay_summary = QLabel('Выберите результат: здесь появятся оценка и полосы трёх генотипов.')
        self.assay_summary.setWordWrap(True)
        self.assay_summary.setTextFormat(Qt.TextFormat.PlainText)
        bottom_layout.addWidget(self.assay_summary)

        filter_layout = QHBoxLayout()
        filter_layout.addWidget(QLabel("Фермент:"))
        self.filter_enzyme_combo = QComboBox()
        self.filter_enzyme_combo.addItem("Все ферменты")
        filter_layout.addWidget(self.filter_enzyme_combo)

        filter_layout.addWidget(QLabel("Паттерн:"))
        self.filter_pattern_combo = QComboBox()
        self.filter_pattern_combo.addItem("Все паттерны")
        filter_layout.addWidget(self.filter_pattern_combo)

        self.open_vis_btn = QPushButton("Открыть визуализацию")
        self.open_vis_btn.setEnabled(False)
        filter_layout.addWidget(self.open_vis_btn)

        filter_layout.addStretch(1)
        bottom_layout.addLayout(filter_layout)

        # Таблица результатов
        self.result_table = QTableWidget()
        self.result_table.setSortingEnabled(True)
        header = self.result_table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        header.setStretchLastSection(True)

        bottom_layout.addWidget(self.result_table)

        main_layout.addLayout(top_layout, stretch=4)
        main_layout.addLayout(bottom_layout, stretch=3)

        central.setLayout(main_layout)
        self.setCentralWidget(central)

        # сигналы
        self.run_btn.clicked.connect(self.run_calc)
        self.cancel_btn.clicked.connect(self.cancel_calc)
        self.export_csv_btn.clicked.connect(self.export_csv)
        self.export_xlsx_btn.clicked.connect(self.export_excel)
        self.filter_enzyme_combo.currentIndexChanged.connect(self.apply_filters)
        self.filter_pattern_combo.currentIndexChanged.connect(self.apply_filters)
        self.result_table.itemSelectionChanged.connect(self.on_table_selection_changed)
        self.open_vis_btn.clicked.connect(self.open_visualization_window)

        self.load_genomes_to_combo()
        self.load_names_to_combo()

        # состояние для визуализации; окно визуализации создаётся по кнопке
        self.last_params: dict = {}
        self.run_report: dict = {}
        self.fa_view: Optional[Fasta] = None
        self.col_index: Dict[str, int] = {}
        self._current_vis_payload: Optional[dict] = None
        self._current_vis_error: str = ""
        self.vis_dialog: Optional[RFLPVisualizationDialog] = None

    def load_genomes_to_combo(self):
        base_dir = get_base_dir()
        folder = os.path.join(base_dir, "genoms")
        self.fasta_combo.clear()
        if not os.path.isdir(folder):
            self.log_msg(f"Папка с геномами не найдена: {folder}", "error")
            return
        files = [f for f in os.listdir(folder) if f.endswith((".fa", ".fasta"))]
        if not files:
            self.log_msg(f"FASTA файлы не найдены в {folder}", "error")
            return
        for f in files:
            self.fasta_combo.addItem(f)
        self.log_msg(f"Найдено {len(files)} FASTA файлов в папке genoms")

    def load_names_to_combo(self):
        base_dir = get_base_dir()
        folder = os.path.join(base_dir, "names")
        self.names_combo.clear()
        if not os.path.isdir(folder):
            self.log_msg(f"Папка с names не найдена: {folder}", "error")
            return
        files = [f for f in os.listdir(folder) if f.endswith(".txt")]
        if not files:
            self.log_msg(f"Файлы names не найдены в {folder}", "error")
        for f in files:
            self.names_combo.addItem(f)
        self.log_msg(f"Загружено {len(files)} файлов names")

    def update_names_selection(self, index):
        fasta_name = self.fasta_combo.currentText()
        base_name = os.path.splitext(fasta_name)[0]
        found = False
        for i in range(self.names_combo.count()):
            names_name = self.names_combo.itemText(i)
            if os.path.splitext(names_name)[0] == base_name:
                self.names_combo.setCurrentIndex(i)
                self.log_msg(
                    f"Выбран файл соответствий: {names_name} для генома {fasta_name}"
                )
                found = True
                break
        if not found and self.names_combo.count() > 0:
            self.log_msg(
                f"Файл соответствий для {fasta_name} не найден среди файлов names",
                "warning",
            )

    def select_file(self, line_edit):
        try:
            path, _ = QFileDialog.getOpenFileName(
                self,
                "Выбрать файл",
                "",
                "Текстовые файлы (*.txt *.tsv *.csv);;Все файлы (*)",
            )
            if path:
                line_edit.setText(path)
                self.log_msg(f"Выбран файл: {path}")
        except Exception:
            tb = traceback.format_exc()
            self.log_msg(f"[ERROR] Ошибка при выборе файла:\n{tb}", "error")

    def _update_design_controls(self):
        self.joint_design_chk.setEnabled(self.use_primer3_chk.isChecked())
        joint = self.use_primer3_chk.isChecked() and self.joint_design_chk.isChecked()
        self.primer_pairs_spin.setEnabled(joint)
        self.top_results_spin.setEnabled(joint)

    def collect_params(self) -> dict:
        base_dir = get_base_dir()
        fasta_selected = self.fasta_combo.currentText()
        names_selected = self.names_combo.currentText()
        return {
            "flank": int(self.flank_spin.value()),
            "min_frag": int(self.min_frag_spin.value()),
            "max_frag": int(self.max_frag_spin.value()),
            "delta": int(self.delta_spin.value()),
            "gain_loss_only": self.gloss_chk.isChecked(),
            "max_cuts": int(self.max_cuts_spin.value()),
            "suppliers": [
                cb.text().split(":", 1)[1].strip()
                for cb in self.supplier_checks
                if cb.isChecked()
            ],
            "fasta": os.path.join(base_dir, "genoms", fasta_selected),
            "input_file": self.input_file.text(),
            "names2_file": (os.path.join(base_dir, "names", names_selected) if names_selected else ""),
            "use_primer3": self.use_primer3_chk.isChecked() and _PRIMER3_AVAILABLE,
            "prod_min": int(self.prod_min_spin.value()),
            "prod_max": int(self.prod_max_spin.value()),
            "tm_min": float(self.tm_min_spin.value()),
            "tm_max": float(self.tm_max_spin.value()),
            "joint_design": self.joint_design_chk.isChecked(),
            "primer_pairs": self.primer_pairs_spin.value(),
            "top_results": self.top_results_spin.value(),
            "gel_resolution_bp": self.gel_bp_spin.value(),
            "gel_resolution_pct": self.gel_pct_spin.value(),
            "gel_stress_factor": self.gel_stress_spin.value(),
            "skip_norm": self.skip_norm_chk.isChecked(),
            "no_save_norm": self.no_save_norm_chk.isChecked(),
            "commercial_only": self.commercial_chk.isChecked(),
            "group_isoschizomers": self.group_iso_chk.isChecked(),
            "snp_end_distance": self.end_distance_spin.value(),
            "exclude_primer_sites": self.primer_site_chk.isChecked(),
            "mask_repeats": self.mask_chk.isChecked(),
        }

    def run_calc(self):
        try:
            if not self.fasta_combo.currentText():
                QMessageBox.warning(self, "Ошибка", "Не выбран FASTA файл")
                return
            params = self.collect_params()
            errors = core_impl.validate_params(params)
            if errors:
                QMessageBox.warning(self, "Проверьте параметры", "\n".join(f"• {e}" for e in errors))
                return

            # сохраняем параметры последнего запуска (нужны для визуализации)
            self.last_params = dict(params)
            self.fa_view = None
            self.col_index = {}
            if getattr(self, "vis_dialog", None) is not None:
                self.vis_dialog.panel.show_empty()
            for k in ("flank", "min_frag", "max_frag", "fasta", "input_file"):
                if not params[k]:
                    QMessageBox.warning(self, "Ошибка", f"Не заполнено обязательное поле: {k}")
                    return

            self.status_log.clear()
            self.progress.setRange(0, 0)
            self.run_btn.setDisabled(True)
            self.cancel_btn.setEnabled(True)

            self.worker = RFLPCalcThread(params)
            self.worker.status.connect(self.log_status)
            self.worker.progress.connect(self._on_progress)
            self.worker.finished.connect(self.show_results)
            self.worker.start()
        except Exception:
            tb = traceback.format_exc()
            self.log_msg(f"[ERROR] Exception in run_calc:\n{tb}", "error")

    def cancel_calc(self):
        if getattr(self, "worker", None) is not None and self.worker.isRunning():
            self.worker.requestInterruption()
            self.cancel_btn.setEnabled(False)
            self.log_msg("[INFO] Запрошена отмена анализа…")

    def _on_progress(self, done: int, total: int):
        self.progress.setRange(0, max(1, total))
        self.progress.setValue(done)
        self.progress.setFormat(f"{done} / {total}")

    def log_status(self, msg):
        self.log_msg(msg)

    def log_msg(self, msg, level="info"):
        if level == "info":
            logger.info(msg)
        elif level == "warning":
            logger.warning(msg)
        elif level == "error":
            logger.error(msg)
        elif level == "debug":
            logger.debug(msg)
        self.status_log.append(msg)
        QApplication.processEvents()

    def show_results(self, header, rows):
        try:
            self.assay_summary.setText('Выберите результат: оценка > 1 означает различимость в модели полос.')
            self.progress.setRange(0, 1)
            self.run_btn.setEnabled(True)
            self.cancel_btn.setEnabled(False)
            self.out_header = header
            self.out_rows = rows
            if getattr(self, "worker", None) is not None:
                self.run_report = dict(getattr(self.worker, "report", {}) or {})

            # индекс колонок по имени (нужен для визуализации)
            self.col_index = {str(name): i for i, name in enumerate(header)}

            # подготавливаем FASTA для визуализации (ленивое открытие)
            self.fa_view = None
            fasta_path = (self.last_params or {}).get("fasta")
            if fasta_path:
                try:
                    self.fa_view = Fasta(fasta_path, as_raw=True, sequence_always_upper=True)
                except TypeError:
                    # fallback для старых версий pyfaidx
                    self.fa_view = Fasta(fasta_path, as_raw=True)
                except Exception as e:
                    self.log_msg(f"[WARN] Не удалось открыть FASTA для визуализации: {fasta_path}: {e}", "warning")
                    self.fa_view = None

            if getattr(self, "vis_dialog", None) is not None:
                self.vis_dialog.panel.show_empty()
            self.result_table.clear()
            self.result_table.setRowCount(0)
            self.result_table.setColumnCount(0)
            self.filter_enzyme_combo.clear()
            self.filter_enzyme_combo.addItem("Все ферменты")
            self.filter_pattern_combo.clear()
            self.filter_pattern_combo.addItem("Все паттерны")
            if not rows:
                return

            self.result_table.setColumnCount(len(header))
            self.result_table.setHorizontalHeaderLabels([LABELS.get(name, name) for name in header])
            for c, name in enumerate(header):
                self.result_table.horizontalHeaderItem(c).setToolTip(LABELS.get(name, name))
            table_header = self.result_table.horizontalHeader()
            for position, name in enumerate(('variant', 'enzyme', 'assay_rank',
                                              'genotype_quality', 'worst_margin',
                                              'genotype_margin', 'primer_pair_index')):
                if name in self.col_index:
                    logical = self.col_index[name]
                    table_header.moveSection(table_header.visualIndex(logical), position)
                    self.result_table.setColumnWidth(logical, 150 if position > 2 else 120)
            self.result_table.setRowCount(len(rows))

            self.result_table.setSortingEnabled(False)
            for r, row in enumerate(rows):
                for c, val in enumerate(row):
                    display = f'{val:.4f}'.rstrip('0').rstrip('.') if isinstance(val, float) else str(val)
                    item = NumericTableItem('' if val is None else display)
                    item.setData(Qt.ItemDataRole.UserRole, val)
                    item.setToolTip(f'{LABELS.get(header[c], header[c])}: {val if val is not None else "—"}')
                    self.result_table.setItem(r, c, item)

                pattern_val = str(rows[r][4]) if len(rows[r]) > 4 else ""
                if pattern_val == "no_enzyme_found":
                    for c in range(len(header)):
                        item = self.result_table.item(r, c)
                        if item:
                            item.setForeground(QtGui.QColor("#777777"))

            self.result_table.setSortingEnabled(True)
            self.populate_filters_from_results()
        except Exception:
            tb = traceback.format_exc()
            self.log_msg(f"[ERROR] Exception in show_results:\n{tb}", "error")

    def populate_filters_from_results(self):
        enzymes = set()
        patterns = set()
        for row in getattr(self, "out_rows", []):
            if len(row) > 2 and row[2]:
                enzymes.add(str(row[2]))
            if len(row) > 4 and row[4]:
                patterns.add(str(row[4]))

        current_enzyme = self.filter_enzyme_combo.currentText()
        current_pattern = self.filter_pattern_combo.currentText()

        self.filter_enzyme_combo.blockSignals(True)
        self.filter_pattern_combo.blockSignals(True)

        self.filter_enzyme_combo.clear()
        self.filter_enzyme_combo.addItem("Все ферменты")
        for e in sorted(enzymes):
            self.filter_enzyme_combo.addItem(e)

        self.filter_pattern_combo.clear()
        self.filter_pattern_combo.addItem("Все паттерны")
        for p in sorted(patterns):
            self.filter_pattern_combo.addItem(p)

        idx_e = self.filter_enzyme_combo.findText(current_enzyme)
        if idx_e != -1:
            self.filter_enzyme_combo.setCurrentIndex(idx_e)

        idx_p = self.filter_pattern_combo.findText(current_pattern)
        if idx_p != -1:
            self.filter_pattern_combo.setCurrentIndex(idx_p)

        self.filter_enzyme_combo.blockSignals(False)
        self.filter_pattern_combo.blockSignals(False)

        self.apply_filters()

    def apply_filters(self):
        if not hasattr(self, "out_rows"):
            return
        sel_enzyme = self.filter_enzyme_combo.currentText()
        sel_pattern = self.filter_pattern_combo.currentText()

        # Read values from the visible table order, because sorting changes
        # rows without changing out_rows.
        enzyme_col = self.col_index.get("enzyme", 2)
        pattern_col = self.col_index.get("pattern", 4)
        for r in range(self.result_table.rowCount()):
            enzyme_item = self.result_table.item(r, enzyme_col)
            pattern_item = self.result_table.item(r, pattern_col)
            enzyme = enzyme_item.text() if enzyme_item else ""
            pattern = pattern_item.text() if pattern_item else ""

            if sel_enzyme != "Все ферменты" and enzyme != sel_enzyme:
                self.result_table.setRowHidden(r, True)
                continue
            if sel_pattern != "Все паттерны" and pattern != sel_pattern:
                self.result_table.setRowHidden(r, True)
                continue
            self.result_table.setRowHidden(r, False)

    def on_table_selection_changed(self):
        """Обновляет состояние визуализации при выборе строки (данные готовятся для отдельного окна)."""
        try:
            selected = self.result_table.selectedItems()
            if not selected:
                self.assay_summary.setText('Выберите результат: оценка > 1 означает различимость в модели полос.')
                self._current_vis_payload = None
                self._current_vis_error = ""
                self.open_vis_btn.setEnabled(False)
                if self.vis_dialog is not None and self.vis_dialog.isVisible():
                    self.vis_dialog.panel.show_empty()
                return

            row = selected[0].row()
            self.assay_summary.setText(self._table_cell(row, 'reason'))
            if 'genotype_quality' in self.col_index:
                quality = self._table_cell(row, 'genotype_quality')
                if quality:
                    self.assay_summary.setText(
                        f"{self._table_cell(row, 'reason')}. "
                        f"Оценка: {self._table_cell(row, 'genotype_margin')}; "
                        f"при ухудшении: {self._table_cell(row, 'worst_margin')}.\n"
                        f"Полосы REF/REF: {self._table_cell(row, 'bands_ref_ref') or 'нет'} | "
                        f"REF/ALT: {self._table_cell(row, 'bands_ref_alt') or 'нет'} | "
                        f"ALT/ALT: {self._table_cell(row, 'bands_alt_alt') or 'нет'}")
            status_col = self.col_index.get("status")
            row_status = self.result_table.item(row, status_col).text() if status_col is not None and self.result_table.item(row, status_col) else ""
            if row_status in ("no_primers", "primer3_error", "no_discriminating_assay", "analysis_error"):
                self._current_vis_payload = None
                self._current_vis_error = self._table_cell(row, "reason") or "Визуализация недоступна: праймеры не подобраны."
                self.open_vis_btn.setEnabled(False)
                return
            payload, err = self._build_vis_payload_for_row(row)
            self._current_vis_payload = payload
            self._current_vis_error = err or ""
            self.open_vis_btn.setEnabled(payload is not None)

            # Если окно визуализации уже открыто — обновляем сразу
            if self.vis_dialog is not None and self.vis_dialog.isVisible():
                self.vis_dialog.show_payload(payload, message_if_empty=err or "")
        except Exception:
            tb = traceback.format_exc()
            self.log_msg(f"[ERROR] Ошибка при подготовке визуализации\n{tb}", "error")

    def _table_cell(self, row: int, col_name: str) -> str:
        idx = self.col_index.get(col_name)
        if idx is None:
            return ""
        item = self.result_table.item(row, idx)
        return item.text() if item is not None else ""

    @staticmethod
    def _safe_float(s: object) -> Optional[float]:
        try:
            if s is None:
                return None
            t = str(s).strip()
            if t == "" or t.lower() == "none":
                return None
            return float(t)
        except Exception:
            return None

    def _build_vis_payload_for_row(self, row: int) -> Tuple[Optional[dict], str]:
        """Собирает данные для визуализации (для отдельного окна).

        Возвращает (payload, error_message). Если payload=None — визуализация невозможна.
        """
        fasta_path = (self.last_params or {}).get("fasta")
        if not fasta_path:
            return None, "Не задан FASTA для визуализации. Запустите анализ заново."

        if self.fa_view is None:
            try:
                self.fa_view = Fasta(fasta_path, as_raw=True, sequence_always_upper=True)
            except TypeError:
                self.fa_view = Fasta(fasta_path, as_raw=True)
            except Exception as e:
                return None, f"Не удалось открыть FASTA: {fasta_path}\n{e}"

        mapped_id = self._table_cell(row, "mapped_id").strip()
        amp_start = _safe_int(self._table_cell(row, "amplicon_start"))
        amp_end = _safe_int(self._table_cell(row, "amplicon_end"))
        snp_offset = _safe_int(self._table_cell(row, "snp_offset_in_amplicon"))

        primer_left = self._table_cell(row, "primers_left").strip()
        primer_right = self._table_cell(row, "primers_right").strip()
        tm_left = self._safe_float(self._table_cell(row, "tm_left"))
        tm_right = self._safe_float(self._table_cell(row, "tm_right"))

        p_left_start = _safe_int(self._table_cell(row, "primer_left_start"))
        p_left_len = _safe_int(self._table_cell(row, "primer_left_len"))
        p_right_start = _safe_int(self._table_cell(row, "primer_right_start"))
        p_right_len = _safe_int(self._table_cell(row, "primer_right_len"))

        if not mapped_id or amp_start is None or amp_end is None:
            return None, "В выбранной строке отсутствуют данные для визуализации (mapped_id / amplicon_start / amplicon_end)."

        if amp_start < 1 or amp_end < amp_start:
            return None, f"Некорректные координаты ампликона: {mapped_id}:{amp_start}-{amp_end}"

        try:
            # pyfaidx: end в срезе является end-exclusive, поэтому подаём amp_end как inclusive
            amplicon_seq = str(self.fa_view[mapped_id][amp_start - 1 : amp_end])
        except Exception as e:
            return None, f"Не удалось получить последовательность из FASTA: {mapped_id}:{amp_start}-{amp_end}\n{e}"

        variant_str = self._table_cell(row, "variant").strip()
        enzyme_name = self._table_cell(row, "enzyme").strip()
        pattern = self._table_cell(row, "pattern").strip()
        frags_ref = self._table_cell(row, "frags_ref").strip()
        frags_alt = self._table_cell(row, "frags_alt").strip()

        payload = dict(
            mapped_id=mapped_id,
            genome_start=int(amp_start),
            genome_end=int(amp_end),
            amplicon_seq=amplicon_seq,
            snp_offset=int(snp_offset) if snp_offset is not None else None,
            primer_left=primer_left,
            primer_right=primer_right,
            tm_left=tm_left,
            tm_right=tm_right,
            primer_left_start=p_left_start,
            primer_left_len=p_left_len,
            primer_right_start=p_right_start,
            primer_right_len=p_right_len,
            variant_str=variant_str,
            enzyme_name=enzyme_name,
            pattern=pattern,
            frags_ref=frags_ref,
            frags_alt=frags_alt,
        )

        return payload, ""

    def open_visualization_window(self):
        """Открывает отдельное окно визуализации (адаптивное)."""
        if self.vis_dialog is None:
            self.vis_dialog = RFLPVisualizationDialog(self)

        # Показать окно и поднять наверх
        self.vis_dialog.show()
        self.vis_dialog.raise_()
        self.vis_dialog.activateWindow()

        if self._current_vis_payload is not None:
            self.vis_dialog.show_payload(self._current_vis_payload)
        else:
            msg = self._current_vis_error or "Выберите строку в таблице результатов."
            self.vis_dialog.show_payload(None, message_if_empty=msg)


    def _ensure_results_dir(self) -> str:
        base_dir = get_base_dir()
        results_dir = os.path.join(base_dir, "results")
        if not os.path.isdir(results_dir):
            os.makedirs(results_dir, exist_ok=True)
        return results_dir

    def _make_timestamp_name(self, ext: str) -> str:
        ts = datetime.now().strftime("%Y%m%d%H%M%S")
        return f"rflp_candidates_{ts}.{ext}"

    def _collect_table_data(self) -> Tuple[list, list]:
        rows = self.result_table.rowCount()
        cols = self.result_table.columnCount()
        if rows == 0 or cols == 0:
            return [], []

        header = list(self.out_header)
        data = []
        any_visible = False
        for r in range(rows):
            if self.result_table.isRowHidden(r):
                continue
            any_visible = True
            data.append([
                self.result_table.item(r, c).data(Qt.ItemDataRole.UserRole)
                if self.result_table.item(r, c) is not None else None
                for c in range(cols)
            ])

        if not any_visible:
            return [], []
        return header, data

    def export_csv(self):
        try:
            header, data = self._collect_table_data()
            if not data:
                QMessageBox.information(self, "Информация", "Нет видимых строк для экспорта")
                return
            results_dir = self._ensure_results_dir()
            filename = self._make_timestamp_name("csv")
            path = os.path.join(results_dir, filename)
            safe_export_csv(path, header, data)
            self.log_msg(f"[EXPORT] CSV сохранён: {path}")
            QMessageBox.information(self, "Экспорт CSV", f"Результаты сохранены:\n{path}")
        except Exception:
            tb = traceback.format_exc()
            self.log_msg(f"[ERROR] Exception in export_csv:\n{tb}", "error")

    def export_excel(self):
        try:
            header, data = self._collect_table_data()
            if not data:
                QMessageBox.information(self, "Информация", "Нет видимых строк для экспорта")
                return
            results_dir = self._ensure_results_dir()
            filename = self._make_timestamp_name("xlsx")
            path = os.path.join(results_dir, filename)

            safe_export_excel(path, header, data, metadata={
                "app_version": core_impl.APP_VERSION,
                "scoring_model": core_impl.MODEL_VERSION,
                "run_timestamp": datetime.now().isoformat(timespec="seconds"),
                "genome_path": (self.last_params or {}).get("fasta", ""),
                "parameters": self.last_params,
                "status_summary": (self.run_report or {}).get("summary", {"rows": len(data)}),
            }, rejected=(self.run_report or {}).get("rejected", []))
            self.log_msg(f"[EXPORT] Excel сохранён: {path}")
            QMessageBox.information(self, "Экспорт Excel", f"Результаты сохранены:\n{path}")
        except Exception:
            tb = traceback.format_exc()
            self.log_msg(f"[ERROR] Exception in export_excel:\n{tb}", "error")


def configure_system_theme(app: QApplication):
    """Use the Windows/system color scheme and refresh when it changes.

    The previous stylesheet forced white widgets, which made all controls
    unreadable when Windows was using a dark theme.  We set semantic palette
    roles once and generate the stylesheet from the detected palette instead.
    """
    state = {"updating": False}

    def refresh(*_):
        if state["updating"]:
            return
        state["updating"] = True
        try:
            hints = app.styleHints()
            scheme = getattr(hints, "colorScheme", None)
            scheme_value = scheme() if callable(scheme) else None
            scheme_name = str(scheme_value).lower()
            palette = app.palette()
            dark = "dark" in scheme_name or palette.color(QtGui.QPalette.ColorRole.Window).lightness() < 128

            if dark:
                colors = dict(window="#20252b", base="#171b20", alternate="#20262d", text="#f1f5f9",
                               muted="#a9b7c5", border="#4a5968", button="#29323b", hover="#344452",
                               selected="#315c70", selected_text="#ffffff", header="#2a3540",
                               disabled="#71808e", progress="#35b7a6")
            else:
                colors = dict(window="#f4f7fb", base="#ffffff", alternate="#f3f7fa", text="#22354b",
                               muted="#61758b", border="#c6d4df", button="#ffffff", hover="#e7f2f4",
                               selected="#cdebe7", selected_text="#123c40", header="#e8f0f5",
                               disabled="#8998a5", progress="#2c9c91")

            pal = QtGui.QPalette(palette)
            roles = {
                QtGui.QPalette.ColorRole.Window: colors["window"],
                QtGui.QPalette.ColorRole.Base: colors["base"],
                QtGui.QPalette.ColorRole.AlternateBase: colors["alternate"],
                QtGui.QPalette.ColorRole.Text: colors["text"],
                QtGui.QPalette.ColorRole.WindowText: colors["text"],
                QtGui.QPalette.ColorRole.Button: colors["button"],
                QtGui.QPalette.ColorRole.ButtonText: colors["text"],
                QtGui.QPalette.ColorRole.Highlight: colors["selected"],
                QtGui.QPalette.ColorRole.HighlightedText: colors["selected_text"],
                QtGui.QPalette.ColorRole.PlaceholderText: colors["muted"],
                QtGui.QPalette.ColorRole.ToolTipBase: colors["header"],
                QtGui.QPalette.ColorRole.ToolTipText: colors["text"],
            }
            for role, value in roles.items():
                pal.setColor(role, QtGui.QColor(value))
            pal.setColor(QtGui.QPalette.ColorGroup.Disabled, QtGui.QPalette.ColorRole.Text, QtGui.QColor(colors["disabled"]))
            pal.setColor(QtGui.QPalette.ColorGroup.Disabled, QtGui.QPalette.ColorRole.ButtonText, QtGui.QColor(colors["disabled"]))
            app.setPalette(pal)
            app.setStyleSheet(f"""
                QMainWindow, QDialog {{ background: {colors['window']}; color: {colors['text']}; }}
                QWidget {{ color: {colors['text']}; }}
                QGroupBox {{ font-weight: 600; border: 1px solid {colors['border']}; border-radius: 7px; margin-top: 10px; padding: 8px; }}
                QGroupBox::title {{ subcontrol-origin: margin; left: 10px; padding: 0 4px; color: {colors['text']}; }}
                QPushButton {{ padding: 6px 10px; border: 1px solid {colors['border']}; border-radius: 5px; background: {colors['button']}; color: {colors['text']}; }}
                QPushButton:hover {{ background: {colors['hover']}; border-color: {colors['progress']}; }}
                QPushButton:disabled {{ color: {colors['disabled']}; background: {colors['window']}; }}
                QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {{ padding: 5px; border: 1px solid {colors['border']}; border-radius: 4px; background: {colors['base']}; color: {colors['text']}; selection-background-color: {colors['selected']}; selection-color: {colors['selected_text']}; }}
                QTableWidget {{ background: {colors['base']}; alternate-background-color: {colors['alternate']}; color: {colors['text']}; gridline-color: {colors['border']}; selection-background-color: {colors['selected']}; selection-color: {colors['selected_text']}; }}
                QHeaderView::section {{ background: {colors['header']}; color: {colors['text']}; padding: 6px; border: none; border-bottom: 1px solid {colors['border']}; font-weight: 600; }}
                QTextEdit {{ background: {colors['base']}; color: {colors['text']}; border: 1px solid {colors['border']}; border-radius: 5px; }}
                QProgressBar {{ border: none; border-radius: 4px; background: {colors['header']}; color: {colors['text']}; text-align: center; }}
                QProgressBar::chunk {{ background: {colors['progress']}; border-radius: 4px; }}
                QToolTip {{ background: {colors['header']}; color: {colors['text']}; border: 1px solid {colors['border']}; }}
            """)
        finally:
            state["updating"] = False

    app.setStyle("Fusion")
    refresh()
    # Windows emits paletteChanged when the user changes the system theme.
    app.paletteChanged.connect(refresh)


def main():
    if not _PRIMER3_AVAILABLE:
        logger.warning(f"primer3 недоступен, ошибка импорта: {_PRIMER3_IMPORT_ERROR}")
    app = QApplication(sys.argv)
    configure_system_theme(app)
    win = RFLPPickerGUI()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
