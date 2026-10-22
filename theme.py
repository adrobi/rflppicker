"""Applies the palette to the application: light, dark or following Windows (spec §9)."""
from PyQt6.QtCore import QObject, Qt, pyqtSignal
from PyQt6.QtGui import QColor, QFont, QFontDatabase, QPalette
from PyQt6.QtWidgets import QApplication

import palette

MODES = ('system', 'light', 'dark')
MODE_TITLES = {'system': 'Системная', 'light': 'Светлая', 'dark': 'Тёмная'}


def stylesheet(c: dict) -> str:
    """Quiet flat style: 28 px controls, 4 px corners on controls, 8 px on panels, 1 px lines."""
    return f"""
        QMainWindow, QDialog {{ background: {c['window']}; color: {c['text']}; }}
        QWidget {{ color: {c['text']}; }}
        QFrame[role="panel"] {{ background: {c['panel']}; border: 1px solid {c['line']}; border-radius: 8px; }}
        QLabel[role="caption"] {{ color: {c['text2']}; font-size: 8.5pt; }}
        QLabel[role="title"] {{ font-size: 11pt; font-weight: 500; }}
        QLabel[role="error"] {{ color: {c['error_fg']}; }}
        QLabel[role="muted"] {{ color: {c['text2']}; }}
        QPushButton, QToolButton {{ min-height: 26px; padding: 0 10px; border: 1px solid {c['line']};
            border-radius: 4px; background: {c['panel']}; color: {c['text']}; }}
        QPushButton:hover, QToolButton:hover {{ border-color: {c['accent']}; }}
        QPushButton:disabled, QToolButton:disabled {{ color: {c['text2']}; background: {c['window']}; }}
        QPushButton[role="primary"] {{ background: {c['accent']}; color: {c['on_accent']}; border-color: {c['accent']};
            font-weight: 500; }}
        QPushButton[role="link"], QToolButton[role="link"] {{ border: none; background: transparent; color: {c['accent']};
            padding: 0; text-align: left; }}
        QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {{ min-height: 26px; padding: 0 6px;
            border: 1px solid {c['line']}; border-radius: 4px; background: {c['panel']}; color: {c['text']};
            selection-background-color: {c['selected_bg']}; selection-color: {c['selected_text']}; }}
        QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus {{ border-color: {c['accent']}; }}
        QTableView, QListWidget, QPlainTextEdit, QTextEdit {{ background: {c['panel']}; color: {c['text']};
            border: 1px solid {c['line']}; border-radius: 8px; gridline-color: {c['line']};
            selection-background-color: {c['selected_bg']}; selection-color: {c['selected_text']}; }}
        QHeaderView::section {{ background: {c['panel']}; color: {c['text2']}; padding: 0 6px; min-height: 28px;
            border: none; border-bottom: 1px solid {c['line']}; }}
        QTabWidget::pane {{ border: none; }}
        QTabBar::tab {{ padding: 4px 12px; border: none; border-bottom: 2px solid transparent;
            color: {c['text2']}; background: transparent; }}
        QTabBar::tab:selected {{ color: {c['accent']}; border-bottom-color: {c['accent']}; }}
        QProgressBar {{ border: none; border-radius: 4px; background: {c['line']}; max-height: 6px; }}
        QProgressBar::chunk {{ background: {c['accent']}; border-radius: 4px; }}
        QGroupBox {{ border: none; margin-top: 14px; font-size: 8.5pt; color: {c['text2']}; }}
        QGroupBox::title {{ subcontrol-origin: margin; left: 0; color: {c['text2']}; }}
        QScrollArea {{ background: transparent; border: none; }}
        QMenu {{ background: {c['panel']}; border: 1px solid {c['line']}; }}
        QMenu::item:selected {{ background: {c['selected_bg']}; color: {c['selected_text']}; }}
        QToolTip {{ background: {c['panel']}; color: {c['text']}; border: 1px solid {c['line']}; }}
        QStatusBar {{ background: {c['window']}; color: {c['text2']}; }}
    """


def mono_font(size: float = 9) -> QFont:
    family = 'Cascadia Mono' if 'Cascadia Mono' in QFontDatabase.families() else 'Consolas'
    font = QFont(family)
    font.setPointSizeF(size)
    return font


class ThemeController(QObject):
    """Owns the application look; ``changed`` carries the active palette to the widgets."""
    changed = pyqtSignal(dict)

    def __init__(self, app: QApplication, mode: str = 'system'):
        super().__init__(app)
        self._app = app
        self._mode = 'system'
        self._colors = palette.colors(False)
        # Set the style once: replacing it later deletes the QStyle that widgets
        # with their own style sheets still reference (segfault on Linux).
        # (After a style sheet is set, app.style() is Qt's wrapper, so its name does not tell.)
        if not app.property('rflp_fusion_style'):
            app.setStyle('Fusion')
            app.setProperty('rflp_fusion_style', True)
        app.styleHints().colorSchemeChanged.connect(self._on_system_scheme)
        self.set_mode(mode)

    @property
    def mode(self) -> str:
        return self._mode

    @property
    def colors(self) -> dict:
        return dict(self._colors)

    @property
    def dark(self) -> bool:
        return self._colors['window'] == palette.DARK['window']

    def _system_dark(self) -> bool:
        return self._app.styleHints().colorScheme() == Qt.ColorScheme.Dark

    def _on_system_scheme(self, *_):
        if self._mode == 'system':
            self.set_mode('system')

    def set_mode(self, mode: str) -> None:
        if mode not in MODES:
            raise ValueError(f'Unknown theme mode: {mode}')
        self._mode = mode
        dark = self._system_dark() if mode == 'system' else mode == 'dark'
        c = self._colors = palette.colors(dark)
        app = self._app
        pal = QPalette()
        roles = {
            QPalette.ColorRole.Window: c['window'], QPalette.ColorRole.Base: c['panel'],
            QPalette.ColorRole.AlternateBase: c['window'], QPalette.ColorRole.Text: c['text'],
            QPalette.ColorRole.WindowText: c['text'], QPalette.ColorRole.Button: c['panel'],
            QPalette.ColorRole.ButtonText: c['text'], QPalette.ColorRole.Highlight: c['selected_bg'],
            QPalette.ColorRole.HighlightedText: c['selected_text'], QPalette.ColorRole.PlaceholderText: c['text2'],
            QPalette.ColorRole.ToolTipBase: c['panel'], QPalette.ColorRole.ToolTipText: c['text'],
            QPalette.ColorRole.Link: c['accent'], QPalette.ColorRole.Mid: c['line'],
        }
        for role, value in roles.items():
            pal.setColor(role, QColor(value))
        for role in (QPalette.ColorRole.Text, QPalette.ColorRole.ButtonText, QPalette.ColorRole.WindowText):
            pal.setColor(QPalette.ColorGroup.Disabled, role, QColor(c['text2']))
        app.setPalette(pal)
        font = QFont('Segoe UI')
        font.setPointSizeF(9)
        app.setFont(font)
        app.setStyleSheet(stylesheet(c))
        self.changed.emit(dict(c))
