"""Price tab for the currently displayed digest scheme."""
import http.client
import re
import time
from dataclasses import asdict
from decimal import Decimal
from html import escape

from PyQt6.QtCore import QDateTime, QLocale, QObject, QRunnable, Qt, QThreadPool, pyqtSignal
from PyQt6.QtWidgets import QHBoxLayout, QLabel, QPushButton, QScrollArea, QVBoxLayout, QWidget

from app_settings import AppSettings
from enzyme_prices import fetch_prices, is_product_link, sources_for
from i18n import language, t


class _PriceResult(QObject):
    finished = pyqtSignal(str, str, object)


class _PriceRequest(QRunnable):
    def __init__(self, enzyme, supplier):
        super().__init__()
        self.enzyme = enzyme
        self.supplier = supplier
        self.result = _PriceResult()

    def run(self):
        try:
            quotes = fetch_prices(self.supplier, self.enzyme)
        except (OSError, http.client.HTTPException, ValueError, LookupError, TypeError, ArithmeticError):
            data = None
        else:
            data = {'time': int(time.time()),
                    'quotes': [dict(asdict(quote), amount=str(quote.amount)) for quote in quotes]}
        self.result.finished.emit(self.enzyme, self.supplier, data)


class PricesPanel(QWidget):
    def __init__(self, parent=None, *, settings=None):
        super().__init__(parent)
        self._settings = settings or AppSettings()
        self._enzymes = ()
        self._pending = set()
        self._failed = set()

        layout = QVBoxLayout(self)
        buttons = QHBoxLayout()
        self.refresh_button = QPushButton(t('Обновить цены'))
        self.refresh_button.clicked.connect(self.refresh)
        buttons.addWidget(self.refresh_button)
        buttons.addStretch()
        layout.addLayout(buttons)

        self.details = QLabel()
        self.details.setWordWrap(True)
        self.details.setTextFormat(Qt.TextFormat.RichText)
        self.details.setTextInteractionFlags(Qt.TextInteractionFlag.TextBrowserInteraction)
        self.details.setOpenExternalLinks(True)
        self.details.setAlignment(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(self.details)
        layout.addWidget(scroll, 1)
        self._render()

    def set_enzymes(self, enzymes):
        """Return whether the digest changed, so the active tab can refresh it."""
        updated = tuple(dict.fromkeys(name for name in enzymes if name and re.fullmatch(r'[A-Za-z0-9]+', name)))
        changed = updated != self._enzymes
        self._enzymes = updated
        self._render()
        return changed

    def _keys(self):
        return [(name, supplier) for name in self._enzymes for supplier in sources_for(name)]

    @staticmethod
    def _cache_key(enzyme, supplier):
        return f'prices/{supplier}/{enzyme}'

    def refresh(self):
        for enzyme, supplier in self._keys():
            key = (enzyme, supplier)
            if key in self._pending:
                continue
            self._pending.add(key)
            self._failed.discard(key)
            request = _PriceRequest(enzyme, supplier)
            request.result.finished.connect(self._receive, Qt.ConnectionType.QueuedConnection)
            QThreadPool.globalInstance().start(request)
        self._render()

    def _receive(self, enzyme, supplier, data):
        key = (enzyme, supplier)
        self._pending.discard(key)
        if data is None:
            self._failed.add(key)
        else:
            self._settings.set_value(self._cache_key(enzyme, supplier), data)
        self._render()

    def _render(self):
        keys = self._keys()
        self.refresh_button.setText(t('Обновить цены'))
        self.refresh_button.setEnabled(bool(keys) and not any(key in self._pending for key in keys))
        if not keys:
            self.details.clear()
            return
        locale = QLocale('ru_RU' if language() == 'ru' else 'en_US')
        sections = []
        for enzyme, supplier in keys:
            key = (enzyme, supplier)
            cached = self._settings.value(self._cache_key(*key), {})
            quotes = cached.get('quotes', [])
            timestamp = cached.get('time')
            provider = t('Клонинг Фасилити') if supplier == 'Клонинг Фасилити' else supplier
            lines = [f'<b>{escape(enzyme)} · {escape(provider)}</b>']
            for quote in quotes:
                amount = Decimal(quote['amount'])
                whole = amount == amount.to_integral_value()
                formatted = locale.toString(float(amount), 'f', 0 if whole else 2)
                currency = '₽' if quote['currency'] == 'RUB' else quote['currency']
                pack = f" · {escape(quote['pack'])}" if quote['pack'] else ''
                lines.append(f"{escape(quote['sku'])}{pack} · <b>{formatted} {escape(currency)}</b>")
            if quotes and is_product_link(supplier, quotes[0]['url']):
                lines.append(f'<a href="{escape(quotes[0]["url"])}">{t("Карточка производителя")}</a>')
            if timestamp is not None:
                if not quotes:
                    lines.append(t('Публичная цена не найдена'))
                checked = locale.toString(QDateTime.fromSecsSinceEpoch(timestamp), QLocale.FormatType.ShortFormat)
                lines.append(t('Проверено {time}', time=checked))
            if key in self._pending:
                lines.append(t('Проверяем цены…'))
            elif key in self._failed:
                lines.append(t('Не удалось обновить'))
            elif timestamp is None:
                lines.append(t('Цена ещё не проверена'))
            sections.append('<br>'.join(lines))
        self.details.setText('<br><br>'.join(sections))
