import tempfile
import unittest
from decimal import Decimal
from pathlib import Path
from unittest.mock import patch

from PyQt6.QtCore import QSettings, QThreadPool
from PyQt6.QtWidgets import QApplication

from app_settings import AppSettings
from enzyme_prices import Quote
from prices_panel import PricesPanel


class PricesPanelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def test_refresh_failure_retains_the_last_price_and_timestamp(self):
        quote = Quote('E013', '500 е.а.', Decimal(6500), 'RUB', 'https://sibenzyme.com/product/acs-i/')
        with tempfile.TemporaryDirectory() as tmp:
            settings = AppSettings(QSettings(str(Path(tmp) / 'settings.ini'), QSettings.Format.IniFormat))
            panel = PricesPanel(settings=settings)
            panel.set_enzymes(('AcsI',))
            with patch('prices_panel.fetch_prices', return_value=(quote,)):
                panel.refresh()
                self.assertTrue(QThreadPool.globalInstance().waitForDone(3000))
                self.app.processEvents()
            self.assertIn('E013', panel.details.text())
            cached = settings.value('prices/SibEnzyme/AcsI')
            self.assertIsInstance(cached['time'], int)
            with patch('prices_panel.fetch_prices', side_effect=TimeoutError):
                panel.refresh()
                self.assertTrue(QThreadPool.globalInstance().waitForDone(3000))
                self.app.processEvents()
            self.assertIn('E013', panel.details.text())
            self.assertIn('Не удалось обновить', panel.details.text())
            self.assertEqual(settings.value('prices/SibEnzyme/AcsI'), cached)
            another = PricesPanel(settings=settings)
            another.set_enzymes(('AcsI',))
            self.assertIn('E013', another.details.text())
            another.set_enzymes(('AciI',))
            self.assertNotIn('E013', another.details.text())
