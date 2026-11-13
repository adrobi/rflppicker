import os
import tempfile
import time
import unittest
from unittest.mock import patch

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
from PyQt6.QtCore import QSettings
from PyQt6.QtWidgets import QApplication, QComboBox

import genome_download
import launch_panel
from app_settings import AppSettings
from genome_dialog import GenomeDownloadDialog
from tests.helpers import dispose, write_text


def finish(dialog):
    """Let the download thread end and its signals arrive."""
    dialog.wait_for_thread()
    QApplication.processEvents()


class GenomeDownloadDialogTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self):
        self._dir = tempfile.TemporaryDirectory()
        self.addCleanup(self._dir.cleanup)
        self.calls = []

    def dialog(self, download):
        dialog = GenomeDownloadDialog(self._dir.name, download=download)
        self.addCleanup(dispose, dialog)
        return dialog

    def ok(self, accession, base_dir, progress_cb=None, status_cb=None, cancel_cb=None):
        self.calls.append((accession, base_dir))
        status_cb('Загрузка с NCBI…')
        progress_cb(5_000_000_000, 10_000_000_000)        # more than 32 bits
        return os.path.join(base_dir, 'genoms', 'x.fa'), os.path.join(base_dir, 'names', 'x.txt')

    def test_presets_and_another_assembly(self):
        dialog = self.dialog(self.ok)
        texts = [dialog.assembly_combo.itemText(i) for i in range(dialog.assembly_combo.count())]
        self.assertEqual(len(texts), len(genome_download.PRESETS) + 1)
        self.assertEqual(texts[1], 'Овца — ARS-UI_Ramb_v2.0 (GCF_016772045.1, 0,8 ГБ)')
        self.assertEqual(texts[-1], 'Другая сборка…')
        self.assertFalse(dialog.accession_edit.isEnabled())
        dialog.assembly_combo.setCurrentIndex(dialog.assembly_combo.count() - 1)
        self.assertTrue(dialog.accession_edit.isEnabled())

    def test_download_of_a_preset(self):
        dialog = self.dialog(self.ok)
        dialog.assembly_combo.setCurrentIndex(1)
        dialog.start()
        finish(dialog)
        self.assertEqual(self.calls, [('GCF_016772045.1', self._dir.name)])
        self.assertEqual(dialog.fasta_path, os.path.join(self._dir.name, 'genoms', 'x.fa'))
        self.assertEqual(dialog.result(), dialog.DialogCode.Accepted)
        self.assertEqual(dialog.progress.value(), 500)                 # per mille

    def test_typed_accession(self):
        dialog = self.dialog(self.ok)
        dialog.assembly_combo.setCurrentIndex(dialog.assembly_combo.count() - 1)
        dialog.accession_edit.setText(' gcf_000003195.3 ')
        dialog.start()
        finish(dialog)
        self.assertEqual(self.calls[0][0], 'gcf_000003195.3')

    def test_enter_starts_the_download(self):
        dialog = self.dialog(self.ok)
        self.assertTrue(dialog.download_button.isDefault())
        self.assertFalse(dialog.cancel_button.autoDefault())

    def test_error_is_shown_and_the_dialog_stays(self):
        def fail(*_args, **_kwargs):
            raise genome_download.GenomeDownloadError('Сборка GCF_1 не найдена на NCBI')

        dialog = self.dialog(fail)
        dialog.start()
        finish(dialog)
        self.assertEqual(dialog.error_label.text(), 'Сборка GCF_1 не найдена на NCBI')
        self.assertFalse(dialog.error_label.isHidden())
        self.assertEqual(dialog.fasta_path, '')
        self.assertTrue(dialog.download_button.isEnabled())
        self.assertNotEqual(dialog.result(), dialog.DialogCode.Accepted)

    def test_cancel_stops_the_download(self):
        seen = []

        def slow(accession, base_dir, progress_cb=None, status_cb=None, cancel_cb=None):
            while not cancel_cb():                 # until the dialog asks to stop
                time.sleep(0.01)
            seen.append(True)
            raise genome_download.DownloadCancelled()

        dialog = self.dialog(slow)
        dialog.start()
        dialog.reject()                            # «Отмена», Esc or the window cross during a download
        self.assertNotEqual(dialog.result(), dialog.DialogCode.Accepted)
        finish(dialog)
        self.assertEqual(seen, [True])
        self.assertEqual(dialog.fasta_path, '')
        self.assertTrue(dialog.error_label.isHidden())


class LaunchPanelDownloadTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def setUp(self):
        self._dir = tempfile.TemporaryDirectory()
        self.addCleanup(self._dir.cleanup)
        self.tmp = self._dir.name
        patcher = patch('launch_panel.get_base_dir', return_value=self.tmp)
        patcher.start()
        self.addCleanup(patcher.stop)
        self.settings = AppSettings(QSettings(os.path.join(self.tmp, 's.ini'), QSettings.Format.IniFormat))
        self.panel = launch_panel.LaunchPanel(self.settings)
        self.addCleanup(dispose, self.panel)
        self.combo = self.panel.findChild(QComboBox, 'genome_combo')

    def test_no_genome_is_chosen_in_an_empty_list(self):
        self.assertEqual([self.combo.itemText(i) for i in range(self.combo.count())],
                         ['Скачать геном…', 'Выбрать файл…'])
        self.assertEqual(self.combo.currentIndex(), -1)
        self.assertEqual(self.panel.genome_path(), '')

    def test_downloaded_genome_is_selected(self):
        fasta = write_text(self.tmp, 'Toy.fa', '>1\nACGT\n')

        class Downloaded:
            fasta_path = fasta

            def __init__(self, *_args, **_kwargs):
                pass

            def exec(self):
                return True

            def deleteLater(self):
                pass

        with patch.object(launch_panel, 'GenomeDownloadDialog', Downloaded):
            self.panel._on_genome_activated(self.combo.findText('Скачать геном…'))
        self.assertEqual(self.panel.genome_path(), fasta)
        self.assertEqual(self.settings.recent_genomes(), [fasta])

    def test_closed_download_dialog_keeps_the_previous_genome(self):
        class Closed:
            fasta_path = ''

            def __init__(self, *_args, **_kwargs):
                pass

            def exec(self):
                return False

            def deleteLater(self):
                pass

        with patch.object(launch_panel, 'GenomeDownloadDialog', Closed):
            self.panel._on_genome_activated(self.combo.findText('Скачать геном…'))
        self.assertEqual(self.panel.genome_path(), '')
        self.assertEqual(self.combo.currentIndex(), -1)


if __name__ == '__main__':
    unittest.main()
