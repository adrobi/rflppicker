import os
import unittest
from unittest.mock import patch

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
from PyQt6.QtGui import QPalette
from PyQt6.QtWidgets import QApplication

import icons
import palette
from theme import MODE_TITLES, MODES, ThemeController, mono_font, stylesheet


class ThemeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QApplication.instance() or QApplication([])

    def test_modes_switch_palette(self):
        theme = ThemeController(self.app, 'light')
        self.assertEqual((theme.mode, theme.dark), ('light', False))
        self.assertEqual(self.app.palette().color(QPalette.ColorRole.Window).name().upper(), '#F3F3F3')
        seen = []
        theme.changed.connect(seen.append)
        theme.set_mode('dark')
        self.assertEqual(seen[-1]['accent'], '#85B7EB')
        self.assertEqual(theme.colors, palette.colors(True))
        self.assertIn('#202020', self.app.styleSheet().upper())
        self.assertEqual(self.app.font().family(), 'Segoe UI')

    def test_style_is_set_once(self):
        # Replacing the QStyle while widgets hold style-sheet wrappers of it crashed on Linux.
        ThemeController(self.app, 'light')
        with patch.object(self.app, 'setStyle') as set_style:
            theme = ThemeController(self.app, 'dark')
            theme.set_mode('light')
        set_style.assert_not_called()

    def test_unknown_mode_is_rejected(self):
        theme = ThemeController(self.app, 'light')
        with self.assertRaises(ValueError):
            theme.set_mode('sepia')

    def test_mode_titles(self):
        self.assertEqual(MODES, ('system', 'light', 'dark'))
        self.assertEqual(MODE_TITLES['dark'], 'Тёмная')

    def test_stylesheet_uses_the_palette(self):
        qss = stylesheet(palette.colors(False))
        self.assertIn('#185FA5', qss.upper())
        self.assertIn('min-height: 26px', qss)          # 28 px with the 1 px border

    def test_mono_font(self):
        self.assertIn(mono_font().family(), ('Cascadia Mono', 'Consolas'))

    def test_palettes_have_the_same_keys(self):
        self.assertEqual(set(palette.colors(False)), set(palette.colors(True)))

    def test_icon_is_drawn_in_the_requested_color(self):
        image = icons.icon('star', '#185FA5', filled=True).pixmap(16, 16).toImage()
        colors = {image.pixelColor(x, y).name().upper() for x in range(16) for y in range(16)
                  if image.pixelColor(x, y).alpha() == 255}
        self.assertIn('#185FA5', colors)

    def test_every_icon_renders(self):
        for name in icons.NAMES:
            self.assertFalse(icons.icon(name, '#1F1F1F').pixmap(16, 16).isNull(), name)


if __name__ == '__main__':
    unittest.main()
