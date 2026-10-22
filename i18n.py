"""Interface language: Russian sources, English translations (i18n spec §3).

Every visible text is written in Russian and passed through t(); in English
the Russian source is the key into translations_en.EN. Texts that the core
writes in Russian are translated at display time by core_text(). The language
is chosen once at startup, before any window is built. Plain Python, no Qt.
"""
import math
import re
from typing import Optional

import translations_en

LANGS = ('ru', 'en')
_language = 'ru'


def set_language(code: str) -> None:
    global _language
    if code not in LANGS:
        raise ValueError(f'Unsupported language: {code!r}')
    _language = code


def language() -> str:
    return _language


def N_(text: str) -> str:
    """Mark a constant as translatable; it is translated with t() when shown."""
    return text


def t(text: str, **values) -> str:
    if _language == 'en':
        text = translations_en.EN.get(text, text)
        if isinstance(text, tuple):   # (one, other): English plural chosen by the count value
            text = text[0] if values.get('count') == 1 else text[1]
    return text.format(**values) if values else text


def core_text(text: str) -> str:
    """Translate a text written by the core; parts it does not know stay as they are."""
    if _language == 'ru' or not text:
        return text
    for pattern, replacement in translations_en.CORE:
        text = re.sub(pattern, replacement, text)
    return text


def fmt(value, digits: int = 1) -> str:
    """A number with the decimal separator of the language; whole numbers stay whole."""
    if value is None or value == '':
        return ''
    if isinstance(value, int) or (isinstance(value, float) and not math.isfinite(value)):
        return str(value)
    text = f'{float(value):.{digits}f}'
    return text.replace('.', ',') if _language == 'ru' else text


def system_language(ui_languages: list) -> str:
    """'ru' when the first interface language of the system is Russian, otherwise 'en'."""
    return 'ru' if ui_languages and str(ui_languages[0]).lower().startswith('ru') else 'en'


def locale_language(name: Optional[str]) -> str:
    """Language from a locale name such as 'ru_RU' or 'Russian_Russia' (command line)."""
    name = (name or '').lower()
    return 'ru' if name.startswith(('ru', 'russian')) else 'en'


def resolve_language(setting, detected: str) -> str:
    """The saved choice if it is a known language, otherwise the detected one."""
    return setting if setting in LANGS else detected
