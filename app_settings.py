"""What the program remembers between launches (spec §5.3).

Values are stored as JSON strings: the INI backend of QSettings would
otherwise return booleans and numbers as text.
"""
import json
import os
from typing import Any, List, Optional

from PyQt6.QtCore import QSettings

from rflp_core import RECOMMENDED_PARAMS

FILE_KEYS = ('names2_file', 'background_vcf')   # kept as files, not as parameters
MAX_RECENT_GENOMES = 5


def _same_type(value: Any, default: Any) -> bool:
    if isinstance(default, bool) or isinstance(value, bool):
        return isinstance(value, bool) and isinstance(default, bool)
    if isinstance(default, float):
        return isinstance(value, (int, float))
    if isinstance(default, list):
        return isinstance(value, list) and all(isinstance(item, str) for item in value)
    return isinstance(value, type(default))


class AppSettings:
    def __init__(self, qsettings: Optional[QSettings] = None):
        self._qs = qsettings or QSettings(QSettings.Format.IniFormat, QSettings.Scope.UserScope,
                                          'RFLP Picker', 'RFLP Picker')

    def value(self, key: str, default: Any = None) -> Any:
        raw = self._qs.value(key)
        if raw is None:
            return default
        try:
            return json.loads(raw)
        except (TypeError, ValueError):
            return default

    def set_value(self, key: str, value: Any) -> None:
        self._qs.setValue(key, json.dumps(value, ensure_ascii=False))

    def params(self) -> dict:
        """Recommended parameters overlaid with stored values of the right type."""
        params = {k: v for k, v in RECOMMENDED_PARAMS.items() if k not in FILE_KEYS}
        stored = self.value('params', {})
        if not isinstance(stored, dict):
            return params
        for key, default in params.items():
            if key in stored and _same_type(stored[key], default):
                value = stored[key]
                params[key] = float(value) if isinstance(default, float) else value
        return params

    def set_params(self, params: dict) -> None:
        keys = [k for k in RECOMMENDED_PARAMS if k not in FILE_KEYS]
        self.set_value('params', {k: params[k] for k in keys if k in params})

    def recent_genomes(self) -> List[str]:
        stored = self.value('recent_genomes', [])
        return [p for p in stored if isinstance(p, str)] if isinstance(stored, list) else []

    def add_recent_genome(self, path: str) -> None:
        key = os.path.normcase(path)
        others = [p for p in self.recent_genomes() if os.path.normcase(p) != key]
        self.set_value('recent_genomes', [path] + others[:MAX_RECENT_GENOMES - 1])

    def _names_map(self) -> dict:
        stored = self.value('names_for_genome', {})
        return stored if isinstance(stored, dict) else {}

    def names_for_genome(self, genome: str) -> Optional[str]:
        """None — not chosen by hand; '' — «Без соответствий»; otherwise the file."""
        return self._names_map().get(os.path.normcase(genome))

    def set_names_for_genome(self, genome: str, names: str) -> None:
        mapping = self._names_map()
        mapping[os.path.normcase(genome)] = names
        self.set_value('names_for_genome', mapping)
