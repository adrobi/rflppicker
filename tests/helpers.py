"""Shared helpers for the test suite."""
import os

import rflp_core

HEADER = rflp_core.BASE_HEADER + rflp_core.ASSAY_HEADER + rflp_core.V12_HEADER


def core_row(**values) -> list:
    """A core result row in HEADER order; unspecified columns are None."""
    unknown = set(values) - set(HEADER)
    if unknown:
        raise KeyError(f'Unknown columns: {sorted(unknown)}')
    return [values.get(name) for name in HEADER]


def write_text(directory: str, name: str, text: str) -> str:
    """Write UTF-8 text to directory/name and return the path."""
    path = os.path.join(directory, name)
    with open(path, 'w', encoding='utf-8', newline='') as handle:
        handle.write(text)
    return path
