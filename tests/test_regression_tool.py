import importlib.util
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location(
    'regression_tool', Path(__file__).resolve().parent.parent / 'tools' / 'regression.py')
regression = importlib.util.module_from_spec(spec)
spec.loader.exec_module(regression)


class CompatibilityTests(unittest.TestCase):
    base = dict(header=['a', 'b'], rows=[[1, 'x'], [2.5, None]], error=None,
                report=dict(summary=dict(rows=2), rejected=[dict(line=3)]))

    def candidate(self, **changes):
        cand = dict(header=['a', 'b', 'new'], rows=[[1, 'x', 'n'], [2.5, None, '']], error=None,
                    report=dict(summary=dict(rows=2, extra=1), rejected=[dict(line=3, more=1)], added=[]))
        cand.update(changes)
        return cand

    def test_appended_columns_and_fields_are_compatible(self):
        self.assertTrue(regression.compatible(self.base, self.candidate()))

    def test_changed_values_rows_or_types_are_not(self):
        for change in (dict(rows=[[1, 'y', 'n'], [2.5, None, '']]), dict(rows=[[1, 'x', 'n']]),
                       dict(rows=[[1.0, 'x', 'n'], [2.5, None, '']]), dict(header=['b', 'a', 'new']),
                       dict(report=dict(summary=dict(rows=3), rejected=[dict(line=3)])),
                       dict(report=dict(summary=dict(rows=2), rejected=[])), dict(error='boom')):
            self.assertFalse(regression.compatible(self.base, self.candidate(**change)), change)


if __name__ == '__main__':
    unittest.main()
