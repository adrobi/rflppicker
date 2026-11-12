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


class ProfileTests(unittest.TestCase):
    def test_base_parameters_are_the_program_defaults(self):
        import rflp_core
        stated = {key: value for key, value in regression.BASE_PARAMS.items()
                  if key not in regression.REMOVED_PARAMS + ('top_results',)}
        self.assertEqual(stated, {key: rflp_core.DEFAULT_PARAMS.get(key) for key in stated})
        self.assertEqual(set(stated) | {'top_results'}, set(rflp_core.DEFAULT_PARAMS))

    def test_checks_off_matches_the_test_helper(self):
        from tests.helpers import CHECKS_OFF
        self.assertEqual(regression.CHECKS_OFF, CHECKS_OFF)

    def test_removed_parameters_do_not_count_as_a_difference(self):
        # 1.4.0 echoes joint_design in the report, later versions drop the key.
        report = dict(parameters=dict(flank=300, joint_design=True), summary=dict(rows=1))
        self.assertEqual(regression.comparable_report(report),
                         dict(parameters=dict(flank=300), summary=dict(rows=1)))
        self.assertEqual(report['parameters'], dict(flank=300, joint_design=True))   # the input is not changed
        self.assertEqual(regression.comparable_report({}), {})
