import importlib.util
import argparse
from pathlib import Path
import unittest

FILE = Path(__file__).resolve().parents[3] / "scripts/analysis/pnjl_bifurcation_spike/compare_production.py"
SPEC = importlib.util.spec_from_file_location("palc_production", FILE)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def row():
    return dict(T=120.0, xi=0.0, coordinate=1.5, converged=True, residual=1e-9,
                solution=[-1., -1., -2., .2, .2], masses=[1., 1., 2.], mu_vec=[1., 1., 1.],
                pressure=20.0, rho_norm=1.5, entropy=1.0, energy=3.0, density_error=1e-10, constraint_error=1e-10)


class ComparisonTests(unittest.TestCase):
    def test_rejected_roots_do_not_count_as_parity(self):
        bad = dict(row(), masses=[3., 1., 2.])
        result = MODULE.compare_rows([bad], [bad])
        self.assertFalse(result["same_admissible_results"])
        self.assertEqual(len(result["coverage"]), 1)

    def test_large_state_difference_is_visible(self):
        changed = dict(row(), mu_vec=[1.02, 1., 1.])
        result = MODULE.compare_rows([row()], [changed])
        self.assertEqual(len(result["mismatches"]), 1)

    def test_missing_points_fail_completeness(self):
        self.assertFalse(MODULE.compare_rows([row()], [])["same_admissible_results"])

    def test_small_roundoff_passes(self):
        changed = dict(row(), pressure=20.0+1e-9)
        self.assertTrue(MODULE.compare_rows([row()], [changed])["same_admissible_results"])

    def test_empty_or_duplicate_workloads_do_not_pass(self):
        self.assertFalse(MODULE.compare_rows([], [])["same_admissible_results"])
        self.assertFalse(MODULE.compare_rows([row(), row()], [row()])["same_admissible_results"])

    def test_nonfinite_thermodynamics_is_rejected(self):
        self.assertFalse(MODULE.admissible(dict(row(), pressure=float("nan"))))

    def test_configured_grid_must_be_complete(self):
        data = dict(config=dict(temperatures=[120.0], anisotropies=[0.0], grid=[1.5, 2.0]))
        self.assertFalse(MODULE.complete_round(data, dict(rows=[row()])))

    def test_no_methods_cannot_pass_promotion(self):
        self.assertTrue(MODULE.promotion_checks({}, 2))

    def test_case_selection_preserves_order_and_rejects_duplicates(self):
        self.assertEqual(MODULE.parse_cases("tmu:legacy,rho:palc"), (("tmu", "legacy"), ("rho", "palc")))
        self.assertEqual(len(MODULE.parse_cases("all")), 9)
        for value in ("rho:other", "rho:palc,rho:palc", ""):
            with self.assertRaises(argparse.ArgumentTypeError):
                MODULE.parse_cases(value)

    def test_subset_check_does_not_require_unrequested_work(self):
        method = dict(repetitions=2, complete_grids=True, source_sha="a"*40,
                      repeat_checks=[dict(same_admissible_results=True)], accepted=1, total=1, config={})
        summary = dict(tmu=dict(methods=dict(legacy=method)))
        self.assertEqual(MODULE.promotion_checks(summary, 2, (("tmu", "legacy"),)), [])
        self.assertTrue(MODULE.promotion_checks(summary, 2))
        method["accepted"] = 0
        self.assertTrue(MODULE.promotion_checks(summary, 2, (("tmu", "legacy"),)))


if __name__ == "__main__":
    unittest.main()
