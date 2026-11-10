import importlib.util
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


if __name__ == "__main__":
    unittest.main()
