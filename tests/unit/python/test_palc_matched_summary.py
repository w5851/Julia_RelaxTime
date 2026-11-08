import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location("palc_summary", ROOT / "scripts/analysis/pnjl_bifurcation_spike/summarize_matched_pilot.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class MatchedPoints(unittest.TestCase):
    def point(self, rho, pressure=1.0, accepted=True):
        return dict(rho=rho, pressure=pressure, solution=[rho, 0.0], accepted=accepted)

    def test_nearby_density_is_not_same_target(self):
        result = MODULE.compare_points([self.point(0.35)], [self.point(0.351)])
        self.assertFalse(result["passed"])
        self.assertEqual(result["points"][0]["reason"], "missing_or_duplicate_target")

    def test_dense_grid_compares_only_shared_targets(self):
        self.assertTrue(MODULE.compare_points([self.point(0.35)],
            [self.point(0.325), self.point(0.35), self.point(0.375)])["passed"])

    def test_failed_or_wrong_branch_cannot_pass(self):
        for other in (self.point(0.35, accepted=False), self.point(0.35, pressure=2.0)):
            self.assertFalse(MODULE.compare_points([self.point(0.35)], [other])["passed"])


if __name__ == "__main__":
    unittest.main()
