import importlib.util
import contextlib
import io
import json
from pathlib import Path
import tempfile
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

    def summary_fixture(self, directory, *, wrong_new=False, wrong_legacy=False, truncate=False):
        methods = ["palc", "seed", "dense", "legacy_seed"]
        (directory / "processes.json").write_text(json.dumps(dict(expected_methods=methods)))
        targets = [0.3+i*0.01 for i in range(21)]
        for method in methods:
            grid = [0.3+i*0.005 for i in range(41)] if method == "dense" else targets
            pressure = 2.0 if (wrong_new and method == "seed") or (wrong_legacy and method == "legacy_seed") else 1.0
            points = [self.point(rho, pressure) for rho in grid]
            if truncate and method == "dense":
                points = points[:-1]
            reps = [dict(repetition=i, phase="first_call" if i == 1 else "warm",
                wall_seconds=1.0 if method == "palc" else 2.0, work={},
                result=dict(points=points, errors=[])) for i in range(1, 7)]
            payload = dict(config=dict(targets=targets), source_sha="same", load_seconds=0.1,
                expected_repetitions=6, repetitions=reps)
            (directory / f"{method}.json").write_text(json.dumps(payload))
        with contextlib.redirect_stdout(io.StringIO()):
            return MODULE.summarize(directory)

    def test_expanded_grid_and_legacy_difference(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = self.summary_fixture(Path(tmp), wrong_legacy=True)
        self.assertTrue(result["parity_passed"])
        self.assertFalse(result["checks"]["seed_vs_legacy_seed"]["passed"])
        self.assertEqual(result["seed_over_palc_warm_speedup"], 2.0)

    def test_new_branch_difference_blocks_speedup(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = self.summary_fixture(Path(tmp), wrong_new=True)
        self.assertFalse(result["parity_passed"])
        self.assertIsNone(result["seed_over_palc_warm_speedup"])

    def test_incomplete_dense_grid_blocks_acceptance(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = self.summary_fixture(Path(tmp), truncate=True)
        self.assertFalse(result["complete"])
        self.assertIsNone(result["seed_over_palc_warm_speedup"])


if __name__ == "__main__":
    unittest.main()
