"""The qualification collector must not report incomplete evidence as passed."""
import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location(
    "transport_joint_summary", ROOT / "scripts/analysis/relaxtime/summarize_transport_joint_validation.py")
SUMMARY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SUMMARY)


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


class JointSummaryTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.joint = self.root / "joint"
        self.joint.mkdir()
        self.config = {"candidate": {"p": 48}, "reference": {"p": 128}}
        self.sources = [{"path": "src/relaxtime/kernel.jl", "sha256": "test-numerical-version"}]
        (self.joint / "manifest.json").write_text(json.dumps({
            "source_files": self.sources, "parameters": {"configurations": self.config}}))
        (self.joint / "status.json").write_text(json.dumps({
            "complete": True, "comparisons": 18, "joint_target_passed": True}))
        write_csv(self.joint / "inputs/review_points.csv", [{"point_id": "a"}, {"point_id": "b"}])
        write_csv(self.joint / "comparisons.csv", [dict(
            point_id=point, observable=name, candidate=1.001, reference=1.0,
            relative_change=0.001, passed=True) for point in ("a", "b") for name in sorted(SUMMARY.OBSERVABLES)])
        self.timing = self.root / "timing"
        self.timing.mkdir()
        (self.timing / "manifest.json").write_text(json.dumps({"source_files": self.sources,
            "parameters": {"candidate": self.config["candidate"], "points": ["a", "b"], "repeats": 2}}))
        write_csv(self.timing / "inputs/review_points.csv", [{"point_id": "a"}, {"point_id": "b"}])
        self.timings = [dict(point_id=point, method=method, repeat=repeat,
                             total_seconds=3 + repeat, cpu_seconds=2 + repeat)
                        for point in ("a", "b") for method in ("none", "grouped", "thermal_weighted")
                        for repeat in (1, 2)]
        self.write_timings()

    def write_timings(self):
        write_csv(self.timing / "timings.csv", self.timings)
        write_csv(self.timing / "observables.csv", [dict(
            point_id=row["point_id"], method=row["method"], repeat=row["repeat"],
            **{name: 1.001 for name in SUMMARY.OBSERVABLES}) for row in self.timings])
        (self.timing / "status.json").write_text(json.dumps({
            "complete": True, "measurements": len(self.timings), "points": 2}))

    def summarize(self, expected=2, benchmark=None):
        return SUMMARY.summarize([self.joint], self.root / "out", expected, benchmark)

    def test_complete_evidence_and_medians(self):
        result = self.summarize(benchmark=self.timing)
        self.assertTrue(result["joint_target_passed"])
        self.assertAlmostEqual(result["point_maxima"][0]["relative_change"], 0.001)
        self.assertEqual(len(result["complete_timing"]["points"]), 6)
        self.assertEqual(result["complete_timing"]["points"][0]["median_wall_seconds"], 4.5)
        self.assertEqual(result["complete_timing"]["candidate_replay_max_relative"], 0.0)

    def test_missing_joint_point_is_not_passed(self):
        self.assertFalse(self.summarize(expected=3)["joint_target_passed"])

    def test_missing_entire_timing_point_is_rejected(self):
        self.timings = [row for row in self.timings if row["point_id"] == "a"]
        self.write_timings()
        with self.assertRaisesRegex(ValueError, "advertised point"):
            self.summarize(benchmark=self.timing)

    def test_repeated_measurement_does_not_replace_missing_repeat(self):
        self.timings[1]["repeat"] = 1
        self.write_timings()
        with self.assertRaisesRegex(ValueError, "repetitions"):
            self.summarize(benchmark=self.timing)

    def test_nonfinite_timing_is_rejected(self):
        self.timings[0]["cpu_seconds"] = "NaN"
        self.write_timings()
        with self.assertRaisesRegex(ValueError, "Invalid cpu_seconds"):
            self.summarize(benchmark=self.timing)

    def test_changed_numerical_code_is_rejected(self):
        path = self.timing / "manifest.json"
        manifest = json.loads(path.read_text())
        manifest["source_files"][0]["sha256"] = "different-numerical-version"
        path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError, "different numerical sources"):
            self.summarize(benchmark=self.timing)

    def test_timing_must_measure_the_qualified_numerical_result(self):
        path = self.timing / "observables.csv"
        rows = SUMMARY.read_rows(path)
        row = next(row for row in rows if row["method"] == "thermal_weighted")
        row["eta_over_s"] = 1.03
        write_csv(path, rows)
        with self.assertRaisesRegex(ValueError, "one-percent target"):
            self.summarize(benchmark=self.timing)


if __name__ == "__main__":
    unittest.main()
