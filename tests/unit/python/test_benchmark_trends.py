from __future__ import annotations

import copy
import importlib.util
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location("benchmark_trends", ROOT / "scripts/dev/benchmark_trends.py")
trends = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(trends)


def raw_reports(directory: Path):
    environment = {key: "fixture" for key in trends.ENVIRONMENT_KEYS}
    environment.update(julia_threads=1, blas_threads=1, cpu_threads=2, julia_version="1.12.5")
    metadata = {
        "provenance": dict(commit="a" * 40, branch="main", source_dirty=False,
                           repository="owner/repo", run_id="30", run_attempt="1", event="schedule"),
        "environment": environment,
    }
    common = dict(schema_version="pnjl_benchmark_raw_v1", generated_at="2026-10-04T02:00:00Z",
                  metadata=metadata, samples=10, methodology=dict(phase="warm", evals=1, seconds=5))
    single = dict(common, config=dict(T_mev=50, p_num=24), benchmarks=[
        dict(label="solver", params="method=:trust_region",
             stats=dict(samples=10, median=100, min=90, memory_bytes=1024, allocs=20))])
    scan = dict(common, tmu_config=dict(T=[50, 100], mu=[0]), trho_config=dict(T=[80], rho=[0, 1]),
                benchmarks=[dict(label="scan", n_points=2, convergence_rate=1.0,
                                 stats=dict(samples=5, per_point_median_ms=80, memory_mb=1, allocs=50))])
    for name, report in (("single_point_benchmark.json", single), ("scan_benchmark.json", scan)):
        trends.write_json(directory / name, report)
    return single, scan


class BenchmarkTrendsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        raw_reports(self.directory)
        self.current = trends.build_snapshot(self.directory, "success")
        self.baseline = copy.deepcopy(self.current)
        self.baseline["generated_at"] = "2026-09-28T02:00:00Z"
        self.baseline["provenance"].update(run_id="20", commit="b" * 40)

    def test_snapshot_records_units_config_samples_and_provenance(self):
        self.assertTrue(self.current["baseline_eligible"])
        self.assertEqual(len(self.current["metrics"]), 8)
        self.assertEqual(self.current["metrics"]["scan / scan / memory"]["value"], 1024**2)
        self.assertEqual(self.current["config"]["single"]["config"]["p_num"], 24)
        self.assertEqual(self.current["samples"]["scan / scan"], 5)
        self.assertEqual(self.current["environment"]["julia_threads"], 1)

    def test_first_run_writes_complete_reports_and_action_summary(self):
        summary_path = self.directory / "action-summary.md"
        with patch.dict(os.environ, {"GITHUB_STEP_SUMMARY": str(summary_path)}):
            self.assertEqual(trends.main(["--results-dir", str(self.directory),
                                          "--threshold-status", "success"]), 0)
        self.assertEqual(trends.read_json(self.directory / "comparison.json")["status"], "missing_baseline")
        self.assertTrue((self.directory / trends.SNAPSHOT).is_file())
        self.assertIn("仅记录", summary_path.read_text(encoding="utf-8"))

    def test_changes_warn_without_failing_and_zero_is_defined(self):
        self.current["metrics"]["single / solver / median"]["value"] = 130
        self.current["metrics"]["single / solver / allocations"]["value"] = 0
        self.baseline["metrics"]["single / solver / allocations"]["value"] = 0
        self.baseline["metrics"]["single / solver / memory"]["value"] = 0
        self.current["metrics"]["scan / scan / convergence"]["value"] = 95
        result = trends.compare(self.current, self.baseline)
        self.assertEqual(result["status"], "compared")
        self.assertEqual(result["warning_count"], 3)
        rows = {row["metric"]: row for row in result["rows"]}
        self.assertAlmostEqual(rows["single / solver / median"]["delta"], 30)
        self.assertEqual(rows["single / solver / allocations"]["delta"], 0)
        self.assertIsNone(rows["single / solver / memory"]["delta"])
        self.assertEqual(rows["scan / scan / convergence"]["delta_unit"], "pp")
        baseline_path = self.directory / "previous.json"
        trends.write_json(baseline_path, self.baseline)
        self.assertEqual(trends.main(["--results-dir", str(self.directory), "--baseline", str(baseline_path),
                                     "--threshold-status", "success"]), 0)

    def test_expired_incompatible_and_invalid_baselines_are_distinct(self):
        old = copy.deepcopy(self.baseline)
        old["generated_at"] = "2026-01-01T00:00:00Z"
        self.assertEqual(trends.compare(self.current, old)["status"], "expired_baseline")
        for field in ("julia_threads", "cpu_model", "manifest_sha256", "workload_sha256", "pnjl_profile", "physics_profile"):
            changed = copy.deepcopy(self.baseline)
            changed["environment"][field] = "different"
            self.assertEqual(trends.compare(self.current, changed)["status"], "incompatible_baseline")
        changed = copy.deepcopy(self.baseline)
        changed["config"]["single"]["config"]["p_num"] = 48
        self.assertEqual(trends.compare(self.current, changed)["status"], "incompatible_baseline")
        for value in (0, -1, float("nan"), float("inf")):
            invalid = copy.deepcopy(self.baseline)
            invalid["metrics"]["single / solver / median"]["value"] = value
            self.assertEqual(trends.compare(self.current, invalid)["status"], "invalid_baseline")
        future = copy.deepcopy(self.baseline)
        future["generated_at"] = "2026-10-05T00:00:00Z"
        self.assertEqual(trends.compare(self.current, future)["status"], "invalid_baseline")

    def test_bad_or_partial_current_reports_fail(self):
        single, scan = raw_reports(self.directory)
        scan["metadata"] = copy.deepcopy(scan["metadata"])
        scan["metadata"]["provenance"]["commit"] = "c" * 40
        trends.write_json(self.directory / "scan_benchmark.json", scan)
        self.assertEqual(trends.main(["--results-dir", str(self.directory)]), 1)
        self.assertEqual(trends.read_json(self.directory / "comparison.json")["status"], "invalid_current")
        raw_reports(self.directory)
        (self.directory / "scan_benchmark.json").unlink()
        self.assertEqual(trends.main(["--results-dir", str(self.directory)]), 1)
        for value in (float("nan"), float("inf"), -1, 0):
            single, _ = raw_reports(self.directory)
            single["benchmarks"][0]["stats"]["median"] = value
            (self.directory / "single_point_benchmark.json").write_text(json.dumps(single), encoding="utf-8")
            with self.assertRaises(ValueError):
                trends.build_snapshot(self.directory, "success")

    def test_dirty_pr_or_failed_runs_cannot_become_main_baselines(self):
        for changes in (dict(branch="feature", event="pull_request"), dict(source_dirty=True), dict(run_id="")):
            single, scan = raw_reports(self.directory)
            single["metadata"]["provenance"].update(changes)
            scan["metadata"]["provenance"].update(changes)
            trends.write_json(self.directory / "single_point_benchmark.json", single)
            trends.write_json(self.directory / "scan_benchmark.json", scan)
            self.assertFalse(trends.build_snapshot(self.directory, "success")["baseline_eligible"])
        raw_reports(self.directory)
        self.assertFalse(trends.build_snapshot(self.directory, "failure")["baseline_eligible"])

    def test_github_uses_latest_compatible_successful_main_run(self):
        def run(run_id):
            return dict(id=run_id, head_branch="main", conclusion="success", event="schedule",
                        created_at="2026-09-28T01:00:00Z", head_sha="b" * 40,
                        html_url=f"https://github.com/owner/repo/actions/runs/{run_id}")

        runs = [run(29), run(28), run(27), run(26), run(20)]
        runs[0]["event"] = "pull_request"
        runs[1]["conclusion"] = "failure"
        runs[2]["head_branch"] = "feature"
        downloads = []

        def fake_gh(*args):
            if "workflows/" in args[-1]:
                return json.dumps(dict(workflow_runs=runs))
            if args[0] == "api":
                return json.dumps(dict(artifacts=[dict(name="pnjl-benchmark", expired=False)]))
            run_id = args[2]
            downloads.append(run_id)
            snapshot = copy.deepcopy(self.baseline)
            snapshot["provenance"]["run_id"] = run_id
            if run_id == "26":
                snapshot["config"]["single"]["config"]["p_num"] = 48
            trends.write_json(Path(args[-1]) / trends.SNAPSHOT, snapshot)
            return ""

        with patch.object(trends, "gh", side_effect=fake_gh):
            baseline, result = trends.github_baseline(self.current, "owner/repo", 35)
        self.assertEqual(downloads, ["26", "20"])
        self.assertEqual(baseline["provenance"]["run_id"], "20")
        self.assertEqual(result["status"], "compared")
        self.assertEqual(result["rejected_candidates"][0]["status"], "incompatible_baseline")

    def test_unavailable_and_expired_artifacts_are_nonfatal(self):
        with patch.object(trends, "gh", side_effect=RuntimeError("unavailable")):
            baseline, result = trends.github_baseline(self.current, "owner/repo", 35)
        self.assertIsNone(baseline)
        self.assertEqual(result["status"], "baseline_unavailable")
        run = dict(id=20, head_branch="main", conclusion="success", event="push",
                   created_at="2026-09-28T00:00:00Z", head_sha="b" * 40, html_url="https://example.test/run")
        responses = [json.dumps(dict(workflow_runs=[run])),
                     json.dumps(dict(artifacts=[dict(name="pnjl-benchmark", expired=True)]))]
        with patch.object(trends, "gh", side_effect=responses):
            _, result = trends.github_baseline(self.current, "owner/repo", 35)
        self.assertEqual(result["status"], "expired_baseline")


if __name__ == "__main__":
    unittest.main()
