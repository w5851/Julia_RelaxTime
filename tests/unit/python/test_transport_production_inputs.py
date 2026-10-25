import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location(
    "production_inputs", ROOT / "scripts/analysis/relaxtime/prepare_transport_production_inputs.py")
inputs = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(inputs)


class ProductionInputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        root = Path(self.temp.name)
        self.phase, self.replay, self.output = root / "phase", root / "replay", root / "output"
        self.phase.mkdir()
        self.replay.mkdir()
        self.checkout = root / "checkout"
        self.numerical_source = self.checkout / "src" / "kernel.jl"
        self.numerical_source.parent.mkdir(parents=True)
        self.numerical_source.write_text("qualified numerical implementation\n")
        self.source_entry = dict(path="src/kernel.jl", sha256=hashlib.sha256(self.numerical_source.read_bytes()).hexdigest())
        for folder in (self.phase, self.replay):
            snapshot = folder / "source" / "src" / "kernel.jl"
            snapshot.parent.mkdir(parents=True)
            snapshot.write_bytes(self.numerical_source.read_bytes())
        (self.replay / "inputs").mkdir()
        points = [(0, 1, -.22), (0, 1, .37), (450, 1.1, .35), (900, 1, -.04),
                  (0, 1.2, -.5), (900, 1, -.003), (900, 1, .003)]
        self.plan = [dict(T_MeV=200 * alpha, muB_MeV=mu, alpha_T=alpha, xi=xi,
                          coexistence_side="quark_side" if (mu, xi) == (900, -.003) else "none")
                     for mu, alpha, xi in points]
        self.anchors = [dict(muB_MeV=mu, certified="true", refinement_delta_MeV=.0001) for mu in (0, 450, 900)]
        self.background = [dict(row, distinct_branches="false", refined_distinct_branches="false",
                                stable_phase="quark", refined_stable_phase="quark", passed="true",
                                max_mass_relative_change=.0001, entropy_relative_change=.0002)
                           for row in self.plan]
        self.comparisons = [dict(point_id=f"point-{i}", observable=name, production=1.0,
                                 qualified_candidate=1.0, qualified_reference=1.0,
                                 replay_relative_change=0, reference_relative_change=0, passed="true")
                            for i in range(12) for name in sorted(inputs.OBSERVABLES)]
        self.targets = [dict(point_id=r["point_id"], observable=r["observable"], candidate=1.0, reference=1.0)
                        for r in self.comparisons]
        self.save()

    def write_csv(self, path, rows):
        with path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    def save(self):
        self.write_csv(self.phase / "sampling_plan.csv", self.plan)
        self.write_csv(self.phase / "anchors.csv", self.anchors)
        self.write_csv(self.phase / "background.csv", self.background)
        self.write_csv(self.replay / "comparisons.csv", self.comparisons)
        reference_path = self.replay / "inputs" / "qualified_comparisons.csv"
        self.write_csv(reference_path, self.targets)
        (self.phase / "status.json").write_text(json.dumps(dict(
            complete=True, phase_plan_passed=True, points=len(self.plan), completed_points=len(self.plan), failures=0)))
        (self.phase / "manifest.json").write_text(json.dumps(dict(source_files=[self.source_entry], parameters=dict(
            profile=dict(id="thermal_weighted_v1", source_sha256="b" * 64),
            reference_comparisons_sha256=hashlib.sha256(reference_path.read_bytes()).hexdigest()))))
        (self.replay / "manifest.json").write_text((self.phase / "manifest.json").read_text())
        (self.replay / "status.json").write_text(json.dumps(dict(
            complete=True, production_replay_passed=True, points=12, completed_points=12,
            comparisons=len(self.comparisons), failures=0)))

    def prepare(self, sha=None):
        return inputs.prepare(self.phase, [self.replay], self.output, "new_case", sha, shard_points=2,
                              source_root=self.checkout)

    def fixture_commit(self, *, add_sources=True):
        # All Git mutations here are confined to the disposable test checkout.
        subprocess.run(["git", "init", "-q", str(self.checkout)], check=True, capture_output=True)
        subprocess.run(["git", "-C", str(self.checkout), "config", "core.autocrlf", "true"],
                       check=True, capture_output=True)
        if add_sources:
            subprocess.run(["git", "-C", str(self.checkout), "add", "src"], check=True, capture_output=True)
        subprocess.run(["git", "-C", str(self.checkout), "-c", "user.name=Transport test fixture",
                        "-c", "user.email=transport-test@example.invalid", "-c", "commit.gpgsign=false",
                        "commit", "--allow-empty", "-qm", "fixture"], check=True, capture_output=True)
        return subprocess.check_output(["git", "-C", str(self.checkout), "rev-parse", "HEAD"],
                                       text=True).strip()

    def test_six_point_pilot_and_full_plan_are_preserved(self):
        result = self.prepare()
        self.assertEqual(sum(j["expected_points"] for j in result["pilot_jobs"]), 6)
        self.assertEqual(sum(j["expected_points"] for j in result["provisional_full_jobs"]), len(self.plan))
        self.assertFalse(result["pilot_input_ready"])
        self.assertFalse(result["full_dispatch_ready"])
        self.assertTrue(all(j["inputs"]["sampling_profile"] == "thermal_weighted_v1" for j in result["pilot_jobs"]))

    def test_commit_only_enables_pilot_input(self):
        sha = self.fixture_commit()
        (self.checkout / "README.md").write_text("unrelated documentation\n")
        result = self.prepare(sha)
        self.assertTrue(result["pilot_input_ready"])
        self.assertFalse(result["full_dispatch_ready"])

    def test_commit_with_different_numerical_content_is_rejected(self):
        qualified = self.numerical_source.read_bytes()
        self.numerical_source.write_text("old numerical implementation\n")
        sha = self.fixture_commit()
        self.numerical_source.write_bytes(qualified)
        with self.assertRaisesRegex(ValueError, "differs from qualified numerical code"):
            self.prepare(sha)
        self.assertFalse(self.output.exists())

    def test_uncommitted_numerical_files_are_not_dispatch_ready(self):
        sha = self.fixture_commit(add_sources=False)
        with self.assertRaisesRegex(ValueError, "numerical file set"):
            self.prepare(sha)
        self.assertFalse(self.output.exists())

    def test_unknown_commit_is_rejected(self):
        self.fixture_commit()
        with self.assertRaisesRegex(ValueError, "unavailable"):
            self.prepare("f" * 40)
        self.assertFalse(self.output.exists())

    def test_case_name_cannot_address_the_parent_results_directory(self):
        for name in (".", "..", "../previous_case"):
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, "invalid case name"):
                inputs.prepare(self.phase, [self.replay], self.output, name, source_root=self.checkout)
        self.assertFalse(self.output.exists())

    def test_passed_flag_cannot_hide_numeric_failure(self):
        self.comparisons[0]["qualified_reference"] = .5
        self.save()
        with self.assertRaisesRegex(ValueError, "numerical acceptance"):
            self.prepare()

    def test_stale_recorded_error_is_rejected(self):
        self.comparisons[0]["replay_relative_change"] = .001
        self.save()
        with self.assertRaisesRegex(ValueError, "actual values"):
            self.prepare()

    def test_missing_or_renamed_observable_is_rejected(self):
        self.comparisons[0]["observable"] = "eta"
        self.save()
        with self.assertRaisesRegex(ValueError, "invalid replay"):
            self.prepare()

    def test_nonfinite_values_are_rejected(self):
        self.comparisons[0]["production"] = "NaN"
        self.save()
        with self.assertRaisesRegex(ValueError, "observable"):
            self.prepare()

    def test_failed_anchor_or_background_is_not_promoted(self):
        self.anchors[0]["refinement_delta_MeV"] = .2
        self.save()
        with self.assertRaisesRegex(ValueError, "anchor"):
            self.prepare()
        self.anchors[0]["refinement_delta_MeV"] = 0
        self.background[0]["entropy_relative_change"] = .02
        self.save()
        with self.assertRaisesRegex(ValueError, "background numerical"):
            self.prepare()

    def test_incomplete_background_and_duplicate_plan_are_rejected(self):
        self.background.pop()
        self.save()
        with self.assertRaisesRegex(ValueError, "cover the plan"):
            self.prepare()
        self.plan.append(dict(self.plan[0]))
        self.save()
        with self.assertRaisesRegex(ValueError, "duplicate phase plan"):
            self.prepare()

    def test_existing_output_is_preserved(self):
        self.output.mkdir()
        sentinel = self.output / "keep.txt"
        sentinel.write_text("existing evidence")
        with self.assertRaisesRegex(ValueError, "new empty"):
            self.prepare()
        self.assertEqual(sentinel.read_text(), "existing evidence")

    def test_different_profile_cannot_be_combined_with_phase_qualification(self):
        manifest = json.loads((self.replay / "manifest.json").read_text())
        manifest["parameters"]["profile"]["source_sha256"] = "c" * 64
        (self.replay / "manifest.json").write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError, "different sampling profiles"):
            self.prepare()

    def test_renamed_points_cannot_replace_the_qualified_set(self):
        for row in self.comparisons:
            row["point_id"] = "unqualified-" + row["point_id"]
        self.save()
        with self.assertRaisesRegex(ValueError, "qualified reference"):
            self.prepare()

    def test_changed_recorded_targets_cannot_hide_a_different_comparison(self):
        for key in ("production", "qualified_candidate", "qualified_reference"):
            self.comparisons[0][key] = 1.001
        self.save()
        with self.assertRaisesRegex(ValueError, "qualified reference"):
            self.prepare()

    def test_changed_numerical_source_requires_requalification(self):
        self.numerical_source.write_text("changed numerical implementation\n")
        with self.assertRaisesRegex(ValueError, "differs from current code"):
            self.prepare()

    def test_changed_snapshot_is_rejected(self):
        (self.replay / "source" / "src" / "kernel.jl").write_text("changed snapshot\n")
        with self.assertRaisesRegex(ValueError, "source snapshot changed"):
            self.prepare()

    def test_missing_source_identity_is_rejected(self):
        manifest = json.loads((self.replay / "manifest.json").read_text())
        manifest["source_files"] = []
        (self.replay / "manifest.json").write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError, "incomplete numerical source"):
            self.prepare()


if __name__ == "__main__":
    unittest.main()
