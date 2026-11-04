import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location("transport_artifact", ROOT / "scripts/relaxtime/check_transport_production_artifact.py")
artifact = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(artifact)


class ArtifactCheckTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.profile = self.folder / "profile.toml"
        self.profile.write_text('algorithm="thermal_weighted"\n', encoding="utf-8")
        self.sha = "a" * 40
        self.point = dict(T_MeV="200", muB_MeV="0", xi="0")
        self.rows = [dict(self.point, run_id="second", converged="true",
                          **{k: "1.0" for k in artifact.OBSERVABLES})]
        self.diagnostics = [dict(self.point, channel=k, run_id="second", rate="1", estimated_error="0.0001", converged="true")
                            for k in sorted(artifact.RATE_KEYS)]
        self.cfg = {"dry_run": False, "compute_bulk": True, "propagator_xi_policy": "match_thermo",
                    "sampling_profile": {"algorithm": "thermal_weighted", "id": "thermal_weighted_v1",
                    "source_sha256": hashlib.sha256(self.profile.read_bytes()).hexdigest(),
                    "sampling": {"atol": 1e-12, "rtol": 1e-3}}}
        self.manifest = {"git_commit": self.sha, "summary": {"error_count": 0, "success_count": 1, "skipped_count": 0}}
        self.source = self.folder / "numerical.jl"
        self.source.write_text("numerical implementation\n", encoding="utf-8")
        contract = {"identity": {"profile_id": "thermal_weighted_v1",
                    "profile_sha256": self.cfg["sampling_profile"]["source_sha256"],
                    "propagator_xi_policy": "match_thermo", "compute_bulk": True,
                    "sources": [{"path": "numerical.jl", "sha256": hashlib.sha256(self.source.read_bytes()).hexdigest()}]}}
        (self.folder / "phase_guided_transport_scan_sampling_contract.json").write_text(json.dumps(contract), encoding="utf-8")
        self.write_csv("sampling_plan.csv", [self.point])
        self.write_csv("failed_points.csv", [], fields=list(self.point))
        self.save()

    def write_csv(self, name, rows, fields=None):
        with (self.folder / name).open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields or list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    def save(self):
        (self.folder / "run_manifest.json").write_text(json.dumps(self.manifest), encoding="utf-8")
        (self.folder / "effective_config.json").write_text(json.dumps({"options": self.cfg}), encoding="utf-8")
        self.write_csv("phase_guided_transport_scan.csv", self.rows, fields=list(self.rows[0]) if self.rows else [*self.point, "run_id", *artifact.OBSERVABLES])
        self.write_csv("phase_guided_transport_scan_sampling.csv", self.diagnostics)

    def check(self, tier="production"):
        return artifact.check_artifact(self.folder, self.profile, self.sha, tier, source_root=self.folder)

    def test_complete_artifact(self):
        self.assertEqual(self.check()["qualified_channels"], 17)

    def test_dry_run_does_not_claim_numerical_completion(self):
        self.cfg["dry_run"] = True
        self.manifest["summary"] = {"dry_run": True}
        self.save()
        self.assertEqual(self.check("dry-run"), {"dry_run": True, "planned_points": 1})
        with self.assertRaises(ValueError):
            self.check()

    def test_missing_final_point(self):
        self.rows = []
        self.save()
        with self.assertRaisesRegex(ValueError, "point set"):
            self.check()

    def test_duplicate_channel(self):
        self.diagnostics.append(dict(self.diagnostics[0]))
        self.save()
        with self.assertRaisesRegex(ValueError, "channel diagnostics"):
            self.check()

    def test_exhausted_budget(self):
        self.diagnostics[0]["estimated_error"] = "0.01"
        self.save()
        with self.assertRaisesRegex(ValueError, "channel budget"):
            self.check()

    def test_nonfinite_final_observable(self):
        self.rows[0]["zeta_over_s"] = "nan"
        self.save()
        with self.assertRaisesRegex(ValueError, "invalid final observable"):
            self.check()

    def test_success_summary_cannot_hide_unconverged_result(self):
        for value in ("false", "", None):
            with self.subTest(converged=value):
                if value is None:
                    self.rows[0].pop("converged", None)
                else:
                    self.rows[0]["converged"] = value
                self.save()
                with self.assertRaisesRegex(ValueError, "unconverged result"):
                    self.check()

    def test_resume_keeps_failed_attempt_evidence(self):
        self.diagnostics.append(dict(self.diagnostics[0], run_id="first", converged="false"))
        self.write_csv("failed_points.csv", [self.point])
        self.manifest["summary"] = {"error_count": 0, "success_count": 0, "skipped_count": 1}
        self.save()
        self.assertEqual(self.check()["completed_points"], 1)

    def test_source_and_profile_identity(self):
        self.manifest["git_commit"] = "b" * 40
        self.save()
        with self.assertRaisesRegex(ValueError, "commit"):
            self.check()
        self.manifest["git_commit"] = self.sha
        self.save()
        self.profile.write_text('algorithm="other"\n', encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "sampling profile"):
            self.check()

    def test_same_commit_cannot_hide_changed_numerical_source(self):
        self.source.write_text("changed implementation\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "numerical source"):
            self.check()

    def test_physical_configuration_must_match_the_result_contract(self):
        self.cfg["propagator_xi_policy"] = "isotropic"
        self.save()
        with self.assertRaisesRegex(ValueError, "physical contract"):
            self.check()


if __name__ == "__main__":
    unittest.main()
