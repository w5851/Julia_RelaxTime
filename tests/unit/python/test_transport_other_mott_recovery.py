"""Recovery preserves completed work and never changes the precision contract."""
import copy
import importlib.util
import os
from pathlib import Path
import unittest

SCRIPT = (Path(os.environ["MOTT_RECOVERY_SCRIPT"]) if "MOTT_RECOVERY_SCRIPT" in os.environ
          else Path(__file__).resolve().parents[3]/
          "scripts/analysis/relaxtime/resume_transport_other_mott_pipeline.py")
spec = importlib.util.spec_from_file_location("mott_recovery", SCRIPT)
recovery = importlib.util.module_from_spec(spec)
spec.loader.exec_module(recovery)


class RecoveryPlanTests(unittest.TestCase):
    def jobs(self):
        return [dict(name="tighter_"+c+"_v1", curve_id=c, point_ids=points) for c, points in
                (("other01", ["pilot", "a"]), ("other03", ["b"]),
                 ("other05", ["other05_p02", "c"]), ("other07", ["d"]))]

    def test_keeps_completed_and_running_and_retries_failed_once(self):
        states = dict(tighter_other01_v1=dict(status="completed_diagnostic"),
            tighter_other03_v1=dict(status="running", point_ids=["b"]),
            tighter_other05_v1=dict(status="failed"))
        result = recovery.reference_recovery_plan(self.jobs(), {"pilot"}, states, {"tighter_other01_v1": ["a"]})
        self.assertEqual(result["reused"], ["tighter_pilot_other01_v1", "tighter_other01_v1"])
        self.assertEqual(result["external"], [dict(name="tighter_other03_v1", point_ids=["b"])])
        self.assertEqual(result["jobs"][0]["point_ids"], ["other05_p02"])
        self.assertEqual(result["jobs"][1]["curve_id"], "other07")
        pending = [p for j in result["jobs"] for p in j["point_ids"]]
        self.assertCountEqual(pending, ["other05_p02", "c", "d"])
        self.assertEqual(len(pending), len(set(pending)))

    def test_incomplete_success_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "coverage mismatch"):
            recovery.reference_recovery_plan(self.jobs(), {"pilot"},
                {"tighter_other01_v1": dict(status="completed_diagnostic")}, {"tighter_other01_v1": []})

    def test_partial_failed_output_needs_explicit_salvage(self):
        with self.assertRaisesRegex(ValueError, "explicit salvage"):
            recovery.reference_recovery_plan(self.jobs(), {"pilot"},
                {"tighter_other05_v1": dict(status="failed")}, {"tighter_other05_v1": ["c"]})

    def test_duplicate_plan_ids_are_rejected(self):
        jobs = self.jobs()
        jobs[-1]["point_ids"].append("c")
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            recovery.reference_recovery_plan(jobs, {"pilot"}, {}, {})

    def test_larger_caps_preserve_numerical_configuration(self):
        original = dict(id="old", algorithm="thermal_weighted", sampling=dict(
            rtol=1e-6, atol=1e-12, seed_order=48, max_evaluations=80000, max_nodes=24000))
        updated = copy.deepcopy(original)
        updated["id"] = "new"
        updated["sampling"].update(max_evaluations=320000, max_nodes=96000)
        recovery.check_budget_only_change(original, updated)
        for field, value in (("rtol", 1e-5), ("seed_order", 64), ("max_nodes", 12000)):
            bad = copy.deepcopy(updated)
            bad["sampling"][field] = value
            with self.assertRaises(ValueError):
                recovery.check_budget_only_change(original, bad)


if __name__ == "__main__":
    unittest.main()
