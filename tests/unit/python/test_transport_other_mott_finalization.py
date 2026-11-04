"""Fast guards for the separate complete-reference postprocessing stage."""
from pathlib import Path
import importlib.util
import unittest

candidates = [parent / "scripts/analysis/relaxtime/finalize_transport_other_mott_acceptance.py"
    for parent in Path(__file__).resolve().parents]
candidates.append(Path(__file__).resolve().parent / "finalize_transport_other_mott_acceptance_v1.py")
path = next((candidate for candidate in candidates if candidate.exists()), None)
if path is None:
    raise FileNotFoundError("Other-Mott finalization script was not found")
spec = importlib.util.spec_from_file_location("other_mott_finalization", path)
finalization = importlib.util.module_from_spec(spec)
spec.loader.exec_module(finalization)


class FinalizationGuardTests(unittest.TestCase):
    def ready(self):
        return dict(status="ready_for_agent_review", pending_jobs=0, active_jobs=[], external_jobs=[],
            failed_runs=[], reference_runs=["tighter_a", "tighter_b"], finished_runs=["tighter_a", "tighter_b"])

    def complete(self):
        return dict(status="completed_diagnostic_review", main_grid_complete=True, points=195, curves=13,
            channels=3315, source_files_checked=340, source_all_match=True, reference_evaluations=175,
            expected_reference_evaluations=175, reference_plan_complete=True, whole_task_complete=False,
            missing_reference_ids=[], missing_xi_intervals=[])

    def test_running_queue_waits(self):
        self.assertFalse(finalization.ready_state(dict(status="running_reference_recovery")))
        self.assertTrue(finalization.ready_state(self.ready()))

    def test_failure_does_not_start_audit(self):
        with self.assertRaises(ValueError):
            finalization.ready_state(dict(status="failed"))
        state = self.ready()
        state["failed_runs"] = ["tighter_a"]
        with self.assertRaises(ValueError):
            finalization.ready_state(state)

    def test_ready_requires_idle_unique_matching_runs(self):
        for change in (dict(pending_jobs=1), dict(active_jobs=[dict(pid=123)]),
                       dict(external_jobs=[dict(name="tighter_c")]),
                       dict(reference_runs=["tighter_a", "tighter_a"]), dict(finished_runs=["tighter_a"])):
            with self.subTest(change=change), self.assertRaises(ValueError):
                finalization.ready_state(dict(self.ready(), **change))

    def test_every_reference_and_split_is_required(self):
        finalization.validate_complete_summary(self.complete(), 175)
        for change in (dict(reference_evaluations=174), dict(reference_plan_complete=False),
                       dict(missing_reference_ids=["p_missing"]), dict(missing_xi_intervals=[dict(left_id="p1")])):
            with self.subTest(change=change), self.assertRaises(ValueError):
                finalization.validate_complete_summary(dict(self.complete(), **change), 175)

    def test_complete_grid_does_not_claim_scientific_completion(self):
        with self.assertRaises(ValueError):
            finalization.validate_complete_summary(dict(self.complete(), whole_task_complete=True), 175)


if __name__ == "__main__":
    unittest.main()
