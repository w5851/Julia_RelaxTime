from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts/analysis/relaxtime"))
try:
    import audit_transport_grid_background as audit
finally:
    sys.path.pop(0)


def background(phase="hadron", distinct="true"):
    return dict(T_MeV="120", m_u="2", m_d="2", m_s="3", entropy="0.1",
                distinct_branches=distinct, stable_phase=phase)


def result(phase="hadron"):
    return dict(T_MeV="120", muB_MeV="900", xi="0", m_u="2", m_d="2", m_s="3",
                s_fm3inv="0.1", phase_curr=phase, **{key: "1" for key in audit.OBSERVABLES})


class TransportGridBackgroundTests(unittest.TestCase):
    def test_replay_keeps_existing_mass_entropy_and_distinct_branch_contract(self):
        ref, row = background(), result()
        self.assertTrue(audit.compare_background(row, ref, "B")["passed"])
        self.assertFalse(audit.compare_background(dict(row, m_u="2.03"), ref, "B")["passed"])
        self.assertFalse(audit.compare_background(dict(row, s_fm3inv="0.102"), ref, "B")["passed"])
        self.assertFalse(audit.compare_background(result("quark"), ref, "B")["passed"])
        self.assertTrue(audit.compare_background(result("quark"), background(distinct="false"), "B")["passed"])

    def test_partial_review_does_not_bridge_missing_xi_points(self):
        keys = [("B", 900.0, 120.0, xi) for xi in (-0.01, 0.0, 0.01)]
        refs = {key: background() for key in keys}
        partial = {keys[0]: result(), keys[-1]: result()}
        self.assertEqual(audit.neighbouring_changes(partial, refs, refs), [])
        full = {key: result() for key in keys}
        self.assertEqual(len(audit.neighbouring_changes(full, refs, refs)), 18)

    def test_known_boundary_changes_are_retained_and_identified(self):
        left, right = ("B", 900.0, 120.0, -0.13), ("B", 900.0, 120.0, -0.12)
        refs = {left: background("hadron"), right: background("quark")}
        rows = {left: result(), right: dict(result("quark"), tau_u="4")}
        changes = audit.neighbouring_changes(rows, refs, refs)
        self.assertEqual(len(changes), 9)
        selected = [change for change in changes if change["selected_for_review"]]
        self.assertEqual(len(selected), 1)
        self.assertTrue(selected[0]["reference_first_order_boundary"])
        self.assertEqual(selected[0]["relative_to_larger"], 0.75)


if __name__ == "__main__":
    unittest.main()
