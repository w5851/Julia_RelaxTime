import copy
import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location(
    "mode_b_inputs", ROOT / "scripts/analysis/relaxtime/prepare_transport_mode_b_qualification.py")
inputs = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(inputs)


class ModeBQualificationTests(unittest.TestCase):
    def setUp(self):
        self.plan = [dict(T_MeV=t, muB_MeV=mu, xi=xi, mode="mode_b_fixed_T_sparse_muB")
                     for t, mu, xi in sorted(inputs.expected_keys())]
        self.background = [dict(row, passed="true", coexistence_resolved="true", m_u=1, m_d=1, m_s=2,
                                entropy=.2, max_mass_relative_change=1e-6, entropy_relative_change=1e-5,
                                distinct_branches="false", refined_distinct_branches="false",
                                stable_phase="hadron", refined_stable_phase="hadron") for row in self.plan]
        self.status = dict(complete=True, phase_plan_passed=True, cold_seed_audits_passed=True,
                           failures=0, points=909, completed_points=909)

    def test_all_curves_extremes_and_isotropic_controls_are_selected(self):
        inputs.validate_background(self.plan, self.background, self.status)
        selected = inputs.select_points(self.plan, [])
        self.assertEqual(len(selected), 15)
        self.assertEqual({p[:2] for p in selected}, {(t, mu) for t in (120, 160, 200) for mu in (0, 450, 900)})
        self.assertTrue(all((t, mu, -.5) in selected for t in (120, 160, 200) for mu in (0, 450, 900)))
        self.assertEqual(sum(p[2] == 0 for p in selected), 3)

    def test_phase_bracket_endpoints_are_added_without_moving_them(self):
        boundary = dict(T_MeV=120, muB_MeV=900, xi=-.123, refined_xi=-.12299,
                        bracket_lower_xi=-.13, bracket_upper_xi=-.12)
        selected = inputs.select_points(self.plan, [boundary])
        self.assertEqual(len(selected), 17)
        self.assertIn((120, 900, -.13), selected)
        self.assertIn((120, 900, -.12), selected)
        with self.assertRaisesRegex(ValueError, "duplicated"):
            inputs.select_points(self.plan, [boundary, boundary])
        with self.assertRaisesRegex(ValueError, "boundary bracket"):
            inputs.select_points(self.plan, [dict(boundary, refined_xi=.1)])

    def test_missing_and_duplicate_grid_points_are_rejected(self):
        for plan in (self.plan[:-1], self.plan[:-1] + self.plan[:1]):
            with self.subTest(length=len(plan)), self.assertRaisesRegex(ValueError, "909-point"):
                inputs.validate_background(plan, self.background, self.status)
        selected_plan = [row for row in self.plan if inputs.key(row) != (120, 0, -.5)]
        with self.assertRaisesRegex(ValueError, "point is missing"):
            inputs.select_points(selected_plan, [])

    def test_unresolved_or_failed_rows_are_not_hidden_by_status(self):
        for change in (dict(passed="false"), dict(coexistence_resolved="false"),
                       dict(entropy=float("nan")), dict(max_mass_relative_change=.02),
                       dict(distinct_branches="true", refined_distinct_branches="true", refined_stable_phase="quark")):
            background = copy.deepcopy(self.background)
            background[0].update(change)
            with self.subTest(change=change), self.assertRaises(ValueError):
                inputs.validate_background(self.plan, background, self.status)

    def test_incomplete_audit_status_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "incomplete"):
            inputs.validate_background(self.plan, self.background, dict(self.status, cold_seed_audits_passed=False))

    def test_stable_seed_uses_omega_and_exact_qualified_masses(self):
        high = dict(m_u=2, m_d=2, m_s=3, omega=-20, Phi=.1, Phibar=.2)
        low = dict(m_u=.1, m_d=.1, m_s=1, omega=-21, Phi=.3, Phibar=.4)
        self.assertEqual(inputs.stable_trace(high, low, low), low)
        self.assertEqual(inputs.stable_trace(low, high, low), low)
        with self.assertRaisesRegex(ValueError, "qualified background masses"):
            inputs.stable_trace(high, low, high)
        with self.assertRaisesRegex(ValueError, "nonfinite"):
            inputs.stable_trace(high, dict(low, Phi=float("nan")), low)


if __name__ == "__main__":
    unittest.main()
