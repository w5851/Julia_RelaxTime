import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[3]
SPEC = importlib.util.spec_from_file_location("transport_full_grid", ROOT / "scripts/analysis/relaxtime/review_transport_full_grid.py")
grid = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(grid)


def job(xis="-0.003,0.003", temperature=125.736):
    return dict(expected_points=len(xis.split(",")), expected_T_MeV=temperature,
                inputs=dict(muB_list="900", alpha_t_list="1", t_list="125.736", xi_list=xis))


def row(xi):
    return dict(T_MeV="125.736", muB_MeV="900", alpha_T="1", xi=str(xi),
                quality_flag="false", quality_reason="ok", quality_metric="1", tau_u="1", tau_ubar="1",
                residual_norm="1e-10")


class TransportFullGridTests(unittest.TestCase):
    def setUp(self):
        self.job = job()
        self.entry = dict(mode="A", job=self.job)
        self.bound = dict(A=dict(full_points=2, provisional_full_jobs=[self.job]))
        self.expected = grid.coverage(self.bound, [self.entry])

    def test_reused_and_new_shards_cover_each_special_boundary_point_once(self):
        entries = [dict(mode="A", job=job("-0.003")), dict(mode="A", job=job("0.003"))]
        self.assertEqual(grid.coverage(self.bound, entries), self.expected)
        with self.assertRaisesRegex(ValueError, "overlap"):
            grid.coverage(self.bound, entries + [entries[0]])
        with self.assertRaisesRegex(ValueError, "cover"):
            grid.coverage(self.bound, entries[:1])
        self.assertEqual(grid.coverage(self.bound, entries[:1], partial=True), self.expected)

    def test_equal_counts_do_not_hide_wrong_grid_or_anchor(self):
        for changed in (job("-0.01,0.01"), job(temperature=125.7)):
            with self.subTest(job=changed), self.assertRaisesRegex(ValueError, "changes"):
                grid.coverage(self.bound, [dict(mode="A", job=changed)])

    def test_modes_have_separate_curve_identity(self):
        bound = dict(self.bound, B=dict(full_points=2, provisional_full_jobs=[self.job]))
        expected = grid.coverage(bound, [self.entry, dict(mode="B", job=self.job)])
        self.assertEqual(len(expected), 4)

    def test_extra_diagnostics_or_duplicate_actual_rows_are_rejected(self):
        rows = [row(-0.003), row(0.003)]
        diagnostics = [{}] * 34
        keys = grid.verify_rows(rows, diagnostics, self.entry, self.expected, set())
        self.assertEqual(keys, set(self.expected))
        with self.assertRaisesRegex(ValueError, "diagnostics"):
            grid.verify_rows(rows, diagnostics + [{}], self.entry, self.expected, set())
        with self.assertRaisesRegex(ValueError, "duplicate"):
            grid.verify_rows(rows, diagnostics, self.entry, self.expected, keys)
        with self.assertRaisesRegex(ValueError, "duplicate"):
            grid.verify_rows([rows[0], rows[0]], diagnostics, self.entry, self.expected, set())

    def test_background_and_quality_failures_do_not_certify_a_grid(self):
        for changes in (dict(residual_norm="nan"), dict(residual_norm="1e-7"),
                        dict(quality_flag="true"), dict(T_MeV="125.73"), dict(alpha_T="1.1")):
            rows = [dict(row(-0.003), **changes), row(0.003)]
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                grid.verify_rows(rows, [{}] * 34, self.entry, self.expected, set())

    def test_finite_high_density_asymmetry_flag_is_preserved_as_a_heuristic(self):
        flagged = dict(row(-0.003), tau_u="7", quality_flag="true",
                       quality_reason="tau_u_ubar_ratio_high", quality_metric="7")
        keys = grid.verify_rows([flagged, row(0.003)], [{}] * 34, self.entry, self.expected, set())
        self.assertEqual(keys, set(self.expected))
        for change in (dict(quality_flag="false"), dict(quality_metric="6.5"), dict(tau_u="nan")):
            with self.subTest(change=change), self.assertRaises(ValueError):
                grid.verify_quality_flag(dict(flagged, **change))


if __name__ == "__main__":
    unittest.main()
