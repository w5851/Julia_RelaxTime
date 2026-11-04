"""Partial precision coverage must remain explicit in diagnostic acceptance."""
import importlib.util
import math
import os
from pathlib import Path
import unittest

SCRIPT = (Path(os.environ["MOTT_ACCEPTANCE_SCRIPT"]) if "MOTT_ACCEPTANCE_SCRIPT" in os.environ
          else Path(__file__).resolve().parents[3]/"scripts/analysis/relaxtime/review_transport_other_mott_acceptance.py")
spec = importlib.util.spec_from_file_location("mott_acceptance", SCRIPT)
review = importlib.util.module_from_spec(spec)
spec.loader.exec_module(review)


class AcceptanceTests(unittest.TestCase):
    def test_partial_xi_edges_are_reported_and_not_fabricated(self):
        edges = [dict(curve_id="a", left_id="l", midpoint_id="m", right_id="r"),
                 dict(curve_id="b", left_id="x", midpoint_id="y", right_id="z")]
        refs = {name: dict(xi=x, y=math.exp(x*x), curve_id="a") for name, x in (("l",0.),("m",.5),("r",1.))}
        rows, missing = review.xi_comparisons(edges, refs, ["y"], lambda a,b,dx: math.log(b/a)/dx)
        self.assertEqual(len(rows), 1)
        self.assertEqual(missing[0]["missing_ids"], ["x", "y", "z"])
        self.assertAlmostEqual(rows[0]["coarse_log_slope"], 1.)
        self.assertAlmostEqual(rows[0]["refined_left_log_slope"], .5)
        self.assertAlmostEqual(rows[0]["refined_right_log_slope"], 1.5)

    def test_cross_curve_refinement_is_rejected(self):
        refs = {name: dict(xi=x, y=1., curve_id=curve) for name, x, curve in
                (("l",0.,"a"),("m",.5,"b"),("r",1.,"a"))}
        with self.assertRaisesRegex(ValueError, "geometry"):
            review.xi_comparisons([dict(curve_id="a",left_id="l",midpoint_id="m",right_id="r")],
                                 refs, ["y"], lambda a,b,dx: 0.)

    def test_missing_and_shifted_grid_points_are_rejected(self):
        point = dict(point_id="p", curve_id="a", window="w", meson="K", xi=0., xi_mott=0.,
                     T_MeV=200., muB_MeV=0., background_residual=0.)
        self.assertEqual(review.check_grid([], [point], False), ["p"])
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            review.check_grid([], [point], True)
        with self.assertRaisesRegex(ValueError, "coordinate"):
            review.check_grid([dict(point, xi=.01)], [point], True)
        with self.assertRaisesRegex(ValueError, "coordinate"):
            review.check_grid([dict(point, xi=math.nan)], [point], True)


if __name__ == "__main__":
    unittest.main()
