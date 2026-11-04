"""Side-local differentiation: known powers, stencil integrity, and noise gain."""
from pathlib import Path
import importlib.util
import math
import os
import unittest

configured = os.environ.get("SIGMA_ASYM_REVIEW_SCRIPT")
script = Path(configured) if configured else Path(__file__).resolve().parents[3]/ \
    "scripts/analysis/relaxtime/review_sigma_mott_asymptotics.py"
spec=importlib.util.spec_from_file_location("sigma_asym",script)
review=importlib.util.module_from_spec(spec)
spec.loader.exec_module(review)


def synthetic(function,steps,error=0.,factors=(.75,.875,1.,1.125,1.25)):
    coords=sorted({side*h*f for side in (-1,1) for h in steps for f in factors})
    return [dict(point_id=str(i),delta_xi=x,mixed_square_rate_fm2=function(x),quadrature_error_fm2=error)
            for i,x in enumerate(coords)]


class SideCurvatureTests(unittest.TestCase):
    def test_quadratic_both_sides(self):
        steps=[.01,.005]
        for row in review.derive(synthetic(lambda x:2+3*x+7*x*x,steps),steps):
            self.assertAlmostEqual(row["curvature"],14.,delta=1e-8)

    def test_three_halves_only_on_resonant_side(self):
        steps=[.01,.005,.0025]
        result=review.derive(synthetic(lambda x:max(-x,0.)**1.5,steps),steps)
        for row in result:
            if row["side"]=="left":
                # Fixed relative steps retain a known O((step/h)^4) bias.
                # The independently evaluated power stencil is 0.749986170506,
                # not exactly the pointwise coefficient 3/4.
                self.assertAlmostEqual(row["scaled_curvature"],.749986170506,delta=1e-10)
                if row["effective_exponent"] is not None:
                    self.assertAlmostEqual(row["effective_exponent"],.5,delta=1e-10)
            else:
                self.assertEqual(row["curvature"],0.)

    def test_missing_point_rejected(self):
        rows=synthetic(lambda x:x*x,[.01])
        rows=[r for r in rows if r["delta_xi"]!=.01]
        with self.assertRaises(ValueError):
            review.second_difference(rows,.01,1,.25)

    def test_point_error_amplifies_as_inverse_step_squared(self):
        error=1e-13
        rows=synthetic(lambda x:x*x,[.01],error)
        value=review.second_difference(rows,.01,1,.125)
        self.assertAlmostEqual(value["quadrature_defect"],4*error/(.00125**2),delta=1e-20)

    def test_logarithmic_correction_has_zero_power_exponent(self):
        steps=[2.**(-j) for j in (14,15,16)]
        result=review.derive(synthetic(lambda x:x*x*math.log(abs(x)),steps),steps)
        for row in result:
            if row["effective_exponent"] is not None:
                self.assertLess(row["effective_exponent"],.15)
                self.assertGreater(row["effective_exponent"],0.)

    def test_smaller_stencil_reduces_known_power_bias(self):
        steps=[.01,.005]
        fn=lambda x:max(-x,0.)**1.5
        coarse=review.derive(synthetic(fn,steps),steps)
        fine=review.derive(synthetic(fn,steps,factors=(.875,.9375,1.,1.0625,1.125)),
            steps,fractions=(.125,.0625))
        for a,b in zip(coarse,fine):
            if a["side"]=="left":
                self.assertLess(abs(b["scaled_curvature"]-.75),abs(a["scaled_curvature"]-.75)/14)
                self.assertAlmostEqual(b["scaled_curvature"],.75,delta=1e-6)

    def test_step_that_crosses_mott_is_rejected(self):
        with self.assertRaises(ValueError):
            review.derive([], [.01], fractions=(1.2,.6))


if __name__=="__main__":
    unittest.main()
