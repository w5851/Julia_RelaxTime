"""Exact finite-difference attribution and cancellation cases (run remotely)."""
import importlib.util
import math
from pathlib import Path
import unittest
import os

SCRIPT = (Path(os.environ["MOTT_SLOPE_REVIEW_SCRIPT"]) if "MOTT_SLOPE_REVIEW_SCRIPT" in os.environ
          else Path(__file__).resolve().parents[3] / "scripts/analysis/relaxtime/review_transport_other_mott_slopes.py")
spec = importlib.util.spec_from_file_location("other_mott_slopes", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class SlopeAttributionTests(unittest.TestCase):
    def test_precision_comparison_keeps_native_edge_with_extra_reference_node(self):
        edge = dict(curve_id="toy",point_left="a",point_right="b",observable="tau_s",
                    xi_left=0.,xi_right=.002,interval_width=.002,log_slope=2.)
        refs = {"a":dict(xi=0.,tau_s=1.),"m":dict(xi=.001,tau_s=math.exp(.001)),
                "b":dict(xi=.002,tau_s=math.exp(.0044))}
        checks = module.compare_native_slopes([edge],refs)
        self.assertEqual(len(checks),1)
        self.assertEqual((checks[0]["point_left"],checks[0]["point_right"]),("a","b"))
        self.assertAlmostEqual(checks[0]["tighter_log_slope"],2.2,places=10)
        self.assertAlmostEqual(checks[0]["absolute_difference"],.2,places=10)
        # A sparse reference missing a native endpoint must not substitute a wider edge.
        del refs["b"]
        self.assertEqual(module.compare_native_slopes([edge],refs),[])

    def test_precision_comparison_rejects_shifted_endpoint(self):
        edge = dict(curve_id="toy",point_left="a",point_right="b",observable="tau_s",
                    xi_left=0.,xi_right=.002,interval_width=.002,log_slope=2.)
        with self.assertRaises(ValueError):
            module.compare_native_slopes([edge],{"a":dict(xi=0.,tau_s=1.),"b":dict(xi=.0021,tau_s=1.01)})

    def test_logmean_near_equal_and_reversal(self):
        for a, b in ((3., 3.), (3., 3.+1e-12), (.1, 10.), (10., .1)):
            value = module.logmean(a, b)
            self.assertTrue(min(a, b) <= value <= max(a, b))
            self.assertAlmostEqual(value, module.logmean(b, a), places=12)
        self.assertAlmostEqual(module.log_slope(3, 3, .001), 0)
        with self.assertRaises(ValueError):
            module.log_slope(1, 2, 0)

    def test_product_competing_changes(self):
        # The density increases while the rate decreases; the product is fixed.
        rate_term, density_term = module.product_difference(2, 4, 6, 3)
        self.assertEqual(rate_term, -9)
        self.assertEqual(density_term, 9)
        self.assertEqual(rate_term + density_term, 0)

    def test_full_species_and_channel_closure(self):
        observables, rates, weights, channels = [], [], [], []
        names = ["usbar_to_usbar", "subar_to_subar"] + [f"process_{k}" for k in range(15)]
        for index, x in enumerate((-.002, 0., .002)):
            pid = f"p{index}"
            row = dict(point_id=pid, curve_id="toy", meson="K", window="toy", T_MeV=200,
                       muB_MeV=0, xi=x, xi_mott=0)
            for k, name in enumerate(names):
                rate = (1+k/10)*(1+(k%3)*x)
                rates.append(dict(point_id=pid, channel=name, estimated_error=rate*1e-7))
            for j, sp in enumerate(module.SPECIES):
                total = 0
                for k, name in enumerate(names):
                    density = (1+j/10)*(1+(k%4)*x)
                    rate = (1+k/10)*(1+(k%3)*x)
                    contribution = density*rate
                    total += contribution
                    channels.append(dict(point_id=pid, species=sp, channel=name, rate_fm2=rate,
                        collision_density_weight_inv_fm3=density,
                        inverse_tau_contribution_inv_fm=contribution))
                row["tauinv_"+sp] = total
                row["tau_"+sp] = 1/total
                item = dict(point_id=pid, species=sp)
                for k, final in enumerate(module.FINAL):
                    item[final+"_coefficient"] = (1+j/4)*(1+(k+1)*x)
                weights.append(item)
            for final in module.FINAL:
                row[final] = sum(w[final+"_coefficient"]*row["tau_"+w["species"]]
                                 for w in weights if w["point_id"] == pid)
            observables.append(row)
        result, error = module.analyze(observables, rates, weights, channels)
        self.assertLess(error, 1e-10)
        self.assertEqual(len(result["slopes"]), 18)
        self.assertEqual(len(result["channel_derivatives"]), 204)
        self.assertEqual(len(result["final_species_derivatives"]), 36)
        a, b = observables[:2]
        for row in result["slopes"][:9]:
            expected = math.log(b[row["observable"]]/a[row["observable"]])/.002
            self.assertAlmostEqual(row["log_slope"], expected, places=10)
            self.assertGreater(row["estimated_s_interpolation_slope_error"], 0)


if __name__ == "__main__":
    unittest.main()
