"""Distinguish finite-window curvature from a persistent first-derivative jump."""
from pathlib import Path
import importlib.util
import math
import os
import unittest
from unittest.mock import patch

configured_script = os.environ.get("SIGMA_REVIEW_SCRIPT")
script = (Path(configured_script) if configured_script else
          Path(__file__).resolve().parents[3]/"scripts/analysis/relaxtime/review_sigma_mott_refinement.py")
spec = importlib.util.spec_from_file_location("sigma_review", script)
review = importlib.util.module_from_spec(spec)
spec.loader.exec_module(review)

configured_finalizer = os.environ.get("SIGMA_FINALIZE_SCRIPT")
finalizer_script = (Path(configured_finalizer) if configured_finalizer else
                   script.with_name("finalize_sigma_mott_checks.py"))
finalizer_spec = importlib.util.spec_from_file_location("sigma_finalizer", finalizer_script)
finalizer = importlib.util.module_from_spec(finalizer_spec)
finalizer_spec.loader.exec_module(finalizer)


def points(log_function):
    center = 0.015327819995582102
    result = {}
    for offset, tag in ((-2e-4,"m200"),(-1e-4,"m100"),(-5e-5,"m050"),(-2e-5,"m020"),
                        (0.,"c000"),(2e-5,"p020"),(5e-5,"p050"),(1e-4,"p100"),(2e-4,"p200")):
        result["sig_"+tag] = dict(xi=center+offset,
            **{key: math.exp(log_function(offset)) for key in review.OBSERVABLES})
    return result


class SigmaCenterDiagnosticsTests(unittest.TestCase):
    def test_smooth_curvature_jump_shrinks_with_width(self):
        rows = review.centered_slopes("smooth", points(lambda x: 3*x+7*x*x))
        self.assertEqual(len(rows), 4*len(review.OBSERVABLES))
        for row in rows:
            self.assertAlmostEqual(row["jump"], 14*row["h"], delta=2e-11)

    def test_true_log_slope_jump_survives_refinement(self):
        rows = review.centered_slopes("cusp", points(lambda x: 3*x+0.7*abs(x)))
        for row in rows:
            self.assertAlmostEqual(row["left_slope"], 2.3, delta=2e-11)
            self.assertAlmostEqual(row["right_slope"], 3.7, delta=2e-11)
            self.assertAlmostEqual(row["jump"], 1.4, delta=2e-11)

    def test_partial_grid_does_not_invent_endpoints(self):
        data = points(lambda x: x*x)
        del data["sig_p020"]
        rows = review.centered_slopes("partial", data)
        self.assertEqual({r["h"] for r in rows}, {2e-4,1e-4,5e-5})

    def test_local_rate_replacement_recomputes_inverse_relaxation_time(self):
        channels = ("uubar_to_uubar", "uubar_to_ddbar")
        row = dict(point_id="sig_c000",xi=0.015)
        for sp in review.SPECIES:
            row["tauinv_"+sp] = 2.0
            row["tau_"+sp] = 0.5
        weights = [dict(species=sp, **{key+"_coefficient": 1.0 for key in review.FINAL}) for sp in review.SPECIES]
        contributions = [dict(species=sp,channel=c,collision_density_weight_inv_fm3=0.25)
                         for sp in review.SPECIES for c in channels]
        source = dict(observable=row,rates=[dict(channel=c,rate=1.) for c in channels],
                      weights=weights,contributions=contributions)
        bands = [dict(point_id="sig_c000",channel=c,x_right_inv_fm2=1e-3,
                      continuous_full_rate_fm2=0.5,discrete_full_rate_fm2=0.3) for c in channels]
        result, changes = review.hybrid_points({"sig_c000":source},bands,1e-3)
        self.assertEqual(len(changes), 2)
        for sp in review.SPECIES:
            self.assertAlmostEqual(result["sig_c000"]["tau_"+sp], 1/2.1)
            self.assertAlmostEqual(result["sig_c000"]["tauinv_"+sp], 2.1)
        for key in review.FINAL:
            self.assertAlmostEqual(result["sig_c000"][key], 6/2.1)

    def test_precision_control_never_substitutes_a_different_tier_center(self):
        coarse = points(lambda x: 3*x+7*x*x)
        inner = {key: row.copy() for key, row in coarse.items()
                 if key in {"sig_m020", "sig_c000", "sig_p020"}}
        # A common multiplicative precision shift must not change a secant.
        # Inserting only this center into the coarse grid would create 1/h bias.
        for row in inner.values():
            for key in review.OBSERVABLES:
                row[key] *= math.exp(1e-5)
        with patch.object(review, "hybrid_points", side_effect=lambda grid, *_: (grid, [])):
            rows, _ = review.thermal_controls({"s6": coarse, "s7": inner}, [], (0.001,))
        self.assertEqual({r["h"] for r in rows if "_s7_" in r["variant"]}, {2e-5})
        self.assertEqual({r["h"] for r in rows if "_s6_" in r["variant"]}, {2e-4, 1e-4, 5e-5, 2e-5})
        for row in rows:
            self.assertAlmostEqual(row["jump"], 14*row["h"], delta=3e-11)

    def test_combined_counterfactual_rebuilds_tau_and_terminal_values(self):
        ids = ("sig_m020", "sig_c000", "sig_p020")
        channels = sorted(finalizer.CUTOFF_CHANNELS)
        hybrid, sources = {}, {}
        for ident in ids:
            hybrid[ident] = {**{"tauinv_"+sp: 2.1 for sp in review.SPECIES},
                             **{"tau_"+sp: 1/2.1 for sp in review.SPECIES}}
            sources[ident] = dict(rates=[dict(channel=c,rate=1.) for c in channels],
                contributions=[dict(species=sp,channel=c,collision_density_weight_inv_fm3=0.25)
                               for sp in review.SPECIES for c in channels],
                weights=[dict(species=sp,**{key+"_coefficient":1. for key in review.FINAL})
                         for sp in review.SPECIES])
        deltas = {"left":{c:0.2 for c in channels},"right":{c:-0.1 for c in channels}}
        result = finalizer.apply_fixed_membership(review,hybrid,sources,ids,deltas)
        for ident, inverse in ((ids[0],2.2),(ids[2],2.05)):
            for sp in review.SPECIES:
                self.assertAlmostEqual(result[ident]["tauinv_"+sp],inverse)
                self.assertAlmostEqual(result[ident]["tau_"+sp],1/inverse)
            for key in review.FINAL:
                self.assertAlmostEqual(result[ident][key],6/inverse)
        self.assertEqual(result[ids[1]],hybrid[ids[1]])
        self.assertEqual(hybrid[ids[0]]["tauinv_u"],2.1)


if __name__ == "__main__":
    unittest.main()
