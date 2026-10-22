"""Physics-coordinate and fixed-observation checks for the exploratory band."""
import copy
import math
from pathlib import Path
import sys
import tempfile
import tomllib
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from scripts.analysis.relaxtime.freezeout_parameter_envelope import (
    DEFAULT_CATALOG, source_parameters, pooled_bounds, interior_bounds,
    trajectory, trajectory_valid, in_box, initial_candidates,
    assess_models, complete_curve_band, select_refined_curves, remap_peak, verify_sources,
)
from scripts.plotting.plot_manifest import sha256_file


class ParameterEnvelopeTests(unittest.TestCase):
    def setUp(self):
        self.catalog = tomllib.loads(DEFAULT_CATALOG.read_text(encoding="utf-8"))

    def test_source_directory_is_repository_relative_and_rejects_changed_bytes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source_dir = root / "inputs"
            source_dir.mkdir()
            source = source_dir / "source.pdf"
            source.write_bytes(b"retained source fixture")
            catalog = dict(source_root="inputs", sources=[dict(pdf=source.name, sha256=sha256_file(source))])
            with patch("scripts.analysis.relaxtime.freezeout_parameter_envelope.ROOT", root):
                relative = verify_sources(catalog)
                catalog["source_root"] = str(source_dir)
                self.assertEqual(relative, verify_sources(catalog))
                source.write_bytes(b"changed source fixture")
                with self.assertRaisesRegex(ValueError, "source PDF hash mismatch"):
                    verify_sources(catalog)

    def test_normalized_quartic_is_the_same_temperature(self):
        for source in self.catalog["sources"]:
            p = source_parameters(source)
            if source["kind"] != "normalized_quartic":
                continue
            t0, k2, k4, d, e = source["values"]
            for energy in (4.87, 7.7, 19.6, 200.):
                actual_t, actual_mu = trajectory(p, energy)
                mu = d / (1 + e*energy)
                self.assertAlmostEqual(actual_mu, mu, places=14)
                self.assertAlmostEqual(actual_t, t0*(1-k2*(mu/t0)**2-k4*(mu/t0)**4), places=14)

    def test_pooled_bounds_retain_marginals_and_partial_freezeout_constraint(self):
        bounds = pooled_bounds(self.catalog)
        self.assertAlmostEqual(bounds[0][0], .143)
        self.assertAlmostEqual(bounds[0][1], .168)
        self.assertEqual(bounds[1][0], 0.)
        self.assertAlmostEqual(bounds[1][1], .155)
        self.assertAlmostEqual(bounds[2][0], .032)
        self.assertAlmostEqual(bounds[2][1], .004/.1539**3)
        self.assertAlmostEqual(bounds[3][1], 1.336)
        self.assertAlmostEqual(bounds[4][1], .281)
        # 144+-6 MeV in Borsanyi is a finite-energy temperature range, not Tf0.
        self.assertGreater(bounds[0][0], .138)

    def test_coordinate_mixing_is_permitted_without_claiming_a_source_fit(self):
        candidates, bounds, inner, rejected = initial_candidates(self.catalog)
        repeated = initial_candidates(self.catalog)
        self.assertEqual((candidates, bounds, inner, rejected), repeated)
        self.assertEqual(sum(c["role"] == "interior_sample" for c in candidates), 32)
        self.assertEqual(sum(c["role"] == "source_anchor" for c in candidates), 8)
        self.assertEqual(sum(c["role"] == "explicit_hybrid" for c in candidates), 6)
        self.assertEqual(len({c["name"] for c in candidates}), len(candidates))
        for c in candidates:
            self.assertTrue(in_box(c["values"], bounds))
            self.assertTrue(trajectory_valid(c["values"], [4.87, 200.]))
            if c["role"] == "interior_sample":
                self.assertTrue(in_box(c["values"], inner))
            ts = [trajectory(c["values"], e)[0] for e in (4.87, 7.7, 27., 200.)]
            self.assertEqual(ts, sorted(ts))
        self.assertAlmostEqual(inner[0][0], .1455)
        with self.assertRaises(ValueError):
            interior_bounds(bounds, .5)

    def test_invalid_combination_can_lie_inside_coordinate_box(self):
        p = [.145, .15, 1.09, 1.33, .147]
        self.assertTrue(in_box(p, pooled_bounds(self.catalog)))
        self.assertFalse(trajectory_valid(p, [4.87, 200.]))
        self.assertFalse(trajectory_valid([math.nan, .1, .1, 1., .2], [4.87, 200.]))

    def test_energy_mapping_preserves_the_same_thermodynamic_peak(self):
        p = [.147,.1,.55,.84,.223]
        mapped = remap_peak(p,11.5,7.7,.273)
        self.assertEqual(mapped[:3],p[:3])
        t0,mu0=trajectory(p,11.5)
        t1,mu1=trajectory(mapped,7.7)
        self.assertAlmostEqual(mu0,mu1,places=14)
        self.assertAlmostEqual(t0,t1,places=14)
        self.assertTrue(in_box(mapped,pooled_bounds(self.catalog)))
        with self.assertRaises(ValueError):
            remap_peak(p,11.5,0,.273)


class FixedObservationBandTests(unittest.TestCase):
    def fixture(self):
        experiments = [dict(energy=7.7, ratio=.22, primary_auau=True),
            dict(energy=200., ratio=.16, primary_auau=True)]
        models = {n: {r["energy"]: dict(ratio=r["ratio"]+offset, passed=True) for r in experiments}
            for n, offset in (("a", .02), ("b", .04))}
        candidates = {n: dict(role="interior_sample") for n in models}
        return models, experiments, candidates

    def test_failure_cannot_improve_score_by_deleting_an_observation(self):
        models, data, candidates = self.fixture()
        models["a"][200.] = dict(ratio=math.nan, passed=False)
        result = assess_models(models, data, candidates, [7.7, 200.], 1.1)
        self.assertEqual(result["ranking"], ["b"])
        self.assertIsNone(result["profiles"]["a"]["rmse"])
        self.assertEqual(result["profiles"]["a"]["failed_required_energies"], [200.])
        self.assertEqual(result["profiles"]["b"]["count"], 2)
        self.assertAlmostEqual(result["profiles"]["b"]["rmse"], .04)

    def test_theory_grid_density_does_not_weight_the_fit(self):
        models, data, candidates = self.fixture()
        result = assess_models(models, data, candidates, [7.7, 200.], 1.1)
        for points in models.values():
            points[10.] = dict(ratio=9., passed=True)
        updated = assess_models(models, data, candidates, [7.7, 200.], 1.1)
        self.assertEqual(result["ranking"], updated["ranking"])
        for name in models:
            for key in ("rmse", "mae", "mean_signed_bias", "count", "eight_point_rmse", "high_energy"):
                self.assertEqual(result["profiles"][name][key], updated["profiles"][name][key])

    def test_horn_diagnostic_reports_shape_separately_from_primary_score(self):
        models, data, candidates = self.fixture()
        data.extend([dict(energy=4.87, ratio=.20, primary_auau=True),
            dict(energy=11.5, ratio=.20, primary_auau=True)])
        models["a"].update({4.87:dict(ratio=.19,passed=True), 11.5:dict(ratio=.21,passed=True)})
        result = assess_models(models, data, candidates, [7.7, 200.], 1.1)
        a = result["profiles"]["a"]
        self.assertAlmostEqual(a["rmse"], .02)
        self.assertAlmostEqual(a["horn_three_point_difference_rmse"], math.sqrt(.001/2))
        self.assertTrue(a["horn_rises_then_falls"])
        self.assertEqual(a["sampled_horn"]["energy_GeV"], 7.7)
        self.assertFalse(a["sampled_horn"]["peak_at_grid_boundary"])

    def test_band_uses_fixed_complete_curves_and_records_extremizers(self):
        models, _, _ = self.fixture()
        models["a"][200.]["ratio"] = .23  # boundaries can cross
        rows, complete = complete_curve_band(models, ["a", "b"], [7.7, 200.])
        self.assertEqual(complete, ["a", "b"])
        self.assertEqual(rows[0]["minimum_profile"], "a")
        self.assertEqual(rows[1]["minimum_profile"], "b")
        models["a"][200.].update(passed=False, ratio=math.nan)
        rows, complete = complete_curve_band(models, ["a", "b"], [7.7, 200.])
        self.assertEqual(complete, ["b"])
        self.assertEqual(rows[0]["minimum"], rows[0]["maximum"])
        self.assertAlmostEqual(rows[0]["minimum"], .26)

    def test_refined_selection_retains_the_rmse_leader_and_horn_alternative(self):
        result = dict(ranking=["a","b","c","d"], near_best_profiles=["a","b","c"],
            profiles={"a":dict(horn_rises_then_falls=False),
                "b":dict(horn_rises_then_falls=True,horn_three_point_difference_rmse=.04),
                "c":dict(horn_rises_then_falls=True,horn_three_point_difference_rmse=.01),
                "d":dict(horn_rises_then_falls=True,horn_three_point_difference_rmse=0.)})
        self.assertEqual(select_refined_curves(result), ["a","c","b"])
        self.assertEqual(select_refined_curves(result,1), ["a"])
        with self.assertRaises(ValueError):
            select_refined_curves(result,0)


if __name__ == "__main__":
    unittest.main()
