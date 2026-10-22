"""Scientific data selection and comparison tests, independent of the solver."""
import copy
import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from scripts.analysis.relaxtime.prepare_freezeout_experiment_data import prepare, errors_quadrature, DEFAULT_INPUT
from scripts.analysis.relaxtime.plot_freezeout_comparison import compare, compare_resolutions


class FreezeoutExperimentTests(unittest.TestCase):
    def test_published_charge_and_centrality_selection(self):
        rows = {r["id"]: r for r in prepare(DEFAULT_INPUT)}
        self.assertEqual(len(rows), 22)
        self.assertEqual(sum(r["primary_auau"] for r in rows.values()), 14)
        self.assertEqual(rows["STAR_BES_7.7"]["ratio"], .2227)
        self.assertNotEqual(rows["STAR_BES_7.7"]["ratio"], .0770)  # negative charge
        self.assertEqual(rows["STAR_200"]["ratio"], .159)
        self.assertEqual(rows["STAR_130"]["centrality"], "0-6%")
        self.assertEqual(rows["NA49_7.62"]["ratio"], .253)
        self.assertEqual(rows["NA49_7.62"]["system"], "Pb+Pb")
        self.assertEqual(rows["AGS_3.32"]["system"], "Au+Au")
        self.assertFalse(rows["ALICE_5020"]["primary_auau"])
        self.assertEqual(rows["ALICE_5020"]["ratio"], .16085)
        self.assertFalse(rows["HADES_2.4"]["primary_auau"])
        self.assertAlmostEqual(rows["STAR_14.5"]["ratio"], 26.7/141)
        self.assertTrue(all("pp" not in r["system"] for r in rows.values()))

    def test_asymmetric_uncertainty_and_subcomponent_exclusion(self):
        lo, hi = errors_quadrature([{"symerror": ".01"},
            {"asymerror": {"minus": "-.021", "plus": ".016"}}])
        self.assertAlmostEqual(lo, math.hypot(.01,.021))
        self.assertAlmostEqual(hi, math.hypot(.01,.016))
        self.assertEqual(errors_quadrature([{"label":"syst.","symerror":.013},
            {"label":"syst. uncorr.","symerror":.006}], exclude=("syst. uncorr.",)),(.013,.013))


class ComparisonTests(unittest.TestCase):
    def fixture(self):
        data = [dict(id=identifier,energy=e,ratio=r,primary_auau=au,auxiliary_horn=not au,
            system="Au+Au" if au else "Pb+Pb") for identifier,e,r,au in (
            ("STAR_BES_7.7",7.7,.2227,True),("STAR_200",200.,.159,True),
            ("NA49_7.62",7.62,.253,False),("NA49_17.27",17.27,.179,False))]
        models = {name:{r["energy"]:dict(energy=r["energy"],ratio=r["ratio"]+offset,
            passed=True,T_MeV=150.,muB_MeV=200.) for r in data} for name,offset in (("a",.1),("b",.02))}
        return models,data

    def test_equal_observation_weight_not_grid_density(self):
        models,data = self.fixture()
        result,_ = compare(models,data)
        self.assertEqual(result["best_default_on_primary_metric"],"b")
        self.assertAlmostEqual(result["profiles"]["a"]["AuAu_primary"]["rmse"],.1)
        for name in models:
            for energy in [5+i*.01 for i in range(100)]:
                models[name][energy] = dict(energy=energy,ratio=.3,passed=True,T_MeV=150.,muB_MeV=200.)
        dense,_ = compare(models,data)
        self.assertEqual(result["profiles"]["a"]["AuAu_primary"],dense["profiles"]["a"]["AuAu_primary"])

    def test_failures_use_common_observations_without_zero_fill(self):
        models,data = self.fixture()
        models["a"][7.7].update(passed=False,ratio=math.nan)
        result,residuals = compare(models,data)
        self.assertEqual(result["common_primary_points"],1)
        self.assertEqual(result["excluded_failed_observations"][0]["id"],"STAR_BES_7.7")
        self.assertTrue(all(r["sqrt_s_NN_GeV"]!=7.7 for r in residuals))
        self.assertEqual(result["profiles"]["b"]["AuAu_primary"]["count"],1)

    def test_missing_energy_requires_direct_calculation(self):
        models,data = self.fixture()
        del models["b"][7.7]
        with self.assertRaisesRegex(ValueError,"directly; no interpolation"):
            compare(models,data)

    def test_resolution_differences_retain_failure_changes(self):
        reference, _ = self.fixture()
        refined = copy.deepcopy(reference)
        refined["a"][7.7]["ratio"] += .01
        refined["a"][200.].update(passed=False, ratio=math.nan)
        refined["a"][9.5] = dict(energy=9.5,ratio=.4,passed=True,T_MeV=150.,muB_MeV=200.)
        result = compare_resolutions(refined,reference)["profiles"]
        self.assertEqual(result["a"]["common_valid_points"],3)
        self.assertEqual(result["a"]["added_refined_energies_GeV"],[9.5])
        self.assertAlmostEqual(result["a"]["max_absolute_difference"]["absolute_ratio_difference"],.01)
        self.assertEqual(result["a"]["failure_comparison"],
            [dict(energy_GeV=200.,reference_passed=True,refined_passed=False)])
        self.assertEqual(result["b"]["max_absolute_difference"]["absolute_ratio_difference"],0.)

    def test_resolution_comparison_rejects_physics_or_support_change(self):
        reference, _ = self.fixture()
        refined = copy.deepcopy(reference)
        refined["a"][7.7]["T_MeV"] += 1.
        with self.assertRaisesRegex(ValueError,"thermodynamic coordinates changed"):
            compare_resolutions(refined,reference)
        refined = copy.deepcopy(reference)
        del refined["a"][7.7]
        with self.assertRaisesRegex(ValueError,"missing reference energies"):
            compare_resolutions(refined,reference)


if __name__ == "__main__":
    unittest.main()
