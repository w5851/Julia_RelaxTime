import copy
import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[3]
spec = importlib.util.spec_from_file_location("constraint_comparison",
    ROOT/"scripts/analysis/pnjl_bifurcation_spike/compare_constraints.py")
C = importlib.util.module_from_spec(spec)
spec.loader.exec_module(C)


class ConstraintComparisonTests(unittest.TestCase):
    def record(self):
        return dict(config=dict(T_values=[0.5]), rows=[dict(coordinate=0.5, converged=True,
            residual=1e-10, constraint_error=1e-11, solution=[-1,-1,-2,0.1,0.2,1,1,1],
            mu_vec=[1,1,1], masses=[1,1,2], pressure=2.0, rho_norm=0.5, entropy=0.5, energy=1.0)])

    def test_invalid_and_missing_roots_are_visible(self):
        left = self.record()
        right = copy.deepcopy(left)
        right["rows"][0]["converged"] = False
        self.assertFalse(C.compare(left,right)["same"])
        self.assertEqual(C.compare(left,right)["coverage"],[0.5])
        right["rows"] = []
        self.assertFalse(C.compare(left,right)["same"])

    def test_charge_conjugate_degeneracy_does_not_count_as_same_state(self):
        left = self.record()
        right = copy.deepcopy(left)
        right["rows"][0]["mu_vec"] = [-1,-1,-1]
        right["rows"][0]["solution"][5:] = [-1,-1,-1]
        right["rows"][0]["rho_norm"] = -0.5
        check = C.compare(left,right)
        self.assertFalse(check["same"])
        self.assertEqual(check["maxima"]["pressure"],0.0)
        self.assertEqual(check["maxima"]["mu_vec"],2.0)

    def test_duplicate_and_unequal_config_rejected(self):
        left = self.record()
        self.assertTrue(C.compare(left,left)["same"])
        right = copy.deepcopy(left)
        right["rows"].append(right["rows"][0])
        self.assertFalse(C.compare(left,right)["same"])
        right = copy.deepcopy(left)
        right["config"]["p_num"] = 48
        self.assertFalse(C.compare(left,right)["same"])

    def test_hard_constraints_and_vector_dimensions(self):
        for field,value in (("masses",[2,2,1]),("mu_vec",[1,1]),("masses",[None,1,2])):
            row = self.record()["rows"][0]
            row[field] = value
            self.assertFalse(C.admissible(row))


if __name__ == "__main__":
    unittest.main()
