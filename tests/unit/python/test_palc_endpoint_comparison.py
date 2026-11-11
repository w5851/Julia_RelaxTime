import importlib.util
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = ROOT / "scripts/analysis/pnjl_bifurcation_spike"
sys.path.insert(0, str(SCRIPT_DIR))
try:
    spec = importlib.util.spec_from_file_location("endpoint_comparison", SCRIPT_DIR / "compare_endpoints.py")
    C = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(C)
finally:
    sys.path.pop(0)


def point(rho):
    return dict(coordinate=1.0, converged=True, residual=1e-10, constraint_error=0.0,
                solution=[-rho,-rho,-2.0,0.1,0.2], mu_vec=[1.0]*3, masses=[1.0,1.0,2.0],
                pressure=2.0 if rho==1.0 else 1.0, rho_norm=rho, entropy=0.5, energy=1.0)


class EndpointCoverageTests(unittest.TestCase):
    def test_selected_agreement_does_not_hide_missing_middle_root(self):
        reference = dict(config={"T_fm": 0.6}, rows=[point(1.0)],
                         all_roots=[dict(coordinate=1.0,roots=[point(1.0),point(2.0),point(3.0)])])
        candidate = dict(config=reference["config"],rows=reference["rows"],
                         all_roots=[dict(coordinate=1.0,roots=[point(1.0),point(3.0)])])
        self.assertTrue(C.compare(reference,candidate)["same"])
        coverage = C.root_coverage(reference,candidate)
        self.assertFalse(coverage["same"])
        self.assertEqual(coverage["missing"],[dict(coordinate=1.0,rho_norm=2.0)])

    def test_candidate_order_does_not_change_coverage(self):
        a = dict(all_roots=[dict(coordinate=1.0,roots=[point(1.0),point(2.0)])])
        b = dict(all_roots=[dict(coordinate=1.0,roots=[point(2.0),point(1.0)])])
        self.assertTrue(C.root_coverage(a,b)["same"])

    def test_invalid_root_and_missing_coordinate_are_visible(self):
        root = point(1.0)
        root["converged"] = False
        a = dict(all_roots=[dict(coordinate=1.0,roots=[point(1.0)])])
        b = dict(all_roots=[dict(coordinate=1.0,roots=[root])])
        self.assertFalse(C.root_coverage(a,b)["same"])
        self.assertEqual(len(C.root_coverage(a,dict(all_roots=[]))["missing"]),1)


if __name__ == "__main__":
    unittest.main()
