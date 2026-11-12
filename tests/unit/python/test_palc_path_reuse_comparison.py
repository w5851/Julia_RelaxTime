import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts/analysis/pnjl_bifurcation_spike"))
try:
    import compare_path_reuse as C
finally:
    sys.path.pop(0)


def point(rho, pressure=1.0):
    return dict(coordinate=1.0, converged=True, residual=1e-10, constraint_error=0.0,
        solution=[-rho, -rho, -2.0, 0.1, 0.2], mu_vec=[1.0]*3, masses=[1.0, 1.0, 2.0],
        pressure=pressure, rho_norm=rho, entropy=0.5, energy=1.0)


class PathReuseComparisonTests(unittest.TestCase):
    def test_historical_entropy_tie_requires_both_conjugate_roots_and_roundoff_pressure(self):
        positive = point(1.0,pressure=20.0)
        positive["solution"] = [-1.0,-1.0,-2.0,0.1,0.2,1.0,1.0,1.0]
        negative = dict(positive,solution=[-1.0,-1.0,-2.0,0.2,0.1,-1.0,-1.0,-1.0],
                        mu_vec=[-1.0]*3,rho_norm=-1.0)
        left = dict(config=dict(constraint_mode="fixed_entropy"),rows=[positive],
                    all_roots=[dict(coordinate=1.0,roots=[positive,negative])])
        right = dict(left, rows=[negative])
        def ties(a,b):
            return C.entropy_conjugate_ties(a,b,C.compare(a,b),C.compare_all_roots(a,b))
        self.assertIsNotNone(ties(left,right))
        self.assertFalse(C.compare(left,right)["same"])
        shifted = dict(negative,pressure=20.0+1e-8)
        self.assertIsNone(ties(left,dict(right,rows=[shifted],
            all_roots=[dict(coordinate=1.0,roots=[positive,shifted])])) )
        self.assertIsNone(ties(left,dict(right,all_roots=[dict(coordinate=1.0,roots=[negative])])))
        self.assertIsNone(ties(dict(left,config=dict(constraint_mode="fixed_sigma")),
                              dict(right,config=dict(constraint_mode="fixed_sigma"))))

    def test_completion_rejects_missing_case_and_duplicate_run(self):
        plan = [dict(name="a", methods=["baseline","default"], repetitions=2)]
        records = [dict(name="a",method=m,repetition=r) for m in plan[0]["methods"] for r in (1,2)]
        self.assertTrue(C.check_completion(dict(expected_runs=plan,records=records))["complete"])
        self.assertFalse(C.check_completion(dict(expected_runs=plan,records=records[:-1]))["complete"])
        self.assertFalse(C.check_completion(dict(expected_runs=plan,records=records+records[:1]))["complete"])
        self.assertFalse(C.check_completion(dict(records=records))["complete"])

    def test_invalid_thermodynamics_fail_without_crashing(self):
        reference = dict(all_roots=[dict(coordinate=1.0, roots=[point(1.0)])])
        invalid = point(1.0)
        invalid["pressure"] = None
        candidate = dict(all_roots=[dict(coordinate=1.0, roots=[invalid])])
        self.assertFalse(C.compare_all_roots(reference, candidate)["same"])

    def test_metastable_pressure_drift_is_checked(self):
        reference = dict(all_roots=[dict(coordinate=1.0, roots=[point(1.0), point(2.0)])])
        changed = dict(all_roots=[dict(coordinate=1.0, roots=[point(1.0), point(2.0, 1.00001)])])
        result = C.compare_all_roots(reference, changed)
        self.assertFalse(result["same"])
        self.assertGreater(result["maxima"]["pressure"], 1e-6)

    def test_duplicate_root_does_not_pass_coverage(self):
        reference = dict(all_roots=[dict(coordinate=1.0, roots=[point(1.0)])])
        candidate = dict(all_roots=[dict(coordinate=1.0, roots=[point(1.0),point(1.0)])])
        self.assertFalse(C.compare_all_roots(reference,candidate)["same"])

    def test_missing_middle_root_fails_summary(self):
        rows = [point(1.0)]
        baseline = dict(name="fold", method="baseline", repetition=2, seconds=0.1,
            allocated_bytes=1024, gc_seconds=0.0, config={}, rows=rows,
            all_roots=[dict(coordinate=1.0, roots=rows + [point(2.0)])], work={}, diagnostics={})
        changed = dict(baseline, method="combined", all_roots=[dict(coordinate=1.0, roots=rows)])
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory)/"input.json"
            destination = Path(directory)/"result.json"
            source.write_text(json.dumps(dict(records=[baseline, changed])), encoding="utf-8")
            self.assertFalse(C.summarize(source, destination))
            self.assertFalse(json.loads(destination.read_text())["all_checks_passed"])


if __name__ == "__main__":
    unittest.main()
