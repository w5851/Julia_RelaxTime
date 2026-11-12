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
    def test_metastable_pressure_drift_is_checked(self):
        reference = dict(all_roots=[dict(coordinate=1.0, roots=[point(1.0), point(2.0)])])
        changed = dict(all_roots=[dict(coordinate=1.0, roots=[point(1.0), point(2.0, 1.00001)])])
        result = C.compare_all_roots(reference, changed)
        self.assertFalse(result["same"])
        self.assertGreater(result["maxima"]["pressure"], 1e-6)

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
