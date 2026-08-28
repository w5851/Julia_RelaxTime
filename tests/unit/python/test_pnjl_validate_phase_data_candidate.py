from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "pnjl" / "validate_phase_data.py"
SPEC = importlib.util.spec_from_file_location("pnjl_validate_phase_data_candidate", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_candidate_schema_validation_is_explicit_and_solver_free() -> None:
    root = ROOT / "data" / "reference" / "pnjl" / "issue130_phase_reference_v1"
    issues, statistics = MODULE.validate_candidate_reference(root, "strict")
    assert statistics["boundary"]["rows"] == 7162
    assert statistics["crossover"]["rows"] == 1343
    assert statistics["cep"]["rows"] == 93
    assert statistics["spinodals"]["rows"] == 6886
    assert sum(row["unresolved"] for row in statistics.values()) == 4074
    assert not any(issue.severity == "error" for issue in issues)
