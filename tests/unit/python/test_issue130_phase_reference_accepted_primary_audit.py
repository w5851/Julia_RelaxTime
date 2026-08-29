from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path


ROOT = Path(__file__).parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "pnjl" / "audit_issue130_phase_reference_accepted_primary.py"


def _load_script():
    spec = importlib.util.spec_from_file_location("issue130_accepted_primary_audit_test", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_accepted_primary_audit_separates_runtime_and_path_gates(tmp_path: Path) -> None:
    audit = _load_script()
    output = tmp_path / "audit"
    manifest = audit.build_audit(ROOT, output)

    decision = json.loads((output / "decision.json").read_text(encoding="utf-8"))
    assert manifest["schema_version"] == "pnjl_issue130_phase_reference_legacy_audit_v3"
    assert decision["accepted_primary_valid"] is True
    assert decision["strict_explicit_valid"] is True
    assert decision["runtime_legacy_fallback_rows"] == 0
    assert decision["runtime_legacy_rollback_enabled"] is False
    assert decision["accepted_primary_runtime_ready"] is True
    assert decision["path_retirement_ready"] is False
    assert decision["physical_deletion_eligible"] is False
    assert decision["active_path_contract_count"] > 0
    assert decision["verdict"] == "accepted_primary_runtime_ready_path_retirement_pending"


def test_accepted_primary_audit_records_neighborhood_not_exact_substitution(tmp_path: Path) -> None:
    audit = _load_script()
    output = tmp_path / "audit"
    audit.build_audit(ROOT, output)
    rows = (output / "tables" / "coverage.csv").read_text(encoding="utf-8").splitlines()
    assert any(line.startswith("boundary,12537,7162,5375,0,3091,48,14,34,0,34,") for line in rows)
    assert any(line.startswith("crossover,3135,1343,1792,0,1343,336,21,190,125,190,") for line in rows)
    assert any(line.startswith("cep,161,90,71,0,90,11,11,0,0,0,") for line in rows)
    assert any(line.startswith("spinodals,11989,6886,5103,0,6886,57,26,31,0,31,") for line in rows)
    assert (output / "tables" / "neighbor_coverage.csv").is_file()
