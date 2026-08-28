from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "pnjl" / "audit_issue130_phase_reference_candidate_only_consumers.py"
SPEC = importlib.util.spec_from_file_location("issue130_candidate_only_consumers", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_real_candidate_only_consumer_audit_is_solver_free(tmp_path: Path) -> None:
    manifest = MODULE.build_audit(ROOT, tmp_path / "audit")
    assert manifest["verdict"] == "candidate_only_contract_supported"
    assert manifest["solver_called"] is False
    assert manifest["reference_write"] is False
    assert manifest["runtime_consumption"] is False
    decision = MODULE.read_json(tmp_path / "audit" / "decision.json")
    assert decision["legacy_rollback_contract_supported"] is True
    assert decision["dynamic_request_key_coverage_complete"] is False
    assert decision["physical_deletion_eligible"] is False
    matrix = (tmp_path / "audit" / "dynamic_request_matrix.csv").read_text(encoding="utf-8")
    assert "cep,0.5,False,0,True,fallback_required" in matrix
