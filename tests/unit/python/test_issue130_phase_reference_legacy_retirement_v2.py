from __future__ import annotations

import json
from pathlib import Path

from scripts.analysis.pnjl.audit_issue130_phase_reference_legacy_retirement_v2 import build_audit


ROOT = Path(__file__).resolve().parents[3]


def test_post_acceptance_audit_keeps_runtime_and_downstream_boundaries(tmp_path: Path) -> None:
    manifest = build_audit(ROOT, tmp_path / "legacy_audit_v2")
    assert manifest["schema_version"] == "pnjl_issue130_phase_reference_legacy_retirement_audit_v2"
    assert manifest["downstream_default_layer"] == "accepted"
    assert manifest["runtime_reference_layer"] == "strict"
    assert manifest["accepted_promotion_status"] == "accepted_for_downstream"
    assert manifest["accepted_runtime_consumption"] is False
    assert manifest["accepted_reference_write"] is False
    assert manifest["solver_called"] is False
    assert manifest["reference_write"] is False
    assert manifest["runtime_consumption"] is False

    output = tmp_path / "legacy_audit_v2"
    decision = json.loads((output / "decision.json").read_text(encoding="utf-8"))
    assert decision["schema_version"] == "pnjl_issue130_phase_reference_legacy_retirement_decision_v2"
    assert decision["accepted_downstream_default"] is True
    assert decision["legacy_fallback_key_count"] > 0
    assert decision["active_consumer_blocker_count"] == 25
    assert decision["unknown_active_reference_count"] == 0
    assert decision["physical_deletion_eligible"] is False

    claims = json.loads((output / "tables" / "claim_ledger.json").read_text(encoding="utf-8"))
    accepted_claim = next(item for item in claims if item["claim_id"] == "accepted_downstream_default")
    assert accepted_claim["status"] == "supported"
    assert "not runtime input" in accepted_claim["boundary"]
