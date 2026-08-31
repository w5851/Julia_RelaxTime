from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

from scripts.analysis.pnjl.audit_issue130_phase_reference_legacy_retirement_v2 import build_audit


ROOT = Path(__file__).resolve().parents[3]
RECOVERY_REF = "9aa4c313901ca0c91e851f58514e3df9aa124df4"


def _recovery_snapshot(tmp_path: Path) -> Path:
    root = tmp_path / "legacy_snapshot"
    for name in (
        "README.md",
        "RETIREMENT_MANIFEST.json",
        "boundary.csv",
        "cep.csv",
        "crossover_dense.csv",
        "crossover_dense.meta.json",
        "phase_reference_dense_manifest.json",
        "spinodals.csv",
    ):
        payload = subprocess.check_output(
            ["git", "show", f"{RECOVERY_REF}:data/reference/pnjl/legacy_phase_reference_v1/{name}"],
            cwd=ROOT,
        )
        root.mkdir(parents=True, exist_ok=True)
        (root / name).write_bytes(payload)
    return root


def test_post_acceptance_audit_keeps_runtime_and_downstream_boundaries(tmp_path: Path) -> None:
    manifest = build_audit(ROOT, tmp_path / "legacy_audit_v2", legacy_root=_recovery_snapshot(tmp_path))
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
    # The retired RS pilot definition used to contribute one active blocker.
    # Other historical legacy locators may still grow, so only the semantic
    # lower bound changes from 25 to 24.
    assert decision["active_consumer_blocker_count"] >= 24
    assert decision["unknown_active_reference_count"] == 0
    assert decision["physical_deletion_eligible"] is False

    with (output / "tables" / "consumer_matrix.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        consumers = list(csv.DictReader(handle))
    retired_path = (
        "docs/analysis/governance/diagnostic_workflow_retirement_wave1_v1/"
        "definitions/relaxtime-issue130-rs-numerical-pilot-v1.yml"
    )
    retired = next(row for row in consumers if row["path"] == retired_path)
    assert retired["active"] == "False"
    assert retired["retirement_blocker"] == "False"
    assert retired["decision"] == "not_an_active_runtime_consumer"
    assert not any(
        row["path"] == ".github/workflows/relaxtime-issue130-rs-numerical-pilot-v1.yml"
        for row in consumers
    )

    claims = json.loads((output / "tables" / "claim_ledger.json").read_text(encoding="utf-8"))
    accepted_claim = next(item for item in claims if item["claim_id"] == "accepted_downstream_default")
    assert accepted_claim["status"] == "supported"
    assert "not runtime input" in accepted_claim["boundary"]
