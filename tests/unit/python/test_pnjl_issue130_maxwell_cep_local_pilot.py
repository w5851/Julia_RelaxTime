import csv
from pathlib import Path

from scripts.analysis.collect_pnjl_issue130_maxwell_cep_local_pilot import target_rows


ROOT = Path(__file__).resolve().parents[3]
TARGET_LIST = ROOT / "docs" / "analysis" / "pnjl" / "issue130_endpoint_refinement_preflight_v1" / "maxwell_local" / "tables" / "target_list.csv"
WORKFLOW = ROOT / ".github" / "workflows" / "pnjl-issue130-maxwell-cep-local-pilot.yml"
RUNNER = ROOT / "scripts" / "analysis" / "pnjl_issue130_maxwell_cep_local_pilot.jl"
COLLECTOR = ROOT / "scripts" / "analysis" / "collect_pnjl_issue130_maxwell_cep_local_pilot.py"


def test_preflight_contains_exactly_eleven_authorized_maxwell_targets():
    rows = target_rows(TARGET_LIST)
    assert len(rows) == 11
    assert all(row["target_kind"] == "maxwell_fixed_xi_T" for row in rows.values())
    assert all(row["pilot_selection"] == "pilot_candidate" for row in rows.values())


def test_workflow_matrix_matches_authorized_target_list():
    text = WORKFLOW.read_text(encoding="utf-8")
    with TARGET_LIST.open(newline="", encoding="utf-8") as handle:
        target_ids = [
            row["target_id"]
            for row in csv.DictReader(handle)
            if row.get("pilot_selection") == "pilot_candidate"
        ]
    assert len(target_ids) == 11
    assert all(f"- {target_id}" in text for target_id in target_ids)
    assert "pnjl_issue130_maxwell_cep_local_pilot_v1" in text
    assert "reference_write" in text or "reference_write" in RUNNER.read_text(encoding="utf-8")
    assert "promote_reference" not in text


def test_runner_and_collector_preserve_diagnostic_boundaries():
    runner = RUNNER.read_text(encoding="utf-8")
    collector = COLLECTOR.read_text(encoding="utf-8")
    assert "maxwell_construction" in runner or "strict_candidate" in runner
    assert "reference_write\" => false" in runner
    assert "oracle_labels_consumed\" => false" in runner
    assert "pilot_candidate" in collector
    assert "solver_or_curve_failure" in collector
    assert "MAX_TARGETED = 12" in runner
