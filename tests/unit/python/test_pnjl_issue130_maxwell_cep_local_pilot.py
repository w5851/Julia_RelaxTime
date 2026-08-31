import csv
import json
from pathlib import Path

from scripts.analysis.collect_pnjl_issue130_maxwell_cep_local_pilot import (
    target_rows,
    validate_target,
)


ROOT = Path(__file__).resolve().parents[3]
TARGET_LIST = ROOT / "docs" / "analysis" / "pnjl" / "issue130_endpoint_refinement_preflight_v1" / "maxwell_local" / "tables" / "target_list.csv"
WORKFLOW = ROOT / "docs" / "analysis" / "governance" / "diagnostic_workflow_retirement_wave2_v1" / "definitions" / "pnjl-issue130-maxwell-cep-local-pilot.yml"
EXPANSION_WORKFLOW = ROOT / ".github" / "workflows" / "pnjl-issue130-maxwell-cep-local-expansion.yml"
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


def test_expansion_selection_is_explicit_and_disjoint_from_pilot():
    rows = target_rows(TARGET_LIST, "input_incomplete", expected_count=276)
    assert len(rows) == 276
    assert all(row["pilot_selection"] == "input_incomplete" for row in rows.values())
    text = EXPANSION_WORKFLOW.read_text(encoding="utf-8")
    assert "pnjl_issue130_maxwell_cep_local_expansion_v1" in text
    assert "TARGET_SELECTION: input_incomplete" in text
    assert "EXPECTED_TARGET_COUNT: \"276\"" in text
    assert 'echo "matrix_left={\\"target_id\\":$left_values}"' in text
    assert 'echo "matrix_right={\\"target_id\\":$right_values}"' in text
    assert "matrix: ${{ fromJSON(needs.prepare.outputs.matrix_left) }}" in text
    assert "matrix: ${{ fromJSON(needs.prepare.outputs.matrix_right) }}" in text
    assert "numerical_left:" in text and "numerical_right:" in text
    assert "needs: [prepare, numerical_left, numerical_right]" in text
    assert "matrix: ${{ fromJSON(needs.prepare.outputs.matrix) }}" not in text
    runner = RUNNER.read_text(encoding="utf-8")
    assert "reference_write" in runner
    assert "oracle_labels_consumed" in runner


def test_target_selection_rejects_duplicate_ids(tmp_path):
    path = tmp_path / "targets.csv"
    path.write_text(
        "target_id,pilot_selection,target_kind\n"
        "duplicate,input_incomplete,maxwell_fixed_xi_T\n"
        "duplicate,input_incomplete,maxwell_fixed_xi_T\n",
        encoding="utf-8",
    )
    try:
        target_rows(path, "input_incomplete")
    except ValueError as error:
        assert "duplicate target_id" in str(error)
    else:
        raise AssertionError("duplicate target IDs must be rejected")


def test_runner_and_collector_preserve_diagnostic_boundaries():
    runner = RUNNER.read_text(encoding="utf-8")
    collector = COLLECTOR.read_text(encoding="utf-8")
    assert "maxwell_construction" in runner or "strict_candidate" in runner
    assert "reference_write\" => false" in runner
    assert "oracle_labels_consumed\" => false" in runner
    assert "pilot_candidate" in collector
    assert "solver_or_curve_failure" in collector
    assert "MAX_TARGETED = 12" in runner
    assert "--source-workflow-sha" in collector
    assert "source_workflow_sha" in collector
    assert '"solver_called": args.run_mode == "numerical"' in collector
    assert "args.run_mode == \"numerical\"" in collector
    assert "headSha" in WORKFLOW.read_text(encoding="utf-8")
    assert "SOURCE_WORKFLOW_SHA" in WORKFLOW.read_text(encoding="utf-8")
    assert 'echo "SOURCE_RUN_ID=$SOURCE_RUN_ID" >> "$GITHUB_ENV"' in EXPANSION_WORKFLOW.read_text(encoding="utf-8")


def _write_target_artifact(
    directory: Path,
    *,
    summary_identity: bool = False,
    null_summary_identity: bool = False,
    bad_provenance: bool = False,
):
    directory.mkdir()
    calculation_sha = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
    workflow_sha = "51c93cf0111415b35bb199376c358782c0f5a2f4"
    summary = {
        "schema_version": "pnjl_issue130_maxwell_cep_local_pilot_v1",
        "target_id": "synthetic_target",
        "xi": 0.0,
        "T_MeV": 123.0,
        "verdict": "candidate_and_geometry_feasible",
        "final_status": "first_order",
        "final_reason": "unique_three_crossing_candidate",
        "final_candidate_count": 1,
        "final_candidate_mu_MeV": 304.0,
        "final_area_residual": 1e-6,
        "final_geometry_converged": True,
        "targeted_additions": 1,
        "curve_points": 2,
        "finite_and_converged": True,
    }
    if summary_identity:
        summary.update(calculation_sha=calculation_sha, workflow_head_sha=workflow_sha)
    if null_summary_identity:
        summary.update(calculation_sha=None, workflow_head_sha=None)
    provenance = {
        "schema_version": "pnjl_issue130_maxwell_cep_local_pilot_v1",
        "target_id": "synthetic_target",
        "calculation_sha": calculation_sha,
        "workflow_head_sha": "bad" if bad_provenance else workflow_sha,
        "solver_called": True,
        "reference_write": False,
        "oracle_labels_consumed": False,
    }
    manifest = {
        "schema_version": "pnjl_issue130_maxwell_cep_local_pilot_v1",
        "target_id": "synthetic_target",
        "calculation_sha": calculation_sha,
        "workflow_head_sha": workflow_sha,
        "reference_write": False,
        "oracle_labels_consumed": False,
    }
    (directory / "target_summary.json").write_text(json.dumps(summary), encoding="utf-8")
    (directory / "provenance.json").write_text(json.dumps(provenance), encoding="utf-8")
    (directory / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    (directory / "curve_points.csv").write_text(
        "target_id,xi,T_MeV,rho,muq_MeV,residual_norm,converged,finite\n"
        "synthetic_target,0.0,123.0,0.0,300.0,1e-8,true,true\n"
        "synthetic_target,0.0,123.0,1.0,304.0,1e-8,true,true\n",
        encoding="utf-8",
    )
    (directory / "slice_metrics.csv").write_text(
        "target_id,targeted_additions\nsynthetic_target,1\n", encoding="utf-8"
    )
    (directory / "policy_frontier.csv").write_text(
        "target_id\nsynthetic_target\n", encoding="utf-8"
    )
    (directory / "method_costs.csv").write_text(
        "target_id\nsynthetic_target\n", encoding="utf-8"
    )


def test_missing_summary_identity_uses_verified_sources(tmp_path):
    directory = tmp_path / "target"
    _write_target_artifact(directory)
    summary, _, errors, _, fallbacks = validate_target(
        "synthetic_target",
        {"xi": "0.0", "T_MeV": "123.0", "reason": "synthetic"},
        directory,
        "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48",
        "51c93cf0111415b35bb199376c358782c0f5a2f4",
    )
    assert errors == []
    assert fallbacks == ["calculation_sha", "workflow_head_sha"]
    assert summary["identity_fallback_fields"] == "calculation_sha;workflow_head_sha"


def test_null_summary_identity_uses_verified_sources(tmp_path):
    directory = tmp_path / "target"
    _write_target_artifact(directory, null_summary_identity=True)
    _, _, errors, _, fallbacks = validate_target(
        "synthetic_target",
        {"xi": "0.0", "T_MeV": "123.0", "reason": "synthetic"},
        directory,
        "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48",
        "51c93cf0111415b35bb199376c358782c0f5a2f4",
    )
    assert errors == []
    assert fallbacks == ["calculation_sha", "workflow_head_sha"]


def test_present_summary_identity_is_checked_strictly(tmp_path):
    directory = tmp_path / "target"
    _write_target_artifact(directory, summary_identity=True)
    summary_path = directory / "target_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["calculation_sha"] = "0" * 40
    summary_path.write_text(json.dumps(summary), encoding="utf-8")
    _, _, errors, _, fallbacks = validate_target(
        "synthetic_target",
        {"xi": "0.0", "T_MeV": "123.0", "reason": "synthetic"},
        directory,
        "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48",
        "51c93cf0111415b35bb199376c358782c0f5a2f4",
    )
    assert any("summary calculation_sha mismatch" in error for error in errors)
    assert fallbacks == []


def test_missing_summary_identity_does_not_mask_source_mismatch(tmp_path):
    directory = tmp_path / "target"
    _write_target_artifact(directory, bad_provenance=True)
    _, _, errors, _, fallbacks = validate_target(
        "synthetic_target",
        {"xi": "0.0", "T_MeV": "123.0", "reason": "synthetic"},
        directory,
        "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48",
        "51c93cf0111415b35bb199376c358782c0f5a2f4",
    )
    assert any("provenance workflow_head_sha mismatch" in error for error in errors)
    assert any("summary workflow_head_sha missing" in error for error in errors)
    assert fallbacks == ["calculation_sha"]
