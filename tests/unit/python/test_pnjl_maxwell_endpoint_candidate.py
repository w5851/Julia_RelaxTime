import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
REPLAY = ROOT / "scripts" / "analysis" / "pnjl_maxwell_endpoint_candidate_feasibility.jl"
ENDPOINT = ROOT / "scripts" / "analysis" / "pnjl_maxwell_endpoint_refinement.jl"
COLLECTOR = ROOT / "scripts" / "analysis" / "collect_pnjl_maxwell_endpoint_refinement.py"
WORKFLOW = ROOT / "docs" / "analysis" / "governance" / "diagnostic_workflow_retirement_wave2_v1" / "definitions" / "pnjl-maxwell-endpoint-refinement.yml"
EVIDENCE = ROOT / "docs" / "analysis" / "pnjl" / "cep_maxwell" / "maxwell_contracts" / "pnjl_maxwell_endpoint_candidate_feasibility_v1"


def test_solver_free_replay_contract_isolated_from_equilibrium_solver():
    text = REPLAY.read_text(encoding="utf-8")
    assert "maxwell_construction" in text
    assert "crossing_count_not_three" in text
    assert "unique_three_crossing_candidate" in text
    assert "run_production_phase_pipeline" not in text
    assert "solver_called" in text
    assert "30809754119" in text
    assert "30805637032" in text
    assert "3217bed3635574f00c04cbee75e843b4c49451db" in text


def test_endpoint_runner_is_actions_only_and_positive_bracket_bound():
    text = ENDPOINT.read_text(encoding="utf-8")
    assert "new_rho_point_session" in text
    assert "positive_rho_bracket" in text
    assert "TARGETED_CAPS = (4, 6, 8, 10, 12)" in text
    assert "0.0001953125" in text
    assert "reference_write" in text


def test_endpoint_workflow_is_single_anchor_and_sha_bound():
    text = WORKFLOW.read_text(encoding="utf-8")
    assert "calculation_ref" in text
    assert "workflow-head-sha" in text
    assert "pnjl_maxwell_endpoint_refinement.jl" in text
    assert "xi=-0.5 T=5" in text
    assert "promote_reference" not in text


def test_candidate_evidence_is_solver_free_and_gate_bound():
    manifest = json.loads((EVIDENCE / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["schema_version"] == "pnjl_maxwell_endpoint_candidate_feasibility_v1"
    assert manifest["verdict"] == "feasible_candidate"
    assert manifest["solver_called"] is False
    assert manifest["required_anchor_gate"] is True
    assert manifest["unique_candidate_gate"] is True
    assert manifest["promotion_boundary"].startswith("diagnostic_only")
    summary = list(csv.DictReader((EVIDENCE / "candidate_summary.csv").open(encoding="utf-8")))
    required = {
        ("independent_oracle", "-0.5", "5.0"),
        ("independent_oracle", "-0.5", "20.0"),
        ("independent_oracle", "0.0", "5.0"),
    }
    rows = {(row["method"], row["xi"], row["T_MeV"]): row for row in summary}
    assert required <= rows.keys()
    assert all(rows[key]["status"] == "first_order" for key in required)
    assert (EVIDENCE / "figures" / "plot_manifest.json").is_file()
