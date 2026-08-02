import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "pnjl_maxwell_tolerance_feasibility.jl"
EVIDENCE = ROOT / "docs" / "analysis" / "pnjl_maxwell_tolerance_contract_feasibility_v1"


def test_feasibility_runner_is_solver_free_and_contract_bound():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "maxwell_construction" in text
    assert "_compare_phase_geometry" in text
    assert "run_production_phase_pipeline" not in text
    assert "solve_pnjl_point" not in text
    assert "solver_called" in text
    assert "30676440627" in text
    assert "fde2b929b60575f1daacb84a1b9b8ff6e3b0a0cc" in text


def test_feasibility_evidence_schema_and_gate_boundary():
    manifest = json.loads((EVIDENCE / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["schema_version"] == "pnjl_maxwell_tolerance_contract_feasibility_v1"
    assert manifest["solver_called"] is False
    assert manifest["strict_solver_tolerance"] == 5e-6
    assert manifest["outer_geometry_gate"]["maxwell_area"] == 5e-5
    assert manifest["point_count"] == 5
    assert manifest["promotion_boundary"].startswith("diagnostic_only")

    point_header = (EVIDENCE / "point_results.csv").read_text(encoding="utf-8").splitlines()[0]
    assert "strict_solver_tol" in point_header
    assert "strict_maxwell_area_gate" in point_header
    assert "verdict" in point_header
