import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
EVIDENCE = ROOT / "docs" / "analysis" / "pnjl" / "cep_maxwell" / "maxwell_contracts" / "pnjl_maxwell_endpoint_candidate_feasibility_v1"


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
