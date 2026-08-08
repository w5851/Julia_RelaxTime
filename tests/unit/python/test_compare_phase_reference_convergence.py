from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "pnjl" / "compare_phase_reference_convergence.py"


def load_module():
    spec = importlib.util.spec_from_file_location("compare_phase_reference_convergence", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_ambiguous_cep_compares_phase_endpoints_and_not_empty_cep_point():
    module = load_module()
    candidate = [
        {
            "xi": "0",
            "result_status": "ambiguous",
            "T_CEP_MeV": "",
            "muq_CEP_MeV": "",
            "muB_CEP_MeV": "",
            "T_bracket_low_MeV": "",
            "T_bracket_high_MeV": "",
            "bracket_width_T_MeV": "0.0625",
            "T_last_first_order_MeV": "130",
            "T_first_monotone_MeV": "131",
            "ambiguity_width_T_MeV": "1",
        }
    ]
    reference = [dict(candidate[0], T_last_first_order_MeV="130.01")]
    rows = module.compare_cep_artifact(candidate, reference)
    summary = module.summarize_comparison(rows)
    assert any(row["metric"] == "T_last_first_order_MeV" and row["matched_count"] == 1 for row in summary)
    assert any(row["metric"] == "T_CEP_MeV" and row["not_applicable_count"] == 1 for row in summary)
    assert not any(row["metric"] == "T_CEP_MeV" and row["non_numeric_count"] for row in summary)


def test_adaptive_xi_is_reported_separately_from_missing_match():
    module = load_module()
    candidate = [{"xi": "0.025", "T_CEP_MeV": "", "muq_CEP_MeV": ""}]
    reference = []
    rows = module.compare_cep_artifact(candidate, reference)
    assert rows[0]["match_status"] == "adaptive_xi_added"
    summary = module.summarize_comparison(rows)
    assert summary[0]["adaptive_xi_count"] == 1
    assert summary[0]["missing_count"] == 0
