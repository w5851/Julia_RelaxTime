from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "pnjl" / "compare_phase_reference_convergence.py"
BUILDER = ROOT / "scripts" / "analysis" / "pnjl" / "build_c2_convergence_audit.py"


def load_module():
    spec = importlib.util.spec_from_file_location("compare_phase_reference_convergence", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_builder():
    spec = importlib.util.spec_from_file_location("build_c2_convergence_audit", BUILDER)
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


def _phase_rows():
    return {
        "boundary": {"rows": []},
        "spinodals": {"rows": []},
        "cep": {"rows": []},
        "crossover": {"rows": []},
    }


def test_public_anchor_table_distinguishes_first_order_regression_from_monotone():
    module = load_module()
    candidate = _phase_rows()
    reference = _phase_rows()
    reference["boundary"]["rows"] = [
        {
            "xi": "-0.35",
            "T_MeV": "51",
            "mu_transition_MeV": "300",
            "rho_hadron": "0.1",
            "rho_quark": "2.0",
            "area_residual": "1e-6",
        }
    ]
    table = module.build_public_anchor_state_table(candidate, reference)
    row = next(item for item in table if item["xi"] == -0.35 and item["T_MeV"] == 51)
    assert row["reference_state"] == module.STATE_FIRST_ORDER
    assert row["candidate_state"] == module.STATE_MONOTONE
    assert row["classification_regression"] == "true"

    candidate["grid_convergence"] = {
        "rows": [{"xi": "-0.35", "T_MeV": "51", "axis": "rho"}]
    }
    table = module.build_public_anchor_state_table(candidate, reference)
    row = next(item for item in table if item["xi"] == -0.35 and item["T_MeV"] == 51)
    assert row["candidate_state"] == module.STATE_AMBIGUOUS


def test_crossover_comparison_uses_physical_mu_union_and_interpolation():
    module = load_module()
    candidate = [
        {"xi": "0.2875", "mu_MeV": "0", "T_crossover_MeV": "180", "rho": "0", "derivative": "10", "converged": "true"},
        {"xi": "0.2875", "mu_MeV": "100", "T_crossover_MeV": "170", "rho": "1", "derivative": "20", "converged": "true"},
    ]
    reference = [
        {"xi": "0.2875", "mu_MeV": "0", "T_crossover_MeV": "180", "rho": "0", "derivative": "10", "converged": "true"},
        {"xi": "0.2875", "mu_MeV": "50", "T_crossover_MeV": "175", "rho": "0.5", "derivative": "15", "converged": "true"},
        {"xi": "0.2875", "mu_MeV": "100", "T_crossover_MeV": "170", "rho": "1", "derivative": "20", "converged": "true"},
    ]
    rows = module.compare_crossover_artifact(candidate, reference)
    assert {row["match_key"] for row in rows if row["metric"] == "rho"} == {
        "0.2875|0",
        "0.2875|50",
        "0.2875|100",
    }
    assert all(row["match_status"] in {"matched", "interpolated"} for row in rows)


def test_cep_gate_separates_endpoint_resolution_from_ambiguity_width():
    module = load_module()
    candidate = _phase_rows()
    reference = _phase_rows()
    row = {
        "xi": "0.225",
        "result_status": "ambiguous",
        "T_bracket_low_MeV": "130",
        "T_bracket_high_MeV": "130.25",
        "bracket_width_T_MeV": "0.25",
        "ambiguity_width_T_MeV": "0.125",
        "temperature_resolution_target_MeV": "0.0625",
    }
    candidate["cep"]["rows"] = [row]
    reference["cep"]["rows"] = [dict(row, bracket_width_T_MeV="0.0625", ambiguity_width_T_MeV="0.0625")]
    gate = module.build_cep_gate_table(candidate, reference)[0]
    assert gate["candidate_endpoint_resolution_MeV"] == 0.0625
    assert gate["candidate_ambiguity_width_T_MeV"] == 0.125
    assert gate["candidate_bracket_pass"] is False


def test_verdict_precedence_keeps_secondary_failures():
    module = load_module()
    anchor = [{"reference_state": module.STATE_FIRST_ORDER, "candidate_state": module.STATE_AMBIGUOUS, "classification_regression": "true"}]
    cep = [{"candidate_bracket_pass": False, "reference_bracket_pass": True}]
    geometry = [{"all_pass": False}]
    verdict, secondary = module.determine_verdict(
        bad_inventory=[{"artifact": "boundary"}],
        anchor_rows=anchor,
        cep_rows=cep,
        geometry_rows=geometry,
        comparison_rows=[],
    )
    assert verdict == "artifact_invalid"
    assert secondary == ["artifact_invalid", "classification_regression", "cep_bracket_unresolved", "first_order_geometry_unstable"]


def test_comparison_exception_uses_existing_geometry_and_crossover_gates():
    module = load_module()
    builder = load_builder()
    geometry_row = {
        "artifact": "boundary",
        "match_status": "matched",
        "metric": "area_residual",
        "abs_diff": "0.00006",
        "rel_diff": "0.1",
    }
    crossover_row = {
        "artifact": "crossover",
        "match_status": "interpolated",
        "metric": "derivative",
        "abs_diff": "1",
        "rel_diff": "0.026",
    }
    assert builder.comparison_exception(geometry_row, module) == "geometry_gate_exceeded"
    assert builder.comparison_exception(crossover_row, module) == "crossover_response_gate_exceeded"
    assert builder.comparison_exception({"match_status": "missing_in_candidate"}, module) == "missing_in_candidate"
