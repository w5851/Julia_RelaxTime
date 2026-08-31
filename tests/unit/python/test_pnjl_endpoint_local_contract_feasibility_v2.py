import csv
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "pnjl_endpoint_local_feasibility_v2.py"
WORKFLOW = ROOT / "docs" / "analysis" / "governance" / "diagnostic_workflow_retirement_wave2_v1" / "definitions" / "pnjl-maxwell-endpoint-production-shadow.yml"


def load_module():
    spec = importlib.util.spec_from_file_location("endpoint_local_feasibility_v2", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _curve_rows(xi: float, temperature: float, values: list[tuple[float, float]], method: str = "production_hybrid"):
    return [
        {
            "xi": str(xi), "T_MeV": str(temperature), "method": method,
            "rho_level": str(index % 4), "rho": str(rho), "muq_MeV": str(mu),
            "finite": "true", "converged": "true",
        }
        for index, (rho, mu) in enumerate(values)
    ]


def _slice(xi: float, temperature: float, *, method: str = "production_hybrid", status: str = "ambiguous_near_critical"):
    return {
        "xi": str(xi), "T_MeV": str(temperature), "method": method,
        "result_status": status, "maxwell_candidate_count": "1",
        "maxwell_crossing_count": "3", "maxwell_endpoint_dependent": "true",
        "mu_transition_MeV": "330.0", "rho_hadron": "0.0015", "rho_quark": "2.5",
        "rho_spinodal_quark": "2.0", "area_residual": "1e-6",
        "geometry_converged": "true", "position_error_MeV": "0.01",
        "density_error": "0.001", "finite_and_converged": "true",
    }


def test_right_crossing_bracket_does_not_require_spinodal_overshoot():
    module = load_module()
    points = [(0.0, 300.0), (0.1, 340.0), (1.0, 320.0), (2.0, 310.0), (2.5, 330.0), (3.0, 335.0)]
    brackets = module._crossing_brackets(points, 330.0, minimum_rho=2.0)
    assert (2.0, 2.5) in brackets


def test_route_uses_stage_b_for_midpoint_selection_and_produces_internal_certificate():
    module = load_module()
    values = [
        (0.0, 300.0), (0.0125, 335.0), (0.5, 350.0), (1.0, 320.0),
        (1.5, 315.0), (2.0, 310.0), (2.5, 330.0), (3.0, 340.0),
    ]
    curves = _curve_rows(-0.5, 20.0, values)
    deep_values = [(0.0, 300.0), (0.0015625, 331.0), (0.003125, 334.0), *values[1:]]
    deep_curves = _curve_rows(-0.5, 20.0, deep_values, method="independent_oracle")
    source = _slice(-0.5, 20.0)
    deep = _slice(-0.5, 20.0, method="independent_oracle", status="confirmed_first_order")
    deep["rho_hadron"] = "0.0010"
    route, trace, brackets = module._route(curves, [source], deep_curves, [deep], -0.5, 20.0)
    assert route["status"] == "feasible"
    assert route["certificate_kind"] == "endpoint_local_geometry_first_order"
    assert route["right_crossing_bracketed"] is True
    assert len(trace) == module.MAX_TARGETED
    assert all(row["route_decision_source"] == "post_sampled_curve_candidate" for row in trace)
    assert brackets[1]["bracket_kind"] == "right_maxwell"


def test_workflow_has_authenticated_replay_and_explicit_run_mode():
    text = WORKFLOW.read_text(encoding="utf-8")
    assert "aggregate_replay" in text
    assert "GH_TOKEN: ${{ github.token }}" in text
    assert "--run-mode aggregate_replay" in text
    assert "deep_oracle_run_id" in text
