import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "pnjl_cep_hybrid_stagec_offline_feasibility.py"


def load_module():
    spec = importlib.util.spec_from_file_location("stagec_offline", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_candidate_runs_require_stable_plus_minus_plus_topology():
    module = load_module()
    monotone = [(0.0, 0.0), (1.0, 1.0), (2.0, 2.0), (3.0, 3.0)]
    assert module._candidate_runs(monotone) == []
    s_shape = [(0.0, 0.0), (1.0, 1.0), (2.0, 0.0), (3.0, -1.0), (4.0, 0.0), (5.0, 1.0)]
    candidates = module._candidate_runs(s_shape)
    assert len(candidates) == 1
    assert candidates[0]["rho_low"] == 1.0
    assert candidates[0]["rho_high"] == 3.0


def test_multiple_candidates_are_not_silently_collapsed():
    module = load_module()
    first = {"rho_low": 0.5, "rho_high": 1.0, "negative_secants": 2, "slope_before": 1.0, "slope_after": 1.0}
    second = {"rho_low": 2.0, "rho_high": 2.5, "negative_secants": 3, "slope_before": 1.0, "slope_after": 1.0}
    merged = module._merge_candidates([first], [second])
    assert len(merged) == 2


def test_local_selection_is_deterministic_and_capped():
    module = load_module()
    points = [(index / 16.0, float(index)) for index in range(65)]
    selected_a = module._select_local_points(points, [0.5, 2.0, 3.5], 12)
    selected_b = module._select_local_points(points, [0.5, 2.0, 3.5], 12)
    assert selected_a == selected_b
    assert len(selected_a) == 12
    assert len({rho for rho, _ in selected_a}) == 12


def test_classification_does_not_confirm_first_order_without_geometry():
    module = load_module()
    points = [(0.0, 0.0), (1.0, 1.0), (2.0, 0.0), (3.0, -1.0), (4.0, 0.0), (5.0, 1.0)]
    row = {
        "result_status": "confirmed_first_order",
        "geometry_converged": "false",
        "mu_transition_MeV": "NaN",
        "rho_hadron": "NaN",
        "rho_quark": "NaN",
        "mu_spinodal_hadron_MeV": "NaN",
        "mu_spinodal_quark_MeV": "NaN",
        "rho_spinodal_hadron": "NaN",
        "rho_spinodal_quark": "NaN",
        "area_residual": "NaN",
    }
    status, reason, candidates = module._classify(points, row)
    assert status == "ambiguous_near_critical"
    assert reason == "maxwell_or_geometry_unresolved"
    assert len(candidates) == 1
