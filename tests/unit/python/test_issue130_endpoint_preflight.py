import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "pnjl" / "issue130_endpoint_preflight.py"


def module():
    spec = importlib.util.spec_from_file_location("issue130_endpoint_preflight", SCRIPT)
    assert spec is not None and spec.loader is not None
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


def gap(xi="0.0", cep="300.0", last="280.0"):
    return {
        "xi": xi,
        "muq_CEP_proxy_MeV": cep,
        "T_CEP_bracket_low_MeV": "120.0",
        "T_CEP_bracket_high_MeV": "120.1",
        "last_retained_muq_MeV": last,
        "last_retained_T_MeV": "121.0",
        "native_mu_gap_MeV": "20.0",
    }


def test_crossover_targets_are_strictly_below_cep():
    targets, summary = module().build_crossover([gap()])
    assert summary["target_count"] == 2
    assert all(280.0 < row["target_mu_MeV"] < 300.0 for row in targets)


def test_invalid_crossover_interval_stops_preflight():
    with pytest.raises(ValueError, match="no valid crossover target interval"):
        module().build_crossover([gap(cep="280.0")])


def test_missing_maxwell_boundary_is_not_pilot_target():
    maxwell = [{
        "xi": "0.0", "T_MeV": "119.0", "mu_MeV": "290.0",
        "grid_status": "rho:rho_geometry_not_converged", "grid_unresolved": "True",
    }]
    grid = [{
        "axis": "rho", "xi": "0.0", "T_MeV": "118.0", "reason": "rho_geometry_not_converged",
        "boundary_row_present": "False",
    }]
    targets, summary = module().build_maxwell([gap()], maxwell, grid)
    assert summary["candidate_count"] == 2
    assert sum(row["pilot_selection"] == "pilot_candidate" for row in targets) == 1
    assert sum(row["pilot_selection"] == "input_incomplete" for row in targets) == 1
