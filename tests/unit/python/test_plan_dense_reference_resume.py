import importlib.util
import json
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_PATH = PROJECT_ROOT / "scripts" / "pnjl" / "plan_dense_reference_resume.py"
TAG = "resume-test"
SHA = "a" * 40


def load_module():
    spec = importlib.util.spec_from_file_location("dense_reference_resume", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def base_config():
    return {
        "tag": TAG,
        "mode": "production",
        "compute_crossover": True,
        "crossover_method": "peak",
        "crossover_variable": "phi_u",
        "adaptive_xi": True,
        "adaptive_temperature": True,
        "rho_geometry_convergence": True,
        "crossover_only": False,
        "crossover_mu0_only": False,
        "T_min_MeV": 1.0,
        "T_max_MeV": 240.0,
        "T_step_MeV": 5.0,
        "temperature_max_refine_level": 3,
        "temperature_position_tol_MeV": 0.05,
        "temperature_density_tol": 0.005,
        "temperature_maxwell_area_tol": 5e-5,
        "rho_min": 0.0,
        "rho_max": 4.0,
        "rho_step": 0.05,
        "rho_position_tol_MeV": 0.025,
        "rho_density_tol": 0.0025,
        "rho_maxwell_area_tol": 5e-5,
        "cep_tol_MeV": 0.05,
        "p_num": 24,
        "t_num": 8,
        "thermo_quadrature_policy": "rs_reduced_adaptive",
        "thermo_quadrature_rtol": 1e-10,
        "thermo_quadrature_atol": 1e-12,
        "thermo_quadrature_maxevals": 20_000_000,
        "iterations": 80,
        "crossover_n_mu": 16,
        "crossover_mu_max_MeV": 450.0,
        "crossover_T_max_MeV": 240.0,
        "xi_max_refine_level": 2,
        "xi_position_tol_MeV": 0.05,
        "xi_density_tol": 0.005,
        "xi_maxwell_area_tol": 5e-5,
        "xi_response_rtol": 0.025,
        "xi_values": [],
        "requested_xi_values": [],
    }


def write_source_manifests(root: Path, values, *, config=None, commits=None):
    config = config or base_config()
    commits = commits or [SHA] * len(values)
    for index, (xi, commit) in enumerate(zip(values, commits)):
        payload = {
            "generator": {"git_commit": commit},
            "config": {**config, "xi_values": [xi], "requested_xi_values": [xi]},
        }
        path = root / f"shard-{index}" / f"phase_reference_{TAG}_manifest.json"
        path.parent.mkdir(parents=True)
        path.write_text(json.dumps(payload), encoding="utf-8")


def test_builds_plan_for_anchors_and_level1_midpoints(tmp_path: Path):
    module = load_module()
    write_source_manifests(tmp_path, [-0.5, -0.25, 0.0, 0.25, 0.5])

    plan = module.build_resume_plan(tmp_path, TAG, SHA, "-0.5,0,0.5")

    assert plan["source_shard_count"] == 5
    assert plan["expected_xi_csv"] == "-0.5,0,0.5"
    assert plan["initial_intervals"] == [["-0.5", "0"], ["0", "0.5"]]
    assert plan["common_args"][:6] == ["--T-min", "1", "--T-max", "240", "--T-step", "5"]
    assert "--p-num" in plan["common_args"]
    assert "--thermo-quadrature-policy" in plan["common_args"]


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"source_calculation_sha": "b" * 40}, "calculation commits differ"),
        ({"values": [-0.5, 0.0, 0.5]}, "not an initial staged grid"),
    ],
)
def test_rejects_sha_or_grid_mismatch(tmp_path: Path, kwargs, message):
    module = load_module()
    values = kwargs.get("values", [-0.5, -0.25, 0.0, 0.25, 0.5])
    write_source_manifests(tmp_path, values)
    with pytest.raises(ValueError, match=message):
        module.build_resume_plan(tmp_path, TAG, kwargs.get("source_calculation_sha", SHA), "-0.5,0,0.5")


def test_rejects_configuration_mismatch(tmp_path: Path):
    module = load_module()
    config = base_config()
    write_source_manifests(tmp_path, [-0.5, -0.25, 0.0, 0.25, 0.5], config=config)
    second = tmp_path / "shard-1" / f"phase_reference_{TAG}_manifest.json"
    payload = json.loads(second.read_text(encoding="utf-8"))
    payload["config"]["p_num"] = 32
    second.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(ValueError, match="configuration mismatch"):
        module.build_resume_plan(tmp_path, TAG, SHA, "-0.5,0,0.5")


def test_rejects_duplicate_xi(tmp_path: Path):
    module = load_module()
    write_source_manifests(tmp_path, [-0.5, -0.25, 0.0, 0.25, 0.25])

    with pytest.raises(ValueError, match="duplicate xi"):
        module.build_resume_plan(tmp_path, TAG, SHA, "-0.5,0,0.5")


def test_rejects_single_anchor(tmp_path: Path):
    module = load_module()
    write_source_manifests(tmp_path, [0.0])

    with pytest.raises(ValueError, match="at least two anchors"):
        module.build_resume_plan(tmp_path, TAG, SHA, "0")
