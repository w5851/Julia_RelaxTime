import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
COLLECTOR = ROOT / "scripts" / "analysis" / "collect_pnjl_cep_narrow_pilot_v2.py"
FREEZER = ROOT / "scripts" / "analysis" / "freeze_pnjl_cep_narrow_pilot_v2_windows.py"
PLOTTER = ROOT / "scripts" / "analysis" / "plot_pnjl_cep_narrow_pilot_v2.py"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_job(root: Path, xi: float, method: str, stage: str, *, oracle: bool = False) -> None:
    job = root / f"{stage}-{method}-{xi}"
    job.mkdir()
    low = 130.0 + (0.1 if oracle else 0.05)
    high = 131.0 + (0.1 if oracle else 0.05)
    summary = {
        "schema_version": "cep_narrow_pilot_v2",
        "xi": xi,
        "method": method,
        "stage": stage,
        "calculation_sha": "a" * 40,
        "finite_and_converged_final": True,
    }
    (job / "job_summary.json").write_text(json.dumps(summary), encoding="utf-8")
    (job / "curve_points.csv").write_text(
        "xi,method,T_MeV,rho_level,rho,muq_MeV,converged,finite\n"
        f"{xi},{method},130,1,1,295,true,true\n",
        encoding="utf-8",
    )
    (job / "slice_metrics.csv").write_text(
        "xi,method,T_MeV,coarse_status,fine_status,coarse_reason,fine_reason,result_status,targeted_additions,wall_seconds\n"
        f"{xi},{method},130,valid,valid,ok,ok,confirmed_first_order,0,1\n"
        f"{xi},{method},131,invalid,invalid,no_s_shape,no_s_shape,confirmed_monotone,0,1\n",
        encoding="utf-8",
    )
    (job / "method_costs.csv").write_text(
        "xi,method,unique_solves,equilibrium_requests,residual_calls,jacobian_calls,newton_iterations,runner_seconds,fallback_count,retry_count\n"
        f"{xi},{method},100,100,200,100,300,10,0,0\n",
        encoding="utf-8",
    )
    (job / "cep_accuracy.csv").write_text(
        "xi,method,result_status,fine_T_last_first_order_MeV,fine_muq_last_first_order_MeV,fine_T_first_monotone_MeV,"
        "coarse_T_last_first_order_MeV,coarse_muq_last_first_order_MeV,coarse_T_first_monotone_MeV,"
        "T_last_first_order_MeV,muq_last_first_order_MeV,T_first_monotone_MeV\n"
        f"{xi},{method},ambiguous,{low},295,{high},{low},295,{high},{low},295,{high}\n",
        encoding="utf-8",
    )
    if stage == "discovery":
        proposal = {
            "schema_version": "cep_narrow_pilot_v2_window_proposal",
            "xi": xi,
            "T_min_MeV": 120.0,
            "T_max_MeV": 140.0,
            "required_T_anchors": [126.0, 130.0, 134.0],
            "calculation_sha": "a" * 40,
            "cep": {"result_status": "ambiguous", "T_last_first_order_MeV": low, "T_first_monotone_MeV": high},
        }
        (job / "window_proposal.json").write_text(json.dumps(proposal), encoding="utf-8")


def test_freeze_windows_requires_three_xi(tmp_path):
    module = load_module(FREEZER, "freeze_v2")
    for xi in (-0.5, 0.0, 0.5):
        job = tmp_path / f"job-{xi}"
        job.mkdir()
        proposal = {
            "schema_version": "cep_narrow_pilot_v2_window_proposal",
            "xi": xi,
            "T_min_MeV": 120.0,
            "T_max_MeV": 140.0,
            "required_T_anchors": [130.0],
            "calculation_sha": "a" * 40,
            "cep": {"result_status": "ambiguous"},
        }
        (job / "window_proposal.json").write_text(json.dumps(proposal), encoding="utf-8")
    payload = module.freeze(tmp_path, tmp_path / "windows.json")
    assert payload["schema_version"] == "cep_narrow_pilot_v2_validation_windows"
    assert len(payload["windows"]) == 3
    assert payload["calculation_sha"] == "a" * 40


def test_collector_v2_accepts_complete_diagnostic_matrix(tmp_path):
    module = load_module(COLLECTOR, "collector_v2")
    for xi in sorted(module.XIS):
        write_job(tmp_path, xi, "rho_support_cascade", "discovery")
        write_job(tmp_path, xi, "c2_dense_baseline", "validation")
        write_job(tmp_path, xi, "high_resolution_oracle", "validation", oracle=True)
    gate = module.collect(tmp_path, tmp_path / "aggregate")
    assert gate["status"] == "pass"
    assert (tmp_path / "aggregate" / "manifest.json").is_file()
    assert (tmp_path / "aggregate" / "claim_ledger.csv").is_file()


def test_collector_v2_reports_contract_failure_for_missing_validation(tmp_path):
    module = load_module(COLLECTOR, "collector_v2_missing")
    for xi in sorted(module.XIS):
        write_job(tmp_path, xi, "rho_support_cascade", "discovery")
    gate = module.collect(tmp_path, tmp_path / "aggregate")
    assert gate["status"] == "workflow_failure"
    assert any("missing validation jobs" in error for error in gate["workflow_contract_errors"])


def test_plotter_v2_smoke(tmp_path):
    collector = load_module(COLLECTOR, "collector_v2_plot")
    for xi in sorted(collector.XIS):
        write_job(tmp_path, xi, "rho_support_cascade", "discovery")
        write_job(tmp_path, xi, "c2_dense_baseline", "validation")
        write_job(tmp_path, xi, "high_resolution_oracle", "validation", oracle=True)
    aggregate = tmp_path / "aggregate"
    collector.collect(tmp_path, aggregate)
    plotter = load_module(PLOTTER, "plotter_v2")
    plotter.plot(aggregate, aggregate / "figures")
    assert (aggregate / "figures" / "cep_intervals.png").is_file()
    assert (aggregate / "figures" / "plot_manifest.json").is_file()
