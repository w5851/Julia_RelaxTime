import csv
import hashlib
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
COLLECTOR = ROOT / "scripts" / "analysis" / "collect_pnjl_cep_hybrid_production_shadow.py"
WORKFLOW = ROOT / ".github" / "workflows" / "pnjl-cep-hybrid-production-shadow.yml"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _status(t: float):
    if t in (5.0, 20.0):
        return "confirmed_first_order"
    if t in (60.0, 100.0):
        return "ambiguous_near_critical"
    return "confirmed_monotone"


def write_job(root: Path, xi: float, method: str, *, cost: int = 100):
    collector = load_module(COLLECTOR, f"hybrid_collector_constants_{xi}_{method}")
    job = root / f"job-{method}-{xi}"
    job.mkdir()
    anchors = collector.ANCHORS[xi]
    sha = "a" * 40
    with (job / "curve_points.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["xi", "method", "T_MeV", "rho_level", "rho", "muq_MeV", "converged", "finite"])
        writer.writeheader()
        for index, t in enumerate(anchors):
            writer.writerow({"xi": xi, "method": method, "T_MeV": t, "rho_level": 0, "rho": 1.0, "muq_MeV": 300.0 + index, "converged": "true", "finite": "true"})
    with (job / "slice_metrics.csv").open("w", newline="", encoding="utf-8") as handle:
        fields = ["xi", "method", "T_MeV", "stage_a_status", "stage_b_status", "stage_c_status", "stage_used", "upgrade_reason", "result_status", "geometry_converged", "position_error_MeV", "density_error", "maxwell_area_gate", "area_residual", "mu_transition_MeV", "rho_hadron", "rho_quark", "mu_spinodal_hadron_MeV", "mu_spinodal_quark_MeV", "rho_spinodal_hadron", "rho_spinodal_quark", "support_low", "support_high", "support_point_count", "targeted_additions", "solver_failure_count"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for t in anchors:
            status = _status(t)
            first = status == "confirmed_first_order"
            writer.writerow({
                "xi": xi, "method": method, "T_MeV": t,
                "stage_a_status": status, "stage_b_status": "not_run", "stage_c_status": "not_run",
                "stage_used": "stage_a", "upgrade_reason": "synthetic", "result_status": status,
                "geometry_converged": "true", "position_error_MeV": 0.01, "density_error": 0.001,
                "maxwell_area_gate": 0.00001, "area_residual": 0.00001,
                "mu_transition_MeV": 300.0 if first else "", "rho_hadron": 1.0 if first else "",
                "rho_quark": 2.0 if first else "", "mu_spinodal_hadron_MeV": 301.0 if first else "",
                "mu_spinodal_quark_MeV": 299.0 if first else "", "rho_spinodal_hadron": 0.9 if first else "",
                "rho_spinodal_quark": 2.1 if first else "", "support_low": "", "support_high": "",
                "support_point_count": 0, "targeted_additions": 9 if method == "production_hybrid" else 0,
                "solver_failure_count": 0,
            })
    (job / "method_costs.csv").write_text(
        "xi,method,unique_solves,equilibrium_requests,fixedrho_requests,residual_calls,jacobian_calls,newton_iterations,runner_seconds,fallback_count,retry_count\n"
        f"{xi},{method},{cost},{cost},{cost},100,50,120,10,0,0\n", encoding="utf-8")
    (job / "cep_accuracy.csv").write_text(
        "xi,method,anchor_T_MeV,result_status,T_last_first_order_MeV,muq_last_first_order_MeV,T_first_monotone_MeV,ambiguity_width_T_MeV,temperature_resolution_target_MeV\n"
        + "".join(f"{xi},{method},{t},{_status(t)},20,300,120,100,0.125\n" for t in anchors), encoding="utf-8")
    hashes = {name: hashlib.sha256((job / name).read_bytes()).hexdigest() for name in ("curve_points.csv", "slice_metrics.csv", "method_costs.csv", "cep_accuracy.csv")}
    (job / "job_summary.json").write_text(json.dumps({
        "schema_version": "cep_cascade_production_shadow_v2", "xi": xi, "method": method,
        "calculation_sha": sha, "workflow_head_sha": sha, "anchors": list(anchors),
        "finite_and_converged_final": True, "curve_file_sha256": hashes,
        "provenance": {"calculation_sha": sha, "reference_write": False, "anchor_run_count": len(anchors)},
    }), encoding="utf-8")


def _matrix(tmp_path: Path):
    for xi in (-0.5, 0.0, 0.5):
        write_job(tmp_path, xi, "production_hybrid", cost=90)
        write_job(tmp_path, xi, "memoized_dense", cost=100)
        write_job(tmp_path, xi, "independent_oracle", cost=110)


def test_hybrid_collector_accepts_complete_matrix(tmp_path):
    module = load_module(COLLECTOR, "hybrid_collector_complete")
    _matrix(tmp_path)
    gate = module.collect(tmp_path, tmp_path / "aggregate", None, "a" * 40)
    assert gate["verdict"] == "full_hybrid_candidate"
    assert (tmp_path / "aggregate" / "geometry_accuracy.csv").is_file()
    assert (tmp_path / "aggregate" / "manifest.json").is_file()


def test_hybrid_collector_rejects_confirmation_on_oracle_ambiguous(tmp_path):
    module = load_module(COLLECTOR, "hybrid_collector_unsupported")
    _matrix(tmp_path)
    oracle = tmp_path / "job-independent_oracle--0.5" / "slice_metrics.csv"
    text = oracle.read_text(encoding="utf-8").replace(
        "stage_a,synthetic,confirmed_first_order,true",
        "stage_a,synthetic,ambiguous_near_critical,true",
        1,
    )
    oracle.write_text(text, encoding="utf-8")
    summary_path = oracle.parent / "job_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["curve_file_sha256"]["slice_metrics.csv"] = hashlib.sha256(oracle.read_bytes()).hexdigest()
    summary_path.write_text(json.dumps(summary), encoding="utf-8")
    gate = module.collect(tmp_path, tmp_path / "aggregate", None, "a" * 40)
    assert gate["verdict"] == "hybrid_integration_failed"
    assert gate["classification_errors"]


def test_hybrid_workflow_contract_has_immutable_matrix_and_replay():
    module = load_module(COLLECTOR, "hybrid_collector_workflow")
    text = WORKFLOW.read_text(encoding="utf-8")
    assert "calculation_ref" in text
    assert "production_hybrid" in text and "independent_oracle" in text
    assert "timeout-minutes: 180" in text
    assert "aggregate_replay" in text and "source_run_id" in text
    assert "matplotlib==3.9.2" in text
    assert len(module.XIS) * len(module.METHODS) == 9
