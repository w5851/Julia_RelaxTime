import csv
import hashlib
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
COLLECTOR = ROOT / "scripts" / "analysis" / "collect_pnjl_cep_production_shadow.py"
WORKFLOW = ROOT / ".github" / "workflows" / "pnjl-cep-production-shadow.yml"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _status(t: float):
    if t in (5.0, 20.0):
        return "confirmed_first_order", "valid", "valid", "ok", "ok"
    if t in (60.0,):
        return "ambiguous_near_critical", "weak_s_shape", "invalid", "maxwell_failed", "no_s_shape"
    return "confirmed_monotone", "invalid", "invalid", "no_s_shape", "no_s_shape"


def write_job(root: Path, xi: float, method: str, *, cost: int = 100):
    collector = load_module(COLLECTOR, f"collector_constants_{xi}_{method}")
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
        fields = ["xi", "method", "T_MeV", "coarse_status", "fine_status", "coarse_reason", "fine_reason", "result_status", "geometry_converged", "position_error_MeV", "density_error", "maxwell_area_gate", "targeted_additions", "solver_failure_count"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for t in anchors:
            status, coarse, fine, coarse_reason, fine_reason = _status(t)
            writer.writerow({"xi": xi, "method": method, "T_MeV": t, "coarse_status": coarse, "fine_status": fine, "coarse_reason": coarse_reason, "fine_reason": fine_reason, "result_status": status, "geometry_converged": "true", "position_error_MeV": 0.01, "density_error": 0.001, "maxwell_area_gate": 0.00001, "targeted_additions": 9 if method == "production_cascade" else 0, "solver_failure_count": 0})

    (job / "method_costs.csv").write_text(
        "xi,method,unique_solves,equilibrium_requests,fixedrho_requests,residual_calls,jacobian_calls,newton_iterations,runner_seconds,fallback_count,retry_count\n"
        f"{xi},{method},{cost},{cost},{cost},100,50,120,10,0,0\n",
        encoding="utf-8",
    )
    (job / "cep_accuracy.csv").write_text(
        "xi,method,anchor_T_MeV,result_status,T_last_first_order_MeV,muq_last_first_order_MeV,T_first_monotone_MeV,ambiguity_width_T_MeV,temperature_resolution_target_MeV\n"
        + "".join(f"{xi},{method},{t},ambiguous,20,300,60,40,0.125\n" for t in anchors),
        encoding="utf-8",
    )
    hashes = {name: hashlib.sha256((job / name).read_bytes()).hexdigest() for name in ("curve_points.csv", "slice_metrics.csv", "method_costs.csv", "cep_accuracy.csv")}
    (job / "job_summary.json").write_text(json.dumps({
        "schema_version": "cep_cascade_production_shadow_v1",
        "xi": xi,
        "method": method,
        "calculation_sha": sha,
        "workflow_head_sha": sha,
        "anchors": list(anchors),
        "finite_and_converged_final": True,
        "curve_file_sha256": hashes,
        "provenance": {"calculation_sha": sha, "reference_write": False, "anchor_run_count": len(anchors)},
    }), encoding="utf-8")


def test_shadow_collector_accepts_complete_matrix(tmp_path):
    module = load_module(COLLECTOR, "collector_complete")
    for xi in sorted(module.XIS):
        write_job(tmp_path, xi, "production_cascade", cost=90)
        write_job(tmp_path, xi, "memoized_dense", cost=100)
        write_job(tmp_path, xi, "independent_oracle", cost=110)
    gate = module.collect(tmp_path, tmp_path / "aggregate", None, "a" * 40)
    assert gate["verdict"] == "full_cascade_candidate"
    assert (tmp_path / "aggregate" / "manifest.json").is_file()
    assert (tmp_path / "aggregate" / "geometry_accuracy.csv").is_file()


def test_shadow_collector_rejects_duplicate_curve_key(tmp_path):
    module = load_module(COLLECTOR, "collector_duplicate")
    write_job(tmp_path, -0.5, "production_cascade")
    curve = tmp_path / "job-production_cascade--0.5" / "curve_points.csv"
    with curve.open("a", encoding="utf-8") as handle:
        handle.write("-0.5,production_cascade,5,0,1,300,true,true\n")
    # Recompute the declared hash so this test isolates duplicate-key validation.
    summary_path = curve.parent / "job_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["curve_file_sha256"]["curve_points.csv"] = hashlib.sha256(curve.read_bytes()).hexdigest()
    summary_path.write_text(json.dumps(summary), encoding="utf-8")
    gate = module.collect(tmp_path, tmp_path / "aggregate", None, "a" * 40)
    assert gate["verdict"] == "workflow_failure"
    assert any("duplicate curve point key" in error for error in gate["workflow_contract_errors"])


def test_shadow_collector_marks_low_temperature_cascade_failure_hybrid(tmp_path):
    module = load_module(COLLECTOR, "collector_hybrid")
    for xi in sorted(module.XIS):
        write_job(tmp_path, xi, "production_cascade", cost=90)
        write_job(tmp_path, xi, "memoized_dense", cost=100)
        write_job(tmp_path, xi, "independent_oracle", cost=110)
    job = tmp_path / "job-production_cascade--0.5"
    path = job / "slice_metrics.csv"
    text = path.read_text(encoding="utf-8").replace("-0.5,production_cascade,5.0,valid,valid,ok,ok,confirmed_first_order,true", "-0.5,production_cascade,5.0,valid,valid,ok,ok,confirmed_first_order,false")
    path.write_text(text, encoding="utf-8")
    summary_path = job / "job_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["curve_file_sha256"]["slice_metrics.csv"] = hashlib.sha256(path.read_bytes()).hexdigest()
    summary_path.write_text(json.dumps(summary), encoding="utf-8")
    gate = module.collect(tmp_path, tmp_path / "aggregate", None, "a" * 40)
    assert gate["verdict"] == "hybrid_required"


def test_shadow_workflow_has_immutable_ref_and_nine_jobs():
    module = load_module(COLLECTOR, "collector_workflow_contract")
    text = WORKFLOW.read_text(encoding="utf-8")
    assert "calculation_ref" in text
    assert "actions: read" in text
    assert "timeout-minutes: 180" in text
    assert len(module.XIS) * len(module.METHODS) == 9
