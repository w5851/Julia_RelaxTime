import csv
import hashlib
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
COLLECTOR = ROOT / "scripts" / "analysis" / "collect_pnjl_cep_deep_oracle.py"
WORKFLOW = ROOT / ".github" / "workflows" / "pnjl-cep-deep-oracle.yml"
SHA = "a" * 40
WORKFLOW_SHA = "b" * 40
POINTS = [(-0.5, 20.0), (0.0, 5.0), (0.0, 20.0), (0.5, 5.0), (0.5, 20.0)]


def load_module(path: Path):
    spec = importlib.util.spec_from_file_location("deep_oracle_collector", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_job(root: Path, xi: float, temperature: float):
    job = root / f"job-xi-{xi}-T-{temperature}"
    job.mkdir()
    curve = job / "curve_points.csv"
    with curve.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["xi", "method", "T_MeV", "rho_level", "rho", "muq_MeV", "residual_norm", "converged", "finite"],
        )
        writer.writeheader()
        writer.writerow({
            "xi": xi,
            "method": "independent_oracle",
            "T_MeV": temperature,
            "rho_level": 0,
            "rho": 1.0,
            "muq_MeV": 300.0,
            "residual_norm": 1e-10,
            "converged": "true",
            "finite": "true",
        })
    (job / "slice_metrics.csv").write_text("xi,T_MeV,result_status\n" f"{xi},{temperature},confirmed_first_order\n", encoding="utf-8")
    (job / "method_costs.csv").write_text("xi,T_MeV,unique_solves\n" f"{xi},{temperature},100\n", encoding="utf-8")
    (job / "cep_accuracy.csv").write_text("xi,T_MeV,result_status\n" f"{xi},{temperature},confirmed_first_order\n", encoding="utf-8")
    (job / "job_summary.json").write_text(json.dumps({
        "schema_version": "cep_deep_oracle_v1",
        "xi": xi,
        "temperature_MeV": temperature,
        "calculation_sha": SHA,
        "workflow_head_sha": WORKFLOW_SHA,
        "rho_coarse_step": 0.003125,
        "rho_fine_step": 0.0015625,
        "curve_file_sha256": {"curve_points.csv": hashlib.sha256(curve.read_bytes()).hexdigest()},
    }), encoding="utf-8")


def test_deep_oracle_collector_requires_five_points(tmp_path):
    module = load_module(COLLECTOR)
    for xi, temperature in POINTS:
        write_job(tmp_path, xi, temperature)
    result = module.collect(tmp_path, tmp_path / "aggregate", SHA, WORKFLOW_SHA, "test")
    assert result["status"] == "complete"
    assert result["job_count"] == 5
    assert (tmp_path / "aggregate" / "manifest.json").is_file()
    assert (tmp_path / "aggregate" / "claim_ledger.csv").is_file()


def test_deep_oracle_workflow_is_fixed_scope_and_provenance_bound():
    text = WORKFLOW.read_text(encoding="utf-8")
    runner = (ROOT / "scripts" / "analysis" / "pnjl_cep_deep_oracle.jl").read_text(encoding="utf-8")
    assert "calculation_ref" in text
    assert "workflow-head-sha" in text
    assert "0.003125" in runner
    assert "0.0015625" in runner
    assert text.count("temperature:") == 5
    assert "pnjl_cep_deep_oracle.jl" in text
    assert "reference_write" not in text
