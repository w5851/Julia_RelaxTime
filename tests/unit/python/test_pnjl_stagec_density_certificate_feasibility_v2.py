import csv
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
EVALUATOR = ROOT / "scripts" / "analysis" / "pnjl_stagec_density_certificate_feasibility_v2.jl"
PLOTTER = ROOT / "scripts" / "analysis" / "plot_pnjl_stagec_density_certificate_feasibility_v2.py"
WORKFLOW = ROOT / ".github" / "workflows" / "pnjl-stagec-density-certificate-feasibility-v2.yml"


def test_v2_evaluator_is_solver_free_and_versioned():
    text = EVALUATOR.read_text(encoding="utf-8")
    assert 'SCHEMA_VERSION = "cep_hybrid_stagec_density_certificate_feasibility_v2"' in text
    assert 'const SOURCE_RUN_ID = "31296511813"' in text
    assert 'const SOURCE_CALCULATION_SHA = "ffa816df0a145f73d7490db1ed9ff10c92e017a4"' in text
    assert "solver_called" in text
    assert "run_production_phase_pipeline" not in text
    assert "candidate_point_index.csv" in text
    assert "Write the manifest last" in text


def test_v2_workflow_is_aggregate_only_and_uses_dual_sha():
    text = WORKFLOW.read_text(encoding="utf-8")
    assert "options: [aggregate_replay]" in text
    assert "source_run_id" in text and "source_calculation_sha" in text
    assert "source_postprocess_sha" in text
    assert 'EXPECTED_SOURCE_RUN_ID: "31296511813"' in text
    assert 'EXPECTED_CALCULATION_SHA: "ffa816df0a145f73d7490db1ed9ff10c92e017a4"' in text
    assert "numeric_jobs" in text and "test \"$numeric_jobs\" = 30" in text
    assert "numeric_failures" in text and "test \"$numeric_failures\" = 0" in text
    assert "git worktree add --detach calculation" in text
    assert "--expected-calculation-sha" in text
    assert "--producer-head-sha" in text
    assert "solver_called" in text
    assert "matplotlib==3.9.2" in text


def test_v2_plotter_emits_figures_and_refreshes_manifest(tmp_path):
    spec = importlib.util.spec_from_file_location("stagec_v2_plotter", PLOTTER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    source = tmp_path / "source" / "job"
    source.mkdir(parents=True)
    rows = [
        {"method": "independent_oracle", "xi": "-0.35", "T_MeV": "51.0", "rho_level": "0", "rho": "0", "muq_MeV": "300"},
        {"method": "independent_oracle", "xi": "-0.35", "T_MeV": "51.0", "rho_level": "0", "rho": "1", "muq_MeV": "301"},
        {"method": "independent_oracle", "xi": "-0.35", "T_MeV": "51.0", "rho_level": "0", "rho": "2", "muq_MeV": "300"},
        {"method": "independent_oracle", "xi": "-0.35", "T_MeV": "51.0", "rho_level": "0", "rho": "3", "muq_MeV": "302"},
    ]
    with (source / "curve_points.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    with (source / "slice_metrics.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["method", "xi", "T_MeV", "mu_transition_MeV"])
        writer.writeheader()
        writer.writerow({"method": "independent_oracle", "xi": "-0.35", "T_MeV": "51.0", "mu_transition_MeV": "300.5"})
    output = tmp_path / "output"
    output.mkdir()
    (output / "manifest.json").write_text(json.dumps({"schema_version": "cep_hybrid_stagec_density_certificate_feasibility_v2"}), encoding="utf-8")
    (output / "route_comparison.csv").write_text("route,cap,xi,T_MeV,maxwell_area_gate\nstage_b_features_v1,12,-0.35,51.0,1e-6\n", encoding="utf-8")
    (output / "selected_point_index.csv").write_text("xi,T_MeV,route,cap,rank,feature,rho,muq_MeV\n", encoding="utf-8")
    result = module.render(tmp_path / "source", output, output)
    assert len(result["figures"]) == 1
    figure = output / result["figures"][0]["path"]
    assert figure.is_file() and figure.stat().st_size > 0
    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    assert "plot_manifest_sha256" in manifest
