import importlib.util
import csv
import json
import shutil
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "pnjl_c2_limited_feasibility.py"
WORKFLOW = ROOT / "docs" / "analysis" / "governance" / "diagnostic_workflow_retirement_wave2_v1" / "definitions" / "pnjl-c2-limited-feasibility.yml"


def load_module():
    spec = importlib.util.spec_from_file_location("c2_limited", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_schema_and_frozen_matrix_are_explicit():
    module = load_module()
    assert module.SCHEMA_VERSION == "pnjl_c2_limited_feasibility_v1"
    assert len(module.XI_GRID) == 10
    assert module.RHO_COUNT == 1281
    assert sum(len(temperatures) for temperatures in module.DENSITY_ANCHORS_BY_XI.values()) == 15
    assert len(module.ALL_ANCHORS) == 15
    assert module.ROUTES == (
        "stage_b_features_v1",
        "balanced_density_features_v2",
        "geometry_feedback_v2",
    )
    assert module.CAPS == (12, 16, 24)
    assert module.CROSSOVER_XI == 0.2875
    assert len(module.CROSSOVER_MU) == 5


def test_scope_plan_freezes_cep_and_crossover_contract(tmp_path):
    module = load_module()
    audit = tmp_path / "audit" / "tables"
    audit.mkdir(parents=True)
    with (audit / "c1_vs_c2_cep_gates.csv").open("w", encoding="utf-8") as handle:
        handle.write("xi,bracket\n" + "\n".join(f"{index},0.2" for index in range(17)))
    cep = module.build_scope_plan("cep", audit.parent, tmp_path / "cep", module.CALCULATION_SHA, "123")
    assert cep["frozen_bracket_count"] == 17
    assert cep["max_new_temperature_slices"] == 18
    cross = module.build_scope_plan("crossover", audit.parent, tmp_path / "cross", module.CALCULATION_SHA, "123")
    assert cross["diagnostic_levels"]["4"]["n_scan"] == 80
    assert cross["gate"]["response_relative"] == 0.025


def test_workflow_is_density_first_and_cost_stopped():
    text = WORKFLOW.read_text(encoding="utf-8")
    assert "scope=density|cep|crossover" not in text  # choice inputs are explicit below
    assert "options: [density, cep, crossover]" in text
    assert "options: [numerical, aggregate_replay]" in text
    assert "issue130_c2_density_feasibility_20260812" in text
    assert "source_run_id" in text and "source_calculation_sha" in text
    assert "failed_job_keys" in text and "rerun_failed_only" in text
    assert "density_numerical" in text and "max-parallel: 10" in text
    assert "C2 limited density xi=" in text
    assert "150" in text
    assert "scope_plan" in text and "crossover" in text
    assert "recover_source" in text
    assert "recovered-source-artifacts" in text


def test_workflow_yaml_has_density_and_deferred_scope_contracts():
    workflow = yaml.safe_load(WORKFLOW.read_text(encoding="utf-8"))
    assert "jobs" in workflow
    assert {"density_numerical", "density_aggregate", "scope_plan"}.issubset(workflow["jobs"])
    assert "scope" in workflow[True]["workflow_dispatch"]["inputs"]
    assert workflow[True]["workflow_dispatch"]["inputs"]["scope"]["default"] == "density"


def test_density_aggregate_pins_plotting_dependency():
    workflow = yaml.safe_load(WORKFLOW.read_text(encoding="utf-8"))
    steps = workflow["jobs"]["density_aggregate"]["steps"]
    setup = next(step for step in steps if step.get("name") == "Setup Python plotting environment")
    install = next(step for step in steps if step.get("name") == "Install fixed plotting dependency")
    assert setup["uses"] == "actions/setup-python@v5"
    assert setup["with"]["python-version"] == "3.11"
    assert install["run"] == "python -m pip install matplotlib==3.9.2"


def test_aggregate_replay_uses_explicit_repository_for_gh_calls():
    text = WORKFLOW.read_text(encoding="utf-8")
    assert 'gh run view "$SOURCE_RUN_ID" --repo "$GITHUB_REPOSITORY"' in text
    assert 'gh run download "$SOURCE_RUN_ID" --repo "$GITHUB_REPOSITORY"' in text


def test_recovery_overlay_restores_rows_from_production_eval(tmp_path):
    module = load_module()
    source = tmp_path / "source"
    artifact = source / "c2-limited-density-xi--0.35"
    eval_dir = artifact / "anchors" / "T_51p0" / "production_eval"
    eval_dir.mkdir(parents=True)
    pool = artifact / "fine_pool.csv"
    fieldnames = ["xi", "T_MeV", "rho", "muq_MeV", "residual_norm", "iterations",
                  "converged", "finite", "sampling_role", "rho_level", "calculation_sha"]
    with pool.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for rho_index in range(100, module.RHO_COUNT):
            rho = rho_index * module.RHO_FINE_STEP
            writer.writerow({"xi": "-0.35", "T_MeV": "51.0", "rho": str(rho),
                             "muq_MeV": "300.0", "residual_norm": "0.0", "iterations": "2",
                             "converged": "true", "finite": "true",
                             "sampling_role": "uniform_nested_fine_pool", "rho_level": "0",
                             "calculation_sha": module.CALCULATION_SHA})
    eval_path = eval_dir / "prod_eval_T51p000000_memoized.csv"
    eval_fields = ["T_MeV", "rho", "xi", "mu_avg_MeV", "residual_norm", "iterations", "converged"]
    with eval_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=eval_fields)
        writer.writeheader()
        for rho_index in range(module.RHO_COUNT):
                writer.writerow({"T_MeV": "51.0", "rho": str(rho_index * module.RHO_FINE_STEP),
                             "xi": "-0.35", "mu_avg_MeV": "300.0", "residual_norm": "0.0",
                             "iterations": "2", "converged": "true"})
    slices = artifact / "slice_metrics.csv"
    slices.write_text("failed_points\n0\n", encoding="utf-8")
    costs = artifact / "method_costs.csv"
    costs.write_text("runner_seconds\n1\n", encoding="utf-8")
    manifest = {
        "schema_version": module.JOB_SCHEMA, "scope": "density", "xi": -0.35,
        "source_run_id": "123", "calculation_sha": module.CALCULATION_SHA,
        "postprocess_sha": "postprocess", "solver_called": True,
        "files": {name: module.sha256(artifact / name)
                  for name in ("fine_pool.csv", "slice_metrics.csv", "method_costs.csv")},
    }
    (artifact / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    output = tmp_path / "recovered"
    # Build a one-shard copy to test the actual helper without requiring ten
    # full fixtures in this focused unit test.
    shutil.copytree(source, output)
    recovered = module._recover_manifest(
        output / artifact.name / "manifest.json", module.CALCULATION_SHA,
        "postprocess", "123", "recovery-head",
    )
    assert recovered["recovered_rows"] == 100
    assert recovered["final_rows"] == module.RHO_COUNT
    overlay_manifest = json.loads((output / artifact.name / "manifest.json").read_text(encoding="utf-8"))
    assert overlay_manifest["recovery"]["method"] == module.RECOVERY_METHOD
    assert overlay_manifest["recovery"]["solver_called"] is False
    assert overlay_manifest["recovery"]["recovery_postprocess_sha"] == "recovery-head"


def test_aggregate_validator_rejects_solver_called(tmp_path):
    module = load_module()
    aggregate = tmp_path / "aggregate"
    aggregate.mkdir()
    manifest = {
        "schema_version": module.SCHEMA_VERSION,
        "source_calculation_sha": module.CALCULATION_SHA,
        "source_run_id": "123",
        "solver_called": True,
    }
    (aggregate / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    for name in ("route_comparison.csv", "component_geometry.csv", "selected_point_index.csv",
                 "candidate_point_index.csv", "cost_frontier.csv", "claim_ledger.csv"):
        (aggregate / name).write_text("route,cap\n", encoding="utf-8")
    try:
        module.validate_aggregate(aggregate, module.CALCULATION_SHA, "123")
    except ValueError as error:
        assert "solver_called" in str(error)
    else:
        raise AssertionError("solver_called=true must be rejected")


def test_plot_manifest_is_written_with_provenance(tmp_path):
    module = load_module()
    aggregate = tmp_path / "aggregate"
    aggregate.mkdir()
    (aggregate / "route_comparison.csv").write_text(
        "route,cap,xi,T_MeV,density_error\n"
        "stage_b_features_v1,12,0.0,51.0,0.0\n"
        "balanced_density_features_v2,12,0.0,51.0,0.0\n"
        "geometry_feedback_v2,12,0.0,51.0,0.0\n",
        encoding="utf-8",
    )
    for name in ("component_geometry.csv", "selected_point_index.csv", "cost_frontier.csv"):
        (aggregate / name).write_text("field\nvalue\n", encoding="utf-8")
    figures = module.render_representative_plot(aggregate)
    plot_manifest = json.loads((aggregate / "plot_manifest.json").read_text(encoding="utf-8"))
    assert plot_manifest["schema_version"] == module.PLOT_SCHEMA_VERSION
    assert plot_manifest["generator_sha256"]
    assert isinstance(figures, list)
    for figure in figures:
        assert (aggregate / figure["path"]).is_file()
