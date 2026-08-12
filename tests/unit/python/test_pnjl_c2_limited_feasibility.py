import importlib.util
import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "pnjl_c2_limited_feasibility.py"
WORKFLOW = ROOT / ".github" / "workflows" / "pnjl-c2-limited-feasibility.yml"


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
