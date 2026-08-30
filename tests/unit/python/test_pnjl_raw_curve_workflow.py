from pathlib import Path

import yaml


REPO_ROOT = Path(__file__).resolve().parents[3]
WORKFLOW = REPO_ROOT / ".github" / "workflows" / "pnjl-raw-curve-archive.yml"
SHARD_WORKFLOW = REPO_ROOT / ".github" / "workflows" / "pnjl-raw-curve-archive-shard.yml"


def test_raw_only_workflow_is_independent_and_strict():
    text = WORKFLOW.read_text(encoding="utf-8")
    shard_text = SHARD_WORKFLOW.read_text(encoding="utf-8")
    payload = yaml.load(text, Loader=yaml.BaseLoader)
    shard_payload = yaml.load(shard_text, Loader=yaml.BaseLoader)
    assert payload["jobs"]["raw_shards"]["uses"] == "./.github/workflows/pnjl-raw-curve-archive-shard.yml"
    assert shard_payload["jobs"]["raw-curve-shard"]["runs-on"] == "ubuntu-latest"
    assert "phase_grid_convergence" in text
    assert "axis=rho" in text or "build-coverage" in text
    assert "no Cartesian completion" in text or "build-coverage" in text
    assert "--require-full-domain" in text
    assert "pnjl-raw-curve-archive-" in text
    assert "independent_oracle" in (REPO_ROOT / "scripts" / "analysis" / "pnjl_raw_curve_archive.py").read_text(encoding="utf-8")


def test_raw_runner_is_equilibrium_only():
    runner = (REPO_ROOT / "scripts" / "analysis" / "pnjl_raw_curve_archive_oracle_job.jl").read_text(encoding="utf-8")
    assert "new_rho_point_session" in runner
    assert "rho_point!" in runner
    assert "0.003125" in runner
    assert "rs_reduced_adaptive" in runner
    assert "p_num=24" in runner
    assert "iterations=80" in runner
    assert "run_production_phase_pipeline" not in runner


def test_existing_dense_workflows_are_not_raw_retention_entrypoints():
    for relative in (
        ".github/workflows/pnjl-dense-reference.yml",
        ".github/workflows/pnjl-dense-reference-shard.yml",
        ".github/workflows/pnjl-dense-reference-resume.yml",
    ):
        text = (REPO_ROOT / relative).read_text(encoding="utf-8")
        assert "retain_raw_curves" not in text
        assert "raw_curve_method" not in text
