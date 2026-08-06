import importlib.util
import json
from pathlib import Path

import pytest
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_PATH = PROJECT_ROOT / "scripts" / "pnjl" / "plan_dense_reference_action.py"
WORKFLOW_PATH = PROJECT_ROOT / ".github" / "workflows" / "pnjl-dense-reference.yml"
ASSESS_WORKFLOW_PATH = PROJECT_ROOT / ".github" / "workflows" / "pnjl-dense-reference-assess-xi.yml"
SHARD_WORKFLOW_PATH = PROJECT_ROOT / ".github" / "workflows" / "pnjl-dense-reference-shard.yml"
REPLAY_WORKFLOW_PATH = PROJECT_ROOT / ".github" / "workflows" / "pnjl-dense-reference-replay.yml"
RESUME_WORKFLOW_PATH = PROJECT_ROOT / ".github" / "workflows" / "pnjl-dense-reference-resume.yml"


def load_module():
    spec = importlib.util.spec_from_file_location("dense_action_plan", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def full_reference_inputs(**updates):
    payload = {
        "xi_values": "",
        "xi_min": "-0.5",
        "xi_max": "0.5",
        "xi_step": "0.05",
        "adaptive_xi": True,
        "crossover_only": False,
        "adaptive_T": True,
        "rho_geometry_convergence": True,
        "overwrite": True,
        "advanced_config_json": "{}",
    }
    payload.update(updates)
    return payload


def test_plans_one_xi_initial_shards_and_defers_recursive_levels():
    module = load_module()
    plan = module.build_action_plan(full_reference_inputs())
    entries = plan["initial_matrix"]["include"]

    assert plan["xi_refine_levels"] == 2
    assert len(entries) == 41
    assert len({entry["xi"] for entry in entries}) == 41
    assert all("," not in entry["xi"] for entry in entries)
    assert [entry["stage"] for entry in entries].count("anchor") == 21
    assert [entry["stage"] for entry in entries].count("level1") == 20
    assert len(plan["initial_intervals"]) == 20
    assert "--adaptive-xi" in plan["common_args"]


def test_explicit_endpoint_is_preserved_and_advanced_controls_win():
    module = load_module()
    plan = module.build_action_plan(
        full_reference_inputs(
            xi_min="-0.5",
            xi_max="0.52",
            xi_step="0.2",
            advanced_config_json=json.dumps({"xi_refine_levels": 3, "xi_position_tol": 0.025}),
        )
    )

    assert plan["expected_xi_csv"].endswith("0.52")
    assert plan["xi_refine_levels"] == 3
    assert plan["xi_position_tol"] == "0.025"
    position_indices = [i for i, token in enumerate(plan["common_args"]) if token == "--xi-position-tol"]
    assert plan["common_args"][position_indices[-1] + 1] == "0.025000000000000001"


@pytest.mark.parametrize(
    "updates, message",
    [
        ({"crossover_only": True}, "adaptive_xi requires crossover_only=false"),
        ({"advanced_config_json": '{"xi_refine_levels":4}'}, "supports at most 3 levels"),
    ],
)
def test_rejects_unstageable_adaptive_plans(updates, message):
    module = load_module()
    with pytest.raises(ValueError, match=message):
        module.build_action_plan(full_reference_inputs(**updates))


@pytest.mark.parametrize(
    "calculation_ref",
    ["main", "fc5d50a5", "G" * 40, "FC5D50A5D22FE79F7CC52EF2ECBDB79376E256E1"],
)
def test_rejects_mutable_or_noncanonical_calculation_ref(calculation_ref):
    module = load_module()
    with pytest.raises(ValueError, match="immutable lowercase 40-character Git SHA"):
        module.build_action_plan(full_reference_inputs(calculation_ref=calculation_ref))


def test_accepts_full_calculation_sha():
    module = load_module()
    plan = module.build_action_plan(
        full_reference_inputs(calculation_ref="fc5d50a5d22fe79f7cc52ef2ecbdb79376e256e1")
    )
    assert plan["initial_shard_count"] == 41


def test_hybrid_policy_is_planned_with_four_rho_levels():
    module = load_module()
    plan = module.build_action_plan(
        full_reference_inputs(
            advanced_config_json=json.dumps({
                "rho_refinement_policy": "rho_support_hybrid",
                "rho_refine_levels": 4,
                "rho_support_fine_step": 0.025,
                "rho_support_targeted_cap": 12,
            })
        )
    )
    assert "rho_support_hybrid" in plan["common_args"]
    assert plan["common_args"][plan["common_args"].index("--rho-refine-levels") + 1] == "4"
    assert "--rho-hybrid-endpoint-policy" in plan["common_args"]
    endpoint_index = plan["common_args"].index("--rho-hybrid-endpoint-policy")
    assert plan["common_args"][endpoint_index + 1] == "three_crossing_endpoint_local_v2"


def test_hybrid_endpoint_policy_can_be_explicitly_overridden_for_replay():
    module = load_module()
    plan = module.build_action_plan(
        full_reference_inputs(
            advanced_config_json=json.dumps({
                "rho_refinement_policy": "rho_support_hybrid",
                "rho_refine_levels": 4,
                "rho_hybrid_endpoint_policy": "bounded_zero_density_v1",
            })
        )
    )
    endpoint_index = plan["common_args"].index("--rho-hybrid-endpoint-policy")
    assert plan["common_args"][endpoint_index + 1] == "bounded_zero_density_v1"


def test_hybrid_policy_rejects_wrong_refinement_level():
    module = load_module()
    with pytest.raises(ValueError, match="rho_support_hybrid requires rho_refine_levels=4"):
        module.build_action_plan(
            full_reference_inputs(
                advanced_config_json='{"rho_refinement_policy":"rho_support_hybrid","rho_refine_levels":3}'
            )
        )


def test_workflow_uses_reusable_one_xi_and_assessment_jobs():
    workflow_text = WORKFLOW_PATH.read_text(encoding="utf-8")
    assess_workflow_text = ASSESS_WORKFLOW_PATH.read_text(encoding="utf-8")
    workflow = yaml.load(workflow_text, Loader=yaml.BaseLoader)
    jobs = workflow["jobs"]
    assert jobs["generate_initial_shards"]["uses"] == "./.github/workflows/pnjl-dense-reference-shard.yml"
    assert jobs["assess_xi_level1"]["uses"] == "./.github/workflows/pnjl-dense-reference-assess-xi.yml"
    assert "xi_list" not in jobs["generate_initial_shards"]["with"]
    assert "os.environ['ADVANCED_CONFIG_JSON']" in workflow_text
    assert "advanced_config_json: `${{ inputs" not in workflow_text
    assert '"--expected-xi-list=${{ needs.plan.outputs.expected_xi_csv }}"' in workflow_text
    assert '"--expected-xi-list=${{ inputs.expected_xi_csv }}"' in assess_workflow_text
    for job_name in ("generate_initial_shards", "generate_xi_level2", "generate_xi_level3"):
        assert jobs[job_name]["with"]["calculation_ref"] == "${{ inputs.calculation_ref }}"
    for job_name in ("assess_xi_level1", "assess_xi_level2", "assess_xi_level3"):
        assert jobs[job_name]["with"]["source_calculation_sha"] == (
            "${{ inputs.calculation_ref || github.sha }}"
        )
    assert (
        '--expected-calculation-git-commit "${{ inputs.calculation_ref || github.sha }}"'
        in workflow_text
    )
    assert '--postprocess-git-commit "${GITHUB_SHA}"' in workflow_text
    assert "actions: read" in workflow_text
    assert "GH_TOKEN: ${{ github.token }}" in workflow_text
    assert 'gh run download "${{ github.run_id }}"' in workflow_text
    assert "GH_TOKEN: ${{ github.token }}" in assess_workflow_text
    assert 'gh run download "${{ github.run_id }}"' in assess_workflow_text
    assert "next_xi_count" in assess_workflow_text


def test_resume_workflow_reuses_source_shards_and_keeps_calculation_sha():
    workflow_text = RESUME_WORKFLOW_PATH.read_text(encoding="utf-8")
    shard_workflow_text = SHARD_WORKFLOW_PATH.read_text(encoding="utf-8")
    workflow = yaml.load(workflow_text, Loader=yaml.BaseLoader)
    inputs = workflow["on"]["workflow_dispatch"]["inputs"]
    assert set(inputs) == {"source_run_id", "source_calculation_sha", "tag", "expected_xi_list"}
    assert 'gh run download "${{ inputs.source_run_id }}"' in workflow_text
    assert "actions/download-artifact@v4" not in workflow_text
    assert "calculation_ref: ${{ inputs.source_calculation_sha }}" in workflow_text
    assert "source and resumed shards" in workflow_text
    assert "diagnostic candidate only; no reference promotion" in workflow_text
    assert "include_current_run_shards: ${{ fromJSON(needs.assess_xi_level2.outputs.next_xi_count) > 0 }}" in workflow_text
    assert "runner.temp" not in workflow_text
    assert "calculation_ref:" in shard_workflow_text
    assert "ref: ${{ inputs.calculation_ref || github.sha }}" in shard_workflow_text


def test_postprocess_replay_uses_cross_run_artifacts_and_dual_provenance():
    workflow_text = REPLAY_WORKFLOW_PATH.read_text(encoding="utf-8")
    workflow = yaml.load(workflow_text, Loader=yaml.BaseLoader)
    inputs = workflow["on"]["workflow_dispatch"]["inputs"]
    assert {"source_run_id", "source_calculation_sha", "tag", "expected_xi_list"} <= set(inputs)
    assert {
        "final_assessment_level",
        "final_assessment_intervals_json",
        "final_assessment_position_tol",
        "final_assessment_density_tol",
        "final_assessment_maxwell_area_tol",
        "final_assessment_response_rtol",
    } <= set(inputs)
    assert 'gh run download "${{ inputs.source_run_id }}"' in workflow_text
    assert 'gh run download "${{ github.run_id }}"' in workflow_text
    assert "GH_TOKEN: ${{ github.token }}" in workflow_text
    assert "actions/download-artifact@v4" not in workflow_text
    jobs = workflow["jobs"]
    assert jobs["assess_final"]["uses"] == "./.github/workflows/pnjl-dense-reference-assess-xi.yml"
    assert jobs["assess_final"]["with"]["source_calculation_sha"] == "${{ inputs.source_calculation_sha }}"
    assert jobs["assess_final"]["with"]["include_current_run_shards"] == "false"
    assert "needs.assess_final.result == 'success'" in workflow_text
    assert "replayed final xi convergence record" in workflow_text
    assert '--expected-calculation-git-commit "${{ inputs.source_calculation_sha }}"' in workflow_text
    assert '--postprocess-git-commit "${GITHUB_SHA}"' in workflow_text
    assert "diagnostic replay only; no reference promotion" in workflow_text


def test_dense_reference_cross_job_downloads_use_paginated_gh_cli():
    for path in (
        WORKFLOW_PATH,
        ASSESS_WORKFLOW_PATH,
        REPLAY_WORKFLOW_PATH,
        RESUME_WORKFLOW_PATH,
    ):
        workflow_text = path.read_text(encoding="utf-8")
        assert "actions/download-artifact@v4" not in workflow_text
        assert "gh run download" in workflow_text
        assert "GH_TOKEN: ${{ github.token }}" in workflow_text
        assert "GH_REPO: ${{ github.repository }}" in workflow_text
