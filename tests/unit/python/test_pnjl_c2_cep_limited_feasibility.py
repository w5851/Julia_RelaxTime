import csv
import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "pnjl_c2_cep_limited_feasibility.py"
WORKFLOW = ROOT / ".github" / "workflows" / "pnjl-c2-cep-limited-feasibility.yml"


def load_module():
    spec = importlib.util.spec_from_file_location("c2_cep_limited", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_v2_schema_and_workflow_contract():
    module = load_module()
    assert module.SCHEMA_VERSION == "pnjl_c2_cep_limited_feasibility_v2"
    assert module.JOB_SCHEMA == "pnjl_c2_cep_limited_feasibility_job_v2"
    text = WORKFLOW.read_text(encoding="utf-8")
    assert "pnjl_c2_cep_limited_feasibility_v2" in text
    assert "issue130_c2_cep_limited_feasibility_v2_20260812" in text
    assert "solver_called=false" not in text.lower()


def test_numerical_job_materializes_hashed_frozen_brackets():
    text = WORKFLOW.read_text(encoding="utf-8")
    relative_path = "docs/analysis/pnjl/c2_limited_feasibility_v1/cep_failures.csv"
    assert f"bracket_rel='{relative_path}'" in text
    assert "bracket_sha='f8775f3cb7e0457722a523b30eb89485239bdd1730e88f372eeae4e2bc21f3d2'" in text
    assert 'test -f "$bracket_rel"' in text
    assert 'cp "$bracket_rel" "calculation/$bracket_rel"' in text
    assert 'sha256sum "calculation/$bracket_rel"' in text


def test_hybrid_local_step_is_nested_below_fine_step():
    job = (ROOT / "scripts" / "analysis" / "pnjl_c2_cep_limited_feasibility_job.jl").read_text(
        encoding="utf-8"
    )
    assert "const RHO_FINE_STEP = 0.003125" in job
    assert "const HYBRID_LOCAL_STEP = RHO_FINE_STEP / 2" in job
    assert "local_step=HYBRID_LOCAL_STEP" in job
    assert "hybrid_local_step_contract" in job


def _write_costs(module, path: Path, *, valid: bool = True):
    path.write_text(
        "xi,method,unique_solves,point_requests,cache_hits,failed_points,runner_seconds\n"
        f"0.0,hybrid,10,12,2,0,1.0\n"
        f"0.0,oracle,20,20,0,{0 if valid else 1},2.0\n",
        encoding="utf-8",
    )


def test_cost_contract_requires_two_methods_and_reconciliation(tmp_path):
    module = load_module()
    costs = tmp_path / "method_costs.csv"
    _write_costs(module, costs)
    rows, seconds, unique = module._validate_costs(costs, 0.0)
    assert {row["method"] for row in rows} == {"hybrid", "oracle"}
    assert seconds == 3.0
    assert unique == 30

    costs.write_text(
        "xi,method,unique_solves,point_requests,cache_hits,failed_points,runner_seconds\n"
        "0.0,hybrid,10,11,2,0,1.0\n"
        "0.0,oracle,20,20,0,0,2.0\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="reconciliation"):
        module._validate_costs(costs, 0.0)


def test_source_validator_rejects_oracle_routing_leak(tmp_path):
    module = load_module()
    artifact = tmp_path / "xi_0"
    artifact.mkdir()
    pool = artifact / "fine_pool.csv"
    pool.write_text(
        "xi,T_MeV,rho,muq_MeV,residual_norm,converged,finite\n"
        "0.0,1.0,0.0,1.0,0.0,true,true\n",
        encoding="utf-8",
    )
    (artifact / "slice_metrics.csv").write_text(
        "T_MeV,hybrid_status,oracle_labels_used_for_routing\n1.0,ambiguous_near_critical,true\n",
        encoding="utf-8",
    )
    costs = artifact / "method_costs.csv"
    _write_costs(module, costs)
    manifest = {
        "schema_version": module.JOB_SCHEMA,
        "scope": "cep",
        "xi": 0.0,
        "source_run_id": "123",
        "calculation_sha": module.CALCULATION_SHA,
        "postprocess_sha": "a" * 40,
        "solver_called": True,
        "provenance": {
            "reference_write": False,
            "route_before_oracle_gate": True,
            "oracle_labels_used_for_routing": True,
        },
        "files": {name: module.sha256(artifact / name)
                  for name in ("fine_pool.csv", "slice_metrics.csv", "method_costs.csv")},
    }
    (artifact / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="oracle routing leakage"):
        module._validate_route_provenance(manifest)


def test_aggregate_validator_requires_cost_reconciliation(tmp_path):
    module = load_module()
    aggregate = tmp_path / "aggregate"
    aggregate.mkdir()
    (aggregate / "manifest.json").write_text(json.dumps({
        "schema_version": module.SCHEMA_VERSION,
        "source_calculation_sha": module.CALCULATION_SHA,
        "source_run_id": "123",
        "solver_called": False,
        "oracle_labels_used_for_routing": False,
        "verdict": "cep_ambiguity_inconclusive",
    }), encoding="utf-8")
    rows = "xi,T_last_first_order_MeV,T_first_monotone_MeV,bracket_width_T_MeV\n" + \
        "\n".join(f"{xi},1,2,1" for xi in module.TARGET_XI) + "\n"
    (aggregate / "cep_bracket_results.csv").write_text(rows, encoding="utf-8")
    for name in ("temperature_states.csv", "method_costs.csv", "route_trace.csv", "claim_ledger.csv"):
        (aggregate / name).write_text("field\nvalue\n", encoding="utf-8")
    (aggregate / "cost_frontier.csv").write_text(
        "runner_minutes,unique_solves,point_requests,cache_hits,point_request_reconciliation\n"
        "1,10,11,0,false\n", encoding="utf-8")
    manifest = json.loads((aggregate / "manifest.json").read_text(encoding="utf-8"))
    manifest["files"] = {
        name: module.sha256(aggregate / name)
        for name in ("cep_bracket_results.csv", "temperature_states.csv", "method_costs.csv",
                     "cost_frontier.csv", "route_trace.csv", "claim_ledger.csv")
    }
    (aggregate / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="reconciliation"):
        module.validate_aggregate(aggregate, module.CALCULATION_SHA, "123")
