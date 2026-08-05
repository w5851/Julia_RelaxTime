import csv
import hashlib
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "build_pnjl_maxwell_endpoint_limit_contract.py"
EVIDENCE = ROOT / "docs" / "analysis" / "pnjl_maxwell_endpoint_limit_contract_v1"


def load_contract_module():
    spec = importlib.util.spec_from_file_location("endpoint_limit_contract", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


CONTRACT = load_contract_module()


def valid_inputs():
    metrics = []
    for level in range(CONTRACT.EXPECTED_TARGETED_CAP + 1):
        upper = CONTRACT.EXPECTED_BASE_STEP / (2**level)
        metrics.append(
            {
                "level": str(level),
                "targeted_additions": str(level),
                "status": "first_order",
                "reason": "unique_three_crossing_candidate",
                "candidate_count": "1",
                "candidate_mu": str(331.5 + level / 1000),
                "candidate_area": "1e-7",
                "rho_hadron": str(upper / 2),
                "rho_quark": str(2.8 + level / 10000),
                "endpoint_dependent": "true",
                "bracket_low": "0.0",
                "bracket_high": str(upper),
                "bracket_width": str(upper),
                "positive_rho_bracket": "false",
                "position_error_MeV": str(0.01 / (level + 1)),
                "density_error": str(0.001 / (level + 1)),
                "geometry_converged": "true" if level else "false",
                "unique_solves": str(1281 + level),
                "solver_failures": "0",
            }
        )

    curves = []
    for index in range(1281 + CONTRACT.EXPECTED_TARGETED_CAP):
        curves.append(
            {
                "xi": str(CONTRACT.EXPECTED_XI),
                "T_MeV": str(CONTRACT.EXPECTED_T_MEV),
                "rho": str(index / 10000),
                "muq_MeV": str(330 + index / 10000),
                "pressure_fm4": "1.0",
                "residual_norm": "1e-10",
                "converged": "true",
                "finite": "true",
                "calculation_sha": CONTRACT.SOURCE_CALCULATION_SHA,
                "workflow_head_sha": CONTRACT.SOURCE_WORKFLOW_HEAD_SHA,
            }
        )

    costs = [
        {
            "targeted_additions": str(CONTRACT.EXPECTED_TARGETED_CAP),
            "unique_solves": str(len(curves)),
            "point_requests": str(len(curves)),
            "fixedrho_requests": str(len(curves)),
            "cache_hits": "0",
            "nonconverged_attempts": "0",
            "fallback_count": "0",
            "governed_rescue_count": "0",
            "retry_count": "0",
            "exception_count": "0",
        }
    ]
    source_manifest = {
        "schema_version": CONTRACT.SOURCE_SCHEMA,
        "verdict": CONTRACT.SOURCE_VERDICT,
        "calculation_sha": CONTRACT.SOURCE_CALCULATION_SHA,
        "workflow_head_sha": CONTRACT.SOURCE_WORKFLOW_HEAD_SHA,
        "reference_write": False,
        "base_rho_step": CONTRACT.EXPECTED_BASE_STEP,
    }
    return metrics, costs, curves, source_manifest


def evaluate(mutator=None):
    metrics, costs, curves, source_manifest = valid_inputs()
    if mutator is not None:
        mutator(metrics, costs, curves, source_manifest)
    return CONTRACT.evaluate_endpoint_limit_contract(
        metrics, costs, curves, source_manifest
    )


def test_endpoint_limit_contract_accepts_converged_zero_bound():
    certificate, errors = evaluate()
    assert errors == []
    assert certificate["contract_pass"] is True
    assert certificate["physical_classification"] == "confirmed_first_order"
    assert certificate["certificate_kind"] == "endpoint_limited_first_order"
    assert certificate["rho_hadron_representative"] == 0.0
    assert certificate["rho_hadron_upper_bound"] == 0.003125 / (2**12)
    assert certificate["strict_positive_rho_claimed"] is False


def test_multiple_candidate_is_inconclusive():
    def mutate(metrics, *_):
        metrics[4]["candidate_count"] = "2"

    certificate, errors = evaluate(mutate)
    assert certificate["contract_pass"] is False
    assert certificate["unique_candidate_gate"] is False
    assert "unique three-crossing candidate gate failed" in errors


def test_non_bisecting_endpoint_trace_is_inconclusive():
    def mutate(metrics, *_):
        metrics[6]["bracket_high"] = metrics[5]["bracket_high"]
        metrics[6]["bracket_width"] = metrics[5]["bracket_width"]

    certificate, errors = evaluate(mutate)
    assert certificate["contract_pass"] is False
    assert certificate["bracket_halving_gate"] is False
    assert "endpoint bracket does not shrink deterministically by bisection" in errors


def test_area_or_solver_failure_cannot_produce_certificate():
    def mutate(metrics, costs, *_):
        metrics[8]["candidate_area"] = "6e-6"
        metrics[8]["solver_failures"] = "1"
        costs[0]["nonconverged_attempts"] = "1"

    certificate, errors = evaluate(mutate)
    assert certificate["contract_pass"] is False
    assert certificate["area_solver_gate"] is False
    assert certificate["trace_accounting_gate"] is False
    assert certificate["cost_gate"] is False
    assert "candidate area exceeds the internal solver tolerance" in errors


def test_committed_evidence_schema_provenance_and_finite_tables():
    manifest = json.loads((EVIDENCE / "manifest.json").read_text(encoding="utf-8"))
    provenance = json.loads((EVIDENCE / "provenance.json").read_text(encoding="utf-8"))
    policy = json.loads((EVIDENCE / "selected_policy.json").read_text(encoding="utf-8"))

    assert manifest["schema_version"] == CONTRACT.SCHEMA_VERSION
    assert manifest["verdict"] == "endpoint_limited_first_order_candidate"
    assert manifest["physical_classification"] == "confirmed_first_order"
    assert manifest["certificate_kind"] == "endpoint_limited_first_order"
    assert manifest["contract_errors"] == []
    assert manifest["solver_called"] is False
    assert manifest["reference_write"] is False
    assert provenance["source_run_id"] == CONTRACT.SOURCE_RUN_ID
    assert provenance["source_calculation_sha"] == CONTRACT.SOURCE_CALCULATION_SHA
    assert provenance["source_workflow_head_sha"] == CONTRACT.SOURCE_WORKFLOW_HEAD_SHA
    assert provenance["source_job_manifest_sha256"] == CONTRACT.SOURCE_JOB_MANIFEST_SHA256
    assert (
        provenance["source_aggregate_manifest_sha256"]
        == CONTRACT.SOURCE_AGGREGATE_MANIFEST_SHA256
    )
    assert provenance["author_decision"]["status"] == "accepted"
    assert policy["rho_hadron_representative"] == 0.0
    assert policy["strict_positive_rho_claimed"] is False

    for relative_path, expected_hash in manifest["files"].items():
        payload = (EVIDENCE / relative_path).read_bytes()
        assert hashlib.sha256(payload).hexdigest() == expected_hash

    invalid = {"nan", "inf", "+inf", "-inf"}
    for path in (EVIDENCE / "tables").glob("*.csv"):
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.reader(handle):
                assert all(cell.strip().lower() not in invalid for cell in row)
