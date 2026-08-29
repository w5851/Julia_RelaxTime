from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "pnjl" / "audit_issue130_phase_reference_legacy_retirement.py"
SPEC = importlib.util.spec_from_file_location("issue130_phase_reference_legacy_retirement", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_semantic_keys_are_decimal_normalized_and_reject_nonfinite() -> None:
    assert MODULE.key_tuple("boundary", {"xi": "0.10000000000000001", "T_MeV": "50.0"}) == (
        0.1,
        50.0,
    )
    with pytest.raises(ValueError):
        MODULE.semantic_float("nan")
    with pytest.raises(ValueError):
        MODULE.semantic_float("inf")


def test_coverage_uses_semantic_keys_and_keeps_uncertified_candidate_out() -> None:
    candidate = {
        "boundary": [
            {"xi": 0.0, "T_MeV": 50.0, "certified": True},
            {"xi": 0.0, "T_MeV": 60.0, "certified": False},
        ],
    }
    legacy = {
        "boundary": [
            {"xi": "0.0", "T_MeV": "50.00000000000001"},
            {"xi": "0.0", "T_MeV": "60.0"},
        ],
    }
    coverage, fallback = MODULE.coverage_matrix(candidate, legacy)
    row = next(item for item in coverage if item["table"] == "boundary")
    assert row["legacy_certified_overlap"] == 1
    assert row["legacy_fallback_required"] == 1
    assert fallback == [
        {
            "table": "boundary",
            "legacy_row": 3,
            "key": '{"T_MeV":60.0,"xi":0.0}',
            "reason": "candidate_key_absent_or_uncertified",
            "disposition": "legacy_fallback_required",
        }
    ]


def test_active_consumer_and_historical_classification() -> None:
    assert MODULE.classify_path("scripts/relaxtime/run_gap_transport_scan.jl") == "active_repository_contract"
    assert MODULE.classify_path("docs/analysis/pnjl/example/README.md") == "historical_evidence"
    assert MODULE.occurrence_role(
        "scripts/relaxtime/phase_reference_adapter.jl", "load_default_phase_reference_runtime"
    ) == "fallback_or_rollback"
    assert MODULE.occurrence_role("docs/analysis/pnjl/example/README.md", "legacy reference") == "legacy_path_contract"


def _candidate_info() -> dict[str, object]:
    return {
        "root_manifest": {
            "calculation_sha": MODULE.CALCULATION_SHA,
            "runtime_consumption": False,
            "source_run_id": "1",
            "replay_run_id": "2",
        },
        "strict_manifest": {
            "calculation_sha": MODULE.CALCULATION_SHA,
            "solver_called": False,
            "reference_write": False,
        },
    }


def _coverage(*, fallback: int, uncertified: int = 0) -> list[dict[str, object]]:
    return [
        {
            "table": "boundary",
            "candidate_rows": 1,
            "candidate_certified_rows": 1,
            "candidate_uncertified_rows": uncertified,
            "candidate_duplicate_keys": 0,
            "legacy_rows": fallback,
            "legacy_duplicate_keys": 0,
            "legacy_certified_overlap": 0,
            "legacy_fallback_required": fallback,
        }
    ]


def test_decision_precedence_retains_snapshot_until_coverage_and_consumers_are_clear() -> None:
    legacy = {"integrity_pass": True}
    inconclusive = MODULE.compute_decision(_candidate_info(), legacy, _coverage(fallback=1), [])
    assert inconclusive["verdict"] == "retirement_inconclusive"
    assert inconclusive["physical_deletion_eligible"] is False
    assert "candidate_or_legacy_key_coverage_incomplete" in inconclusive["stop_reasons"]

    eligible = MODULE.compute_decision(_candidate_info(), legacy, _coverage(fallback=0), [])
    assert eligible["verdict"] == "legacy_physical_deletion_candidate"
    assert eligible["physical_deletion_eligible"] is True

    # An uncertified candidate-only row is reported, but it is not a legacy
    # dependency unless the same semantic key is needed from the snapshot.
    candidate_diag = MODULE.compute_decision(
        _candidate_info(), legacy, _coverage(fallback=0, uncertified=1), []
    )
    assert candidate_diag["coverage_complete"] is True
    assert candidate_diag["physical_deletion_eligible"] is True

    blocker = MODULE.compute_decision(
        _candidate_info(),
        legacy,
        _coverage(fallback=0),
        [{"retirement_blocker": True, "kind": "runtime_adapter"}],
    )
    assert blocker["verdict"] == "retirement_inconclusive"
    assert blocker["physical_deletion_eligible"] is False


def test_candidate_provenance_invalid_stops_before_deletion() -> None:
    info = _candidate_info()
    info["strict_manifest"] = {
        "calculation_sha": MODULE.CALCULATION_SHA,
        "solver_called": True,
        "reference_write": False,
    }
    decision = MODULE.compute_decision(info, {"integrity_pass": True}, _coverage(fallback=0), [])
    assert decision["verdict"] == "legacy_audit_input_invalid"
    assert decision["physical_deletion_eligible"] is False
