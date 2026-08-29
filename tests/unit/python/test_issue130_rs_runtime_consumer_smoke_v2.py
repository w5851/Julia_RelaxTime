from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "pnjl" / "build_issue130_rs_runtime_consumer_smoke_v2.py"
SPEC = importlib.util.spec_from_file_location("issue130_rs_consumer_smoke_v2", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _point(mode: str, *, phase_T: float = 198.2) -> dict[str, object]:
    return {
        "mode": mode,
        "source_kind": "legacy" if mode == "legacy" else "candidate",
        "boundary_rows": 5 if mode == "legacy" else 49,
        "boundary_xi_used": 0.0,
        "crossover_xi_used": 0.0,
        "first_order_xi_used": 0.0,
        "mode_a_phase_kind": "crossover",
        "mode_a_phase_T_MeV": phase_T,
        "crossover_T_MeV": phase_T,
        "first_order_T_MeV": 130.0 if mode == "legacy" else 131.0,
        "solver_called": False,
    }


def test_v2_constants_pin_merge_and_solver_free_contract() -> None:
    assert MODULE.MERGE_SHA == "3b19246fb911be4a2efa75fbe14fcb9a793ca256"
    assert len(MODULE.CALCULATION_SHA) == 40
    assert MODULE.SCHEMA.endswith("_v2")


def test_parity_rows_keep_candidate_legacy_numeric_delta() -> None:
    payload = {"consumer_points": [_point("runtime"), _point("diagnostic", phase_T=198.18), _point("legacy", phase_T=198.1)]}
    rows = MODULE._parity_rows(payload)
    phase_row = next(row for row in rows if row["field"] == "mode_a_phase_T_MeV")
    assert phase_row["runtime_minus_legacy"] == pytest.approx(0.1)
    assert phase_row["verdict"] == "mode_sensitive"


def test_source_smoke_rejects_solver_invocation() -> None:
    payload = {
        "schema_version": "pnjl_issue130_rs_runtime_consumer_smoke_v1",
        "repo_head": MODULE.MERGE_SHA,
        "solver_called": True,
        "sources": {},
    }
    with pytest.raises(ValueError, match="solver_called"):
        MODULE._validate_source_smoke(payload, "candidate-manifest")


def test_source_smoke_accepts_author_accepted_fallback_view() -> None:
    diagnostics = {
        "runtime_view": "certified_candidate_with_accepted_then_legacy_fallback",
        "fallback_enabled": True,
        "fallback_order": "strict_candidate>accepted_downstream>legacy_snapshot",
        "candidate_manifest_sha256": "candidate-manifest",
        "accepted_manifest_sha256": "a" * 64,
        "accepted_layer_manifest_sha256": "b" * 64,
        "accepted_fallback_row_counts": {"boundary": 2},
    }
    payload = {
        "schema_version": "pnjl_issue130_rs_runtime_consumer_smoke_v1",
        "repo_head": MODULE.MERGE_SHA,
        "solver_called": False,
        "sources": {
            "runtime": {"source_kind": "candidate", "runtime_enabled": True, "diagnostics": diagnostics},
            "diagnostic": {"source_kind": "candidate", "runtime_enabled": False, "diagnostics": {"runtime_view": "diagnostic_all_rows"}},
            "legacy": {"source_kind": "legacy", "runtime_enabled": True, "diagnostics": {"runtime_view": "legacy"}},
        },
        "consumer_points": [_point("runtime"), _point("diagnostic"), _point("legacy")],
    }
    MODULE._validate_source_smoke(payload, "candidate-manifest")
