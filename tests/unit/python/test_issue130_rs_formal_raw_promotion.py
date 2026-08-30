"""Contract checks for the Issue #130 RS raw promotion metadata."""

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
CASE = "first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2"
CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"


def test_registry_promotes_both_modes_as_raw_only():
    registry = json.loads(
        (ROOT / "data/outputs/results/relaxtime/transport/phase_guided/production_registry.json").read_text(
            encoding="utf-8"
        )
    )
    entries = [entry for entry in registry["entries"] if entry["case_slug"] == CASE]
    assert {entry["mode"] for entry in entries} == {
        "mode_a_fixed_muB_phase_scaled",
        "mode_b_fixed_T_sparse_muB",
    }
    assert all(entry["status"] == "approved" for entry in entries)
    assert all(entry["manuscript_eligible"] is False for entry in entries)
    assert all(entry["source_commit"] == CALCULATION_SHA for entry in entries)


def test_promotion_manifest_keeps_numerical_and_runtime_boundaries():
    manifest = json.loads(
        (ROOT / "docs/analysis/relaxtime/issue130_rs_transport_formal_raw_promotion_v1/manifest.json").read_text(
            encoding="utf-8"
        )
    )
    assert manifest["package_status"] == "approved_formal_raw"
    assert manifest["numerical_audit_status"] == "diagnostic_only"
    assert manifest["promotion_solver_called"] is False
    assert manifest["source_solver_called"] is True
    assert manifest["runtime_default_switch"] is False
    assert manifest["analysis_default_switch"] is True
    assert manifest["manuscript_eligible"] is False
    assert len(manifest["entries"]) == 2


def test_analysis_entrypoints_use_versioned_prod_v2_defaults():
    jump = (ROOT / "scripts/analysis/relaxtime/build_phase_guided_transport_xi001_jump_analysis.py").read_text(
        encoding="utf-8"
    )
    mechanism = (ROOT / "scripts/analysis/relaxtime/phase_guided_p128_mechanism_scan.jl").read_text(
        encoding="utf-8"
    )
    assert f'CASE = "{CASE}"' in jump
    assert f'const CASE = "{CASE}"' in mechanism
    assert "phase_guided_transport_p128_xi001_analysis_v2" in jump
    assert "phase_guided_transport_p128_xi001_analysis_v2" in mechanism
