from __future__ import annotations

import json
from pathlib import Path

from scripts.analysis.pnjl import audit_issue130_phase_reference_physical_deletion as deletion

ROOT = Path(__file__).parents[3]
REFERENCE_ROOT = ROOT / "data" / "reference" / "pnjl"


def test_retirement_snapshot_is_deleted_but_recoverable() -> None:
    package = ROOT / deletion.PACKAGE_ROOT
    manifest = json.loads((package / deletion.MANIFEST_NAME).read_text(encoding="utf-8"))
    assert manifest["schema_version"] == deletion.SCHEMA_VERSION
    assert manifest["physical_deletion_applied_in_branch"] is True
    assert manifest["fallback_available_after_merge"] is False
    assert manifest["rollback_available_after_merge"] is False
    rows = deletion._read_allowlist(package / deletion.ALLOWLIST_NAME)
    assert len(rows) == 8
    assert all(not deletion._repo_path(ROOT, row["path"]).exists() for row in rows)
    assert all(deletion._git_tree_files(ROOT, deletion.RECOVERY_REF, row["path"]) for row in rows)


def test_dense_legacy_files_are_retired_from_canonical_root() -> None:
    retired = {
        "boundary.csv",
        "cep.csv",
        "crossover_dense.csv",
        "crossover_dense.meta.json",
        "phase_reference_dense_manifest.json",
        "spinodals.csv",
    }
    assert all(not (REFERENCE_ROOT / name).exists() for name in retired)
    assert not (REFERENCE_ROOT / "legacy_phase_reference_v1").exists()

    # This fixed-point table is a separate historical input, not the dense
    # Issue #130 legacy bundle retired by this change.
    assert (REFERENCE_ROOT / "crossover.csv").is_file()


def test_runtime_consumers_default_to_accepted_and_do_not_fallback_to_legacy() -> None:
    consumers = [
        ROOT / "scripts" / "relaxtime" / "phase_reference_adapter.jl",
        ROOT / "scripts" / "relaxtime" / "run_gap_transport_scan.jl",
        ROOT / "scripts" / "relaxtime" / "gap_transport_scan_phase_equilibrium.jl",
        ROOT / "scripts" / "relaxtime" / "generate_xi_smoothness_params.jl",
        ROOT / "scripts" / "relaxtime" / "xi_smoothness_sampling_lib.jl",
        ROOT / "scripts" / "pnjl" / "plot_phase_diagram.py",
        ROOT / "scripts" / "pnjl" / "validate_phase_data.py",
        ROOT / "scripts" / "dev" / "export_pnjl_chi_b_taylordiff_baseline.jl",
    ]
    for path in consumers:
        text = path.read_text(encoding="utf-8")
        assert "issue130_phase_reference_v2" in text, path
        assert "legacy_phase_reference_v1" not in text, path
    adapter = (ROOT / "scripts" / "relaxtime" / "phase_reference_adapter.jl").read_text(encoding="utf-8")
    runtime_scan = (ROOT / "scripts" / "relaxtime" / "run_gap_transport_scan.jl").read_text(encoding="utf-8")
    assert "load_phase_reference_runtime_with_fallback" not in adapter
    assert "phase-reference-mode legacy" not in runtime_scan


def test_validation_tool_defaults_to_accepted_and_keeps_explicit_path_override() -> None:
    validator = (ROOT / "scripts" / "pnjl" / "validate_phase_data.py").read_text(encoding="utf-8")
    assert "issue130_phase_reference_v2" in validator
    assert "--phase-reference-root" in validator
    assert "legacy_phase_reference_v1" not in validator


def test_validation_adapter_maps_accepted_tables_without_writing_files() -> None:
    from scripts.pnjl.validate_phase_data import (
        DEFAULT_PHASE_REFERENCE_LAYER,
        DEFAULT_PHASE_REFERENCE_ROOT,
        load_candidate_validation_data,
    )

    assert DEFAULT_PHASE_REFERENCE_LAYER == "accepted"
    data = load_candidate_validation_data(DEFAULT_PHASE_REFERENCE_ROOT)
    boundary_headers, boundary_rows = data["boundary"]
    spinodal_headers, spinodal_rows = data["spinodals"]
    crossover_headers, crossover_rows = data["crossover"]
    assert {"xi", "T_MeV", "mu_transition_MeV"} <= set(boundary_headers)
    assert {"xi", "T_MeV", "mu_spinodal_hadron_MeV"} <= set(spinodal_headers)
    assert {"xi", "mu_MeV", "T_crossover_chiral_MeV"} <= set(crossover_headers)
    assert len(boundary_rows) == 12537
    assert len(spinodal_rows) == 11989
    assert len(crossover_rows) == 3135
