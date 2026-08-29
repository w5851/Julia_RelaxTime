from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).parents[3]
REFERENCE_ROOT = ROOT / "data" / "reference" / "pnjl"
SNAPSHOT_ROOT = REFERENCE_ROOT / "legacy_phase_reference_v1"


def test_retirement_snapshot_manifest_matches_bytes() -> None:
    manifest = json.loads((SNAPSHOT_ROOT / "RETIREMENT_MANIFEST.json").read_text(encoding="utf-8"))
    assert manifest["schema_version"] == "pnjl_legacy_phase_reference_retirement_v1"
    assert manifest["status"] == "retired_canonical_snapshot"
    assert manifest["runtime_canonical_source"].endswith("issue130_phase_reference_v1/strict")
    assert manifest["fallback_and_rollback"] == "explicit_snapshot_only"
    assert manifest["canonical_root_status"] == "dense_legacy_paths_absent"
    assert manifest["solver_called"] is False

    for record in manifest["files"]:
        assert record["source_path"] == f"data/reference/pnjl/{record['path']}"
        path = SNAPSHOT_ROOT / record["path"]
        payload = path.read_bytes()
        assert len(payload) == record["bytes"]
        assert hashlib.sha256(payload).hexdigest() == record["sha256"]


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
    assert all((SNAPSHOT_ROOT / name).is_file() for name in retired)

    # This fixed-point table is a separate historical input, not the dense
    # Issue #130 legacy bundle retired by this change.
    assert (REFERENCE_ROOT / "crossover.csv").is_file()


def test_runtime_consumers_default_to_accepted_and_do_not_fallback_to_legacy() -> None:
    consumers = [
        ROOT / "scripts" / "relaxtime" / "phase_reference_adapter.jl",
        ROOT / "scripts" / "relaxtime" / "run_gap_transport_scan.jl",
        ROOT / "scripts" / "pnjl" / "plot_phase_diagram.py",
    ]
    for path in consumers:
        text = path.read_text(encoding="utf-8")
        assert "issue130_phase_reference_v2" in text, path
    adapter = (ROOT / "scripts" / "relaxtime" / "phase_reference_adapter.jl").read_text(encoding="utf-8")
    runtime_scan = (ROOT / "scripts" / "relaxtime" / "run_gap_transport_scan.jl").read_text(encoding="utf-8")
    assert "load_phase_reference_runtime_with_fallback" not in adapter
    assert "phase-reference-mode legacy" not in runtime_scan


def test_historical_validation_tool_remains_explicitly_snapshot_bound() -> None:
    validator = (ROOT / "scripts" / "pnjl" / "validate_phase_data.py").read_text(encoding="utf-8")
    assert "legacy_phase_reference_v1" in validator
