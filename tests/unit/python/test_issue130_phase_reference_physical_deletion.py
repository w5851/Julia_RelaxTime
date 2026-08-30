from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).parents[3]
SCRIPT = ROOT / "scripts/analysis/pnjl/audit_issue130_phase_reference_physical_deletion.py"
SPEC = importlib.util.spec_from_file_location("issue130_pnjl_phase_reference_physical_deletion", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_pnjl_physical_deletion_is_solver_free_and_exactly_allowlisted() -> None:
    report = MODULE.validate(ROOT)
    assert report["verdict"] == "physical_deletion_proposal_valid"
    assert report["solver_called"] is False
    assert report["production_write"] is False
    assert report["manifest"] == {
        "rows": 8,
        "deleted_files": 8,
        "deleted_bytes": 50749,
    }
    assert report["consumers"]["active_consumer_count"] == len(MODULE.ACTIVE_CONSUMERS)


def test_candidate_and_separate_crossover_reference_remain() -> None:
    candidate = ROOT / MODULE.CANDIDATE_ROOT
    assert (candidate / "manifest.json").is_file()
    assert (candidate / "accepted" / "manifest.json").is_file()
    assert (candidate / "strict" / "manifest.json").is_file()
    assert (ROOT / MODULE.FIXED_CROSSOVER).is_file()
    assert not (ROOT / MODULE.LEGACY_ROOT).exists()


def test_deleted_snapshot_is_recoverable_from_merge_ref() -> None:
    package = ROOT / MODULE.PACKAGE_ROOT
    rows = MODULE._read_allowlist(package / MODULE.ALLOWLIST_NAME)
    for row in rows:
        assert not MODULE._repo_path(ROOT, row["path"]).exists()
        assert MODULE._git_tree_files(ROOT, MODULE.RECOVERY_REF, row["path"]) == [row["path"]]
