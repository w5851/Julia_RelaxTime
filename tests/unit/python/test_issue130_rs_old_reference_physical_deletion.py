from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts/analysis/relaxtime/audit_issue130_rs_old_reference_physical_deletion.py"
SPEC = importlib.util.spec_from_file_location("issue130_rs_old_reference_physical_deletion", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_physical_deletion_proposal_is_solver_free_and_allowlisted() -> None:
    report = MODULE.validate(ROOT)
    assert report["verdict"] == "physical_deletion_proposal_valid"
    assert report["solver_called"] is False
    assert report["production_write"] is False
    assert report["manifest"] == {
        "rows": 6,
        "deleted_files": 112,
        "deleted_bytes": 41673063,
    }


def test_current_prod_v2_and_pnjl_candidate_reference_remain() -> None:
    for mode in MODULE.MODES:
        assert (
            ROOT / MODULE.RESULT_ROOT / mode / MODULE.CURRENT_CASE / "manifest.json"
        ).is_file()
        assert (
            ROOT / MODULE.FIGURE_ROOT / mode / MODULE.CURRENT_CASE / "plot_manifest.json"
        ).is_file()
    assert (ROOT / "data/reference/pnjl/issue130_phase_reference_v2").is_dir()


def test_deleted_snapshot_paths_are_absent_but_recoverable_from_merge() -> None:
    package = ROOT / MODULE.PACKAGE_ROOT
    rows = MODULE.read_allowlist(package / MODULE.ALLOWLIST_NAME)
    for row in rows:
        assert not MODULE._repo_path(ROOT, row["path"]).exists()
        assert row["recovery_ref"] == MODULE.RECOVERY_REF
        assert MODULE._git_tree_files(ROOT, MODULE.RECOVERY_REF, row["path"])


def test_registry_closes_rs_fallback_without_touching_pnjl_fallback() -> None:
    result = MODULE.validate_registry(
        ROOT,
        (MODULE.PACKAGE_ROOT / MODULE.MANIFEST_NAME).as_posix(),
    )
    assert result["legacy_entries"] == 2
    assert result["current_prod_v2_untouched"] is True
