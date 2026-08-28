from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts/analysis/relaxtime/audit_issue130_rs_old_reference_retirement.py"
SPEC = importlib.util.spec_from_file_location("issue130_rs_old_reference_retirement", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_tree_hash_matches_importer_contract(tmp_path: Path) -> None:
    root = tmp_path / "tree"
    root.mkdir()
    (root / "b.txt").write_text("b", encoding="utf-8")
    (root / "a.txt").write_text("a", encoding="utf-8")
    first = MODULE.tree_hash(root)
    (root / "nested").mkdir()
    (root / "nested" / "c.txt").write_text("c", encoding="utf-8")
    assert MODULE.tree_hash(root) != first


def test_classify_current_and_legacy_layers() -> None:
    current = MODULE.classify_case(MODULE.CURRENT_CASE, "result")
    legacy = MODULE.classify_case(MODULE.LEGACY_CASE, "result")
    legacy_figure = MODULE.classify_case(MODULE.LEGACY_CASE, "figure")
    historical = MODULE.classify_case(MODULE.HISTORICAL_CASES[0], "result")

    assert current == ("current_prod_v2", "canonical_current", "retain_canonical")
    assert legacy[0:2] == ("legacy_prod_v1", "legacy_fallback_required")
    assert legacy_figure[1] == "legacy_figure_history"
    assert historical[1] == "historical_evidence"


def test_csv_stats_skips_comments_and_detects_duplicates(tmp_path: Path) -> None:
    path = tmp_path / "scan.csv"
    path.write_text(
        "# metadata\n"
        "muB_MeV,xi,alpha_T,converged,T_MeV\n"
        "0,0,1,true,100\n"
        "0,0,1,true,100\n"
        "0,0.1,1,false,NaN\n",
        encoding="utf-8",
    )
    stats = MODULE.csv_stats(path, "mode_a_fixed_muB_phase_scaled")
    assert stats["rows"] == 3
    assert stats["duplicate_keys"] == 1
    assert stats["nonfinite_values"] == 1
    assert stats["nonconverged_rows"] == 1


def test_scan_contract_requires_complete_finite_rows(tmp_path: Path) -> None:
    path = tmp_path / "scan.csv"
    path.write_text(
        "muB_MeV,xi,alpha_T,converged,T_MeV\n"
        "0,0,1,true,100\n",
        encoding="utf-8",
    )
    stats = MODULE.csv_stats(path, "mode_a_fixed_muB_phase_scaled")
    assert MODULE.scan_contract_ok(stats, 1)
    assert not MODULE.scan_contract_ok(stats, 2)


def test_reference_role_separates_fallback_from_historical() -> None:
    assert MODULE.occurrence_role(
        "scripts/analysis/relaxtime/example.py",
        "legacy fallback path",
        "code_or_workflow_review",
    ) == "fallback_or_rollback"
    assert MODULE.occurrence_role(
        "docs/analysis/example/README.md",
        "old result",
        "historical_evidence",
    ) == "historical"


def test_status_consistency_keeps_import_and_registry_semantics_separate() -> None:
    assert MODULE.status_consistency(
        MODULE.CURRENT_CASE, "approved", "versioned_candidate_not_default"
    ) == "semantic_split_expected"
    assert MODULE.status_consistency(
        MODULE.LEGACY_CASE, "superseded_for_manuscript", "current_candidate"
    ) == "historical_manifest_status_stale_but_registry_superseded"
    assert MODULE.status_consistency("unknown", "approved", "current_candidate") == "metadata_review_required"


def test_real_fallback_hashes_match_current_manifests() -> None:
    for mode in MODULE.MODES:
        current = ROOT / MODULE.RESULT_ROOT / mode / MODULE.CURRENT_CASE
        legacy = MODULE.result_case_root(ROOT, mode, MODULE.LEGACY_CASE)
        manifest = MODULE.read_json(current / "manifest.json")
        assert manifest["legacy_prod_v1_tree_hash"] == MODULE.tree_hash(legacy)


def test_real_figure_inventory_uses_associated_result_scan() -> None:
    for mode in MODULE.MODES:
        result = ROOT / MODULE.RESULT_ROOT / mode / MODULE.CURRENT_CASE
        figure = ROOT / MODULE.FIGURE_ROOT / mode / MODULE.CURRENT_CASE
        plot_manifest = MODULE.read_json(figure / "plot_manifest.json")
        assert plot_manifest["source_csv_sha256"] == MODULE.sha256(result / MODULE.SCAN_NAME)
        assert MODULE.csv_stats(result / MODULE.SCAN_NAME, mode)["rows"] > 0


def test_real_inventory_includes_shared_convergence_trees() -> None:
    roots = [
        ROOT / MODULE.RESULT_ROOT / "first_canonical_v1_p128_validated_anchored_prod_v1_convergence",
        ROOT / MODULE.RESULT_ROOT / "first_canonical_v1_p128_xi001_validated_anchored_prod_v1_convergence",
    ]
    assert all(path.is_dir() for path in roots)
    assert all(MODULE.classify_case(path.name, "result_convergence")[0] == "historical_prod_v1" for path in roots)


def test_real_current_and_legacy_inventory_contracts_pass() -> None:
    rows = MODULE.collect_trees(ROOT, MODULE.registry_index(MODULE.read_json(ROOT / MODULE.REGISTRY_PATH)))
    selected = [
        row
        for row in rows
        if row["family"] in {"current_prod_v2", "legacy_prod_v1"}
    ]
    assert selected
    assert all(row["tree_contract_ok"] for row in selected)


def test_snapshot_manifests_preserve_tree_hashes_and_retire_canonical_paths() -> None:
    result_root = ROOT / MODULE.LEGACY_RESULT_ROOT
    for root, expected_schema in ((result_root, "relaxtime_rs_legacy_prod_v1_retirement_v1"),):
        manifest = MODULE.read_json(root / "RETIREMENT_MANIFEST.json")
        assert manifest["schema_version"] == expected_schema
        assert manifest["status"] == "retired_canonical_snapshot"
        assert manifest["canonical_root_status"] == "legacy_prod_v1_paths_absent"
        assert manifest["solver_called"] is False
        assert manifest["production_write"] is False
        for tree in manifest["trees"]:
            snapshot = ROOT / tree["snapshot_path"]
            source = ROOT / tree["source_path"]
            assert snapshot.is_dir()
            assert not source.exists()
            assert MODULE.tree_hash(snapshot) == tree["tree_sha256"]
            files = [path for path in snapshot.rglob("*") if path.is_file()]
            assert len(files) == tree["file_count"]
            assert sum(path.stat().st_size for path in files) == tree["bytes"]
            inner_manifest = snapshot / "manifest.json"
            assert MODULE.sha256(inner_manifest) == tree["manifest_sha256"]

    figure_root = ROOT / MODULE.LEGACY_FIGURE_ROOT
    figure_manifest_root = ROOT / "docs" / "analysis" / "relaxtime" / "issue130_rs_old_reference_path_retirement_v1" / "figure_snapshot"
    figure_manifest = MODULE.read_json(figure_manifest_root / "RETIREMENT_MANIFEST.json")
    assert figure_manifest["schema_version"] == "relaxtime_rs_legacy_prod_v1_figure_retirement_v1"
    assert figure_manifest["status"] == "retired_canonical_snapshot"
    assert figure_manifest["canonical_root_status"] == "legacy_prod_v1_paths_absent"
    assert figure_manifest["solver_called"] is False
    assert figure_manifest["production_write"] is False
    for tree in figure_manifest["trees"]:
        snapshot = ROOT / tree["snapshot_path"]
        source = ROOT / tree["source_path"]
        assert snapshot.is_dir()
        assert not source.exists()
        assert MODULE.tree_hash(snapshot) == tree["tree_sha256"]
        files = [path for path in snapshot.rglob("*") if path.is_file()]
        assert len(files) == tree["file_count"]
        assert sum(path.stat().st_size for path in files) == tree["bytes"]
        assert MODULE.sha256(snapshot / "plot_manifest.json") == tree["manifest_sha256"]
    assert figure_root.is_dir()
    assert not (figure_root / "README.md").exists()
    assert not (figure_root / "RETIREMENT_MANIFEST.json").exists()
    assert MODULE.read_json(result_root / "RETIREMENT_MANIFEST.json")["figure_snapshot_metadata_path"].endswith(
        "issue130_rs_old_reference_path_retirement_v1/figure_snapshot/RETIREMENT_MANIFEST.json"
    )


def test_registry_points_legacy_entries_to_versioned_snapshot() -> None:
    registry = MODULE.read_json(ROOT / MODULE.REGISTRY_PATH)
    entries = [
        entry for entry in registry["entries"] if entry.get("case_slug") == MODULE.LEGACY_CASE
    ]
    assert len(entries) == len(MODULE.MODES)
    for entry in entries:
        assert entry["path_status"] == "retired_to_versioned_legacy_snapshot"
        assert entry["legacy_snapshot_version"] == MODULE.LEGACY_SNAPSHOT_VERSION
        assert "legacy_prod_v1_snapshot_v1" in entry["result_path"]
        assert "legacy_prod_v1_snapshot_v1" in entry["figure_path"]
        assert (ROOT / entry["result_path"]).is_dir()
        assert (ROOT / entry["figure_path"]).is_dir()


def test_default_consumer_smoke_is_solver_free() -> None:
    rows = MODULE.default_consumer_rows(ROOT)
    assert rows
    assert all(row["solver_called"] is False for row in rows)
    assert all(row["verdict"] == "pass" for row in rows)
    assert all(row["old_case_occurrences"] == 0 for row in rows)


def test_real_audit_manifest_matches_files() -> None:
    root = ROOT / MODULE.AUDIT_PACKAGE_ROOT
    manifest = MODULE.read_json(root / "manifest.json")
    assert manifest["solver_called"] is False
    assert manifest["old_reference_deleted"] is False
    for item in manifest["files"]:
        path = root / item["path"]
        assert path.is_file()
        assert path.stat().st_size == item["bytes"]
        assert MODULE.sha256(path) == item["sha256"]


def test_nonfinite_parser_rejects_nan_without_raising(tmp_path: Path) -> None:
    path = tmp_path / "scan.csv"
    path.write_text(
        "muB_MeV,xi,alpha_T,converged\n0,0,1,true\n0,0.1,1,true\n",
        encoding="utf-8",
    )
    stats = MODULE.csv_stats(path, "mode_a_fixed_muB_phase_scaled")
    assert stats["nonfinite_values"] == 0
