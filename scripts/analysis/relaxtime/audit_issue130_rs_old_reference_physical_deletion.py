#!/usr/bin/env python3
"""Validate the Issue #130 RS legacy-snapshot deletion proposal.

This is a solver-free governance check.  It verifies the exact deletion
allowlist, confirms that every target is absent from the working tree, checks
that the path-retirement merge still contains the recorded recovery trees,
and ensures that current ``prod_v2`` consumers remain available.  It never
deletes files and never calls Julia or a numerical solver.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path, PureWindowsPath
from typing import Any


SCHEMA_VERSION = "issue130_rs_old_reference_physical_deletion_v1"
PACKAGE_ROOT = Path("docs/analysis/relaxtime/issue130_rs_old_reference_physical_deletion_v1")
MANIFEST_NAME = "deletion_manifest.json"
ALLOWLIST_NAME = "deletion_allowlist.csv"
CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
RECOVERY_REF = "74b53b47ebcca2b292cee72f70a70a84b0d2eea5"
LEGACY_SNAPSHOT_VERSION = "legacy_prod_v1_snapshot_v1"
LEGACY_CASE = "first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1"
CURRENT_CASE = "first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2"
MODES = ("mode_a_fixed_muB_phase_scaled", "mode_b_fixed_T_sparse_muB")
RESULT_ROOT = Path("data/outputs/results/relaxtime/transport/phase_guided")
FIGURE_ROOT = Path("data/outputs/figures/relaxtime/transport/phase_guided")
ACTIVE_RS_CONSUMERS = (
    "scripts/analysis/relaxtime/build_phase_guided_transport_xi001_jump_analysis.py",
    "scripts/analysis/relaxtime/phase_guided_p128_mechanism_scan.jl",
    "scripts/relaxtime/run_phase_guided_transport_scan.jl",
    "scripts/relaxtime/phase_reference_adapter.jl",
)


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _contract_relative_sort_key(relative: str) -> tuple[str, str]:
    """Match the case-insensitive Windows ``Path`` ordering used at import.

    The historical result/figure manifests were created on Windows with
    ``sorted(Path(...))``.  ``git ls-tree`` uses bytewise path ordering,
    which differs for names such as ``README.md`` and
    ``channel_diagnostics.csv``.  Keep the recorded tree hashes reproducible
    across the local Windows checkout and Linux Actions runners by making the
    ordering explicit here.
    """

    normalized = PureWindowsPath(relative).as_posix()
    return normalized.casefold(), normalized


def tree_hash_from_git(repo_root: Path, ref: str, prefix: str) -> tuple[int, int, str]:
    names = _git_tree_files(repo_root, ref, prefix)
    digest = hashlib.sha256()
    total_bytes = 0
    prefix_clean = prefix.rstrip("/")
    relative_names = []
    for name in names:
        if name.startswith(prefix_clean + "/"):
            relative = name[len(prefix_clean) + 1 :]
            full_name = name
        else:
            # Be tolerant of Git implementations that return paths relative
            # to a directory pathspec rather than repository-relative paths.
            relative = name
            full_name = f"{prefix_clean}/{name}"
        relative_names.append((relative, full_name))
    for relative, full_name in sorted(
        relative_names, key=lambda item: _contract_relative_sort_key(item[0])
    ):
        payload = _git_show_bytes(repo_root, ref, full_name)
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(sha256_bytes(payload).encode("ascii"))
        digest.update(b"\n")
        total_bytes += len(payload)
    return len(names), total_bytes, digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def read_allowlist(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    required = {
        "path_type", "layer", "mode", "case_slug", "path", "file_count", "bytes",
        "tree_sha256", "manifest_sha256", "pre_delete_sha256", "recovery_ref",
        "expected_after", "status", "reason",
    }
    if not rows:
        raise ValueError(f"empty deletion allowlist: {path}")
    missing = sorted(required - set(rows[0]))
    if missing:
        raise ValueError(f"deletion allowlist missing fields: {missing}")
    return [{str(key): str(value or "") for key, value in row.items()} for row in rows]


def _git_bytes(repo_root: Path, *args: str) -> bytes:
    try:
        return subprocess.check_output(["git", *args], cwd=repo_root, stderr=subprocess.DEVNULL)
    except (OSError, subprocess.CalledProcessError) as exc:
        raise ValueError(f"git command failed: git {' '.join(args)}") from exc


def _git_show_bytes(repo_root: Path, ref: str, path: str) -> bytes:
    return _git_bytes(repo_root, "show", f"{ref}:{path}")


def _git_tree_files(repo_root: Path, ref: str, prefix: str) -> list[str]:
    raw = _git_bytes(repo_root, "ls-tree", "-r", "--name-only", ref, "--", prefix)
    return [line for line in raw.decode("utf-8").splitlines() if line]


def _parse_int(value: str, field: str, path: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid {field} for {path}: {value!r}") from exc


def _repo_path(repo_root: Path, relative: str) -> Path:
    path = (repo_root / Path(relative.replace("/", "\\"))).resolve()
    root = repo_root.resolve()
    if path != root and root not in path.parents:
        raise ValueError(f"allowlist path escapes repository: {relative}")
    return path


def validate_allowlist(repo_root: Path, manifest: dict[str, Any], rows: list[dict[str, str]]) -> dict[str, Any]:
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("deletion manifest schema mismatch")
    if manifest.get("path_retirement_merge_sha") != RECOVERY_REF:
        raise ValueError("unexpected recovery ref")
    if manifest.get("calculation_sha") != CALCULATION_SHA:
        raise ValueError("calculation SHA mismatch")
    if manifest.get("physical_deletion_applied_in_branch") is not True:
        raise ValueError("branch does not record applied physical deletion")
    if manifest.get("merge_authorization_required") is not True:
        raise ValueError("merge authorization gate is missing")
    if manifest.get("fallback_available_after_merge") is not False:
        raise ValueError("fallback-after-merge boundary is not explicit")
    if manifest.get("rollback_available_after_merge") is not False:
        raise ValueError("rollback-after-merge boundary is not explicit")

    expected_scope = {
        "data/outputs/results/relaxtime/transport/phase_guided/legacy_prod_v1_snapshot_v1",
        "data/outputs/figures/relaxtime/transport/phase_guided/legacy_prod_v1_snapshot_v1",
    }
    if set(manifest.get("scope", [])) != expected_scope:
        raise ValueError("deletion scope drifted")
    if len(rows) != 6:
        raise ValueError(f"expected six allowlist rows, got {len(rows)}")

    deleted_files = 0
    deleted_bytes = 0
    directory_rows = 0
    for row in rows:
        path = row["path"]
        if row["expected_after"] != "absent" or row["status"] != "proposed":
            raise ValueError(f"invalid deletion state for {path}")
        if row["recovery_ref"] != RECOVERY_REF:
            raise ValueError(f"recovery ref mismatch for {path}")
        working = _repo_path(repo_root, path)
        if working.exists():
            raise ValueError(f"allowlisted deletion target still exists: {path}")

        expected_count = _parse_int(row["file_count"], "file_count", path)
        expected_bytes = _parse_int(row["bytes"], "bytes", path)
        if row["path_type"] == "file":
            if expected_count != 1:
                raise ValueError(f"file row must have file_count=1: {path}")
            names = _git_tree_files(repo_root, RECOVERY_REF, path)
            if names != [path]:
                raise ValueError(f"recovery file missing or ambiguous: {path}")
            payload = _git_show_bytes(repo_root, RECOVERY_REF, path)
            if len(payload) != expected_bytes or sha256_bytes(payload) != row["pre_delete_sha256"]:
                raise ValueError(f"recovery file hash/size mismatch: {path}")
            deleted_files += 1
            deleted_bytes += expected_bytes
        elif row["path_type"] == "directory":
            directory_rows += 1
            names = _git_tree_files(repo_root, RECOVERY_REF, path)
            count, total, tree_sha = tree_hash_from_git(repo_root, RECOVERY_REF, path)
            if count != expected_count or total != expected_bytes or tree_sha != row["tree_sha256"]:
                raise ValueError(f"recovery tree hash/count mismatch: {path}")
            deleted_files += count
            deleted_bytes += total
        else:
            raise ValueError(f"unsupported allowlist path_type for {path}: {row['path_type']}")

    if directory_rows != 4:
        raise ValueError(f"expected four directory rows, got {directory_rows}")
    if deleted_files != int(manifest.get("deleted_file_count", -1)):
        raise ValueError("deleted file count does not match manifest")
    if deleted_bytes != int(manifest.get("deleted_bytes", -1)):
        raise ValueError("deleted byte count does not match manifest")
    return {"rows": len(rows), "deleted_files": deleted_files, "deleted_bytes": deleted_bytes}


def validate_registry(repo_root: Path, manifest_rel: str) -> dict[str, Any]:
    registry_path = repo_root / RESULT_ROOT / "production_registry.json"
    registry = read_json(registry_path)
    entries = [
        entry for entry in registry.get("entries", [])
        if isinstance(entry, dict) and entry.get("case_slug") == LEGACY_CASE
    ]
    if {entry.get("mode") for entry in entries} != set(MODES):
        raise ValueError("registry legacy entries do not cover both modes")
    for entry in entries:
        if entry.get("path_status") != "physically_deleted":
            raise ValueError(f"registry path status is not physically_deleted: {entry.get('mode')}")
        if entry.get("fallback_available") is not False or entry.get("rollback_available") is not False:
            raise ValueError(f"registry fallback boundary is not closed: {entry.get('mode')}")
        if entry.get("physical_deletion_manifest_path") != manifest_rel:
            raise ValueError(f"registry deletion manifest pointer mismatch: {entry.get('mode')}")
        if LEGACY_SNAPSHOT_VERSION not in str(entry.get("result_path", "")):
            raise ValueError("registry must retain the deleted path as historical locator")
    return {"legacy_entries": len(entries), "current_prod_v2_untouched": True}


def validate_current_consumers(repo_root: Path) -> dict[str, Any]:
    old_tokens = (LEGACY_CASE, LEGACY_SNAPSHOT_VERSION)
    checked = []
    for relative in ACTIVE_RS_CONSUMERS:
        path = repo_root / Path(relative.replace("/", "\\"))
        if not path.is_file():
            raise ValueError(f"missing active consumer: {relative}")
        text = path.read_text(encoding="utf-8")
        if any(token in text for token in old_tokens):
            raise ValueError(f"active consumer retains deleted RS snapshot token: {relative}")
        checked.append(relative)
    for mode in MODES:
        current_result = repo_root / RESULT_ROOT / mode / CURRENT_CASE
        current_figure = repo_root / FIGURE_ROOT / mode / CURRENT_CASE
        if not (current_result / "manifest.json").is_file():
            raise ValueError(f"current prod_v2 result missing: {current_result}")
        if not (current_figure / "plot_manifest.json").is_file():
            raise ValueError(f"current prod_v2 figure missing: {current_figure}")
    return {"active_consumers_checked": checked, "current_modes_checked": list(MODES)}


def validate(repo_root: Path) -> dict[str, Any]:
    package = repo_root / PACKAGE_ROOT
    manifest_path = package / MANIFEST_NAME
    allowlist_path = package / ALLOWLIST_NAME
    manifest = read_json(manifest_path)
    rows = read_allowlist(allowlist_path)
    manifest_rel = (PACKAGE_ROOT / MANIFEST_NAME).as_posix()
    return {
        "schema_version": SCHEMA_VERSION,
        "verdict": "physical_deletion_proposal_valid",
        "manifest": validate_allowlist(repo_root, manifest, rows),
        "registry": validate_registry(repo_root, manifest_rel),
        "consumers": validate_current_consumers(repo_root),
        "solver_called": False,
        "production_write": False,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    report = validate(args.repo_root.resolve())
    print(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
