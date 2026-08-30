#!/usr/bin/env python3
"""Validate the Issue #130 PNJL legacy phase-reference deletion proposal.

This is a solver-free governance validator.  The deletion itself is applied
in the proposal branch; this script checks that the exact allowlist is absent,
that the pre-delete snapshot is recoverable from the path-retirement merge,
and that the accepted/strict candidate remains intact.  It never restores,
deletes, rewrites, or numerically evaluates a phase-reference file.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path, PureWindowsPath
from typing import Any


SCHEMA_VERSION = "issue130_pnjl_phase_reference_physical_deletion_v1"
PACKAGE_ROOT = Path(
    "docs/analysis/pnjl/phase_reference/issue130_phase_reference_physical_deletion_v1"
)
MANIFEST_NAME = "deletion_manifest.json"
ALLOWLIST_NAME = "deletion_allowlist.csv"
LEGACY_ROOT = Path("data/reference/pnjl/legacy_phase_reference_v1")
CANDIDATE_ROOT = Path("data/reference/pnjl/issue130_phase_reference_v2")
ACCEPTED_MANIFEST = CANDIDATE_ROOT / "accepted" / "manifest.json"
STRICT_MANIFEST = CANDIDATE_ROOT / "strict" / "manifest.json"
FIXED_CROSSOVER = Path("data/reference/pnjl/crossover.csv")
CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
# PR #284 is the last merge containing the complete legacy snapshot.
RECOVERY_REF = "9aa4c313901ca0c91e851f58514e3df9aa124df4"
LEGACY_FILES = (
    "boundary.csv",
    "cep.csv",
    "crossover_dense.csv",
    "crossover_dense.meta.json",
    "phase_reference_dense_manifest.json",
    "README.md",
    "RETIREMENT_MANIFEST.json",
    "spinodals.csv",
)

# These are the executable default consumers whose source contracts must not
# depend on the deleted snapshot.  Historical audit/replay scripts are not in
# this list: they remain usable only when an external recovery tree is given.
ACTIVE_CONSUMERS = (
    "scripts/relaxtime/phase_reference_adapter.jl",
    "scripts/relaxtime/run_gap_transport_scan.jl",
    "scripts/relaxtime/gap_transport_scan_phase_equilibrium.jl",
    "scripts/relaxtime/generate_xi_smoothness_params.jl",
    "scripts/relaxtime/xi_smoothness_sampling_lib.jl",
    "scripts/pnjl/plot_phase_diagram.py",
    "scripts/pnjl/validate_phase_data.py",
    "scripts/dev/export_pnjl_chi_b_taylordiff_baseline.jl",
)


def _git_bytes(repo_root: Path, *args: str) -> bytes:
    try:
        return subprocess.check_output(["git", *args], cwd=repo_root, stderr=subprocess.DEVNULL)
    except (OSError, subprocess.CalledProcessError) as exc:
        raise ValueError(f"git command failed: git {' '.join(args)}") from exc


def _git_tree_files(repo_root: Path, ref: str, prefix: str) -> list[str]:
    raw = _git_bytes(repo_root, "ls-tree", "-r", "--name-only", ref, "--", prefix)
    return [line for line in raw.decode("utf-8").splitlines() if line]


def _git_show_bytes(repo_root: Path, ref: str, path: str) -> bytes:
    return _git_bytes(repo_root, "show", f"{ref}:{path}")


def _sort_key(relative: str) -> tuple[str, str]:
    normalized = PureWindowsPath(relative).as_posix()
    return normalized.casefold(), normalized


def _tree_hash_from_git(repo_root: Path, ref: str, prefix: str) -> tuple[int, int, str]:
    names = _git_tree_files(repo_root, ref, prefix)
    clean = prefix.rstrip("/")
    entries: list[tuple[str, str]] = []
    for name in names:
        if name.startswith(clean + "/"):
            entries.append((name[len(clean) + 1 :], name))
        else:
            entries.append((name, f"{clean}/{name}"))
    digest = hashlib.sha256()
    total_bytes = 0
    for relative, full_name in sorted(entries, key=lambda item: _sort_key(item[0])):
        payload = _git_show_bytes(repo_root, ref, full_name)
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(hashlib.sha256(payload).hexdigest().encode("ascii"))
        digest.update(b"\n")
        total_bytes += len(payload)
    return len(entries), total_bytes, digest.hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _read_allowlist(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    required = {
        "path_type",
        "path",
        "bytes",
        "pre_delete_sha256",
        "recovery_ref",
        "expected_after",
        "status",
        "reason",
    }
    if not rows:
        raise ValueError(f"empty deletion allowlist: {path}")
    missing = sorted(required - set(rows[0]))
    if missing:
        raise ValueError(f"deletion allowlist missing fields: {missing}")
    return [{str(key): str(value or "") for key, value in row.items()} for row in rows]


def _repo_path(repo_root: Path, relative: str) -> Path:
    path = (repo_root / Path(relative.replace("/", "\\"))).resolve()
    root = repo_root.resolve()
    if path != root and root not in path.parents:
        raise ValueError(f"allowlist path escapes repository: {relative}")
    return path


def _parse_int(value: str, field: str, path: str) -> int:
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid {field} for {path}: {value!r}") from exc


def validate_allowlist(repo_root: Path, manifest: dict[str, Any], rows: list[dict[str, str]]) -> dict[str, Any]:
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("deletion manifest schema mismatch")
    if manifest.get("recovery_ref") != RECOVERY_REF:
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

    expected_paths = {f"{LEGACY_ROOT.as_posix()}/{name}" for name in LEGACY_FILES}
    if {row["path"] for row in rows} != expected_paths:
        raise ValueError("deletion allowlist paths drifted")
    if len(rows) != len(LEGACY_FILES):
        raise ValueError(f"expected {len(LEGACY_FILES)} allowlist rows, got {len(rows)}")

    deleted_bytes = 0
    for row in rows:
        path = row["path"]
        if row["path_type"] != "file":
            raise ValueError(f"PNJL legacy allowlist only permits files: {path}")
        if row["expected_after"] != "absent" or row["status"] != "proposed":
            raise ValueError(f"invalid deletion state for {path}")
        if row["recovery_ref"] != RECOVERY_REF:
            raise ValueError(f"recovery ref mismatch for {path}")
        if _repo_path(repo_root, path).exists():
            raise ValueError(f"allowlisted deletion target still exists: {path}")
        expected_bytes = _parse_int(row["bytes"], "bytes", path)
        names = _git_tree_files(repo_root, RECOVERY_REF, path)
        if names != [path]:
            raise ValueError(f"recovery file missing or ambiguous: {path}")
        payload = _git_show_bytes(repo_root, RECOVERY_REF, path)
        if len(payload) != expected_bytes or hashlib.sha256(payload).hexdigest() != row["pre_delete_sha256"]:
            raise ValueError(f"recovery file hash/size mismatch: {path}")
        deleted_bytes += expected_bytes

    if repo_root.joinpath(*LEGACY_ROOT.parts).exists():
        raise ValueError("legacy snapshot directory still exists after deletion")
    expected_count, expected_total, expected_tree = _tree_hash_from_git(repo_root, RECOVERY_REF, LEGACY_ROOT.as_posix())
    if expected_count != len(LEGACY_FILES):
        raise ValueError("recovery snapshot file count drifted")
    if expected_total != int(manifest.get("pre_delete_bytes", -1)):
        raise ValueError("pre-delete byte total mismatch")
    if expected_tree != manifest.get("pre_delete_tree_sha256"):
        raise ValueError("pre-delete tree hash mismatch")
    if deleted_bytes != int(manifest.get("deleted_bytes", -1)):
        raise ValueError("deleted byte total mismatch")
    return {"rows": len(rows), "deleted_files": len(rows), "deleted_bytes": deleted_bytes}


def validate_candidate(repo_root: Path) -> dict[str, Any]:
    candidate_manifest_path = repo_root / CANDIDATE_ROOT / "manifest.json"
    accepted_manifest_path = repo_root / ACCEPTED_MANIFEST
    strict_manifest_path = repo_root / STRICT_MANIFEST
    for path in (candidate_manifest_path, accepted_manifest_path, strict_manifest_path):
        if not path.is_file():
            raise ValueError(f"candidate manifest missing: {path}")
    candidate = _read_json(candidate_manifest_path)
    accepted = _read_json(accepted_manifest_path)
    strict = _read_json(strict_manifest_path)
    if candidate.get("calculation_sha") != CALCULATION_SHA:
        raise ValueError("candidate calculation SHA mismatch")
    if accepted.get("promotion_status") != "accepted_for_downstream":
        raise ValueError("accepted layer is not promoted for downstream")
    if candidate.get("downstream_default_layer") != "accepted":
        raise ValueError("accepted is not the declared downstream default")
    if strict.get("layer") != "strict":
        raise ValueError("strict layer manifest mismatch")
    # ``source_manifest_sha256`` refers to the imported v1 source manifest and
    # is intentionally not the hash of this v2 manifest.  Validate only that
    # the immutable source provenance is present here.
    if not candidate.get("source_manifest_sha256"):
        raise ValueError("candidate source manifest provenance is missing")
    if not (repo_root / FIXED_CROSSOVER).is_file():
        raise ValueError("separate crossover.csv reference disappeared")
    return {
        "candidate_manifest_sha256": _sha256_file(candidate_manifest_path),
        "accepted_manifest_sha256": _sha256_file(accepted_manifest_path),
        "strict_manifest_sha256": _sha256_file(strict_manifest_path),
        "accepted_rows": accepted.get("rows", {}),
        "strict_rows": strict.get("rows", {}),
        "fixed_crossover_sha256": _sha256_file(repo_root / FIXED_CROSSOVER),
    }


def validate_consumers(repo_root: Path) -> dict[str, Any]:
    checked: list[str] = []
    for relative in ACTIVE_CONSUMERS:
        path = repo_root / Path(relative.replace("/", "\\"))
        if not path.is_file():
            raise ValueError(f"active consumer missing: {relative}")
        text = path.read_text(encoding="utf-8")
        if "legacy_phase_reference_v1" in text:
            raise ValueError(f"active consumer retains deleted legacy token: {relative}")
        checked.append(relative)
    return {"active_consumers_checked": checked, "active_consumer_count": len(checked)}


def validate(repo_root: Path) -> dict[str, Any]:
    package = repo_root / PACKAGE_ROOT
    manifest_path = package / MANIFEST_NAME
    allowlist_path = package / ALLOWLIST_NAME
    manifest = _read_json(manifest_path)
    rows = _read_allowlist(allowlist_path)
    return {
        "schema_version": SCHEMA_VERSION,
        "verdict": "physical_deletion_proposal_valid",
        "manifest": validate_allowlist(repo_root, manifest, rows),
        "candidate": validate_candidate(repo_root),
        "consumers": validate_consumers(repo_root),
        "solver_called": False,
        "production_write": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, required=True)
    args = parser.parse_args()
    report = validate(args.repo_root.resolve())
    print(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
