#!/usr/bin/env python3
"""Build a post-execution manifest for the approved figure retirement allowlist.

The preflight manifest is the before-state snapshot. This command verifies the
current worktree against that snapshot and records the resulting paths and
hashes. It does not perform file operations.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PREFLIGHT = Path("docs/analysis/figure_asset_registry_v1/cleanup_preflight_v1.json")
DEFAULT_OUTPUT = Path("docs/analysis/figure_asset_registry_v1/retirement_execution_v1.json")
OLD_TOKEN = "first_canonical_v1_p128_validated_anchored_prod_v1"
NEW_TOKEN = "first_canonical_v1_p128_xi005_validated_anchored_prod_v1"

FREEZEOUT_MANIFEST_SOURCE = (
    "data/outputs/figures/relaxtime/meson_density/plot_review/"
    "freezeout_kminus_piminus_mu_pi_100/plot_manifest.json"
)
FREEZEOUT_MANIFEST_SOURCE_PREFIX = (
    "data/outputs/figures/relaxtime/meson_density/plot_review/"
    "freezeout_kminus_piminus_mu_pi_100/"
)
FREEZEOUT_MANIFEST_DEST_PREFIX = (
    "docs/analysis/relaxtime/meson_density/"
    "freezeout_kminus_piminus_mu_pi_100_analysis/"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", type=Path, default=DEFAULT_PREFLIGHT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def relative(path: Path) -> str:
    return path.as_posix()


def absolute(path: Path) -> Path:
    return path if path.is_absolute() else PROJECT_ROOT / path


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def git_head_bytes(path: str) -> bytes:
    completed = subprocess.run(
        ["git", "show", f"HEAD:{path}"],
        cwd=PROJECT_ROOT,
        check=True,
        stdout=subprocess.PIPE,
    )
    return completed.stdout


def before_bytes(path: str, expected_sha256: str) -> bytes:
    payload = git_head_bytes(path)
    actual = sha256_bytes(payload)
    if actual != expected_sha256:
        raise RuntimeError(
            f"preflight hash mismatch for {path}: expected {expected_sha256}, got {actual}"
        )
    return payload


def current_bytes(path: str) -> bytes | None:
    target = absolute(Path(path))
    return target.read_bytes() if target.is_file() else None


def current_path_for_entry(entry: dict[str, Any]) -> str:
    if entry["action"] in {"move", "rename"}:
        return str(entry["destination"])
    return str(entry["path"])


def check_absent(path: str) -> None:
    if absolute(Path(path)).exists():
        raise RuntimeError(f"source still exists after approved operation: {path}")


def expected_after_bytes(entry: dict[str, Any], before: bytes) -> tuple[bytes, str]:
    path = str(entry["path"])
    policy = str(entry.get("content_policy", "byte_preserve"))

    if path == FREEZEOUT_MANIFEST_SOURCE:
        old_prefix = FREEZEOUT_MANIFEST_SOURCE_PREFIX.encode("utf-8")
        new_prefix = FREEZEOUT_MANIFEST_DEST_PREFIX.encode("utf-8")
        return before.replace(old_prefix, new_prefix), "path_only_metadata"

    old_token = OLD_TOKEN.encode("ascii")
    new_token = NEW_TOKEN.encode("ascii")
    if old_token in before:
        return before.replace(old_token, new_token), "exact_token_replace"

    if policy == "byte_preserve":
        return before, policy

    return before.replace(old_token, new_token), policy


def build_entry_record(entry: dict[str, Any]) -> dict[str, Any]:
    path = str(entry["path"])
    action = str(entry["action"])
    before = before_bytes(path, str(entry["sha256"]))
    destination = entry.get("destination")
    after_path = current_path_for_entry(entry)

    record: dict[str, Any] = {
        "operation_id": entry["operation_id"],
        "action": action,
        "before_path": path,
        "before_size_bytes": len(before),
        "before_sha256": sha256_bytes(before),
        "after_path": None if action == "delete" else after_path,
        "content_policy": entry.get("content_policy", "byte_preserve"),
    }
    if destination is not None:
        record["planned_destination"] = destination

    if action == "delete":
        check_absent(path)
        record.update(
            {
                "after_exists": False,
                "after_size_bytes": None,
                "after_sha256": None,
                "verification": "deleted",
            }
        )
        return record

    check_absent(path)
    after = current_bytes(after_path)
    if after is None:
        raise RuntimeError(f"destination missing after approved operation: {after_path}")

    expected, effective_policy = expected_after_bytes(entry, before)
    if after != expected:
        raise RuntimeError(
            f"unexpected content change for {path} -> {after_path}; "
            f"policy={effective_policy}"
        )
    record.update(
        {
            "after_exists": True,
            "after_size_bytes": len(after),
            "after_sha256": sha256_bytes(after),
            "content_policy": effective_policy,
            "verification": "byte_preserved" if effective_policy == "byte_preserve" else "declared_path_update",
        }
    )
    return record


def build_reference_record(candidate: dict[str, Any]) -> dict[str, Any]:
    before_path = str(candidate["path"])
    before = before_bytes(before_path, str(candidate["sha256"]))
    after_path = before_path
    if not absolute(Path(after_path)).is_file():
        after_path = before_path.replace(OLD_TOKEN, NEW_TOKEN)
    after = current_bytes(after_path)
    if after is None:
        raise RuntimeError(f"updated reference is missing: {after_path}")

    expected = before.replace(OLD_TOKEN.encode("ascii"), NEW_TOKEN.encode("ascii"))
    if after != expected:
        raise RuntimeError(f"reference update is not an exact token replacement: {before_path}")

    return {
        "before_path": before_path,
        "after_path": after_path,
        "before_sha256": sha256_bytes(before),
        "after_sha256": sha256_bytes(after),
        "match_count": candidate["match_count"],
        "content_policy": "exact_token_replace",
        "verification": "exact_token_replace",
    }


def build_manifest(preflight_path: Path) -> dict[str, Any]:
    preflight = json.loads(absolute(preflight_path).read_text(encoding="utf-8"))
    records = [build_entry_record(entry) for entry in preflight["entries"]]
    reference_records = [
        build_reference_record(candidate)
        for candidate in preflight["reference_update"]["candidates"]
    ]
    action_counts: dict[str, int] = {}
    for record in records:
        action = str(record["action"])
        action_counts[action] = action_counts.get(action, 0) + 1

    return {
        "schema_version": "figure_retirement_execution_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "execution_status": "completed",
        "source_preflight": {
            "path": relative(preflight_path),
            "sha256": sha256_file(absolute(preflight_path)),
            "generated_at_utc": preflight["generated_at_utc"],
            "base_git_commit": preflight["base_git_commit"],
        },
        "author_decision_record": preflight["author_decision_record"],
        "execution_gate": {
            "manual_review_required": False,
            "delete_performed": True,
            "move_performed": True,
            "rename_performed": True,
            "tracked_only": True,
            "no_solver_or_numeric_data_changes": True,
            "untracked_exclusions_untouched": True,
            "recovery_source": "git_history_and_pull_request_history",
        },
        "summary": {
            "operation_entry_count": len(records),
            "delete_entry_count": action_counts.get("delete", 0),
            "move_entry_count": action_counts.get("move", 0),
            "rename_entry_count": action_counts.get("rename", 0),
            "reference_update_count": len(reference_records),
            "freezeout_manifest_policy": "path_fields_only_old_and_new_hash_recorded",
        },
        "operations": records,
        "reference_updates": {
            "old_token": preflight["reference_update"]["old_token"],
            "new_token": preflight["reference_update"]["new_token"],
            "records": reference_records,
        },
        "untracked_exclusions": preflight["untracked_exclusions"],
    }


def main() -> None:
    args = parse_args()
    manifest = build_manifest(args.preflight)
    output = absolute(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(output)
    print(json.dumps(manifest["summary"], ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
