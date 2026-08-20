#!/usr/bin/env python3
"""Audit analysis-package output metadata without changing evidence files."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Any


PACKAGE_SPECS = (
    (
        "phase_reference_current_state_freeze_v1",
        "docs/analysis/pnjl/phase_reference/phase_reference_current_state_freeze_v1",
        "a564d3288ba35913172894c3ccc47db6f90d68d7",
        "docs/analysis/pnjl/phase_reference_current_state_freeze_v1",
    ),
    (
        "phase_reference_limited_evidence_audit_v1",
        "docs/analysis/pnjl/phase_reference/phase_reference_limited_evidence_audit_v1",
        "a564d3288ba35913172894c3ccc47db6f90d68d7",
        "docs/analysis/pnjl/phase_reference_limited_evidence_audit_v1",
    ),
    (
        "phase_reference_manual_overlay_promotion_audit_v1",
        "docs/analysis/pnjl/phase_reference/phase_reference_manual_overlay_promotion_audit_v1",
        "a564d3288ba35913172894c3ccc47db6f90d68d7",
        "docs/analysis/pnjl/phase_reference_manual_overlay_promotion_audit_v1",
    ),
    (
        "phase_guided_transport_p128_xi001_analysis",
        "docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_p128_xi001_analysis",
        "5f75efab086bccac95eb07e84b7c5c15c8f011ba",
        "docs/analysis/relaxtime/phase_guided_transport_p128_xi001_analysis",
    ),
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def run_git(root: Path, *args: str) -> str:
    result = subprocess.run(
        ["git", *args],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def git_blob_bytes(root: Path, commit: str, relative_path: str) -> bytes | None:
    try:
        return subprocess.check_output(
            ["git", "cat-file", "blob", f"{commit}:{relative_path}"],
            cwd=root,
        )
    except subprocess.CalledProcessError:
        return None


def normalize_manifest_path(raw: str, package_rel: PurePosixPath) -> str:
    normalized = raw.replace("\\", "/").lstrip("./")
    prefixes = [package_rel.as_posix().rstrip("/")]
    # The transport package retained its pre-migration root in the manifest.
    if len(package_rel.parts) >= 2:
        prefixes.append(
            (package_rel.parent.parent / package_rel.name).as_posix().rstrip("/")
        )
    for prefix in prefixes:
        prefix = prefix + "/"
        if normalized.startswith(prefix):
            return normalized[len(prefix) :]
    return normalized


def metadata_record(
    path: Path, relative_to: Path, *, include_text: bool = False
) -> dict[str, Any]:
    relative = path.relative_to(relative_to).as_posix()
    if not path.is_file():
        return {"path": relative, "exists": False}
    raw = path.read_bytes()
    record = {
        "path": relative,
        "exists": True,
        "bytes": len(raw),
        "sha256": sha256_bytes(raw),
    }
    if include_text:
        record["content_utf8"] = raw.decode("utf-8")
    return record


def _normalized_text(data: bytes) -> str | None:
    try:
        return data.decode("utf-8").replace("\r\n", "\n").replace("\r", "\n")
    except UnicodeDecodeError:
        return None


def classify_mismatch(
    root: Path,
    package: Path,
    entry: dict[str, Any],
    inspected: dict[str, Any],
    baseline_commit: str,
    baseline_rel: PurePosixPath,
) -> dict[str, Any]:
    result = {
        "raw_path": entry.get("raw_path"),
        "path": entry.get("path"),
        "declared_sha256": entry.get("declared_sha256"),
        "actual_sha256": inspected.get("actual_sha256"),
    }
    if inspected.get("status") != "mismatch" or not inspected.get("exists"):
        result["classification"] = "author_check"
        result["reason"] = "the entry is not an existing hash mismatch"
        return result

    relative = PurePosixPath(str(entry["path"]))
    baseline_path = baseline_rel.joinpath(*relative.parts)
    baseline_data = git_blob_bytes(root, baseline_commit, baseline_path.as_posix())
    current_data = package.joinpath(*relative.parts).read_bytes()
    result["baseline"] = {
        "commit": baseline_commit,
        "path": baseline_path.as_posix(),
    }
    if baseline_data is None:
        result["classification"] = "author_check"
        result["reason"] = "historical baseline file is unavailable"
        return result

    result["baseline"]["bytes"] = len(baseline_data)
    result["baseline"]["sha256"] = sha256_bytes(baseline_data)
    if baseline_data == current_data:
        result["classification"] = "stale_metadata"
        result["reason"] = (
            "the current file is byte-identical to the first committed artifact, "
            "while the declared hash was already stale at that baseline"
        )
        return result

    baseline_text = _normalized_text(baseline_data)
    current_text = _normalized_text(current_data)
    if baseline_text is not None and current_text is not None:
        json_equivalent = False
        try:
            json_equivalent = json.loads(baseline_text) == json.loads(current_text)
        except json.JSONDecodeError:
            pass
        if baseline_text == current_text or json_equivalent:
            result["classification"] = "newline_or_serialization_drift"
            result["reason"] = (
                "the historical baseline differs only by text line endings or "
                "JSON serialization while retaining the same parsed content"
            )
            return result

    result["classification"] = "content_change"
    result["reason"] = (
        "the historical baseline and current file differ in content; do not "
        "force a metadata rewrite without author review"
    )
    return result


def load_declared_outputs(
    manifest: dict[str, Any], package_rel: PurePosixPath
) -> list[dict[str, Any]]:
    if isinstance(manifest.get("output_files"), dict):
        return [
            {
                "raw_path": str(path),
                "path": normalize_manifest_path(str(path), package_rel),
                "declared_sha256": str(value).lower(),
            }
            for path, value in manifest["output_files"].items()
        ]
    if isinstance(manifest.get("outputs"), list):
        entries = []
        for item in manifest["outputs"]:
            if not isinstance(item, dict) or "path" not in item:
                entries.append({"raw_entry": item, "invalid": True})
                continue
            entry = {
                "raw_path": str(item["path"]),
                "path": normalize_manifest_path(str(item["path"]), package_rel),
                "declared_sha256": str(item.get("sha256", "")).lower(),
            }
            if "bytes" in item:
                entry["declared_bytes"] = item["bytes"]
            entries.append(entry)
        return entries
    return []


def parse_checksums(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        return []
    entries = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split(maxsplit=1)
        if len(parts) != 2 or len(parts[0]) != 64:
            entries.append({"line": line_number, "raw": line, "invalid": True})
            continue
        entries.append(
            {
                "line": line_number,
                "path": parts[1].strip().replace("\\", "/"),
                "declared_sha256": parts[0].lower(),
            }
        )
    return entries


def audit_package(
    root: Path,
    name: str,
    package_rel_text: str,
    baseline_commit: str,
    baseline_rel_text: str,
) -> dict[str, Any]:
    package_rel = PurePosixPath(package_rel_text)
    baseline_rel = PurePosixPath(baseline_rel_text)
    package = root.joinpath(*package_rel.parts)
    manifest_path = package / "manifest.json"
    checksum_path = package / "checksums.sha256"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    output_entries = load_declared_outputs(manifest, package_rel)
    checksum_entries = parse_checksums(checksum_path)

    declared_paths = {
        entry["path"]
        for entry in output_entries
        if not entry.get("invalid") and "path" in entry
    }
    actual_files = sorted(
        path.relative_to(package).as_posix()
        for path in package.rglob("*")
        if path.is_file()
    )
    actual_output_files = [
        path for path in actual_files if path not in {"manifest.json", "checksums.sha256"}
    ]

    def inspect_entry(entry: dict[str, Any]) -> dict[str, Any]:
        if entry.get("invalid"):
            return dict(entry, status="invalid")
        relative = entry["path"]
        path = package.joinpath(*PurePosixPath(relative).parts)
        result = dict(entry)
        result["exists"] = path.is_file()
        if path.is_file():
            result["actual_bytes"] = path.stat().st_size
            result["actual_sha256"] = sha256_file(path)
            result["status"] = (
                "match" if result["actual_sha256"] == entry["declared_sha256"] else "mismatch"
            )
            if "declared_bytes" in entry:
                result["byte_count_status"] = (
                    "match" if result["actual_bytes"] == entry["declared_bytes"] else "mismatch"
                )
        else:
            result["status"] = "missing"
        return result

    output_audit = [inspect_entry(entry) for entry in output_entries]
    checksum_audit = [inspect_entry(entry) for entry in checksum_entries]
    expected_files = declared_paths | {
        entry["path"]
        for entry in checksum_entries
        if not entry.get("invalid") and "path" in entry
    }
    expected_files.add("manifest.json")
    if checksum_path.is_file():
        expected_files.add("checksums.sha256")
    extra_files = sorted(set(actual_files) - expected_files)
    missing_files = sorted(expected_files - set(actual_files))

    return {
        "name": name,
        "root": package_rel_text,
        "metadata_before": {
            "manifest": metadata_record(manifest_path, root, include_text=True),
            "checksums": metadata_record(checksum_path, root, include_text=True),
        },
        "historical_baseline": {
            "commit": baseline_commit,
            "root": baseline_rel_text,
        },
        "file_inventory": [
            metadata_record(package / relative, package) for relative in actual_files
        ],
        "manifest_output_audit": output_audit,
        "checksums_audit": checksum_audit,
        "declared_output_count": len(output_entries),
        "manifest_mismatch_count": sum(
            entry.get("status") == "mismatch" for entry in output_audit
        ),
        "manifest_missing_count": sum(
            entry.get("status") == "missing" for entry in output_audit
        ),
        "checksum_entry_count": len(checksum_entries),
        "checksum_mismatch_count": sum(
            entry.get("status") == "mismatch" for entry in checksum_audit
        ),
        "checksum_missing_count": sum(
            entry.get("status") == "missing" for entry in checksum_audit
        ),
        "extra_output_files": extra_files,
        "missing_output_files": missing_files,
        "mismatch_classifications": [
            classify_mismatch(
                root,
                package,
                entry,
                inspected,
                baseline_commit,
                baseline_rel,
            )
            for entry, inspected in [
                *zip(output_entries, output_audit),
                *zip(checksum_entries, checksum_audit),
            ]
            if inspected.get("status") == "mismatch"
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="write the JSON audit report")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="return non-zero when any declared hash/path/file-set check fails",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    packages = [
        audit_package(root, name, package_rel, baseline_commit, baseline_rel)
        for name, package_rel, baseline_commit, baseline_rel in PACKAGE_SPECS
    ]
    report: dict[str, Any] = {
        "schema": "docs_analysis_metadata_audit_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "repo": {
            "head": run_git(root, "rev-parse", "HEAD"),
            "branch": run_git(root, "branch", "--show-current"),
            "status_porcelain": run_git(root, "status", "--porcelain"),
        },
        "packages": packages,
    }
    report_text = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(report_text, encoding="utf-8", newline="\n")
    else:
        sys.stdout.write(report_text)

    failures = []
    for package in packages:
        if package["missing_output_files"] or package["extra_output_files"]:
            failures.append(package["name"] + ":file-set")
        if package["manifest_missing_count"] or package["checksum_missing_count"]:
            failures.append(package["name"] + ":missing")
        if package["manifest_mismatch_count"] or package["checksum_mismatch_count"]:
            failures.append(package["name"] + ":hash")
    if args.strict and failures:
        print("metadata audit failed: " + ", ".join(failures), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
