"""Verify or restore the retained chemical-freezeout review evidence.

Raw point records, source snapshots and search responses are stored in ZIP archives with
repository-relative member paths. Existing files with different bytes are
rejected before restoration; frozen records are never rewritten.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import zipfile

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MANIFEST = Path(
    "data/outputs/results/relaxtime/meson_density/freezeout_review_20261022/archive_manifest.json"
)


def safe_path(root: Path, name: str) -> Path:
    relative = PurePosixPath(name)
    if not name or relative.is_absolute() or ".." in relative.parts or "\\" in name or ":" in name:
        raise ValueError(f"invalid repository-relative evidence path: {name}")
    target = (root / relative).resolve()
    if target == root or not target.is_relative_to(root):
        raise ValueError(f"evidence path escapes repository: {name}")
    return target


def check_record(path: Path, record: dict) -> None:
    data = path.read_bytes()
    if len(data) != record["bytes"] or hashlib.sha256(data).hexdigest() != record["sha256"]:
        raise ValueError(f"evidence bytes or SHA-256 mismatch: {path}")


def verify(manifest_path: Path, root: Path) -> tuple[dict, dict]:
    root = root.resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != "freezeout_review_evidence_v1":
        raise ValueError("unsupported freezeout evidence manifest schema")
    for record in manifest["retained_files"]:
        check_record(safe_path(root, record["path"]), record)
    members = set()
    for record in manifest["archives"]:
        path = safe_path(root, record["path"])
        check_record(path, record)
        with zipfile.ZipFile(path) as archive:
            entries = archive.infolist()
            if len(entries) != record["files"]:
                raise ValueError(f"archive member count mismatch: {path}")
            if sum(entry.file_size for entry in entries) != record["uncompressed_bytes"]:
                raise ValueError(f"archive uncompressed size mismatch: {path}")
            for entry in entries:
                target = safe_path(root, entry.filename)
                if entry.is_dir() or target in members:
                    raise ValueError(f"duplicate or directory archive member: {entry.filename}")
                members.add(target)
            if archive.testzip() is not None:
                raise ValueError(f"archive CRC mismatch: {path}")
    result = dict(retained_files_checked=len(manifest["retained_files"]),
                  archives_checked=len(manifest["archives"]), archive_members_checked=len(members))
    return manifest, result


def restore(manifest_path: Path, root: Path) -> dict:
    root = root.resolve()
    manifest, result = verify(manifest_path, root)
    # Check every conflict before writing any restored file.
    for record in manifest["archives"]:
        with zipfile.ZipFile(safe_path(root, record["path"])) as archive:
            for entry in archive.infolist():
                target = safe_path(root, entry.filename)
                if target.exists() and (not target.is_file() or target.read_bytes() != archive.read(entry)):
                    raise ValueError(f"restore would replace existing evidence: {target}")
    written = 0
    for record in manifest["archives"]:
        with zipfile.ZipFile(safe_path(root, record["path"])) as archive:
            for entry in archive.infolist():
                target = safe_path(root, entry.filename)
                if target.exists():
                    continue
                target.parent.mkdir(parents=True, exist_ok=True)
                with target.open("xb") as handle:
                    handle.write(archive.read(entry))
                written += 1
    return dict(result, restored_files=written)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("verify", "restore"))
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    args = parser.parse_args()
    root = args.root.resolve()
    manifest = args.manifest if args.manifest.is_absolute() else root / args.manifest
    result = restore(manifest, root) if args.action == "restore" else verify(manifest, root)[1]
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
