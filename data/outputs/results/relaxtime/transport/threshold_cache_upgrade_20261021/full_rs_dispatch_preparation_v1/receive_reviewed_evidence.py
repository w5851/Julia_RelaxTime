"""Receive this batch's reviewed archives using byte and metadata checks only."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import stat
import zipfile


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def checked_path(base, relative):
    base = base.resolve()
    path = (base / relative.replace("\\", "/")).resolve()
    if not path.is_relative_to(base) or path == base:
        raise ValueError(f"path leaves the received directory: {relative}")
    return path


def receive(folder):
    archive = folder / "returned_evidence.zip"
    declared = read(folder / "returned_archive_manifest.json")
    if digest(archive) != declared["sha256"] or archive.stat().st_size != declared["bytes"]:
        raise ValueError("received archive identity differs from the remote declaration")
    target = folder / "returned"
    if not target.resolve().is_relative_to(folder.resolve()):
        raise ValueError("received directory resolves outside this case")
    with zipfile.ZipFile(archive) as bundle:
        members = bundle.infolist()
        destinations = [checked_path(target, member.filename) for member in members]
        if (len(members) != declared["files"] or
                len({os.path.normcase(str(p)) for p in destinations}) != len(members) or
                any(m.is_dir() or stat.S_ISLNK(m.external_attr >> 16) for m in members)):
            raise ValueError("invalid or duplicate archive members")
        for member, path in zip(members, destinations):
            content = bundle.read(member)
            if path.exists():
                if not path.is_file() or path.read_bytes() != content:
                    raise ValueError(f"existing received file differs: {member.filename}")
            else:
                path.parent.mkdir(parents=True, exist_ok=True)
                with path.open("xb") as stream:
                    stream.write(content)
        actual = {os.path.normcase(str(p.resolve())) for p in target.rglob("*") if p.is_file()}
        if actual != {os.path.normcase(str(p)) for p in destinations}:
            raise ValueError("unexpected files in the received directory")
    # These are pre-run local identities, not values reconstructed from outputs.
    copied = {}
    for name in ("package_manifest.json", "control_manifest.json", "input_manifest.json"):
        before, after = folder / name, target / name
        if before.is_file() and after.is_file():
            if before.read_bytes() != after.read_bytes():
                raise ValueError(f"pre-run manifest changed: {name}")
            copied[name] = digest(after)
    control_path = folder / "control_manifest.json"
    if control_path.is_file():
        for relative, checksum in read(control_path)["files"].items():
            before = checked_path(folder, relative)
            if digest(before) != checksum:
                raise ValueError(f"frozen local control changed: {relative}")
            after = checked_path(target, relative)
            if after.is_file() and digest(after) != checksum:
                raise ValueError(f"returned control changed: {relative}")
    source_count = input_count = 0
    # Match the numerical runs reviewed by this batch; copied input manifests
    # are evidence files, not additional producers with their own source tree.
    for stage in ("cheap", "kernel", "outer", "joint", "production"):
        path = target / stage / "manifest.json"
        if not path.is_file():
            continue
        manifest = read(path)
        sources = manifest.get("source_files")
        if not isinstance(sources, list):
            continue
        for item in sources:
            frozen = checked_path(path.parent / "source", item["path"])
            if digest(frozen) != item["sha256"]:
                raise ValueError(f"source snapshot changed: {frozen}")
            source_count += 1
        for item in manifest.get("input_files", []):
            frozen = checked_path(path.parent, item["snapshot"])
            if digest(frozen) != item["sha256"]:
                raise ValueError(f"input snapshot changed: {frozen}")
            input_count += 1
        identity = manifest.get("source_archive")
        if identity:
            frozen = checked_path(path.parent, identity["package_manifest_snapshot"])
            if digest(frozen) != identity["package_manifest_sha256"]:
                raise ValueError("source package manifest snapshot changed")
    review_path = target / "reviewed_result.json"
    review = read(review_path)
    if "source_snapshots_verified" in review and source_count != review["source_snapshots_verified"]:
        raise ValueError("source snapshot count differs from the remote review")
    reviewers = [p.name for p in target.glob("*.py") if digest(p) == review["reviewer_sha256"]]
    if len(reviewers) != 1:
        raise ValueError("the actual remote reviewer is missing or ambiguous")
    receipt = dict(received_utc=datetime.now(timezone.utc).isoformat(), archive=declared,
                   extracted_files_verified=declared["files"],
                   pre_run_manifests_verified=copied,
                   source_snapshots_verified=source_count, input_snapshots_verified=input_count,
                   remote_review_sha256=digest(review_path), remote_reviewer=reviewers[0],
                   receiver_sha256=digest(Path(__file__)),
                   check_scope="Archive bytes, frozen snapshots and metadata only; numerical review ran on quote-remote")
    receipt_path = folder / "local_receipt.json"
    if receipt_path.exists():
        previous = read(receipt_path)
        for name in ("archive", "remote_review_sha256", "source_snapshots_verified", "input_snapshots_verified"):
            if previous[name] != receipt[name]:
                raise ValueError("existing receipt differs; preserve it and investigate")
        return previous
    with receipt_path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(receipt, stream, indent=2)
        stream.write("\n")
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_directory", type=Path)
    args = parser.parse_args()
    print(json.dumps(receive(args.case_directory.resolve()), indent=2))
