#!/usr/bin/env python3
"""Verify and extract an eta diagnostic return; local metadata checks only."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path, PurePosixPath
import stat
import zipfile


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def safe_name(name):
    path = PurePosixPath(name)
    if not name or path.is_absolute() or ".." in path.parts or "\\" in name or ":" in name:
        raise ValueError(f"Unsafe archive path: {name}")
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    args = parser.parse_args()
    case = args.case.resolve()
    returned = case/"returned"
    receipt_path = case/"local_return_verification.json"
    if returned.exists() or receipt_path.exists():
        raise FileExistsError("Return or verification already exists; retain existing evidence")
    transfer = read(case/"archive_receipt.json")
    archive = case.joinpath(*safe_name(transfer["archive"]).parts)
    assert archive.stat().st_size == transfer["bytes"]
    assert sha(archive) == transfer["sha256"]
    with zipfile.ZipFile(archive) as bundle:
        entries = bundle.infolist()
        assert len(entries) == len({e.filename for e in entries}), "Duplicate archive entries"
        for entry in entries:
            safe_name(entry.filename)
            assert not entry.is_dir() and not stat.S_ISLNK(entry.external_attr >> 16)
        manifest_bytes = bundle.read("package_manifest.json")
        assert hashlib.sha256(manifest_bytes).hexdigest() == transfer["package_manifest_sha256"]
        manifest = json.loads(manifest_bytes)
        expected = {r["path"]: r for r in manifest["files"]}
        assert len(expected) == len(manifest["files"]) == transfer["manifest_files"]
        assert {e.filename for e in entries} == set(expected)|{"package_manifest.json"}
        for name, record in expected.items():
            data = bundle.read(name)
            assert len(data) == record["bytes"], name
            assert hashlib.sha256(data).hexdigest() == record["sha256"], name
        returned.mkdir()
        for entry in entries:
            target = returned.joinpath(*safe_name(entry.filename).parts)
            assert target.resolve().is_relative_to(returned)
            target.parent.mkdir(parents=True, exist_ok=True)
            with target.open("xb") as stream:
                stream.write(bundle.read(entry.filename))
    for name, record in expected.items():
        target = returned.joinpath(*safe_name(name).parts)
        assert target.stat().st_size == record["bytes"] and sha(target) == record["sha256"], name
    assert sha(returned/"package_manifest.json") == transfer["package_manifest_sha256"]
    validation = read(returned/"validation.json")
    assert validation["status"] == "completed_diagnostic"
    assert validation["calculation_sha"] == "2e5c411c716d5cbc524c3fc6cfb335fe05b39444"
    assert validation["points"] == 51 and validation["channels"] == 867
    assert validation["comparison_points"] == 18
    assert validation["source_files_checked"] == 340 and validation["all_source_hashes_match"]
    for name in ("locate_v3.status.json", "transport_queue_v2.status.json", "peak_refinement_v1.status.json"):
        stage = read(returned/name)
        assert stage["status"] == "finished", name
        assert all(w["exit_code"] == 0 for w in stage["workers"]), name
    finish_path = case/"finish_eta_v1.status.json"
    assert read(finish_path)["status"] == "finished"
    review = read(returned/"review_final_v1/review.json")
    assert review["status"] == "completed_diagnostic_review"
    assert sha(returned/"review_final_v1/review.json") == validation["review_sha256"]
    for name, digest in review["output_sha256"].items():
        assert sha(returned/"review_final_v1"/name) == digest, name
    receipt = dict(status="verified_return",received_utc=datetime.now(timezone.utc).isoformat(),
        archive=transfer,files_verified=len(expected),
        final_status_sha256=sha(finish_path),validation_sha256=sha(returned/"validation.json"),
        final_review_sha256=sha(returned/"review_final_v1/review.json"),
        fine_points=validation["points"],reference_evaluations=validation["comparison_points"],
        local_scope="Archive integrity, metadata and source receipts; no local numerical analysis or solve",
        scientific_convergence_review="separate")
    receipt_path.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(receipt,indent=2))


if __name__ == "__main__":
    main()
