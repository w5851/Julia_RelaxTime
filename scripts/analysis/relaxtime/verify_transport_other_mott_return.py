#!/usr/bin/env python3
"""Local archive/hash/metadata verification only; no physics calculation."""
from pathlib import Path
import argparse
import hashlib
import json
import stat
import zipfile
from verify_transport_eta_return import safe_name, sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("package", type=Path)
    args = parser.parse_args()
    root = args.package.resolve()
    returned = root / "returned"
    if returned.exists():
        raise FileExistsError(returned)
    receipt = json.loads((root / "archive_receipt.json").read_text(encoding="utf-8-sig"))
    archive = root.joinpath(*safe_name(receipt["archive"]).parts)
    assert archive.stat().st_size == receipt["bytes"] and sha(archive) == receipt["sha256"]
    with zipfile.ZipFile(archive) as bundle:
        entries = bundle.infolist()
        assert len(entries) == len({r.filename for r in entries})
        for entry in entries:
            safe_name(entry.filename)
            assert not entry.is_dir() and not stat.S_ISLNK(entry.external_attr >> 16)
        raw = bundle.read("package_manifest.json")
        assert hashlib.sha256(raw).hexdigest() == receipt["package_manifest_sha256"]
        manifest = json.loads(raw)
        expected = {r["path"]: r for r in manifest["files"]}
        assert len(expected) == receipt["manifest_files"]
        assert {r.filename for r in entries} == set(expected) | {"package_manifest.json"}
        for name, record in expected.items():
            content = bundle.read(name)
            assert len(content) == record["bytes"] and hashlib.sha256(content).hexdigest() == record["sha256"], name
        returned.mkdir()
        for entry in entries:
            path = returned.joinpath(*safe_name(entry.filename).parts)
            assert path.resolve().is_relative_to(returned)
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("xb") as stream:
                stream.write(bundle.read(entry.filename))
    for name, record in expected.items():
        path = returned.joinpath(*safe_name(name).parts)
        assert path.stat().st_size == record["bytes"] and sha(path) == record["sha256"]
    validation = json.loads((returned / "phase_validation.json").read_text(encoding="utf-8"))
    assert validation["status"] == "completed_diagnostic_phases"
    assert validation["source_files_checked"] == 340 and validation["all_source_hashes_match"]
    result = dict(status="verified_phase_return", files_verified=len(expected), archive=receipt,
                  completed_runs=validation["completed_runs"], whole_task_complete=validation["whole_task_complete"],
                  local_scope="archive, file identity and metadata only")
    (root / "local_return_verification.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
