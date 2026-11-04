"""Verify the retained pilot bundle, review on W, and return a batch plan."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import traceback
import zipfile


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2)
        stream.write("\n")


def child(root, name):
    path = (root / name).resolve()
    require(path != root and path.is_relative_to(root), "path leaves the review case")
    return path


def verify_files(root, files):
    for name, checksum in files.items():
        require(digest(child(root, name)) == checksum, "input changed: " + name)


def unpack(root, archive, manifest):
    require(digest(archive) == manifest["sha256"], "runner bundle checksum differs")
    require(archive.stat().st_size == manifest["bytes"], "runner bundle size differs")
    with zipfile.ZipFile(archive) as bundle:
        members = [m for m in bundle.infolist() if not m.is_dir()]
        names = [m.filename for m in members]
        destinations = [child(root, name) for name in names]
        require(len(set(destinations)) == len(names), "duplicate runner bundle destinations")
        require(set(names) == set(manifest["files"]), "runner bundle file list differs")
        for member, path in zip(members, destinations):
            require(not stat.S_ISLNK(member.external_attr >> 16), "runner bundle symlink")
            data = bundle.read(member)
            require(hashlib.sha256(data).hexdigest() == manifest["files"][member.filename],
                    "runner bundle member differs: " + member.filename)
            if path.exists():
                require(path.is_file() and path.read_bytes() == data,
                        "retained runner evidence differs: " + member.filename)
            else:
                path.parent.mkdir(parents=True, exist_ok=True)
                with path.open("xb") as stream:
                    stream.write(data)
    return len(members)


def run(root):
    require(os.environ.get("COMPUTERNAME", "").upper() == "W", "review must run on quote-remote W")
    require(root.is_relative_to(Path(r"D:\CodexCompute\Julia_RelaxTime").resolve()),
            "unexpected remote workspace")
    execution = root / "execution_v2"
    execution.mkdir()
    status = dict(host=os.environ["COMPUTERNAME"], passed=False,
                  calculation_sha=None, reviewer_exit_code=None, planner_exit_code=None,
                  started_utc=datetime.now(timezone.utc).isoformat(),
                  controller_sha256=digest(Path(__file__)), numerical_calculation_performed=False)
    try:
        contract = read(root / "pilot_review_execution_manifest_v2.json")
        verify_files(root, contract["files"])
        foundation = read(root / "review_foundation_manifest.json")
        require(digest(root / "review_foundation.zip") == foundation["sha256"] ==
                contract["foundation_sha256"], "foundation changed")
        control = read(root / "review_control_manifest.json")
        verify_files(root, control["files"])
        status["calculation_sha"] = control["calculation_sha"]
        count = unpack(root, root / "runner_results.zip", read(root / "runner_results_manifest.json"))
        require(read(root / "collection_status.json")["ready_for_remote_review"],
                "runner collection is incomplete")
        module_spec = importlib.util.spec_from_file_location("frozen_pilot_review", root / "review_transport_runner_pilot_v2.py")
        reviewer = importlib.util.module_from_spec(module_spec)
        module_spec.loader.exec_module(reviewer)
        source = read(root / "source_manifest.json")
        reviewer.verified_zip(root / "calculation-source.zip", root / "source", source["sha256"])
        write(execution / "input_receipt.json", dict(inputs_verified=True, source_files=source["files"],
              bundle_files=count, execution_manifest_sha256=digest(root / "pilot_review_execution_manifest_v2.json"),
              runner_bundle_sha256=digest(root / "runner_results.zip")))
        for name, script, arguments in (
                ("tests", "reviewer_check_v2/tests/unit/python/test_transport_runner_pilot.py", []),
                ("reviewer", "review_transport_runner_pilot_v2.py", [str(root)]),
                ("planner", "plan_remaining_batches.py", [str(root)])):
            with (execution / (name + ".stdout.log")).open("xb") as stdout, \
                 (execution / (name + ".stderr.log")).open("xb") as stderr:
                result = subprocess.run([sys.executable, "-B", str(root / script), *arguments],
                                        cwd=root, stdout=stdout, stderr=stderr, check=False)
            status[name + "_exit_code"] = result.returncode
            require(result.returncode == 0, name + " failed; inspect retained logs")
        verify_files(root, contract["files"])
        verify_files(root, control["files"])
        reviewer.verified_zip(root / "calculation-source.zip", root / "source", source["sha256"])
        require("Ran 7 tests" in (execution / "tests.stderr.log").read_text(encoding="utf-8"),
                "reviewer test count differs")
        status["reviewer_tests_passed"] = 7
        status["passed"] = True
        status["inputs_unchanged_after_review"] = True
    except Exception:
        status["error"] = traceback.format_exc()
    status["finished_utc"] = datetime.now(timezone.utc).isoformat()
    write(execution / "remote_status.json", status)
    selected = [p for p in execution.iterdir() if p.is_file()]
    selected += [root / name for name in ("pilot_review_execution_manifest_v2.json",
                 "runner_results_manifest.json", "plan_remaining_batches.py", "run_remote_pilot_review_v2.py")]
    if (root / "reviewed_pilot.json").is_file():
        selected.append(root / "reviewed_pilot.json")
    if (root / "remaining_batches").is_dir():
        selected += [p for p in (root / "remaining_batches").rglob("*") if p.is_file()]
    selected += [p for p in (root / "execution_v1").iterdir() if p.is_file()]
    selected += [p for p in (root / "reviewer_check_v2").rglob("*") if p.is_file()]
    selected += [root / name for name in ("review_transport_runner_pilot_v2.py", "revision_manifest_v2.json",
                 "runner-pilot-review-return-v1.zip", "runner-pilot-review-return-v1.manifest.json")]
    selected = sorted(set(p for p in selected if p.is_file()))
    files = {p.relative_to(root).as_posix(): dict(sha256=digest(p), bytes=p.stat().st_size) for p in selected}
    manifest = execution / "package_manifest.json"
    write(manifest, dict(files=files, passed=status["passed"]))
    selected.append(manifest)
    archive = root / "runner-pilot-review-return-v2.zip"
    with zipfile.ZipFile(archive, "x", compression=zipfile.ZIP_DEFLATED) as bundle:
        for path in selected:
            bundle.write(path, path.relative_to(root).as_posix())
    receipt = dict(passed=status["passed"], sha256=digest(archive), bytes=archive.stat().st_size,
                   files=len(selected), remote_archive=str(archive))
    write(root / "runner-pilot-review-return-v2.manifest.json", receipt)
    print(json.dumps(receipt, indent=2))
    return 0 if status["passed"] else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    sys.exit(run(args.directory.resolve()))
