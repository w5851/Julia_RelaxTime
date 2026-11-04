"""Prepare the reviewed pilot requests and source archive; never dispatch jobs."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import zipfile


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def prepare(root, bound, output):
    binding = read(bound / "binding_manifest.json")
    commit, ref = binding["calculation_sha"], binding["workflow_ref"]
    require(not output.exists(), "use a new pilot evidence directory")
    require(subprocess.check_output(["git", "-C", str(root), "rev-parse", ref], text=True).strip() == commit,
            "local execution ref differs from the calculation commit")
    expected_ref = "refs/tags/" + ref
    remote = subprocess.check_output(["git", "-C", str(root), "ls-remote", "--tags", "origin", expected_ref], text=True).strip()
    require(remote == commit + "\t" + expected_ref, "remote calculation tag is missing or differs")
    verified_files = {}
    for name, checksum in binding["files"].items():
        path = (bound / name).resolve()
        require(path.is_relative_to(bound) and digest(path) == checksum, "bound input changed")
        verified_files[name] = checksum
    workflow = json.loads(subprocess.check_output(["gh", "api", "repos/w5851/Julia_RelaxTime/actions/workflows/291078197"]))
    require(workflow["state"] == "active" and workflow["path"] == binding["workflow"], "workflow is unavailable")
    output.mkdir(parents=True)
    for name in verified_files:
        shutil.copyfile(bound / name, output / name)
    shutil.copyfile(bound / "binding_manifest.json", output / "binding_manifest.json")
    archive = output / "calculation-source.zip"
    selected = ["Project.toml", "Manifest.toml", "src", "config/models/pnjl", "config/physics", "scripts/relaxtime",
                "scripts/analysis/relaxtime/prepare_transport_production_inputs.py",
                "benchmark/relaxtime/transport_runner_timing.jl"]
    subprocess.run(["git", "-C", str(root), "archive", "--format=zip", "--output=" + str(archive), commit, "--", *selected], check=True)
    source = output / "source"
    source.mkdir()
    with zipfile.ZipFile(archive) as bundle:
        for member in bundle.infolist():
            path = (source / member.filename).resolve()
            require(path.is_relative_to(source), "source archive path leaves its directory")
            if member.is_dir():
                path.mkdir(parents=True, exist_ok=True)
            else:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(bundle.read(member))
    source_manifest = dict(calculation_sha=commit, sha256=digest(archive), bytes=archive.stat().st_size,
                           selected_paths=selected, files=sum(p.is_file() for p in source.rglob("*")),
                           scope="Committed numerical sources and artifact-check helpers, without the large historical data archive")
    (output / "source_manifest.json").write_text(json.dumps(source_manifest, indent=2)+"\n", encoding="utf-8")
    files = {"A": "mode_a_runner_inputs.json", "B": "mode_b_runner_inputs.json"}
    planned = []
    for mode, filename in files.items():
        inputs = read(output / filename)
        require(inputs["pilot_input_ready"] and inputs["calculation_sha"] == commit and
                inputs["workflow_sha"] == commit and inputs["workflow_ref"] == ref,
                "input does not identify this execution ref")
        for job in inputs["pilot_jobs"]:
            require(job["inputs"]["calculation_sha"] == commit and job["inputs"]["run_tier"] == "production",
                    "pilot job differs from the calculation")
            request = dict(ref=ref, inputs={k: (str(v).lower() if isinstance(v, bool) else v)
                                           for k, v in job["inputs"].items()})
            name = job["label"] + ".request.json"
            (output / name).write_text(json.dumps(request, indent=2)+"\n", encoding="utf-8")
            planned.append(dict(mode=mode, job=job, request_file=name, request_sha256=digest(output / name)))
    manifest = dict(created_utc=datetime.now(timezone.utc).isoformat(), repository="w5851/Julia_RelaxTime",
                    calculation_sha=commit, workflow_sha=commit, workflow_ref=ref, workflow_id=291078197,
                    workflow_path=binding["workflow"], bound_inputs=files, bound_file_hashes=verified_files,
                    binding_manifest_sha256=digest(output / "binding_manifest.json"), source_archive_sha256=digest(archive),
                    remote_ref_verified=True, remote_ref=expected_ref, planned_runs=planned, runs=[],
                    controller_sha256=digest(Path(__file__)), full_dispatch_ready=False)
    (output / "dispatch_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n", encoding="utf-8")
    return dict(directory=str(output), calculation_sha=commit, workflow_ref=ref, pilot_jobs=len(planned),
                pilot_points=sum(entry["job"]["expected_points"] for entry in planned), actions_dispatched=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("bound_inputs", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    print(json.dumps(prepare(args.root.resolve(), args.bound_inputs.resolve(), args.output.resolve()), indent=2))
