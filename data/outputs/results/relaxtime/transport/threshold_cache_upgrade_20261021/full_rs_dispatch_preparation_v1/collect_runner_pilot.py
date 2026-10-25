"""Download completed pilot evidence with archive checks only; no numerical review."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import stat
import subprocess
import zipfile


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def api(endpoint):
    return json.loads(subprocess.check_output(["gh", "api", endpoint]))


def retain(path, value):
    if path.exists():
        return read(path)
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2)
        stream.write("\n")
    return value


def collect(directory):
    specification = read(directory / "dispatch_manifest.json")
    repository = specification["repository"]
    results = []
    for entry in specification["runs"]:
        run_id = entry["run_id"]
        run = api(f"repos/{repository}/actions/runs/{run_id}")
        require(run["head_sha"] == specification["workflow_sha"] and
                run["run_attempt"] == entry.get("run_attempt", 1), "run identity changed")
        record = dict(run_id=run_id, label=entry["job"]["label"], status=run["status"], conclusion=run["conclusion"])
        results.append(record)
        if run["status"] != "completed":
            continue
        folder = (directory / entry["directory"]).resolve()
        require(folder != directory and folder.is_relative_to(directory), "run directory leaves the case")
        folder.mkdir(parents=True, exist_ok=True)
        retain(folder / "run.json", run)
        retain(folder / "jobs.json", api(f"repos/{repository}/actions/runs/{run_id}/jobs?per_page=100"))
        artifacts = api(f"repos/{repository}/actions/runs/{run_id}/artifacts?per_page=100")
        expected_name = f"transport-{entry['job']['inputs']['case_name']}-{entry['job']['label']}"
        selected = [a for a in artifacts["artifacts"] if a["name"] == expected_name]
        require(len(selected) <= 1, "duplicate pilot artifacts")
        if not selected:
            record["artifact_available"] = False
            continue
        artifact = retain(folder / "artifact.json", selected[0])
        require(not artifact["expired"] and artifact["workflow_run"]["id"] == run_id,
                "artifact expired or belongs to another run")
        archive = folder / "artifact.zip"
        if not archive.exists():
            partial = folder / "artifact.download"
            require(not partial.exists(), "a partial download remains; inspect it before retrying")
            with partial.open("xb") as stream:
                subprocess.run(["gh", "api", f"repos/{repository}/actions/artifacts/{artifact['id']}/zip"],
                               stdout=stream, check=True)
            partial.rename(archive)
        checksum = digest(archive)
        if artifact.get("digest"):
            require(artifact["digest"] == "sha256:" + checksum, "GitHub artifact digest mismatch")
        unpacked = folder / "artifact"
        with zipfile.ZipFile(archive) as bundle:
            members = [m for m in bundle.infolist() if not m.is_dir()]
            destinations = [(unpacked / m.filename).resolve() for m in members]
            require(len({os.path.normcase(str(p)) for p in destinations}) == len(members),
                    "duplicate archive destinations")
            require(all(p != unpacked and p.is_relative_to(unpacked) for p in destinations),
                    "archive member leaves its directory")
            for member, path in zip(members, destinations):
                require(not stat.S_ISLNK(member.external_attr >> 16), "archive symlink")
                data = bundle.read(member)
                if path.exists():
                    require(path.is_file() and path.read_bytes() == data, "existing extracted file differs")
                else:
                    path.parent.mkdir(parents=True, exist_ok=True)
                    with path.open("xb") as stream:
                        stream.write(data)
            require({p.resolve() for p in unpacked.rglob("*") if p.is_file()} == set(destinations),
                    "unexpected files in extracted artifact")
        manifest = dict(sha256=checksum, bytes=archive.stat().st_size, files=len(members), artifact_id=artifact["id"],
                        github_digest=artifact.get("digest"), downloaded_utc=datetime.now(timezone.utc).isoformat())
        previous = retain(folder / "archive_manifest.json", manifest)
        require(all(previous[key] == manifest[key] for key in ("sha256", "bytes", "files", "artifact_id")),
                "retained artifact identity changed")
        record["artifact_available"] = True
        record["archive_sha256"] = checksum
    ready = (len(results) == len(specification["planned_runs"]) and bool(results) and
             all(row["status"] == "completed" and row["conclusion"] == "success" and
                 row.get("artifact_available") for row in results))
    status = dict(runs=results, ready_for_remote_review=ready, numerical_review_performed=False,
                  collector_sha256=digest(Path(__file__)), updated_utc=datetime.now(timezone.utc).isoformat())
    (directory / "collection_status.json").write_text(json.dumps(status, indent=2)+"\n", encoding="utf-8")
    return status


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    print(json.dumps(collect(args.directory.resolve()), indent=2))
