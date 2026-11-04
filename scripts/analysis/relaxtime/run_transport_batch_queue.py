"""Resume a bounded Actions queue and collect immutable transport artifacts.

This controller performs network and file-identity operations only. Numerical
acceptance is a separate operation on the compute host. A saved dispatch attempt
is never automatically posted twice, even when the first HTTP reply was lost.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import time
import zipfile


def require(condition, message):
    if not condition:
        raise ValueError(message)


def now():
    return datetime.now(timezone.utc).isoformat()


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def retain(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        return read(path)
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(data, stream, indent=2, ensure_ascii=False)
        stream.write("\n")
    return data


def child(root, name):
    root = root.resolve()
    target = (root / name).resolve()
    require(target != root and target.is_relative_to(root), "path leaves case directory")
    return target


def api(endpoint):
    # Only reads can be retried automatically.
    for attempt in range(3):
        try:
            return json.loads(subprocess.check_output(["gh", "api", endpoint], timeout=60))
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
            if attempt == 2:
                raise
            time.sleep(2 ** attempt)


def title(entry):
    return f"transport {entry['job']['inputs']['case_name']} / {entry['job']['label']} / production"


def list_runs(spec, get=api):
    route = f"repos/{spec['repository']}/actions/workflows/{spec['workflow_id']}/runs"
    records = {}
    page = 1
    while True:
        data = get(route + f"?event=workflow_dispatch&head_sha={spec['workflow_sha']}&per_page=100&page={page}")
        for run in data["workflow_runs"]:
            records[run["id"]] = run
        if len(data["workflow_runs"]) < 100:
            break
        page += 1
    return list(records.values())


def recover_entry(entry, spec, runs):
    matches = [run for run in runs if run["display_title"] == title(entry)]
    require(len(matches) <= 1, "duplicate workflow runs for " + entry["job"]["label"])
    if not matches:
        return None
    run = matches[0]
    require(run["head_sha"] == spec["workflow_sha"] and run["head_branch"] == spec["workflow_ref"]
            and run["event"] == "workflow_dispatch" and run["run_attempt"] == 1,
            "workflow identity or attempt differs from frozen plan")
    return dict(entry, run_id=run["id"], run_attempt=run["run_attempt"], run_url=run["html_url"],
                directory=f"runs/{run['id']}")


def dispatch(directory, entry, spec):
    request = child(directory, entry["request_file"])
    require(digest(request) == entry["request_sha256"], "request bytes changed")
    payload = read(request)
    # Prepared GitHub API requests encode boolean inputs as lowercase strings.
    inputs = {key: str(value).lower() if isinstance(value, bool) else value
              for key, value in entry["job"]["inputs"].items()}
    require(payload == dict(ref=spec["workflow_ref"], inputs=inputs),
            "request differs from frozen job")
    path = directory / "attempts" / (entry["job"]["label"] + ".json")
    require(not path.exists(), "dispatch attempt already exists; recover it instead of posting again")
    attempt = dict(label=entry["job"]["label"], request_sha256=entry["request_sha256"],
                   started_utc=now(), request_may_have_been_accepted=True)
    retain(path, attempt)
    print("DISPATCH " + entry["job"]["label"], flush=True)
    route = f"repos/{spec['repository']}/actions/workflows/{spec['workflow_id']}/dispatches"
    try:
        result = subprocess.run(["gh", "api", "--method", "POST", route, "--input", str(request)],
                                capture_output=True, encoding="utf-8", errors="replace", timeout=60)
        attempt.update(returncode=result.returncode, stdout=result.stdout, stderr=result.stderr, response_utc=now())
    except subprocess.TimeoutExpired:
        attempt.update(response_utc=now(), transport_timeout=True)
    write(path, attempt)
    # Recovery handles both an HTTP success and an uncertain network outcome.
    return attempt


def extract_verified(archive, unpacked):
    unpacked = unpacked.resolve()
    with zipfile.ZipFile(archive) as bundle:
        members = [member for member in bundle.infolist() if not member.is_dir()]
        destinations = [child(unpacked, member.filename) for member in members]
        require(len({os.path.normcase(str(p)) for p in destinations}) == len(members),
                "duplicate archive destinations")
        for member, destination in zip(members, destinations):
            require(not stat.S_ISLNK(member.external_attr >> 16), "archive symlink")
            data = bundle.read(member)
            if destination.exists():
                require(destination.is_file() and destination.read_bytes() == data,
                        "extracted artifact differs from archive")
            else:
                destination.parent.mkdir(parents=True, exist_ok=True)
                with destination.open("xb") as stream:
                    stream.write(data)
        require({p.resolve() for p in unpacked.rglob("*") if p.is_file()} == set(destinations),
                "unexpected extracted files")
    return len(members)


def collect_run(directory, entry, spec, run):
    folder = child(directory, entry["directory"])
    require(run["id"] == entry["run_id"] and run["status"] == "completed" and
            run["head_sha"] == spec["workflow_sha"] and run["run_attempt"] == entry["run_attempt"],
            "collection run identity changed")
    retain(folder / "run.json", run)
    repository, run_id = spec["repository"], entry["run_id"]
    if not (folder / "jobs.json").exists():
        retain(folder / "jobs.json", api(f"repos/{repository}/actions/runs/{run_id}/jobs?per_page=100"))
    if (folder / "archive_manifest.json").exists():
        archive = read(folder / "archive_manifest.json")
        require(digest(folder / "artifact.zip") == archive["sha256"], "retained archive changed")
        return dict(artifact_available=True, archive_sha256=archive["sha256"])
    if (folder / "artifact.json").exists():
        artifact = read(folder / "artifact.json")
    else:
        artifacts = api(f"repos/{repository}/actions/runs/{run_id}/artifacts?per_page=100")
        expected = f"transport-{entry['job']['inputs']['case_name']}-{entry['job']['label']}"
        selected = [artifact for artifact in artifacts["artifacts"] if artifact["name"] == expected]
        require(len(selected) <= 1, "duplicate artifacts")
        if not selected:
            retain(folder / "no_artifact_observed.json", dict(observed_utc=now(), artifacts=artifacts))
            return dict(artifact_available=False)
        artifact = retain(folder / "artifact.json", selected[0])
    require(not artifact["expired"] and artifact["workflow_run"]["id"] == run_id and
            artifact["workflow_run"]["head_sha"] == spec["workflow_sha"], "artifact identity differs")
    archive = folder / "artifact.zip"
    if not archive.exists():
        # Interrupted downloads are retained; a new attempt never overwrites them.
        number = 1
        while (folder / f"artifact.download.{number}").exists():
            number += 1
        partial = folder / f"artifact.download.{number}"
        with partial.open("xb") as stream:
            subprocess.run(["gh", "api", f"repos/{repository}/actions/artifacts/{artifact['id']}/zip"],
                           stdout=stream, check=True, timeout=180)
        require(artifact.get("digest") == "sha256:" + digest(partial), "GitHub artifact digest mismatch")
        partial.rename(archive)
    checksum = digest(archive)
    require(artifact.get("digest") == "sha256:" + checksum, "GitHub archive digest differs")
    files = extract_verified(archive, folder / "artifact")
    retain(folder / "archive_manifest.json", dict(sha256=checksum, bytes=archive.stat().st_size,
           files=files, artifact_id=artifact["id"], github_digest=artifact["digest"], downloaded_utc=now()))
    print(f"COLLECT {run_id} {entry['job']['label']} {run['conclusion']}", flush=True)
    return dict(artifact_available=True, archive_sha256=checksum)


def prepare(pilot, directory):
    require(not directory.exists(), "new queue directory already exists")
    previous = read(pilot / "dispatch_manifest.json")
    batches = pilot / "returned/remaining_batches"
    manifest = read(batches / "batches_manifest.json")
    require(digest(pilot / "returned/reviewed_pilot.json") == manifest["review_sha256"], "pilot review changed")
    for name, checksum in manifest["files"].items():
        require(digest(child(batches, name)) == checksum, "batch bytes changed: " + name)
    require(previous["calculation_sha"] == previous["workflow_sha"] == manifest["calculation_sha"],
            "calculation identity mismatch")
    ref = api(f"repos/{previous['repository']}/git/ref/tags/{previous['workflow_ref']}")
    require(ref["object"]["type"] == "commit" and ref["object"]["sha"] == previous["workflow_sha"],
            "remote execution tag moved")
    grouped = {}
    for mode in ("A", "B"):
        plan = read(batches / f"mode_{mode.lower()}_batches.json")
        require(plan["full_dispatch_ready"] and plan["runner_cost_qualified"], "batch is not ready")
        require(all(plan[k] == previous[k] for k in ("calculation_sha", "workflow_sha", "workflow_ref")),
                "batch identity differs")
        grouped[mode] = [dict(mode=mode, job=job, request_file="requests/" + job["label"] + ".request.json",
                       request_sha256=manifest["files"][job["label"] + ".request.json"])
                       for job in plan["jobs"]]
    # Interleave modes; the cap limits outstanding jobs, not scientific coverage.
    planned = [entries[index] for index in range(max(map(len, grouped.values())))
               for entries in grouped.values() if index < len(entries)]
    require(len(planned) == 135 and sum(e["job"]["expected_points"] for e in planned) == 1808,
            "full-grid batch count changed")
    require(len({title(entry) for entry in planned}) == len(planned), "duplicate planned job titles")
    directory.mkdir(parents=True)
    (directory / "requests").mkdir()
    (directory / "inputs").mkdir()
    for name in manifest["files"]:
        target = "requests" if name.endswith(".request.json") else "inputs"
        shutil.copyfile(batches / name, directory / target / name)
    for original, relative in ((batches / "batches_manifest.json", "inputs/batches_manifest.json"),
        (pilot / "dispatch_manifest.json", "inputs/pilot_dispatch_manifest.json"),
        (pilot / "returned/reviewed_pilot.json", "inputs/reviewed_pilot.json"),
        (pilot / "source_manifest.json", "source_manifest.json"),
        (pilot / "calculation-source.zip", "calculation-source.zip"),
        (pilot / "mode_a_runner_inputs.json", "inputs/mode_a_runner_inputs.json"),
        (pilot / "mode_b_runner_inputs.json", "inputs/mode_b_runner_inputs.json")):
        shutil.copyfile(original, directory / relative)
    for entry in previous["runs"]:
        shutil.copytree(child(pilot, entry["directory"]), directory / "pilot_runs" / str(entry["run_id"]))
    spec = {key: previous[key] for key in ("repository", "calculation_sha", "workflow_sha", "workflow_ref",
                                         "workflow_id", "workflow_path", "source_archive_sha256")}
    spec.update(created_utc=now(), planned_runs=planned, runs=[], reused_pilot_runs=previous["runs"],
                max_outstanding=20, full_grid_points=1819, remaining_points=1808,
                full_length_pilot_is_prerequisite=False, remote_ref_verified=True,
                numerical_review_host="quote-remote", formal_production_promoted=False,
                controller_sha256=digest(Path(__file__)),
                authorization="User authorized all remaining batches, artifact retrieval and full-grid acceptance; no wait for a full-length pilot; supplemental numerical validation on quote-remote.")
    spec["frozen_files"] = {path.relative_to(directory).as_posix(): digest(path)
                            for path in directory.rglob("*") if path.is_file()}
    retain(directory / "queue_plan.json", spec)
    write(directory / "dispatch_manifest.json", spec)
    return spec


def cycle(directory, spec):
    runs = list_runs(spec)
    recovered = [found for entry in spec["planned_runs"] if (found := recover_entry(entry, spec, runs))]
    by_id = {run["id"]: run for run in runs}
    for old in spec["runs"]:
        require(old in recovered, "previously recorded run disappeared or changed")
    spec["runs"] = recovered
    write(directory / "dispatch_manifest.json", spec)
    results, errors = [], []
    for entry in recovered:
        run = by_id[entry["run_id"]]
        record = dict(run_id=run["id"], label=entry["job"]["label"], status=run["status"],
                      conclusion=run["conclusion"], expected_points=entry["job"]["expected_points"])
        if run["status"] == "completed":
            try:
                record.update(collect_run(directory, entry, spec, run))
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
                errors.append(dict(run_id=run["id"], error=str(error)))
        results.append(record)
    present = {entry["job"]["label"] for entry in recovered}
    outstanding = sum(run["status"] != "completed" for run in results)
    uncertain = []
    for entry in spec["planned_runs"]:
        label = entry["job"]["label"]
        if label in present:
            continue
        attempt_path = directory / "attempts" / (label + ".json")
        if attempt_path.exists():
            uncertain.append(label)
            outstanding += 1
    posted = []
    for entry in spec["planned_runs"]:
        label = entry["job"]["label"]
        if label in present or label in uncertain or outstanding >= spec["max_outstanding"]:
            continue
        dispatch(directory, entry, spec)
        posted.append(label)
        outstanding += 1
        time.sleep(2)
    complete = len(results) == len(spec["planned_runs"]) and all(r["status"] == "completed" for r in results)
    ready = complete and all(r["conclusion"] == "success" and r.get("artifact_available") for r in results)
    status = dict(updated_utc=now(), total_jobs=len(spec["planned_runs"]), recovered_jobs=len(results),
                  completed_jobs=sum(r["status"] == "completed" for r in results),
                  successful_jobs=sum(r["conclusion"] == "success" for r in results),
                  collected_jobs=sum(bool(r.get("artifact_available")) for r in results),
                  successful_points=sum(r["expected_points"] for r in results if r["conclusion"] == "success"),
                  newly_posted_labels=posted, unresolved_attempts=uncertain, runs=results, collection_errors=errors,
                  all_jobs_terminal=complete, ready_for_remote_review=ready, numerical_review_performed=False)
    write(directory / "collection_status.json", status)
    print(json.dumps({k: v for k, v in status.items() if k not in ("runs", "newly_posted_labels")}), flush=True)
    return status


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--prepare-from", type=Path)
    parser.add_argument("--once", action="store_true")
    args = parser.parse_args()
    directory = args.directory.resolve()
    if args.prepare_from:
        spec = prepare(args.prepare_from.resolve(), directory)
        print(json.dumps(dict(prepared_jobs=len(spec["planned_runs"]), directory=str(directory))))
        return
    spec = read(directory / "dispatch_manifest.json")
    plan = read(directory / "queue_plan.json")
    require(all(spec[k] == v for k, v in plan.items() if k != "runs"), "dispatch plan changed")
    for name, checksum in plan["frozen_files"].items():
        require(digest(child(directory, name)) == checksum, "frozen input changed: " + name)
    ref = api(f"repos/{spec['repository']}/git/ref/tags/{spec['workflow_ref']}")
    require(ref["object"]["sha"] == spec["workflow_sha"], "remote execution tag changed")
    # OS file locking releases automatically if the controller is interrupted.
    with (directory / "controller.lock").open("a+b") as lock:
        lock.seek(0)
        if os.name == "nt":
            import msvcrt
            msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
        else:
            import fcntl
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        write(directory / "controller_process.json", dict(pid=os.getpid(), started_utc=now(),
              controller_sha256=digest(Path(__file__))))
        while True:
            try:
                status = cycle(directory, spec)
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
                # Dispatch has its own retained attempt and is never retried
                # here. A failed read/collection can resume at the next cycle.
                write(directory / "network_retry_status.json", dict(updated_utc=now(), error=str(error),
                      recovering=True, controller_pid=os.getpid()))
                print("RECONNECT Actions read/collection: " + str(error), flush=True)
                if args.once:
                    raise
                time.sleep(30)
                continue
            retry_path = directory / "network_retry_status.json"
            if retry_path.exists():
                retry = read(retry_path)
                if retry["recovering"]:
                    retry.update(recovering=False, recovered_utc=now())
                    write(retry_path, retry)
            if args.once or status["ready_for_remote_review"]:
                break
            # A successful run can become visible before its artifact index.
            # Keep collecting until every successful artifact is available.
            if (status["all_jobs_terminal"] and not status["collection_errors"] and
                    any(run["conclusion"] != "success" for run in status["runs"])):
                break
            time.sleep(30)


if __name__ == "__main__":
    main()
