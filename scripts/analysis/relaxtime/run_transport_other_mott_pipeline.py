#!/usr/bin/env python3
"""Unattended quote-remote diagnostic queue, with precision and xi controls."""
from pathlib import Path
from datetime import datetime, timezone
from collections import defaultdict
import argparse
import csv
import ctypes
import hashlib
import json
import os
import subprocess
import sys
import time
import traceback

SOURCE = Path("D:/CodexCompute/trfull-20261026-v1/case/source")
JULIA = Path("D:/CodexCompute/Julia_RelaxTime/transport-preflight-20261023-v1/runtime/julia-1.12.5/bin/julia.exe")
PREREQUISITES = ("locate_v1", "midpoint_screen_v1", "fine_pilot_other01_v1", "tighter_pilot_other01_v1")
REVIEWER = "review_transport_other_mott_slopes_v3.py"
PACKAGER = "package_transport_other_mott_v2.py"
PIPELINE_ID = "pipeline_v2"
SCRIPTS = ("run_transport_other_mott_refinement.jl", "run_transport_mott_probe.jl",
           "transport_mott_thresholds.jl", "transport_mott_kaon_endpoint.jl",
           REVIEWER, "review_transport_other_mott_midpoints.py",
           "prepare_transport_other_mott_convergence.py", PACKAGER,
           Path(__file__).name, "run_pipeline_v2.ps1", "test_transport_other_mott_slopes_v2.py")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_json(path, value):
    temporary = path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    temporary.replace(path)


def now():
    return datetime.now(timezone.utc).isoformat()


def free_gib():
    class MemoryStatus(ctypes.Structure):
        _fields_ = [("length", ctypes.c_ulong), ("load", ctypes.c_ulong),
                    *[(n, ctypes.c_ulonglong) for n in ("total", "available", "total_page", "available_page",
                                                       "total_virtual", "available_virtual", "extended")]]
    status = MemoryStatus()
    status.length = ctypes.sizeof(status)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
        raise OSError("Cannot obtain physical-memory headroom")
    return status.available/1024**3


def prepare(case):
    path = case / (PIPELINE_ID+".config.json")
    if path.exists():
        raise FileExistsError(path)
    config = dict(schema="transport_other_mott_pipeline_v2", host="quote-remote", max_workers=3,
        minimum_launch_free_GiB=3.25, prerequisites=list(PREREQUISITES), prepared_utc=now(),
        calculation_sha="2e5c411c716d5cbc524c3fc6cfb335fe05b39444",
        script_sha256={name: sha(case/name) for name in SCRIPTS},
        input_sha256={name: sha(case/"inputs"/name) for name in
                      ("source_rows.csv", "targets.json", "base.toml", "fine.toml", "tighter.toml")},
        numerical_scope="All 13 selected closed-threshold windows; three final coefficients, six tau, full 17 channels, species weights; automatic s and xi controls.",
        final_state="ready_for_agent_review; no promotion or PDF generation")
    write_json(path, config)
    print(json.dumps(config), flush=True)


def run(case):
    config = read_json(case/(PIPELINE_ID+".config.json"))
    expected = read_json(case/"expected_source_hashes.json")
    lock_path = case/(PIPELINE_ID+".lock")
    lock = lock_path.open("x", encoding="utf-8")
    lock.write(str(os.getpid()))
    lock.flush()
    state = dict(status="waiting_for_prerequisites", started_utc=now(), pid=os.getpid(), jobs=[])

    def checkpoint(**values):
        state.update(values, updated_utc=now())
        write_json(case/(PIPELINE_ID+".status.json"), state)

    def verify_source():
        assert len(expected) == 340
        for name, digest in expected.items():
            assert sha(SOURCE/name) == digest, "Frozen numerical source changed: "+name
        for name, digest in config["script_sha256"].items():
            assert sha(case/name) == digest, "Frozen diagnostic script changed: "+name
        for name, digest in config["input_sha256"].items():
            assert sha(case/"inputs"/name) == digest, "Frozen input changed: "+name

    def python_stage(script, arguments, log):
        with (case/log).open("xb") as stream:
            completed = subprocess.run([sys.executable, str(case/script), *map(str, arguments)],
                stdout=stream, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                creationflags=subprocess.CREATE_NO_WINDOW, check=False)
        if completed.returncode:
            raise RuntimeError(f"{script} failed; retained log {log}")

    def run_jobs(jobs, plan_dir):
        pending, active, finished, failed = list(jobs), [], [], []
        while pending or active:
            for worker in list(active):
                process = worker["process"]
                if process.poll() is None:
                    continue
                for stream in worker["streams"]:
                    stream.close()
                job = worker["job"]
                verify_source()
                status = dict(status="completed_diagnostic" if process.returncode == 0 else "failed",
                    started_utc=worker["started"], finished_utc=now(), exit_code=process.returncode,
                    source_files_checked_before=340, source_files_checked_after=340, source_all_match=True,
                    arguments=worker["arguments"], point_ids=job["point_ids"])
                if process.returncode == 0:
                    summary = read_json(case/job["name"]/"summary.json")
                    assert summary["status"] == "completed_diagnostic" and summary["points"] == len(job["point_ids"])
                    for name, digest in summary["files"].items():
                        assert sha(case/job["name"]/name) == digest
                    finished.append(job["name"])
                else:
                    failed.append(job["name"])
                write_json(case/(job["name"]+".status.json"), status)
                active.remove(worker)
                print("FINISHED", job["name"], "exit", process.returncode, flush=True)
            if failed:
                pending.clear()
            if pending and len(active) < config["max_workers"] and free_gib() >= config["minimum_launch_free_GiB"]:
                verify_source()
                job = pending.pop(0)
                out = case/job["name"]
                if out.exists() or (case/(job["name"]+".status.json")).exists():
                    raise FileExistsError(out)
                env = dict(os.environ, JULIA_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
                    JULIA_PKG_PRECOMPILE_AUTO="0", MOTT_LOCATE_DIR=str(plan_dir), MOTT_POINT_IDS=",".join(job["point_ids"]))
                argv = [str(JULIA), "--startup-file=no", "--threads=1", "--project="+str(SOURCE),
                    str(case/"run_transport_other_mott_refinement.jl"), str(SOURCE), str(case/"inputs"),
                    str(out), "transport", job["curve_id"], job["variant"]]
                streams = ((case/(job["name"]+".stdout.log")).open("xb"),
                           (case/(job["name"]+".stderr.log")).open("xb"))
                process = subprocess.Popen(argv, cwd=SOURCE, env=env, stdout=streams[0], stderr=streams[1],
                    stdin=subprocess.DEVNULL, creationflags=subprocess.CREATE_NO_WINDOW)
                started = now()
                active.append(dict(job=job, process=process, streams=streams, arguments=argv, started=started))
                write_json(case/(job["name"]+".status.json"), dict(status="running", pid=process.pid,
                    started_utc=started, arguments=argv, point_ids=job["point_ids"]))
                print("STARTED", job["name"], "points", len(job["point_ids"]), "pid", process.pid, flush=True)
            checkpoint(active_jobs=[dict(name=w["job"]["name"], pid=w["process"].pid) for w in active],
                       pending_jobs=len(pending), current_stage_finished=finished, current_stage_failed=failed)
            if pending or active:
                time.sleep(10)
        if failed:
            raise RuntimeError("Stopped dispatch after failed workers: "+",".join(failed))
        return finished

    try:
        verify_source()
        deadline = time.time()+48*3600
        while True:
            statuses = {name: read_json(case/(name+".status.json"))
                        if (case/(name+".status.json")).exists() else dict(status="not_started") for name in PREREQUISITES}
            checkpoint(prerequisites={name: value["status"] for name, value in statuses.items()})
            if any(value["status"] == "failed" for value in statuses.values()):
                raise RuntimeError("Prerequisite failed; do not dispatch transport")
            if all(value["status"] == "completed_diagnostic" for value in statuses.values()):
                break
            if time.time() > deadline:
                raise TimeoutError("Prerequisites did not finish within 48 hours")
            time.sleep(15)
        python_stage("review_transport_other_mott_midpoints.py", [case, case/"midpoint_review_v1"], "midpoint_review_v1.log")
        assert not read_json(case/"midpoint_review_v1/summary.json")["new_closed_candidates"], "New closed thresholds need target selection"
        locations = read_json(case/"locate_v1/locations.json")
        targets = read_json(case/"inputs/targets.json")["targets"]
        assert {v["curve_id"] for v in locations} == {v["curve_id"] for v in targets}
        location_checks = []
        for location in locations:
            center = location["levels"]["base"]["xi"]
            shifts = {level: value["xi"]-center for level, value in location["levels"].items()}
            assert max(abs(v) for v in shifts.values()) < 1e-6, location["curve_id"]
            location_checks.append(dict(curve_id=location["curve_id"], meson=location["meson"],
                                        xi=center, shifts=shifts))
        write_json(case/"location_review_v1.json", dict(status="completed_location_controls", locations=location_checks,
            maximum_allowed_shift=1e-6, locator_sha256=sha(case/"locate_v1/locations.json")))
        plan = read_csv(case/"locate_v1/transport_plan.csv")
        assert len(plan) == len({r["point_id"] for r in plan})
        original = {r["point_id"]: r for r in plan}
        pilot = read_csv(case/"fine_pilot_other01_v1/observables.csv")
        for row in pilot:
            assert abs(float(row["xi"])-float(original[row["point_id"]]["xi"])) < 1e-13
        done = {r["point_id"] for r in pilot}
        grouped = defaultdict(list)
        for row in plan:
            if row["point_id"] not in done:
                grouped[row["curve_id"]].append(row["point_id"])
        priority = ("other09", "other07", "other08", "other06", "other01", "other11", "other02",
                    "other04", "other03", "other12", "other13", "other10", "other05")
        jobs = [dict(name="fine_"+curve+"_v1", curve_id=curve, variant="fine", point_ids=grouped[curve])
                for curve in priority if grouped[curve]]
        write_json(case/"fine_queue_plan_v1.json", dict(jobs=jobs, reused_pilot_points=sorted(done), planned_points=len(plan)))
        checkpoint(status="running_fine_transport", jobs=jobs, expected_fine_points=len(plan))
        fine_runs = ["fine_pilot_other01_v1", *run_jobs(jobs, case/"locate_v1")]
        python_stage(REVIEWER, [case, case/"fine_slopes_v1", "--runs", ",".join(fine_runs)], "fine_slopes_v1.log")
        python_stage("prepare_transport_other_mott_convergence.py", [case, case/"convergence_plan_v1"], "convergence_plan_v1.log")
        reference = read_json(case/"convergence_plan_v1/summary.json")
        checkpoint(status="running_precision_and_xi_controls", reference_jobs=reference["jobs"],
                   expected_reference_points=reference["total_evaluations"])
        reference_pilot_ids = {r["point_id"] for r in read_csv(case/"tighter_pilot_other01_v1/observables.csv")}
        for job in reference["jobs"]:
            job["point_ids"] = [p for p in job["point_ids"] if p not in reference_pilot_ids]
        reference_jobs = [j for j in reference["jobs"] if j["point_ids"]]
        ref_runs = ["tighter_pilot_other01_v1", *run_jobs(reference_jobs, case/"convergence_plan_v1")]
        python_stage(REVIEWER, [case, case/"final_slopes_v1", "--runs", ",".join(fine_runs),
            "--comparison-runs", ",".join(ref_runs)], "final_slopes_v1.log")
        verify_source()
        completed = list(dict.fromkeys(["survey_v1", "endpoint_checks_v1", *PREREQUISITES, *fine_runs[1:], *ref_runs]))
        python_stage(PACKAGER, [case, SOURCE, "complete_snapshot_v1", "--runs", ",".join(completed),
            "--reviews", "survey_review_v2,midpoint_review_v1,fine_slopes_v1,convergence_plan_v1,final_slopes_v1"], "complete_snapshot_v1.log")
        checkpoint(status="ready_for_agent_review", finished_utc=now(), fine_runs=fine_runs, reference_runs=ref_runs,
                   final_review="final_slopes_v1/summary.json", package="packages/complete_snapshot_v1/archive_receipt.json")
        print("READY_FOR_AGENT_REVIEW", flush=True)
    except BaseException as error:
        checkpoint(status="failed", error=str(error), traceback=traceback.format_exc())
        raise
    finally:
        lock.close()
        lock_path.unlink()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("--prepare", action="store_true")
    args = parser.parse_args()
    case_root = args.case.resolve()
    prepare(case_root) if args.prepare else run(case_root)
