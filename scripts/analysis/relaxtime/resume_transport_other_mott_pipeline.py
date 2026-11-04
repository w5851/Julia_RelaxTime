#!/usr/bin/env python3
"""Resume missing reference runs on quote-remote without changing tolerances."""
from pathlib import Path
import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import time
import tomllib
import traceback

RECOVERY = "recovery_v1"
PROFILE = "tighter_budget_v1"


def reference_recovery_plan(jobs, pilot_ids, states, completed_ids):
    """Keep completed/running runs; replace only failed or unstarted work."""
    reused, external, pending = ["tighter_pilot_other01_v1"], [], []
    planned = set()
    for job in jobs:
        wanted = set(job["point_ids"])
        if len(wanted) != len(job["point_ids"]) or planned & wanted:
            raise ValueError("Duplicate reference point identities")
        planned.update(wanted)
        needed = [point for point in job["point_ids"] if point not in pilot_ids]
        name = job["name"]
        status = states.get(name, {}).get("status", "not_started")
        if status == "completed_diagnostic":
            if set(completed_ids[name]) != set(needed):
                raise ValueError("Completed reference coverage mismatch: " + name)
            reused.append(name)
        elif status == "running":
            if set(states[name]["point_ids"]) != set(needed):
                raise ValueError("Running reference coverage mismatch: " + name)
            external.append(dict(name=name, point_ids=needed))
        elif status in ("failed", "not_started"):
            if completed_ids.get(name):
                raise ValueError("Partial completed points require explicit salvage: " + name)
            if not needed:
                continue
            if job["curve_id"] == "other05" and "other05_p02" in needed:
                pending.append(dict(name="tighter_other05_budget_probe_v1", curve_id="other05",
                    variant=PROFILE, point_ids=["other05_p02"], priority=0))
                needed.remove("other05_p02")
            if needed:
                priority = {"other07": 1, "other08": 2, "other09": 3, "other06": 4,
                            "other11": 5, "other05": 10}.get(job["curve_id"], 6)
                pending.append(dict(name="tighter_"+job["curve_id"]+"_recovery_v1",
                    curve_id=job["curve_id"], variant=PROFILE, point_ids=needed, priority=priority))
        else:
            raise ValueError("Unsupported reference status: " + status)
    if not set(pilot_ids) <= planned:
        raise ValueError("Pilot points are absent from the reference plan")
    return dict(reused=reused, external=external, jobs=sorted(pending, key=lambda j: (j["priority"], j["name"])),
                expected_reference_ids=sorted(planned))


def check_budget_only_change(original, updated):
    before, after = json.loads(json.dumps(original)), json.loads(json.dumps(updated))
    before.pop("id")
    after.pop("id")
    for name in ("max_evaluations", "max_nodes"):
        if after["sampling"][name] < before["sampling"][name]:
            raise ValueError("Recovery budget must not decrease")
        after["sampling"][name] = before["sampling"][name]
    if before != after:
        raise ValueError("Recovery may change only profile identity and resource caps")


def load_pipeline(case):
    spec = importlib.util.spec_from_file_location("frozen_mott_pipeline", case/"run_transport_other_mott_pipeline_v2.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def prepare(case, base):
    config_path = case/(RECOVERY+".config.json")
    if config_path.exists():
        raise FileExistsError(config_path)
    failure = base.read_json(case/"tighter_other05_v1.status.json")
    assert failure["status"] == "failed" and failure["source_all_match"]
    assert "evaluation_or_seed_budget" in (case/"tighter_other05_v1.stderr.log").read_text(encoding="utf-8-sig")
    original = case/"inputs/tighter.toml"
    updated = case/"inputs"/(PROFILE+".toml")
    if updated.exists():
        raise FileExistsError(updated)
    content = original.read_text(encoding="utf-8")
    content = content.replace('id = "eta_xi_tighter_v1"', 'id = "other_mott_tighter_budget_v1"')
    content = content.replace("max_evaluations = 80000", "max_evaluations = 320000")
    content = content.replace("max_nodes = 24000", "max_nodes = 96000")
    check_budget_only_change(tomllib.loads(original.read_text(encoding="utf-8")), tomllib.loads(content))
    jobs = base.read_json(case/"convergence_plan_v1/summary.json")["jobs"]
    states, completed = {}, {}
    for job in jobs:
        status = case/(job["name"]+".status.json")
        if status.exists():
            states[job["name"]] = base.read_json(status)
        table = case/job["name"]/"observables.csv"
        if table.exists():
            completed[job["name"]] = [r["point_id"] for r in base.read_csv(table)]
    pilot = {r["point_id"] for r in base.read_csv(case/"tighter_pilot_other01_v1/observables.csv")}
    plan = reference_recovery_plan(jobs, pilot, states, completed)
    old_state = base.read_json(case/"pipeline_v2.status.json")
    assert old_state["current_stage_failed"] == ["tighter_other05_v1"]
    assert old_state["pending_jobs"] == 0
    with updated.open("x", encoding="utf-8") as stream:
        stream.write(content)
    config = dict(schema="other_mott_reference_recovery_v1", created_utc=base.now(), **plan,
        max_workers=3, minimum_launch_free_GiB=3.25,
        fine_runs=list(base.read_json(case/"fine_slopes_v1/summary.json")["provenance"]),
        original_failure=failure,
        files_sha256={p.relative_to(case).as_posix(): base.sha(p) for p in (
            case/"pipeline_v2.config.json", case/"convergence_plan_v1/transport_plan.csv",
            case/"convergence_plan_v1/summary.json", updated, case/Path(__file__).name,
            case/"test_transport_other_mott_recovery_v2.py")},
        change="Only max_evaluations 80000->320000 and max_nodes 24000->96000; tolerances and all quadrature settings unchanged.",
        package="packages/complete_recovery_snapshot_v1/archive_receipt.json")
    base.write_json(config_path, config)
    print(json.dumps(config), flush=True)


def run(case, base):
    config = base.read_json(case/(RECOVERY+".config.json"))
    original = base.read_json(case/"pipeline_v2.config.json")
    source = base.read_json(case/"expected_source_hashes.json")
    lock = (case/(RECOVERY+".lock")).open("x", encoding="utf-8")
    lock.write(str(os.getpid()))
    lock.flush()
    state = dict(status="running_reference_recovery", pid=os.getpid(), started_utc=base.now())

    def checkpoint(**changes):
        state.update(changes, updated_utc=base.now())
        base.write_json(case/(RECOVERY+".status.json"), state)

    def verify():
        assert len(source) == 340
        for name, digest in source.items():
            assert base.sha(base.SOURCE/name) == digest, name
        for names, folder in ((original["script_sha256"], case), (original["input_sha256"], case/"inputs"),
                              (config["files_sha256"], case)):
            for name, digest in names.items():
                assert base.sha(folder/name) == digest, name

    def completed(name, point_ids=None):
        status = base.read_json(case/(name+".status.json"))
        assert status["status"] == "completed_diagnostic" and status["source_all_match"]
        assert status["source_files_checked_after"] == 340
        summary = base.read_json(case/name/"summary.json")
        assert summary["status"] == "completed_diagnostic"
        for file, digest in summary["files"].items():
            assert base.sha(case/name/file) == digest
        actual = [r["point_id"] for r in base.read_csv(case/name/"observables.csv")]
        assert len(actual) == len(set(actual))
        if point_ids is not None:
            assert set(actual) == set(point_ids), name
        return actual

    def stage(script, arguments, log):
        with (case/log).open("xb") as stream:
            subprocess.run([sys.executable, str(case/script), *map(str, arguments)], cwd=case,
                stdout=stream, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                creationflags=subprocess.CREATE_NO_WINDOW, check=True)

    active, external = [], list(config["external"])
    pending, references, failed = list(config["jobs"]), list(config["reused"]), []
    try:
        verify()
        for name in references:
            completed(name)
        while pending or active or external:
            for item in list(external):
                status = base.read_json(case/(item["name"]+".status.json"))
                if status["status"] == "running":
                    continue
                if status["status"] == "completed_diagnostic":
                    completed(item["name"], item["point_ids"])
                    references.append(item["name"])
                else:
                    failed.append(item["name"])
                external.remove(item)
            for worker in list(active):
                if worker["process"].poll() is None:
                    continue
                for stream in worker["streams"]:
                    stream.close()
                verify()
                job = worker["job"]
                status = dict(status="completed_diagnostic" if worker["process"].returncode == 0 else "failed",
                    exit_code=worker["process"].returncode, started_utc=worker["started"], finished_utc=base.now(),
                    source_files_checked_before=340, source_files_checked_after=340, source_all_match=True,
                    arguments=worker["argv"], point_ids=job["point_ids"], recovery_config_sha256=base.sha(case/(RECOVERY+".config.json")))
                base.write_json(case/(job["name"]+".status.json"), status)
                if status["status"] == "completed_diagnostic":
                    completed(job["name"], job["point_ids"])
                    references.append(job["name"])
                else:
                    failed.append(job["name"])
                active.remove(worker)
                print("FINISHED", job["name"], status["status"], flush=True)
            if failed:
                pending.clear()
            if pending and len(active)+len(external) < config["max_workers"] and base.free_gib() >= config["minimum_launch_free_GiB"]:
                verify()
                job = pending.pop(0)
                assert not (case/job["name"]).exists() and not (case/(job["name"]+".status.json")).exists()
                env = dict(os.environ, JULIA_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
                    JULIA_PKG_PRECOMPILE_AUTO="0", MOTT_LOCATE_DIR=str(case/"convergence_plan_v1"),
                    MOTT_POINT_IDS=",".join(job["point_ids"]))
                argv = [str(base.JULIA), "--startup-file=no", "--threads=1", "--project="+str(base.SOURCE),
                    str(case/"run_transport_other_mott_refinement.jl"), str(base.SOURCE), str(case/"inputs"),
                    str(case/job["name"]), "transport", job["curve_id"], job["variant"]]
                streams = [(case/(job["name"]+suffix)).open("xb") for suffix in (".stdout.log", ".stderr.log")]
                process = subprocess.Popen(argv, cwd=base.SOURCE, env=env, stdin=subprocess.DEVNULL,
                    stdout=streams[0], stderr=streams[1], creationflags=subprocess.CREATE_NO_WINDOW)
                started = base.now()
                active.append(dict(job=job, process=process, streams=streams, started=started, argv=argv))
                base.write_json(case/(job["name"]+".status.json"), dict(status="running", pid=process.pid,
                    point_ids=job["point_ids"], started_utc=started, arguments=argv))
                print("STARTED", job["name"], "points", len(job["point_ids"]), "pid", process.pid, flush=True)
            checkpoint(active_jobs=[dict(name=w["job"]["name"], pid=w["process"].pid) for w in active],
                external_jobs=external, pending_jobs=len(pending), finished_runs=references, failed_runs=failed)
            if pending or active or external:
                time.sleep(10)
        if failed:
            raise RuntimeError("Reference recovery failed; retained outputs: "+",".join(failed))
        ids = [point for name in references for point in completed(name)]
        assert len(ids) == len(set(ids)) and set(ids) == set(config["expected_reference_ids"])
        fine_ids = [point for name in config["fine_runs"] for point in completed(name)]
        assert len(fine_ids) == len(set(fine_ids)) == 195
        stage(base.REVIEWER, [case, case/"final_slopes_recovery_v1", "--runs", ",".join(config["fine_runs"]),
            "--comparison-runs", ",".join(references)], "final_slopes_recovery_v1.log")
        verify()
        all_runs = list(dict.fromkeys(["survey_v1", "endpoint_checks_v1", "locate_v1", "midpoint_screen_v1",
            *config["fine_runs"], *references]))
        # Include the recovery config and original failed run as evidence, outside accepted runs.
        evidence = case/"recovery_evidence_v1"
        evidence.mkdir(exist_ok=False)
        for path in [case/(RECOVERY+".config.json"), case/"recovery_unit_v2.stdout.log",
                     case/"recovery_unit_v2.stderr.log", *case.glob("tighter_other05_v1.*")]:
            if path.is_file():
                shutil.copyfile(path, evidence/path.name)
        shutil.copytree(case/"tighter_other05_v1", evidence/"failed_tighter_other05_v1")
        stage(base.PACKAGER, [case, base.SOURCE, "complete_recovery_snapshot_v1", "--runs", ",".join(all_runs),
            "--reviews", "survey_review_v2,midpoint_review_v1,fine_slopes_v1,convergence_plan_v1,final_slopes_recovery_v1,pilot_plan_other01_v1,recovery_evidence_v1"],
            "complete_recovery_snapshot_v1.log")
        checkpoint(status="ready_for_agent_review", finished_utc=base.now(), fine_runs=config["fine_runs"],
            reference_runs=references, final_review="final_slopes_recovery_v1/summary.json", package=config["package"])
    except BaseException as error:
        checkpoint(status="failed", error=str(error), traceback=traceback.format_exc())
        raise
    finally:
        lock.close()
        (case/(RECOVERY+".lock")).unlink()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("--prepare", action="store_true")
    args = parser.parse_args()
    case = args.case.resolve()
    pipeline = load_pipeline(case)
    prepare(case, pipeline) if args.prepare else run(case, pipeline)
