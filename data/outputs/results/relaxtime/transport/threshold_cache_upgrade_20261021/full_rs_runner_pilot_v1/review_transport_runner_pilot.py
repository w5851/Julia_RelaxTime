"""Review downloaded transport pilot artifacts and a conservative runner budget.

Run this reviewer on the compute host. Its inputs retain the GitHub API replies,
downloaded ZIPs, dispatch specification and the exact committed source archive.
Whole-job measurements are never relabelled as steady per-point measurements.
"""
from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import re
import stat
import tomllib
import zipfile


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finite(value, name, *, positive=False):
    number = float(value)
    require(math.isfinite(number) and (number > 0 if positive else number >= 0),
            f"invalid {name}")
    return number


def child(root, relative):
    root = root.resolve()
    path = (root / relative).resolve()
    require(path != root and path.is_relative_to(root), "path leaves its evidence directory")
    return path


def verified_zip(archive, unpacked, checksum):
    require(digest(archive) == checksum, "downloaded archive checksum differs")
    unpacked = unpacked.resolve()
    with zipfile.ZipFile(archive) as bundle:
        members = [m for m in bundle.infolist() if not m.is_dir()]
        paths = [child(unpacked, m.filename) for m in members]
        require(len(paths) == len(set(paths)), "duplicate archive destinations")
        for member, path in zip(members, paths):
            require(not stat.S_ISLNK(member.external_attr >> 16), "archive contains a symlink")
            require(path.is_file() and path.read_bytes() == bundle.read(member),
                    f"extracted archive member differs: {member.filename}")
        require({p.resolve() for p in unpacked.rglob("*") if p.is_file()} == set(paths),
                "extracted archive has missing or additional files")


def parse_process_time(text):
    def field(label):
        match = re.search(r"^\s*" + re.escape(label) + r":\s*(\S.*)$", text, re.MULTILINE)
        require(match is not None, f"GNU time field is missing: {label}")
        return match[1].strip()

    parts = field("Elapsed (wall clock) time (h:mm:ss or m:ss)").split(":")
    require(2 <= len(parts) <= 3, "invalid GNU time elapsed format")
    seconds = 0.0
    for part in parts:
        seconds = 60 * seconds + finite(part, "GNU time elapsed component")
    require(seconds > 0 and field("Exit status") == "0", "process did not complete successfully")
    return dict(wall_seconds=seconds,
                user_seconds=finite(field("User time (seconds)"), "user CPU time"),
                system_seconds=finite(field("System time (seconds)"), "system CPU time"),
                peak_rss_bytes=1024 * finite(field("Maximum resident set size (kbytes)"),
                                             "peak resident memory", positive=True))


def check_point_set(rows, job):
    args = job["inputs"]
    temperatures = ([float(job["expected_T_MeV"])] if args["mode"] == "fixed-muB-phase-scaled"
                    else [float(v) for v in args["t_list"].split(",")])
    expected = list(itertools.product(temperatures,
                    [float(v) for v in args["muB_list"].split(",")],
                    [float(v) for v in args["xi_list"].split(",")]))
    actual = [tuple(float(row[k]) for k in ("T_MeV", "muB_MeV", "xi")) for row in rows]
    require(len(expected) == len(set(expected)) == len(actual) == job["expected_points"],
            "pilot point count differs from dispatch")
    require(len(set(actual)) == len(actual), "duplicate pilot point")
    # A freshly solved phase anchor may differ by roundoff across platforms.
    unmatched = list(expected)
    for point in actual:
        matches = [other for other in unmatched if
                   abs(point[0] - other[0]) <= 1e-5 and
                   abs(point[1] - other[1]) <= 1e-10 and
                   abs(point[2] - other[2]) <= 1e-10]
        require(len(matches) == 1, "pilot coordinates differ from dispatch")
        unmatched.remove(matches[0])


def expected_argv(args):
    argv = ["--mode", args["mode"], "--case-name", args["case_name"],
            "--sampling-profile", args["sampling_profile"], "--xi-list", args["xi_list"],
            "--muB-list", args["muB_list"], "--channel-diagnostics"]
    argv += (["--alphaT-list", args["alpha_t_list"]] if args["mode"] == "fixed-muB-phase-scaled"
             else ["--T-list", args["t_list"]])
    return argv


def check_timing(timing, process, args, source_root, run_id, workflow_sha):
    require(timing["entrypoint_returned"] and timing["julia_version"] == "1.12.5",
            "timed entrypoint failed or used another Julia version")
    require(all(timing[k] == 1 for k in ("julia_threads", "blas_threads_at_start",
                                        "blas_threads_after_load", "blas_threads_at_end")),
            "runner did not use one Julia and one BLAS thread")
    require(timing["calculation_sha"] == args["calculation_sha"] and
            timing["workflow_sha"] == workflow_sha and str(timing["github_run_id"]) == str(run_id),
            "timing identity differs from dispatch")
    require(timing["argv"] == expected_argv(args), "actual CLI argv differs from dispatch")
    for field, relative in (("entrypoint_sha256", "scripts/relaxtime/run_phase_guided_transport_scan.jl"),
                            ("timing_script_sha256", "benchmark/relaxtime/transport_runner_timing.jl")):
        require(timing[field] == digest(source_root / relative), "timing source differs from calculation")
    require(timing["compile_time_included_in_elapsed"], "ambiguous compilation timing")
    for stage in ("load", "scan"):
        phase = timing[stage]
        require(phase["compile_timing_available"], "JIT timing is missing")
        for field in ("elapsed_seconds", "gc_seconds", "compile_seconds", "recompile_seconds", "allocated_bytes"):
            finite(phase[field], f"{stage}.{field}")
    phases = timing["load"]["elapsed_seconds"] + timing["scan"]["elapsed_seconds"]
    measured = finite(timing["measured_wrapper_seconds"], "wrapper elapsed")
    initialization = finite(timing["wrapper_initialization_seconds"], "wrapper initialization")
    require(measured >= phases and process["wall_seconds"] + 0.02 >= measured + initialization,
            "inconsistent nested timing intervals")
    memory = finite(timing["total_memory_bytes"], "runner memory", positive=True)
    require(process["peak_rss_bytes"] <= memory, "peak RSS exceeds reported host memory")


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def review_run(directory, entry, specification, source_root):
    job, run_id = entry["job"], entry["run_id"]
    args = job["inputs"]
    require(args["run_tier"] == "production" and args["render_plots"] is False,
            "cost pilot must use the numerical production tier without plotting")
    run = read(directory / "run.json")
    require(run["id"] == run_id and run["event"] == "workflow_dispatch" and
            run["status"] == "completed" and run["conclusion"] == "success" and
            run["head_sha"] == specification["workflow_sha"] and
            run["head_branch"] == specification["workflow_ref"] and
            run["repository"]["full_name"] == specification["repository"] and
            run["path"].split("@")[0] == specification["workflow_path"],
            "GitHub run does not identify the successful requested workflow")
    require(run["run_attempt"] == entry.get("run_attempt", 1), "unrecorded workflow rerun")
    jobs = read(directory / "jobs.json")["jobs"]
    require(len(jobs) == 1 and jobs[0]["name"] == "generate" and
            jobs[0]["status"] == "completed" and jobs[0]["conclusion"] == "success" and
            jobs[0]["run_id"] == run_id, "pilot job did not complete successfully")
    job_record = jobs[0]
    started, completed = (datetime.fromisoformat(job_record[k].replace("Z", "+00:00"))
                          for k in ("started_at", "completed_at"))
    job_seconds = finite((completed - started).total_seconds(), "Actions job elapsed", positive=True)
    artifact = read(directory / "artifact.json")
    require(artifact["name"] == f"transport-{args['case_name']}-{args['shard_label']}" and
            not artifact["expired"] and artifact["workflow_run"]["id"] == run_id and
            artifact["workflow_run"]["head_sha"] == specification["workflow_sha"],
            "artifact metadata differs from requested run")
    archive = read(directory / "archive_manifest.json")
    checksum = archive["sha256"]
    require((directory / "artifact.zip").stat().st_size == archive["bytes"], "download size changed")
    if artifact.get("digest"):
        require(artifact["digest"] == "sha256:" + checksum, "GitHub artifact digest differs")
    unpacked = directory / "artifact"
    verified_zip(directory / "artifact.zip", unpacked, checksum)
    times = list(unpacked.rglob("runtime_timing.toml"))
    manifests = list(unpacked.rglob("run_manifest.json"))
    require(len(times) == len(manifests) == 1, "missing or ambiguous pilot artifact")
    timing = tomllib.loads(times[0].read_text(encoding="utf-8-sig"))
    process = parse_process_time((times[0].parent / "process_time.txt").read_text(encoding="utf-8"))
    check_timing(timing, process, args, source_root, run_id, specification["workflow_sha"])
    require(job_seconds + 2 >= process["wall_seconds"], "Actions elapsed is shorter than its scan process")
    result = manifests[0].parent
    checker = load_module("pilot_artifact", source_root / "scripts/relaxtime/check_transport_production_artifact.py")
    profile = source_root / "config/models/pnjl/transport_sampling" / (args["sampling_profile"] + ".toml")
    checked = checker.check_artifact(result, profile, args["calculation_sha"], "production", source_root)
    rows = checker.read_rows(result / "phase_guided_transport_scan.csv")
    check_point_set(rows, job)
    cfg = read(result / "effective_config.json")["options"]
    summary = read(result / "run_manifest.json")["summary"]
    require(cfg["compute_bulk"] and cfg["propagator_xi_policy"] == "match_thermo" and
            summary["skipped_count"] == 0 and summary["success_count"] == job["expected_points"],
            "pilot changed physics or reused existing results")
    inputs = load_module("pilot_inputs", source_root / "scripts/analysis/relaxtime/prepare_transport_production_inputs.py")
    contract = read(result / "phase_guided_transport_scan_sampling_contract.json")["identity"]
    present = {p.relative_to(source_root).as_posix() for p in source_root.rglob("*")
               if p.is_file() and inputs.is_numerical_source(p.relative_to(source_root).as_posix())}
    require({item["path"] for item in contract["sources"]} == present, "incomplete numerical source identity")
    diagnostics = checker.read_rows(result / "phase_guided_transport_scan_sampling.csv")
    channel_times = []
    for row in rows:
        seconds = sum(finite(d["seconds"], "channel timing") for d in diagnostics
                      if checker.key(d) == checker.key(row) and d["run_id"] == row["run_id"])
        channel_times.append(dict(T_MeV=float(row["T_MeV"]), muB_MeV=float(row["muB_MeV"]),
                                  xi=float(row["xi"]), channel_seconds=seconds))
    return dict(mode=entry["mode"], label=job["label"], run_id=run_id, run_url=run["html_url"],
                calculation_sha=args["calculation_sha"], artifact_sha256=checksum,
                completed_points=checked["completed_points"], qualified_channels=checked["qualified_channels"],
                process=process, job_seconds=job_seconds, load=timing["load"], scan=timing["scan"],
                cpu_name=timing["cpu_name"], cpu_threads=timing["cpu_threads"],
                total_memory_bytes=timing["total_memory_bytes"], channel_times=channel_times,
                full_single_point_process_seconds=process["wall_seconds"] if len(rows) == 1 else None,
                steady_point_seconds=None,
                input_evidence_sha256={name: digest(directory / name) for name in
                    ("run.json", "jobs.json", "artifact.json", "archive_manifest.json")})


def calibrate_cost(runs, *, target_seconds=10800, safety_factor=2.0, proposed_points=20):
    require(runs and target_seconds > 0 and safety_factor >= 1 and proposed_points >= 1,
            "invalid cost planning inputs")
    scan = max(finite(r["scan"]["elapsed_seconds"], "pilot scan", positive=True) for r in runs)
    overhead = max(max(0.0, finite(r["job_seconds"], "pilot job elapsed") -
                       finite(r["scan"]["elapsed_seconds"], "pilot scan")) for r in runs)
    memory_fraction = max(finite(r["process"]["peak_rss_bytes"], "peak RSS", positive=True) /
                          finite(r["total_memory_bytes"], "host memory", positive=True) for r in runs)
    require(memory_fraction <= 0.8, "pilot peak RSS exceeds 80 percent of host memory")
    reserved_overhead = safety_factor * overhead
    per_point_budget = safety_factor * scan
    points = min(proposed_points, math.floor((target_seconds - reserved_overhead) / per_point_budget))
    require(points >= 1, "a one-point shard exceeds the planning budget")
    return dict(max_points_per_shard=points, job_target_seconds=target_seconds,
                safety_factor=safety_factor, reserved_overhead_seconds=reserved_overhead,
                per_point_budget_seconds=per_point_budget, maximum_pilot_scan_seconds=scan,
                maximum_peak_memory_fraction=memory_fraction, steady_point_seconds=None,
                basis="Twice the largest measured whole pilot scan per budgeted point, without division by point count; twice the observed job overhead is reserved once",
                scope="Conservative plan from measured pilots, not a proven upper bound for every grid point")


def review(directory):
    spec_path = directory / "dispatch_manifest.json"
    specification = read(spec_path)
    source_root = directory / "source"
    source = read(directory / "source_manifest.json")
    require(source["calculation_sha"] == specification["calculation_sha"] == specification["workflow_sha"],
            "source archive and workflow use different calculations")
    verified_zip(directory / "calculation-source.zip", source_root, source["sha256"])
    entries = specification["runs"]
    require(entries and len({e["run_id"] for e in entries}) == len(entries), "missing or duplicate pilot runs")
    labels = {(entry["mode"], entry["job"]["label"]) for entry in entries}
    expected = {(mode, job["label"]) for mode, filename in specification["bound_inputs"].items()
                for job in read(child(directory, filename))["pilot_jobs"]}
    require(labels == expected and len(labels) == len(entries), "pilot collection is incomplete")
    for entry in entries:
        bound = read(child(directory, specification["bound_inputs"][entry["mode"]]))
        require(entry["job"] in bound["pilot_jobs"] and
                entry["job"]["inputs"]["calculation_sha"] == specification["calculation_sha"],
                "dispatch job differs from the bound pilot")
    runs = [review_run(child(directory, e["directory"]), e, specification, source_root) for e in entries]
    budgets = {mode: calibrate_cost([run for run in runs if run["mode"] == mode])
               for mode in sorted({run["mode"] for run in runs})}
    report = dict(pilot_artifacts_passed=True, runner_cost_qualified=True, modes=budgets, runs=runs,
                  whole_grid_transport_calculated=False, formal_production_promoted=False,
                  reviewer_sha256=digest(Path(__file__)), dispatch_manifest_sha256=digest(spec_path),
                  source_manifest_sha256=digest(directory / "source_manifest.json"))
    with (directory / "reviewed_pilot.json").open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    print(json.dumps(review(args.directory.resolve()), indent=2))
