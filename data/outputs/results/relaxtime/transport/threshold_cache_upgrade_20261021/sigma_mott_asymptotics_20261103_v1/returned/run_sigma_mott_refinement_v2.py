#!/usr/bin/env python3
"""Isolated quote-remote sigma diagnostics; never updates the frozen source."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import csv
import ctypes
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
import traceback
from types import SimpleNamespace

OLD = Path("D:/CodexCompute/trmott-other-slopes-20261029-v1")
SOURCE = Path("D:/CodexCompute/trfull-20261026-v1/case/source")
JULIA = Path("D:/CodexCompute/Julia_RelaxTime/transport-preflight-20261023-v1/runtime/julia-1.12.5/bin/julia.exe")
CALCULATION_SHA = "2e5c411c716d5cbc524c3fc6cfb335fe05b39444"
DEPENDENCIES = ("run_transport_other_mott_refinement.jl", "run_transport_mott_probe.jl",
                "transport_mott_thresholds.jl", "transport_mott_kaon_endpoint.jl")


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def verify_source(case):
    expected = read_json(case / "expected_source_hashes.json")
    if len(expected) != 340:
        raise ValueError("Unexpected frozen numerical source inventory")
    for name, expected_sha in expected.items():
        if sha(SOURCE / name) != expected_sha:
            raise ValueError("Frozen numerical source changed: " + name)
    return len(expected)


def prepare(case):
    if (case / "contract.json").exists() or (case / "inputs").exists():
        raise FileExistsError("Refusing to replace an existing sigma case")
    inputs = case / "inputs"
    inputs.mkdir()
    for name in DEPENDENCIES:
        shutil.copy2(OLD / name, case / name)
    shutil.copy2(OLD / "expected_source_hashes.json", case / "expected_source_hashes.json")
    for name in ("targets.json", "source_rows.csv", "base.toml", "tighter_budget_v1.toml"):
        shutil.copy2(OLD / "inputs" / name, inputs / name)
    prior = inputs / "previous"
    prior.mkdir()
    for name in ("observables.csv", "rates.csv", "species_weights.csv", "channel_contributions.csv", "summary.json"):
        shutil.copy2(OLD / "tighter_other07_recovery_v1" / name, prior / name)
    shutil.copy2(OLD / "locate_v1" / "locations.json", inputs / "locations.json")
    target = next(t for t in read_json(inputs / "targets.json")["targets"] if t["curve_id"] == "other07")
    location = next(t for t in read_json(inputs / "locations.json") if t["curve_id"] == "other07")
    center = location["levels"]["base"]["xi"]
    offsets = [(-2e-4, "m200"), (-1e-4, "m100"), (-5e-5, "m050"), (-2e-5, "m020"),
               (0.0, "c000"), (2e-5, "p020"), (5e-5, "p050"), (1e-4, "p100"), (2e-4, "p200")]
    plan = [dict(point_id="sig_" + tag, curve_id="other07", meson="sigma", cut="uu",
                 window=target["window"], T_MeV=target["T_MeV"], muB_MeV=target["muB_MeV"],
                 xi=center + delta, xi_mott=center, delta_xi=delta, kind="explicit_center_refinement",
                 original_point_id="") for delta, tag in offsets]
    with (inputs / "transport_plan.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(plan[0]))
        writer.writeheader()
        writer.writerows(plan)
    original = (inputs / "tighter_budget_v1.toml").read_text(encoding="utf-8-sig")
    for tag, tolerance in (("s6", "1.0e-06"), ("s7", "1.0e-07"), ("s8", "1.0e-08")):
        profile = original.replace('id = "other_mott_tighter_budget_v1"', f'id = "sigma_center_{tag}"')
        profile = profile.replace("rtol = 1.0e-06", "rtol = " + tolerance)
        if tag != "s6":
            profile = profile.replace("max_evaluations = 320000", "max_evaluations = 960000")
            profile = profile.replace("max_nodes = 96000", "max_nodes = 280000")
            profile = profile.replace("max_rounds = 50", "max_rounds = 80")
        (inputs / (tag + ".toml")).write_text(profile, encoding="utf-8")
    checked = verify_source(case)
    contract = dict(schema="sigma_mott_center_and_s_checks_v1", prepared_utc=now(), host="quote-remote",
        calculation_sha=CALCULATION_SHA, numerical_source=str(SOURCE), julia=str(JULIA),
        target=target, location=location, plan=plan, source_files_verified=checked,
        input_sha256={p.relative_to(case).as_posix(): sha(p) for p in inputs.rglob("*") if p.is_file()},
        dependency_sha256={name: sha(case / name) for name in DEPENDENCIES},
        old_case=str(OLD), maximum_workers=3, minimum_launch_free_GiB=3.25,
        scope="Diagnostic center xi refinement and near-threshold s spectra/integrals; fixed numerical source.",
        production_kernel_changed=False, publication_or_pdf_promotion=False)
    write_json(case / "contract.json", contract)
    print(json.dumps(dict(status="prepared", case=str(case), center=center, points=len(plan), source_files=checked)), flush=True)


def free_gib():
    class MemoryStatus(ctypes.Structure):
        _fields_ = [("length", ctypes.c_ulong), ("load", ctypes.c_ulong),
                    *[(n, ctypes.c_ulonglong) for n in ("total", "available", "total_page", "available_page",
                                                       "total_virtual", "available_virtual", "extended")]]
    status = MemoryStatus()
    status.length = ctypes.sizeof(status)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):
        raise OSError("Cannot obtain physical-memory headroom")
    return status.available / 1024**3


def job(case, args):
    if not args.name.replace("_", "").isalnum():
        raise ValueError("Invalid run name")
    status_path = case / (args.name + ".status.json")
    with status_path.open("x", encoding="utf-8") as stream:
        json.dump(dict(status="preflight", started_utc=now(), pid=os.getpid()), stream)
    state = dict(status="running", started_utc=now(), pid=os.getpid(), name=args.name,
                 stage=args.stage, profile=args.profile, points=args.points, calculation_sha=CALCULATION_SHA)
    try:
        contract = read_json(case / "contract.json")
        protected = dict(contract["input_sha256"], **contract["dependency_sha256"])
        script = "run_transport_other_mott_refinement.jl" if args.stage == "transport" else args.script
        if Path(script).name != script or not script.endswith(".jl"):
            raise ValueError("Julia diagnostic script must be a case-local basename")
        protected[script] = sha(case / script)
        protected[Path(__file__).name] = sha(Path(__file__))
        def verify():
            n = verify_source(case)
            for name, digest in protected.items():
                if sha(case / name) != digest:
                    raise ValueError("Frozen case file changed: " + name)
            return n
        state["source_files_checked_before"] = verify()
        state["file_sha256"] = protected
        headroom = free_gib()
        if headroom < contract["minimum_launch_free_GiB"]:
            raise RuntimeError(f"Insufficient memory headroom: {headroom:.3f} GiB")
        env = os.environ.copy()
        env.update(JULIA_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
                   JULIA_PKG_PRECOMPILE_AUTO="0", MOTT_POINT_IDS=args.points,
                   MOTT_LOCATE_DIR=str(case / "inputs"))
        command = [str(JULIA), "--startup-file=no", "--threads=1", "--project=" + str(SOURCE),
                   str(case / script), str(SOURCE), str(case / "inputs"), str(case / args.name),
                   args.stage, "other07", args.profile]
        state["command"] = command
        state["launch_free_GiB"] = headroom
        with (case / (args.name + ".stdout.log")).open("xb") as out, (case / (args.name + ".stderr.log")).open("xb") as err:
            proc = subprocess.Popen(command, cwd=SOURCE, env=env, stdin=subprocess.DEVNULL,
                                    stdout=out, stderr=err, creationflags=subprocess.CREATE_NO_WINDOW)
            state["julia_pid"] = proc.pid
            write_json(status_path, state)
            code = proc.wait()
        state["exit_code"] = code
        state["source_files_checked_after"] = verify()
        if code != 0:
            raise RuntimeError("Julia diagnostic exited with code " + str(code))
        if not (case / args.name / "summary.json").is_file():
            raise RuntimeError("Julia returned without a completion summary")
        state.update(status="completed_diagnostic", source_all_match=True, finished_utc=now())
        write_json(status_path, state)
        print(json.dumps(state), flush=True)
    except BaseException as exc:
        state.update(status="failed", error=str(exc), finished_utc=now())
        write_json(status_path, state)
        traceback.print_exc()
        raise


def queue_controls(case):
    """Wait for the two retained workers, then run memory-heavy controls serially."""
    status_path = case / "controls_v1.status.json"
    with status_path.open("x", encoding="utf-8") as stream:
        json.dump(dict(status="waiting_for_center_grid", started_utc=now(), pid=os.getpid()), stream)
    tasks = [dict(name="thermal_audit_v1", stage="spectrum", profile="s6", points="all",
                  script="run_sigma_mott_thermal_audit.jl"),
             *[dict(name=name, stage="transport", profile="s7", points=point,
                    script="run_sigma_mott_spectrum.jl") for name, point in
               (("center_s7_v2", "sig_c000"), ("left_s7_v1", "sig_m020"), ("right_s7_v1", "sig_p020"))]]
    state = dict(status="waiting_for_center_grid", started_utc=now(), pid=os.getpid(),
                 tasks=tasks, completed=[], numerical_workers_after_prerequisites=1)
    try:
        while True:
            statuses = [read_json(case/(name+".status.json"))["status"] for name in ("left_s6_v1", "right_s6_v1")]
            if any(value == "failed" for value in statuses):
                raise RuntimeError("A required center-grid worker failed")
            if all(value == "completed_diagnostic" for value in statuses):
                break
            state.update(updated_utc=now(), prerequisite_statuses=statuses, free_GiB=free_gib())
            write_json(status_path, state)
            time.sleep(20)
        for task in tasks:
            while free_gib() < 5.0:
                state.update(status="waiting_for_memory", next_job=task["name"], updated_utc=now(), free_GiB=free_gib())
                write_json(status_path, state)
                time.sleep(20)
            state.update(status="running_control", current_job=task["name"], updated_utc=now())
            write_json(status_path, state)
            job(case, SimpleNamespace(**task))
            state["completed"].append(task["name"])
        state.update(status="completed_diagnostic", finished_utc=now(), current_job=None)
        write_json(status_path, state)
    except BaseException as exc:
        state.update(status="failed", error=str(exc), finished_utc=now())
        write_json(status_path, state)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "job", "queue"))
    parser.add_argument("case", type=Path)
    parser.add_argument("--name")
    parser.add_argument("--stage", choices=("transport", "spectrum"), default="transport")
    parser.add_argument("--profile", default="s6")
    parser.add_argument("--points", default="all")
    parser.add_argument("--script", default="run_sigma_mott_spectrum.jl")
    args = parser.parse_args()
    if args.action == "prepare":
        prepare(args.case.resolve())
    elif args.action == "queue":
        queue_controls(args.case.resolve())
    else:
        if args.name is None:
            parser.error("job requires --name")
        job(args.case.resolve(), args)


if __name__ == "__main__":
    main()
