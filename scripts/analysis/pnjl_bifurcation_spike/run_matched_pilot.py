"""Run identical cold/warm workloads in separate Julia processes on one runner."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
METHODS = ("palc", "seed", "multi", "production", "dense", "legacy_seed", "seed_reverse")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=("t120_rho035", "t130_rho080"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--timeout-seconds", type=int, default=480)
    parser.add_argument("--profile", choices=("local", "long", "higher"), default="local")
    parser.add_argument("--methods", nargs="+", choices=METHODS, default=list(METHODS))
    parser.add_argument("--repetitions", type=int, default=6)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    source = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    records = []
    started = time.monotonic()
    environment = dict(os.environ, PALC_PILOT_PROFILE=args.profile, PALC_PILOT_REPETITIONS=str(args.repetitions))
    for method in args.methods:
        destination = args.output / f"{method}.json"
        command = ["julia", "--startup-file=no", "--project=scripts/analysis/pnjl_bifurcation_spike",
                   "scripts/analysis/pnjl_bifurcation_spike/run_matched_pilot.jl",
                   method, args.case, str(destination.resolve())]
        method_started = time.monotonic()
        print(f"Starting {args.case} {method}", flush=True)
        try:
            result = subprocess.run(command, cwd=ROOT, env=environment, timeout=args.timeout_seconds, check=False)
            status, code = "finished", result.returncode
        except subprocess.TimeoutExpired:
            status, code = "timeout", None
        records.append(dict(method=method, status=status, returncode=code,
                            process_seconds=time.monotonic() - method_started))
        (args.output / "processes.json").write_text(json.dumps(dict(
            schema="palc_matched_processes_v2", case=args.case, source_sha=source,
            profile=args.profile, expected_methods=args.methods, expected_repetitions=args.repetitions,
            actions_run_id=os.environ.get("GITHUB_RUN_ID"),
            actions_run_attempt=os.environ.get("GITHUB_RUN_ATTEMPT"),
            timeout_seconds=args.timeout_seconds, processes=records,
            runner_compute_seconds=time.monotonic() - started), indent=2) + "\n", encoding="utf-8")
    if any(r["status"] != "finished" or r["returncode"] != 0 for r in records):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
