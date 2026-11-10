"""Bounded production-entry comparison on one runner, with a shared time budget."""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
STRATEGIES = ("legacy", "single_pass", "palc")
ROUTES = ("asym", "rho", "tmu")
ALL_CASES = tuple((route, strategy) for strategy in ("palc", "single_pass", "legacy") for route in ROUTES)


def parse_cases(value):
    if value == "all":
        return ALL_CASES
    cases = tuple(tuple(item.strip().split(":")) for item in value.split(","))
    if not cases or any(case not in ALL_CASES for case in cases) or len(set(cases)) != len(cases):
        raise argparse.ArgumentTypeError("cases must be unique route:strategy pairs or all")
    return cases


def admissible(row):
    finite = lambda values: all(x is not None and math.isfinite(x) for x in values)
    return (row.get("converged", False) and row.get("residual") is not None
            and math.isfinite(row["residual"]) and row["residual"] <= 1e-6
            and finite(row["solution"] + row["mu_vec"] + row["masses"])
            and finite([row[k] for k in ("pressure", "rho_norm", "entropy", "energy")])
            and row["density_error"] is not None and row["density_error"] <= 1e-6
            and row["constraint_error"] is not None and row["constraint_error"] <= 1e-6
            and min(row["masses"]) > 0 and row["masses"][2] > row["masses"][0]
            and all(0 <= x <= 1 for x in row["solution"][3:5]))


def compare_rows(left, right):
    by_key = lambda rows: {(r["T"], r["xi"], r["coordinate"]): r for r in rows}
    a, b = by_key(left), by_key(right)
    mismatches, coverage = [], []
    maxima = dict(solution=0.0, mu_vec=0.0, pressure=0.0, rho_norm=0.0, entropy=0.0, energy=0.0)
    for key in sorted(set(a) | set(b)):
        if key not in a or key not in b or not (admissible(a[key]) and admissible(b[key])):
            coverage.append(dict(key=key, left=admissible(a[key]) if key in a else False,
                                 right=admissible(b[key]) if key in b else False))
            continue
        deltas = {k: (max(abs(x-y) for x,y in zip(a[key][k],b[key][k])) if isinstance(a[key][k],list)
                      else abs(a[key][k]-b[key][k])) for k in maxima}
        for name,value in deltas.items():
            maxima[name] = max(maxima[name],value)
        if any(value > (1e-4 if name in ("solution","mu_vec") else 1e-6) for name,value in deltas.items()):
            mismatches.append(dict(key=key,deltas=deltas,left_pressure=a[key]["pressure"],right_pressure=b[key]["pressure"]))
    duplicates = len(a) != len(left) or len(b) != len(right)
    return dict(total=len(set(a)|set(b)), coverage=coverage, mismatches=mismatches,
                duplicate_keys=duplicates, maxima=maxima,
                same_admissible_results=bool(a) and not duplicates and not coverage and not mismatches)


def complete_round(data, iteration):
    config = data["config"]
    expected = {(t, xi, x) for t in config["temperatures"] for xi in config["anisotropies"] for x in config["grid"]}
    actual = {(r["T"], r["xi"], r["coordinate"]) for r in iteration["rows"]}
    return expected == actual and len(iteration["rows"]) == len(expected)


def promotion_checks(summary, repetitions, cases=ALL_CASES):
    problems = []
    source_shas = set()
    for route in ROUTES:
        group = summary.get(route, {})
        methods = group.get("methods", {})
        required = [strategy for candidate_route, strategy in cases if candidate_route == route]
        for strategy in required:
            method = methods.get(strategy)
            if not method or method["repetitions"] != repetitions or not method["complete_grids"]:
                problems.append(f"{route}/{strategy}: incomplete workload")
                continue
            source_shas.add(method["source_sha"])
            if method["accepted"] != method["total"]:
                problems.append(f"{route}/{strategy}: rejected points")
            if not all(check["same_admissible_results"] for check in method["repeat_checks"]):
                problems.append(f"{route}/{strategy}: repeated results differ or contain rejected points")
        if {"single_pass", "palc"}.issubset(required) and not group.get("comparisons", {}).get("single_pass_vs_palc", {}).get("same_admissible_results", False):
            problems.append(f"{route}: single_pass/PALC parity failed")
        if required and all(s in methods for s in required) and any(methods[s]["config"] != methods[required[0]]["config"] for s in required):
            problems.append(f"{route}: unequal configurations")
    if len(source_shas) != 1:
        problems.append("source SHA is missing or inconsistent")
    return problems


def summarize(output):
    result, lines = {}, ["# Production PALC comparison", "", "| Route | Strategy | First call/s | Warm median/s | Admissible points |", "|---|---|---:|---:|---:|"]
    for route in ROUTES:
        group = {}
        for strategy in STRATEGIES:
            file = output / f"{route}_{strategy}.json"
            if not file.exists():
                continue
            data = json.loads(file.read_text(encoding="utf-8"))
            rounds = data["rounds"]
            if not rounds:
                continue
            first = rounds[0]
            warm = [r["seconds"] for r in rounds[1:]]
            group[strategy] = dict(first_seconds=first["seconds"], warm_median=statistics.median(warm) if warm else None,
                repetitions=len(rounds), accepted=sum(admissible(r) for r in first["rows"]), total=len(first["rows"]),
                repeat_checks=[compare_rows(first["rows"],r["rows"]) for r in rounds[1:]],
                complete_grids=all(complete_round(data,r) for r in rounds),
                phase_pipeline_check=data.get("phase_pipeline_check"),
                curves=first["curves"], source_sha=data["source_sha"], config=data["config"],
                actions_run_id=data.get("actions_run_id"))
            lines.append(f"|{route}|{strategy}|{first['seconds']:.3f}|{statistics.median(warm) if warm else float('nan'):.3f}|{group[strategy]['accepted']}/{group[strategy]['total']}|")
        comparisons = {}
        for other in ("legacy", "single_pass"):
            if other in group and "palc" in group:
                left = json.loads((output / f"{route}_{other}.json").read_text())["rounds"][0]["rows"]
                right = json.loads((output / f"{route}_palc.json").read_text())["rounds"][0]["rows"]
                comparisons[other+"_vs_palc"] = compare_rows(left,right)
        result[route] = dict(methods=group,comparisons=comparisons)
    (output/"summary.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    for route, group in result.items():
        for name, comparison in group["comparisons"].items():
            lines.append(f"\n- {route} {name}: {len(comparison['coverage'])} coverage differences, {len(comparison['mismatches'])} numerical differences; maxima={comparison['maxima']}.")
    (output/"summary.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    print("\n".join(lines),flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--profile",choices=("quick","medium"),default="medium")
    parser.add_argument("--budget-seconds",type=int,default=2800)
    parser.add_argument("--repetitions",type=int,default=2)
    parser.add_argument("--cases",type=parse_cases,default=ALL_CASES,
                        help="Ordered route:strategy pairs, comma separated; default: all nine")
    parser.add_argument("--summarize-only",action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    if args.summarize_only:
        summarize(args.output)
        return 0
    if not 1 <= args.budget_seconds <= 3000 or not 1 <= args.repetitions <= 3:
        parser.error("budget must be 1..3000 seconds and repetitions 1..3")
    start = time.monotonic()
    processes = []
    for route, strategy in args.cases:
        remaining = args.budget_seconds - (time.monotonic()-start)
        if remaining <= 30:
            processes.append(dict(route=route,strategy=strategy,status="budget_exhausted"))
            continue
        destination = (args.output/f"{route}_{strategy}.json").resolve()
        command = ["julia","--startup-file=no","--project=.",
                   "scripts/analysis/pnjl_bifurcation_spike/run_production_comparison.jl",
                   route,strategy,str(destination),args.profile,str(args.repetitions)]
        print(f"Starting {route}/{strategy}; shared budget {remaining:.0f}s",flush=True)
        began = time.monotonic()
        try:
            completed = subprocess.run(command,cwd=ROOT,timeout=remaining,check=False)
            record = dict(route=route,strategy=strategy,status="finished",returncode=completed.returncode)
        except subprocess.TimeoutExpired:
            record = dict(route=route,strategy=strategy,status="timeout",returncode=None)
        record["process_seconds"] = time.monotonic()-began
        processes.append(record)
        (args.output/"processes.json").write_text(json.dumps(dict(processes=processes,
            seconds=time.monotonic()-start,budget_seconds=args.budget_seconds,
            actions_run_id=os.environ.get("GITHUB_RUN_ID")),indent=2)+"\n",encoding="utf-8")
        if record.get("returncode") not in (0,None):
            summarize(args.output)
            return 1
    summary = summarize(args.output)
    complete = len(processes) == len(args.cases) and all(p.get("returncode")==0 for p in processes)
    # Legacy differences are evidence for review; new methods must agree before default promotion.
    problems = promotion_checks(summary, args.repetitions, args.cases)
    (args.output/"checks.json").write_text(json.dumps(dict(complete=complete,
        complete_nine_method_workload=complete and set(args.cases)==set(ALL_CASES),
        cases=args.cases, problems=problems),indent=2)+"\n",encoding="utf-8")
    return 0 if complete and not problems else 2


if __name__ == "__main__":
    sys.exit(main())
