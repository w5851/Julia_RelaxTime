"""Compact completed Actions batches without retaining repeated state arrays."""
from __future__ import annotations

import argparse
import contextlib
import io
import json
from pathlib import Path

from summarize_matched_pilot import compare_points, summarize


def collect_case(directory):
    batches = []
    raw_batches = []
    reference = None
    environment = None
    for batch_dir in sorted(directory.glob("batch_*")):
        with contextlib.redirect_stdout(io.StringIO()):
            summary = summarize(batch_dir)
        if not summary["parity_passed"]:
            raise ValueError(f"unqualified batch: {batch_dir}")
        processes = json.loads((batch_dir / "processes.json").read_text(encoding="utf-8"))
        methods = {}
        raw_methods = {}
        for method in processes["expected_methods"]:
            data = json.loads((batch_dir / f"{method}.json").read_text(encoding="utf-8"))
            repetitions = data["repetitions"]
            raw_methods[method] = repetitions[-1]["result"]["points"]
            methods[method] = dict(load_seconds=data["load_seconds"],
                wall_seconds=[r["wall_seconds"] for r in repetitions],
                gc_seconds=[r["gc_seconds"] for r in repetitions],
                allocated_bytes=[r["allocated_bytes"] for r in repetitions],
                compile_seconds=[r.get("compile_seconds") for r in repetitions],
                stages=[r["result"].get("stages") for r in repetitions] if method == "palc" else None,
                last_work=repetitions[-1]["work"],
                work_identical=all(r["work"] == repetitions[-1]["work"] for r in repetitions),
                accepted_all=summary["measurements"][method]["accepted_all"])
            if reference is None and method == "palc":
                reference = repetitions[-1]["result"]
                environment = {key:data[key] for key in ("config","julia_version","bifurcationkit_version","threads","blas_threads")}
        checks = {}
        for name, check in summary["checks"].items():
            finite = [p for p in check["points"] if "state_distance" in p]
            checks[name] = dict(passed=check["passed"], count=len(check["points"]),
                max_state_distance=max((p["state_distance"] for p in finite), default=None),
                max_pressure_delta=max((abs(p["pressure_delta"]) for p in finite), default=None),
                differences=[p for p in check["points"] if not p["passed"]])
        batches.append(dict(batch=batch_dir.name, processes=processes,
            parity_passed=summary["parity_passed"], seed_over_palc_warm_speedup=summary["seed_over_palc_warm_speedup"],
            checks=checks, methods=methods))
        raw_batches.append(raw_methods)
    if not batches:
        raise ValueError(f"no batches: {directory}")
    independent = {}
    for index, batch in enumerate(batches[1:], start=1):
        for method in batch["methods"]:
            if method == "legacy_seed":
                continue
            result = compare_points(raw_batches[0][method], raw_batches[index][method])
            if not result["passed"]:
                raise ValueError(f"independent-process mismatch: {directory}/{method}")
            independent[f"{batch['batch']}:{method}"] = dict(passed=True,
                max_state_distance=max(p["state_distance"] for p in result["points"]),
                max_pressure_delta=max(abs(p["pressure_delta"]) for p in result["points"]))
    def rows(points):
        if not points:
            return dict(columns=[], rows=[])
        columns = sorted(points[0])
        return dict(columns=columns, rows=[[p.get(k) for k in columns] for p in points])
    # Keep one qualified numerical path and its branch corrections; equivalent
    # method/repetition states are represented by the actual comparison errors.
    for key in ("points", "anchors", "candidates"):
        reference[key] = rows(reference[key])
    legacy = raw_batches[0].get("legacy_seed", [])
    differences = []
    for point in legacy:
        comparison = compare_points([point], raw_batches[0]["palc"])
        if not comparison["passed"]:
            differences.append(point)
    return dict(environment=environment, batches=batches, numerical_reference=reference,
        legacy_different_points=rows(differences), independent_process_checks=independent)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-dir", action="append", type=Path, required=True)
    parser.add_argument("--candidate-checks", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cases = {}
    for directory in args.case_dir:
        case = collect_case(directory)
        identity = case["batches"][0]["processes"]
        cases[f"{identity['profile']}:{identity['case']}"] = case
    result = dict(schema="palc_single_pass_delivery_v1", cases=cases,
        candidate_checks=json.loads(args.candidate_checks.read_text(encoding="utf-8")),
        interpretation="same-discretization homogeneous root tracking; first-order coexistence requires separate thermodynamic construction")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, separators=(",", ":"), ensure_ascii=False)+"\n", encoding="utf-8")
    print(f"Saved {len(cases)} cases, {args.output.stat().st_size} bytes: {args.output}")


if __name__ == "__main__":
    main()
