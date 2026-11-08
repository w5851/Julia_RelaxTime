"""Compare equal-density outputs and costs; incomplete work never yields a speedup."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import statistics

METHODS = ("palc", "seed", "multi", "production", "dense")


def finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def compare_points(left, right, *, state_tol=1e-4, pressure_tol=1e-6):
    rows = []
    for a in left:
        matches = [b for b in right if abs(a["rho"] - b["rho"]) <= 1e-12]
        if len(matches) != 1:
            rows.append(dict(rho=a["rho"], passed=False, reason="missing_or_duplicate_target"))
            continue
        b = matches[0]
        valid = a["accepted"] and b["accepted"] and finite(a["pressure"]) and finite(b["pressure"])
        valid = valid and all(finite(x) for x in a["solution"] + b["solution"])
        if not valid:
            rows.append(dict(rho=a["rho"], passed=False, reason="unaccepted_point"))
            continue
        distance = math.dist(a["solution"], b["solution"])
        pressure_delta = a["pressure"] - b["pressure"]
        rows.append(dict(rho=a["rho"], passed=distance <= state_tol and abs(pressure_delta) <= pressure_tol,
                         state_distance=distance, pressure_delta=pressure_delta))
    return dict(passed=bool(rows) and all(r["passed"] for r in rows), points=rows)


def summarize(directory):
    directory.mkdir(parents=True, exist_ok=True)
    measurements = {}
    cases = {}
    for method in METHODS:
        path = directory / f"{method}.json"
        if not path.exists():
            measurements[method] = dict(complete=False, reason="no_output")
            continue
        data = json.loads(path.read_text(encoding="utf-8"))
        reps = data["repetitions"]
        warm = [r for r in reps if r["phase"] == "warm"]
        cases[method] = data
        expected_count = 9 if method == "dense" else 5
        measurements[method] = dict(complete=len(reps) == 3, source_sha=data["source_sha"],
            load_seconds=data["load_seconds"], first_call_seconds=reps[0]["wall_seconds"],
            warm_seconds=[r["wall_seconds"] for r in warm],
            warm_median_seconds=statistics.median(r["wall_seconds"] for r in warm) if warm else None,
            accepted_all=all(len(r["result"]["points"]) == expected_count and
                all(p["accepted"] for p in r["result"]["points"]) and not r["result"]["errors"] for r in reps),
            last_work=reps[-1]["work"])
    checks = {}
    for method, data in cases.items():
        reps = data["repetitions"]
        for rep in reps[1:]:
            checks[f"{method}_repeat_{rep['repetition']}"] = compare_points(
                rep["result"]["points"], reps[0]["result"]["points"])
    for other in ("production", "multi", "dense", "palc"):
        if "seed" in cases and other in cases:
            seed = cases["seed"]["repetitions"][-1]["result"]["points"]
            candidate = cases[other]["repetitions"][-1]["result"]["points"]
            checks[f"seed_vs_{other}"] = compare_points(seed, candidate)
    complete = len(cases) == 5 and all(m.get("complete") and m.get("accepted_all") for m in measurements.values())
    identities_match = len({d["source_sha"] for d in cases.values()}) == 1 and len({
        json.dumps(d["config"], sort_keys=True) for d in cases.values()}) == 1
    parity = complete and identities_match and all(check["passed"] for check in checks.values())
    speedup = None
    if parity:
        speedup = measurements["seed"]["warm_median_seconds"] / measurements["palc"]["warm_median_seconds"]
    result = dict(schema="palc_matched_summary_v1", complete=complete, identities_match=identities_match,
        parity_passed=parity, seed_over_palc_warm_speedup=speedup, measurements=measurements, checks=checks,
        interpretation="same_discretization_local_pilot; dense is half-step seed continuation, not an independent physical oracle")
    (directory / "summary.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    lines = ["# PALC matched pilot", "", f"Complete: {complete}; parity: {parity}", "",
             "| Method | Load / s | First call / s | Warm median / s | Accepted |",
             "| --- | ---: | ---: | ---: | --- |"]
    for method, m in measurements.items():
        lines.append(f"| {method} | {m.get('load_seconds')} | {m.get('first_call_seconds')} | "
                     f"{m.get('warm_median_seconds')} | {m.get('accepted_all', False)} |")
    lines += ["", f"Seed / PALC warm cost ratio (only if parity passes): {speedup}", "",
              "The 8 x 4 quadrature is a solver diagnostic, not an integration-convergence claim.",
              "All methods start from fresh seeds on every repetition. Load is separate; first call includes remaining JIT.",
              "PALC includes anchor discovery and correction to the same target densities.",
              "Production parity covers the equilibrium stage when downstream meson steps succeed; no meson densities are run."]
    (directory / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    arguments = parser.parse_args()
    summarize(arguments.directory)
