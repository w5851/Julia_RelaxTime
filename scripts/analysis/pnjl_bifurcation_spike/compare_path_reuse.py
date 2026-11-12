"""Compare all target roots and selected states in a bounded path optimization run."""
import argparse
from collections import Counter
import json
from pathlib import Path
import statistics
import subprocess

from compare_constraints import ROOT, admissible, compare
from compare_endpoints import root_coverage


def compare_all_roots(reference, candidate):
    result = root_coverage(reference, candidate)
    maxima = dict(solution=0.0, pressure=0.0, rho_norm=0.0, entropy=0.0, energy=0.0)
    others = {row["coordinate"]: row["roots"] for row in candidate["all_roots"]}
    differences = []
    for row in reference["all_roots"]:
        if len(row["roots"]) != len(others.get(row["coordinate"], [])):
            differences.append(dict(coordinate=row["coordinate"], root_count=[
                len(row["roots"]), len(others.get(row["coordinate"], []))]))
        for root in row["roots"]:
            if not admissible(root):
                continue
            matches = [r for r in others.get(row["coordinate"], []) if admissible(r)]
            if not matches:
                continue
            match = min(matches, key=lambda other: max(abs(x-y) for x,y in zip(root["solution"],other["solution"])))
            delta = dict(solution=max(abs(x-y) for x,y in zip(root["solution"],match["solution"])))
            delta.update({key: abs(root[key]-match[key]) for key in maxima if key != "solution"})
            for key, value in delta.items():
                maxima[key] = max(maxima[key], value)
            if any(value > (1e-4 if key == "solution" else 1e-6) for key,value in delta.items()):
                differences.append(dict(coordinate=row["coordinate"], deltas=delta))
    result.update(maxima=maxima, differences=differences)
    result["same"] &= not differences
    return result


def check_completion(data):
    plan = data.get("expected_runs")
    if not plan:
        return dict(complete=False, reason="expected run plan missing")
    expected = {(row["name"], method, repetition) for row in plan
                for method in row["methods"] for repetition in range(1, row["repetitions"]+1)}
    actual = Counter((r["name"], r["method"], r["repetition"]) for r in data["records"])
    missing = sorted(expected-actual.keys())
    unexpected = sorted(actual.keys()-expected)
    duplicate = sorted(key for key, count in actual.items() if count != 1)
    return dict(complete=bool(expected) and not missing and not unexpected and not duplicate,
                expected_count=len(expected), completed_count=len(actual),
                missing=missing, unexpected=unexpected, duplicates=duplicate)


def summarize(source, destination):
    data = json.loads(source.read_text(encoding="utf-8"))
    compact = {k: v for k, v in data.items() if k != "records"}
    compact["cases"] = {}
    frozen_path = ROOT / "data/outputs/results/analysis/palc_matched_20261108/anchor_optimization.json"
    frozen = json.loads(frozen_path.read_text(encoding="utf-8")) if frozen_path.exists() else None
    lines = ["# PALC path optimization", "", "| Path | Method | Warm s | Baseline / method | Allocated MiB | All roots and selected |",
             "|---|---|---:|---:|---:|---|"]
    ok = bool(data["records"])
    for name in dict.fromkeys(row["name"] for row in data["records"]):
        records = [r for r in data["records"] if r["name"] == name]
        methods = list(dict.fromkeys(r["method"] for r in records))
        baseline = next((r for r in records if r["method"] == "baseline"), None)
        if baseline is None:
            ok = False
            lines.append(f"\n{name}: baseline missing; comparison incomplete.")
            continue
        base_rounds = [r for r in records if r["method"] == "baseline" and r["repetition"] > 1]
        base_time = statistics.median(r["seconds"] for r in base_rounds) if base_rounds else None
        results = {}
        for method in methods:
            rounds = [r for r in records if r["method"] == method]
            checks = [dict(repetition=r["repetition"], selected=compare(baseline, r),
                           candidates=compare_all_roots(baseline, r)) for r in rounds]
            good = all(c["selected"]["same"] and c["candidates"]["same"] for c in checks)
            ok &= good
            warm = [r for r in rounds if r["repetition"] > 1]
            seconds = statistics.median(r["seconds"] for r in warm) if warm else None
            memory = statistics.median(r["allocated_bytes"] for r in warm) if warm else None
            speedup = base_time / seconds if base_time and seconds else None
            results[method] = dict(timings_seconds=[r["seconds"] for r in rounds],
                allocated_bytes=[r["allocated_bytes"] for r in rounds],
                gc_seconds=[r["gc_seconds"] for r in rounds], warm_median=seconds,
                baseline_over_method=speedup, work=rounds[-1]["work"],
                diagnostics=rounds[-1]["diagnostics"], checks=checks)
            lines.append(f"| {name} | {method} | {seconds:.6f} | {speedup:.3f} | {memory / 2**20:.3f} | {good} |"
                         if seconds is not None and speedup is not None else f"| {name} | {method} | — | — | — | {good} |")
        root_rows = {r["coordinate"]:r["roots"] for r in baseline["all_roots"]}
        selection = []
        for row in baseline["rows"]:
            roots = root_rows.get(row["coordinate"], [])
            # Retain state arrays once; root_index is zero-based in this JSON.
            index = min(range(len(roots)), key=lambda i:max(abs(x-y) for x,y in
                zip(row["solution"],roots[i]["solution"]))) if admissible(row) and roots else None
            selection.append(dict(coordinate=row["coordinate"],root_index=index))
        compact["cases"][name] = dict(config=baseline["config"], reference_selection=selection,
            reference_all_roots=baseline["all_roots"], methods=results)
        if baseline.get("microbenchmark"):
            compact["cases"][name]["microbenchmark"] = baseline["microbenchmark"]
        if frozen is not None:
            for run_name in ("sigma_recovery", "promotion"):
                old = frozen["runs"][run_name]["cases"].get(name)
                if old is None or old["config"] != baseline["config"]:
                    continue
                root_rows = {r["coordinate"]: r["roots"] for r in old["reference_all_roots"]}
                selected_rows = [root_rows[r["coordinate"]][r["root_index"]] for r in old["reference_selection"]]
                reference = dict(config=old["config"], rows=selected_rows,
                    all_roots=old["reference_all_roots"])
                check = dict(source_sha=frozen["runs"][run_name]["source_sha"],
                    selected=compare(reference, baseline), candidates=compare_all_roots(reference, baseline))
                compact["cases"][name]["prior_reference_check"] = check
                ok &= check["selected"]["same"] and check["candidates"]["same"]
                break
    sequence_costs = {}
    for name, case in compact["cases"].items():
        if not name.startswith("seq_"):
            continue
        sequence = name.rsplit("_t", 1)[0]
        bucket = sequence_costs.setdefault(sequence, {})
        for method, result in case["methods"].items():
            if result["warm_median"] is None:
                continue
            totals = bucket.setdefault(method, dict(seconds=0.0, allocated_bytes=0.0,
                root_iterations=0, anchor_solves=0, reused_anchor_solves=0))
            totals["seconds"] += result["warm_median"]
            totals["allocated_bytes"] += statistics.median(result["allocated_bytes"][1:])
            for key in ("root_iterations", "anchor_solves", "reused_anchor_solves"):
                totals[key] += result["work"][key]
    compact["sequence_costs"] = sequence_costs
    if sequence_costs:
        lines += ["", "## Temperature sequence totals", "", "| Sequence | Method | Warm s | Baseline / method |",
                  "|---|---|---:|---:|"]
        for sequence, methods in sequence_costs.items():
            for method, cost in methods.items():
                cost["baseline_over_method"] = methods["baseline"]["seconds"] / cost["seconds"]
                lines.append(f'| {sequence} | {method} | {cost["seconds"]:.6f} | {cost["baseline_over_method"]:.3f} |')
    compact["numeric_checks_passed"] = ok
    compact["completion"] = check_completion(data)
    compact["all_checks_passed"] = ok and compact["completion"]["complete"]
    lines += ["", f"Numerical checks: {ok}; planned run complete: {compact['completion']['complete']}."]
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(compact, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    source.with_name("summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines), flush=True)
    return compact["all_checks_passed"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--cases", default="all")
    parser.add_argument("--methods", default="all")
    parser.add_argument("--budget-seconds", type=int, default=1200)
    parser.add_argument("--summarize-only", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.repetitions <= 5 or not 1 <= args.budget_seconds <= 3000:
        parser.error("repetitions must be 1..5 and budget must be 1..3000 seconds")
    args.output.mkdir(parents=True, exist_ok=True)
    source = args.source or args.output / "path_reuse.json"
    code = 0
    if not args.summarize_only:
        command = ["julia", "--startup-file=no", "--project=.",
            "scripts/analysis/pnjl_bifurcation_spike/run_path_reuse_comparison.jl",
            str(source.resolve()), str(args.repetitions), args.cases, args.methods]
        try:
            code = subprocess.run(command, cwd=ROOT, check=False, timeout=args.budget_seconds).returncode
        except subprocess.TimeoutExpired:
            code = 124
    if source.exists():
        ok = summarize(source, args.output / "path_reuse_compact.json")
        return code or (0 if ok else 1)
    return code or 1


if __name__ == "__main__":
    raise SystemExit(main())
