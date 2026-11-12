"""Compare all target roots and selected states in a bounded path optimization run."""
import argparse
import json
from pathlib import Path
import statistics
import subprocess

from compare_constraints import ROOT, compare
from compare_endpoints import root_coverage


def compare_all_roots(reference, candidate):
    result = root_coverage(reference, candidate)
    maxima = dict(solution=0.0, pressure=0.0, rho_norm=0.0, entropy=0.0, energy=0.0)
    others = {row["coordinate"]: row["roots"] for row in candidate["all_roots"]}
    differences = []
    for row in reference["all_roots"]:
        for root in row["roots"]:
            matches = others.get(row["coordinate"], [])
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


def summarize(source, destination):
    data = json.loads(source.read_text(encoding="utf-8"))
    compact = {k: v for k, v in data.items() if k != "records"}
    compact["cases"] = {}
    lines = ["# PALC path optimization", "", "| Path | Method | Warm s | Baseline / method | Allocated MiB | All roots and selected |",
             "|---|---|---:|---:|---:|---|"]
    ok = True
    for name in dict.fromkeys(row["name"] for row in data["records"]):
        records = [r for r in data["records"] if r["name"] == name]
        methods = list(dict.fromkeys(r["method"] for r in records))
        baseline = next(r for r in records if r["method"] == "baseline")
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
                         if seconds is not None else f"| {name} | {method} | — | — | — | {good} |")
        compact["cases"][name] = dict(config=baseline["config"], reference_rows=baseline["rows"],
            reference_all_roots=baseline["all_roots"], methods=results)
    compact["all_checks_passed"] = ok
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(compact, ensure_ascii=False, separators=(",", ":")) + "\n", encoding="utf-8")
    source.with_name("summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines), flush=True)
    return ok


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--cases", default="all")
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
            str(source.resolve()), str(args.repetitions), args.cases]
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
