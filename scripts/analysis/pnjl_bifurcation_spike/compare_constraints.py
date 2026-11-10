"""Small constraint-path and single-anchor comparison; no full phase production."""
import argparse
import json
import math
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[3]


def admissible(row):
    return (row.get("converged", False)
            and all(row.get(key) is not None and math.isfinite(row[key])
                    for key in ("residual", "constraint_error", "pressure", "rho_norm", "entropy", "energy"))
            and row["residual"] <= 1e-6 and row["constraint_error"] <= 1e-6
            and len(row.get("solution", [])) in (5, 8)
            and len(row.get("mu_vec", [])) == len(row.get("masses", [])) == 3
            and all(x is not None and math.isfinite(x) for key in ("solution", "mu_vec", "masses") for x in row[key])
            and min(row["masses"]) > 0 and row["masses"][2] > row["masses"][0]
            and all(0 <= x <= 1 for x in row["solution"][3:5]))


def compare(left, right):
    a, b = ({row["coordinate"]: row for row in data["rows"]} for data in (left, right))
    duplicate_keys = len(a) != len(left["rows"]) or len(b) != len(right["rows"])
    same_config = left.get("config") == right.get("config")
    coverage, differences = [], []
    maxima = dict(solution=0.0, mu_vec=0.0, pressure=0.0, rho_norm=0.0, entropy=0.0, energy=0.0)
    for key in sorted(set(a) | set(b)):
        if key not in a or key not in b or not (admissible(a[key]) and admissible(b[key])):
            coverage.append(key)
            continue
        deltas = {k: max(abs(x-y) for x,y in zip(a[key][k], b[key][k]))
                  if isinstance(a[key][k], list) else abs(a[key][k]-b[key][k]) for k in maxima}
        for k,v in deltas.items():
            maxima[k] = max(maxima[k], v)
        if any(v > (1e-4 if k in ("solution", "mu_vec") else 1e-6) for k,v in deltas.items()):
            differences.append(dict(coordinate=key, deltas=deltas,
                                    left_mu=a[key]["mu_vec"], right_mu=b[key]["mu_vec"]))
    return dict(same=bool(a) and not duplicate_keys and same_config and not coverage and not differences,
                coverage=coverage, differences=differences, maxima=maxima,
                duplicate_keys=duplicate_keys, same_config=same_config)


def collect_results(source, destination):
    data = json.loads(source.read_text(encoding="utf-8"))
    compact = {k:v for k,v in data.items() if k != "records"}
    compact["cases"] = {}
    for name in dict.fromkeys(r["name"] for r in data["records"]):
        records = [r for r in data["records"] if r["name"] == name]
        methods = {}
        for method in dict.fromkeys(r["method"] for r in records):
            rounds = [r for r in records if r["method"] == method]
            first = rounds[0]
            methods[method] = dict(rows=first["rows"], config=first["config"],
                timings_seconds=[r["seconds"] for r in rounds], work=rounds[-1]["work"],
                diagnostics=first["diagnostics"], candidate_counts=first["candidate_counts"],
                repetitions=len(rounds), repeated_rows_exact=all(r["rows"] == first["rows"] for r in rounds),
                repeat_checks=[compare(first,r) for r in rounds[1:]])
        baseline = "single_seed" if name.startswith("tmu_") else "legacy"
        comparisons = compare(methods[baseline], methods["palc"]) if baseline in methods and "palc" in methods else None
        compact["cases"][name] = dict(methods=methods, comparison=comparisons)
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(json.dumps(compact,ensure_ascii=False,separators=(",",":"))+"\n",encoding="utf-8")
    return compact


def summarize(output):
    source = output / "constraint_comparison.json"
    if not source.exists():
        (output / "summary.md").write_text("No completed constraint results.\n", encoding="utf-8")
        return
    data = json.loads(source.read_text(encoding="utf-8"))
    lines = ["# Constraint PALC comparison", "", "| Path | Method | Round | Seconds | Admissible |",
             "|---|---|---:|---:|---:|"]
    for item in data["records"]:
        good = sum(admissible(row) for row in item["rows"])
        lines.append(f'|{item["name"]}|{item["method"]}|{item["repetition"]}|{item["seconds"]:.3f}|{good}/{len(item["rows"])}|')
    comparisons = {}
    for name in dict.fromkeys(r["name"] for r in data["records"]):
        by_key = {(r["method"], r["repetition"]): r for r in data["records"] if r["name"] == name}
        baseline = "single_seed" if name.startswith("tmu_") else "legacy"
        for repetition in (1,2,3):
            a,b = by_key.get((baseline,repetition)), by_key.get(("palc",repetition))
            if a and b:
                result = compare(a,b)
                result["baseline_over_palc"] = a["seconds"]/b["seconds"]
                comparisons[f"{name}_{repetition}"] = result
                lines.append(f'\n{name} round {repetition}: same={result["same"]}, coverage={len(result["coverage"])}, differences={len(result["differences"])}, time ratio={result["baseline_over_palc"]:.3f}.')
    (output / "summary.json").write_text(json.dumps(dict(source_sha=data["source_sha"],comparisons=comparisons),indent=2)+"\n",encoding="utf-8")
    (output / "summary.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    print("\n".join(lines),flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--budget-seconds",type=int,default=900)
    parser.add_argument("--repetitions",type=int,default=2)
    parser.add_argument("--cases",choices=("all","entropy","sigma","bqs","single_seed"),default="all")
    parser.add_argument("--summarize-only",action="store_true")
    parser.add_argument("--collect",type=Path,help="Compact completed results to one retained JSON without solving")
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    if args.collect:
        collect_results(args.output/"constraint_comparison.json",args.collect)
        return 0
    code=0
    if not args.summarize_only:
        if not 1 <= args.budget_seconds <= 1200 or not 1 <= args.repetitions <= 3:
            parser.error("constraint budget must be 1..1200 s and repetitions 1..3")
        command=["julia","--startup-file=no","--project=.",
                 "scripts/analysis/pnjl_bifurcation_spike/run_constraint_comparison.jl",
                 str((args.output/"constraint_comparison.json").resolve()),str(args.repetitions),args.cases]
        try:
            code=subprocess.run(command,cwd=ROOT,timeout=args.budget_seconds,check=False).returncode
        except subprocess.TimeoutExpired:
            print("Constraint comparison stopped at the shared budget; retain completed cases.",flush=True)
            code=124
    summarize(args.output)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
