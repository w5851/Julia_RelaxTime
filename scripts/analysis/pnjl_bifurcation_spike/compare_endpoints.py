"""Compare two endpoint natural continuation with catalog and single-anchor PALC."""
import argparse
import json
from pathlib import Path
import statistics
import subprocess

from compare_constraints import ROOT, admissible, compare


def root_coverage(reference, candidate):
    left = {r["coordinate"]: r["roots"] for r in reference["all_roots"]}
    right = {r["coordinate"]: r["roots"] for r in candidate["all_roots"]}
    missing, extra, invalid = [], [], []
    def close(a, b):
        return len(a["solution"]) == len(b["solution"]) and max(
            abs(x-y) for x,y in zip(a["solution"], b["solution"])) <= 1e-4
    for coordinate in sorted(set(left) | set(right)):
        a, b = left.get(coordinate, []), right.get(coordinate, [])
        invalid.extend(dict(coordinate=coordinate, side=side) for side, rows in (("reference",a),("candidate",b))
                       for row in rows if not admissible(row))
        missing.extend(dict(coordinate=coordinate, rho_norm=row["rho_norm"])
                       for row in a if not any(close(row, other) for other in b))
        extra.extend(dict(coordinate=coordinate, rho_norm=row["rho_norm"])
                     for row in b if not any(close(row, other) for other in a))
    return dict(same=not (missing or extra or invalid), missing=missing, extra=extra, invalid=invalid)


def summarize(output, collect=None):
    source = output / "endpoint_comparison.json"
    if not source.exists():
        (output/"summary.md").write_text("No completed endpoint results.\n",encoding="utf-8")
        return
    data = json.loads(source.read_text(encoding="utf-8"))
    compact = {k:v for k,v in data.items() if k != "records"}
    compact["cases"] = {}
    lines = ["# Two-endpoint and PALC comparison", "", "| Path | Method | Warm median / s | Selected agreement | Full candidate agreement |", "|---|---|---:|---|---|"]
    for name in dict.fromkeys(r["name"] for r in data["records"]):
        records = [r for r in data["records"] if r["name"] == name]
        methods = {}
        for method in dict.fromkeys(r["method"] for r in records):
            rounds = sorted((r for r in records if r["method"] == method),key=lambda r:r["repetition"])
            first = rounds[0]
            times = [r["seconds"] for r in rounds]
            methods[method] = dict(config=first["config"],rows=first["rows"],all_roots=first["all_roots"],
                timings_seconds=times, warm_median=statistics.median(times[1:]) if len(times)>1 else None,
                work=rounds[-1]["work"],diagnostics=first["diagnostics"],
                repeated_rows_exact=all(r["rows"] == first["rows"] for r in rounds),
                repeated_candidates_exact=all(r["all_roots"] == first["all_roots"] for r in rounds),
                repeat_checks=[dict(selected=compare(first,r),candidates=root_coverage(first,r)) for r in rounds[1:]])
        comparisons = {}
        if "catalog" in methods:
            for method,record in methods.items():
                result = dict(selected=compare(methods["catalog"],record),candidates=root_coverage(methods["catalog"],record))
                comparisons[method] = result
                warm = record["warm_median"]
                lines.append(f'| {name} | {method} | {warm:.6f} | {result["selected"]["same"]} | {result["candidates"]["same"]} |'
                             if warm is not None else f'| {name} | {method} | — | {result["selected"]["same"]} | {result["candidates"]["same"]} |')
        # Keep the reference states once, plus every method's measured cost and
        # comparison. Raw repeated state arrays remain in the Actions artifact.
        reference = methods.get("catalog", next(iter(methods.values())))
        retained_methods = {}
        for method, record in methods.items():
            retained_methods[method] = {k:v for k,v in record.items() if k not in ("rows","all_roots","config")}
            retained_methods[method]["candidate_counts"] = [len(row["roots"]) for row in record["all_roots"]]
        compact["cases"][name] = dict(config=reference["config"],
            reference_rows=reference["rows"],reference_all_roots=reference["all_roots"],
            methods=retained_methods,comparisons=comparisons)
    (output/"summary.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    target = collect or output/"endpoint_compact.json"
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(json.dumps(compact,ensure_ascii=False,separators=(",",":"))+"\n",encoding="utf-8")
    print("\n".join(lines),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--repetitions",type=int,default=3)
    parser.add_argument("--scope",choices=("endpoints","anchors","promotion"),default="endpoints")
    parser.add_argument("--cases",default="all",help="Comma separated path names, or all")
    parser.add_argument("--probe-anchors",action="store_true")
    parser.add_argument("--budget-seconds",type=int,default=900)
    parser.add_argument("--summarize-only",action="store_true")
    parser.add_argument("--collect",type=Path)
    args=parser.parse_args()
    if not 1 <= args.repetitions <= 3 or not 1 <= args.budget_seconds <= 1200:
        parser.error("repetitions must be 1..3 and budget must be 1..1200 seconds")
    args.output.mkdir(parents=True,exist_ok=True)
    code=0
    if not args.summarize_only and args.collect is None:
        cmd=["julia","--startup-file=no","--project=.","scripts/analysis/pnjl_bifurcation_spike/run_endpoint_comparison.jl",
             str((args.output/"endpoint_comparison.json").resolve()),str(args.repetitions),args.scope,
             args.cases,str(args.probe_anchors).lower()]
        try:
            code=subprocess.run(cmd,cwd=ROOT,timeout=args.budget_seconds,check=False).returncode
        except subprocess.TimeoutExpired:
            print("Endpoint comparison stopped at its shared budget.",flush=True)
            code=124
    summarize(args.output,args.collect)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
