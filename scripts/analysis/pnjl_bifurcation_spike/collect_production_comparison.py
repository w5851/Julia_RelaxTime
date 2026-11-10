"""Collect completed production batches; compare values across runs, timings only within runs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from compare_production import ROUTES, STRATEGIES, admissible, compare_rows, complete_round


def collect(directories):
    batches, by_route = [], {route: [] for route in ROUTES}
    for directory in directories:
        process_file = directory / "processes.json"
        processes = json.loads(process_file.read_text(encoding="utf-8")) if process_file.exists() else None
        methods, raw = {}, {}
        for route in ROUTES:
            for strategy in STRATEGIES:
                file = directory / f"{route}_{strategy}.json"
                if not file.exists():
                    continue
                data = json.loads(file.read_text(encoding="utf-8"))
                rounds = data["rounds"]
                if not rounds:
                    continue
                first, last = rounds[0], rounds[-1]
                key = f"{route}:{strategy}"
                methods[key] = dict(source_sha=data["source_sha"], actions_run_id=data["actions_run_id"],
                    config=data["config"], julia_version=data["julia_version"],
                    bifurcationkit_version=data["bifurcationkit_version"],
                    seconds=[r["seconds"] for r in rounds], load_seconds=data["load_seconds"],
                    complete_grids=all(complete_round(data, r) for r in rounds),
                    accepted=[sum(admissible(p) for p in r["rows"]) for r in rounds],
                    total=[len(r["rows"]) for r in rounds],
                    repeat_checks=[compare_rows(first["rows"], r["rows"]) for r in rounds[1:]],
                    last_work=last["work"], last_path_work=last["single_pass_work"],
                    paths=first["paths"], curves=first["curves"],
                    phase_pipeline_check=data.get("phase_pipeline_check"))
                raw[key] = first["rows"]
                by_route[route].append((methods[key], strategy, first["rows"]))
        comparisons = {}
        for route in ROUTES:
            for strategy in ("legacy", "single_pass"):
                left, right = f"{route}:{strategy}", f"{route}:palc"
                if left in raw and right in raw:
                    comparison = compare_rows(raw[left], raw[right])
                    if min(len(methods[left]["seconds"]), len(methods[right]["seconds"])) >= 2:
                        comparison["warm_cost_ratio"] = methods[left]["seconds"][1] / methods[right]["seconds"][1]
                    comparisons[f"{left}_vs_palc"] = comparison
        batches.append(dict(actions_run_id=next(iter(methods.values()))["actions_run_id"] if methods else None,
                            processes=processes, methods=methods, same_runner_comparisons=comparisons))
    routes = {}
    for route, measurements in by_route.items():
        references = [measurement for measurement in measurements if measurement[1] == "palc"]
        if not references:
            continue
        reference, _, rows = references[0]
        comparisons, different_rows = {}, {}
        for method, strategy, other in measurements:
            key = f"{method['actions_run_id']}:{strategy}"
            comparison = compare_rows(rows, other)
            # Early batches did not serialize the fixed iterations=300 field.
            config_a, config_b = dict(reference["config"]), dict(method["config"])
            config_a.setdefault("iterations", 300)
            config_b.setdefault("iterations", 300)
            comparison["same_configuration"] = config_a == config_b
            comparison["source_sha"] = method["source_sha"]
            comparison["timing_comparable"] = reference["actions_run_id"] == method["actions_run_id"]
            comparisons[key] = comparison
            differing_keys = {tuple(item["key"]) for field in ("coverage", "mismatches") for item in comparison[field]}
            if differing_keys:
                different_rows[key] = [r for r in other if (r["T"], r["xi"], r["coordinate"]) in differing_keys]
        columns = sorted(rows[0])
        routes[route] = dict(reference_run_id=reference["actions_run_id"], source_sha=reference["source_sha"],
                             config=reference["config"], columns=columns,
                             rows=[[row.get(name) for name in columns] for row in rows],
                             comparisons=comparisons, different_rows=different_rows)
    return dict(schema="palc_production_delivery_v1", batches=batches, routes=routes,
                timing_contract="Ratios use the same Actions run; different runs are numerical parity checks only.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-dir", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = collect(args.batch_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, separators=(",", ":"), allow_nan=False) + "\n", encoding="utf-8")
    print(f"Saved {len(result['batches'])} batches and {len(result['routes'])} representative paths: {args.output}")


if __name__ == "__main__":
    main()
