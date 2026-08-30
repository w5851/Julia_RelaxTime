"""Validator and aggregate contract for the CEP limited-feasibility scope.

This module is deliberately solver-free.  Numerical artifacts are validated
for complete fine rho pools and provenance; the Julia evaluator owns the
production-parity classification and Maxwell replay.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

SCHEMA_VERSION = "pnjl_c2_cep_limited_feasibility_v2"
JOB_SCHEMA = "pnjl_c2_cep_limited_feasibility_job_v2"
CALCULATION_SHA = "4c9703c3be45b76608ab57d375082e29418bfd05"
TARGET_XI = (-0.45, -0.34375, -0.275, -0.21875, 0.025, 0.125, 0.15,
             0.2, 0.225, 0.35, 0.3625, 0.38125, 0.39375, 0.4, 0.4375,
             0.45, 0.5)
RHO_STEP = 0.003125
RHO_MAX = 4.0
RHO_COUNT = int(round(RHO_MAX / RHO_STEP)) + 1
WIDTH_GATE = 0.1
COST_STOP_MINUTES = 200.0
PLOT_SCHEMA_VERSION = "pnjl_c2_cep_limited_feasibility_plot_v2"


def parse_target_xi(value: str | None) -> tuple[float, ...]:
    """Resolve an explicit CEP subset without changing the default 17-point contract."""
    if value is None or not value.strip():
        return TARGET_XI
    try:
        requested = tuple(float(item.strip()) for item in value.split(",") if item.strip())
    except ValueError as error:
        raise ValueError(f"invalid --target-xi value: {value!r}") from error
    if not requested or len(set(requested)) != len(requested):
        raise ValueError("--target-xi must contain unique xi values")
    if any(xi not in TARGET_XI for xi in requested):
        raise ValueError(f"--target-xi contains an xi outside the frozen matrix: {requested}")
    return tuple(sorted(requested))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def finite(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _manifest_paths(input_dir: Path) -> list[Path]:
    return sorted(input_dir.rglob("manifest.json"))


def _declared_file(manifest: dict[str, Any], manifest_path: Path, name: str) -> Path:
    declared = manifest.get("files", {}).get(name)
    if not declared:
        raise ValueError(f"manifest missing files.{name}: {manifest_path}")
    path = manifest_path.parent / name
    if not path.is_file():
        raise ValueError(f"missing {name}: {manifest_path}")
    if str(declared) != sha256(path):
        raise ValueError(f"hash mismatch for {name}: {manifest_path}")
    return path


def _rho_index(value: str) -> int:
    rho = float(value)
    index = round(rho / RHO_STEP)
    if not 0 <= index <= round(RHO_MAX / RHO_STEP) or not math.isclose(
        rho, index * RHO_STEP, abs_tol=3e-7, rel_tol=0.0
    ):
        raise ValueError(f"rho is off the frozen grid: {value}")
    return index


def _validate_pool(path: Path, xi: float) -> tuple[set[tuple[float, int]], int]:
    rows = read_csv(path)
    keys: set[tuple[float, int]] = set()
    for row in rows:
        row_xi = float(row["xi"])
        T = float(row["T_MeV"])
        key = (round(T, 8), _rho_index(row["rho"]))
        if not math.isclose(row_xi, xi, abs_tol=1e-8, rel_tol=0.0):
            raise ValueError(f"fine-pool xi mismatch: {row}")
        if key in keys:
            raise ValueError(f"duplicate fine-pool key: {key}")
        keys.add(key)
        required = ("xi", "T_MeV", "rho", "muq_MeV", "residual_norm")
        if not all(finite(row.get(field)) for field in required):
            raise ValueError(f"non-finite fine-pool row: {row}")
        if row.get("converged", "").lower() not in {"true", "1"}:
            raise ValueError(f"non-converged fine-pool row: {row}")
        if row.get("finite", "").lower() not in {"true", "1"}:
            raise ValueError(f"non-finite fine-pool row: {row}")
    temperatures = {T for T, _ in keys}
    for T in temperatures:
        count = sum(1 for row_T, _ in keys if row_T == T)
        if count != RHO_COUNT:
            raise ValueError(f"incomplete fine-pool grid xi={xi}, T={T}: {count}")
    return keys, len(rows)


def _validate_costs(path: Path, xi: float) -> tuple[list[dict[str, str]], float, int]:
    rows = read_csv(path)
    if {row.get("method") for row in rows} != {"hybrid", "oracle"}:
        raise ValueError(f"method_costs must contain hybrid and oracle rows: {path}")
    total_seconds = 0.0
    total_unique = 0
    for row in rows:
        if row.get("xi") and not math.isclose(float(row["xi"]), xi, abs_tol=1e-8, rel_tol=0.0):
            raise ValueError(f"method_costs xi mismatch: {path}")
        fields = ("unique_solves", "point_requests", "cache_hits", "failed_points", "runner_seconds")
        if not all(finite(row.get(field)) for field in fields):
            raise ValueError(f"non-finite method cost row: {path}")
        unique = int(float(row["unique_solves"]))
        requests = int(float(row["point_requests"]))
        hits = int(float(row["cache_hits"]))
        failed = int(float(row["failed_points"]))
        if min(unique, requests, hits, failed) < 0 or requests != unique + hits:
            raise ValueError(f"method cost reconciliation failure: {path}")
        if failed != 0:
            raise ValueError(f"method cost records failed points: {path}")
        total_unique += unique
        total_seconds += float(row["runner_seconds"])
    return rows, total_seconds, total_unique


def _validate_route_provenance(manifest: dict[str, Any]) -> None:
    provenance = manifest.get("provenance", {})
    if provenance.get("reference_write") is not False:
        raise ValueError("source artifact must not write reference")
    if provenance.get("route_before_oracle_gate") is not True:
        raise ValueError("source route provenance missing")
    if provenance.get("oracle_labels_used_for_routing") is not False:
        raise ValueError("oracle routing leakage")


def validate_source(input_dir: Path, expected_sha: str, expected_postprocess_sha: str,
                    source_run_id: str,
                    target_xi: tuple[float, ...] = TARGET_XI) -> dict[str, Any]:
    manifests = _manifest_paths(input_dir)
    if len(manifests) != len(target_xi):
        raise ValueError(f"expected {len(target_xi)} CEP manifests, got {len(manifests)}")
    seen: set[float] = set()
    input_files: list[dict[str, str]] = []
    total_rows = 0
    total_seconds = 0.0
    for manifest_path in manifests:
        manifest = read_json(manifest_path)
        if manifest.get("schema_version") != JOB_SCHEMA:
            raise ValueError(f"unexpected job schema: {manifest_path}")
        if manifest.get("scope") != "cep":
            raise ValueError(f"unexpected job scope: {manifest_path}")
        xi = float(manifest.get("xi"))
        if xi not in target_xi or xi in seen:
            raise ValueError(f"duplicate/unexpected xi: {manifest_path}")
        seen.add(xi)
        if str(manifest.get("calculation_sha", "")).lower() != expected_sha.lower():
            raise ValueError(f"calculation SHA mismatch: {manifest_path}")
        if expected_postprocess_sha and manifest.get("postprocess_sha") != expected_postprocess_sha:
            raise ValueError(f"postprocess SHA mismatch: {manifest_path}")
        if str(manifest.get("source_run_id", source_run_id)) != str(source_run_id):
            raise ValueError(f"source run mismatch: {manifest_path}")
        if manifest.get("solver_called") is not True:
            raise ValueError(f"source artifact must record solver_called=true: {manifest_path}")
        try:
            _validate_route_provenance(manifest)
        except ValueError as error:
            raise ValueError(f"{error}: {manifest_path}") from error
        pool = _declared_file(manifest, manifest_path, "fine_pool.csv")
        slices = _declared_file(manifest, manifest_path, "slice_metrics.csv")
        costs = _declared_file(manifest, manifest_path, "method_costs.csv")
        _, rows = _validate_pool(pool, xi)
        slice_rows = read_csv(slices)
        if not slice_rows or len(slice_rows) > 2:
            raise ValueError(f"unexpected slice count for xi={xi}: {len(slice_rows)}")
        planned = {float(row["T_MeV"]) for row in slice_rows}
        pool_temperatures = {float(row["T_MeV"]) for row in read_csv(pool)}
        if planned != pool_temperatures:
            raise ValueError(f"slice/pool temperature mismatch for xi={xi}")
        cost_rows, cost_seconds, _ = _validate_costs(costs, xi)
        total_seconds += cost_seconds
        total_rows += rows
        for path in (pool, slices, costs, manifest_path):
            input_files.append({"path": path.relative_to(input_dir).as_posix(), "sha256": sha256(path)})
    if seen != set(target_xi):
        raise ValueError(f"incomplete xi matrix: {sorted(seen)}")
    return {
        "schema_version": SCHEMA_VERSION,
        "source_run_id": str(source_run_id),
        "source_job_count": len(manifests),
        "target_xi": list(target_xi),
        "source_calculation_sha": expected_sha,
        "source_postprocess_sha": expected_postprocess_sha,
        "solver_called": True,
        "fine_pool_rows": total_rows,
        "runner_minutes": total_seconds / 60.0,
        "input_files": input_files,
    }


def render_plot(aggregate_dir: Path) -> list[dict[str, str]]:
    rows = read_csv(aggregate_dir / "cep_bracket_results.csv")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []
    figure_dir = aggregate_dir / "figures"
    figure_dir.mkdir(exist_ok=True)
    path = figure_dir / "cep_bracket_widths.png"
    rows = sorted(rows, key=lambda row: float(row["xi"]))
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar([float(row["xi"]) for row in rows], [float(row["bracket_width_T_MeV"]) for row in rows], width=0.008)
    ax.axhline(WIDTH_GATE, color="black", linewidth=0.8)
    ax.set_xlabel(r"$\xi$")
    ax.set_ylabel("CEP bracket width [MeV]")
    ax.set_title("C2 limited CEP feasibility")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return [{"path": "figures/cep_bracket_widths.png", "sha256": sha256(path)}]


def validate_aggregate(aggregate_dir: Path, expected_sha: str, source_run_id: str,
                       target_xi: tuple[float, ...] = TARGET_XI,
                       expected_workflow_schema: str = "") -> dict[str, Any]:
    manifest_path = aggregate_dir / "manifest.json"
    manifest = read_json(manifest_path)
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"unexpected aggregate schema: {manifest_path}")
    if manifest.get("source_calculation_sha") != expected_sha:
        raise ValueError("aggregate calculation SHA mismatch")
    if str(manifest.get("source_run_id")) != str(source_run_id):
        raise ValueError("aggregate source run mismatch")
    if expected_workflow_schema and manifest.get("workflow_schema") != expected_workflow_schema:
        raise ValueError("aggregate workflow schema mismatch")
    if manifest.get("solver_called") is not False:
        raise ValueError("aggregate must record solver_called=false")
    if manifest.get("oracle_labels_used_for_routing") is not False:
        raise ValueError("oracle labels must not route CEP points")
    declared_files = manifest.get("files", {})
    if not isinstance(declared_files, dict) or not declared_files:
        raise ValueError("aggregate output file hashes are missing")
    for name, declared_hash in declared_files.items():
        path = aggregate_dir / name
        if not path.is_file():
            raise ValueError(f"aggregate declared file is missing: {name}")
        if str(declared_hash) != sha256(path):
            raise ValueError(f"aggregate output hash mismatch: {name}")
    for name in ("cep_bracket_results.csv", "temperature_states.csv", "method_costs.csv",
                 "cost_frontier.csv", "route_trace.csv", "claim_ledger.csv"):
        if not (aggregate_dir / name).is_file():
            raise ValueError(f"aggregate missing {name}")
    rows = read_csv(aggregate_dir / "cep_bracket_results.csv")
    if len(rows) != len(target_xi) or {float(row["xi"]) for row in rows} != set(target_xi):
        raise ValueError("aggregate CEP bracket matrix mismatch")
    cost_rows = read_csv(aggregate_dir / "cost_frontier.csv")
    if not cost_rows or not finite(cost_rows[0].get("runner_minutes")):
        raise ValueError("aggregate missing runner cost")
    cost = cost_rows[0]
    if cost.get("point_request_reconciliation", "").lower() not in {"true", "1"}:
        raise ValueError("aggregate cost reconciliation failed")
    if not finite(cost.get("unique_solves")) or not finite(cost.get("point_requests")) or not finite(cost.get("cache_hits")):
        raise ValueError("aggregate cost counters are not finite")
    if int(float(cost["point_requests"])) != int(float(cost["unique_solves"])) + int(float(cost["cache_hits"])):
        raise ValueError("aggregate cost counters do not reconcile")
    method_rows = read_csv(aggregate_dir / "method_costs.csv")
    if len(method_rows) != len(target_xi):
        raise ValueError("aggregate method-cost matrix mismatch")
    if not all(row.get("method") == "aggregate" for row in method_rows):
        raise ValueError("aggregate method-cost rows have unexpected method")
    return manifest


def _refresh_aggregate_file_hashes(aggregate_dir: Path, manifest: dict[str, Any]) -> None:
    files: dict[str, str] = {}
    for path in sorted(aggregate_dir.rglob("*")):
        if not path.is_file() or path.name == "manifest.json":
            continue
        files[path.relative_to(aggregate_dir).as_posix()] = sha256(path)
    manifest["files"] = files
    (aggregate_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("source_validate", "aggregate_validate"), required=True)
    parser.add_argument("--input-dir", type=Path)
    parser.add_argument("--aggregate-dir", type=Path)
    parser.add_argument("--expected-calculation-sha", default=CALCULATION_SHA)
    parser.add_argument("--expected-source-postprocess-sha", default="")
    parser.add_argument("--source-run-id", required=True)
    parser.add_argument(
        "--target-xi",
        default="",
        help="comma-separated subset of the frozen CEP xi matrix; default is all 17 points",
    )
    parser.add_argument("--expected-workflow-schema", default="")
    args = parser.parse_args(argv)
    expected = args.expected_calculation_sha.lower()
    if len(expected) != 40 or any(ch not in "0123456789abcdef" for ch in expected):
        raise SystemExit("expected-calculation-sha must be lowercase 40-hex")
    try:
        target_xi = parse_target_xi(args.target_xi)
    except ValueError as error:
        raise SystemExit(str(error)) from error
    if args.mode == "source_validate":
        if args.input_dir is None:
            raise SystemExit("source_validate requires --input-dir")
        result = validate_source(
            args.input_dir,
            expected,
            args.expected_source_postprocess_sha,
            args.source_run_id,
            target_xi,
        )
        print(json.dumps(result, sort_keys=True))
        return 0 if result["runner_minutes"] <= COST_STOP_MINUTES else 2
    if args.aggregate_dir is None:
        raise SystemExit("aggregate_validate requires --aggregate-dir")
    manifest = validate_aggregate(
        args.aggregate_dir,
        expected,
        args.source_run_id,
        target_xi,
        args.expected_workflow_schema,
    )
    figures = render_plot(args.aggregate_dir)
    (args.aggregate_dir / "plot_manifest.json").write_text(json.dumps({
        "schema_version": PLOT_SCHEMA_VERSION, "aggregate_schema_version": SCHEMA_VERSION,
        "workflow_schema": manifest.get("workflow_schema", ""),
        "source_calculation_sha": expected, "source_run_id": str(args.source_run_id),
        "figures": figures,
    }, indent=2) + "\n", encoding="utf-8")
    _refresh_aggregate_file_hashes(args.aggregate_dir, manifest)
    print(json.dumps({"verdict": manifest.get("verdict"), "figures": figures}, sort_keys=True))
    return 0 if manifest.get("verdict") == "cep_feasible_candidate" else 2


if __name__ == "__main__":
    raise SystemExit(main())
