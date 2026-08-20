"""Collector and contract validator for Issue #130 C2 limited feasibility.

The numerical density shards are intentionally small, immutable fine-pool
producers.  This script validates their provenance, creates the scope plans for
CEP/crossover, and validates the solver-free Julia aggregate.  It never calls
the PNJL equilibrium solver.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import subprocess
from pathlib import Path
from typing import Any

SCHEMA_VERSION = "pnjl_c2_limited_feasibility_v1"
JOB_SCHEMA = "pnjl_c2_limited_feasibility_job_v1"
PLOT_SCHEMA_VERSION = "pnjl_c2_limited_feasibility_plot_v1"
CALCULATION_SHA = "4c9703c3be45b76608ab57d375082e29418bfd05"
XI_GRID = (-0.5, -0.35, -0.25, -0.2, -0.15, -0.1, 0.0, 0.3, 0.35, 0.5)
DENSITY_ANCHORS_BY_XI = {
    -0.5: (60.0, 160.0),
    -0.35: (51.0,),
    -0.25: (41.0,),
    -0.2: (41.0,),
    -0.15: (41.0,),
    -0.1: (41.0,),
    0.0: (51.0, 60.0, 145.0),
    0.3: (21.0,),
    0.35: (51.0, 101.0),
    0.5: (60.0, 120.0),
}
RHO_FINE_STEP = 0.003125
RHO_MAX = 4.0
RHO_COUNT = int(round(RHO_MAX / RHO_FINE_STEP)) + 1
DENSITY_ANCHORS = (
    (-0.35, 51.0), (-0.25, 41.0), (-0.2, 41.0), (-0.15, 41.0),
    (-0.1, 41.0), (0.0, 51.0), (0.3, 21.0), (0.35, 51.0),
    (0.35, 101.0),
)
FIRST_ORDER_CONTROLS = ((-0.5, 60.0), (0.0, 60.0), (0.5, 60.0))
MONOTONE_CONTROLS = ((-0.5, 160.0), (0.0, 145.0), (0.5, 120.0))
ALL_ANCHORS = DENSITY_ANCHORS + FIRST_ORDER_CONTROLS + MONOTONE_CONTROLS
ROUTES = ("stage_b_features_v1", "balanced_density_features_v2", "geometry_feedback_v2")
CAPS = (12, 16, 24)
CEP_XIS = (0.05, 0.15, 0.225, 0.35, 0.4, 0.45, 0.5)
CROSSOVER_XI = 0.2875
CROSSOVER_MU = (0.0, 148.7809455570, 272.7650668545, 297.7487748943, 322.3587153735)
RUNNER_COST_STOP_MINUTES = 150.0
RECOVERY_SCHEMA_VERSION = "pnjl_c2_limited_feasibility_recovery_v1"
RECOVERY_METHOD = "production_eval_materialization_v1"


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


def _row_key(row: dict[str, str]) -> tuple[float, float, int]:
    xi = float(row["xi"])
    temperature = float(row["T_MeV"])
    rho = float(row["rho"])
    rho_index = round(rho / RHO_FINE_STEP)
    if not math.isclose(rho, rho_index * RHO_FINE_STEP, abs_tol=3e-7, rel_tol=0.0):
        raise ValueError(f"rho is off the frozen grid: {row}")
    return round(xi, 8), round(temperature, 8), rho_index


def _production_eval_files(artifact_dir: Path) -> list[Path]:
    return sorted(artifact_dir.rglob("production_eval/*memoized*.csv"))


def _recovered_row(row: dict[str, str], expected_sha: str) -> dict[str, str]:
    required = ("xi", "T_MeV", "rho", "mu_avg_MeV", "residual_norm", "iterations", "converged")
    missing = [field for field in required if field not in row]
    if missing:
        raise ValueError(f"production eval missing fields {missing}")
    if not all(finite(row[field]) for field in required[:-1]):
        raise ValueError(f"non-finite production eval row: {row}")
    if row["converged"].lower() not in {"true", "1"}:
        raise ValueError(f"non-converged production eval row: {row}")
    return {
        "xi": row["xi"],
        "T_MeV": row["T_MeV"],
        "rho": row["rho"],
        "muq_MeV": row["mu_avg_MeV"],
        "residual_norm": row["residual_norm"],
        "iterations": row["iterations"],
        "converged": "true",
        "finite": "true",
        "sampling_role": RECOVERY_METHOD,
        "rho_level": "0",
        "calculation_sha": expected_sha,
    }


def _recover_manifest(manifest_path: Path, expected_sha: str,
                      expected_postprocess_sha: str, source_run_id: str,
                      recovery_postprocess_sha: str = "") -> dict[str, Any]:
    artifact_dir = manifest_path.parent
    original_manifest_bytes = manifest_path.read_bytes()
    manifest = read_json(manifest_path)
    if manifest.get("schema_version") != JOB_SCHEMA:
        raise ValueError(f"unexpected job schema: {manifest_path}")
    if manifest.get("scope") != "density":
        raise ValueError(f"unexpected job scope: {manifest_path}")
    if str(manifest.get("calculation_sha", "")).lower() != expected_sha.lower():
        raise ValueError(f"calculation SHA mismatch: {manifest_path}")
    if expected_postprocess_sha and manifest.get("postprocess_sha") != expected_postprocess_sha:
        raise ValueError(f"postprocess SHA mismatch: {manifest_path}")
    if str(manifest.get("source_run_id", source_run_id)) != str(source_run_id):
        raise ValueError(f"source run mismatch: {manifest_path}")
    if manifest.get("solver_called") is not True:
        raise ValueError(f"source artifact must record solver_called=true: {manifest_path}")
    xi = float(manifest["xi"])
    if xi not in XI_GRID:
        raise ValueError(f"unexpected xi in manifest: {manifest_path}")
    pool = _declared_file(manifest, manifest_path, "fine_pool.csv")
    rows = read_csv(pool)
    existing: dict[tuple[float, float, int], dict[str, str]] = {}
    for row in rows:
        key = _row_key(row)
        if key in existing:
            raise ValueError(f"duplicate source fine-pool key: {key}")
        existing[key] = row

    eval_files = _production_eval_files(artifact_dir)
    if not eval_files:
        raise ValueError(f"missing complete production eval source: {manifest_path}")
    production: dict[tuple[float, float, int], dict[str, str]] = {}
    eval_hashes: list[dict[str, str]] = []
    for eval_path in eval_files:
        eval_hashes.append({
            "path": eval_path.relative_to(artifact_dir).as_posix(),
            "sha256": sha256(eval_path),
        })
        for row in read_csv(eval_path):
            recovered = _recovered_row(row, expected_sha)
            key = _row_key(recovered)
            if key[0] != round(xi, 8):
                raise ValueError(f"production eval xi mismatch: {eval_path}")
            if key in production:
                raise ValueError(f"duplicate production eval key: {key}")
            production[key] = recovered

    expected_keys = {
        (round(xi, 8), round(temperature, 8), rho_index)
        for temperature in DENSITY_ANCHORS_BY_XI[xi]
        for rho_index in range(RHO_COUNT)
    }
    missing = sorted(expected_keys - set(existing))
    if not set(production).issubset(expected_keys):
        extra = sorted(set(production) - expected_keys)
        raise ValueError(f"production eval has unexpected keys: {extra[:3]}")
    unrecoverable = [key for key in missing if key not in production]
    if unrecoverable:
        raise ValueError(f"production eval cannot recover fine-pool keys: {unrecoverable[:3]}")
    recovered_rows = [production[key] for key in missing]
    merged_rows = list(existing.values()) + recovered_rows
    merged_rows.sort(key=lambda row: (float(row["T_MeV"]), -float(row["rho"])))
    fieldnames = (
        "xi", "T_MeV", "rho", "muq_MeV", "residual_norm", "iterations",
        "converged", "finite", "sampling_role", "rho_level", "calculation_sha",
    )
    with pool.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows({field: row[field] for field in fieldnames} for row in merged_rows)

    source_pool_sha = str(manifest["files"]["fine_pool.csv"])
    recovery = {
        "schema_version": RECOVERY_SCHEMA_VERSION,
        "method": RECOVERY_METHOD,
        "solver_called": False,
        "source_manifest_sha256": hashlib.sha256(original_manifest_bytes).hexdigest(),
        "source_fine_pool_sha256": source_pool_sha,
        "source_fine_pool_rows": len(rows),
        "recovered_fine_pool_rows": len(recovered_rows),
        "recovery_postprocess_sha": recovery_postprocess_sha,
        "production_eval_files": eval_hashes,
    }
    manifest["recovery"] = recovery
    manifest["files"]["fine_pool.csv"] = sha256(pool)
    manifest["recovery_overlay_sha256"] = sha256(pool)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return {
        "manifest": manifest_path.relative_to(manifest_path.parents[1]).as_posix(),
        "xi": xi,
        "source_rows": len(rows),
        "recovered_rows": len(recovered_rows),
        "final_rows": len(merged_rows),
        "source_manifest_sha256": recovery["source_manifest_sha256"],
        "recovered_fine_pool_sha256": sha256(pool),
    }


def recover_source(input_dir: Path, output_dir: Path, expected_sha: str,
                   expected_postprocess_sha: str, source_run_id: str,
                   recovery_postprocess_sha: str = "") -> dict[str, Any]:
    if not input_dir.is_dir():
        raise ValueError(f"missing source artifact directory: {input_dir}")
    if output_dir.exists():
        raise ValueError(f"recovery output already exists: {output_dir}")
    shutil.copytree(input_dir, output_dir)
    manifests = _manifest_paths(output_dir)
    if len(manifests) != len(XI_GRID):
        raise ValueError(f"expected {len(XI_GRID)} xi manifests, got {len(manifests)}")
    results = [
        _recover_manifest(
            path, expected_sha, expected_postprocess_sha, source_run_id,
            recovery_postprocess_sha,
        )
        for path in manifests
    ]
    if {float(result["xi"]) for result in results} != set(XI_GRID):
        raise ValueError("incomplete xi matrix after recovery")
    overlay = {
        "schema_version": RECOVERY_SCHEMA_VERSION,
        "method": RECOVERY_METHOD,
        "source_run_id": source_run_id,
        "source_calculation_sha": expected_sha,
        "source_postprocess_sha": expected_postprocess_sha,
        "recovery_postprocess_sha": recovery_postprocess_sha,
        "solver_called": False,
        "shards": results,
    }
    (output_dir / "recovery_manifest.json").write_text(
        json.dumps(overlay, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return overlay


def validate_source(input_dir: Path, expected_sha: str, expected_postprocess_sha: str,
                    source_run_id: str) -> dict[str, Any]:
    if not input_dir.is_dir():
        raise ValueError(f"missing source artifact directory: {input_dir}")
    manifests = _manifest_paths(input_dir)
    if len(manifests) != len(XI_GRID):
        raise ValueError(f"expected {len(XI_GRID)} xi manifests, got {len(manifests)}")
    seen_xi: set[float] = set()
    input_files: list[dict[str, str]] = []
    total_seconds = 0.0
    total_rows = 0
    failed_points = 0
    for manifest_path in manifests:
        manifest = read_json(manifest_path)
        if manifest.get("schema_version") != JOB_SCHEMA:
            raise ValueError(f"unexpected job schema: {manifest_path}")
        if manifest.get("scope") != "density":
            raise ValueError(f"unexpected job scope: {manifest_path}")
        if str(manifest.get("calculation_sha", "")).lower() != expected_sha.lower():
            raise ValueError(f"calculation SHA mismatch: {manifest_path}")
        if expected_postprocess_sha and manifest.get("postprocess_sha") != expected_postprocess_sha:
            raise ValueError(f"postprocess SHA mismatch: {manifest_path}")
        xi = float(manifest.get("xi"))
        if xi not in XI_GRID or xi in seen_xi:
            raise ValueError(f"duplicate or unexpected xi: {manifest_path}")
        seen_xi.add(xi)
        if manifest.get("solver_called") is not True:
            raise ValueError("source fine-pool job must record solver_called=true")
        pool = _declared_file(manifest, manifest_path, "fine_pool.csv")
        slices = _declared_file(manifest, manifest_path, "slice_metrics.csv")
        costs = _declared_file(manifest, manifest_path, "method_costs.csv")
        pool_rows = read_csv(pool)
        if not pool_rows:
            raise ValueError(f"empty fine pool: {pool}")
        keys: set[tuple[float, float, int]] = set()
        for row in pool_rows:
            row_xi = float(row["xi"])
            temperature = float(row["T_MeV"])
            rho = float(row["rho"])
            if not math.isclose(row_xi, xi, abs_tol=1e-8, rel_tol=0.0):
                raise ValueError(f"fine-pool xi mismatch: expected {xi}, got {row_xi}")
            if not any(math.isclose(temperature, expected, abs_tol=1e-8, rel_tol=0.0)
                       for expected in DENSITY_ANCHORS_BY_XI[xi]):
                raise ValueError(f"unexpected fine-pool temperature: {(xi, temperature)}")
            rho_index = round(rho / RHO_FINE_STEP)
            if not 0 <= rho_index <= round(RHO_MAX / RHO_FINE_STEP) or not math.isclose(
                rho, rho_index * RHO_FINE_STEP, abs_tol=3e-7, rel_tol=0.0
            ):
                raise ValueError(f"fine-pool rho is off the frozen grid: {(xi, temperature, rho)}")
            key = (round(row_xi, 8), round(temperature, 8), rho_index)
            if key in keys:
                raise ValueError(f"duplicate fine-pool key: {key}")
            keys.add(key)
            if row.get("converged", "").lower() not in {"true", "1"}:
                raise ValueError(f"non-converged fine-pool row: {key}")
            if row.get("finite", "").lower() not in {"true", "1"}:
                raise ValueError(f"non-finite fine-pool row: {key}")
            if not all(finite(row.get(field)) for field in ("xi", "T_MeV", "rho", "muq_MeV", "residual_norm")):
                raise ValueError(f"invalid fine-pool row: {key}")
        expected_keys = {
            (round(xi, 8), round(temperature, 8), rho_index)
            for temperature in DENSITY_ANCHORS_BY_XI[xi]
            for rho_index in range(RHO_COUNT)
        }
        if keys != expected_keys:
            missing = sorted(expected_keys - keys)
            extra = sorted(keys - expected_keys)
            raise ValueError(
                f"incomplete fine-pool grid for xi={xi}: rows={len(keys)}, "
                f"missing={missing[:3]}, extra={extra[:3]}"
            )
        for row in read_csv(slices):
            failed_points += int(float(row.get("failed_points", 0) or 0))
        for row in read_csv(costs):
            total_seconds += float(row.get("runner_seconds", 0) or 0)
        total_rows += len(pool_rows)
        for path in (pool, slices, costs, manifest_path):
            input_files.append({"path": path.relative_to(input_dir).as_posix(), "sha256": sha256(path)})
    if seen_xi != set(XI_GRID):
        raise ValueError(f"incomplete xi matrix: {sorted(seen_xi)}")
    return {
        "schema_version": SCHEMA_VERSION,
        "source_run_id": source_run_id,
        "source_job_count": len(manifests),
        "source_calculation_sha": expected_sha,
        "source_postprocess_sha": expected_postprocess_sha,
        "solver_called": True,
        "fine_pool_rows": total_rows,
        "failed_points": failed_points,
        "runner_minutes": total_seconds / 60.0,
        "input_files": input_files,
    }


def _audit_rows(audit_dir: Path, name: str) -> list[dict[str, str]]:
    path = audit_dir / "tables" / name
    if name == "c1_vs_c2_cep_gates.csv":
        frozen = audit_dir.parent.parent / "c2_followups" / "c2_limited_feasibility_v1" / "cep_failures.csv"
        if frozen.is_file():
            path = frozen
    if not path.is_file():
        raise ValueError(f"missing frozen audit table: {path}")
    return read_csv(path)


def build_scope_plan(scope: str, audit_dir: Path, output_dir: Path,
                     expected_sha: str, source_run_id: str) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    if scope == "cep":
        rows = _audit_rows(audit_dir, "c1_vs_c2_cep_gates.csv")
        if len(rows) != 17:
            raise ValueError(f"expected frozen 17 CEP brackets, got {len(rows)}")
        plan = {
            "scope": "cep", "status": "plan_ready",
            "source_audit": "c2_convergence_audit_v1",
            "frozen_bracket_count": len(rows), "midpoint_policy": "one_midpoint_per_bracket",
            "xi_0p225_extra_bisection": True, "max_new_temperature_slices": 18,
            "node_step_MeV": 0.0625, "window_extension_MeV": 0.25,
            "gate": {"bracket_width_MeV": 0.1, "low": "first_order", "high": "monotone"},
        }
    elif scope == "crossover":
        plan = {
            "scope": "crossover", "status": "plan_ready", "xi": CROSSOVER_XI,
            "mu_MeV": list(CROSSOVER_MU), "integral_configs": {
                "c1": {"rtol": 1e-10, "atol": 1e-12},
                "c2": {"rtol": 1e-8, "atol": 1e-10},
            },
            "diagnostic_levels": {"2": {"n_scan": 20}, "3": {"n_scan": 40}, "4": {"n_scan": 80}},
            "gate": {"T_MeV": 0.05, "density": 0.005, "response_relative": 0.025},
        }
    else:
        raise ValueError(f"unsupported scope plan: {scope}")
    plan.update({"schema_version": SCHEMA_VERSION, "calculation_sha": expected_sha,
                 "source_run_id": source_run_id, "solver_called": False})
    (output_dir / "scope_plan.json").write_text(json.dumps(plan, indent=2) + "\n", encoding="utf-8")
    return plan


def validate_aggregate(aggregate_dir: Path, expected_sha: str, source_run_id: str) -> dict[str, Any]:
    manifest_path = aggregate_dir / "manifest.json"
    manifest = read_json(manifest_path)
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"unexpected aggregate schema: {manifest_path}")
    if manifest.get("source_calculation_sha") != expected_sha:
        raise ValueError("aggregate calculation SHA mismatch")
    if str(manifest.get("source_run_id")) != str(source_run_id):
        raise ValueError("aggregate source run mismatch")
    if manifest.get("solver_called") is not False:
        raise ValueError("aggregate must record solver_called=false")
    for name in ("route_comparison.csv", "component_geometry.csv", "selected_point_index.csv",
                 "candidate_point_index.csv", "cost_frontier.csv", "claim_ledger.csv"):
        if not (aggregate_dir / name).is_file():
            raise ValueError(f"aggregate missing {name}")
    frontier = read_csv(aggregate_dir / "cost_frontier.csv")
    routes = {row.get("route") for row in frontier}
    if routes != set(ROUTES):
        raise ValueError(f"route frontier mismatch: {routes}")
    caps = {int(row["cap"]) for row in frontier}
    if caps != set(CAPS):
        raise ValueError(f"cap frontier mismatch: {caps}")
    for row in frontier:
        for field in ("unique_solves", "targeted_points", "dense_unique_solves"):
            if not finite(row.get(field)):
                raise ValueError(f"invalid cost field {field}: {row}")
    return manifest


def render_representative_plot(aggregate_dir: Path) -> list[dict[str, Any]]:
    """Create a lightweight figure and a provenance sidecar from replay rows.

    The plot is deliberately aggregate-only: no solver curve is read here.  A
    sidecar is written even when matplotlib is unavailable so the evidence
    package still records the attempted renderer and its input hashes.
    """
    input_names = (
        "route_comparison.csv",
        "component_geometry.csv",
        "selected_point_index.csv",
        "cost_frontier.csv",
    )
    input_hashes = {
        name: sha256(aggregate_dir / name)
        for name in input_names
        if (aggregate_dir / name).is_file()
    }
    script_path = Path(__file__).resolve()
    plot_manifest: dict[str, Any] = {
        "schema_version": PLOT_SCHEMA_VERSION,
        "generator": "scripts/analysis/pnjl_c2_limited_feasibility.py",
        "generator_sha256": sha256(script_path),
        "input_file_hashes": input_hashes,
        "figures": [],
        "renderer": "matplotlib",
    }
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        plot_manifest["renderer_status"] = "matplotlib_unavailable"
        (aggregate_dir / "plot_manifest.json").write_text(
            json.dumps(plot_manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return []
    rows = read_csv(aggregate_dir / "route_comparison.csv")
    if not rows:
        plot_manifest["renderer_status"] = "no_route_rows"
        (aggregate_dir / "plot_manifest.json").write_text(
            json.dumps(plot_manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return []
    figure_dir = aggregate_dir / "figures"
    figure_dir.mkdir(exist_ok=True)
    figure_path = figure_dir / "density_route_frontier.png"
    fig, ax = plt.subplots(figsize=(7, 4))
    for route in ROUTES:
        route_rows = [row for row in rows if row.get("route") == route and row.get("cap") == "12"]
        route_rows.sort(key=lambda row: (float(row["T_MeV"]), float(row["xi"])))
        if not route_rows:
            continue
        ax.plot([float(row["T_MeV"]) for row in route_rows],
                [float(row["density_error"]) if finite(row.get("density_error")) else math.nan for row in route_rows],
                marker=".", linewidth=0.8, label=route)
    ax.set_xlabel("T [MeV]")
    ax.set_ylabel("density geometry error")
    ax.set_title("C2 limited Stage-C route comparison (cap 12)")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(figure_path, dpi=150)
    plt.close(fig)
    plot_manifest["renderer_status"] = "rendered"
    plot_manifest["figures"] = [{
        "path": "figures/density_route_frontier.png",
        "sha256": sha256(figure_path),
        "role": "route_frontier_cap12",
    }]
    (aggregate_dir / "plot_manifest.json").write_text(
        json.dumps(plot_manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return plot_manifest["figures"]


def write_manifest_patch(aggregate_dir: Path, figures: list[dict[str, Any]]) -> None:
    path = aggregate_dir / "manifest.json"
    manifest = read_json(path)
    manifest["figures"] = figures
    plot_manifest_path = aggregate_dir / "plot_manifest.json"
    if not plot_manifest_path.is_file():
        raise ValueError("plot renderer did not create plot_manifest.json")
    manifest["plot_manifest_sha256"] = sha256(plot_manifest_path)
    files = dict(manifest.get("files", {}))
    for file in aggregate_dir.rglob("*"):
        if file.is_file() and file.name != "manifest.json":
            files[file.relative_to(aggregate_dir).as_posix()] = sha256(file)
    manifest["files"] = files
    path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scope", choices=("density", "cep", "crossover"), required=True)
    parser.add_argument("--mode", choices=("source_validate", "recover_source", "aggregate_validate", "scope_plan"), required=True)
    parser.add_argument("--input-dir", type=Path)
    parser.add_argument("--aggregate-dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--audit-dir", type=Path)
    parser.add_argument("--expected-calculation-sha", default=CALCULATION_SHA)
    parser.add_argument("--expected-source-postprocess-sha", default="")
    parser.add_argument("--recovery-postprocess-sha", default="")
    parser.add_argument("--source-run-id", required=True)
    args = parser.parse_args(argv)
    expected = args.expected_calculation_sha.lower()
    if len(expected) != 40 or any(ch not in "0123456789abcdef" for ch in expected):
        raise SystemExit("expected-calculation-sha must be lowercase 40-hex")
    if args.mode == "source_validate":
        if args.scope != "density" or args.input_dir is None:
            raise SystemExit("source_validate currently requires density and --input-dir")
        result = validate_source(args.input_dir, expected, args.expected_source_postprocess_sha, args.source_run_id)
        print(json.dumps(result, sort_keys=True))
        return 0 if result["failed_points"] == 0 and result["runner_minutes"] <= RUNNER_COST_STOP_MINUTES else 2
    if args.mode == "recover_source":
        if args.scope != "density" or args.input_dir is None or args.output_dir is None:
            raise SystemExit("recover_source requires density, --input-dir and --output-dir")
        result = recover_source(
            args.input_dir, args.output_dir, expected,
            args.expected_source_postprocess_sha, args.source_run_id,
            args.recovery_postprocess_sha,
        )
        print(json.dumps(result, sort_keys=True))
        return 0
    if args.mode == "scope_plan":
        if args.scope == "density":
            raise SystemExit("density scope requires numerical/source artifacts")
        if args.output_dir is None or args.audit_dir is None:
            raise SystemExit("scope_plan requires --output-dir and --audit-dir")
        print(json.dumps(build_scope_plan(args.scope, args.audit_dir, args.output_dir, expected, args.source_run_id), sort_keys=True))
        return 0
    if args.aggregate_dir is None:
        raise SystemExit("aggregate_validate requires --aggregate-dir")
    manifest = validate_aggregate(args.aggregate_dir, expected, args.source_run_id)
    figures = render_representative_plot(args.aggregate_dir)
    write_manifest_patch(args.aggregate_dir, figures)
    final_manifest = read_json(args.aggregate_dir / "manifest.json")
    plot_manifest_path = args.aggregate_dir / "plot_manifest.json"
    if final_manifest.get("plot_manifest_sha256") != sha256(plot_manifest_path):
        raise ValueError("plot manifest hash was not recorded in aggregate manifest")
    print(json.dumps({"verdict": manifest.get("verdict"), "figures": figures}, sort_keys=True))
    return 0 if manifest.get("verdict") == "density_feasible_candidate" else 2


if __name__ == "__main__":
    raise SystemExit(main())
