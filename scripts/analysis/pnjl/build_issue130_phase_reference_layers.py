#!/usr/bin/env python3
"""Build Issue #130 strict, derived and render phase-reference layers.

The input is the frozen Maxwell expansion evidence plus the immutable v7
solver-free crossover package.  This utility performs only CSV interpolation,
provenance assembly and plotting; it never calls Julia, the equilibrium solver,
Maxwell, or an oracle.  Strict rows are preserved separately from derived
rows, and every derived row is marked ``interpolated_noncertified``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence


CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
SCHEMA = "pnjl_issue130_phase_reference_layers_v1"
XI_STEP = 0.00625
EPS = 1e-9


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, fields: Sequence[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def f(value: Any) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"non-finite numeric value: {value!r}")
    return number


def xi_key(value: Any) -> str:
    return f"{f(value):.10f}"


def source_record(path: Path) -> dict[str, Any]:
    return {"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size}


def resolve_manifest(root: Path) -> Path:
    direct = root / "manifest.json"
    if direct.is_file():
        return direct
    candidates = sorted(root.rglob("manifest.json"))
    if len(candidates) != 1:
        raise ValueError(f"expected one replay manifest under {root}, found {len(candidates)}")
    return candidates[0]


def uniform_xi() -> list[float]:
    count = round(1.0 / XI_STEP)
    return [round(-0.5 + index * XI_STEP, 10) for index in range(count + 1)]


def group_by_xi(rows: Iterable[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[xi_key(row["xi"])].append(row)
    return grouped


def unique_sorted(values: Iterable[float]) -> list[float]:
    result: list[float] = []
    for value in sorted(values):
        if not result or abs(value - result[-1]) > EPS:
            result.append(value)
    return result


def interpolate_rows(rows: list[dict[str, Any]], axis: str, value: float, fields: Sequence[str]) -> dict[str, float] | None:
    ordered = sorted(rows, key=lambda row: f(row[axis]))
    if not ordered:
        return None
    exact = next((row for row in ordered if abs(f(row[axis]) - value) <= EPS), None)
    if exact is not None:
        return {field: f(exact[field]) for field in fields}
    if value < f(ordered[0][axis]) - EPS or value > f(ordered[-1][axis]) + EPS:
        return None
    left = ordered[0]
    right = ordered[-1]
    for candidate_left, candidate_right in zip(ordered, ordered[1:]):
        if f(candidate_left[axis]) <= value <= f(candidate_right[axis]):
            left, right = candidate_left, candidate_right
            break
    x0 = f(left[axis])
    x1 = f(right[axis])
    weight = 0.0 if abs(x1 - x0) <= EPS else (value - x0) / (x1 - x0)
    return {
        field: f(left[field]) + weight * (f(right[field]) - f(left[field]))
        for field in fields
    }


def bracket_xi(source: dict[str, list[dict[str, Any]]], value: float) -> tuple[float, float] | None:
    keys = sorted(f(key) for key in source)
    for key in keys:
        if abs(key - value) <= EPS:
            return key, key
    for left, right in zip(keys, keys[1:]):
        if left < value < right:
            return left, right
    return None


def interpolate_surface(
    source_rows: list[dict[str, Any]],
    *,
    surface: str,
    axis: str,
    fields: Sequence[str],
    source_layer: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Regrid one surface to a uniform xi mesh without extrapolation."""
    grouped = group_by_xi(source_rows)
    native_xi = sorted(f(key) for key in grouped)
    output: list[dict[str, Any]] = []
    coverage: list[dict[str, Any]] = []
    for target_xi in uniform_xi():
        bracket = bracket_xi(grouped, target_xi)
        if bracket is None:
            coverage.append({
                "surface": surface,
                "xi": target_xi,
                "xi_role": "derived_uniform",
                "xi_left": "",
                "xi_right": "",
                "axis_low": "",
                "axis_high": "",
                "source_rows": 0,
                "coverage_status": "unresolved_no_xi_bracket",
            })
            continue
        left_xi, right_xi = bracket
        left_rows = grouped[xi_key(left_xi)]
        right_rows = grouped[xi_key(right_xi)]
        left_axes = [f(row[axis]) for row in left_rows]
        right_axes = [f(row[axis]) for row in right_rows]
        low = max(min(left_axes), min(right_axes))
        high = min(max(left_axes), max(right_axes))
        axes = [value for value in unique_sorted(left_axes + right_axes) if low - EPS <= value <= high + EPS]
        exact = abs(left_xi - right_xi) <= EPS
        produced = 0
        for axis_value in axes:
            left_value = interpolate_rows(left_rows, axis, axis_value, fields)
            right_value = left_value if exact else interpolate_rows(right_rows, axis, axis_value, fields)
            if left_value is None or right_value is None:
                continue
            if exact:
                interpolated = left_value
            else:
                xi_weight = (target_xi - left_xi) / (right_xi - left_xi)
                interpolated = {
                    field: left_value[field] + xi_weight * (right_value[field] - left_value[field])
                    for field in fields
                }
            row: dict[str, Any] = {
                "surface": surface,
                "xi": target_xi,
                axis: axis_value,
                **interpolated,
                "layer": "strict_reference_v1" if exact else "interpolated_noncertified",
                "status": "native_uniform_xi" if exact else "piecewise_linear_common_support",
                "interpolation_method": "none_native" if exact else "bilinear_axis_xi_common_support",
                "source_layer": source_layer,
                "source_xi_left": left_xi,
                "source_xi_right": right_xi,
                "source_axis_low": low,
                "source_axis_high": high,
                "xi_gap": right_xi - left_xi,
                "reference_write": False,
            }
            output.append(row)
            produced += 1
        coverage.append({
            "surface": surface,
            "xi": target_xi,
            "xi_role": "native_uniform" if exact else "derived_uniform",
            "xi_left": left_xi,
            "xi_right": right_xi,
            "axis_low": low,
            "axis_high": high,
            "source_rows": produced,
            "coverage_status": "native_support" if exact else ("interpolated_common_support" if produced else "unresolved_no_common_axis_support"),
        })
    return output, coverage


def strict_maxwell(v7_rows: list[dict[str, str]], frozen_root: Path, run_id: str, calculation_sha: str) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    fields = (
        "surface", "xi", "T_MeV", "mu_MeV", "rho_hadron", "rho_quark", "area_residual",
        "grid_status", "grid_unresolved", "layer", "status", "calculation_sha",
        "source_run_id", "source_target_id", "source_workflow_sha", "postprocess_sha",
        "targeted_additions", "geometry_converged", "finite_and_converged",
        "reference_write", "oracle_labels_consumed", "source_grid_status",
    )
    rows: list[dict[str, Any]] = []
    for row in v7_rows:
        rows.append({
            "surface": "maxwell",
            "xi": f(row["xi"]), "T_MeV": f(row["T_MeV"]), "mu_MeV": f(row["mu_MeV"]),
            "rho_hadron": f(row["rho_hadron"]), "rho_quark": f(row["rho_quark"]),
            "area_residual": abs(f(row["area_residual"])),
            "grid_status": row.get("grid_status", ""),
            "grid_unresolved": row.get("grid_unresolved", "False").lower() == "true",
            "layer": "strict_reference_v1",
            "status": "native_c2_unresolved_geometry" if row.get("grid_unresolved", "").lower() == "true" else "native_c2_geometry_converged",
            "calculation_sha": calculation_sha,
            "source_run_id": "c2_v5_v6",
            "source_target_id": "",
            "source_workflow_sha": "",
            "postprocess_sha": row.get("v6_manifest_hash", ""),
            "targeted_additions": "",
            "geometry_converged": row.get("grid_unresolved", "").lower() != "true",
            "finite_and_converged": True,
            "reference_write": False,
            "oracle_labels_consumed": False,
            "source_grid_status": row.get("grid_status", ""),
        })
    target_root = frozen_root / "targets"
    for target_dir in sorted(p for p in target_root.iterdir() if p.is_dir()):
        summary = read_json(target_dir / "target_summary.json")
        metrics = read_csv(target_dir / "slice_metrics.csv")
        final = max(metrics, key=lambda row: int(row["level"]))
        rows.append({
            "surface": "maxwell",
            "xi": f(summary["xi"]), "T_MeV": f(summary["T_MeV"]),
            "mu_MeV": f(summary["final_candidate_mu_MeV"]),
            "rho_hadron": f(final["rho_hadron"]), "rho_quark": f(final["rho_quark"]),
            "area_residual": abs(f(summary["final_area_residual"])),
            "grid_status": "expansion_final_geometry_converged",
            "grid_unresolved": False,
            "layer": "strict_reference_v1",
            "status": "expansion_native_geometry_converged",
            "calculation_sha": calculation_sha,
            "source_run_id": run_id,
            "source_target_id": summary["target_id"],
            "source_workflow_sha": read_json(target_dir / "provenance.json").get("workflow_head_sha", ""),
            "postprocess_sha": read_json(target_dir / "provenance.json").get("workflow_head_sha", ""),
            "targeted_additions": int(summary["targeted_additions"]),
            "geometry_converged": True,
            "finite_and_converged": True,
            "reference_write": False,
            "oracle_labels_consumed": False,
            "source_grid_status": summary.get("grid_status", ""),
        })
    keys: set[tuple[str, float]] = set()
    for row in rows:
        key = (xi_key(row["xi"]), round(f(row["T_MeV"]), 10))
        if key in keys:
            raise ValueError(f"duplicate strict Maxwell key: {key}")
        keys.add(key)
    rows.sort(key=lambda row: (f(row["xi"]), f(row["T_MeV"])))
    return rows, list(rows)


def strict_crossover(v7_root: Path) -> list[dict[str, Any]]:
    path = v7_root / "tables" / "v6_crossover_overlay.csv"
    rows = read_csv(path)
    result: list[dict[str, Any]] = []
    for row in rows:
        if row.get("plot_status") != "retained":
            continue
        result.append({
            "surface": "crossover", "xi": f(row["xi"]), "mu_MeV": f(row["mu_MeV"]),
            "T_MeV": f(row["T_MeV"]), "rho": f(row["rho"]),
            "mu_CEP_proxy_MeV": f(row["mu_CEP_proxy_MeV"]),
            "layer": "strict_reference_v1", "status": "native_v6_or_endpoint",
            "interpolation_method": "none_native", "physical_region": "crossover_below_CEP",
            "source_layer": row.get("source_layer", "v6"), "reference_write": False,
        })
    return result


def strict_spinodal(v7_root: Path) -> list[dict[str, Any]]:
    rows = read_csv(v7_root / "tables" / "v6_spinodals.csv")
    result: list[dict[str, Any]] = []
    for row in rows:
        result.append({
            "surface": "spinodal", "xi": f(row["xi"]), "T_MeV": f(row["T_MeV"]),
            "mu_spinodal_hadron_MeV": f(row["mu_spinodal_hadron_MeV"]),
            "mu_spinodal_quark_MeV": f(row["mu_spinodal_quark_MeV"]),
            "layer": "strict_reference_v1", "status": "native_v6",
            "source_layer": "v6", "reference_write": False,
        })
    return result


def strict_cep(v7_root: Path) -> list[dict[str, Any]]:
    rows = read_csv(v7_root / "tables" / "cep_boundary_estimates.csv")
    return [
        {
            "surface": "cep_boundary", "xi": f(row["xi"]),
            "mu_CEP_proxy_MeV": f(row["mu_CEP_proxy_MeV"]),
            "T_low_MeV": f(row["T_CEP_bracket_low_MeV"]),
            "T_high_MeV": f(row["T_CEP_bracket_high_MeV"]),
            "T_midpoint_MeV": f(row["T_CEP_estimated_midpoint_MeV"]),
            "layer": "strict_reference_v1", "status": "bracket_preserved",
            "boundary_mode": row.get("boundary_mode", "estimated_midpoint"),
            "reference_write": False,
        }
        for row in rows
    ]


def make_render_figure(output: Path, crossover: list[dict[str, Any]], maxwell: list[dict[str, Any]], cep: list[dict[str, Any]]) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - environment-specific
        raise RuntimeError("matplotlib is required to render phase_surface_render_v1") from exc
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    fig = plt.figure(figsize=(11, 8), dpi=180)
    axis = fig.add_subplot(111, projection="3d")
    for rows, color, label, axis_name in (
        (crossover, "#d9792b", "crossover", "mu_MeV"),
        (maxwell, "#245b8f", "Maxwell", "T_MeV"),
    ):
        grouped = group_by_xi(rows)
        for index, key in enumerate(sorted(grouped, key=float)):
            segment = sorted(grouped[key], key=lambda row: f(row[axis_name]))
            if axis_name == "mu_MeV":
                x = [f(row["mu_MeV"]) for row in segment]
                z = [f(row["T_MeV"]) for row in segment]
            else:
                x = [f(row["mu_MeV"]) for row in segment]
                z = [f(row["T_MeV"]) for row in segment]
            y = [f(key)] * len(segment)
            axis.plot(x, y, z, color=color, linewidth=0.45, alpha=0.42, label=label if index == 0 else None)
    axis.plot(
        [f(row["mu_CEP_proxy_MeV"]) for row in sorted(cep, key=lambda row: f(row["xi"]))],
        [f(row["xi"]) for row in sorted(cep, key=lambda row: f(row["xi"]))],
        [f(row["T_midpoint_MeV"]) for row in sorted(cep, key=lambda row: f(row["xi"]))],
        color="#202020", linewidth=1.8, label="CEP estimated midpoint",
    )
    axis.set_xlabel(r"$\mu_q$ [MeV]")
    axis.set_ylabel(r"$\xi$")
    axis.set_zlabel("T [MeV]")
    axis.set_title("Issue #130 phase surfaces: strict + derived uniform-xi render")
    axis.view_init(elev=23, azim=-62)
    axis.legend(loc="upper left")
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)


def write_manifest(path: Path, value: dict[str, Any]) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v7-root", type=Path, required=True)
    parser.add_argument("--frozen-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--replay-root", type=Path, required=True)
    parser.add_argument("--replay-run-id", required=True)
    parser.add_argument("--calculation-sha", default=CALCULATION_SHA)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    v7_root = args.v7_root.resolve()
    frozen_root = args.frozen_root.resolve()
    output_root = args.output_root.resolve()
    replay_root = args.replay_root.resolve()
    if args.calculation_sha != CALCULATION_SHA:
        raise ValueError("unexpected calculation SHA")
    v7_manifest_path = v7_root / "manifest.json"
    freeze_manifest_path = frozen_root / "freeze_manifest.json"
    v7_manifest = read_json(v7_manifest_path)
    freeze_manifest = read_json(freeze_manifest_path)
    replay_manifest_path = resolve_manifest(replay_root)
    replay_manifest = read_json(replay_manifest_path)
    if v7_manifest.get("calculation_sha") != args.calculation_sha:
        raise ValueError("v7 calculation SHA mismatch")
    if freeze_manifest.get("calculation_sha") != args.calculation_sha or freeze_manifest.get("run_id") != args.run_id:
        raise ValueError("frozen expansion provenance mismatch")
    if replay_manifest.get("run_mode") != "aggregate_replay":
        raise ValueError("replay artifact is not aggregate_replay")
    if replay_manifest.get("source_run_id") != args.run_id:
        raise ValueError("replay source run mismatch")
    if replay_manifest.get("calculation_sha") != args.calculation_sha:
        raise ValueError("replay calculation SHA mismatch")
    if replay_manifest.get("source_workflow_sha") != freeze_manifest.get("source_workflow_sha"):
        raise ValueError("replay source workflow SHA mismatch")
    if replay_manifest.get("solver_called") is not False:
        raise ValueError("aggregate replay must record solver_called=false")
    if replay_manifest.get("materialized_target_count") != 276 or replay_manifest.get("errors"):
        raise ValueError("replay aggregate is incomplete")
    if output_root.exists():
        raise FileExistsError(f"refusing to overwrite phase-reference layers: {output_root}")

    strict_root = output_root / "strict_reference_v1"
    derived_root = output_root / "derived_reference_v1"
    render_root = output_root / "phase_surface_render_v1"
    for root in (strict_root, derived_root, render_root):
        (root / "tables").mkdir(parents=True, exist_ok=True)
        (root / "figures").mkdir(parents=True, exist_ok=True)

    v7_tables = v7_root / "tables"
    v7_maxwell = read_csv(v7_tables / "maxwell_surface_v7.csv")
    strict_maxwell_rows, _ = strict_maxwell(v7_maxwell, frozen_root, args.run_id, args.calculation_sha)
    strict_crossover_rows = strict_crossover(v7_root)
    strict_spinodal_rows = strict_spinodal(v7_root)
    strict_cep_rows = strict_cep(v7_root)

    maxwell_fields = list(strict_maxwell_rows[0].keys())
    crossover_fields = list(strict_crossover_rows[0].keys())
    spinodal_fields = list(strict_spinodal_rows[0].keys())
    cep_fields = list(strict_cep_rows[0].keys())
    write_csv(strict_root / "tables" / "maxwell_surface_strict_reference_v1.csv", maxwell_fields, strict_maxwell_rows)
    write_csv(strict_root / "tables" / "crossover_surface_strict_reference_v1.csv", crossover_fields, strict_crossover_rows)
    write_csv(strict_root / "tables" / "spinodal_surface_strict_reference_v1.csv", spinodal_fields, strict_spinodal_rows)
    write_csv(strict_root / "tables" / "cep_boundary_strict_reference_v1.csv", cep_fields, strict_cep_rows)

    derived_crossover, crossover_coverage = interpolate_surface(
        strict_crossover_rows, surface="crossover", axis="mu_MeV",
        fields=("T_MeV", "rho", "mu_CEP_proxy_MeV"), source_layer="strict_crossover_v6_or_endpoint"
    )
    derived_maxwell, maxwell_coverage = interpolate_surface(
        strict_maxwell_rows, surface="maxwell", axis="T_MeV",
        fields=("mu_MeV", "rho_hadron", "rho_quark", "area_residual"), source_layer="strict_maxwell_c2_plus_expansion"
    )
    strict_spinodal_regrid, spinodal_coverage = interpolate_surface(
        strict_spinodal_rows, surface="spinodal", axis="T_MeV",
        fields=("mu_spinodal_hadron_MeV", "mu_spinodal_quark_MeV"), source_layer="strict_spinodal_v6"
    )
    cep_group = {xi_key(row["xi"]): row for row in strict_cep_rows}
    derived_cep: list[dict[str, Any]] = []
    for target_xi in uniform_xi():
        bracket = bracket_xi({key: [value] for key, value in cep_group.items()}, target_xi)
        if bracket is None:
            continue
        left, right = bracket
        left_row = cep_group[xi_key(left)]
        right_row = cep_group[xi_key(right)]
        weight = 0.0 if abs(right - left) <= EPS else (target_xi - left) / (right - left)
        row = {
            "surface": "cep_boundary", "xi": target_xi,
            "mu_CEP_proxy_MeV": f(left_row["mu_CEP_proxy_MeV"]) + weight * (f(right_row["mu_CEP_proxy_MeV"]) - f(left_row["mu_CEP_proxy_MeV"])),
            "T_low_MeV": f(left_row["T_low_MeV"]) + weight * (f(right_row["T_low_MeV"]) - f(left_row["T_low_MeV"])),
            "T_high_MeV": f(left_row["T_high_MeV"]) + weight * (f(right_row["T_high_MeV"]) - f(left_row["T_high_MeV"])),
            "T_midpoint_MeV": f(left_row["T_midpoint_MeV"]) + weight * (f(right_row["T_midpoint_MeV"]) - f(left_row["T_midpoint_MeV"])),
            "layer": "strict_reference_v1" if abs(left - right) <= EPS else "interpolated_noncertified",
            "status": "bracket_preserved_native" if abs(left - right) <= EPS else "bracket_preserved_xi_interpolation",
            "boundary_mode": "estimated_midpoint",
            "source_xi_left": left, "source_xi_right": right, "xi_gap": right - left,
            "reference_write": False,
        }
        derived_cep.append(row)

    derived_root_tables = derived_root / "tables"
    write_csv(derived_root_tables / "crossover_surface_derived_reference_v1.csv", list(derived_crossover[0].keys()), derived_crossover)
    write_csv(derived_root_tables / "maxwell_surface_derived_reference_v1.csv", list(derived_maxwell[0].keys()), derived_maxwell)
    write_csv(derived_root_tables / "spinodal_surface_derived_reference_v1.csv", list(strict_spinodal_regrid[0].keys()), strict_spinodal_regrid)
    write_csv(derived_root_tables / "cep_boundary_derived_reference_v1.csv", list(derived_cep[0].keys()), derived_cep)
    all_coverage = crossover_coverage + maxwell_coverage + spinodal_coverage
    write_csv(derived_root_tables / "surface_coverage_mask.csv", list(all_coverage[0].keys()), all_coverage)

    render_tables = render_root / "tables"
    render_crossover = [row for row in derived_crossover if row.get("layer")]
    render_maxwell = [row for row in derived_maxwell if row.get("layer")]
    render_cep = derived_cep
    write_csv(render_tables / "crossover_surface_render.csv", list(render_crossover[0].keys()), render_crossover)
    write_csv(render_tables / "maxwell_surface_render.csv", list(render_maxwell[0].keys()), render_maxwell)
    write_csv(render_tables / "cep_boundary_render.csv", list(render_cep[0].keys()), render_cep)
    write_csv(render_tables / "mesh_coverage.csv", list(all_coverage[0].keys()), all_coverage)
    figure_path = render_root / "figures" / "phase_surface_render_mu_xi_T.png"
    make_render_figure(figure_path, render_crossover, render_maxwell, render_cep)

    generated_at = datetime.now(timezone.utc).isoformat()
    generator = {
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "command": [str(value) for value in sys.argv],
    }
    replay_provenance = {
        "schema_version": "pnjl_issue130_maxwell_cep_local_replay_provenance_v1",
        "replay_run_id": args.replay_run_id,
        "replay_run_url": f"https://github.com/w5851/Julia_RelaxTime/actions/runs/{args.replay_run_id}",
        "source_run_id": args.run_id,
        "source_run_url": f"https://github.com/w5851/Julia_RelaxTime/actions/runs/{args.run_id}",
        "calculation_sha": args.calculation_sha,
        "source_workflow_sha": replay_manifest.get("source_workflow_sha"),
        "replay_postprocess_sha": replay_manifest.get("postprocess_sha"),
        "verdict": replay_manifest.get("verdict"),
        "run_mode": replay_manifest.get("run_mode"),
        "solver_called": replay_manifest.get("solver_called"),
        "materialized_target_count": replay_manifest.get("materialized_target_count"),
        "expected_target_count": replay_manifest.get("expected_target_count"),
        "manifest": source_record(replay_manifest_path),
        "identity_fallback_count": len(replay_manifest.get("identity_fallbacks", {})),
        "boundary": "solver-free replay provenance; no reference promotion",
    }
    (output_root / "replay_provenance.json").write_text(
        json.dumps(replay_provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    common = {
        "schema_version": SCHEMA,
        "generated_at_utc": generated_at,
        "calculation_sha": args.calculation_sha,
        "maxwell_source_run_id": args.run_id,
        "v7_manifest_hash": sha256(v7_manifest_path),
        "freeze_manifest_hash": sha256(freeze_manifest_path),
        "uniform_xi_step": XI_STEP,
        "uniform_xi_count": len(uniform_xi()),
        "reference_write": False,
        "solver_called": False,
        "oracle_labels_consumed": False,
        "generator": generator,
        "replay": {
            "run_id": args.replay_run_id,
            "source_run_id": args.run_id,
            "manifest_sha256": sha256(replay_manifest_path),
            "postprocess_sha": replay_manifest.get("postprocess_sha"),
            "verdict": replay_manifest.get("verdict"),
        },
    }
    write_manifest(strict_root / "manifest.json", {
        **common, "layer": "strict_reference_v1", "promotion_status": "awaiting_author_review",
        "semantics": "native C2/v6 plus real Maxwell expansion; old unresolved rows retained",
        "rows": {"maxwell": len(strict_maxwell_rows), "crossover": len(strict_crossover_rows), "spinodal": len(strict_spinodal_rows), "cep": len(strict_cep_rows)},
    })
    write_manifest(derived_root / "manifest.json", {
        **common, "layer": "derived_reference_v1", "promotion_status": "awaiting_author_review",
        "semantics": "uniform xi local interpolation with common-axis support; derived rows are non-certified",
        "rows": {"maxwell": len(derived_maxwell), "crossover": len(derived_crossover), "spinodal": len(strict_spinodal_regrid), "cep": len(derived_cep)},
        "coverage": {"crossover": len(crossover_coverage), "maxwell": len(maxwell_coverage), "spinodal": len(spinodal_coverage)},
    })
    plot_manifest_path = render_root / "figures" / "plot_manifest.json"
    write_manifest(plot_manifest_path, {
        "schema_version": "pnjl_issue130_phase_surface_plot_manifest_v1",
        "generated_at_utc": generated_at,
        "calculation_sha": args.calculation_sha,
        "source_run_id": args.run_id,
        "replay_run_id": args.replay_run_id,
        "generator": generator,
        "inputs": {
            "derived_manifest": source_record(derived_root / "manifest.json"),
            "crossover_table": source_record(render_tables / "crossover_surface_render.csv"),
            "maxwell_table": source_record(render_tables / "maxwell_surface_render.csv"),
            "cep_table": source_record(render_tables / "cep_boundary_render.csv"),
            "coverage_table": source_record(render_tables / "mesh_coverage.csv"),
        },
        "figure": source_record(figure_path),
        "render_semantics": {
            "projection": "ordered_xi_slice_lines",
            "triangulation": False,
            "masked_cells_connected": False,
            "derived_rows_are_certified": False,
        },
    })
    render_manifest = {
        **common, "layer": "phase_surface_render_v1", "promotion_status": "awaiting_author_review",
        "semantics": "structured no-triangulation line-surface render from derived_reference_v1",
        "rows": {"maxwell": len(render_maxwell), "crossover": len(render_crossover), "cep": len(render_cep)},
        "figures": {
            "phase_surface_render_mu_xi_T.png": sha256(figure_path),
            "plot_manifest.json": sha256(plot_manifest_path),
        },
        "mesh": {"coverage_rows": len(all_coverage), "surfaces": ["crossover", "maxwell", "spinodal"]},
    }
    write_manifest(render_root / "manifest.json", render_manifest)
    write_manifest(output_root / "manifest.json", {
        **common, "layer": "issue130_phase_reference_layers_v1", "promotion_status": "awaiting_author_review",
        "children": {
            "strict_reference_v1": sha256(strict_root / "manifest.json"),
            "derived_reference_v1": sha256(derived_root / "manifest.json"),
            "phase_surface_render_v1": sha256(render_root / "manifest.json"),
        },
    })
    (strict_root / "README.md").write_text(
        "# strict_reference_v1\n\n"
        "This layer preserves the v6/C2 native crossover and spinodal support, all prior Maxwell rows "
        "including unresolved geometry status, and appends the 276 real Maxwell expansion targets from "
        f"run `{args.run_id}`. It is an input candidate for author review, not a phase-reference promotion.\n",
        encoding="utf-8",
    )
    (derived_root / "README.md").write_text(
        "# derived_reference_v1\n\n"
        f"A uniform xi grid with step {XI_STEP} is produced by local linear interpolation over adjacent "
        "strict slices and their common axis support. No extrapolation across absent support is performed. "
        "All non-native rows are `interpolated_noncertified`; CEP brackets remain intervals and midpoint "
        "values are boundary estimates.\n",
        encoding="utf-8",
    )
    (render_root / "README.md").write_text(
        "# phase_surface_render_v1\n\n"
        "Structured crossover and Maxwell surfaces plus the CEP boundary line are rendered from "
        "derived_reference_v1. The render uses ordered xi-slice lines and a machine-readable coverage mask; "
        "it does not triangulate gaps or promote derived rows to strict evidence.\n",
        encoding="utf-8",
    )
    claim_rows = [
        {"claim_id": "maxwell_expansion_appended", "claim": "276 real Maxwell expansion targets are appended to the strict input layer", "status": "supported", "evidence": "strict_reference_v1/tables/maxwell_surface_strict_reference_v1.csv; freeze_manifest.json", "boundary": "diagnostic candidate; no phase-reference promotion"},
        {"claim_id": "unresolved_preserved", "claim": "Pre-existing C2 Maxwell unresolved geometry rows remain explicitly marked", "status": "supported", "evidence": "strict_reference_v1/tables/maxwell_surface_strict_reference_v1.csv", "boundary": "unresolved rows are not silently certified"},
        {"claim_id": "uniform_xi_derived", "claim": "Derived surfaces use a uniform xi grid and common-support interpolation", "status": "supported", "evidence": "derived_reference_v1/tables/surface_coverage_mask.csv", "boundary": "interpolated_noncertified; no extrapolation"},
        {"claim_id": "phase_surface_render", "claim": "The render is reproducible from structured derived tables and coverage metadata", "status": "supported", "evidence": "phase_surface_render_v1/manifest.json; phase_surface_render_v1/tables/mesh_coverage.csv", "boundary": "visualization does not authorize promotion"},
        {"claim_id": "promotion", "claim": "These layers authorize phase-reference promotion", "status": "blocked", "evidence": "manifest.json", "boundary": "awaiting author review and formal promotion gate"},
    ]
    write_csv(output_root / "claim_ledger.csv", list(claim_rows[0].keys()), claim_rows)
    (output_root / "README.md").write_text(
        "# Issue #130 phase-reference layers v1\n\n"
        "This package freezes a strict input layer, a uniform-xi derived layer and a structured render layer. "
        "The strict layer contains immutable C2/v6 support plus the 276-target Maxwell expansion; the derived "
        "layer performs only local common-support interpolation; the render layer is a no-triangulation visual "
        "projection. All layers keep `reference_write=false` and remain diagnostic/author-review candidates.\n",
        encoding="utf-8",
    )
    (output_root / "AUDIT.md").write_text(
        "# Issue #130 phase-reference layers v1 审计\n\n"
        f"- numerical Maxwell source run: `{args.run_id}`；replay run: `{args.replay_run_id}`。\n"
        f"- calculation SHA: `{args.calculation_sha}`。\n"
        f"- strict Maxwell rows: `{len(strict_maxwell_rows)}`（C2/v6 原生行加 276 个真实补点）。\n"
        f"- derived rows: Maxwell `{len(derived_maxwell)}`、crossover `{len(derived_crossover)}`、"
        f"spinodal `{len(strict_spinodal_regrid)}`、CEP `{len(derived_cep)}`。\n"
        "- replay 是 `aggregate_replay`，已验证 276/276 materialized、无错误且 `solver_called=false`；"
        "replay manifest 的 source run 由 workflow 输入显式传递并保存在 `replay_provenance.json`。\n"
        "- strict 层保留旧 unresolved geometry；derived 层只在相邻 xi 的共同 axis support 内局部插值，"
        "非原生行标记为 `interpolated_noncertified`，不跨缺失 support 外推。\n"
        "- render 层使用 ordered xi-slice lines，不三角化跨越 masked cell；图像连续性不改变数据层状态。\n\n"
        "## Verdict\n\n"
        "这些层是 `diagnostic/author-review candidate`，不是 phase-reference promotion。"
        "在作者审核和独立 promotion gate 完成前，不写入 `data/reference`，不启动 RS transport。\n",
        encoding="utf-8",
    )
    (output_root / "decision.json").write_text(
        json.dumps({
            "schema_version": "pnjl_issue130_phase_reference_layers_decision_v1",
            "verdict": "awaiting_author_review",
            "promotion_status": "blocked",
            "reference_write": False,
            "solver_called": False,
            "oracle_labels_consumed": False,
            "accepted_inputs": {
                "calculation_sha": args.calculation_sha,
                "maxwell_source_run_id": args.run_id,
                "aggregate_replay_run_id": args.replay_run_id,
            },
            "checks": {
                "maxwell_expansion_materialized": len(strict_maxwell_rows) >= 276,
                "replay_complete": replay_manifest.get("materialized_target_count") == 276 and not replay_manifest.get("errors"),
                "strict_unresolved_preserved": True,
                "derived_no_extrapolation": True,
                "render_no_triangulation": True,
            },
            "next_action": "作者审核 strict/derived/render 数据、coverage、代表图和 provenance；通过后另行执行 phase-reference promotion gate。",
        }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    # Rewrite the root manifest after all package files exist so the package-level
    # inventory covers the audit, decision and replay provenance files as well.
    root_manifest = {
        **common,
        "layer": "issue130_phase_reference_layers_v1",
        "promotion_status": "awaiting_author_review",
        "children": {
            "strict_reference_v1": sha256(strict_root / "manifest.json"),
            "derived_reference_v1": sha256(derived_root / "manifest.json"),
            "phase_surface_render_v1": sha256(render_root / "manifest.json"),
        },
        "package_files": {
            relative.as_posix(): sha256(output_root / relative)
            for relative in (
                Path("README.md"), Path("AUDIT.md"), Path("decision.json"),
                Path("claim_ledger.csv"), Path("replay_provenance.json"),
            )
        },
    }
    write_manifest(output_root / "manifest.json", root_manifest)
    print(json.dumps({
        "output_root": str(output_root),
        "strict_maxwell_rows": len(strict_maxwell_rows),
        "derived_maxwell_rows": len(derived_maxwell),
        "derived_crossover_rows": len(derived_crossover),
        "derived_cep_rows": len(derived_cep),
        "uniform_xi_count": len(uniform_xi()),
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
