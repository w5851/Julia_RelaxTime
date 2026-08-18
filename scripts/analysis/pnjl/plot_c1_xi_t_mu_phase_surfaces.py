#!/usr/bin/env python3
"""Plot diagnostic PNJL (xi, T, mu) phase surfaces from a C1 artifact.

The script is deliberately solver-free.  It consumes already materialized
boundary and crossover CSV files, validates their keys and finite values, and
writes a traceable diagnostic package under ``docs/analysis``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.lines import Line2D


BOUNDARY_COLUMNS = ("xi", "T_MeV", "mu_transition_MeV", "rho_hadron", "rho_quark", "converged")
SPINODAL_COLUMNS = (
    "xi",
    "T_MeV",
    "mu_spinodal_hadron_MeV",
    "mu_spinodal_quark_MeV",
    "rho_spinodal_hadron",
    "rho_spinodal_quark",
)
CROSSOVER_COLUMNS = ("xi", "mu_MeV", "T_crossover_MeV", "rho", "method", "converged")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--boundary", required=True, type=Path)
    parser.add_argument("--spinodals", required=True, type=Path)
    parser.add_argument("--crossover", required=True, type=Path)
    parser.add_argument("--grid-convergence", required=True, type=Path)
    parser.add_argument("--source-manifest", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument(
        "--orientation",
        choices=("xi_T_mu", "mu_xi_T"),
        default="xi_T_mu",
        help="coordinate order for the rendered axes (default: xi,T,mu)",
    )
    return parser.parse_args()


def fail(message: str) -> None:
    raise SystemExit(f"[c1-phase-surfaces] {message}")


def read_csv(path: Path, required: Iterable[str]) -> list[dict[str, str]]:
    if not path.is_file():
        fail(f"missing input: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        missing = sorted(set(required) - fields)
        if missing:
            fail(f"{path.name} missing columns: {missing}")
        return list(reader)


def finite(row: dict[str, str], field: str, path: Path, row_number: int) -> float:
    try:
        value = float(row[field])
    except (KeyError, TypeError, ValueError) as exc:
        fail(f"{path.name}:{row_number} non-numeric {field}: {row.get(field)!r}")
        raise AssertionError from exc
    if not math.isfinite(value):
        fail(f"{path.name}:{row_number} non-finite {field}: {value}")
    return value


def bool_value(row: dict[str, str], field: str, path: Path, row_number: int) -> bool:
    value = row.get(field, "").strip().lower()
    if value not in {"true", "false"}:
        fail(f"{path.name}:{row_number} invalid boolean {field}: {value!r}")
    return value == "true"


def rounded_key(values: Iterable[float]) -> tuple[float, ...]:
    return tuple(round(value, 10) for value in values)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def repo_head() -> str | None:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip()


def validate_inputs(args: argparse.Namespace) -> tuple[list[dict], list[dict], list[dict], list[dict], dict]:
    boundary = read_csv(args.boundary, BOUNDARY_COLUMNS)
    spinodals = read_csv(args.spinodals, SPINODAL_COLUMNS)
    crossover = read_csv(args.crossover, CROSSOVER_COLUMNS)
    grid = read_csv(args.grid_convergence, ("axis", "xi", "T_MeV", "converged", "reason"))

    boundary_keys: set[tuple[float, ...]] = set()
    for number, row in enumerate(boundary, 2):
        key = rounded_key((finite(row, "xi", args.boundary, number), finite(row, "T_MeV", args.boundary, number)))
        if key in boundary_keys:
            fail(f"duplicate boundary key at row {number}: {key}")
        boundary_keys.add(key)
        for field in BOUNDARY_COLUMNS[2:5]:
            finite(row, field, args.boundary, number)
        if not bool_value(row, "converged", args.boundary, number):
            fail(f"boundary contains non-converged row {number}; this plot only uses converged rows")

    spinodal_keys: set[tuple[float, ...]] = set()
    for number, row in enumerate(spinodals, 2):
        key = rounded_key((finite(row, "xi", args.spinodals, number), finite(row, "T_MeV", args.spinodals, number)))
        if key in spinodal_keys:
            fail(f"duplicate spinodal key at row {number}: {key}")
        spinodal_keys.add(key)
        for field in SPINODAL_COLUMNS[2:]:
            finite(row, field, args.spinodals, number)

    crossover_keys: set[tuple[float, ...]] = set()
    for number, row in enumerate(crossover, 2):
        key = rounded_key((finite(row, "xi", args.crossover, number), finite(row, "mu_MeV", args.crossover, number)))
        if key in crossover_keys:
            fail(f"duplicate crossover key at row {number}: {key}")
        crossover_keys.add(key)
        for field in ("T_crossover_MeV", "rho"):
            finite(row, field, args.crossover, number)
        if not bool_value(row, "converged", args.crossover, number):
            fail(f"crossover contains non-converged row {number}; this plot only uses converged rows")

    for number, row in enumerate(grid, 2):
        finite(row, "xi", args.grid_convergence, number)
        # xi-axis interpolation records intentionally have no physical T
        # coordinate; rho and T records must carry both coordinates.
        if row["axis"].strip().lower() != "xi":
            finite(row, "T_MeV", args.grid_convergence, number)
        bool_value(row, "converged", args.grid_convergence, number)

    try:
        source_manifest = json.loads(args.source_manifest.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        fail(f"cannot read source manifest {args.source_manifest}: {exc}")

    expected = {
        "boundary": len(boundary),
        "spinodals": len(spinodals),
        "crossover": len(crossover),
        "grid_convergence": len(grid),
    }
    manifest_artifacts = source_manifest.get("artifacts", {})
    for name, count in expected.items():
        declared = manifest_artifacts.get(name, {}).get("row_count")
        if declared is not None and int(declared) != count:
            fail(f"manifest row-count mismatch for {name}: declared={declared}, actual={count}")

    return boundary, spinodals, crossover, grid, source_manifest


def set_triangle_mask(triangulation: mtri.Triangulation, x: np.ndarray, y: np.ndarray, max_dx: float, max_dy: float) -> None:
    triangles = triangulation.triangles
    dx = np.ptp(x[triangles], axis=1)
    dy = np.ptp(y[triangles], axis=1)
    triangulation.set_mask((dx > max_dx) | (dy > max_dy))


def surface_triangles(x: np.ndarray, y: np.ndarray, max_dx: float, max_dy: float) -> np.ndarray:
    triangulation = mtri.Triangulation(x, y)
    set_triangle_mask(triangulation, x, y, max_dx, max_dy)
    mask = triangulation.mask
    return triangulation.triangles if mask is None else triangulation.triangles[~mask]


def audit_crossover_vs_maxwell(
    boundary: list[dict],
    spinodals: list[dict],
    crossover: list[dict],
) -> tuple[list[dict], dict[str, int], list[dict]]:
    """Classify crossover samples relative to the Maxwell/spinodal curves.

    This is a geometric audit only.  It does not reinterpret the crossover
    detector or promote a response peak to a phase boundary.
    """
    boundary_by_xi: dict[float, list[tuple[float, float]]] = defaultdict(list)
    spinodal_by_xi: dict[float, list[tuple[float, float, float]]] = defaultdict(list)
    crossover_by_xi: dict[float, list[dict]] = defaultdict(list)
    for row in boundary:
        xi = round(float(row["xi"]), 10)
        boundary_by_xi[xi].append((float(row["T_MeV"]), float(row["mu_transition_MeV"])))
    for row in spinodals:
        xi = round(float(row["xi"]), 10)
        mu_a = float(row["mu_spinodal_hadron_MeV"])
        mu_b = float(row["mu_spinodal_quark_MeV"])
        spinodal_by_xi[xi].append((float(row["T_MeV"]), min(mu_a, mu_b), max(mu_a, mu_b)))
    for row in crossover:
        crossover_by_xi[round(float(row["xi"]), 10)].append(row)

    audit_rows: list[dict] = []
    relation_counts: dict[str, int] = defaultdict(int)
    endpoint_rows: list[dict] = []
    for xi, rows in sorted(crossover_by_xi.items()):
        boundary_rows = sorted(boundary_by_xi[xi])
        spinodal_rows = sorted(spinodal_by_xi[xi])
        boundary_T = np.asarray([row[0] for row in boundary_rows])
        boundary_mu = np.asarray([row[1] for row in boundary_rows])
        spinodal_T = np.asarray([row[0] for row in spinodal_rows])
        spinodal_low = np.asarray([row[1] for row in spinodal_rows])
        spinodal_high = np.asarray([row[2] for row in spinodal_rows])
        T_maxwell_max = float(boundary_T.max())
        T_crossover_min = min(float(row["T_crossover_MeV"]) for row in rows)
        endpoint_rows.append(
            {
                "xi": f"{xi:.12g}",
                "T_maxwell_max_MeV": f"{T_maxwell_max:.12g}",
                "T_crossover_min_MeV": f"{T_crossover_min:.12g}",
                "delta_T_crossover_minus_maxwell_MeV": f"{T_crossover_min - T_maxwell_max:.12g}",
                "overlap_in_T": str(T_crossover_min <= T_maxwell_max).lower(),
            }
        )
        for row in rows:
            mu = float(row["mu_MeV"])
            T = float(row["T_crossover_MeV"])
            if T < boundary_T.min() or T > boundary_T.max():
                relation = "above_maxwell_T_range"
                mu_maxwell = ""
                mu_low = ""
                mu_high = ""
                maxwell_side = "outside_boundary_T_range"
            else:
                mu_maxwell_value = float(np.interp(T, boundary_T, boundary_mu))
                mu_low_value = float(np.interp(T, spinodal_T, spinodal_low))
                mu_high_value = float(np.interp(T, spinodal_T, spinodal_high))
                mu_maxwell = f"{mu_maxwell_value:.12g}"
                mu_low = f"{mu_low_value:.12g}"
                mu_high = f"{mu_high_value:.12g}"
                if mu < mu_low_value:
                    relation = "stable_hadron_side"
                elif mu > mu_high_value:
                    relation = "stable_quark_side"
                else:
                    relation = "inside_spinodal_band"
                if abs(mu - mu_maxwell_value) <= 0.5:
                    maxwell_side = "near_maxwell"
                elif mu < mu_maxwell_value:
                    maxwell_side = "below_maxwell_mu"
                else:
                    maxwell_side = "above_maxwell_mu"
            relation_counts[relation] += 1
            audit_rows.append(
                {
                    "xi": f"{xi:.12g}",
                    "mu_crossover_MeV": f"{mu:.12g}",
                    "T_crossover_MeV": f"{T:.12g}",
                    "mu_maxwell_interp_MeV": mu_maxwell,
                    "mu_spinodal_low_interp_MeV": mu_low,
                    "mu_spinodal_high_interp_MeV": mu_high,
                    "relation_to_spinodals": relation,
                    "relation_to_maxwell": maxwell_side,
                }
            )
    return audit_rows, dict(sorted(relation_counts.items())), endpoint_rows


def plot_surface(
    boundary: list[dict],
    crossover: list[dict],
    grid: list[dict],
    output_path: Path,
    orientation: str,
) -> dict:
    bx = np.asarray([float(row["xi"]) for row in boundary])
    by = np.asarray([float(row["T_MeV"]) for row in boundary])
    bz = np.asarray([float(row["mu_transition_MeV"]) for row in boundary])
    cx = np.asarray([float(row["xi"]) for row in crossover])
    cy = np.asarray([float(row["T_crossover_MeV"]) for row in crossover])
    cz = np.asarray([float(row["mu_MeV"]) for row in crossover])

    if orientation == "mu_xi_T":
        boundary_x, boundary_y, boundary_z = bz, bx, by
        crossover_x, crossover_y, crossover_z = cz, cx, cy
        x_label, y_label, z_label = r"$\mu$ (MeV)", r"$\xi$", r"$T$ (MeV)"
        title = "C1 diagnostic phase surfaces in $(\\mu, \\xi, T)$"
        box_aspect = (1.55, 1.0, 1.35)
    else:
        boundary_x, boundary_y, boundary_z = bx, by, bz
        crossover_x, crossover_y, crossover_z = cx, cy, cz
        x_label, y_label, z_label = r"$\xi$", r"$T$ (MeV)", r"$\mu$ (MeV)"
        title = "C1 diagnostic phase surfaces in $(\\xi, T, \\mu)$"
        box_aspect = (1.0, 1.55, 1.35)

    figure = plt.figure(figsize=(12, 8.5))
    axis = figure.add_subplot(111, projection="3d")

    # Maxwell is sampled on (xi, T), while crossover is sampled on (xi, mu).
    # The resulting vertices are always rendered in the requested (xi, T, mu)
    # coordinate order; masks prevent triangulation across large data gaps.
    boundary_triangles = surface_triangles(bx, by, max_dx=0.11, max_dy=12.0)
    crossover_parameter_triangles = surface_triangles(cx, cz, max_dx=0.11, max_dy=55.0)
    axis.plot_trisurf(
        boundary_x,
        boundary_y,
        boundary_z,
        triangles=boundary_triangles,
        color="#2f6db0",
        alpha=0.55,
        linewidth=0.15,
        edgecolor="#24527f",
        shade=True,
    )
    axis.plot_trisurf(
        crossover_x,
        crossover_y,
        crossover_z,
        triangles=crossover_parameter_triangles,
        color="#d9792b",
        alpha=0.38,
        linewidth=0.12,
        edgecolor="#9a4b11",
        shade=True,
    )

    # Thin generator curves make the sparse/adaptive sampling visible and keep
    # the surface interpretation auditable.
    for xi in sorted({round(float(row["xi"]), 10) for row in boundary}):
        rows = sorted((row for row in boundary if round(float(row["xi"]), 10) == xi), key=lambda row: float(row["T_MeV"]))
        if orientation == "mu_xi_T":
            axis.plot([float(row["mu_transition_MeV"]) for row in rows], [xi] * len(rows), [float(row["T_MeV"]) for row in rows], color="#1f4f86", alpha=0.22, linewidth=0.45)
        else:
            axis.plot([xi] * len(rows), [float(row["T_MeV"]) for row in rows], [float(row["mu_transition_MeV"]) for row in rows], color="#1f4f86", alpha=0.22, linewidth=0.45)
    for xi in sorted({round(float(row["xi"]), 10) for row in crossover}):
        rows = sorted((row for row in crossover if round(float(row["xi"]), 10) == xi), key=lambda row: float(row["mu_MeV"]))
        if orientation == "mu_xi_T":
            axis.plot([float(row["mu_MeV"]) for row in rows], [xi] * len(rows), [float(row["T_crossover_MeV"]) for row in rows], color="#9a4b11", alpha=0.18, linewidth=0.35)
        else:
            axis.plot([xi] * len(rows), [float(row["T_crossover_MeV"]) for row in rows], [float(row["mu_MeV"]) for row in rows], color="#9a4b11", alpha=0.18, linewidth=0.35)

    legend = [
        Line2D([0], [0], color="#2f6db0", lw=8, alpha=0.55, label="Maxwell coexistence surface"),
        Line2D([0], [0], color="#d9792b", lw=8, alpha=0.45, label="crossover pseudocritical surface"),
    ]
    axis.legend(handles=legend, loc="upper left", bbox_to_anchor=(0.02, 0.98), framealpha=0.9)
    axis.set_xlabel(x_label, labelpad=10)
    axis.set_ylabel(y_label, labelpad=10)
    axis.set_zlabel(z_label, labelpad=10)
    axis.set_title(title)
    axis.view_init(elev=25, azim=-62)
    try:
        axis.set_box_aspect(box_aspect)
    except AttributeError:
        pass
    figure.text(
        0.015,
        0.015,
        f"C1 diagnostic-only | Maxwell points={len(boundary)} | crossover points={len(crossover)} | grid convergence records={len(grid)}; unresolved={sum(row['converged'].strip().lower() == 'false' for row in grid)}",
        fontsize=8,
        color="#444444",
    )
    figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(figure)
    return {
        "boundary_triangle_count": int(len(boundary_triangles)),
        "crossover_triangle_count": int(len(crossover_parameter_triangles)),
        "path": str(output_path),
        "orientation": orientation,
    }


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    boundary, spinodals, crossover, grid, source_manifest = validate_inputs(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = args.out_dir / "figures"
    tables_dir = args.out_dir / "tables"
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    unresolved_by_axis: dict[str, int] = defaultdict(int)
    for row in grid:
        if row["converged"].strip().lower() == "false":
            unresolved_by_axis[row["axis"]] += 1
    summary_rows = [
        {
            "surface": "maxwell_coexistence",
            "input_rows": len(boundary),
            "xi_count": len({round(float(row["xi"]), 10) for row in boundary}),
            "T_min_MeV": min(float(row["T_MeV"]) for row in boundary),
            "T_max_MeV": max(float(row["T_MeV"]) for row in boundary),
            "mu_min_MeV": min(float(row["mu_transition_MeV"]) for row in boundary),
            "mu_max_MeV": max(float(row["mu_transition_MeV"]) for row in boundary),
            "converged_rows": sum(bool_value(row, "converged", args.boundary, number) for number, row in enumerate(boundary, 2)),
        },
        {
            "surface": "crossover_pseudocritical",
            "input_rows": len(crossover),
            "xi_count": len({round(float(row["xi"]), 10) for row in crossover}),
            "T_min_MeV": min(float(row["T_crossover_MeV"]) for row in crossover),
            "T_max_MeV": max(float(row["T_crossover_MeV"]) for row in crossover),
            "mu_min_MeV": min(float(row["mu_MeV"]) for row in crossover),
            "mu_max_MeV": max(float(row["mu_MeV"]) for row in crossover),
            "converged_rows": sum(bool_value(row, "converged", args.crossover, number) for number, row in enumerate(crossover, 2)),
        },
    ]
    boundary_by_xi: dict[float, list[dict]] = defaultdict(list)
    for row in boundary:
        boundary_by_xi[round(float(row["xi"]), 10)].append(row)
    maxwell_support_rows = []
    for xi, rows in sorted(boundary_by_xi.items()):
        top = max(rows, key=lambda row: float(row["T_MeV"]))
        bottom = min(rows, key=lambda row: float(row["T_MeV"]))
        maxwell_support_rows.append(
            {
                "xi": f"{xi:.12g}",
                "boundary_row_count": len(rows),
                "T_min_MeV": f"{float(bottom['T_MeV']):.12g}",
                "T_max_MeV": f"{float(top['T_MeV']):.12g}",
                "mu_at_T_max_MeV": f"{float(top['mu_transition_MeV']):.12g}",
                "mu_at_T_min_MeV": f"{float(bottom['mu_transition_MeV']):.12g}",
            }
        )
    write_csv(
        tables_dir / "surface_point_summary.csv",
        ["surface", "input_rows", "xi_count", "T_min_MeV", "T_max_MeV", "mu_min_MeV", "mu_max_MeV", "converged_rows"],
        summary_rows,
    )
    write_csv(
        tables_dir / "maxwell_surface_support_audit.csv",
        ["xi", "boundary_row_count", "T_min_MeV", "T_max_MeV", "mu_at_T_max_MeV", "mu_at_T_min_MeV"],
        maxwell_support_rows,
    )

    crossover_audit_rows, crossover_relation_counts, endpoint_rows = audit_crossover_vs_maxwell(boundary, spinodals, crossover)
    write_csv(
        tables_dir / "crossover_vs_maxwell_audit.csv",
        [
            "xi",
            "mu_crossover_MeV",
            "T_crossover_MeV",
            "mu_maxwell_interp_MeV",
            "mu_spinodal_low_interp_MeV",
            "mu_spinodal_high_interp_MeV",
            "relation_to_spinodals",
            "relation_to_maxwell",
        ],
        crossover_audit_rows,
    )
    write_csv(
        tables_dir / "crossover_maxwell_endpoint_separation.csv",
        ["xi", "T_maxwell_max_MeV", "T_crossover_min_MeV", "delta_T_crossover_minus_maxwell_MeV", "overlap_in_T"],
        endpoint_rows,
    )

    figure_name = "c1_mu_xi_T_phase_surfaces.png" if args.orientation == "mu_xi_T" else "c1_xi_T_mu_phase_surfaces.png"
    figure_path = figures_dir / figure_name
    figure_info = plot_surface(boundary, crossover, grid, figure_path, args.orientation)
    source_paths = {
        "boundary": args.boundary,
        "spinodals": args.spinodals,
        "crossover": args.crossover,
        "grid_convergence": args.grid_convergence,
        "source_manifest": args.source_manifest,
    }
    script_path = Path(__file__).resolve()
    figure_info.update({
        "sha256": sha256(figure_path),
        "generator": str(script_path),
        "generator_sha256": sha256(script_path),
    })
    claim_rows = [
        {
            "claim_id": "OBS-MAXWELL-SURFACE",
            "claim_type": "observation",
            "status": "candidate",
            "claim": "C1 converged boundary rows define a diagnostic Maxwell coexistence surface in (xi,T,mu).",
            "evidence_file": args.boundary.name,
            "evidence_field": "xi,T_MeV,mu_transition_MeV,converged",
        },
        {
            "claim_id": "OBS-CROSSOVER-SURFACE",
            "claim_type": "observation",
            "status": "candidate",
            "claim": "C1 converged peak-crossover rows define a diagnostic pseudocritical surface in (xi,T,mu).",
            "evidence_file": args.crossover.name,
            "evidence_field": "xi,mu_MeV,T_crossover_MeV,converged",
        },
        {
            "claim_id": "LIMIT-C1-CONVERGENCE",
            "claim_type": "boundary",
            "status": "author_check",
            "claim": "The surfaces are not a formal phase-reference claim because the C1 grid convergence artifact contains unresolved records and CEP is ambiguous.",
            "evidence_file": args.grid_convergence.name,
            "evidence_field": "axis,converged,reason",
        },
    ]
    write_csv(
        tables_dir / "claim_ledger.csv",
        ["claim_id", "claim_type", "status", "claim", "evidence_file", "evidence_field"],
        claim_rows,
    )
    manifest = {
        "schema_version": "pnjl_c1_phase_surfaces_v2",
        "status": "diagnostic_only",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "repo_head": repo_head(),
        "source_manifest": str(args.source_manifest),
        "source_provenance": source_manifest.get("provenance", {}),
        "input_files": {name: {"path": str(path), "sha256": sha256(path)} for name, path in source_paths.items()},
        "row_counts": {"boundary": len(boundary), "spinodals": len(spinodals), "crossover": len(crossover), "grid_convergence": len(grid)},
        "unresolved_grid_records_by_axis": dict(sorted(unresolved_by_axis.items())),
        "surface_summary": summary_rows,
        "maxwell_surface_support_audit": maxwell_support_rows,
        "crossover_relation_counts": crossover_relation_counts,
        "crossover_maxwell_endpoint_separation": endpoint_rows,
        "figure": figure_info,
        "scope": {
            "axes": ["mu_MeV", "xi", "T_MeV"] if args.orientation == "mu_xi_T" else ["xi", "T_MeV", "mu_MeV"],
            "maxwell_surface": "boundary mu_transition_MeV over xi,T",
            "crossover_surface": "crossover T_crossover_MeV over xi,mu",
            "no_solver_called": True,
            "cep_surface": "not plotted; C1 CEP rows are ambiguous",
            "unresolved_regions": "not interpolated into a physical claim",
        },
    }
    (args.out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (figures_dir / "plot_manifest.json").write_text(
        json.dumps(
            {
                "schema_version": "pnjl_plot_manifest_v2",
                "generator": str(script_path),
                "generator_sha256": sha256(script_path),
        "input_files": {name: {"path": str(path), "sha256": sha256(path)} for name, path in source_paths.items()},
                "figures": [figure_info],
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )
    readme = f"""# C1 `({args.orientation})` 诊断相面

本目录基于 C1 artifact `31762201725`（calculation SHA `{source_manifest.get('provenance', {}).get('calculation_git_commit', 'unknown')}`）生成，不调用 solver。

图中蓝色半透明曲面是 `boundary_*.csv` 的 Maxwell 共存化学势，橙色半透明曲面是 `crossover_*.csv` 的 peak crossover。坐标顺序为 `({"mu, xi, T" if args.orientation == "mu_xi_T" else "xi, T, mu"})`。

本图只是诊断可视化，不是 phase-reference 晋升结果。C1 的 `phase_grid_convergence` 共 {len(grid)} 条记录，其中 {sum(unresolved_by_axis.values())} 条未闭合；这些区域没有被强行插值成物理结论。C1 CEP 仍为 ambiguous，因此未绘制 CEP 曲线。

`crossover` 是对每个采样 `(xi,mu)` 独立在温区内寻找 `peak`，原始 CSV 没有用 Maxwell 面做物理遮罩。因此 `tables/crossover_vs_maxwell_audit.csv` 只做几何关系审计，不把一阶区域内的响应峰自动解释成第二张物理 crossover 面。`tables/crossover_maxwell_endpoint_separation.csv` 给出每个 ξ 的 Maxwell 最高温度与 crossover 最低温度差；它用于检查两张面是否在同一温度域重叠。

输入审计、哈希、三角剖分数量、派生点表、Maxwell 支持上边界、crossover 关系审计、端点温度分离和 claim ledger 见 `manifest.json`、`figures/plot_manifest.json`、`tables/surface_point_summary.csv`、`tables/maxwell_surface_support_audit.csv`、`tables/crossover_vs_maxwell_audit.csv`、`tables/crossover_maxwell_endpoint_separation.csv` 和 `tables/claim_ledger.csv`。
"""
    (args.out_dir / "README.md").write_text(readme, encoding="utf-8")
    print(json.dumps({"out_dir": str(args.out_dir), "figure": figure_info, "unresolved_grid_records_by_axis": dict(unresolved_by_axis)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
