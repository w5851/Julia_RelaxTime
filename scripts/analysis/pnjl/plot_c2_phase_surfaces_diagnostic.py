#!/usr/bin/env python3
"""Render a C2 diagnostic phase-surface view without calling the solver.

The figure is intended for global author review.  It renders finite C2
boundary/crossover support, projects each CEP temperature bracket at the last
first-order chemical potential, and overlays C2 grid-convergence status on the
Maxwell support.  No unresolved record is promoted or overwritten.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.lines import Line2D


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--c2-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument(
        "--mode",
        choices=("detailed", "visual_closed", "visual_lines"),
        default="detailed",
        help="detailed keeps audit overlays; visual_closed hides status overlays; visual_lines disables triangulation",
    )
    parser.add_argument("--max-dxi", type=float, default=0.07)
    parser.add_argument("--max-dT", type=float, default=8.0)
    parser.add_argument("--max-dmuq", type=float, default=55.0)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"CSV has no rows: {path}")
    return rows


def write_csv(path: Path, fields: Iterable[str], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields))
        writer.writeheader()
        writer.writerows(rows)


def finite(row: dict[str, str], field: str) -> float:
    value = float(row[field])
    if not math.isfinite(value):
        raise ValueError(f"non-finite {field}: {value}")
    return value


def is_true(value: str | bool) -> bool:
    return str(value).strip().lower() == "true"


def optional_finite(row: dict[str, str], field: str) -> float | None:
    """Return a finite optional numeric diagnostic, preserving blank cells."""
    raw = row.get(field, "")
    if raw is None or not str(raw).strip():
        return None
    try:
        value = float(raw)
    except ValueError as exc:
        raise ValueError(f"invalid numeric {field}: {raw!r}") from exc
    if not math.isfinite(value):
        raise ValueError(f"non-finite optional {field}: {value}")
    return value


def key(xi: float, T: float) -> tuple[float, float]:
    return round(xi, 10), round(T, 10)


def masked_triangles(x: np.ndarray, y: np.ndarray, max_dx: float, max_dy: float) -> np.ndarray:
    triangulation = mtri.Triangulation(x, y)
    triangles = triangulation.triangles
    dx = np.ptp(x[triangles], axis=1)
    dy = np.ptp(y[triangles], axis=1)
    return triangles[(dx <= max_dx) & (dy <= max_dy)]


def main() -> int:
    args = parse_args()
    c2_root = args.c2_root.resolve()
    reference = c2_root / "reference"
    output_root = args.output_root.resolve()
    figures = output_root / "figures"
    tables = output_root / "tables"
    figures.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)

    boundary_path = next(reference.glob("boundary_*.csv"))
    spinodal_path = next(reference.glob("spinodals_*.csv"))
    crossover_path = next(reference.glob("crossover_*.csv"))
    cep_path = next(reference.glob("cep_*.csv"))
    grid_path = next(reference.glob("phase_grid_convergence_*.csv"))
    source_manifest_path = next(reference.glob("phase_reference_*_manifest.json"))

    boundary = read_csv(boundary_path)
    spinodals = read_csv(spinodal_path)
    crossover = read_csv(crossover_path)
    cep = read_csv(cep_path)
    grid = read_csv(grid_path)
    source_manifest = json.loads(source_manifest_path.read_text(encoding="utf-8"))
    calculation_sha = str(source_manifest.get("calculation_git_commit") or source_manifest.get("git_commit") or "")
    if not calculation_sha:
        validation_path = c2_root / "validation_report.json"
        validation = json.loads(validation_path.read_text(encoding="utf-8")) if validation_path.is_file() else {}
        calculation_sha = str(validation.get("calculation_git_commit") or validation.get("git_commit") or "")

    required_boundary = {"xi", "T_MeV", "mu_transition_MeV", "rho_hadron", "rho_quark", "area_residual", "converged"}
    required_crossover = {"xi", "mu_MeV", "T_crossover_MeV", "rho", "converged"}
    required_cep = {"xi", "T_bracket_low_MeV", "T_bracket_high_MeV", "muq_last_first_order_MeV", "result_status"}
    if not required_boundary.issubset(boundary[0]):
        raise ValueError("boundary schema is incomplete")
    if not required_crossover.issubset(crossover[0]):
        raise ValueError("crossover schema is incomplete")
    if not required_cep.issubset(cep[0]):
        raise ValueError("CEP schema is incomplete")

    boundary_keys: set[tuple[float, float]] = set()
    boundary_plot: list[dict[str, object]] = []
    for row in boundary:
        xi = finite(row, "xi")
        T = finite(row, "T_MeV")
        item_key = key(xi, T)
        if item_key in boundary_keys:
            raise ValueError(f"duplicate boundary key: {item_key}")
        boundary_keys.add(item_key)
        if not is_true(row["converged"]):
            continue
        boundary_plot.append({
            "xi": xi,
            "T_MeV": T,
            # The phase-reference CSV stores the quark chemical potential
            # directly.  Keep the same mu_q convention for every surface.
            "mu_MeV": finite(row, "mu_transition_MeV"),
            "rho_hadron": finite(row, "rho_hadron"),
            "rho_quark": finite(row, "rho_quark"),
            "area_residual": finite(row, "area_residual"),
            "key": item_key,
        })

    cep_mu_by_xi: dict[float, float] = {}
    cep_mu_source_by_xi: dict[float, str] = {}
    cep_bracket_by_xi: dict[float, tuple[float, float]] = {}
    for row in cep:
        xi = finite(row, "xi")
        cep_key = round(xi, 10)
        if cep_key in cep_mu_by_xi:
            raise ValueError(f"duplicate CEP xi: {cep_key}")
        confirmed_mu = row.get("muq_CEP_MeV", "").strip()
        if confirmed_mu:
            try:
                confirmed_value = float(confirmed_mu)
            except ValueError as exc:
                raise ValueError(f"invalid muq_CEP_MeV for xi={xi}: {confirmed_mu!r}") from exc
            if not math.isfinite(confirmed_value):
                raise ValueError(f"non-finite muq_CEP_MeV for xi={xi}: {confirmed_value}")
            cep_mu_by_xi[cep_key] = confirmed_value
            cep_mu_source_by_xi[cep_key] = "muq_CEP_MeV"
        else:
            cep_mu_by_xi[cep_key] = finite(row, "muq_last_first_order_MeV")
            cep_mu_source_by_xi[cep_key] = "muq_last_first_order_MeV"
        low = finite(row, "T_bracket_low_MeV")
        high = finite(row, "T_bracket_high_MeV")
        cep_bracket_by_xi[cep_key] = (min(low, high), max(low, high))

    crossover_plot: list[dict[str, object]] = []
    crossover_excluded: list[dict[str, object]] = []
    crossover_filter_rows: list[dict[str, object]] = []
    crossover_support_by_xi: defaultdict[float, list[tuple[float, float]]] = defaultdict(list)
    crossover_keys: set[tuple[float, float]] = set()
    for row in crossover:
        xi = finite(row, "xi")
        muq = finite(row, "mu_MeV")
        item_key = (round(xi, 10), round(muq, 10))
        if item_key in crossover_keys:
            raise ValueError(f"duplicate crossover key: {item_key}")
        crossover_keys.add(item_key)
        if not is_true(row["converged"]):
            continue
        T = finite(row, "T_crossover_MeV")
        rho = finite(row, "rho")
        crossover_support_by_xi[round(xi, 10)].append((muq, T))
        cep_mu = cep_mu_by_xi[round(xi, 10)]
        if muq > cep_mu:
            crossover_excluded.append({
                "xi": xi,
                "mu_MeV": muq,
                "T_MeV": T,
                "rho": rho,
                "cep_mu_proxy_MeV": cep_mu,
                "cep_mu_source": cep_mu_source_by_xi[round(xi, 10)],
                "filter_status": "excluded_response_peak_above_CEP",
            })
            crossover_filter_rows.append({
                "xi": xi,
                "mu_MeV": muq,
                "T_MeV": T,
                "rho": rho,
                "cep_mu_proxy_MeV": cep_mu,
                "cep_mu_source": cep_mu_source_by_xi[round(xi, 10)],
                "physical_status": "excluded_response_peak_above_CEP",
                "reason": "mu_gt_mu_CEP_first_order_region_has_Maxwell_only",
            })
            continue
        crossover_plot.append({
            "xi": xi,
            "muq_MeV": muq,
            "mu_MeV": muq,
            "T_MeV": T,
            "rho": rho,
        })
        crossover_filter_rows.append({
            "xi": xi,
            "mu_MeV": muq,
            "T_MeV": T,
            "rho": rho,
            "cep_mu_proxy_MeV": cep_mu,
            "cep_mu_source": cep_mu_source_by_xi[round(xi, 10)],
            "physical_status": "retained_crossover",
            "reason": "mu_le_mu_CEP_proxy",
        })

    # Keep the native sampling gap visible in the audit package.  This is a
    # source-resolution diagnostic, not a license to connect the gap.
    crossover_sampling_gap_rows: list[dict[str, object]] = []
    for xi_key in sorted(crossover_support_by_xi):
        values = sorted(crossover_support_by_xi[xi_key])
        cep_mu = cep_mu_by_xi[xi_key]
        retained = [item for item in values if item[0] <= cep_mu]
        excluded = [item for item in values if item[0] > cep_mu]
        last_retained = retained[-1] if retained else (None, None)
        first_excluded = excluded[0] if excluded else (None, None)
        low, high = cep_bracket_by_xi.get(xi_key, (None, None))
        native_gap = (
            first_excluded[0] - last_retained[0]
            if first_excluded[0] is not None and last_retained[0] is not None
            else None
        )
        crossover_sampling_gap_rows.append({
            "xi": xi_key,
            "muq_CEP_proxy_MeV": cep_mu,
            "T_CEP_bracket_low_MeV": low,
            "T_CEP_bracket_high_MeV": high,
            "last_retained_muq_MeV": last_retained[0],
            "last_retained_T_MeV": last_retained[1],
            "first_excluded_muq_MeV": first_excluded[0],
            "first_excluded_T_MeV": first_excluded[1],
            "native_mu_gap_MeV": native_gap,
            "interpretation": "native crossover sampling gap across CEP proxy; no implicit connector",
        })

    grid_by_coordinate: defaultdict[tuple[float, float], set[str]] = defaultdict(set)
    xi_unresolved: defaultdict[float, set[str]] = defaultdict(set)
    unresolved_grid_rows = 0
    for row in grid:
        xi = finite(row, "xi")
        axis = row["axis"].strip().lower()
        if is_true(row["converged"]):
            continue
        unresolved_grid_rows += 1
        reason = row.get("reason", "unknown")
        if axis == "xi":
            xi_unresolved[round(xi, 10)].add(f"xi:{reason}")
            continue
        T = finite(row, "T_MeV")
        grid_by_coordinate[key(xi, T)].add(f"{axis}:{reason}")

    # An unresolved grid row can still carry useful refinement/geometry
    # diagnostics, but it is never a boundary point unless the source
    # boundary CSV contains a converged row.  Keep this distinction explicit
    # for the no-triangulation review package.
    unresolved_diagnostic_fields = (
        "left",
        "right",
        "midpoint",
        "position_error_MeV",
        "density_error",
        "maxwell_area",
    )
    unresolved_diagnostic_rows: list[dict[str, object]] = []
    for row in grid:
        if is_true(row["converged"]):
            continue
        xi_value = finite(row, "xi")
        T_value = optional_finite(row, "T_MeV")
        values = {field: optional_finite(row, field) for field in unresolved_diagnostic_fields}
        present = [field for field, value in values.items() if value is not None]
        if len(present) == len(unresolved_diagnostic_fields):
            summary_status = "complete_candidate_diagnostics"
        elif present:
            summary_status = "partial_candidate_diagnostics"
        else:
            summary_status = "no_candidate_diagnostics"
        unresolved_diagnostic_rows.append({
            "axis": row["axis"],
            "xi": xi_value,
            "T_MeV": T_value,
            "level": row["level"],
            "reason": row.get("reason", ""),
            "bracket_semantics": f"{row['axis'].strip().lower()}_refinement_level",
            "summary_status": summary_status,
            "boundary_row_present": key(xi_value, T_value) in boundary_keys if T_value is not None else False,
            "candidate_diagnostic_field_count": len(present),
            "candidate_diagnostic_fields": ",".join(present),
            **values,
        })
    unresolved_diagnostic_fields_csv = (
        "axis",
        "xi",
        "T_MeV",
        "level",
        "reason",
        "bracket_semantics",
        "summary_status",
        "boundary_row_present",
        "candidate_diagnostic_field_count",
        "candidate_diagnostic_fields",
        *unresolved_diagnostic_fields,
    )
    write_csv(tables / "grid_unresolved_diagnostics.csv", unresolved_diagnostic_fields_csv, unresolved_diagnostic_rows)

    status_rows: list[dict[str, object]] = []
    for item in boundary_plot:
        reasons = set(grid_by_coordinate.get(item["key"], set()))
        reasons.update(xi_unresolved.get(item["xi"], set()))
        status = "automatic_grid_pass" if not reasons else "|".join(sorted(reasons))
        status_rows.append({
            "xi": item["xi"],
            "T_MeV": item["T_MeV"],
            "mu_MeV": item["mu_MeV"],
            "rho_hadron": item["rho_hadron"],
            "rho_quark": item["rho_quark"],
            "area_residual": item["area_residual"],
            "grid_status": status,
            "grid_unresolved": status != "automatic_grid_pass",
        })
    write_csv(tables / "maxwell_surface_point_status.csv", status_rows[0].keys(), status_rows)
    write_csv(tables / "crossover_physical_filter.csv", crossover_filter_rows[0].keys(), crossover_filter_rows)
    write_csv(tables / "crossover_cep_sampling_gap.csv", crossover_sampling_gap_rows[0].keys(), crossover_sampling_gap_rows)
    if crossover_excluded:
        write_csv(tables / "crossover_excluded_response_peaks.csv", crossover_excluded[0].keys(), crossover_excluded)

    status_counts = Counter(row["grid_status"] for row in status_rows)
    surface_summary = [
        {"surface": "maxwell_coexistence", "rows": len(boundary_plot), "source_rows": len(boundary_plot), "xi_count": len({row["xi"] for row in boundary_plot}), "grid_unresolved_rows": sum(row["grid_unresolved"] for row in status_rows)},
        {"surface": "crossover_pseudocritical", "rows": len(crossover_plot), "source_rows": len(crossover_plot) + len(crossover_excluded), "xi_count": len({row["xi"] for row in crossover_plot}), "grid_unresolved_rows": "not_directly_joined_to_(xi,T)_rho_grid"},
        {"surface": "cep_bracket_projection", "rows": len(cep), "source_rows": len(cep), "xi_count": len({finite(row, "xi") for row in cep}), "grid_unresolved_rows": "bracket_status_preserved"},
    ]
    write_csv(tables / "surface_summary.csv", surface_summary[0].keys(), surface_summary)
    maxwell_max_T_by_xi = {
        xi: max(float(row["T_MeV"]) for row in boundary_plot if row["xi"] == xi)
        for xi in sorted({row["xi"] for row in boundary_plot})
    }
    crossover_min_T_by_xi = {
        xi: min(float(row["T_MeV"]) for row in crossover_plot if row["xi"] == xi)
        for xi in sorted({row["xi"] for row in crossover_plot})
    }
    endpoint_separation_rows = []
    for xi in sorted(set(maxwell_max_T_by_xi) & set(crossover_min_T_by_xi)):
        maxwell_T = maxwell_max_T_by_xi[xi]
        crossover_T = crossover_min_T_by_xi[xi]
        endpoint_separation_rows.append({
            "xi": xi,
            "T_maxwell_max_MeV": maxwell_T,
            "T_crossover_min_MeV": crossover_T,
            "delta_T_crossover_minus_maxwell_MeV": crossover_T - maxwell_T,
            "overlap_in_T": crossover_T <= maxwell_T,
            "interpretation": "diagnostic endpoint separation; not a phase certificate",
        })
    write_csv(
        tables / "crossover_maxwell_endpoint_separation.csv",
        endpoint_separation_rows[0].keys(),
        endpoint_separation_rows,
    )
    write_csv(tables / "maxwell_grid_status_counts.csv", ("grid_status", "row_count", "interpretation"), [
        {"grid_status": status, "row_count": count, "interpretation": "C2 diagnostic status overlay; no production override"}
        for status, count in sorted(status_counts.items())
    ])

    is_visual_closed = args.mode == "visual_closed"
    is_visual_lines = args.mode == "visual_lines"
    hide_diagnostic_overlays = is_visual_closed or is_visual_lines
    schema_version = "pnjl_c2_phase_surfaces_diagnostic_v4" if is_visual_closed else ("pnjl_c2_phase_surfaces_diagnostic_v5" if is_visual_lines else "pnjl_c2_phase_surfaces_diagnostic_v3")
    figure_filename = (
        "c2_phase_surfaces_mu_xi_T_visual_closed.png"
        if is_visual_closed
        else ("c2_phase_surfaces_mu_xi_T_no_triangulation.png" if is_visual_lines else "c2_phase_surfaces_mu_xi_T_diagnostic.png")
    )

    bx = np.asarray([row["xi"] for row in boundary_plot], dtype=float)
    by = np.asarray([row["T_MeV"] for row in boundary_plot], dtype=float)
    bz = np.asarray([row["mu_MeV"] for row in boundary_plot], dtype=float)
    cx = np.asarray([row["xi"] for row in crossover_plot], dtype=float)
    cq = np.asarray([row["muq_MeV"] for row in crossover_plot], dtype=float)
    cz = np.asarray([row["T_MeV"] for row in crossover_plot], dtype=float)
    boundary_triangle_max_dT = max(args.max_dT, 16.0) if is_visual_closed else args.max_dT
    if is_visual_lines:
        boundary_triangles = np.empty((0, 3), dtype=np.int32)
        crossover_triangles = np.empty((0, 3), dtype=np.int32)
    else:
        boundary_triangles = masked_triangles(bx, by, args.max_dxi, boundary_triangle_max_dT)
        crossover_triangles = masked_triangles(cx, cq, args.max_dxi, args.max_dmuq)

    figure = plt.figure(figsize=(12.6, 8.8))
    axis = figure.add_subplot(111, projection="3d")
    if not is_visual_lines:
        if len(boundary_triangles):
            axis.plot_trisurf(bz, bx, by, triangles=boundary_triangles, color="#2f6db0", alpha=0.52, linewidth=0.12, edgecolor="#24527f", shade=True)
        if len(crossover_triangles):
            axis.plot_trisurf(cq, cx, cz, triangles=crossover_triangles, color="#d9792b", alpha=0.35, linewidth=0.10, edgecolor="#9a4b11", shade=True)

    def plot_native_segments(rows: list[dict[str, object]], *, coordinate: str, max_gap: float, color: str, alpha: float, linewidth: float) -> None:
        """Plot only adjacent native support points; never bridge a sampling gap."""
        if not rows:
            return
        ordered = sorted(rows, key=lambda row: float(row[coordinate]))
        segment: list[dict[str, object]] = [ordered[0]]
        for row in ordered[1:]:
            previous = segment[-1]
            if float(row[coordinate]) - float(previous[coordinate]) <= max_gap:
                segment.append(row)
                continue
            if len(segment) >= 2:
                axis.plot(
                    [item["mu_MeV"] for item in segment],
                    [float(segment[0]["xi"])] * len(segment),
                    [item["T_MeV"] for item in segment],
                    color=color,
                    alpha=alpha,
                    linewidth=linewidth,
                )
            segment = [row]
        if len(segment) >= 2:
            axis.plot(
                [item["mu_MeV"] for item in segment],
                [float(segment[0]["xi"])] * len(segment),
                [item["T_MeV"] for item in segment],
                color=color,
                alpha=alpha,
                linewidth=linewidth,
            )
        if is_visual_lines:
            axis.scatter(
                [item["mu_MeV"] for item in ordered],
                [float(ordered[0]["xi"])] * len(ordered),
                [item["T_MeV"] for item in ordered],
                color=color,
                s=4,
                alpha=alpha,
                depthshade=False,
            )

    # Thin generator lines show the actual adaptive support without filling gaps.
    for xi in sorted({row["xi"] for row in boundary_plot}):
        rows = list(row for row in boundary_plot if row["xi"] == xi)
        if is_visual_lines:
            plot_native_segments(rows, coordinate="T_MeV", max_gap=args.max_dT, color="#1f4f86", alpha=0.82, linewidth=0.65)
        else:
            rows = sorted(rows, key=lambda row: row["T_MeV"])
            axis.plot([row["mu_MeV"] for row in rows], [xi] * len(rows), [row["T_MeV"] for row in rows], color="#1f4f86", alpha=0.18, linewidth=0.35)
    for xi in sorted({row["xi"] for row in crossover_plot}):
        rows = list(row for row in crossover_plot if row["xi"] == xi)
        if is_visual_lines:
            plot_native_segments(rows, coordinate="muq_MeV", max_gap=args.max_dmuq, color="#9a4b11", alpha=0.78, linewidth=0.55)
        else:
            rows = sorted(rows, key=lambda row: row["muq_MeV"])
            axis.plot([row["mu_MeV"] for row in rows], [xi] * len(rows), [row["T_MeV"] for row in rows], color="#9a4b11", alpha=0.16, linewidth=0.30)
    if crossover_excluded and not hide_diagnostic_overlays:
        axis.scatter(
            [row["mu_MeV"] for row in crossover_excluded],
            [row["xi"] for row in crossover_excluded],
            [row["T_MeV"] for row in crossover_excluded],
            color="#666666",
            marker="x",
            s=8,
            alpha=0.34,
            depthshade=False,
            label="response peak excluded for $\\mu_q>\\mu_{CEP}$",
        )

    marker_styles = [
        ("rho", "#c43d3d", "rho geometry unresolved"),
        ("temperature", "#7b4fa3", "temperature interpolation unresolved"),
        ("xi", "#2f855a", "xi interpolation unresolved"),
    ]
    if not hide_diagnostic_overlays:
        for token, color, label in marker_styles:
            selected = [row for row in status_rows if token in str(row["grid_status"])]
            if selected:
                axis.scatter([row["mu_MeV"] for row in selected], [row["xi"] for row in selected], [row["T_MeV"] for row in selected], color=color, s=3, alpha=0.28, depthshade=False, label=label)

    cep_projection_count = 0
    for row in cep:
        low = finite(row, "T_bracket_low_MeV")
        high = finite(row, "T_bracket_high_MeV")
        mu = finite(row, "muq_last_first_order_MeV")
        xi = finite(row, "xi")
        if high < low:
            low, high = high, low
        midpoint = 0.5 * (low + high)
        axis.plot([mu, mu], [xi, xi], [low, high], color="#333333", linestyle="--", linewidth=0.65, alpha=0.55)
        axis.scatter([mu], [xi], [midpoint], color="#333333", marker="D", s=10, facecolors="none", depthshade=False)
        cep_projection_count += 1

    legend = [
        Line2D([0], [0], color="#2f6db0", lw=2.0 if is_visual_lines else 8, alpha=0.82 if is_visual_lines else 0.52, label="Maxwell native support" if is_visual_lines else "Maxwell coexistence surface"),
        Line2D([0], [0], color="#d9792b", lw=2.0 if is_visual_lines else 8, alpha=0.78 if is_visual_lines else 0.35, label="crossover native support" if is_visual_lines else "crossover pseudocritical surface"),
        Line2D([0], [0], color="#333333", lw=1.2, linestyle="--", marker="D", markerfacecolor="none", label="CEP bracket projection / midpoint"),
    ]
    if not hide_diagnostic_overlays:
        legend.insert(2, Line2D([0], [0], color="#666666", lw=0, marker="x", alpha=0.65, label="response peak excluded above CEP"))
    axis.legend(handles=legend, loc="upper left", bbox_to_anchor=(0.01, 0.98), framealpha=0.92, fontsize=8)
    axis.set_xlabel(r"$\mu_q$ [MeV]", labelpad=10)
    axis.set_ylabel(r"$\xi$", labelpad=10)
    axis.set_zlabel(r"$T$ [MeV]", labelpad=10)
    title = (
        r"C2 phase curves (no triangulation): $\mu_q$–$\xi$–$T$"
        if is_visual_lines
        else (r"C2 phase surfaces (visual closure only): $\mu_q$–$\xi$–$T$" if is_visual_closed else r"C2 diagnostic phase surfaces: $\mu_q$–$\xi$–$T$")
    )
    axis.set_title(title, pad=12)
    axis.view_init(elev=25, azim=-62)
    try:
        axis.set_box_aspect((1.55, 1.0, 1.25))
    except AttributeError:
        pass
    axis.text2D(
        0.01,
        0.01,
        f"{'no-triangulation diagnostic' if is_visual_lines else ('visualization-only' if is_visual_closed else 'diagnostic-only')} | Maxwell={len(boundary_plot)} | crossover={len(crossover_plot)} | excluded peaks={len(crossover_excluded)} | CEP brackets={cep_projection_count} | unresolved grid={unresolved_grid_rows}; {'native ordered support only' if is_visual_lines else ('finite unresolved rows rendered blue' if is_visual_closed else 'raw support retained')}",
        transform=figure.transFigure,
        fontsize=7,
        color="#444444",
    )
    figure.tight_layout()
    figure_path = figures / figure_filename
    figure.savefig(figure_path, dpi=220, bbox_inches="tight")
    plt.close(figure)

    script_path = Path(__file__).resolve()
    script_digest = sha256(script_path)
    plot_manifest = {
        "schema_version": schema_version,
        "figure_mode": "visualization_only_closed_surface" if is_visual_closed else ("diagnostic_no_triangulation" if is_visual_lines else "diagnostic_global_author_review"),
        "orientation": "x=mu_q, y=xi, z=T",
        "calculation_sha": calculation_sha,
        "solver_called": False,
        "reference_write": False,
        "generator": {"script": str(script_path), "script_sha256": script_digest},
        "surface_policy": {
            "maxwell": "all finite converged boundary rows; unresolved grid statuses rendered with the same blue surface color" if is_visual_closed else ("ordered native boundary support rendered as blue lines; no triangulation" if is_visual_lines else "all finite converged boundary rows; unresolved grid statuses overlaid"),
            "crossover": "finite converged response peaks with mu_q <= the per-xi CEP chemical-potential proxy",
            "excluded_response_peaks": "source crossover rows with mu_q above the per-xi last-first-order CEP proxy; retained in tables but not drawn" if hide_diagnostic_overlays else "source crossover rows with mu_q above the per-xi last-first-order CEP proxy; retained as gray audit points",
            "cep_mu_source": "prefer muq_CEP_MeV when finite; otherwise muq_last_first_order_MeV",
            "cep": "explicit T bracket projected at muq_last_first_order; midpoint is not a CEP",
            "interpolation": "disabled; native ordered support segments only, with gaps left blank" if is_visual_lines else "triangulation only within explicit max gaps; no gap filling",
            "closure_visualization_only": is_visual_closed,
            "triangulation": "disabled" if is_visual_lines else "enabled",
            "unresolved_rendering": "same blue Maxwell surface color; no status markers; no synthetic rows" if is_visual_closed else ("native line support only; unresolved rows omitted" if is_visual_lines else "status-specific markers"),
        },
        "triangle_limits": {"max_dxi": args.max_dxi, "max_dT": boundary_triangle_max_dT, "max_dmuq": args.max_dmuq},
        "visual_closure": {
            "enabled": is_visual_closed,
            "source_rows_only": True,
            "max_boundary_dT_MeV": boundary_triangle_max_dT,
            "note": "no triangulation; native line-gap limits only" if is_visual_lines else "explicit display-only triangulation limit; no synthetic rows and no production certificate effect",
        },
        "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in {
            "boundary": boundary_path,
            "spinodals": spinodal_path,
            "crossover": crossover_path,
            "cep": cep_path,
            "grid_convergence": grid_path,
            "source_manifest": source_manifest_path,
        }.items()},
        "output": {"path": str(figure_path), "sha256": sha256(figure_path)},
    }
    (figures / "plot_manifest.json").write_text(json.dumps(plot_manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    claim_rows = [
        {"claim_id": "global_surface_topology", "claim": "Current C2 finite boundary and crossover outputs form a globally inspectable diagnostic surface", "status": "supported", "evidence": f"figures/{figure_filename}", "boundary": "visual closure is display-only" if is_visual_closed else ("line-only native support view; no surface fill or gap connector" if is_visual_lines else "diagnostic surface only")},
        {"claim_id": "unresolved_status_visible", "claim": "C2 grid unresolved records are visible as status markers rather than silently discarded", "status": "supported" if not hide_diagnostic_overlays else "not_plotted_but_retained_in_table", "evidence": "tables/maxwell_surface_point_status.csv", "boundary": "display mode does not close certificates"},
        {"claim_id": "cep_bracket_semantics", "claim": "CEP information is represented as bracket projections and midpoints, not a certified single-valued CEP line", "status": "supported", "evidence": "figures/plot_manifest.json", "boundary": "mu_q projection uses last first-order point"},
        {"claim_id": "response_peak_physical_filter", "claim": "Response peaks above the per-xi CEP chemical-potential proxy are not plotted as physical crossover", "status": "supported", "evidence": "tables/crossover_physical_filter.csv", "boundary": "CEP proxy remains bracket-derived and is not a certified CEP line"},
        {"claim_id": "crossover_sampling_gap", "claim": "Near-CEP crossover gaps are retained as native sampling gaps and are not implicitly connected", "status": "supported", "evidence": "tables/crossover_cep_sampling_gap.csv", "boundary": "gap size reflects discrete response-peak sampling and CEP proxy filtering"},
        {"claim_id": "crossover_endpoint_separation", "claim": "In every shared xi slice, the lowest retained crossover temperature is above the highest retained Maxwell temperature", "status": "supported", "evidence": "tables/crossover_maxwell_endpoint_separation.csv", "boundary": "discrete endpoint separation only; not a proof of CEP closure"},
        {"claim_id": "phase_reference_promotion", "claim": "The global figure alone authorizes phase-reference promotion", "status": "blocked", "evidence": "tables/maxwell_grid_status_counts.csv", "boundary": "local Maxwell/geometry gates remain unresolved"},
        {"claim_id": "unresolved_maxwell_diagnostic_semantics", "claim": "Unresolved grid rows may retain partial Maxwell diagnostics but are not boundary points", "status": "supported", "evidence": "tables/grid_unresolved_diagnostics.csv", "boundary": "only finite converged boundary rows are plotted as Maxwell support"},
    ]
    write_csv(tables / "claim_ledger.csv", claim_rows[0].keys(), claim_rows)

    decision = {
        "schema_version": schema_version,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "calculation_sha": calculation_sha,
        "solver_called": False,
        "reference_write": False,
        "verdict": "visualization_only_closed_surface_ready" if is_visual_closed else ("diagnostic_no_triangulation_ready" if is_visual_lines else "global_diagnostic_surface_ready_for_author_review"),
        "promotion_effect": "none",
        "counts": {
            "boundary_rows_plotted": len(boundary_plot),
            "crossover_rows_plotted": len(crossover_plot),
            "crossover_rows_excluded": len(crossover_excluded),
            "cep_brackets_projected": cep_projection_count,
            "unresolved_grid_rows": unresolved_grid_rows,
            "maxwell_grid_status_counts": dict(status_counts),
            "unresolved_diagnostic_summary_counts": dict(Counter(row["summary_status"] for row in unresolved_diagnostic_rows)),
        },
    }
    (output_root / "decision.json").write_text(json.dumps(decision, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    mode_readme = (
        "本 v4 版本为 `visualization-only`：所有有限/converged Maxwell boundary 行（包括自动 geometry/interpolation 未闭合行）统一以蓝色绘制，隐藏状态颜色和高于 CEP 的灰色响应峰。它不生成合成数值行、不改变源曲线；仅使用 manifest 明确记录的 16 MeV display-only 三角化上限改善视觉连续性，不把未闭合证书升级为闭合。\n"
        if is_visual_closed else
        "本 v5 模式完全禁用三角化：Maxwell 与 crossover 只绘制原生有序 support 的相邻线段，超过显式采样间隔门限的 gap 保持为空，不填面、不跨 gap、不生成合成点。该图用于区分数据 support 缺失与三角化造成的视觉空洞，不替代逐点 Maxwell/geometry 证书。\n"
        if is_visual_lines else
        "本包用于作者从全局相图检查相结构拓扑，不替代逐点 Maxwell/geometry 证书。C2 grid 未闭合状态以诊断 marker 叠加显示。\n"
    )
    (output_root / "README.md").write_text(
        ("# C2 diagnostic phase surfaces v4 (visual closed)\n\n" if is_visual_closed else ("# C2 diagnostic phase curves v5 (no triangulation)\n\n" if is_visual_lines else "# C2 diagnostic phase surfaces v3\n\n"))
        + "本包的 x 轴为夸克化学势 `mu_q`，y 轴为 `xi`，z 轴为 `T`。Maxwell 与 crossover 的物理筛选遵循同一条规则：同一 `(xi, mu_q)` 不同时赋予两种物理含义；在一阶区 `mu_q > mu_CEP` 的偏导响应峰只作为诊断数据，不绘制为 crossover。\n\n"
        + mode_readme
        + "`tables/crossover_physical_filter.csv` 记录响应峰筛选，`tables/crossover_cep_sampling_gap.csv` 记录每个 xi 切片最后保留点到首个被排除点的原生 mu 间隔；大 xi 处的视觉缺口主要由该离散间隔造成，不能跨 gap 隐式补线。`tables/maxwell_surface_point_status.csv` 保留原始 grid 状态；`tables/grid_unresolved_diagnostics.csv` 记录 unresolved 行是否保留细化层级/几何诊断字段，以及是否存在同坐标的 converged boundary 行。这里的 `left/right/midpoint` 是细化层级元数据，不是 Maxwell 化学势 bracket。CEP 仍保留 bracket 和 midpoint 语义，midpoint 不是 confirmed CEP。\n\n"
        + f"主图见 `figures/{figure_filename}`，证据表见 `tables/claim_ledger.csv`。verdict 为 `{('visualization_only_closed_surface_ready' if is_visual_closed else ('diagnostic_no_triangulation_ready' if is_visual_lines else 'global_diagnostic_surface_ready_for_author_review'))}`，不产生 reference promotion。\n",
        encoding="utf-8",
    )

    output_files = {}
    for path in sorted(output_root.rglob("*")):
        if path.is_file() and path.name != "manifest.json":
            output_files[path.relative_to(output_root).as_posix()] = sha256(path)
    manifest = {
        "schema_version": schema_version,
        "generated_at_utc": decision["generated_at_utc"],
        "calculation_sha": calculation_sha,
        "solver_called": False,
        "reference_write": False,
        "generator": {"script": str(script_path), "script_sha256": script_digest},
        "verdict": decision["verdict"],
        "inputs": plot_manifest["inputs"],
        "counts": decision["counts"],
        "output_files": output_files,
    }
    (output_root / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"verdict": decision["verdict"], "output_root": str(output_root), "counts": decision["counts"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
