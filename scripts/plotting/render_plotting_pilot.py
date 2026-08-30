#!/usr/bin/env python3
"""Render three read-only plotting-contract pilot cases.

The pilot consumes existing CSV/JSON artifacts only. It creates new sibling
directories, refuses to overwrite an existing case, and never calls Julia or
any solver. The output is intentionally a review package, not a final SOP.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import sys
from typing import Any

_SCRIPT_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_SCRIPT_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_PROJECT_ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.tri import Triangulation
import numpy as np

from scripts.plotting.plot_manifest import (
    PROJECT_ROOT,
    build_manifest,
    generator_record,
    input_record,
    output_record,
    runtime_record,
    write_manifest,
)
from scripts.plotting.plot_style import (
    PlotProfile,
    configure_matplotlib,
    figure_size_inches,
    load_profile,
    phase_style,
)


PHASE_REFERENCE = PROJECT_ROOT / "data" / "outputs" / "results" / "pnjl" / "phase_diagram" / "figure4_phase_diagram_prod_v1" / "reference"
MESON_CSV = PROJECT_ROOT / "data" / "outputs" / "results" / "relaxtime" / "meson_mass" / "path_scan" / "freezeout" / "default_baseline_freezeout_xi0_loggrid_1to200_n30.csv"
C1_ROOT = Path(r"D:\Desktop\Julia_RelaxTime_issue130_artifacts\c1_complete_acceptance_31762201725_attempt2\reference")
C1_FILE_NAMES = {
    "boundary": "boundary_issue130_stagec_density_v2_c1_integral_tight_20260813.csv",
    "crossover": "crossover_issue130_stagec_density_v2_c1_integral_tight_20260813.csv",
    "grid_convergence": "phase_grid_convergence_issue130_stagec_density_v2_c1_integral_tight_20260813.csv",
    "cep": "cep_issue130_stagec_density_v2_c1_integral_tight_20260813.csv",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(row for row in handle if row.strip() and not row.startswith("#")))


def finite(row: dict[str, str], field: str) -> float:
    value = float(row[field])
    if not math.isfinite(value):
        raise ValueError(f"non-finite {field}: {value}")
    return value


def is_true(value: str) -> bool:
    return value.strip().lower() == "true"


def exact_xi(row: dict[str, str], target: float) -> bool:
    return math.isclose(float(row["xi"]), target, rel_tol=0.0, abs_tol=1.0e-10)


def refuse_existing(path: Path) -> Path:
    if path.exists():
        raise FileExistsError(f"refusing to overwrite existing pilot case: {path}")
    path.mkdir(parents=True, exist_ok=False)
    return path


def save_outputs(
    figure: Any,
    *,
    output_dir: Path,
    stem: str,
    profile: PlotProfile,
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for fmt in profile.formats:
        output = output_dir / f"{stem}.{fmt}"
        figure.savefig(output, format=fmt, dpi=profile.dpi, bbox_inches="tight", pad_inches=0.08)
        records.append(
            output_record(
                output,
                fmt=fmt,
                dpi=profile.dpi,
                vector=fmt in {"svg", "pdf"},
            )
        )
    plt.close(figure)
    return records


def generator(command: str, font_info: dict[str, str]) -> dict[str, Any]:
    import matplotlib as mpl

    return generator_record(
        Path(__file__),
        command=command,
        runtime=runtime_record({"matplotlib": mpl.__version__, "font": font_info}),
    )


def phase_case(output_root: Path, command: str, suffix: str = "") -> Path:
    profile = load_profile("audit_v1")
    font_info = configure_matplotlib(profile)
    out_dir = refuse_existing(
        output_root
        / "pnjl"
        / "phase_diagram"
        / f"figure4_phase_diagram__plotv1__audit{suffix}"
    )
    boundary_path = PHASE_REFERENCE / "boundary_figure4_phase_diagram_prod_v1_c1_p24t8.csv"
    crossover_path = PHASE_REFERENCE / "crossover_figure4_phase_diagram_prod_v1_c1_p24t8.csv"
    spinodal_path = PHASE_REFERENCE / "spinodals_figure4_phase_diagram_prod_v1_c1_p24t8.csv"
    cep_path = PHASE_REFERENCE / "cep_figure4_phase_diagram_prod_v1_c1_p24t8.csv"
    boundary = read_csv(boundary_path)
    crossover = read_csv(crossover_path)
    spinodals = read_csv(spinodal_path)
    cep = read_csv(cep_path)
    xi_values = [-0.5, -0.25, 0.0, 0.25, 0.5]
    colors = profile.colors
    figure, axis = plt.subplots(figsize=figure_size_inches(profile, "double_column"))
    series: list[dict[str, Any]] = []
    for color_index, xi in enumerate(xi_values):
        color = colors[color_index % len(colors)]
        b_rows = sorted((row for row in boundary if exact_xi(row, xi)), key=lambda row: finite(row, "T_MeV"))
        c_rows = sorted((row for row in crossover if exact_xi(row, xi)), key=lambda row: finite(row, "T_crossover_MeV"))
        s_rows = sorted((row for row in spinodals if exact_xi(row, xi)), key=lambda row: finite(row, "T_MeV"))
        e_rows = [row for row in cep if exact_xi(row, xi)]
        if not b_rows or not c_rows or not s_rows or len(e_rows) != 1:
            raise ValueError(f"phase pilot support is incomplete at xi={xi}")

        first_order_style = phase_style(profile, "first_order")
        axis.plot(
            [3.0 * finite(row, "mu_transition_MeV") for row in b_rows],
            [finite(row, "T_MeV") for row in b_rows],
            color=color,
            linestyle=first_order_style["linestyle"],
            alpha=float(first_order_style["alpha"]),
            label=f"xi={xi:g}",
        )
        series.append(
            {
                "series_id": f"first_order_xi_{xi:g}",
                "state": "first_order",
                "x_field": "mu_transition_MeV",
                "y_field": "T_MeV",
                "support_rule": "raw converged boundary support rows at exact xi; sorted by T_MeV",
                "mask_rule": "nonfinite rows rejected; no row insertion",
                "row_count": len(b_rows),
            }
        )

        crossover_style = phase_style(profile, "crossover")
        axis.plot(
            [3.0 * finite(row, "mu_MeV") for row in c_rows],
            [finite(row, "T_crossover_MeV") for row in c_rows],
            color=color,
            linestyle=crossover_style["linestyle"],
            marker=crossover_style["marker"] if profile.draw_support_points else "",
            markevery=max(1, len(c_rows) // 9) if profile.draw_support_points else None,
            alpha=float(crossover_style["alpha"]),
        )
        series.append(
            {
                "series_id": f"crossover_xi_{xi:g}",
                "state": "crossover",
                "x_field": "mu_MeV",
                "y_field": "T_crossover_MeV",
                "support_rule": "raw converged crossover support rows at exact xi; sorted by T_crossover_MeV",
                "mask_rule": "nonfinite rows rejected; no interpolation",
                "row_count": len(c_rows),
            }
        )

        spinodal_style = phase_style(profile, "spinodal")
        for branch, field in (("hadron", "mu_spinodal_hadron_MeV"), ("quark", "mu_spinodal_quark_MeV")):
            axis.plot(
                [3.0 * finite(row, field) for row in s_rows],
                [finite(row, "T_MeV") for row in s_rows],
                color=color,
                linestyle=spinodal_style["linestyle"],
                alpha=float(spinodal_style["alpha"]),
            )
            series.append(
                {
                    "series_id": f"spinodal_{branch}_xi_{xi:g}",
                    "state": "spinodal",
                    "x_field": field,
                    "y_field": "T_MeV",
                    "support_rule": "raw spinodal support rows at exact xi; sorted by T_MeV",
                    "mask_rule": "nonfinite rows rejected; no interpolation",
                    "row_count": len(s_rows),
                }
            )

        endpoint = e_rows[0]
        cep_style = phase_style(profile, "cep_confirmed")
        axis.plot(
            [finite(endpoint, "muB_CEP_MeV")],
            [finite(endpoint, "T_CEP_MeV")],
            color=color,
            linestyle=cep_style["linestyle"],
            marker=cep_style["marker"],
            markersize=7.0,
            markerfacecolor=color,
            markeredgecolor=color,
            alpha=float(cep_style["alpha"]),
        )
        series.append(
            {
                "series_id": f"cep_confirmed_xi_{xi:g}",
                "state": "cep_confirmed",
                "x_field": "muB_CEP_MeV",
                "y_field": "T_CEP_MeV",
                "support_rule": "raw CEP reference row at exact xi",
                "mask_rule": "empty or nonfinite endpoint rejected",
                "row_count": 1,
            }
        )

    state_handles = [
        Line2D([0], [0], color="#222222", linestyle="-", label="first-order / coexistence"),
        Line2D([0], [0], color="#222222", linestyle="--", label="crossover"),
        Line2D([0], [0], color="#222222", linestyle=":", alpha=0.55, label="spinodal"),
        Line2D([0], [0], color="#222222", marker="o", linestyle="None", markersize=7, label="confirmed CEP"),
    ]
    state_legend = axis.legend(handles=state_handles, loc="upper right", fontsize=7)
    axis.add_artist(state_legend)
    axis.legend(title="xi", loc="lower left", fontsize=7, title_fontsize=7, frameon=False)
    axis.set_xlabel(r"$\mu_B$ [MeV]")
    axis.set_ylabel(r"$T$ [MeV]")
    axis.set_xlim(left=0.0)
    axis.set_ylim(bottom=0.0)
    axis.text(
        0.01,
        0.01,
        "audit pilot | raw support only | mu_B = 3 mu_q where applicable",
        transform=axis.transAxes,
        fontsize=5.5,
        color="#444444",
    )
    outputs = save_outputs(figure, output_dir=out_dir, stem="figure4_phase_diagram_TmuB_audit", profile=profile)
    inputs = [
        input_record(boundary_path, role="calculation_result", schema="figure4_boundary_v1", units={"T_MeV": "MeV", "mu_transition_MeV": "MeV"}),
        input_record(crossover_path, role="calculation_result", schema="figure4_crossover_v1", units={"T_crossover_MeV": "MeV", "mu_MeV": "MeV"}),
        input_record(spinodal_path, role="calculation_result", schema="figure4_spinodal_v1", units={"T_MeV": "MeV", "mu_spinodal_*_MeV": "MeV"}),
        input_record(cep_path, role="calculation_result", schema="figure4_cep_v1", units={"T_CEP_MeV": "MeV", "muB_CEP_MeV": "MeV"}),
    ]
    manifest = build_manifest(
        asset_id="pnjl.figure4.tmuB.audit.pilot_v1",
        figure_family="pnjl_phase_diagram",
        case_slug=f"figure4_phase_diagram__plotv1__audit{suffix}",
        figure_mode="audit",
        semantic_status="audit_phase_state_support",
        style_profile=profile.profile_id,
        publication_scope="internal_review",
        generator=generator(command, font_info),
        inputs=inputs,
        axes=[
            {"field": "muB_MeV", "source_field": "mu_transition_MeV|mu_MeV|mu_spinodal_*_MeV|muB_CEP_MeV", "source_unit": "MeV", "display_unit": "MeV", "transform": "mu_B=3*mu_q for boundary/crossover/spinodal; CEP field is direct"},
            {"field": "T_MeV", "source_field": "T_MeV|T_crossover_MeV|T_CEP_MeV", "source_unit": "MeV", "display_unit": "MeV", "transform": "identity"},
        ],
        series=series,
        outputs=outputs,
        selection_rule="selected xi in {-0.5,-0.25,0,0.25,0.5}; raw support rows only",
        interpolation_policy="none",
        connector_policy="forbidden",
        missing_value_policy="reject nonfinite rows and retain support gaps as gaps",
        validation={"finite": True, "duplicate_keys": True, "support": True, "strict_gate": False, "selected_xi_count": len(xi_values)},
        rendering={"column": "double_column", "figure_size_inches": list(figure_size_inches(profile, "double_column")), "bbox_inches": "tight", "pad_inches": 0.08, "requested_dpi": profile.dpi, "support_visual_policy": profile.support_visual_policy},
    )
    write_manifest(out_dir / "plot_manifest.json", manifest)
    return out_dir


def meson_case(output_root: Path, command: str, suffix: str = "") -> Path:
    profile = load_profile("strict_origin_like_v1")
    font_info = configure_matplotlib(profile)
    out_dir = refuse_existing(
        output_root
        / "relaxtime"
        / "meson_mass"
        / "path_scan"
        / "freezeout"
        / f"default_baseline_freezeout_xi0__plotv1__strict{suffix}"
    )
    rows = read_csv(MESON_CSV)
    valid = [
        row
        for row in rows
        if is_true(row["equilibrium_converged"])
        and finite(row, "T_MeV") > 0.0
        and all(math.isfinite(float(row[field])) for field in ("sqrt_s_NN_GeV", "M_pi", "M_K"))
    ]
    valid.sort(key=lambda row: finite(row, "path_point_index"))
    if not valid:
        raise ValueError("meson strict pilot has no valid rows")
    x = [finite(row, "sqrt_s_NN_GeV") for row in valid]
    if len(x) != len(set(x)):
        raise ValueError("meson strict pilot has duplicate x support")
    segments: list[list[dict[str, str]]] = []
    for row in valid:
        path_index = int(finite(row, "path_point_index"))
        if not segments or path_index != int(finite(segments[-1][-1], "path_point_index")) + 1:
            segments.append([])
        segments[-1].append(row)
    figure, axis = plt.subplots(figsize=figure_size_inches(profile, "single_column"))
    support_marker = "o" if profile.draw_support_points else ""
    for index, segment in enumerate(segments):
        axis.plot(
            [finite(row, "sqrt_s_NN_GeV") for row in segment],
            [finite(row, "M_pi") for row in segment],
            color=profile.colors[0],
            marker=support_marker,
            label=r"$\pi$" if index == 0 else None,
        )
        axis.plot(
            [finite(row, "sqrt_s_NN_GeV") for row in segment],
            [finite(row, "M_K") for row in segment],
            color=profile.colors[1],
            marker=support_marker,
            label=r"$K$" if index == 0 else None,
        )
    axis.set_xlabel(r"$\sqrt{s_{NN}}$ [GeV]")
    axis.set_ylabel(r"mass [fm$^{-1}$]")
    axis.set_xscale("log")
    axis.legend(frameon=False, loc="best")
    outputs = save_outputs(figure, output_dir=out_dir, stem="freezeout_meson_mass_pi_K_strict", profile=profile)
    manifest = build_manifest(
        asset_id="relaxtime.meson_mass.freezeout.pi_K.strict.pilot_v1",
        figure_family="meson_mass_path_scan",
        case_slug=f"default_baseline_freezeout_xi0__plotv1__strict{suffix}",
        figure_mode="strict",
        semantic_status="confirmed_finite_model_support_only",
        style_profile=profile.profile_id,
        publication_scope="main_text_candidate",
        generator=generator(command, font_info),
        inputs=[input_record(MESON_CSV, role="calculation_result", schema="freezeout_meson_mass_scan_v1", units={"sqrt_s_NN_GeV": "GeV", "T_MeV": "MeV", "muB_MeV": "MeV", "M_pi": "fm^-1", "M_K": "fm^-1"})],
        axes=[
            {"field": "sqrt_s_NN_GeV", "source_field": "sqrt_s_NN_GeV", "source_unit": "GeV", "display_unit": "GeV", "transform": "identity; logarithmic scale"},
            {"field": "mass_fm_inv", "source_field": "M_pi|M_K", "source_unit": "fm^-1", "display_unit": "fm^-1", "transform": "identity"},
        ],
        series=[
            {"series_id": "model_pi", "state": "model_support", "x_field": "sqrt_s_NN_GeV", "y_field": "M_pi", "support_rule": "equilibrium_converged=true and T_MeV>0; raw path rows", "mask_rule": "nonfinite rows excluded; raw-support polyline only; no imputation or connector rows", "row_count": len(valid)},
            {"series_id": "model_K", "state": "model_support", "x_field": "sqrt_s_NN_GeV", "y_field": "M_K", "support_rule": "equilibrium_converged=true and T_MeV>0; raw path rows", "mask_rule": "nonfinite rows excluded; raw-support polyline only; no imputation or connector rows", "row_count": len(valid)},
        ],
        outputs=outputs,
        selection_rule="equilibrium_converged=true and T_MeV>0; source path order sorted by path_point_index; each contiguous support segment rendered separately",
        interpolation_policy="none",
        connector_policy="forbidden",
        missing_value_policy="exclude invalid rows and preserve gaps as visual segment breaks",
        validation={"finite": True, "duplicate_keys": True, "support": True, "strict_gate": True, "input_rows": len(rows), "plotted_rows": len(valid), "excluded_rows": len(rows) - len(valid), "support_segments": len(segments)},
        rendering={"column": "single_column", "figure_size_inches": list(figure_size_inches(profile, "single_column")), "bbox_inches": "tight", "pad_inches": 0.08, "requested_dpi": profile.dpi, "support_visual_policy": profile.support_visual_policy},
    )
    write_manifest(out_dir / "plot_manifest.json", manifest)
    return out_dir


def c1_case(analysis_root: Path, c1_root: Path, command: str, suffix: str = "") -> Path:
    profile = load_profile("candidate_origin_like_v1")
    font_info = configure_matplotlib(profile)
    out_dir = refuse_existing(analysis_root / f"plotting_pilot_c1_phase_surface__plotv1__estimated_midpoint{suffix}")
    paths = {name: c1_root / filename for name, filename in C1_FILE_NAMES.items()}
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(f"missing C1 diagnostic input: {path}")
    boundary = [row for row in read_csv(paths["boundary"]) if is_true(row["converged"])]
    crossover = [row for row in read_csv(paths["crossover"]) if is_true(row["converged"])]
    grid = read_csv(paths["grid_convergence"])
    cep = read_csv(paths["cep"])
    bx = np.asarray([finite(row, "xi") for row in boundary])
    by = np.asarray([finite(row, "T_MeV") for row in boundary])
    bz = np.asarray([finite(row, "mu_transition_MeV") for row in boundary])
    cx = np.asarray([finite(row, "xi") for row in crossover])
    cy = np.asarray([finite(row, "T_crossover_MeV") for row in crossover])
    cz = np.asarray([finite(row, "mu_MeV") for row in crossover])
    triangulation = Triangulation(bx, by)
    triangle_values = triangulation.triangles
    dx = np.ptp(bx[triangle_values], axis=1)
    dy = np.ptp(by[triangle_values], axis=1)
    triangulation.set_mask((dx > 0.04) | (dy > 8.0))
    usable = triangulation.triangles if triangulation.mask is None else triangulation.triangles[~triangulation.mask]

    figure = plt.figure(figsize=figure_size_inches(profile, "double_column"))
    axis = figure.add_subplot(111, projection="3d")
    if len(usable):
        axis.plot_trisurf(bx, by, bz, triangles=usable, color=profile.colors[0], alpha=0.33, linewidth=0.12)
    if profile.draw_support_points:
        axis.scatter(cx, cy, cz, color=profile.colors[1], s=7, alpha=0.55, label="crossover support")
    else:
        crossover_groups: dict[float, list[tuple[float, float]]] = {}
        for row in crossover:
            xi = finite(row, "xi")
            crossover_groups.setdefault(xi, []).append(
                (finite(row, "T_crossover_MeV"), finite(row, "mu_MeV"))
            )
        crossover_style = phase_style(profile, "crossover")
        for xi, values in sorted(crossover_groups.items()):
            values.sort(key=lambda item: item[0])
            axis.plot(
                [xi] * len(values),
                [item[0] for item in values],
                [item[1] for item in values],
                color=profile.colors[1],
                linestyle=crossover_style["linestyle"],
                linewidth=0.4,
                alpha=0.32,
            )
    midpoint_rows = sorted(cep, key=lambda item: finite(item, "xi"))
    if profile.draw_support_points:
        for row in midpoint_rows:
            xi = finite(row, "xi")
            low = finite(row, "T_bracket_low_MeV")
            high = finite(row, "T_bracket_high_MeV")
            z = finite(row, "muq_last_first_order_MeV")
            midpoint = 0.5 * (low + high)
            axis.plot([xi, xi], [low, high], [z, z], color="#222222", linestyle="--", linewidth=0.7, alpha=0.7)
            axis.scatter([xi], [midpoint], [z], color="#222222", marker="D", s=14, facecolors="none")
    elif midpoint_rows:
        xis = [finite(row, "xi") for row in midpoint_rows]
        lows = [finite(row, "T_bracket_low_MeV") for row in midpoint_rows]
        highs = [finite(row, "T_bracket_high_MeV") for row in midpoint_rows]
        midpoints = [0.5 * (low + high) for low, high in zip(lows, highs)]
        zs = [finite(row, "muq_last_first_order_MeV") for row in midpoint_rows]
        for values, linestyle, linewidth, alpha in (
            (lows, "--", 1.15, 0.8),
            (highs, "--", 1.15, 0.8),
            (midpoints, ":", 1.35, 0.8),
        ):
            axis.plot(xis, values, zs, color="#222222", linestyle=linestyle, linewidth=linewidth, alpha=alpha)
    axis.set_xlabel(r"$\xi$")
    axis.set_ylabel(r"$T$ [MeV]")
    axis.set_zlabel(r"$\mu$ [MeV]")
    axis.set_title("C1 phase-surface estimated-midpoint review", fontsize=9, pad=8)
    axis.text2D(
        0.01,
        0.01,
        f"estimated midpoint | internal review only | unresolved grid records retained in manifest ({sum(not is_true(row['converged']) for row in grid)})",
        transform=figure.transFigure,
        fontsize=6.5,
        color="#444444",
    )
    outputs = save_outputs(figure, output_dir=out_dir, stem="c1_phase_surface_estimated_midpoint", profile=profile)
    unresolved = sum(not is_true(row["converged"]) for row in grid)
    manifest = build_manifest(
        asset_id="pnjl.c1.phase_surface.estimated_midpoint.pilot_v1",
        figure_family="pnjl_c1_phase_surface_diagnostic",
        case_slug=f"plotting_pilot_c1_phase_surface__plotv1__estimated_midpoint{suffix}",
        figure_mode="estimated_midpoint",
        semantic_status="diagnostic_only_bracket; midpoint_not_physical_CEP",
        style_profile=profile.profile_id,
        publication_scope="supplement_or_internal_review",
        generator=generator(command, font_info),
        inputs=[
            input_record(paths["boundary"], role="diagnostic_calculation_result", schema="c1_boundary_v2", units={"xi": "dimensionless", "T_MeV": "MeV", "mu_transition_MeV": "MeV"}),
            input_record(paths["crossover"], role="diagnostic_calculation_result", schema="c1_crossover_v2", units={"xi": "dimensionless", "T_crossover_MeV": "MeV", "mu_MeV": "MeV"}),
            input_record(paths["grid_convergence"], role="diagnostic_validation", schema="c1_grid_convergence_v2", units={"T_MeV": "MeV", "xi": "dimensionless"}),
            input_record(paths["cep"], role="diagnostic_bracket", schema="c1_cep_bracket_v2", units={"xi": "dimensionless", "T_bracket_*_MeV": "MeV", "muq_last_first_order_MeV": "MeV"}),
        ],
        axes=[
            {"field": "xi", "source_field": "xi", "source_unit": "dimensionless", "display_unit": "dimensionless", "transform": "identity"},
            {"field": "T_MeV", "source_field": "T_MeV|T_crossover_MeV|T_bracket_*_MeV", "source_unit": "MeV", "display_unit": "MeV", "transform": "identity"},
            {"field": "mu_MeV", "source_field": "mu_transition_MeV|mu_MeV|muq_last_first_order_MeV", "source_unit": "MeV", "display_unit": "MeV", "transform": "identity; quark-chemical-potential convention retained"},
        ],
        series=[
            {"series_id": "c1_maxwell_surface", "state": "diagnostic_surface", "x_field": "xi", "y_field": "T_MeV", "z_field": "mu_transition_MeV", "support_rule": "converged=true boundary rows only", "mask_rule": "triangles masked when dx>0.04 or dT>8 MeV; no gap filling", "row_count": len(boundary)},
            {"series_id": "c1_crossover_support", "state": "diagnostic_surface", "x_field": "xi", "y_field": "T_crossover_MeV", "z_field": "mu_MeV", "support_rule": "converged=true crossover rows only; grouped by exact xi and sorted by T", "mask_rule": "raw-support polylines or audit scatter according to profile; no surface interpolation or gap filling", "row_count": len(crossover)},
            {"series_id": "c1_cep_bracket_midpoint", "state": "estimated_midpoint", "x_field": "xi", "y_field": "T_bracket_low_MeV..T_bracket_high_MeV", "z_field": "muq_last_first_order_MeV", "support_rule": "explicit bracket endpoints from C1 artifact", "mask_rule": "midpoint is 0.5*(low+high) and is not a CEP", "row_count": len(midpoint_rows)},
        ],
        outputs=outputs,
        selection_rule="converged=true rows for surfaces; all explicit CEP bracket rows for diagnostic line policy",
        interpolation_policy="none",
        connector_policy="forbidden",
        missing_value_policy="retain unresolved/support gaps as masked or omitted diagnostic regions",
        validation={"finite": True, "duplicate_keys": True, "support": False, "strict_gate": False, "unresolved_grid_records": unresolved, "midpoint_is_physical_CEP": False, "triangle_count_after_mask": int(len(usable))},
        rendering={"column": "double_column", "figure_size_inches": list(figure_size_inches(profile, "double_column")), "bbox_inches": "tight", "pad_inches": 0.08, "requested_dpi": profile.dpi, "support_visual_policy": profile.support_visual_policy},
        calculation_sha="3c5f6b3c9bd535cff7657364dadb2efc31f2ea48",
        postprocess_sha="fd359e792a89beb5ab12349bba761dc58ee16761",
        source_run_id="31762201725",
    )
    write_manifest(out_dir / "plot_manifest.json", manifest)
    return out_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=PROJECT_ROOT / "data" / "outputs" / "figures")
    parser.add_argument("--analysis-root", type=Path, default=PROJECT_ROOT / "docs" / "analysis")
    parser.add_argument("--c1-root", type=Path, default=C1_ROOT)
    parser.add_argument("--only", choices=("all", "phase", "meson", "c1"), default="all")
    parser.add_argument("--suffix", default="", help="suffix for a new sibling pilot case, e.g. __pilot_v2")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    command = "python " + " ".join(str(item) for item in sys.argv)
    outputs: list[Path] = []
    if args.only in {"all", "phase"}:
        outputs.append(phase_case(args.output_root.resolve(), command, args.suffix))
    if args.only in {"all", "meson"}:
        outputs.append(meson_case(args.output_root.resolve(), command, args.suffix))
    if args.only in {"all", "c1"}:
        outputs.append(c1_case(args.analysis_root.resolve(), args.c1_root.resolve(), command, args.suffix))
    for output in outputs:
        print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
