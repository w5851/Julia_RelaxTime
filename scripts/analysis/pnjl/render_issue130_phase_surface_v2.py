#!/usr/bin/env python3
"""Render the Issue #130 phase-surface figure-layer v2.

This is a style-only display variant.  It reuses the same derived tables and
the same display-only CEP/low-temperature policies as the current v9 render,
but uses an inner-origin Cartesian 3-D axis triad, marks the actual crossover
intersection with the mu_q=0 plane outside the legend, and makes the CEP
boundary more prominent.  It never writes any strict/derived numerical table.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.analysis.pnjl.render_issue130_phase_surface_v9 import (
    CALCULATION_SHA,
    SOURCE_POSTPROCESS_SHA,
    SOURCE_REPLAY_ID,
    SOURCE_RUN_ID,
    append_cep_boundary_display,
    build_surface_quads,
    extend_maxwell_groups_to_temperature_floor,
    file_record,
    group_by_xi,
    input_record,
    number,
    read_csv,
)


DEFAULT_INPUT_ROOT = Path("data/reference/pnjl/issue130_phase_reference_v1")
DEFAULT_OUTPUT_ROOT = Path(
    "data/outputs/figures/pnjl/phase_reference/issue130_phase_reference_v2"
)
SOURCE_RENDER_ROOT = Path(
    "docs/analysis/pnjl/phase_reference/issue130_phase_reference_layers_v1/"
    "phase_surface_render_v9"
)
STYLE_PROFILE = "inner_origin_cartesian_3d_mu0_intersection_v2"
CEP_LINE_COLOR = "#202020"
CEP_LINE_WIDTH = 3.2
CEP_HALO_WIDTH = 5.2
MU_ZERO_COLOR = "#6b2e2e"
MU_ZERO_LINE_WIDTH = 1.4
DISPLAY_T_FLOOR_MEV = 0.0


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def render_style_metadata() -> dict[str, Any]:
    """Return the explicit, data-independent v2 rendering contract."""

    return {
        "style_profile": STYLE_PROFILE,
        "coordinate_grid": "no_pane_or_grid_inner_axis_triad",
        "coordinate_axes": "inner_origin_cartesian_axes_with_manual_ticks_and_labels",
        "projection": "orthographic",
        "mu_zero_marker": {
            "enabled": True,
            "legend": False,
            "style": "thin_dashed_crossover_plane_intersection",
            "color": MU_ZERO_COLOR,
            "linewidth": MU_ZERO_LINE_WIDTH,
            "annotation": r"crossover $\cap\,(\mu_q=0)$",
        },
        "cartesian_scaffold": {
            "origin": {"mu_q": 0.0, "xi": 0.0, "T": 0.0},
            "planes": [],
            "pane_grid": False,
            "axis_color": "#202020",
            "axis_linewidth": 1.6,
        },
        "cep_boundary": {
            "color": CEP_LINE_COLOR,
            "linewidth": CEP_LINE_WIDTH,
            "halo_linewidth": CEP_HALO_WIDTH,
            "label": "CEP boundary (estimated midpoint)",
        },
        "formats": ["png"],
        "pdf_emitted": False,
    }


def _configure_axis_only(axis: Any) -> None:
    """Hide mplot3d's boundary axes, panes and grids.

    The visible coordinate frame is drawn explicitly at the physical origin;
    leaving mplot3d's default boundary frame enabled would add a second,
    misleading set of axes on the box edges.
    """

    axis.grid(False)
    for axis_object in (axis.xaxis, axis.yaxis, axis.zaxis):
        # mplot3d stores pane/grid styling in _axinfo.  This is the stable
        # rendering hook across the Matplotlib versions used by this project.
        axis_object._axinfo["grid"]["linewidth"] = 0.0
        axis_object._axinfo["grid"]["color"] = (1.0, 1.0, 1.0, 0.0)
        axis_object._axinfo["axisline"]["linewidth"] = 1.2
        axis_object._axinfo["axisline"]["color"] = "#202020"
        axis_object.line.set_color("#202020")
        axis_object.line.set_linewidth(1.2)
        axis_object.pane.set_facecolor((1.0, 1.0, 1.0, 0.0))
        axis_object.pane.set_edgecolor((1.0, 1.0, 1.0, 0.0))
        axis_object.set_visible(False)
    # mplot3d can still draw a boundary z-axis through its private axis
    # artist even after individual Axis objects are hidden.  Disable that
    # artist as a whole; manual inner-axis lines/text remain ordinary artists.
    axis.set_axis_off()


def _draw_inner_cartesian_axes(axis: Any) -> None:
    """Draw inner-origin Cartesian axes without panes or background grids."""

    mu_ticks = [0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0]
    xi_ticks = [-0.5, -0.25, 0.0, 0.25, 0.5]
    temperature_ticks = [0.0, 25.0, 50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0]
    # The three axes meet at the physical origin (mu_q=0, xi=0, T=0).
    axis_style = {"color": "#202020", "linewidth": 1.6, "alpha": 0.95}
    axis.plot([0.0, 400.0], [0.0, 0.0], [0.0, 0.0], **axis_style)
    axis.plot([0.0, 0.0], [-0.5, 0.5], [0.0, 0.0], **axis_style)
    axis.plot([0.0, 0.0], [0.0, 0.0], [0.0, 225.0], **axis_style)

    # Manual ticks keep the frame visible after the default boundary axes are
    # hidden.  Their cross-segment lengths are small in the corresponding
    # physical coordinate so they do not form a replacement pane grid.
    for mu in mu_ticks:
        axis.plot([mu, mu], [-0.014, 0.014], [0.0, 0.0], **axis_style)
        if mu != 0.0:
            axis.text(mu, -0.045, -3.5, f"{mu:g}", ha="center", va="top", fontsize=9)
    for xi in xi_ticks:
        axis.plot([-6.0, 6.0], [xi, xi], [0.0, 0.0], **axis_style)
        if xi != 0.0:
            axis.text(-12.0, xi, -3.5, f"{xi:g}", ha="right", va="center", fontsize=9)
    for temperature in temperature_ticks:
        axis.plot([-6.0, 6.0], [0.0, 0.0], [temperature, temperature], **axis_style)
        if temperature != 0.0:
            axis.text(-12.0, 0.0, temperature, f"{temperature:g}", ha="right", va="center", fontsize=9)

    axis.text(-10.0, -0.025, -3.5, "0", ha="right", va="top", fontsize=9)
    axis.text(405.0, 0.0, 0.0, r"$\mu_q$ [MeV]", ha="left", va="center", fontsize=12)
    axis.text(0.0, 0.54, 0.0, r"$\xi$", ha="center", va="center", fontsize=12)
    axis.text(0.0, 0.0, 230.0, "T [MeV]", ha="center", va="bottom", fontsize=12)


def render_figure(
    output_png: Path,
    crossover_quads: list[list[tuple[float, float, float]]],
    maxwell_quads: list[list[tuple[float, float, float]]],
    cep_rows: list[dict[str, str]],
    mu_zero_crossover_rows: list[dict[str, str]],
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    except ImportError as exc:  # pragma: no cover - environment-specific
        raise RuntimeError("matplotlib is required for the v2 phase render") from exc

    fig = plt.figure(figsize=(11, 8), dpi=600)
    axis = fig.add_subplot(111, projection="3d")
    # Use the same non-perspective projection as a textbook Cartesian frame:
    # the three axis directions remain orthogonal instead of converging toward
    # a perspective vanishing point.
    axis.set_proj_type("ortho")
    _configure_axis_only(axis)
    _draw_inner_cartesian_axes(axis)

    crossover_collection = Poly3DCollection(
        crossover_quads,
        facecolors="#d9792b",
        edgecolors="none",
        linewidths=0.0,
        alpha=0.60,
    )
    maxwell_collection = Poly3DCollection(
        maxwell_quads,
        facecolors="#245b8f",
        edgecolors="none",
        linewidths=0.0,
        alpha=0.60,
    )
    axis.add_collection3d(crossover_collection)
    axis.add_collection3d(maxwell_collection)

    ordered_cep = sorted(cep_rows, key=lambda row: number(row, "xi"))
    cep_mu = [number(row, "mu_CEP_proxy_MeV") for row in ordered_cep]
    cep_xi = [number(row, "xi") for row in ordered_cep]
    cep_T = [number(row, "T_midpoint_MeV") for row in ordered_cep]
    # A light under-stroke keeps the boundary legible where it crosses either
    # translucent surface without introducing another physical surface.
    axis.plot(cep_mu, cep_xi, cep_T, color="white", linewidth=CEP_HALO_WIDTH, alpha=0.9)
    axis.plot(
        cep_mu,
        cep_xi,
        cep_T,
        color=CEP_LINE_COLOR,
        linewidth=CEP_LINE_WIDTH,
        label="CEP boundary (estimated midpoint)",
    )
    # Mark the actual crossover/μq=0 intersection in the figure itself, but
    # keep it out of the legend because it is a derived reference intersection
    # rather than a separate physical phase surface.
    ordered_mu_zero = sorted(mu_zero_crossover_rows, key=lambda row: number(row, "xi"))
    axis.plot(
        [number(row, "mu_MeV") for row in ordered_mu_zero],
        [number(row, "xi") for row in ordered_mu_zero],
        [number(row, "T_MeV") for row in ordered_mu_zero],
        color=MU_ZERO_COLOR,
        linewidth=MU_ZERO_LINE_WIDTH,
        linestyle=(0, (4, 2)),
        alpha=0.95,
    )
    if ordered_mu_zero:
        label_row = ordered_mu_zero[len(ordered_mu_zero) // 2]
        axis.text(
            number(label_row, "mu_MeV") + 5.0,
            number(label_row, "xi"),
            number(label_row, "T_MeV") + 3.0,
            r"crossover $\cap\,(\mu_q=0)$",
            color=MU_ZERO_COLOR,
            fontsize=10,
        )

    legend_handles = [
        Line2D([0], [0], color="#d9792b", linewidth=5, alpha=0.60, label="crossover"),
        Line2D(
            [0],
            [0],
            color="#245b8f",
            linewidth=5,
            alpha=0.60,
            label="first-order coexistence (Maxwell)",
        ),
        Line2D(
            [0],
            [0],
            color=CEP_LINE_COLOR,
            linewidth=CEP_LINE_WIDTH,
            label="CEP boundary (estimated midpoint)",
        ),
    ]
    axis.legend(handles=legend_handles, loc="upper left", frameon=True)
    axis.set_title("PNJL phase structure: crossover and first-order coexistence")
    axis.set_xlim(0.0, 400.0)
    axis.set_ylim(-0.5, 0.5)
    axis.set_zlim(0.0, 225.0)
    axis.view_init(elev=23, azim=-62)
    fig.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=600, bbox_inches="tight")
    plt.close(fig)


def build_plot_manifest(
    *,
    root: Path,
    output_png: Path,
    generator: Path,
    inputs: list[dict[str, Any]],
) -> dict[str, Any]:
    return {
        "schema_version": "plot_manifest_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "asset_id": "pnjl_issue130_figure4_phase_diagram_prod_v2",
        "figure_family": "pnjl_phase_diagram",
        "case_slug": "figure4_phase_diagram_prod_v2",
        "figure_version": "v2",
        "asset_kind": "figure",
        "figure_mode": "style_only_display_variant",
        "semantic_status": "visualization_only_style_variant",
        "style_profile": STYLE_PROFILE,
        "publication_scope": "external_display_candidate",
        "generator": {
            **file_record(generator, root),
            "role": "render_generator",
            "command": "render_issue130_phase_surface_v2.py",
            "runtime": {
                "python": platform.python_version(),
                "platform": platform.platform(),
                "executable": sys.executable,
            },
        },
        "calculation_sha": CALCULATION_SHA,
        "postprocess_sha": SOURCE_POSTPROCESS_SHA,
        "source_run_id": SOURCE_RUN_ID,
        "replay_run_id": SOURCE_REPLAY_ID,
        "source_render_layer": "phase_surface_render_v9",
        "inputs": inputs,
        "axes": [
            {"field": "mu_q", "source_unit": "MeV", "display_unit": "MeV", "transform": "identity"},
            {"field": "xi", "source_unit": "dimensionless", "display_unit": "dimensionless", "transform": "identity"},
            {"field": "T", "source_unit": "MeV", "display_unit": "MeV", "transform": "identity"},
        ],
        "rendering": {
            **render_style_metadata(),
            "projection": "orthographic_cartesian_surface",
            "triangulation": False,
            "filled_surface": True,
            "figure_size_inches": [11.0, 8.0],
            "boundary_display": "CEP estimated midpoint with endpoint-constrained crossover closure",
            "low_temperature_display": "Maxwell display-only linear continuation to T=0 from first two native T rows",
        },
        "data_contract": {
            "strict_and_derived_tables_unchanged": True,
            "numeric_source_rows_unchanged": True,
            "solver_called": False,
            "public_alias_replaced": False,
            "mu_zero_marker_added": True,
        },
        "outputs": [
            {**file_record(output_png, root), "format": "png", "dpi": 600, "vector": False},
        ],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--grid-points", type=int, default=160)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.repo_root.resolve()
    input_root = (root / args.input_root).resolve()
    output_root = (root / args.output_root).resolve()
    if output_root.exists():
        raise FileExistsError(f"refusing to overwrite v2 evidence root: {output_root}")
    if args.grid_points < 8:
        raise ValueError("--grid-points must be at least 8")

    crossover_path = input_root / "derived/tables/crossover_surface_derived_reference_v1.csv"
    maxwell_path = input_root / "derived/tables/maxwell_surface_derived_reference_v1.csv"
    cep_path = input_root / "derived/tables/cep_boundary_derived_reference_v1.csv"
    coverage_path = input_root / "derived/tables/surface_coverage_mask.csv"
    for path in (crossover_path, maxwell_path, cep_path, coverage_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    crossover_rows = read_csv(crossover_path)
    maxwell_rows = read_csv(maxwell_path)
    cep_rows = read_csv(cep_path)
    read_csv(coverage_path)
    mu_zero_crossover_rows = [
        row for row in crossover_rows if abs(number(row, "mu_MeV")) <= 1.0e-10
    ]
    if len(mu_zero_crossover_rows) < 2:
        raise ValueError("crossover input lacks a usable mu_q=0 intersection line")
    crossover_groups, _closure_rows = append_cep_boundary_display(group_by_xi(crossover_rows), cep_rows)
    maxwell_groups, maxwell_extrapolation_rows = extend_maxwell_groups_to_temperature_floor(
        group_by_xi(maxwell_rows),
        floor_T_MeV=DISPLAY_T_FLOOR_MEV,
    )
    crossover_quads, crossover_summary = build_surface_quads(
        crossover_groups,
        axis="mu_MeV",
        value_field="T_MeV",
        grid_points=args.grid_points,
        max_gap=26.0,
    )
    maxwell_quads, maxwell_summary = build_surface_quads(
        maxwell_groups,
        axis="T_MeV",
        value_field="mu_MeV",
        grid_points=args.grid_points,
        max_gap=18.0,
    )

    figure_root = output_root
    output_png = figure_root / "phase_surface_render_mu_xi_T.png"
    render_figure(
        output_png,
        crossover_quads,
        maxwell_quads,
        cep_rows,
        mu_zero_crossover_rows,
    )

    source_inputs = [
        input_record(crossover_path, root, "derived_crossover_table"),
        input_record(maxwell_path, root, "derived_maxwell_table"),
        input_record(cep_path, root, "derived_cep_boundary_table"),
        input_record(coverage_path, root, "derived_coverage_mask"),
    ]
    plot_manifest_path = figure_root / "plot_manifest.json"
    summary_rows = [
        {
            "surface": "crossover",
            **crossover_summary,
            "display_policy": "cep_boundary_constrained_common_support",
        },
        {
            "surface": "first_order_coexistence_maxwell",
            **maxwell_summary,
            "display_policy": "bounded_display_extrapolation_to_T0",
            "display_T_floor_MeV": DISPLAY_T_FLOOR_MEV,
            "extrapolated_xi_count": sum(bool(row["synthetic_row"]) for row in maxwell_extrapolation_rows),
            "max_extrapolation_gap_MeV": max(
                float(row["extension_gap_MeV"]) for row in maxwell_extrapolation_rows
            ),
        },
    ]
    source_render_manifest = root / SOURCE_RENDER_ROOT / "manifest.json"
    plot_manifest = build_plot_manifest(
        root=root,
        output_png=output_png,
        generator=Path(__file__).resolve(),
        inputs=source_inputs,
    )
    plot_manifest.update(
        {
            "source_render_layer": "phase_surface_render_v9",
            "source_render_manifest": file_record(source_render_manifest, root) if source_render_manifest.is_file() else None,
            "source_display_policy": "v9_endpoint_and_low_temperature_display_replay",
            "solver_called": False,
            "strict_and_derived_tables_unchanged": True,
            "numeric_source_rows_unchanged": True,
            "public_alias_replaced": False,
            "mu_zero_marker_added": True,
            "mu_zero_marker_semantics": "crossover_surface_intersection_with_mu_q_zero_plane",
            "mu_zero_intersection_row_count": len(mu_zero_crossover_rows),
            "formal_output_root": "data/outputs/figures/pnjl/phase_reference/issue130_phase_reference_v2",
            "surface_summary": summary_rows,
        }
    )
    _write_json(plot_manifest_path, plot_manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
