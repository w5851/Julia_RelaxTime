#!/usr/bin/env python3
"""Render the Issue #130 Figure 4 phase surfaces as a v3 style candidate.

The v3 renderer is display-only.  It reads the immutable derived crossover,
Maxwell and CEP tables, retains the v9 endpoint/low-temperature display
policies, and writes a versioned figure package.  It does not call the solver,
write strict/derived rows, or replace the phase-reference public alias.

Compared with v2, v3 improves presentation only: a balanced orthographic view,
explicit inner-origin labels for all three physical coordinates, an external
legend, a neutral directional-light shading of the two surfaces, and a
high-contrast CEP/``crossover intersect mu_q=0`` line treatment.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.analysis.pnjl.render_issue130_phase_surface_v2 import (
    CALCULATION_SHA,
    DEFAULT_INPUT_ROOT,
    SOURCE_POSTPROCESS_SHA,
    SOURCE_REPLAY_ID,
    SOURCE_RUN_ID,
    SOURCE_RENDER_ROOT,
    append_cep_boundary_display,
    build_surface_quads,
    extend_maxwell_groups_to_temperature_floor,
    file_record,
    group_by_xi,
    input_record,
    number,
    read_csv,
)


DEFAULT_OUTPUT_ROOT = Path("data/outputs/figures/pnjl/phase_reference/issue130_phase_reference_v3")
STYLE_PROFILE = "balanced_orthographic_inner_axes_neutral_light_v3"
CEP_LINE_COLOR = "#f1bf3a"
CEP_LINE_WIDTH = 2.8
CEP_HALO_WIDTH = 5.8
MU_ZERO_COLOR = "#7a1f2e"
MU_ZERO_LINE_WIDTH = 1.6
AXIS_COLOR = "#222222"
DISPLAY_T_FLOOR_MEV = 0.0


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def render_style_metadata() -> dict[str, Any]:
    return {
        "style_profile": STYLE_PROFILE,
        "coordinate_grid": "no_pane_or_grid_inner_axis_triad",
        "coordinate_axes": "inner_origin_cartesian_axes_with_explicit_mu_xi_T_labels",
        "projection": "orthographic",
        "view": {"elevation_deg": 20.0, "azimuth_deg": -55.0},
        "surface_shading": "neutral_directional_light_geometry_only",
        "mu_zero_marker": {
            "enabled": True,
            "legend": False,
            "style": "white_halo_dashed_crossover_plane_intersection",
            "color": MU_ZERO_COLOR,
            "linewidth": MU_ZERO_LINE_WIDTH,
            "annotation": r"crossover $\cap\,(\mu_q=0)$",
        },
        "cartesian_scaffold": {
            "origin": {"mu_q": 0.0, "xi": 0.0, "T": 0.0},
            "planes": [],
            "pane_grid": False,
            "axis_color": AXIS_COLOR,
            "axis_linewidth": 1.5,
        },
        "axis_labels": {
            "mu_q": r"$\mu_q$ [MeV]",
            "xi": r"$\xi$ (dimensionless)",
            "T": "T [MeV]",
        },
        "legend": {"placement": "compact_figure_top_center_outside_data", "columns": 3},
        "formats": ["png"],
        "pdf_emitted": False,
        "cep_boundary": {
            "color": CEP_LINE_COLOR,
            "linewidth": CEP_LINE_WIDTH,
            "halo_linewidth": CEP_HALO_WIDTH,
            "label": "CEP boundary (estimated midpoint)",
        },
    }


def _configure_axis(axis: Any) -> None:
    """Remove every default mplot3d pane/grid/boundary axis."""

    axis.grid(False)
    for axis_object in (axis.xaxis, axis.yaxis, axis.zaxis):
        axis_object._axinfo["grid"]["linewidth"] = 0.0
        axis_object._axinfo["grid"]["color"] = (1.0, 1.0, 1.0, 0.0)
        axis_object.pane.set_facecolor((1.0, 1.0, 1.0, 0.0))
        axis_object.pane.set_edgecolor((1.0, 1.0, 1.0, 0.0))
        axis_object.set_visible(False)
    axis.set_axis_off()


def _label_box() -> dict[str, Any]:
    return {
        "facecolor": "white",
        "edgecolor": "none",
        "alpha": 0.84,
        "pad": 1.5,
    }


def _draw_inner_axes(axis: Any) -> None:
    """Draw only the physical Cartesian triad and its manual tick labels."""

    mu_ticks = [0.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0]
    xi_ticks = [-0.5, -0.25, 0.0, 0.25, 0.5]
    temperature_ticks = [0.0, 25.0, 50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0, 225.0]
    axis_style = {"color": AXIS_COLOR, "linewidth": 1.5, "alpha": 0.95}

    # All three axes meet at (mu_q, xi, T)=(0, 0, 0).
    axis.plot([0.0, 400.0], [0.0, 0.0], [0.0, 0.0], **axis_style)
    axis.plot([0.0, 0.0], [-0.5, 0.5], [0.0, 0.0], **axis_style)
    axis.plot([0.0, 0.0], [0.0, 0.0], [0.0, 225.0], **axis_style)

    for mu in mu_ticks:
        axis.plot([mu, mu], [-0.013, 0.013], [0.0, 0.0], **axis_style)
        if mu != 0.0:
            axis.text(mu, -0.047, -5.0, f"{mu:g}", ha="center", va="top", fontsize=8.5)
    for xi in xi_ticks:
        axis.plot([-5.5, 5.5], [xi, xi], [0.0, 0.0], **axis_style)
        if xi != 0.0:
            axis.text(-11.5, xi, -4.0, f"{xi:g}", ha="right", va="center", fontsize=8.5)
    for temperature in temperature_ticks:
        axis.plot([-5.5, 5.5], [0.0, 0.0], [temperature, temperature], **axis_style)
        if temperature != 0.0:
            axis.text(-12.0, 0.0, temperature, f"{temperature:g}", ha="right", va="center", fontsize=8.5)

    axis.text(-10.5, -0.035, -5.0, "0", ha="right", va="top", fontsize=8.5)
    axis.text(
        411.0,
        0.0,
        0.0,
        r"$\mu_q$ [MeV]",
        ha="left",
        va="center",
        fontsize=12,
    )
    # Put the dimensionless axis label at the positive xi end, outside the
    # origin cluster, so it cannot be mistaken for a fourth physical axis.
    axis.text(
        0.0,
        0.57,
        -2.0,
        r"$\xi$ (dimensionless)",
        ha="left",
        va="top",
        fontsize=11,
    )
    axis.text(
        0.0,
        0.0,
        233.0,
        "T [MeV]",
        ha="center",
        va="bottom",
        fontsize=12,
    )


def _shade_quads(
    quads: Iterable[list[tuple[float, float, float]]],
    base_color: str,
    *,
    alpha: float = 0.72,
) -> list[tuple[float, float, float, float]]:
    """Apply geometry-only Lambert-like shading to each display quad.

    Coordinates are normalized before computing normals so the MeV scale of
    μq/T does not swamp the dimensionless ξ direction.  The light changes only
    the visual tone; it does not encode a physical observable or alter rows.
    """

    import numpy as np
    from matplotlib.colors import to_rgb

    base = np.asarray(to_rgb(base_color), dtype=float)
    light = np.asarray([0.52, -0.38, 0.76], dtype=float)
    light /= np.linalg.norm(light)
    colors: list[tuple[float, float, float, float]] = []
    for quad in quads:
        points = np.asarray(
            [[x / 400.0, y / 1.0, z / 225.0] for x, y, z in quad[:3]],
            dtype=float,
        )
        normal = np.cross(points[1] - points[0], points[2] - points[0])
        norm = np.linalg.norm(normal)
        cosine = abs(float(np.dot(normal / norm, light))) if norm > 0.0 else 0.0
        intensity = 0.76 + 0.24 * cosine
        rgb = np.clip(base * intensity, 0.0, 1.0)
        colors.append((float(rgb[0]), float(rgb[1]), float(rgb[2]), alpha))
    return colors


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
        raise RuntimeError("matplotlib is required for the v3 phase render") from exc

    fig = plt.figure(figsize=(12.0, 8.2), dpi=600)
    axis = fig.add_subplot(111, projection="3d")
    axis.set_proj_type("ortho")
    _configure_axis(axis)
    _draw_inner_axes(axis)

    axis.add_collection3d(
        Poly3DCollection(
            crossover_quads,
            facecolors=_shade_quads(crossover_quads, "#d9792b", alpha=0.80),
            edgecolors="none",
            linewidths=0.0,
        )
    )
    axis.add_collection3d(
        Poly3DCollection(
            maxwell_quads,
            facecolors=_shade_quads(maxwell_quads, "#245b8f", alpha=0.80),
            edgecolors="none",
            linewidths=0.0,
        )
    )

    ordered_cep = sorted(cep_rows, key=lambda row: number(row, "xi"))
    cep_mu = [number(row, "mu_CEP_proxy_MeV") for row in ordered_cep]
    cep_xi = [number(row, "xi") for row in ordered_cep]
    cep_T = [number(row, "T_midpoint_MeV") for row in ordered_cep]
    axis.plot(cep_mu, cep_xi, cep_T, color="white", linewidth=CEP_HALO_WIDTH, alpha=0.95)
    axis.plot(cep_mu, cep_xi, cep_T, color=CEP_LINE_COLOR, linewidth=CEP_LINE_WIDTH)

    ordered_mu_zero = sorted(mu_zero_crossover_rows, key=lambda row: number(row, "xi"))
    mu_zero_mu = [number(row, "mu_MeV") for row in ordered_mu_zero]
    mu_zero_xi = [number(row, "xi") for row in ordered_mu_zero]
    mu_zero_T = [number(row, "T_MeV") for row in ordered_mu_zero]
    axis.plot(mu_zero_mu, mu_zero_xi, mu_zero_T, color="white", linewidth=3.8, alpha=0.95)
    axis.plot(
        mu_zero_mu,
        mu_zero_xi,
        mu_zero_T,
        color=MU_ZERO_COLOR,
        linewidth=MU_ZERO_LINE_WIDTH,
        linestyle=(0, (3, 2)),
        alpha=0.98,
    )
    if ordered_mu_zero:
        # Use a fixed 2-D callout in the reserved white area rather than
        # placing text on the orange surface.  The maroon dashed trace itself
        # remains the physical location cue.
        axis.text2D(
            0.18,
            0.72,
            r"$\mu_q=0$ crossover trace",
            transform=axis.transAxes,
            color=MU_ZERO_COLOR,
            fontsize=9.5,
            bbox=_label_box(),
        )

    legend_handles = [
        Line2D([0], [0], color="#d9792b", linewidth=6, alpha=0.78, label="crossover"),
        Line2D(
            [0],
            [0],
            color="#245b8f",
            linewidth=6,
            alpha=0.78,
            label="first-order (Maxwell)",
        ),
        Line2D([0], [0], color=CEP_LINE_COLOR, linewidth=CEP_LINE_WIDTH, label="CEP boundary (estimated midpoint)"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.935),
        ncol=3,
        frameon=True,
        framealpha=0.88,
        fontsize=9.5,
        borderpad=0.35,
        handlelength=2.4,
        columnspacing=1.4,
    )
    fig.suptitle(
        "PNJL phase structure: crossover and first-order coexistence",
        fontsize=13.5,
        y=0.975,
    )
    axis.set_xlim(-18.0, 420.0)
    axis.set_ylim(-0.62, 0.62)
    axis.set_zlim(-8.0, 238.0)
    axis.set_box_aspect((1.75, 1.0, 1.22))
    axis.view_init(elev=20.0, azim=-55.0)
    # Reserve only a compact header; the 3-D artist otherwise leaves a large
    # empty band between the legend and the physical surfaces.
    fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=0.84)
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
        "asset_id": "pnjl_issue130_phase_reference_v3",
        "figure_family": "pnjl_phase_reference",
        "case_slug": "issue130_phase_reference_v3",
        "figure_version": "v3",
        "asset_kind": "figure",
        "figure_mode": "style_only_display_variant",
        "semantic_status": "visualization_only_style_variant",
        "style_profile": STYLE_PROFILE,
        "publication_scope": "external_display_candidate",
        "generator": {
            **file_record(generator, root),
            "role": "render_generator",
            "command": "render_issue130_phase_surface_v3.py",
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
            "triangulation": False,
            "filled_surface": True,
            "figure_size_inches": [12.0, 8.2],
            "boundary_display": "CEP estimated midpoint with endpoint-constrained crossover closure",
            "low_temperature_display": "Maxwell display-only linear continuation to T=0 from first two native T rows",
        },
        "data_contract": {
            "strict_and_derived_tables_unchanged": True,
            "numeric_source_rows_unchanged": True,
            "solver_called": False,
            "public_alias_replaced": False,
            "mu_zero_marker_added": True,
            "mu_zero_marker_semantics": "crossover_surface_intersection_with_mu_q_zero_plane",
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
        raise FileExistsError(f"refusing to overwrite v3 evidence root: {output_root}")
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
    mu_zero_crossover_rows = [row for row in crossover_rows if abs(number(row, "mu_MeV")) <= 1.0e-10]
    if len(mu_zero_crossover_rows) < 2:
        raise ValueError("crossover input lacks a usable mu_q=0 intersection line")

    crossover_groups, _closure_rows = append_cep_boundary_display(group_by_xi(crossover_rows), cep_rows)
    maxwell_groups, _maxwell_trace = extend_maxwell_groups_to_temperature_floor(
        group_by_xi(maxwell_rows), floor_T_MeV=DISPLAY_T_FLOOR_MEV
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

    output_png = output_root / "phase_surface_render_mu_xi_T.png"
    render_figure(output_png, crossover_quads, maxwell_quads, cep_rows, mu_zero_crossover_rows)

    inputs = [
        input_record(crossover_path, root, "derived_crossover_table"),
        input_record(maxwell_path, root, "derived_maxwell_table"),
        input_record(cep_path, root, "derived_cep_boundary_table"),
        input_record(coverage_path, root, "derived_coverage_mask"),
    ]
    manifest = build_plot_manifest(
        root=root,
        output_png=output_png,
        generator=Path(__file__).resolve(),
        inputs=inputs,
    )
    manifest["surface_summary"] = [
        {"surface": "crossover", **crossover_summary, "display_policy": "cep_boundary_constrained_common_support"},
        {"surface": "first_order_coexistence_maxwell", **maxwell_summary, "display_policy": "bounded_display_extrapolation_to_T0"},
    ]
    manifest["source_render_manifest"] = (
        file_record(root / SOURCE_RENDER_ROOT / "manifest.json", root)
        if (root / SOURCE_RENDER_ROOT / "manifest.json").is_file()
        else None
    )
    _write_json(output_root / "plot_manifest.json", manifest)
    _write_json(
        output_root / "decision.json",
        {
            "verdict": "display_candidate",
            "figure_version": "v3",
            "reference_write": False,
            "runtime_consumption": False,
            "solver_called": False,
            "strict_or_derived_csv_modified": False,
            "public_alias_replaced": False,
            "supersedes": "data/outputs/figures/pnjl/phase_reference/issue130_phase_reference_v2/",
        },
    )
    (output_root / "README.md").write_text(
        "# Issue #130 phase-reference Figure 4 v3\n\n"
        "Display-only style candidate generated from the immutable derived tables. "
        "The three physical axes are mu_q [MeV], xi (dimensionless), and T [MeV], "
        "meeting at the inner origin. Pane/grid artists are disabled. Neutral "
        "directional shading changes visual tone only; it is not a data encoding. "
        "The crossover intersection with mu_q=0 is shown with a white halo and "
        "maroon dashed core, outside the legend. No strict/derived data or public "
        "alias is modified.\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "output_root": str(output_root),
                "png": str(output_png),
                "solver_called": False,
                "strict_and_derived_tables_unchanged": True,
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
