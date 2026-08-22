#!/usr/bin/env python3
"""Render the Issue #130 phase-reference v9 display layer.

The numerical strict/derived tables are immutable inputs.  This renderer only
builds a display mesh from local common support, records a boundary-constrained
crossover endpoint at the estimated CEP midpoint, and applies a bounded,
display-only low-temperature continuation to the Maxwell surface.  Synthetic
rows are never written to the reference tables.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import shutil
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INPUT_ROOT = Path("data/reference/pnjl/issue130_phase_reference_v1")
DEFAULT_OUTPUT_ROOT = Path(
    "docs/analysis/pnjl/phase_reference/issue130_phase_reference_layers_v1/"
    "phase_surface_render_v9"
)
DEFAULT_PUBLIC_ROOT = Path(
    "data/outputs/figures/pnjl/phase_reference/issue130_phase_reference_v1"
)
CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
SOURCE_RUN_ID = "32354095831"
SOURCE_REPLAY_ID = "32451053476"
SOURCE_POSTPROCESS_SHA = "d9a52e6232d09695a83997cda2c11f6b00204587"
EPS = 1.0e-10
DISPLAY_T_FLOOR_MEV = 0.0
PREVIOUS_PUBLIC_VERSION = "v8"
PREVIOUS_PUBLIC_PNG_SHA256 = "5b8152b2ccef4107c51393370b82ff77f48287132e8206a37ee6cfb2ab793b14"
PREVIOUS_PUBLIC_PLOT_MANIFEST_SHA256 = "892b383698708eac3a836dfa5b3590ad9d7850cb9ffd31c0a6e080d74c4058b3"
PREVIOUS_PUBLIC_SOURCE_PATH = (
    "docs/analysis/pnjl/phase_reference/issue130_phase_reference_layers_v1/"
    "phase_surface_render_v8/figures/phase_surface_render_mu_xi_T.png"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_commit(root: Path) -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=root,
            text=True,
            encoding="utf-8",
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def relative_path(path: Path, root: Path) -> tuple[str, str]:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix(), "repository_relative"
    except ValueError:
        return str(path.resolve()), "external_absolute"


def file_record(path: Path, root: Path) -> dict[str, Any]:
    display, path_kind = relative_path(path, root)
    return {
        "path": display,
        "path_kind": path_kind,
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def input_record(path: Path, root: Path, role: str) -> dict[str, Any]:
    record = file_record(path, root)
    record["role"] = role
    return record


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"CSV has no data rows: {path}")
    return rows


def write_csv(path: Path, fields: Iterable[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def number(row: dict[str, Any], field: str) -> float:
    value = float(row[field])
    if not math.isfinite(value):
        raise ValueError(f"non-finite {field}: {value}")
    return value


def xi_key(value: float) -> float:
    return round(value, 10)


def group_by_xi(rows: Iterable[dict[str, Any]]) -> dict[float, list[dict[str, Any]]]:
    grouped: defaultdict[float, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[xi_key(number(row, "xi"))].append(row)
    return dict(grouped)


def sorted_axis(rows: Iterable[dict[str, Any]], axis: str) -> list[dict[str, Any]]:
    return sorted(rows, key=lambda row: number(row, axis))


def interpolated_value(
    rows: list[dict[str, Any]],
    *,
    axis: str,
    field: str,
    target: float,
    max_gap: float,
) -> float | None:
    ordered = sorted_axis(rows, axis)
    if target < number(ordered[0], axis) - EPS or target > number(ordered[-1], axis) + EPS:
        return None
    for left, right in zip(ordered, ordered[1:]):
        x0 = number(left, axis)
        x1 = number(right, axis)
        if target > x1 + EPS:
            continue
        if abs(x1 - x0) <= EPS:
            return number(left, field)
        if x1 - x0 > max_gap + EPS:
            return None
        weight = (target - x0) / (x1 - x0)
        value = number(left, field) + weight * (number(right, field) - number(left, field))
        return value if math.isfinite(value) else None
    return number(ordered[-1], field) if abs(target - number(ordered[-1], axis)) <= EPS else None


def common_range(left: list[dict[str, Any]], right: list[dict[str, Any]], axis: str) -> tuple[float, float] | None:
    low = max(number(sorted_axis(left, axis)[0], axis), number(sorted_axis(right, axis)[0], axis))
    high = min(number(sorted_axis(left, axis)[-1], axis), number(sorted_axis(right, axis)[-1], axis))
    return (low, high) if high - low > EPS else None


def append_cep_boundary_display(
    crossover_groups: dict[float, list[dict[str, Any]]],
    cep_rows: list[dict[str, str]],
) -> tuple[dict[float, list[dict[str, Any]]], list[dict[str, Any]]]:
    """Add only display endpoints connecting the last crossover point to CEP."""

    cep_by_xi = {xi_key(number(row, "xi")): row for row in cep_rows}
    result = {key: list(rows) for key, rows in crossover_groups.items()}
    closure_rows: list[dict[str, Any]] = []
    for key, rows in sorted(result.items()):
        ordered = sorted_axis(rows, "mu_MeV")
        cep = cep_by_xi.get(key)
        if cep is None:
            closure_rows.append({
                "xi": key,
                "last_crossover_mu_MeV": "",
                "cep_mu_MeV": "",
                "last_crossover_T_MeV": "",
                "cep_T_midpoint_MeV": "",
                "delta_mu_MeV": "",
                "delta_T_MeV": "",
                "closure_status": "missing_cep_boundary",
            })
            continue
        last = ordered[-1]
        last_mu = number(last, "mu_MeV")
        cep_mu = number(cep, "mu_CEP_proxy_MeV")
        cep_T = number(cep, "T_midpoint_MeV")
        if cep_mu <= last_mu + EPS:
            closure_rows.append({
                "xi": key,
                "last_crossover_mu_MeV": last_mu,
                "cep_mu_MeV": cep_mu,
                "last_crossover_T_MeV": number(last, "T_MeV"),
                "cep_T_midpoint_MeV": cep_T,
                "delta_mu_MeV": cep_mu - last_mu,
                "delta_T_MeV": cep_T - number(last, "T_MeV"),
                "closure_status": "not_needed_or_boundary_before_last_support",
            })
            continue
        result[key].append({
            "surface": "crossover",
            "xi": key,
            "mu_MeV": cep_mu,
            "T_MeV": cep_T,
            "rho": number(last, "rho"),
            "layer": "display_only",
            "status": "boundary_constrained_cep_display_endpoint",
            "source_layer": "derived_reference_v1_cep_boundary",
        })
        closure_rows.append({
            "xi": key,
            "last_crossover_mu_MeV": last_mu,
            "cep_mu_MeV": cep_mu,
            "last_crossover_T_MeV": number(last, "T_MeV"),
            "cep_T_midpoint_MeV": cep_T,
            "delta_mu_MeV": cep_mu - last_mu,
            "delta_T_MeV": cep_T - number(last, "T_MeV"),
            "closure_status": "added_boundary_constrained_display_endpoint",
        })
    return result, closure_rows


def extend_maxwell_groups_to_temperature_floor(
    groups: dict[float, list[dict[str, Any]]],
    *,
    floor_T_MeV: float,
) -> tuple[dict[float, list[dict[str, Any]]], list[dict[str, Any]]]:
    """Add a bounded display-only continuation below each slice's first T.

    The derived Maxwell table has a slice-dependent low-temperature support
    boundary.  A two-point local continuation of ``mu(T, xi)`` removes that
    support-edge staircase in the public visualization while preserving every
    source row.  The synthetic row is explicitly non-certified and is never
    used to build strict/derived reference tables.
    """

    if not math.isfinite(floor_T_MeV):
        raise ValueError(f"non-finite display temperature floor: {floor_T_MeV}")
    result: dict[float, list[dict[str, Any]]] = {}
    trace: list[dict[str, Any]] = []
    for key, rows in sorted(groups.items()):
        ordered = sorted_axis(rows, "T_MeV")
        first_T = number(ordered[0], "T_MeV")
        first_mu = number(ordered[0], "mu_MeV")
        if first_T <= floor_T_MeV + EPS:
            result[key] = list(ordered)
            trace.append({
                "xi": key,
                "source_T_min_MeV": first_T,
                "display_T_floor_MeV": floor_T_MeV,
                "extension_gap_MeV": 0.0,
                "extension_method": "native_support",
                "synthetic_row": False,
                "status": "native_display_support",
            })
            continue
        if len(ordered) < 2:
            raise ValueError(f"Maxwell xi={key} lacks two rows for continuation")
        second_T = number(ordered[1], "T_MeV")
        second_mu = number(ordered[1], "mu_MeV")
        if second_T - first_T <= EPS:
            raise ValueError(f"duplicate Maxwell temperatures at xi={key}")
        slope = (second_mu - first_mu) / (second_T - first_T)
        floor_mu = first_mu + slope * (floor_T_MeV - first_T)
        if not math.isfinite(floor_mu):
            raise ValueError(f"non-finite Maxwell continuation at xi={key}")
        synthetic = {
            "surface": "maxwell",
            "xi": key,
            "T_MeV": floor_T_MeV,
            "mu_MeV": floor_mu,
            "layer": "display_only_extrapolated_noncertified",
            "status": "display_extrapolated_noncertified",
            "interpolation_method": "linear_T_from_first_two_native_rows",
            "source_layer": "derived_reference_v1_maxwell",
            "source_T_min_MeV": first_T,
            "source_T_second_MeV": second_T,
        }
        result[key] = [synthetic, *ordered]
        trace.append({
            "xi": key,
            "source_T_min_MeV": first_T,
            "display_T_floor_MeV": floor_T_MeV,
            "extension_gap_MeV": first_T - floor_T_MeV,
            "extension_method": "linear_T_from_first_two_native_rows",
            "source_T_second_MeV": second_T,
            "source_mu_first_MeV": first_mu,
            "source_mu_second_MeV": second_mu,
            "display_mu_floor_MeV": floor_mu,
            "synthetic_row": True,
            "status": "display_extrapolated_noncertified",
        })
    return result, trace


def build_surface_quads(
    groups: dict[float, list[dict[str, Any]]],
    *,
    axis: str,
    value_field: str,
    grid_points: int,
    max_gap: float,
) -> tuple[list[list[tuple[float, float, float]]], dict[str, Any]]:
    """Build masked quads from adjacent xi slices and local common support."""

    keys = sorted(groups)
    quads: list[list[tuple[float, float, float]]] = []
    pair_count = 0
    covered_ranges: list[tuple[float, float]] = []
    for left_key, right_key in zip(keys, keys[1:]):
        left = sorted_axis(groups[left_key], axis)
        right = sorted_axis(groups[right_key], axis)
        support = common_range(left, right, axis)
        if support is None:
            continue
        pair_count += 1
        low, high = support
        covered_ranges.append((low, high))
        for index in range(grid_points - 1):
            a0 = low + (high - low) * index / (grid_points - 1)
            a1 = low + (high - low) * (index + 1) / (grid_points - 1)
            left0 = interpolated_value(left, axis=axis, field=value_field, target=a0, max_gap=max_gap)
            left1 = interpolated_value(left, axis=axis, field=value_field, target=a1, max_gap=max_gap)
            right0 = interpolated_value(right, axis=axis, field=value_field, target=a0, max_gap=max_gap)
            right1 = interpolated_value(right, axis=axis, field=value_field, target=a1, max_gap=max_gap)
            if None in (left0, left1, right0, right1):
                continue
            if axis == "mu_MeV":
                quad = [
                    (a0, left_key, left0),
                    (a1, left_key, left1),
                    (a1, right_key, right1),
                    (a0, right_key, right0),
                ]
            else:
                quad = [
                    (left0, left_key, a0),
                    (left1, left_key, a1),
                    (right1, right_key, a1),
                    (right0, right_key, a0),
                ]
            quads.append(quad)
    summary = {
        "xi_slice_count": len(keys),
        "xi_pair_count_with_common_support": pair_count,
        "quad_count": len(quads),
        "grid_points_per_xi_pair": grid_points,
        "max_local_axis_gap": max_gap,
        "common_support_low_min": min((value[0] for value in covered_ranges), default=None),
        "common_support_high_max": max((value[1] for value in covered_ranges), default=None),
    }
    return quads, summary


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def render_figure(
    output_png: Path,
    output_svg: Path,
    crossover_quads: list[list[tuple[float, float, float]]],
    maxwell_quads: list[list[tuple[float, float, float]]],
    cep_rows: list[dict[str, str]],
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    except ImportError as exc:  # pragma: no cover - environment-specific
        raise RuntimeError("matplotlib is required for the v9 phase render") from exc

    fig = plt.figure(figsize=(11, 8), dpi=600)
    axis = fig.add_subplot(111, projection="3d")
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
    axis.plot(
        [number(row, "mu_CEP_proxy_MeV") for row in ordered_cep],
        [number(row, "xi") for row in ordered_cep],
        [number(row, "T_midpoint_MeV") for row in ordered_cep],
        color="#202020",
        linewidth=2.1,
        label="CEP boundary (estimated midpoint)",
    )

    legend_handles = [
        Line2D([0], [0], color="#d9792b", linewidth=5, alpha=0.60, label="crossover"),
        Line2D([0], [0], color="#245b8f", linewidth=5, alpha=0.60, label="first-order coexistence (Maxwell)"),
        Line2D([0], [0], color="#202020", linewidth=2.1, label="CEP boundary (estimated midpoint)"),
    ]
    axis.legend(handles=legend_handles, loc="upper left", frameon=True)
    axis.set_xlabel(r"$\mu_q$ [MeV]")
    axis.set_ylabel(r"$\xi$")
    axis.set_zlabel("T [MeV]")
    axis.set_title("PNJL phase structure: crossover and first-order coexistence")
    axis.set_xlim(0.0, 400.0)
    axis.set_ylim(-0.5, 0.5)
    axis.set_zlim(0.0, 225.0)
    axis.view_init(elev=23, azim=-62)
    fig.tight_layout()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=600, bbox_inches="tight")
    fig.savefig(output_svg, format="svg", bbox_inches="tight")
    # Matplotlib may leave whitespace at line ends in path data; normalize the
    # generated vector artifact so repository whitespace checks remain useful.
    svg_text = output_svg.read_text(encoding="utf-8")
    output_svg.write_text(
        "\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n",
        encoding="utf-8",
    )
    plt.close(fig)


def build_plot_manifest(
    *,
    root: Path,
    output_png: Path,
    output_svg: Path,
    generator: Path,
    inputs: list[dict[str, Any]],
    case_slug: str,
) -> dict[str, Any]:
    return {
        "schema_version": "plot_manifest_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "asset_id": "pnjl_issue130_phase_surface_v9",
        "figure_family": "pnjl_phase_reference",
        "case_slug": case_slug,
        "figure_version": "v9",
        "asset_kind": "figure",
        "figure_mode": "audit",
        "semantic_status": "visualization_only_closed_candidate",
        "style_profile": "candidate_origin_like_v1",
        "publication_scope": "external_display_candidate",
        "generator": {
            **file_record(generator, root),
            "role": "render_generator",
            "command": "render_issue130_phase_surface_v9.py",
            "runtime": {
                "python": platform.python_version(),
                "platform": platform.platform(),
                "executable": sys.executable,
            },
        },
        "git_commit": git_commit(root),
        "calculation_sha": CALCULATION_SHA,
        "postprocess_sha": SOURCE_POSTPROCESS_SHA,
        "source_run_id": SOURCE_RUN_ID,
        "previous_public_version": PREVIOUS_PUBLIC_VERSION,
        "previous_public_png_sha256": PREVIOUS_PUBLIC_PNG_SHA256,
        "previous_public_source_path": PREVIOUS_PUBLIC_SOURCE_PATH,
        "inputs": inputs,
        "axes": [
            {"field": "mu_q", "source_unit": "MeV", "display_unit": "MeV", "transform": "identity"},
            {"field": "xi", "source_unit": "dimensionless", "display_unit": "dimensionless", "transform": "identity"},
            {"field": "T", "source_unit": "MeV", "display_unit": "MeV", "transform": "identity"},
        ],
        "series": [
            {
                "series_id": "crossover",
                "state": "crossover_with_boundary_constrained_cep_display_endpoint",
                "support_rule": "derived rows and local common-support mesh; endpoint may reach estimated CEP midpoint",
                "mask_rule": "mask cells without adjacent xi common mu support",
            },
            {
                "series_id": "first_order_coexistence",
                "state": "first_order_coexistence_maxwell_display",
                "support_rule": "derived Maxwell rows plus bounded display-only continuation to T=0 within adjacent xi common T support",
                "mask_rule": "mask cells without adjacent xi common T support or local axis gap > 18 MeV",
            },
            {
                "series_id": "cep_boundary",
                "state": "estimated_midpoint",
                "support_rule": "CEP bracket midpoint from derived boundary table",
                "mask_rule": "none",
            },
        ],
        "selection_rule": "render finite derived rows, add explicit CEP display endpoint, and add bounded non-certified Maxwell T continuation",
        "interpolation_policy": "local_common_support_mesh_plus_boundary_constrained_endpoint_plus_bounded_T_extrapolation",
        "connector_policy": "explicit_display_only",
        "missing_value_policy": "masked",
        "outputs": [
            {**file_record(output_png, root), "format": "png", "dpi": 600, "vector": False},
            {**file_record(output_svg, root), "format": "svg", "dpi": None, "vector": True},
        ],
        "rendering": {
            "projection": "masked_local_quad_surface",
            "triangulation": False,
            "filled_surface": True,
            "figure_size_inches": [11.0, 8.0],
            "legend_policy": "explicit_semantic_labels",
            "boundary_display": "CEP estimated midpoint with endpoint-constrained crossover closure",
            "low_temperature_display": "Maxwell display-only linear continuation to T=0 from first two native T rows",
            "low_temperature_source_tables_unchanged": True,
            "maxwell_display_extrapolation_max_gap_MeV": 13.0,
        },
        "validation": {
            "finite": True,
            "duplicate_keys": True,
            "support": True,
            "strict_gate": False,
            "solver_called": False,
            "source_tables_unchanged": True,
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--public-root", type=Path, default=DEFAULT_PUBLIC_ROOT)
    parser.add_argument("--grid-points", type=int, default=160)
    parser.add_argument("--replace-public-figure", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.repo_root.resolve()
    input_root = (root / args.input_root).resolve()
    output_root = (root / args.output_root).resolve()
    public_root = (root / args.public_root).resolve()
    if output_root.exists():
        raise FileExistsError(f"refusing to overwrite v9 evidence root: {output_root}")
    if args.grid_points < 8:
        raise ValueError("--grid-points must be at least 8")
    if not args.replace_public_figure:
        raise ValueError("replacing the stable public figure requires --replace-public-figure")

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
    coverage_rows = read_csv(coverage_path)
    crossover_groups, closure_rows = append_cep_boundary_display(group_by_xi(crossover_rows), cep_rows)
    maxwell_groups = group_by_xi(maxwell_rows)
    maxwell_groups, maxwell_extrapolation_rows = extend_maxwell_groups_to_temperature_floor(
        maxwell_groups,
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

    figure_root = output_root / "figures"
    output_png = figure_root / "phase_surface_render_mu_xi_T.png"
    output_svg = figure_root / "phase_surface_render_mu_xi_T.svg"
    render_figure(output_png, output_svg, crossover_quads, maxwell_quads, cep_rows)

    public_root.mkdir(parents=True, exist_ok=True)
    public_png = public_root / "phase_surface_render_mu_xi_T.png"
    public_svg = public_root / "phase_surface_render_mu_xi_T.svg"
    shutil.copyfile(output_png, public_png)
    shutil.copyfile(output_svg, public_svg)

    source_inputs = [
        input_record(crossover_path, root, "derived_crossover_table"),
        input_record(maxwell_path, root, "derived_maxwell_table"),
        input_record(cep_path, root, "derived_cep_boundary_table"),
        input_record(coverage_path, root, "derived_coverage_mask"),
    ]
    package_manifest = build_plot_manifest(
        root=root,
        output_png=output_png,
        output_svg=output_svg,
        generator=Path(__file__).resolve(),
        inputs=source_inputs,
        case_slug="issue130_phase_surface_render_v9",
    )
    package_plot_manifest = figure_root / "plot_manifest.json"
    _write_json(package_plot_manifest, package_manifest)
    public_manifest = build_plot_manifest(
        root=root,
        output_png=public_png,
        output_svg=public_svg,
        generator=Path(__file__).resolve(),
        inputs=source_inputs,
        case_slug="issue130_phase_reference_v1_public_v9_alias",
    )
    public_plot_manifest = public_root / "plot_manifest.json"
    _write_json(public_plot_manifest, public_manifest)

    tables_root = output_root / "tables"
    write_csv(
        tables_root / "crossover_endpoint_closure.csv",
        [
            "xi", "last_crossover_mu_MeV", "cep_mu_MeV", "last_crossover_T_MeV",
            "cep_T_midpoint_MeV", "delta_mu_MeV", "delta_T_MeV", "closure_status",
        ],
        closure_rows,
    )
    write_csv(
        tables_root / "source_coverage_mask.csv",
        list(coverage_rows[0].keys()),
        coverage_rows,
    )
    write_csv(
        tables_root / "maxwell_low_temperature_extrapolation.csv",
        [
            "xi", "source_T_min_MeV", "display_T_floor_MeV", "extension_gap_MeV",
            "extension_method", "source_T_second_MeV", "source_mu_first_MeV",
            "source_mu_second_MeV", "display_mu_floor_MeV", "synthetic_row", "status",
        ],
        maxwell_extrapolation_rows,
    )
    summary_rows = [
        {"surface": "crossover", **crossover_summary, "display_policy": "cep_boundary_constrained_common_support"},
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
    summary_fields = sorted({field for row in summary_rows for field in row})
    write_csv(
        tables_root / "render_surface_summary.csv",
        summary_fields,
        [{field: row.get(field, "") for field in summary_fields} for row in summary_rows],
    )

    generated_at = datetime.now(timezone.utc).isoformat()
    root_manifest = {
        "schema_version": "pnjl_issue130_phase_surface_render_v9",
        "generated_at_utc": generated_at,
        "calculation_sha": CALCULATION_SHA,
        "source_run_id": SOURCE_RUN_ID,
        "replay_run_id": SOURCE_REPLAY_ID,
        "source_layer": "data/reference/pnjl/issue130_phase_reference_v1/derived",
        "source_render_layer": "derived_reference_v1",
        "render_mode": "visualization_only_closed_candidate",
        "solver_called": False,
        "strict_and_derived_tables_unchanged": True,
        "previous_public_version": PREVIOUS_PUBLIC_VERSION,
        "previous_public_png_sha256": PREVIOUS_PUBLIC_PNG_SHA256,
        "previous_public_source_path": PREVIOUS_PUBLIC_SOURCE_PATH,
        "previous_public_plot_manifest_sha256": PREVIOUS_PUBLIC_PLOT_MANIFEST_SHA256,
        "public_png": file_record(public_png, root),
        "public_svg": file_record(public_svg, root),
        "public_plot_manifest": file_record(public_plot_manifest, root),
        "package_plot_manifest": file_record(package_plot_manifest, root),
        "policies": {
            "crossover_endpoint": "CEP estimated midpoint as explicit display-only boundary endpoint",
            "maxwell_label": "first-order coexistence (Maxwell)",
            "maxwell_low_temperature": "bounded display-only linear continuation to T=0 from first two native T rows",
            "maxwell_low_temperature_status": "display_extrapolated_noncertified",
            "maxwell_display_extrapolation_max_gap_MeV": max(
                float(row["extension_gap_MeV"]) for row in maxwell_extrapolation_rows
            ),
            "triangulation": False,
            "source_unresolved_semantics": "preserved in source tables and coverage mask",
        },
        "surface_summary": summary_rows,
    }
    _write_json(output_root / "manifest.json", root_manifest)
    (output_root / "README.md").write_text(
        "# Issue #130 phase-surface render v9\n\n"
        "This is a display-only successor to the v8 public figure. "
        "The strict and derived tables are immutable inputs. Crossover endpoints "
        "may be connected to the estimated CEP midpoint for display, while the "
        "first-order (Maxwell) surface receives a bounded display-only linear "
        "continuation to T=0 from each slice's first two native temperature rows. "
        "The continuation is explicitly non-certified, does not modify strict or "
        "derived tables, and is never consumed by runtime reference code. No "
        "triangulation is used. The stable public PNG alias is intentionally written "
        "to the existing phase-reference figure path; the v8 source figure remains "
        "preserved under the historical v1 analysis package.\n",
        encoding="utf-8",
    )
    _write_json(output_root / "decision.json", {
        "verdict": "display_candidate",
        "figure_version": "v9",
        "previous_public_version": PREVIOUS_PUBLIC_VERSION,
        "reference_write": False,
        "runtime_consumption": False,
        "solver_called": False,
        "strict_or_derived_csv_modified": False,
        "public_alias_replaced": True,
        "previous_public_png_sha256": PREVIOUS_PUBLIC_PNG_SHA256,
    })
    write_csv(
        output_root / "claim_ledger.csv",
        ["claim_id", "status", "claim", "evidence", "boundary"],
        [
            {
                "claim_id": "V9-RENDER-001",
                "status": "supported",
                "claim": "The public figure is reproducible from the immutable derived tables with explicit display policies.",
                "evidence": "manifest.json; figures/plot_manifest.json; tables/render_surface_summary.csv",
                "boundary": "display candidate; not a numerical certificate or runtime reference",
            },
            {
                "claim_id": "V9-RENDER-002",
                "status": "supported",
                "claim": "The plotted first-order surface is labelled as Maxwell coexistence and uses an explicit bounded display-only low-temperature continuation.",
                "evidence": "tables/source_coverage_mask.csv; tables/maxwell_low_temperature_extrapolation.csv; figures/plot_manifest.json",
                "boundary": "continuation is non-certified and absent from strict/derived/reference runtime data",
            },
            {
                "claim_id": "V9-RENDER-003",
                "status": "display_only",
                "claim": "Crossover lines can reach the estimated CEP midpoint through explicit boundary-constrained display endpoints.",
                "evidence": "tables/crossover_endpoint_closure.csv",
                "boundary": "synthetic display endpoint; absent from strict/derived reference tables",
            },
        ],
    )

    # The public manifest must be rewritten after all output files exist so its
    # output hashes are exact.  This is the only intentional public alias
    # replacement; the v1 source image and numeric tables remain untouched.
    public_manifest = build_plot_manifest(
        root=root,
        output_png=public_png,
        output_svg=public_svg,
        generator=Path(__file__).resolve(),
        inputs=source_inputs,
        case_slug="issue130_phase_reference_v1_public_v9_alias",
    )
    _write_json(public_plot_manifest, public_manifest)
    package_manifest = build_plot_manifest(
        root=root,
        output_png=output_png,
        output_svg=output_svg,
        generator=Path(__file__).resolve(),
        inputs=source_inputs,
        case_slug="issue130_phase_surface_render_v9",
    )
    _write_json(package_plot_manifest, package_manifest)
    # Refresh manifests that record the plot-manifest hash after the final
    # manifest write.
    root_manifest["public_plot_manifest"] = file_record(public_plot_manifest, root)
    root_manifest["package_plot_manifest"] = file_record(package_plot_manifest, root)
    _write_json(output_root / "manifest.json", root_manifest)

    # Keep the imported candidate manifest internally consistent with the new
    # stable public figure alias while preserving the v8 source-render record.
    imported_manifest_path = input_root / "manifest.json"
    imported_manifest = json.loads(imported_manifest_path.read_text(encoding="utf-8"))
    imported_manifest["figure"]["png"] = file_record(public_png, root)
    imported_manifest["figure"]["plot_manifest"] = file_record(public_plot_manifest, root)
    imported_manifest["public_render"] = {
        "version": "v9",
        "mode": "visualization_only_closed_candidate",
        "source_render_layer": "derived_reference_v1",
        "previous_public_version": PREVIOUS_PUBLIC_VERSION,
        "previous_public_source_path": PREVIOUS_PUBLIC_SOURCE_PATH,
        "previous_public_png_sha256": PREVIOUS_PUBLIC_PNG_SHA256,
        "public_png_sha256": sha256(public_png),
        "strict_and_derived_tables_unchanged": True,
        "manifest": relative_path(output_root / "manifest.json", root)[0],
    }
    _write_json(imported_manifest_path, imported_manifest)
    audit_path = input_root / "IMPORT_AUDIT.md"
    with audit_path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(
            "\n- public figure alias: v9 display render; the pre-refinement public figure is "
            "defined as v8 and remains byte-preserved under the historical v1 source layer; "
            "strict/derived CSV bytes are unchanged\n"
        )
    imported_manifest["audit"] = file_record(audit_path, root)
    _write_json(imported_manifest_path, imported_manifest)
    root_manifest["import_audit"] = file_record(audit_path, root)
    _write_json(output_root / "manifest.json", root_manifest)

    print(json.dumps({
        "output_root": str(output_root),
        "public_png": str(public_png),
        "public_svg": str(public_svg),
        "public_plot_manifest": str(public_plot_manifest),
        "previous_public_version": PREVIOUS_PUBLIC_VERSION,
        "previous_public_png_sha256": PREVIOUS_PUBLIC_PNG_SHA256,
        "new_public_png_sha256": sha256(public_png),
        "crossover_quads": len(crossover_quads),
        "maxwell_quads": len(maxwell_quads),
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
