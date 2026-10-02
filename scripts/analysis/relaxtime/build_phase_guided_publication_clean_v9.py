#!/usr/bin/env python3
"""Build the PNG-only ``publication_clean_v9`` figure review layer.

v9 is a solver-free display derivative of the accepted v5 point table.  It
keeps the v5 values, phase gaps, adjustment records, and raw provenance
unchanged.  The review layer changes only the figure layout and style:

* one reviewed in-axes legend per single chart;
* one compact shared legend in a reviewed sparse host panel of each native
  composite, including the first-order endpoint marker key;
* three labeled x-axis major ticks and about five linear y-axis major ticks;
* one linear minor tick between adjacent linear major ticks;
* a small y margin so curves do not touch the frame;
* a local log-y scale for positive high-range panels, with the complete
  ``mu_B=900 MeV`` column of the four-row relaxation-time composite forced to
  log-y for column consistency.

Only 600 dpi PNG files are emitted in this stage.  Vector delivery is a later
author-approved stage and is deliberately not generated here.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import datetime as dt
import importlib.util
import json
import math
from pathlib import Path
import shutil
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, LogLocator, MaxNLocator

from scripts.plotting.plot_manifest import (
    build_manifest,
    generator_record,
    input_record,
    runtime_record,
    sha256_file,
    write_manifest,
)
from scripts.plotting.plot_quality import export_figure
from scripts.plotting.plot_style import configure_axis_ticks, configure_matplotlib, load_profile
from scripts.plotting.validate_plot_artifact import validate_manifest


V8_SCRIPT = ROOT / "scripts/analysis/relaxtime/build_phase_guided_publication_clean_v8.py"
V8_SPEC = importlib.util.spec_from_file_location("publication_clean_v8_inputs_for_v9", V8_SCRIPT)
if V8_SPEC is None or V8_SPEC.loader is None:
    raise RuntimeError(f"cannot load v8 input reader: {V8_SCRIPT}")
V8 = importlib.util.module_from_spec(V8_SPEC)
V8_SPEC.loader.exec_module(V8)
V7 = V8.V7

PARENT = V8.PARENT
V2 = PARENT.V2
TRANSPORT_ANALYSIS_ROOT = PARENT.TRANSPORT_ANALYSIS_ROOT
V5_ROOT = PARENT.V5_ROOT
V5_TABLE_ROOT = PARENT.V5_TABLE_ROOT
V5_MANIFEST = PARENT.V5_MANIFEST
V5_PLOT_MANIFEST = PARENT.V5_PLOT_MANIFEST
V5_FIGURE_LAYER_MANIFEST = PARENT.V5_FIGURE_LAYER_MANIFEST
CURRENT_POINTER = PARENT.CURRENT_POINTER

ANALYSIS_ROOT = TRANSPORT_ANALYSIS_ROOT / "phase_guided_transport_publication_clean_v9_png_review"
FIGURE_ROOT = ROOT / "data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v9_png_review"
TABLE_ROOT = ANALYSIS_ROOT / "tables"
STYLE_MANIFEST = ANALYSIS_ROOT / "v9_display_style_manifest.json"
CAPTION_HANDOFF = ANALYSIS_ROOT / "caption_handoff.md"

OBSERVABLES = list(V8.OBSERVABLES)
UNITS = dict(V8.UNITS)
LABELS = dict(V8.LABELS)
COMPOSITES = dict(V8.COMPOSITES)
PANELS = list(V8.PANELS)

ENDPOINT_LABELS = dict(V8.ENDPOINT_LABELS)
ENDPOINT_MARKERS = dict(V8.ENDPOINT_MARKERS)
X_MAJOR_TICKS = (-0.5, 0.0, 0.5)
X_MAJOR_TICK_LABELS = {
    -0.5: "-0.5",
    0.0: "0",
    0.5: "0.5",
}
LINEAR_MINOR_SUBDIVISIONS = 2
LOG_Y_DYNAMIC_RANGE_THRESHOLD = 20.0
Y_MARGIN = 0.06
FORCED_LOG_COMPOSITE = "figure1_relaxation_times_comparison"
FORCED_LOG_PANEL = "muB900.0"
COMPOSITE_LEGEND_HOSTS = {
    "figure1_relaxation_times_comparison": "row0_col1",
    "figure2_transport_coefficients_comparison": "row1_col1",
}


def series_label(mode: str, row: dict[str, str]) -> str:
    if mode == "mode_a":
        alpha = float(row["plot_series"].removeprefix("alpha"))
        return rf"$\alpha_T = {alpha:.1f}$"
    return rf"$\mu_B = {float(row['muB_MeV']):.0f}\;\mathrm{{MeV}}$"


def panel_title(mode: str, panel: str) -> str:
    if mode == "mode_a":
        return rf"$\mu_B = {float(panel.removeprefix('muB')):.0f}\;\mathrm{{MeV}}$"
    return rf"$T = {float(panel.removeprefix('T')):.0f}\;\mathrm{{MeV}}$"


def format_x_tick(value: float, _position: int) -> str:
    for tick, label in X_MAJOR_TICK_LABELS.items():
        if math.isclose(value, tick, rel_tol=0.0, abs_tol=1.0e-9):
            return label
    return ""


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def group_inputs(points: list[dict[str, str]], gaps: list[dict[str, Any]]) -> tuple[dict, dict]:
    grouped = defaultdict(list)
    gap_map = defaultdict(list)
    for row in points:
        grouped[(row["mode_key"], row["plot_panel"], row["observable"], row["plot_series"])].append(row)
    for rows in grouped.values():
        rows.sort(key=lambda row: float(row["xi"]))
    for gap in gaps:
        gap_map[(gap["mode_key"], gap["plot_panel"], gap["plot_series"])].append(gap)
    return grouped, gap_map


def axis_spec_for_rows(
    rows: list[dict[str, Any]],
    gaps: list[dict[str, Any]],
    *,
    force_log: bool = False,
) -> dict[str, Any]:
    """Choose the v9 display transform without changing any point values."""

    values = [float(row["clean_value"]) for row in rows]
    if not values or not all(math.isfinite(value) for value in values):
        raise ValueError("cannot choose an axis scale for empty or non-finite data")
    data_min = min(values)
    data_max = max(values)
    ratio = data_max / data_min if data_min != 0.0 else math.inf
    parent = V2.figure_axis_spec(rows, gaps)
    common = {
        "data_min": data_min,
        "data_max": data_max,
        "dynamic_range_ratio": ratio,
        "first_order_gap_present": bool(gaps),
        "threshold": LOG_Y_DYNAMIC_RANGE_THRESHOLD,
        "parent_axis_scale": parent["axis_scale"],
        "parent_axis_scale_reason": parent["axis_scale_reason"],
    }
    if force_log and data_min <= 0.0:
        raise ValueError("forced log-y panel contains a non-positive display value")
    if force_log:
        return {
            **common,
            "axis_scale": "log",
            "axis_scale_reason": "forced_log_for_complete_muB900_relaxation_column",
        }
    if data_min <= 0.0:
        return {**common, "axis_scale": "linear", "axis_scale_reason": "non_positive_display_value"}
    if ratio < LOG_Y_DYNAMIC_RANGE_THRESHOLD:
        return {
            **common,
            "axis_scale": "linear",
            "axis_scale_reason": "positive_range_ratio_below_v9_local_log_threshold",
        }
    return {
        **common,
        "axis_scale": "log",
        "axis_scale_reason": "positive_range_ratio_ge_20_for_v9_local_log_readability",
    }


def render_panel(
    ax: Any,
    mode: str,
    panel: str,
    observable: str,
    grouped: dict,
    gap_map: dict,
    profile: Any,
    *,
    force_log: bool = False,
) -> dict[str, Any]:
    names = sorted(key[3] for key in grouped if key[:3] == (mode, panel, observable))
    rows_all: list[dict[str, str]] = []
    gaps_all: list[dict[str, Any]] = []
    series_records: list[dict[str, Any]] = []
    endpoint_records: list[dict[str, Any]] = []
    for index, name in enumerate(names):
        rows = grouped[(mode, panel, observable, name)]
        gaps = gap_map.get((mode, panel, name), [])
        segments = V2.split_curve_segments(rows, gaps)
        color = profile.colors[index]
        linestyle = profile.data["palette"]["parameter_linestyles"][index]
        for segment in segments:
            ax.plot(
                [float(row["xi"]) for row in segment],
                [float(row["clean_value"]) for row in segment],
                color=color,
                linestyle=linestyle,
            )
        point_index = {V2.canonical_xi(row["xi"]): row for row in rows}
        for gap in gaps:
            for role in ("left", "right"):
                endpoint = gap["endpoint_rows"][role]
                phase = endpoint["phase_curr"]
                if phase not in ENDPOINT_MARKERS:
                    raise ValueError(f"unresolved first-order endpoint: {gap['boundary_id']}")
                point = point_index[V2.canonical_xi(endpoint["xi"])]
                ax.scatter(
                    [float(point["xi"])],
                    [float(point["clean_value"])],
                    marker=ENDPOINT_MARKERS[phase],
                    s=profile.data["marker_size_pt"] ** 2,
                    facecolor="white",
                    edgecolor=color,
                    linewidth=0.8,
                    zorder=6,
                )
                endpoint_records.append(
                    {
                        "boundary_id": gap["boundary_id"],
                        "role": role,
                        "phase": phase,
                        "xi": point["xi"],
                        "clean_value": point["clean_value"],
                        "color": color,
                        "series": name,
                    }
                )
        series_records.append(
            {
                "series_id": f"{mode}.{panel}.{observable}.{name}",
                "state": "author_accepted_display_derivative",
                "label": series_label(mode, rows[0]),
                "support_rule": "v5 frozen clean_value rows; no new values or support points",
                "mask_rule": "preserve audited v2 first-order gap; split each gap into separate segments",
                "row_count": len(rows),
                "segment_count": len(segments),
                "color": color,
                "linestyle": linestyle,
            }
        )
        rows_all.extend(rows)
        gaps_all.extend(gaps)

    axis_spec = axis_spec_for_rows(rows_all, gaps_all, force_log=force_log)
    ax.set_yscale(axis_spec["axis_scale"])
    ax.set_xlim(-0.52, 0.52)
    ax.set_xticks(X_MAJOR_TICKS)
    ax.xaxis.set_major_formatter(FuncFormatter(format_x_tick))
    if axis_spec["axis_scale"] == "linear":
        ax.yaxis.set_major_locator(
            MaxNLocator(
                nbins=5,
                min_n_ticks=4,
                steps=[1, 2, 2.5, 5, 10],
            )
        )
    else:
        ax.yaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0,)))
    ax.margins(y=Y_MARGIN)
    configure_axis_ticks(ax, profile)
    if observable == "sigma_over_T":
        from matplotlib.ticker import FormatStrFormatter

        ax.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))
    # Three explicit x labels keep the repeated bottom row legible in the
    # thesis-like composite layout.
    ax.tick_params(axis="x", which="major", pad=6)
    ax.tick_params(axis="y", which="major", pad=3)
    ax.set_ylabel(LABELS[observable])
    ax.yaxis.labelpad = 3
    return {
        "mode_key": mode,
        "plot_panel": panel,
        "observable": observable,
        "series": series_records,
        "endpoints": endpoint_records,
        **axis_spec,
    }


def endpoint_colors(specs: list[dict[str, Any]]) -> list[str]:
    colors = []
    for spec in specs:
        colors.extend(record["color"] for record in spec["endpoints"])
    return sorted(set(colors))


def legend_handles(
    specs: list[dict[str, Any]],
    profile: Any,
    *,
    include_endpoints: bool,
    line_break_endpoint_labels: bool = False,
) -> list[Any]:
    first = specs[0]
    handles = [
        Line2D(
            [],
            [],
            color=item["color"],
            linestyle=item["linestyle"],
            linewidth=profile.data["line_width_pt"],
            label=item["label"],
        )
        for item in first["series"]
    ]
    if include_endpoints:
        colors = endpoint_colors(specs)
        marker_color = colors[0] if len(colors) == 1 else "#000000"
        for phase, label in ENDPOINT_LABELS.items():
            if line_break_endpoint_labels:
                label = label.replace(" (", "\n(", 1)
            handles.append(
                Line2D(
                    [],
                    [],
                    linestyle="None",
                    marker=ENDPOINT_MARKERS[phase],
                    markerfacecolor="white",
                    markeredgecolor=marker_color,
                    markersize=profile.data["marker_size_pt"],
                    label=label,
                )
            )
    return handles


def apply_in_axes_legend(ax: Any, handles: list[Any], profile: Any, *, composite: bool = False) -> Any:
    ncol = 1 if composite else (2 if len(handles) >= 4 else max(1, len(handles)))
    legend_kwargs = {
        "loc": "upper left" if composite else "best",
        "ncol": ncol,
    }
    if composite:
        legend_kwargs["bbox_to_anchor"] = (0.02, 0.98)
    legend = ax.legend(
        handles=handles,
        **legend_kwargs,
        fontsize=profile.data["legend_font_size_pt"],
        frameon=False,
        handlelength=1.25 if composite else 1.8,
        columnspacing=0.5,
        handletextpad=0.35,
        labelspacing=0.12 if composite else 0.25,
        borderpad=0.1 if composite else 0.25,
    )
    legend.set_zorder(10)
    return legend


def write_chart(
    figure: Any,
    stem: Path,
    specs: list[dict[str, Any]],
    profile: Any,
    inputs: list[dict],
    font: dict,
    command: str,
    *,
    composite: bool,
) -> tuple[dict[str, Any], Path]:
    outputs, quality = export_figure(figure, stem, profile, formats=("png",))
    axes_records = []
    for spec in specs:
        axes_records.extend(
            [
                {
                    "field": "xi",
                    "source_unit": "dimensionless",
                    "display_unit": "dimensionless",
                    "label": r"$\xi$",
                    "transform": "identity",
                    "panel": spec["plot_panel"],
                },
                {
                    "field": spec["observable"],
                    "source_unit": UNITS[spec["observable"]],
                    "display_unit": UNITS[spec["observable"]],
                    "label": LABELS[spec["observable"]],
                    "transform": spec["axis_scale"],
                    "panel": spec["plot_panel"],
                },
            ]
        )
    legend_policy = "best_in_axes_reviewed"
    legend_host = COMPOSITE_LEGEND_HOSTS.get(stem.name) if composite else None
    manifest = build_manifest(
        asset_id=f"relaxtime.phase_guided.v9.png_review.{relative(stem)}",
        figure_family="phase_guided_transport",
        case_slug=stem.name,
        figure_mode="audit",
        semantic_status="author_review_display_derivative",
        style_profile=profile.profile_id,
        publication_scope="internal_review",
        generator=generator_record(
            Path(__file__),
            command=command,
            runtime=runtime_record({"matplotlib": matplotlib.__version__, "font": font}),
        ),
        inputs=inputs,
        axes=axes_records,
        series=[item for spec in specs for item in spec["series"]],
        outputs=outputs,
        selection_rule="render frozen v5 clean_value table in original xi order; rows/channels/phase gaps unchanged",
        interpolation_policy="inherited v5 display adjustments; no new interpolation or replacement",
        connector_policy="forbidden",
        missing_value_policy="preserve first-order and missing-support gaps",
        validation={"finite": True, "duplicate_keys": True, "support": True, "strict_gate": False},
        rendering={
            "column": "double_column",
            "figure_size_inches": quality["figure_size_inches"],
            "size_override_reason": "native manuscript composite at final physical size" if composite else None,
            "bbox_inches": None,
            "quality": quality,
            "legend_policy": legend_policy,
            "legend_location": f"{legend_host}_upper_left" if composite else "best",
            "legend_host_panel": legend_host,
            "legend_outside": False,
            "legend_contents": "parameter curves plus restored/broken first-order endpoint marker key",
            "legend_marker_color_policy": "endpoint marker colors match the data series; the compact key labels open circles as restored and open squares as broken",
            "case_layout_contract": "thesis_like single in-axes legend in a case-selected sparse host panel; curve visibility reviewed",
            "color_route": "undecided_review",
            "grayscale_review": "required before author acceptance",
            "y_axis_policy": "Figure 1 muB900 complete column forced log-y; other positive high-range panels use local log-y when dynamic range ratio >= 20; linear otherwise",
            "tick_policy": {
                "x_major_ticks": list(X_MAJOR_TICKS),
                "linear_y_major_locator": "MaxNLocator(nbins=5, min_n_ticks=4, steps=[1,2,2.5,5,10])",
                "linear_minor_subdivisions": LINEAR_MINOR_SUBDIVISIONS,
                "log_major_locator": "LogLocator(base=10, subs=[1])",
                "log_minor_locator": "profile log locator subs=[2,3,4,5,6,7,8,9]",
                "sides": ["top", "bottom", "left", "right"],
                "direction": "in",
            },
            "panel_specs": specs,
            "x_tick_label_policy": "three labels at -0.5, 0, and 0.5",
            "sigma_over_T_tick_policy": "linear y labels formatted to three decimal places; leading zero retained",
            "delivery_stage": "png_review",
            "vector_delivery_pending": True,
            "output_formats": ["png"],
        },
        calculation_sha=PARENT.V5.V4.V3.V2.V1.CALCULATION_SHA,
    )
    manifest.update(
        {
            "manuscript_eligible": False,
            "current_publication_layer": False,
            "solver_called": False,
            "canonical_data_modified": False,
            "new_display_values": False,
            "delivery_stage": "png_review",
            "vector_delivery_pending": True,
            "workflow_head_sha": PARENT.V5.V4.V3.V2.V1.WORKFLOW_HEAD_SHA,
            "source_run_id_note": "multiple upstream runs; see frozen v5 per-point run_id and raw manifests",
            "numerical_status": "inherited_author_accepted_display_only",
            "raw_manuscript_eligible": False,
        }
    )
    manifest_path = stem.with_suffix(".plot_manifest.json")
    write_manifest(manifest_path, manifest)
    violations = validate_manifest(manifest_path)
    if violations:
        raise ValueError(f"{stem.name}: " + "; ".join(violations))
    plt.close(figure)
    chart = {
        "stem": relative(stem),
        "manifest": relative(manifest_path),
        "manifest_sha256": sha256_file(manifest_path),
        "kind": "composite" if composite else "single",
        "mode_key": specs[0]["mode_key"],
        "legend_policy": legend_policy,
        "axis_scales": sorted(set(spec["axis_scale"] for spec in specs)),
        "outputs": outputs,
        "minimum_capital_numeral_height_mm": quality["minimum_capital_numeral_height_mm"],
    }
    return chart, manifest_path


def render_single(mode: str, panel: str, observable: str, grouped: dict, gap_map: dict, profile: Any) -> tuple[Any, list[dict]]:
    figure, ax = plt.subplots(figsize=(6.75, 4.6))
    figure.subplots_adjust(left=0.16, right=0.97, bottom=0.14, top=0.95)
    spec = render_panel(
        ax,
        mode,
        panel,
        observable,
        grouped,
        gap_map,
        profile,
        force_log=(
            mode == "mode_a"
            and panel == FORCED_LOG_PANEL
            and observable in COMPOSITES["figure1_relaxation_times_comparison"]
        ),
    )
    ax.set_xlabel(r"$\xi$")
    figure.suptitle(panel_title(mode, panel), y=0.985)
    apply_in_axes_legend(
        ax,
        legend_handles([spec], profile, include_endpoints=bool(spec["endpoints"])),
        profile,
    )
    return figure, [spec]


def render_composite(observables: list[str], grouped: dict, gap_map: dict, profile: Any) -> tuple[Any, list[dict]]:
    nrows = len(observables)
    figure, axes = plt.subplots(nrows, 3, figsize=(6.75, 7.45 if nrows == 4 else 6.15), sharex=True)
    figure.subplots_adjust(left=0.12, right=0.98, bottom=0.10, top=0.925, hspace=0.14, wspace=0.34)
    specs: list[dict[str, Any]] = []
    for row, observable in enumerate(observables):
        for col, panel in enumerate(V7.PANELS):
            ax = axes[row, col]
            force_log = nrows == 4 and panel == FORCED_LOG_PANEL
            spec = render_panel(
                ax,
                "mode_a",
                panel,
                observable,
                grouped,
                gap_map,
                profile,
                force_log=force_log,
            )
            specs.append(spec)
            if col:
                ax.set_ylabel("")
            if row == 0:
                ax.set_title(panel_title("mode_a", panel), pad=4)
            if row == nrows - 1:
                ax.set_xlabel(r"$\xi$")
            ax.text(
                0.02,
                0.96,
                f"({chr(97 + row * 3 + col)})",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=10.5,
            )
    handles = legend_handles(specs, profile, include_endpoints=any(spec["endpoints"] for spec in specs))
    # Case-level host selection keeps the compact key inside a sparse panel
    # without covering a steep blue branch.
    legend_host = axes[0, 1] if nrows == 4 else axes[1, 1]
    apply_in_axes_legend(legend_host, handles, profile, composite=True)
    return figure, specs


def parent_records(package: dict[str, Any], profile: Any) -> list[dict[str, Any]]:
    records = V7.parent_records(package, profile)
    known = {record["path"] for record in records}
    for path, role in ((V8_SCRIPT, "previous_review_renderer"), (Path(__file__), "v9_generator")):
        record = input_record(path, role=role)
        if record["path"] not in known:
            records.append(record)
    return records


def write_caption_handoff(points: list[dict[str, str]], specs: list[dict[str, Any]]) -> Path:
    temperature: dict[tuple[str, str], float] = {}
    for row in points:
        if row["mode_key"] != "mode_a":
            continue
        key = (row["plot_panel"], row["plot_series"])
        value = float(row["T_MeV"])
        if key in temperature and not math.isclose(temperature[key], value, rel_tol=0.0, abs_tol=1.0e-9):
            raise ValueError(f"temperature varies along xi for {key}")
        temperature[key] = value
    table = "\n".join(
        f"| {panel} | " + " | ".join(f"{temperature[(panel, f'alpha{alpha:.1f}')]:.6f}" for alpha in (1, 1.1, 1.2)) + " |"
        for panel in PANELS
    )
    seen: set[tuple[str, str, str]] = set()
    log_rows = []
    for spec in specs:
        key = (spec["mode_key"], spec["plot_panel"], spec["observable"])
        if key in seen or spec["axis_scale"] != "log":
            continue
        seen.add(key)
        log_rows.append(
            f"- `{spec['mode_key']}` / `{spec['plot_panel']}` / `{spec['observable']}`: "
            f"ratio={spec['dynamic_range_ratio']:.6g}, positive data, local log-y."
        )
    log_text = "\n".join(log_rows) if log_rows else "- None; all panels remain linear."
    CAPTION_HANDOFF.parent.mkdir(parents=True, exist_ok=True)
    CAPTION_HANDOFF.write_text(
        r"""# publication_clean_v9 PNG review caption handoff

v9 changes figure layout and display transforms only. It does not add points,
replace values, reconnect phase branches, smooth data, call a solver, or
change the v5 provenance chain.

Recommended caption additions:

> Solid, dashed, and dash-dotted curves correspond to the three parameter
> values shown in the in-axes legend. For Figure 1, these are alpha_T = 1.0,
> 1.1, and 1.2; for Figure 2, the legend retains the corresponding fixed
> $\mu_B$ values. The temperatures associated with $\alpha_T$ are listed
> below. Open circles and open squares denote the chirally restored and
> chirally broken endpoints of the first-order transition, respectively.
> Disconnected
> segments display the two phase branches separately. Panels with a positive
> displayed range ratio of at least 20 use a local logarithmic y axis. In
> Figure 1, the complete $\mu_B=900\,\mathrm{MeV}$ column uses log-y for
> column consistency; remaining panels use a linear y axis unless listed in
> the manifest. The relaxation times are given in $\mathrm{fm}$.

The in-axes legend is a reviewed layout exception for this dense multi-panel
family. It is not a new repository-wide default: the default remains an
external legend unless a case-level contract records and reviews the exception.

Exact frozen temperature mapping (MeV):

| panel | alpha_T=1.0 | alpha_T=1.1 | alpha_T=1.2 |
| --- | ---: | ---: | ---: |
""" + table + "\n\nPanels using local log-y in v9:\n\n" + log_text + "\n\n"
        "The v9 PNG review layer remains manuscript_eligible=false and does not replace the v5/current pointer.\n",
        encoding="utf-8",
    )
    return CAPTION_HANDOFF


def style_manifest_payload(package: dict[str, Any], profile: Any) -> dict[str, Any]:
    return {
        "schema": "publication_clean_v9_png_review_display_style_v1",
        "status": "author_review_required",
        "manuscript_eligible": False,
        "current_publication_layer": False,
        "delivery_stage": "png_review",
        "vector_delivery_pending": True,
        "solver_called": False,
        "canonical_data_modified": False,
        "source_parent_artifact": relative(V5_ROOT),
        "source_parent_manifest_sha256": sha256_file(V5_MANIFEST),
        "source_parent_plot_manifest_sha256": sha256_file(V5_PLOT_MANIFEST),
        "source_parent_figure_layer_manifest": relative(V5_FIGURE_LAYER_MANIFEST),
        "source_parent_figure_layer_manifest_sha256": sha256_file(V5_FIGURE_LAYER_MANIFEST),
        "source_calculation_sha": package.get("calculation_sha"),
        "source_workflow_head_sha": package.get("workflow_head_sha"),
        "generator": relative(Path(__file__)),
        "generator_sha256": sha256_file(Path(__file__)),
        "legend": {
            "policy": "one reviewed in-axes parameter legend per single chart; one compact shared legend with endpoint marker key in a case-selected sparse host panel per composite",
            "font_size_pt": profile.data["legend_font_size_pt"],
            "endpoint_labels": ENDPOINT_LABELS,
            "endpoint_key_location": "in-axes legend and caption_handoff.md",
            "marker_color_policy": "match endpoint data series; open circle=restored and open square=broken",
        },
        "ticks": {
            "x_major_ticks": list(X_MAJOR_TICKS),
            "linear_y_major_locator": "MaxNLocator(nbins=5, min_n_ticks=4, steps=[1,2,2.5,5,10])",
            "linear_minor_subdivisions": LINEAR_MINOR_SUBDIVISIONS,
            "log_major_locator": "LogLocator(base=10, subs=[1])",
            "log_minor_locator": "LogLocator(base=10, subs=[2,3,4,5,6,7,8,9])",
            "sides": ["top", "bottom", "left", "right"],
            "direction": "in",
        },
        "y_axis": {
            "rule": "Figure 1 muB900 complete column forced log-y; otherwise positive data_max/data_min >= 20 uses local log-y",
            "threshold": LOG_Y_DYNAMIC_RANGE_THRESHOLD,
            "margin": Y_MARGIN,
            "values_unchanged": True,
        },
        "figure_count": 72,
        "mode_counts": {"mode_a": 36, "mode_b": 36},
        "formats": ["png"],
        "composite_layout": {
            "x_tick_label_policy": "three labels at -0.5, 0, and 0.5",
            "legend_host_panel": dict(COMPOSITE_LEGEND_HOSTS),
            "figure1_muB900_log_column": True,
            "sigma_over_T_precision": "three decimal places with leading zero retained",
        },
        "profile_sha256": sha256_file(profile.path),
    }


def package_output_paths() -> list[Path]:
    index = FIGURE_ROOT / "plot_manifest.json"
    return [
        ANALYSIS_ROOT / "README.md",
        CAPTION_HANDOFF,
        STYLE_MANIFEST,
        *sorted(TABLE_ROOT.glob("*.csv")),
        index,
        *sorted(FIGURE_ROOT.rglob("*.png")),
    ]


def main() -> None:
    global ANALYSIS_ROOT, FIGURE_ROOT, TABLE_ROOT, STYLE_MANIFEST, CAPTION_HANDOFF

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--png-review", action="store_true", help="generate the PNG-only author-review stage")
    parser.add_argument("--composites-only", action="store_true", help="write a temporary composite preview only")
    parser.add_argument("--output-root", type=Path, default=FIGURE_ROOT)
    args = parser.parse_args()
    if not args.png_review and not args.composites_only:
        raise ValueError("v9 is intentionally PNG-review-only; pass --png-review")
    if args.composites_only and args.output_root.resolve() == FIGURE_ROOT.resolve():
        raise ValueError("composites-only preview requires a separate --output-root")
    if not args.composites_only and (FIGURE_ROOT.exists() or ANALYSIS_ROOT.exists()):
        raise FileExistsError("refusing to overwrite an existing v9 review package")

    points, package, _, _, gaps, table_hashes = PARENT.load_v5_inputs()
    keys = [(row["mode_key"], row["plot_panel"], row["plot_series"], row["observable"], row["xi"]) for row in points]
    if len(keys) != len(set(keys)):
        raise ValueError("duplicate frozen v5 point keys")
    for item in package["outputs"]:
        path = ROOT / item["path"]
        if sha256_file(path) != item["sha256"]:
            raise ValueError(f"parent package hash mismatch: {path}")

    profile = load_profile("candidate_aps_v2")
    font = configure_matplotlib(profile)
    grouped, gap_map = group_inputs(points, gaps)
    if args.composites_only:
        args.output_root.mkdir(parents=True)
        for name, observables in COMPOSITES.items():
            figure, _ = render_composite(observables, grouped, gap_map, profile)
            figure.savefig(args.output_root / f"{name}.png", dpi=150, bbox_inches=None)
            plt.close(figure)
        return

    ANALYSIS_ROOT.mkdir(parents=True)
    shutil.copytree(V5_TABLE_ROOT, TABLE_ROOT, dirs_exist_ok=False)
    for name, expected in table_hashes.items():
        if sha256_file(TABLE_ROOT / name) != expected:
            raise ValueError(f"v9 changed inherited table bytes: {name}")
    inputs = parent_records(package, profile)
    command = "python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v9.py --png-review"
    charts: list[dict[str, Any]] = []
    all_specs: list[dict[str, Any]] = []
    for mode in ("mode_a", "mode_b"):
        panels = sorted({key[1] for key in grouped if key[0] == mode})
        for panel in panels:
            for observable in OBSERVABLES:
                figure, specs = render_single(mode, panel, observable, grouped, gap_map, profile)
                all_specs.extend(specs)
                mode_dir = "mode_a_fixed_muB_phase_scaled" if mode == "mode_a" else "mode_b_fixed_T_sparse_muB"
                stem = FIGURE_ROOT / mode_dir / f"plot_panel={panel}" / f"{observable}_vs_xi"
                chart, _ = write_chart(figure, stem, specs, profile, inputs, font, command, composite=False)
                charts.append(chart)
        print(f"[v9] {mode}: 36 single charts validated", flush=True)
    for name, observables in COMPOSITES.items():
        figure, specs = render_composite(observables, grouped, gap_map, profile)
        all_specs.extend(specs)
        chart, _ = write_chart(figure, FIGURE_ROOT / "composites" / name, specs, profile, inputs, font, command, composite=True)
        charts.append(chart)

    for item in inputs:
        if sha256_file(ROOT / item["path"]) != item["sha256"]:
            raise ValueError(f"input changed during v9 render: {item['path']}")
    for name, expected in table_hashes.items():
        if sha256_file(TABLE_ROOT / name) != expected:
            raise ValueError(f"v9 changed inherited table bytes after render: {name}")

    caption_path = write_caption_handoff(points, all_specs)
    style = style_manifest_payload(package, profile)
    style["generated_at_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
    write_manifest(STYLE_MANIFEST, style)

    index = {
        "schema": "publication_clean_v9_png_review_figure_index_v1",
        "status": "author_review_required",
        "manuscript_eligible": False,
        "current_publication_layer": False,
        "solver_called": False,
        "canonical_data_modified": False,
        "delivery_stage": "png_review",
        "vector_delivery_pending": True,
        "single_figure_count": 72,
        "mode_counts": {"mode_a": 36, "mode_b": 36},
        "composite_figure_count": 2,
        "style_manifest": relative(STYLE_MANIFEST),
        "style_manifest_sha256": sha256_file(STYLE_MANIFEST),
        "parent_manifest": relative(V5_MANIFEST),
        "parent_manifest_sha256": sha256_file(V5_MANIFEST),
        "charts": charts,
    }
    index_path = FIGURE_ROOT / "plot_manifest.json"
    write_manifest(index_path, index)

    (ANALYSIS_ROOT / "README.md").write_text(
        f"""# publication_clean_v9 PNG review layer

v9 is a solver-free display derivative of the accepted v5 point table. It
preserves all v5 CSV values, display-adjustment records, phase gaps, raw
provenance, and the current publication pointer. It emits 72 single PNG
charts (36 mode A and 36 mode B) plus two native composite PNG charts.

The layout follows a reviewed thesis-like exception to the default external
legend policy: each single chart has one in-axes legend, and each composite
has one compact shared legend in a case-selected sparse host panel. The key
labels alpha_T curves and open-circle/open-square first-order endpoints.
Linear panels use three labeled x-axis major ticks, one linear minor tick
between adjacent major ticks, and a MaxNLocator linear y-axis. Figure 1's
complete muB=900 MeV column uses log-y; other positive panels with
data_max/data_min >= 20 use local log-y. The v9 manifest records the parent
scale and the v9 display scale for each panel.

This is a PNG-only author-review stage. It remains
`manuscript_eligible=false`, `current_publication_layer=false`, and
`vector_delivery_pending=true`; it does not replace publication_clean_v5.
No solver, high-rate gate, numerical production, paper-project edit, or new
display value was used.

Parent v5 manifest SHA256: `{sha256_file(V5_MANIFEST)}`
Style manifest: `{relative(STYLE_MANIFEST)}`
Caption handoff: `{relative(caption_path)}`

Reproduce in an absent sibling directory:

```powershell
python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v9.py --png-review
```
""",
        encoding="utf-8",
    )

    package_manifest = {
        "schema": "publication_clean_v9_png_review_package_v1",
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "status": "author_review_required",
        "manuscript_eligible": False,
        "current_publication_layer": False,
        "delivery_stage": "png_review",
        "vector_delivery_pending": True,
        "author_acceptance": None,
        "solver_called": False,
        "canonical_data_modified": False,
        "numerical_status": "inherited_author_accepted_display_only",
        "raw_manuscript_eligible": False,
        "calculation_sha": package.get("calculation_sha"),
        "workflow_head_sha": package.get("workflow_head_sha"),
        "parent_manifest": relative(V5_MANIFEST),
        "parent_manifest_sha256": sha256_file(V5_MANIFEST),
        "inputs": inputs,
        "point_row_count": len(points),
        "inherited_table_hashes": table_hashes,
        "figure_index": relative(index_path),
        "figure_index_sha256": sha256_file(index_path),
        "outputs": [input_record(path, role="v9_display_document_or_table") for path in package_output_paths()],
    }
    write_manifest(ANALYSIS_ROOT / "manifest.json", package_manifest)
    print(
        json.dumps(
            {
                "figures": relative(FIGURE_ROOT),
                "manifest": relative(ANALYSIS_ROOT / "manifest.json"),
                "single_charts": 72,
                "composites": 2,
                "delivery_stage": "png_review",
                "vector_delivery_pending": True,
                "manuscript_eligible": False,
                "log_panel_count": sum(spec["axis_scale"] == "log" for spec in all_specs),
            },
            ensure_ascii=True,
        )
    )


if __name__ == "__main__":
    main()
