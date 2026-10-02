#!/usr/bin/env python3
"""Render the v11 PNG review layer from frozen v5 display values.

Reuse v10's point selection, phase-gap splitting, palette, and physical layout.
Only axes, numeric tick formatting, and endpoint-key typography change.
No solver, new smoothing, PDF delivery, or publication promotion is performed.
"""

from __future__ import annotations

import argparse
import datetime as dt
from decimal import Decimal
import importlib.util
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
from matplotlib.ticker import FixedLocator, FormatStrFormatter, FuncFormatter, LogLocator, MaxNLocator

from scripts.plotting.plot_manifest import (
    build_manifest, generator_record, input_record, runtime_record, sha256_file, write_manifest,
)
from scripts.plotting.plot_quality import export_figure
from scripts.plotting.plot_style import configure_matplotlib, load_profile
from scripts.plotting.validate_plot_artifact import validate_manifest_record
from scripts.plotting.plot_bundle import build_bundle


V10_SCRIPT = ROOT / "scripts/analysis/relaxtime/build_phase_guided_publication_clean_v10.py"
V10_SPEC = importlib.util.spec_from_file_location("publication_clean_v10_inputs_for_v11", V10_SCRIPT)
if V10_SPEC is None or V10_SPEC.loader is None:
    raise RuntimeError(f"cannot load the retained v10 renderer: {V10_SCRIPT}")
V10 = importlib.util.module_from_spec(V10_SPEC)
V10_SPEC.loader.exec_module(V10)
PARENT = V10.PARENT
V2 = V10.V2
V5_TABLE_ROOT = V10.V5_TABLE_ROOT
V5_MANIFEST = V10.V5_MANIFEST
CURRENT_POINTER = V10.CURRENT_POINTER
ANALYSIS_ROOT = V10.TRANSPORT_ANALYSIS_ROOT / "phase_guided_transport_publication_clean_v11_png_review"
FIGURE_ROOT = ROOT / "data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v11_png_review"
TABLE_ROOT = ANALYSIS_ROOT / "tables"
STYLE_MANIFEST = ANALYSIS_ROOT / "v11_display_style_manifest.json"
CAPTION_HANDOFF = ANALYSIS_ROOT / "caption_handoff.md"
OBSERVABLES = V10.OBSERVABLES
LABELS = V10.LABELS
UNITS = V10.UNITS
PANELS = V10.PANELS
COMPOSITES = V10.COMPOSITES
TAU_OBSERVABLES = tuple(name for name in OBSERVABLES if name.startswith("tau_"))
X_MAJOR_TICKS = V10.X_MAJOR_TICKS
LINEAR_MINOR_SUBDIVISIONS = V10.LINEAR_MINOR_SUBDIVISIONS
COMPOSITE_LAYOUT = V10.COMPOSITE_LAYOUT
COMPOSITE_SIZE_IN = V10.COMPOSITE_SIZE_IN
COMPOSITE_FONT_SIZE_PT = V10.COMPOSITE_FONT_SIZE_PT
COMPOSITE_LEGEND_FONT_SIZE_PT = V10.COMPOSITE_LEGEND_FONT_SIZE_PT
COMPOSITE_LEGEND_HOSTS = V10.COMPOSITE_LEGEND_HOSTS
ENDPOINT_LABELS = {"quark": "First-order (restored)", "hadron": "First-order (broken)"}
LOG_MAJOR_SUBS = (1.0, 2.0, 5.0)
LOG_MINOR_SUBS = tuple(range(2, 10))
ENDPOINT_TITLE = "First-order"
LOG_MIN_LABEL_FRACTION = 0.12
COMPOSITE_LEGEND_ANCHOR = (-0.02, 1.04)


def relative(path: Path) -> str:
    return V10.relative(path)


def visible_ticks(axis: Any) -> list[float]:
    low, high = axis.get_ylim()
    return [float(value) for value in axis.get_yticks() if low <= value <= high]


def format_log_tick(value: float, _position: int = 0) -> str:
    if not math.isfinite(value) or value <= 0:
        return ""
    number = Decimal(f"{value:.10g}")
    return format(number.normalize(), "f")


def readable_log_ticks(low: float, high: float, *, minimum: int = 3, maximum: int = 6) -> list[float]:
    """Prefer 1-2-5 ticks, supplement narrow views, and thin wide views."""
    if not (math.isfinite(low) and math.isfinite(high) and 0 < low < high):
        raise ValueError("log tick limits must be finite, positive, and increasing")

    def located(subs: tuple[float, ...]) -> list[float]:
        ticks = LogLocator(base=10, subs=subs, numticks=100).tick_values(low, high)
        return sorted({float(f"{value:.12g}") for value in ticks if low <= value <= high})

    span = math.log(high / low)
    preferred = located(LOG_MAJOR_SUBS)
    if len(preferred) > maximum:
        targets = [preferred[0] * (preferred[-1] / preferred[0]) ** (index / (maximum - 1)) for index in range(maximum)]
        preferred = sorted({min(preferred, key=lambda value: abs(math.log(value / target))) for target in targets})
    ticks = []
    for value in preferred:
        if not ticks or math.log(value / ticks[-1]) / span >= LOG_MIN_LABEL_FRACTION:
            ticks.append(value)
    candidates = located((1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9))
    if len(candidates) < minimum:
        candidates = sorted(set(candidates) | {
            float(f"{value:.12g}")
            for value in MaxNLocator(nbins=5, min_n_ticks=4).tick_values(low, high)
            if low <= value <= high
        })

    def separated(value: float) -> bool:
        return all(abs(math.log(value / tick)) / span >= LOG_MIN_LABEL_FRACTION for tick in ticks)

    # Add edge labels only when the nearest current label leaves a large gap;
    # close pairs such as 0.5/0.6 cannot coexist in the compact physical panel.
    for target, lower in ((low, True), (high, False)):
        if len(ticks) >= maximum:
            break
        edge = min(ticks) if lower and ticks else max(ticks) if ticks else target
        if ticks and abs(math.log(edge / target)) / span <= 0.15:
            continue
        available = [value for value in candidates if separated(value) and (not ticks or (value < edge if lower else value > edge))]
        if available:
            ticks.append(min(available, key=lambda value: abs(math.log(value / target))))
    while len(ticks) < minimum:
        available = [value for value in candidates if separated(value)]
        if not available:
            break
        ticks.append(max(available, key=lambda value: min((abs(math.log(value / tick)) for tick in ticks), default=span)))
    if len(ticks) < 3:
        raise ValueError("log view cannot support at least three readable numeric labels")
    return sorted(ticks)


def decimal_places_for_ticks(ticks: list[float]) -> int:
    if len(ticks) < 2:
        raise ValueError("decimal precision requires at least two major ticks")
    step = min(right - left for left, right in zip(ticks, ticks[1:]) if right > left)
    return max(0, -Decimal(f"{step:.12g}").normalize().as_tuple().exponent)


def capture_tick_spec(axis: Any, spec: dict[str, Any]) -> None:
    values = visible_ticks(axis)
    formatter = axis.yaxis.get_major_formatter()
    spec["y_limits"] = list(axis.get_ylim())
    spec["y_major_ticks"] = values
    spec["y_major_tick_labels"] = [formatter(value, index) for index, value in enumerate(values)]


def render_panel(
    axis: Any, mode: str, panel: str, observable: str, grouped: dict, gap_map: dict, profile: Any,
    *, composite: bool = False,
) -> dict[str, Any]:
    force_log = mode == "mode_a" and observable in TAU_OBSERVABLES
    spec = V10.render_panel(
        axis, mode, panel, observable, grouped, gap_map, profile,
        force_log=force_log, composite=composite,
    )
    if force_log:
        spec["axis_scale_reason"] = "uniform_log_for_mode_a_relaxation_family"
    if spec["axis_scale"] == "log":
        ticks = readable_log_ticks(*axis.get_ylim())
        axis.yaxis.set_major_locator(FixedLocator(ticks))
        axis.yaxis.set_major_formatter(FuncFormatter(format_log_tick))
        spec["y_tick_formatter_policy"] = "ordinary_numeric_log_labels_without_redundant_zeros"
        spec["y_tick_decimal_places"] = None
    else:
        precision = 3 if observable == "sigma_over_T" else decimal_places_for_ticks(visible_ticks(axis))
        axis.yaxis.set_major_formatter(FormatStrFormatter(f"%.{precision}f"))
        spec["y_tick_formatter_policy"] = "fixed_decimal_from_major_tick_spacing"
        spec["y_tick_decimal_places"] = precision
    capture_tick_spec(axis, spec)
    return spec


def apply_row_precision(axes: Any, specs: list[dict[str, Any]], observable: str) -> int | None:
    if all(axis.get_yscale() == "log" for axis in axes):
        return None
    if not all(axis.get_yscale() == "linear" for axis in axes):
        raise ValueError("a v11 composite row cannot mix linear and log axes")
    precision = max(decimal_places_for_ticks(visible_ticks(axis)) for axis in axes)
    if observable == "sigma_over_T":
        precision = 3
    for axis, spec in zip(axes, specs):
        axis.yaxis.set_major_formatter(FormatStrFormatter(f"%.{precision}f"))
        spec["y_tick_decimal_places"] = precision
        spec["y_tick_formatter_policy"] = "row_shared_fixed_decimal_from_finest_major_tick_spacing"
        capture_tick_spec(axis, spec)
    return precision


def legend_handles(specs: list[dict[str, Any]], profile: Any, *, endpoints: bool) -> list[Any]:
    handles = V10.legend_handles(specs, profile, include_endpoints=endpoints)
    if endpoints:
        for handle, label in zip(handles[len(specs[0]["series"]):], ENDPOINT_LABELS.values()):
            handle.set_label(label)
    return handles


def apply_legend(axis: Any, handles: list[Any], profile: Any, *, composite: bool, title: str | None = None) -> Any:
    font_size = COMPOSITE_LEGEND_FONT_SIZE_PT if composite else profile.data["legend_font_size_pt"]
    kwargs = {"loc": "upper left", "bbox_to_anchor": COMPOSITE_LEGEND_ANCHOR} if composite else {"loc": "best"}
    legend = axis.legend(
        handles=handles, **kwargs,
        ncol=1 if composite else (2 if len(handles) >= 4 else max(1, len(handles))),
        fontsize=font_size, title_fontsize=font_size, title=title,
        alignment="left", frameon=False,
        handlelength=1.1 if composite else 1.8,
        columnspacing=0.45 if composite else 0.5,
        handletextpad=0.3 if composite else 0.35,
        labelspacing=0.16 if composite else 0.25,
        borderpad=0.08 if composite else 0.25,
    )
    legend.set_zorder(10)
    return legend


def render_single(mode: str, panel: str, observable: str, grouped: dict, gap_map: dict, profile: Any) -> tuple[Any, list[dict]]:
    figure, axis = plt.subplots(figsize=(6.75, 4.6))
    figure.subplots_adjust(left=0.16, right=0.97, bottom=0.14, top=0.95)
    spec = render_panel(axis, mode, panel, observable, grouped, gap_map, profile)
    axis.set_xlabel(r"$\xi$")
    figure.suptitle(V10.panel_title(mode, panel), y=0.985)
    apply_legend(axis, legend_handles([spec], profile, endpoints=bool(spec["endpoints"])), profile, composite=False)
    return figure, [spec]


def render_composite(observables: list[str], grouped: dict, gap_map: dict, profile: Any) -> tuple[Any, list[dict]]:
    with matplotlib.rc_context({
        "font.size": COMPOSITE_FONT_SIZE_PT, "axes.labelsize": COMPOSITE_FONT_SIZE_PT,
        "axes.titlesize": COMPOSITE_FONT_SIZE_PT, "xtick.labelsize": COMPOSITE_FONT_SIZE_PT,
        "ytick.labelsize": COMPOSITE_FONT_SIZE_PT,
    }):
        case_name = next(name for name, values in COMPOSITES.items() if values == observables)
        figure, axes = plt.subplots(len(observables), 3, figsize=COMPOSITE_SIZE_IN[case_name], sharex=True)
        figure.subplots_adjust(**COMPOSITE_LAYOUT)
        specs = []
        for row, observable in enumerate(observables):
            row_specs = []
            for col, panel in enumerate(PANELS):
                axis = axes[row, col]
                spec = render_panel(axis, "mode_a", panel, observable, grouped, gap_map, profile, composite=True)
                row_specs.append(spec)
                if col:
                    axis.set_ylabel("")
                if row == 0:
                    axis.set_title(V10.panel_title("mode_a", panel), pad=4)
                if row == len(observables) - 1:
                    axis.set_xlabel(r"$\xi$")
                axis.text(0, 1.015, f"({chr(97 + row * 3 + col)})", transform=axis.transAxes,
                          ha="left", va="bottom", fontsize=10)
            apply_row_precision(axes[row], row_specs, observable)
            specs.extend(row_specs)
        handles = legend_handles(specs, profile, endpoints=True)
        handles[3].set_label("restored")
        handles[4].set_label("broken")
        for role, host in COMPOSITE_LEGEND_HOSTS[case_name].items():
            host_row, host_col = (int(value) for value in host.removeprefix("row").split("_col"))
            apply_legend(axes[host_row, host_col], handles[:3] if role == "parameters" else handles[3:],
                         profile, composite=True, title=ENDPOINT_TITLE if role == "endpoints" else None)
        return figure, specs


def write_chart(figure: Any, stem: Path, specs: list[dict], profile: Any,
                inputs: list[dict], font: dict, *, composite: bool) -> dict[str, Any]:
    outputs, quality = export_figure(figure, stem, profile, formats=("png",))
    axes = []
    for spec in specs:
        axes.extend([
            {"field": "xi", "source_unit": "dimensionless", "display_unit": "dimensionless",
             "label": r"$\xi$", "transform": "identity", "panel": spec["plot_panel"]},
            {"field": spec["observable"], "source_unit": UNITS[spec["observable"]],
             "display_unit": UNITS[spec["observable"]], "label": LABELS[spec["observable"]],
             "transform": spec["axis_scale"], "panel": spec["plot_panel"]},
        ])
    policy = "shared_in_panel_reviewed_geometry_checked" if composite else "best_in_axes_reviewed"
    manifest = build_manifest(
        asset_id=f"relaxtime.phase_guided.v11.png_review.{relative(stem)}",
        figure_family="phase_guided_transport", case_slug=stem.name, figure_mode="audit",
        semantic_status="author_review_display_derivative", style_profile=profile.profile_id,
        publication_scope="internal_review",
        generator=generator_record(Path(__file__), command="python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v11.py --png-review",
                                   runtime=runtime_record({"matplotlib": matplotlib.__version__, "font": font})),
        inputs=inputs, axes=axes, series=[item for spec in specs for item in spec["series"]], outputs=outputs,
        selection_rule="frozen v5 clean_value table in original xi order; phase gaps unchanged",
        interpolation_policy="inherited v5 display adjustments; no new interpolation or replacement",
        connector_policy="forbidden", missing_value_policy="preserve first-order and missing-support gaps",
        validation={"finite": True, "duplicate_keys": True, "support": True, "strict_gate": False},
        rendering={
            "column": "double_column", "figure_size_inches": quality["figure_size_inches"],
            "size_override_reason": "native manuscript composite at final physical size" if composite else None,
            "typography_exception": "dense_composite_review_compact_typography" if composite else None,
            "font_overrides_pt": {"labels_titles_ticks": COMPOSITE_FONT_SIZE_PT, "legend": COMPOSITE_LEGEND_FONT_SIZE_PT,
                                  "panel_labels": 10} if composite else None,
            "subplot_layout": dict(COMPOSITE_LAYOUT) if composite else None, "bbox_inches": None,
            "quality": quality, "legend_policy": policy, "legend_outside": False,
            "legend_host_panel": COMPOSITE_LEGEND_HOSTS.get(stem.name) if composite else None,
            "legend_location": {"loc": "upper left", "bbox_to_anchor": list(COMPOSITE_LEGEND_ANCHOR)} if composite else "best",
            "legend_contents": "alpha_T key and restored/broken endpoint key under the First-order heading" if composite else "panel parameter series and any actual First-order endpoint records",
            "legend_title": ENDPOINT_TITLE if composite else None, "legend_alignment": "left",
            "legend_title_weight": "normal", "legend_title_font_size_pt": COMPOSITE_LEGEND_FONT_SIZE_PT if composite else profile.data["legend_font_size_pt"],
            "legend_font_size_pt": COMPOSITE_LEGEND_FONT_SIZE_PT if composite else profile.data["legend_font_size_pt"],
            "legend_scope": {"parameters": "all panels" if composite else "series in this panel",
                             "endpoints": "muB900.0 alpha1.0 first-order endpoint curves only" if composite else "actual endpoints listed in panel_specs only"},
            "legend_marker_color_policy": "open markers match their actual endpoint data series",
            "case_layout_contract": "retained v10 thesis-like layout; new axes and key typography reviewed geometrically",
            "panel_label_policy": "uniform above-frame upper-left labels; no data area occupied",
            "color_route": "undecided_review", "grayscale_review": "required before author acceptance",
            "y_axis_policy": "all mode-A tau panels use log-y; transport composite linear; mode-B scale selection retained from v10",
            "tick_policy": {
                "x_major_ticks": list(X_MAJOR_TICKS), "linear_minor_subdivisions": LINEAR_MINOR_SUBDIVISIONS,
                "linear_y_major_locator": "v10 MaxNLocator; sigma/T exact integer-milli grid",
                "log_major_locator": "1-2-5 LogLocator candidates; supplemental neat values for narrow ranges; FixedLocator",
                "log_major_formatter": "ordinary numbers; redundant zeros removed; no offset",
                "log_minimum_label_separation_fraction": LOG_MIN_LABEL_FRACTION,
                "log_minor_locator": "LogLocator(base=10, subs=[2,3,4,5,6,7,8,9]); duplicate major positions removed",
                "direction": "in", "sides": ["top", "bottom", "left", "right"],
            },
            "panel_specs": specs, "x_tick_label_policy": "three labels at -0.5, 0, and 0.5",
            "linear_decimal_policy": "shared within each observable row, from its finest major tick spacing",
            "sigma_over_T_tick_policy": "exact integer-milli ticks; three decimal places; leading zero retained",
            "delivery_stage": "png_review", "vector_delivery_pending": True, "output_formats": ["png"],
        },
        calculation_sha=PARENT.V5.V4.V3.V2.V1.CALCULATION_SHA,
    )
    manifest.update({
        "manuscript_eligible": False, "current_publication_layer": False, "solver_called": False,
        "canonical_data_modified": False, "new_display_values": False, "delivery_stage": "png_review",
        "vector_delivery_pending": True, "workflow_head_sha": PARENT.V5.V4.V3.V2.V1.WORKFLOW_HEAD_SHA,
        "numerical_status": "inherited_author_accepted_display_only", "raw_manuscript_eligible": False,
    })
    violations = validate_manifest_record(manifest)
    plt.close(figure)
    if violations:
        raise ValueError(f"{stem.name}: " + "; ".join(violations))
    return {"stem": relative(stem), "manifest_record": manifest,
            "kind": "composite" if composite else "single", "mode_key": specs[0]["mode_key"],
            "axis_scales": sorted({spec["axis_scale"] for spec in specs}), "outputs": outputs}


def retained_v10_records() -> list[dict[str, Any]]:
    paths = [*V10.FIGURE_ROOT.rglob("*"), *V10.ANALYSIS_ROOT.rglob("*")]
    return [input_record(path, role="retained_v10_review_artifact") for path in sorted(paths) if path.is_file()]


def parent_records(package: dict, profile: Any) -> list[dict]:
    records = V10.parent_records(package, profile)
    records.extend(retained_v10_records())
    records.append(input_record(Path(__file__), role="v11_generator"))
    return records


def write_caption_handoff(points: list[dict[str, str]]) -> None:
    temperature = {}
    for row in points:
        if row["mode_key"] == "mode_a":
            key = (row["plot_panel"], row["plot_series"])
            value = float(row["T_MeV"])
            if key in temperature and not math.isclose(value, temperature[key], abs_tol=1e-9, rel_tol=0):
                raise ValueError(f"temperature varies along xi for {key}")
            temperature[key] = value
    mapping = "\n".join(
        f"| {panel} | " + " | ".join(f"{temperature[(panel, f'alpha{alpha:.1f}')]:.6f}" for alpha in (1, 1.1, 1.2)) + " |"
        for panel in PANELS
    )
    CAPTION_HANDOFF.write_text(r"""# publication_clean_v11 PNG review caption handoff

Figure identifiers follow this repository, not the reversed numbering in an
external review: Figure 1 is the four-row relaxation-time composite; Figure 2
is the three-row transport-coefficient composite.

Suggested caption text:

> Columns correspond to $\mu_B=0$, 450, and 900 MeV from left to right.
> The color and line-style key in panel (a) applies throughout: solid,
> dashed, and dash-dotted curves correspond to $\alpha_T=1.0$, 1.1, and 1.2.
> For $\mu_B=900\,\mathrm{MeV}$ and $\alpha_T=1.0$, open circles and squares
> mark the endpoints of the chirally restored and chirally broken branches
> at the first-order transition, respectively. Disconnected segments display
> these phase branches separately. Figure 1 rows show $\tau_u$, $\tau_s$,
> $\tau_{\bar u}$, and $\tau_{\bar s}$ in fm, all on logarithmic y axes.
> Figure 2 rows show $\eta/s$, $\zeta/s$, and $\sigma/T$ on linear y axes.
> Vertical ranges differ between panels. The corresponding fixed
> temperatures are listed below.

Do not say both legends apply to all panels: the parameter key is common,
whereas the endpoint key in (c) applies only to the endpoint-bearing curves.
Log-y displays relative changes; apparent curvature is not evidence of
absolute growth-rate saturation. Different y spans also prevent direct
comparison of visual slopes.

The retained v10 physical layout uses 11 pt labels/titles/ticks and 10.5 pt
legends. The First-order title is left-aligned, regular-weight, and the same
size as the parameter and endpoint entries. All new labels require fresh
curve/endpoint-intersection and text-clipping checks. Linear eta/s and zeta/s
rows have two decimal places, sigma/T has three; log labels use ordinary
numbers without redundant zeros. These are display rules, not uncertainty.

Exact frozen temperature mapping (MeV):

| panel | alpha_T=1.0 | alpha_T=1.1 | alpha_T=1.2 |
| --- | ---: | ---: | ---: |
""" + mapping + r"""

This package is PNG-only and review-only. Math scripts below 2 mm retain the
explicit dense_composite_review_compact_typography review exception; this
does not claim full APS submission compliance. v5 values, phase gaps, raw
provenance, and current=v5 remain unchanged. No solver, gate, new smoothing,
paper-project edit, PDF delivery, or manuscript eligibility is authorized.
""", encoding="utf-8")


def style_manifest_payload() -> dict[str, Any]:
    return {
        "schema": "publication_clean_v11_png_review_display_style_v1", "task_classification": "independent",
        "delivery_stage": "png_review", "manuscript_eligible": False, "current_publication_layer": False,
        "vector_delivery_pending": True, "values_unchanged": True,
        "parent_manifest": relative(V5_MANIFEST), "parent_manifest_sha256": sha256_file(V5_MANIFEST),
        "previous_review_manifest": relative(V10.ANALYSIS_ROOT / "manifest.json"),
        "previous_review_manifest_sha256": sha256_file(V10.ANALYSIS_ROOT / "manifest.json"),
        "generator": relative(Path(__file__)), "generator_sha256": sha256_file(Path(__file__)),
        "legend": {"title": ENDPOINT_TITLE, "alignment": "left", "weight": "normal",
                   "endpoint_labels": ENDPOINT_LABELS, "hosts": COMPOSITE_LEGEND_HOSTS,
                   "composite_font_size_pt": COMPOSITE_LEGEND_FONT_SIZE_PT,
                   "title_font_size_pt": COMPOSITE_LEGEND_FONT_SIZE_PT, "composite_anchor": list(COMPOSITE_LEGEND_ANCHOR),
                   "composite_parameters_scope": "all panels", "composite_endpoints_scope": "muB900.0 alpha1.0 first-order endpoint curves",
                   "single_endpoints_scope": "actual endpoints listed in each panel_specs record"},
        "axes": {"relaxation_composite": "all 12 panels log-y", "transport_composite": "all 9 panels linear",
                 "matching_mode_a_tau_singles": "all log-y", "mode_b_scales": "retained from v10",
                 "independent_y_ranges": True, "x_major_ticks": list(X_MAJOR_TICKS),
                 "linear_minor_subdivisions": LINEAR_MINOR_SUBDIVISIONS,
                 "log_major_subs": list(LOG_MAJOR_SUBS), "log_minor_subs": list(LOG_MINOR_SUBS),
                 "log_minimum_label_separation_fraction": LOG_MIN_LABEL_FRACTION,
                 "log_label_policy": "ordinary numeric 1-2-5; neat supplements for narrow views; no redundant zeros",
                 "linear_row_decimal_places": {"eta_over_s": 2, "zeta_over_s": 2, "sigma_over_T": 3}},
        "composite_layout": {"figure_sizes_inches": COMPOSITE_SIZE_IN, "subplot_layout": COMPOSITE_LAYOUT,
                             "labels_titles_ticks_pt": COMPOSITE_FONT_SIZE_PT,
                             "typography_exception": "dense_composite_review_compact_typography"},
        "figure_count": 72, "mode_counts": {"mode_a": 36, "mode_b": 36}, "composite_count": 2, "formats": ["png"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--png-review", action="store_true")
    parser.add_argument("--composites-only", action="store_true", help="temporary preview in a separate absent directory")
    parser.add_argument("--output-root", type=Path)
    args = parser.parse_args()
    if not args.png_review and not args.composites_only:
        raise ValueError("v11 is PNG-review-only; pass --png-review")
    if args.output_root is not None and not args.composites_only:
        raise ValueError("--output-root is reserved for a temporary composites-only preview")
    if args.composites_only and (args.output_root is None or args.output_root.resolve() == FIGURE_ROOT.resolve()):
        raise ValueError("composites-only requires a separate absent --output-root")
    if not args.composites_only and (FIGURE_ROOT.exists() or ANALYSIS_ROOT.exists()):
        raise FileExistsError("refusing to overwrite an existing v11 review package")
    points, package, _, _, gaps, table_hashes = PARENT.load_v5_inputs()
    keys = [(row["mode_key"], row["plot_panel"], row["plot_series"], row["observable"], row["xi"]) for row in points]
    if len(keys) != len(set(keys)):
        raise ValueError("duplicate frozen v5 point keys")
    for record in package["outputs"]:
        if sha256_file(ROOT / record["path"]) != record["sha256"]:
            raise ValueError(f"v5 package hash mismatch: {record['path']}")
    profile = load_profile("candidate_aps_v2")
    font = configure_matplotlib(profile)
    grouped, gap_map = V10.group_inputs(points, gaps)
    if args.composites_only:
        args.output_root.mkdir(parents=True, exist_ok=False)
        for name, observables in COMPOSITES.items():
            figure, _ = render_composite(observables, grouped, gap_map, profile)
            figure.savefig(args.output_root / f"{name}.png", dpi=150)
            plt.close(figure)
        return

    inputs = parent_records(package, profile)
    ANALYSIS_ROOT.mkdir(parents=True, exist_ok=False)
    shutil.copytree(V5_TABLE_ROOT, TABLE_ROOT)
    charts = []
    for mode in ("mode_a", "mode_b"):
        for panel in sorted({key[1] for key in grouped if key[0] == mode}):
            for observable in OBSERVABLES:
                figure, specs = render_single(mode, panel, observable, grouped, gap_map, profile)
                mode_dir = "mode_a_fixed_muB_phase_scaled" if mode == "mode_a" else "mode_b_fixed_T_sparse_muB"
                stem = FIGURE_ROOT / mode_dir / f"plot_panel={panel}" / f"{observable}_vs_xi"
                charts.append(write_chart(figure, stem, specs, profile, inputs, font, composite=False))
        print(f"[v11] {mode}: 36 single charts validated", flush=True)
    for name, observables in COMPOSITES.items():
        figure, specs = render_composite(observables, grouped, gap_map, profile)
        charts.append(write_chart(figure, FIGURE_ROOT / "composites" / name, specs, profile, inputs, font, composite=True))
    for record in inputs:
        if sha256_file(ROOT / record["path"]) != record["sha256"]:
            raise ValueError(f"input changed during render: {record['path']}")
    for name, expected in table_hashes.items():
        if sha256_file(TABLE_ROOT / name) != expected or sha256_file(V5_TABLE_ROOT / name) != expected:
            raise ValueError(f"v11 changed inherited table bytes: {name}")
    write_caption_handoff(points)
    style = style_manifest_payload()
    style["generated_at_utc"] = dt.datetime.now(dt.timezone.utc).isoformat()
    write_manifest(STYLE_MANIFEST, style)
    index_path = FIGURE_ROOT / "plot_manifest.json"
    write_manifest(index_path, build_bundle({
        "schema": "publication_clean_v11_png_review_figure_index_v1", "status": "author_review_required",
        "manuscript_eligible": False, "current_publication_layer": False, "solver_called": False,
        "canonical_data_modified": False, "delivery_stage": "png_review", "vector_delivery_pending": True,
        "single_figure_count": 72, "composite_figure_count": 2, "mode_counts": {"mode_a": 36, "mode_b": 36},
        "style_manifest": relative(STYLE_MANIFEST), "style_manifest_sha256": sha256_file(STYLE_MANIFEST), "charts": charts,
    }, [chart["manifest_record"] for chart in charts]))
    (ANALYSIS_ROOT / "README.md").write_text(f"""# publication_clean_v11 PNG review layer

v11 reuses frozen v5 data and v10's reviewed physical layout. Figure 1 is the
12-panel relaxation-time composite, now uniformly log-y with ordinary numeric
ticks. Figure 2 is the nine-panel transport composite, still linear; eta/s and
zeta/s share two decimal places within each row and sigma/T shares three.
Matching mode-A tau singles use log-y. The First-order endpoint key in (c)
is left-aligned and the same 10.5 pt size as the parameter key in (a).
Caption scope separates the common parameter key from muB900.0 alpha1.0
endpoint-bearing curves. Each panel retains its independent y range.

The package contains 72 single and two composite 600 dpi PNGs. It remains
manuscript_eligible=false, current_publication_layer=false, and
vector_delivery_pending=true. The v10 dense-composite typography review
exception remains explicit; this is not full APS submission compliance.
No v5 values, phase gaps, raw provenance, current pointer, prior v10 file,
solver, numerical gate, smoothing, PDF, or paper-project file was changed.

Manifest and caption: `manifest.json`, `v11_display_style_manifest.json`,
`caption_handoff.md`. Retained v10 artifact hashes are included as inputs.
This independent figure task does not replace the primary task-ledger track.

```powershell
python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v11.py --png-review
```
""", encoding="utf-8")
    outputs = [ANALYSIS_ROOT / "README.md", CAPTION_HANDOFF, STYLE_MANIFEST, index_path,
               *sorted(TABLE_ROOT.glob("*.csv")), *sorted(FIGURE_ROOT.rglob("*.png"))]
    write_manifest(ANALYSIS_ROOT / "manifest.json", {
        "schema": "publication_clean_v11_png_review_package_v1", "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "task_classification": "independent", "status": "author_review_required", "author_acceptance": None,
        "manuscript_eligible": False, "current_publication_layer": False, "solver_called": False,
        "canonical_data_modified": False, "new_display_values": False, "delivery_stage": "png_review",
        "vector_delivery_pending": True, "numerical_status": "inherited_author_accepted_display_only",
        "raw_manuscript_eligible": False, "calculation_sha": package.get("calculation_sha"),
        "workflow_head_sha": package.get("workflow_head_sha"), "parent_manifest": relative(V5_MANIFEST),
        "parent_manifest_sha256": sha256_file(V5_MANIFEST), "inputs": inputs, "point_row_count": len(points),
        "inherited_table_hashes": table_hashes, "figure_index": relative(index_path),
        "figure_index_sha256": sha256_file(index_path),
        "outputs": [input_record(path, role="v11_display_document_or_table") for path in outputs],
    })
    print(f"[v11] 74 PNGs validated: {relative(FIGURE_ROOT)}", flush=True)


if __name__ == "__main__":
    main()
