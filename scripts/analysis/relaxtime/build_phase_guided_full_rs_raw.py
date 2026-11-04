#!/usr/bin/env python3
"""Render the accepted full-RS grid as an unmodified, v13-shaped PNG review."""

from __future__ import annotations

import argparse
from collections import defaultdict
from contextlib import contextmanager
from dataclasses import dataclass
from decimal import Decimal
from itertools import product
import json
import math
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import csv
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FormatStrFormatter, FuncFormatter, LogLocator, MaxNLocator
from PIL import Image, ImageOps

from scripts.plotting.plot_accessibility import accessibility_record
from scripts.plotting.plot_bundle import build_bundle
from scripts.plotting.plot_manifest import (
    build_manifest, generator_record, input_record, output_record, runtime_record,
    sha256_file, write_manifest,
)
from scripts.plotting.plot_quality import (
    export_figure, inspect_export, measure_figure, placement_limits, visible_legends,
)
from scripts.plotting.plot_series import SeriesStyle, SeriesStyleRegistry
from scripts.plotting.plot_style import configure_axis_ticks, load_profile, style_context
from scripts.plotting.validate_plot_artifact import validate_manifest, validate_manifest_record

PRODUCTION_ROOT = ROOT / "data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021/full_rs_full_grid_v1"
DEFAULT_CASE = "full_rs_raw_20261027_png_review"
FIGURE_PARENT = ROOT / "data/outputs/figures/relaxtime/transport/phase_guided"
RESULT_PARENT = ROOT / "data/outputs/results/relaxtime/transport/phase_guided"
CALCULATION_SHA = "2e5c411c716d5cbc524c3fc6cfb335fe05b39444"
MODE_DIRS = {"mode_a": "mode_a_fixed_muB_phase_scaled", "mode_b": "mode_b_fixed_T_sparse_muB"}
PANELS = {"mode_a": ("muB0.0", "muB450.0", "muB900.0"), "mode_b": ("T120.0", "T160.0", "T200.0")}
PARAMETERS = {"mode_a": ("alpha1.0", "alpha1.1", "alpha1.2"), "mode_b": ("muB0.0", "muB450.0", "muB900.0")}
OBSERVABLES = ("tau_u", "tau_d", "tau_s", "tau_ubar", "tau_dbar", "tau_sbar",
               "eta", "zeta", "sigma", "eta_over_s", "zeta_over_s", "sigma_over_T")
UNITS = {**{name: "fm" for name in OBSERVABLES if name.startswith("tau_")},
         "eta": "fm^-3", "zeta": "fm^-3", "sigma": "fm^-1",
         "eta_over_s": "dimensionless", "zeta_over_s": "dimensionless", "sigma_over_T": "dimensionless"}
LABELS = {
    **{f"tau_{q}": rf"$\tau_{{{q}}}\;(\mathrm{{fm}})$" for q in ("u", "d", "s")},
    **{f"tau_{q}bar": rf"$\tau_{{\bar{{{q}}}}}\;(\mathrm{{fm}})$" for q in ("u", "d", "s")},
    "eta": r"$\eta\;(\mathrm{fm}^{-3})$", "zeta": r"$\zeta\;(\mathrm{fm}^{-3})$",
    "sigma": r"$\sigma\;(\mathrm{fm}^{-1})$", "eta_over_s": r"$\eta/s$",
    "zeta_over_s": r"$\zeta/s$", "sigma_over_T": r"$\sigma/T$",
}
COMPOSITES = {
    "figure1_relaxation_times_comparison": ("tau_u", "tau_s", "tau_ubar", "tau_sbar"),
    "figure2_transport_coefficients_comparison": ("eta_over_s", "zeta_over_s", "sigma_over_T"),
}
COMPOSITE_SIZES = {"figure1_relaxation_times_comparison": (6.75, 7.1),
                   "figure2_transport_coefficients_comparison": (6.75, 5.9)}
DETAILS = (("muB0.0", "tau_u", "(a)", (0.57, 0.66)),
           ("muB0.0", "tau_ubar", "(g)", (0.57, 0.66)),
           ("muB450.0", "tau_ubar", "(h)", (0.375, 0.43)),
           ("muB900.0", "tau_ubar", "(i)", (0.27, 0.41)))
WIDTH_IN, LABEL_PT, TICK_PT, LEGEND_PT = 6.75, 13.0, 11.0, 11.0
ENDPOINT_MARKERS = {"quark": "o", "hadron": "s"}
LAYOUT_ERRORS = ("clipped_text", "text_overlap_pairs", "legend_curve_overlap_count",
                 "legend_landmark_overlap_count", "legend_in_axes_overflow_count", "legend_pair_overlap_count")


@dataclass
class RawInputs:
    groups: dict
    gaps: dict
    sources: dict
    evidence_paths: list[Path]
    acceptance: dict


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def display_path(path):
    path = path.resolve()
    return path.relative_to(ROOT).as_posix() if path.is_relative_to(ROOT) else path.as_posix()


def validate_rows(rows, mode):
    """Validate plotted fields; preserve source strings and finite quality flags."""
    required = {*OBSERVABLES, "xi", "T_MeV", "muq_MeV", "muB_MeV", "plot_panel",
                "plot_series", "mode", "converged", "quality_flag", "quality_reason",
                "phase_curr", "phase_structure", "alpha_T", "T_phase_base_MeV"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"{mode}: missing raw source columns")
    groups, seen = defaultdict(list), set()
    for csv_row, raw in enumerate(rows, start=2):
        row = dict(raw, csv_row=csv_row)
        panel, parameter = row["plot_panel"], row["plot_series"]
        if row["mode"] != MODE_DIRS[mode] or panel not in PANELS[mode] or parameter not in PARAMETERS[mode]:
            raise ValueError(f"{mode}: unexpected panel/parameter identity at row {csv_row}")
        values = {name: float(row[name]) for name in (*OBSERVABLES, "xi", "T_MeV", "muq_MeV", "muB_MeV")}
        if not all(math.isfinite(value) for value in values.values()):
            raise ValueError(f"{mode}: non-finite plotted input at row {csv_row}")
        if any(values[name] <= 0 for name in OBSERVABLES):
            raise ValueError(f"{mode}: non-positive transport input at row {csv_row}")
        if row["converged"] != "true":
            raise ValueError(f"{mode}: unconverged row {csv_row}")
        if row["quality_flag"] not in ("true", "false"):
            raise ValueError(f"{mode}: invalid quality flag at row {csv_row}")
        if row["quality_flag"] == "true" and row["quality_reason"] != "tau_u_ubar_ratio_high":
            raise ValueError(f"{mode}: unreviewed quality flag at row {csv_row}")
        if not math.isclose(values["muB_MeV"], 3 * values["muq_MeV"], abs_tol=1e-9):
            raise ValueError(f"{mode}: inconsistent chemical-potential identity")
        if mode == "mode_a":
            if (float(panel[3:]) != values["muB_MeV"] or float(parameter[5:]) != float(row["alpha_T"])
                    or not math.isclose(values["T_MeV"], float(row["alpha_T"]) * float(row["T_phase_base_MeV"]),
                                        rel_tol=1e-12)):
                raise ValueError(f"{mode}: inconsistent phase-scaled temperature identity")
        elif float(panel[1:]) != values["T_MeV"] or float(parameter[3:]) != values["muB_MeV"]:
            raise ValueError(f"{mode}: inconsistent fixed-temperature identity")
        identity = (panel, parameter, values["xi"])
        if identity in seen:
            raise ValueError(f"{mode}: duplicate source point {identity}")
        seen.add(identity)
        groups[(mode, panel, parameter)].append(row)
    for key, points in groups.items():
        points.sort(key=lambda row: float(row["xi"]))
        if len({row["T_MeV"] for row in points}) != 1:
            raise ValueError(f"{key}: changing temperature inside a fixed curve")
    return dict(groups)


def boundary_map(audit, groups):
    """Use only reviewed first-order pairs, never threshold/roughness-based masks."""
    result = defaultdict(set)
    for pair in audit["known_boundary_pairs"]:
        if pair["reference_first_order_boundary"] is not True:
            raise ValueError("boundary has no first-order reference evidence")
        mode = "mode_a" if pair["mode"] == "A" else "mode_b"
        mu = pair["muB_MeV"]
        parameter = pair["curve_parameter"]
        key = ((mode, f"muB{mu:.1f}", f"alpha{parameter:.1f}") if mode == "mode_a"
               else (mode, f"T{parameter:.1f}", f"muB{mu:.1f}"))
        left, right = float(pair["xi_left"]), float(pair["xi_right"])
        rows = groups[key]
        matches = [(a, b) for a, b in zip(rows, rows[1:])
                   if (float(a["xi"]), float(b["xi"])) == (left, right)]
        if len(matches) != 1 or {matches[0][0]["phase_curr"], matches[0][1]["phase_curr"]} != {"quark", "hadron"}:
            raise ValueError(f"{key}: reviewed gap does not match adjacent source branches")
        result[key].add((left, right))
    return {key: sorted(value) for key, value in result.items()}


def load_inputs(production_root=PRODUCTION_ROOT):
    acceptance_path = production_root / "final_review_v1/full_grid_acceptance.json"
    audit_path = production_root / "final_review_v1/background/background_audit.json"
    grid = production_root / "returned/transfer_090/reviews/transfer_090"
    acceptance, audit = read_json(acceptance_path), read_json(audit_path)
    if (acceptance["full_grid_numerical_acceptance_passed"] is not True
            or acceptance["artifact_contract_passed"] is not True
            or acceptance["calculation_sha"] != CALCULATION_SHA
            or acceptance["points"] != 1819):
        raise ValueError("the frozen full-RS grid has not passed its numerical acceptance")
    if sha256_file(audit_path) != acceptance["background_audit_sha256"] or audit["all_backgrounds_passed"] is not True:
        raise ValueError("background acceptance evidence changed")
    review_path = grid / "reviewed_grid.json"
    if sha256_file(review_path) != acceptance["grid_review_sha256"]:
        raise ValueError("grid review evidence changed")
    evidence = [acceptance_path, audit_path, review_path]
    for relative, digest in acceptance["merged_files"].items():
        path = grid / relative
        if sha256_file(path) != digest:
            raise ValueError(f"frozen input hash mismatch: {relative}")
        evidence.append(path)
    groups, sources = {}, {}
    for mode, expected in (("mode_a", 910), ("mode_b", 909)):
        path = grid / mode / "phase_guided_transport_scan.csv"
        with path.open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
        if len(rows) != expected:
            raise ValueError(f"{mode}: expected {expected} source points")
        groups.update(validate_rows(rows, mode))
        sources[mode] = path
    expected_keys = {(mode, panel, parameter) for mode in MODE_DIRS
                     for panel in PANELS[mode] for parameter in PARAMETERS[mode]}
    if set(groups) != expected_keys:
        raise ValueError("incomplete raw panel/parameter coverage")
    flags = sum(row["quality_flag"] == "true" for rows in groups.values() for row in rows)
    if flags != acceptance["quality_flags_retained"]:
        raise ValueError("retained quality flags differ from the accepted grid")
    gaps = boundary_map(audit, groups)
    if sum(map(len, gaps.values())) != 2:
        raise ValueError("expected exactly two reviewed first-order gaps")
    return RawInputs(groups, gaps, sources, evidence, acceptance)


def split_segments(rows, gaps):
    boundaries = set(gaps)
    segments = [[]]
    for index, row in enumerate(rows):
        if index and (float(rows[index - 1]["xi"]), float(row["xi"])) in boundaries:
            segments.append([])
        segments[-1].append(row)
    if len(segments) != len(boundaries) + 1:
        raise ValueError("gap must identify adjacent original points")
    return segments


def make_registry(mode, profile):
    # A parameter identity is shared across facets; panel + observable + identity
    # fully identifies each raw curve in panel_specs/source_curves.
    return SeriesStyleRegistry(f"full_rs_raw_{mode}_parameters_v1", tuple(
        SeriesStyle.from_profile(profile, f"{mode}.{parameter}", index)
        for index, parameter in enumerate(PARAMETERS[mode])))


def panel_title(mode, panel):
    return (rf"$\mu_B={float(panel[3:]):.0f}\;\mathrm{{MeV}}$" if mode == "mode_a"
            else rf"$T={float(panel[1:]):.0f}\;\mathrm{{MeV}}$")


def visible_ticks(axis):
    low, high = axis.get_ylim()
    return [float(value) for value in axis.get_yticks() if low <= value <= high]


def decimal_places(ticks):
    step = min(b - a for a, b in zip(ticks, ticks[1:]))
    return max(0, -Decimal(f"{step:.12g}").normalize().as_tuple().exponent)


def readable_log_ticks(low, high):
    """Keep ordinary 1/2/5-style labels readable in the v13 physical layout."""
    span = math.log(high / low)
    def located(subs):
        return sorted({float(f"{v:.12g}") for v in LogLocator(subs=subs, numticks=100).tick_values(low, high)
                       if low <= v <= high})
    preferred = located((1, 2, 5))
    ticks = []
    for value in preferred:
        if not ticks or math.log(value / ticks[-1]) / span >= 0.15:
            ticks.append(value)
    candidates = located((1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9))
    if len(candidates) < 3:
        candidates += [float(v) for v in MaxNLocator(nbins=5).tick_values(low, high) if low <= v <= high]
    while len(ticks) < 3:
        available = [v for v in candidates if all(abs(math.log(v / t)) / span >= 0.15 for t in ticks)]
        if not available:
            raise ValueError("cannot place three readable log labels")
        ticks.append(max(available, key=lambda v: min((abs(math.log(v / t)) for t in ticks), default=span)))
    return sorted(ticks)


def capture_ticks(axis, spec):
    ticks = visible_ticks(axis)
    spec.update(y_limits=list(axis.get_ylim()), x_limits=list(axis.get_xlim()), y_major_ticks=ticks,
                y_major_tick_labels=[axis.yaxis.get_major_formatter()(value, i) for i, value in enumerate(ticks)])


def render_panel(axis, mode, panel, observable, data, profile, *, composite=False, detail=None):
    registry = make_registry(mode, profile)
    curves, endpoints, values = [], [], []
    for parameter, style in zip(PARAMETERS[mode], registry.styles):
        key = (mode, panel, parameter)
        rows, gaps = data.groups[key], data.gaps.get(key, [])
        segments = split_segments(rows, gaps)
        for segment in segments:
            line, = axis.plot([float(r["xi"]) for r in segment], [float(r[observable]) for r in segment],
                              **style.plot_kwargs())
            line.get_path().should_simplify = False
        for left, right in gaps:
            for role, xi in (("left", left), ("right", right)):
                point = next(r for r in rows if float(r["xi"]) == xi)
                phase, value = point["phase_curr"], float(point[observable])
                axis.scatter([xi], [value], marker=ENDPOINT_MARKERS[phase],
                             s=profile.data["marker_size_pt"] ** 2, facecolor="white",
                             edgecolor=style.color, linewidth=0.8, zorder=6)
                endpoints.append(dict(series=parameter, role=role, phase=phase, xi=xi, value=value,
                                      color=style.color, csv_row=point["csv_row"],
                                      meaning="nearest original sample on this side of the reviewed discontinuity"))
        curves.append(dict(series_id=style.series_id, curve_id=f"{mode}.{panel}.{observable}.{parameter}",
                           parameter=parameter, mode_key=mode, plot_panel=panel, observable=observable,
                           source_csv=display_path(data.sources[mode]),
                           source_column=observable, source_rows=[r["csv_row"] for r in rows],
                           row_count=len(rows), segment_count=len(segments),
                           segment_source_rows=[[r["csv_row"] for r in segment] for segment in segments],
                           quality_flags_retained=sum(r["quality_flag"] == "true" for r in rows),
                           color=style.color, linestyle=style.linestyle, marker=style.marker))
        values.extend(float(row[observable]) for row in rows)
    force_log = mode == "mode_a" and observable.startswith("tau_")
    scale = "log" if force_log or max(values) / min(values) >= 20 else "linear"
    reason = "mode_A_relaxation_family_log" if force_log else "raw_dynamic_range_at_least_20" if scale == "log" else "linear_raw_range"
    if detail is not None:
        scale, reason = "linear", "explicit_v13_comparable_low_value_viewport"
    axis.set_yscale(scale)
    axis.set_xlim(-0.52, 0.52)
    axis.set_xticks([-0.5, 0, 0.5], ["-0.5", "0", "0.5"])
    axis.margins(y=0.06)
    if detail is not None:
        axis.set_ylim(*detail)
    if scale == "log":
        axis.yaxis.set_major_locator(FixedLocator(readable_log_ticks(*axis.get_ylim())))
        axis.yaxis.set_major_formatter(FuncFormatter(lambda value, _: format(Decimal(f"{value:.10g}").normalize(), "f")))
    else:
        axis.yaxis.set_major_locator(MaxNLocator(nbins=4 if composite else 5, min_n_ticks=3,
                                                 steps=[1, 2, 5, 10]))
        precision = max(3 if observable == "sigma_over_T" else 0, decimal_places(visible_ticks(axis)))
        axis.yaxis.set_major_formatter(FormatStrFormatter(f"%.{precision}f"))
    configure_axis_ticks(axis, profile)
    axis.tick_params(axis="y", which="major", pad=3)
    axis.set_ylabel(LABELS[observable], labelpad=3)
    spec = dict(mode_key=mode, plot_panel=panel, observable=observable, series=curves, endpoints=endpoints,
                axis_scale=scale, axis_scale_reason=reason, data_min=min(values), data_max=max(values),
                first_order_gap_present=bool(endpoints))
    capture_ticks(axis, spec)
    return spec


def clear_layout(quality):
    return not any(quality[key] for key in LAYOUT_ERRORS)


def place_legends(figure, specs, profile, *, hosts=(0, 2), single=False, endpoints=True):
    """Retain v13's compact keys and deterministic, measured corner selection."""
    mode = specs[0]["mode_key"]
    registry = make_registry(mode, profile)
    labels = (["1.0", "1.1", "1.2"] if mode == "mode_a" else ["0", "450", "900"])
    title = r"$\alpha_T$" if mode == "mode_a" else r"$\mu_B\;(\mathrm{MeV})$"
    handles = [Line2D([], [], color=s.color, linestyle=s.linestyle,
                      linewidth=profile.data["line_width_pt"], label=label)
               for s, label in zip(registry.styles, labels)]
    keys = [("parameters", title, handles, hosts[0])]
    targets = sorted({(s["mode_key"], s["plot_panel"], e["series"])
                      for s in specs for e in s["endpoints"]})
    if endpoints and targets:
        colors = {e["color"] for s in specs for e in s["endpoints"]}
        color = next(iter(colors)) if len(colors) == 1 else "#000000"
        handles = [Line2D([], [], linestyle="None", marker=marker, markerfacecolor="white",
                          markeredgecolor=color, markersize=profile.data["marker_size_pt"], label=label)
                   for marker, label in (("o", "restored"), ("s", "broken"))]
        keys.append(("endpoints", "First-order", handles, hosts[1]))

    def create(key, location, *, outside=False, anchor=None):
        role, title, group, host = key
        owner = figure if outside else figure.axes[host]
        if not outside and owner.get_legend() is not None:
            owner.add_artist(owner.get_legend())
        legend = owner.legend(
            handles=group, title=title, loc=location, bbox_to_anchor=anchor,
            ncol=len(group) if outside else 1, frameon=False, alignment="left",
            fontsize=LEGEND_PT, title_fontsize=LABEL_PT if role == "parameters" else LEGEND_PT,
            handlelength=1.35 if role == "parameters" else 0.85, handletextpad=0.3,
            columnspacing=0.9, labelspacing=0.10, borderpad=0.05, borderaxespad=0.25,
        )
        record = dict(role=role, placement="outside_axes" if outside else "in_axes",
                      host_axes_index=None if outside else host, location=location, title=title,
                      labels=[handle.get_label() for handle in group],
                      scope="all parameter curves" if role == "parameters" else "listed sampled first-order branch pairs only")
        if role == "parameters":
            record["series_keys"] = list(PARAMETERS[mode])
        else:
            record["applies_to"] = [dict(mode_key=m, plot_panel=p, plot_series=s) for m, p, s in targets]
        return legend, record

    for locations in product(("upper left", "upper right", "lower left", "lower right"), repeat=len(keys)):
        attempt = [create(key, location) for key, location in zip(keys, locations)]
        if clear_layout(measure_figure(figure, intended_width_inches=WIDTH_IN)):
            break
        for legend, _ in attempt:
            legend.remove()
    else:
        height = figure.get_size_inches()[1]
        title_space = 0.28 if single else 0
        reserved = title_space + (1.05 if len(keys) == 2 else 0.65)
        figure.subplots_adjust(top=1 - reserved / height)
        attempt = [create(key, "upper center", outside=True,
                          anchor=(0.55, 1 - (title_space + index * 0.46) / height))
                   for index, key in enumerate(keys)]
        quality = measure_figure(figure, intended_width_inches=WIDTH_IN)
        if not clear_layout(quality):
            raise ValueError("no unobstructed legend layout: " + json.dumps(
                {key: quality[key] for key in LAYOUT_ERRORS if quality[key]}, ensure_ascii=False))
    declarations = {id(legend): record for legend, record in attempt}
    return [{**declarations[id(legend)], "legend_index": index,
             "selection_reason": "first clear declared corner combination" if legend.axes is not None
             else "corners obstruct data or labels; compact top fallback"}
            for index, legend in enumerate(visible_legends(figure))]


@contextmanager
def figure_context(profile, *, compact=False):
    with style_context(profile) as context:
        with matplotlib.rc_context({"axes.labelsize": LABEL_PT, "axes.titlesize": LABEL_PT,
                                    "xtick.labelsize": TICK_PT if compact else LABEL_PT,
                                    "ytick.labelsize": TICK_PT if compact else LABEL_PT,
                                    "path.simplify": False}):
            yield context


def render_single(mode, panel, observable, data, profile):
    with figure_context(profile) as context:
        figure, axis = plt.subplots(figsize=(WIDTH_IN, 4.6))
        context.bind(figure)
        spec = render_panel(axis, mode, panel, observable, data, profile)
        figure.subplots_adjust(left=0.16, right=0.97, bottom=0.14, top=0.91)
        axis.set_xlabel(r"$\xi$")
        figure.suptitle(panel_title(mode, panel), y=0.99)
        placements = place_legends(figure, [spec], profile, hosts=(0, 0), single=True)
        return figure, [spec], placements


def render_composite(name, data, profile):
    observables = COMPOSITES[name]
    with figure_context(profile, compact=True) as context:
        figure, axes = plt.subplots(len(observables), 3, figsize=COMPOSITE_SIZES[name], sharex=True)
        context.bind(figure)
        figure.subplots_adjust(left=0.115, right=0.98, bottom=0.085, top=0.94, hspace=0.20, wspace=0.31)
        specs = []
        for row, observable in enumerate(observables):
            row_specs = []
            for col, panel in enumerate(PANELS["mode_a"]):
                axis = axes[row, col]
                spec = render_panel(axis, "mode_a", panel, observable, data, profile, composite=True)
                row_specs.append(spec)
                if col:
                    axis.set_ylabel("")
                if row == 0:
                    axis.set_title(panel_title("mode_a", panel), pad=5)
                if row == len(observables) - 1:
                    axis.set_xlabel(r"$\xi$")
                axis.text(0, 1.02, f"({chr(97 + row * 3 + col)})", transform=axis.transAxes,
                          ha="left", va="bottom", fontsize=TICK_PT)
            if all(axis.get_yscale() == "linear" for axis in axes[row]):
                precision = max(decimal_places(visible_ticks(axis)) for axis in axes[row])
                precision = max(precision, 3 if observable == "sigma_over_T" else 0)
                for axis, spec in zip(axes[row], row_specs):
                    axis.yaxis.set_major_formatter(FormatStrFormatter(f"%.{precision}f"))
                    capture_ticks(axis, spec)
            specs.extend(row_specs)
        placements = place_legends(figure, specs, profile)
        return figure, specs, placements


def render_detail(data, profile):
    with figure_context(profile, compact=True) as context:
        figure, axes = plt.subplots(2, 2, figsize=(WIDTH_IN, 5.5), sharex=True)
        context.bind(figure)
        figure.subplots_adjust(left=0.12, right=0.97, bottom=0.12, top=0.92, hspace=0.28, wspace=0.32)
        specs = []
        for axis, (panel, observable, source_panel, limits) in zip(axes.flat, DETAILS):
            spec = render_panel(axis, "mode_a", panel, observable, data, profile, composite=True, detail=limits)
            axis.set_title(panel_title("mode_a", panel), pad=5)
            axis.text(0, 1.02, source_panel, transform=axis.transAxes, va="bottom", fontsize=TICK_PT)
            spec.update(source_figure="figure1_relaxation_times_comparison", source_panel=source_panel,
                        viewport_clipping="all original vertices retained; values outside this explicit window are not visible")
            specs.append(spec)
        for axis in axes[-1]:
            axis.set_xlabel(r"$\xi$")
        placements = place_legends(figure, specs, profile, endpoints=False)
        return figure, specs, placements


def verify_raw_vertices(figure, specs, data):
    """Check actual artists, including cropped detail paths, against source rows."""
    vertices = 0
    for axis, spec in zip(figure.axes, specs, strict=True):
        expected = []
        for curve in spec["series"]:
            key = (spec["mode_key"], spec["plot_panel"], curve["parameter"])
            expected.extend((curve["series_id"], segment) for segment in
                            split_segments(data.groups[key], data.gaps.get(key, [])))
        if len(axis.lines) != len(expected):
            raise ValueError("raw curve segment count changed")
        for line, (identity, rows) in zip(axis.lines, expected):
            if (line.get_gid() != identity
                    or list(line.get_xdata()) != [float(row["xi"]) for row in rows]
                    or list(line.get_ydata()) != [float(row[spec["observable"]]) for row in rows]
                    or line.get_path().should_simplify):
                raise ValueError(f"raw vertices changed: {spec['plot_panel']} {spec['observable']}")
            vertices += len(rows)
        if len(axis.collections) != len(spec["endpoints"]):
            raise ValueError("unexpected endpoint artists")
        for artist, endpoint in zip(axis.collections, spec["endpoints"]):
            if artist.get_offsets().tolist() != [[endpoint["xi"], endpoint["value"]]]:
                raise ValueError("sampled branch marker moved away from its original vertex")
        if "viewport_clipping" not in spec:
            low, high = axis.get_ylim()
            if not low <= spec["data_min"] <= spec["data_max"] <= high:
                raise ValueError("full-range view clipped raw values")
    return vertices


def chart_requests():
    for name in COMPOSITES:
        yield Path("composites") / name, "composite", (name,)
    yield Path("details/figure1_low_value_details"), "detail", ()
    for mode in MODE_DIRS:
        for panel in PANELS[mode]:
            for observable in OBSERVABLES:
                yield Path(MODE_DIRS[mode]) / f"plot_panel={panel}" / f"{observable}_vs_xi", "single", (mode, panel, observable)


def render_request(kind, arguments, data, profile):
    return (render_single(*arguments, data, profile) if kind == "single" else
            render_composite(*arguments, data, profile) if kind == "composite" else render_detail(data, profile))


def preview(output_root):
    """Disposable layout preview; final review bundles use build_review."""
    output_root.mkdir(parents=True, exist_ok=False)
    data, profile = load_inputs(), load_profile("candidate_aps_v2")
    checks = []
    for relative, kind, args in list(chart_requests())[:3]:
        figure, specs, placements = render_request(kind, args, data, profile)
        try:
            vertices = verify_raw_vertices(figure, specs, data)
            figure.set_dpi(profile.dpi)
            outputs, quality = export_figure(figure, output_root / relative, profile, formats=["png"])
            check = dict(name=relative.as_posix(), raw_vertices=vertices, placements=placements,
                         glyph_mm=quality["minimum_capital_numeral_height_mm"],
                         **{key: quality[key] for key in LAYOUT_ERRORS})
            checks.append(check)
            print(json.dumps(check), flush=True)
        finally:
            plt.close(figure)
    write_manifest(output_root / "preview_checks.json", {"status": "layout_preview_only", "checks": checks})


def freeze_sources(results, profile):
    paths = [Path(__file__), profile.path,
             ROOT / ".agents/skills/plotting-sop/SKILL.md",
             ROOT / "docs/guides/sop/workflows/figure_production.md",
             ROOT / "docs/api/relaxtime/transport/TransportCoefficients.md",
             ROOT / "docs/analysis/relaxtime/phase_guided_transport/plotting_case_contract.md",
             ROOT / "tests/unit/python/test_phase_guided_full_rs_raw.py",
             *sorted((ROOT / "scripts/plotting").glob("*.py"))]
    records = []
    for path in paths:
        target = results / "source" / path.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        record = input_record(path, role="raw_renderer_or_plotting_contract")
        record["source_snapshot"] = input_record(target, role="immutable_source_copy")
        records.append(record)
    return records


def write_handoff(results, data):
    temperatures = ["| μB (MeV) | αT=1.0 | αT=1.1 | αT=1.2 |", "| ---: | ---: | ---: | ---: |"]
    for panel in PANELS["mode_a"]:
        values = [data.groups[("mode_a", panel, parameter)][0]["T_MeV"] for parameter in PARAMETERS["mode_a"]]
        temperatures.append("| " + panel[3:] + " | " + " | ".join(values) + " |")
    text = r"""# Full-RS 原始 PNG 图注与参数

本包从已验收的 1819 个原始点直接绘制：mode-A 910 点，mode-B 909 点。
12 个量均读取 CSV 同名列；不重新计算比值，不平滑、插值、拟合、删点或替换值。
470 个有限 τ 比值提示点全部保留。相邻原始点以直线相连；关闭路径顶点简化。
图像为本机纯 Python 绘制，没有调用求解器或补做本机数值验证。

图包与 v13 对应：两张 mode-A 主图、四面板低值细节图、72 张 A/B 单图；
另有逐张同尺寸灰度 PNG。只借用历史结构与排版，未读取旧 clean_value 点表、
平滑表、替换表或旧 gap 表。

Figure 1 建议图注：

> Raw full-RS relaxation times versus the anisotropy parameter ξ. Columns show
> μB=0, 450 and 900 MeV; rows show τu, τs, τū and τs̄. All vertical axes are
> logarithmic with independent ranges. Solid, dashed and dash-dotted curves
> denote αT=1.0, 1.1 and 1.2, with T=αT Tphase(μB, ξ=0).
> For μB=900 MeV and αT=1.0, the line is disconnected between the original
> samples at ξ=−0.003 and +0.003. Open circles and squares identify the
> nearest sampled points on the chirally restored and broken sides of the
> reviewed first-order discontinuity; they are not extrapolated coexistence endpoints.
> All original samples are retained.

Figure 2 建议图注：

> Raw full-RS transport ratios η/s, ζ/s and σ/T versus ξ, with the same
> columns and parameter encoding as Figure 1. Vertical axes are linear and
> use independent ranges. The first-order markers have the same limited
> scope as in Figure 1. Ratios are read directly from the accepted source CSV.

低值细节图建议图注：

> Linear low-value views of Figure 1 panels (a), (g), (h) and (i), using
> exactly the same raw vertices. Their vertical windows are 0.57–0.66,
> 0.57–0.66, 0.375–0.43 and 0.27–0.41 fm. The first two windows extend
> 0.01 fm below the v13 windows to retain the new αT=1.2 low-value branch.
> Curves outside these windows continue in the full-range main figure;
> leaving a viewport is not a physical endpoint.

单图：mode-A 的参数图例为 αT；mode-B 固定 T=120/160/200 MeV，
参数图例为 μB=0/450/900 MeV。mode-B 的 T=120 MeV、μB=900 MeV
只在 ξ=−0.13 与 −0.12 两个原始点之间断线，并保留两侧点及相态 marker。
其余相态字段变化（包括 crossover）不触发断线。单图保持全部原始 y 范围；
mode-A 弛豫时间统一 log-y，其余单图仅在正值动态范围至少 20 时选 log-y。
每图具体轴范围、输入行和分段记录在 manifest。

自然单位：τ 为 fm；η、ζ 为 fm⁻³；σ 为 fm⁻¹；η/s、ζ/s、σ/T 和 ξ 无量纲。
μB=3μq。源表 hadron/quark 对应此 PNJL 计算的手征破缺/恢复分支，
不意味着另算了一套强子气体输运。分段依据来自本次 quote-remote 背景审计。

温度映射（原始字符串，显示位数不表示不确定度）：

""" + "\n".join(temperatures) + r"""

版面：原生宽 171.45 mm（6.75 in），600 dpi，普通及数学字体来源逐图记录。
推荐原生宽度查看；逐图 placement_limits 给出最小插入宽度，不默认允许缩成单栏。
彩色与灰度各 75 张。PNG 作者接受、矢量 PDF、论文资格与当前图层晋升分别处理。
本次交付阶段为 png_review，manuscript_eligible=false，current_publication_layer=false。
图像可用于检查局部起伏；排版通过或相邻变化筛查通过不证明所有非物理尖峰已消除。
"""
    path = results / "caption_handoff.md"
    path.write_text(text, encoding="utf-8", newline="\n")
    return path


def write_chart(figure, relative, specs, placements, kind, data, profile, figures, results, inputs, generator, hash_cache):
    vertices = verify_raw_vertices(figure, specs, data)
    try:
        figure.set_dpi(profile.dpi)
        outputs, quality = export_figure(figure, figures / relative, profile, formats=["png"])
        verify_raw_vertices(figure, specs, data)
    finally:
        plt.close(figure)
    outputs[0]["role"] = "color_review"
    color = figures / relative.with_suffix(".png")
    gray = figures / "grayscale_review" / relative.with_suffix(".png")
    gray.parent.mkdir(parents=True, exist_ok=True)
    with Image.open(color) as source:
        ImageOps.grayscale(source.convert("RGB")).save(gray, dpi=(profile.dpi, profile.dpi))
    gray_record = output_record(gray, fmt="png", dpi=profile.dpi, vector=False)
    gray_record.update(role="grayscale_review", inspection=inspect_export(gray),
                       source_color_sha256=outputs[0]["sha256"], conversion="Pillow RGB to L; grayscale only")
    outputs.append(gray_record)
    mode = specs[0]["mode_key"]
    registry = make_registry(mode, profile)
    series = []
    for style in registry.styles:
        curves = [curve for spec in specs for curve in spec["series"] if curve["series_id"] == style.series_id]
        series.append({**style.record(), "state": "accepted_full_grid_raw",
                       "support_rule": "every original source row of each explicitly listed faceted curve",
                       "mask_rule": "no point masks; disconnect only reviewed first-order adjacent pairs",
                       "source_curves": curves, "row_count": sum(curve["row_count"] for curve in curves)})
    axes = []
    for index, spec in enumerate(specs):
        axes.extend([
            dict(field="xi", source_unit="dimensionless", display_unit="dimensionless", label=r"$\xi$",
                 transform="identity", scale="linear", panel_index=index),
            dict(field=spec["observable"], source_unit=UNITS[spec["observable"]], display_unit=UNITS[spec["observable"]],
                 label=LABELS[spec["observable"]], transform="identity", scale=spec["axis_scale"], panel_index=index),
        ])
    manifest = build_manifest(
        asset_id=f"relaxtime.phase_guided.{figures.name}.{relative.as_posix()}",
        figure_family="phase_guided_transport", case_slug=figures.name, figure_mode="audit",
        semantic_status="raw_full_rs_grid_png_review", style_profile=profile.profile_id,
        publication_scope="internal_review", generator=generator, inputs=inputs, axes=axes, series=series, outputs=outputs,
        selection_rule="all original rows for every listed mode/panel/parameter; ascending original xi",
        interpolation_policy="none", connector_policy="forbidden",
        missing_value_policy="fail on missing/nonfinite/duplicate data; preserve reviewed physical discontinuities",
        validation=dict(finite=True, duplicate_keys=True, support=True, strict_gate=False,
                        original_vertices_verified=True, original_vertex_count=vertices),
        calculation_sha=CALCULATION_SHA, source_run_id=data.acceptance["workflow_ref"],
        rendering=dict(
            column="double_column", figure_size_inches=quality["figure_size_inches"],
            size_override_reason="v13 homologous composite/detail canvas" if kind != "single" else None,
            typography_exception=None, quality=quality, placement_limits=placement_limits(quality, profile, outputs),
            legend_policy="declared_geometry_checked", legend_placements=placements,
            legend_outside=all(item["placement"] == "outside_axes" for item in placements),
            panel_specs=specs, color_route="undecided_review",
            series_styles=registry.record([style.series_id for style in registry.styles]),
            series_identity_scope="parameter identity shared across facets; source_curves record full panel/observable identity",
            review_contract="prd_review_v1",
            accessibility=accessibility_record([
                dict(value=style.series_id, color=style.color, noncolor_encoding=f"line style {style.linestyle}")
                for style in registry.styles], backgrounds=["#ffffff"],
                background_scope="opaque white axes; no fills or shaded regions"),
            grayscale_review=dict(status="author_review_required", output=gray_record["path"],
                                  machine_checks="same dimensions and pixel-exact RGB-to-L conversion"),
            path_simplification=False,
            drawing_rule="straight line segments between consecutive original samples; no added support vertices",
            caption_handoff=(results / "caption_handoff.md").relative_to(ROOT).as_posix(),
            delivery_stage="png_review", vector_delivery_pending=True, output_formats=["png"],
        ),
    )
    manifest.update(manuscript_eligible=False, current_publication_layer=False, solver_called=False,
                    canonical_data_modified=False, new_display_values=False, delivery_stage="png_review",
                    vector_delivery_pending=True, formal_production_promoted=data.acceptance["formal_production_promoted"],
                    numerical_status="full_grid_numerical_acceptance_passed")
    violations = validate_manifest_record(manifest, hash_cache=hash_cache)
    if violations:
        raise ValueError(f"{relative}: " + "; ".join(violations))
    return dict(stem=(figures / relative).relative_to(ROOT).as_posix(), kind=kind,
                mode_key=mode, outputs=outputs, manifest_record=manifest)


def build_review(case=DEFAULT_CASE, *, reference_root=None):
    if not case or Path(case).name != case or any(character in case for character in ("/", "\\", ":")):
        raise ValueError("case must be a single new directory name")
    figures, results = FIGURE_PARENT / case, RESULT_PARENT / case
    for path in (figures, results):
        if path.exists():
            raise FileExistsError(f"refusing to overwrite existing case: {path}")
    data, profile = load_inputs(), load_profile("candidate_aps_v2")
    expected = {relative.with_suffix(".png").as_posix() for relative, _, _ in chart_requests()}
    reference_root = reference_root or Path(
        "D:/Desktop/Julia_RelaxTime/data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v13_png_review")
    actual_reference = {path.relative_to(reference_root).as_posix() for path in reference_root.rglob("*.png")
                        if "grayscale_review" not in path.relative_to(reference_root).parts}
    if actual_reference != expected:
        raise ValueError("reference v13 structure differs from the 75 requested color views")
    figures.mkdir(parents=True)
    results.mkdir(parents=True)
    inputs = [input_record(path, role="accepted_raw_grid" if path in data.sources.values() else "numerical_acceptance_evidence")
              for path in data.evidence_paths]
    inputs.extend(freeze_sources(results, profile))
    reference_record = input_record(reference_root / "plot_manifest.json", role="v13_layout_reference_only")
    inputs.append(reference_record)
    caption = write_handoff(results, data)
    inputs.append(input_record(caption, role="raw_figure_caption_and_parameter_mapping"))
    generator = generator_record(Path(__file__),
                                 command=f"python scripts/analysis/relaxtime/build_phase_guided_full_rs_raw.py --png-review --case {case}",
                                 source_snapshot=results / "source" / Path(__file__).relative_to(ROOT),
                                 runtime=runtime_record({"matplotlib": matplotlib.__version__, "Pillow": Image.__version__}))
    charts, hash_cache = [], {}
    try:
        for relative, kind, arguments in chart_requests():
            figure, specs, placements = render_request(kind, arguments, data, profile)
            charts.append(write_chart(figure, relative, specs, placements, kind, data, profile,
                                      figures, results, inputs, generator, hash_cache))
            print(f"[{len(charts):02d}/75] {relative.as_posix()}", flush=True)
        # Re-read hashes after rendering rather than relying only on the validator cache.
        for record in inputs:
            path = Path(record["path"])
            if not path.is_absolute():
                path = ROOT / path
            if sha256_file(path) != record["sha256"]:
                raise ValueError(f"protected input changed while rendering: {path}")
        records = [chart["manifest_record"] for chart in charts]
        package = dict(case_slug=case, single_figure_count=72, composite_figure_count=2, detail_figure_count=1,
                       color_figure_count=75, grayscale_figure_count=75, charts=charts,
                       delivery_stage="png_review", manuscript_eligible=False, current_publication_layer=False,
                       vector_delivery_pending=True, author_acceptance=None,
                       layout_reference=reference_record, raw_input_point_count=1819,
                       numerical_data_postprocessing="none", source_inputs_verified_unchanged=True)
        manifest_path = figures / "plot_manifest.json"
        write_manifest(manifest_path, build_bundle(package, records))
        errors = validate_manifest(manifest_path)
        if errors:
            raise ValueError("; ".join(errors))
        color_names = {path.relative_to(figures).as_posix() for path in figures.rglob("*.png")
                       if "grayscale_review" not in path.relative_to(figures).parts}
        if color_names != expected or len(list(figures.rglob("*.png"))) != 150:
            raise ValueError("final PNG inventory does not match the requested 75 color/gray pairs")
        for chart in charts:
            color, gray = (ROOT / record["path"] for record in chart["outputs"])
            with Image.open(color) as a, Image.open(gray) as b:
                if a.size != b.size or b.mode != "L" or ImageOps.grayscale(a.convert("RGB")).tobytes() != b.tobytes():
                    raise ValueError(f"grayscale conversion mismatch: {gray}")
        quality = [record["rendering"]["quality"] for record in records]
        report = dict(status="passed", figure_count=75, png_count=150, raw_input_points=1819,
                      raw_observables=len(OBSERVABLES), original_vertices_verified=sum(
                          record["validation"]["original_vertex_count"] for record in records),
                      original_data_values_changed=0, removed_source_points=0, new_support_points=0,
                      retained_quality_flags=data.acceptance["quality_flags_retained"],
                      physical_discontinuities=[dict(mode=m, panel=p, parameter=s, pairs=v)
                                                for (m, p, s), v in data.gaps.items()],
                      minimum_glyph_height_mm=min(q["minimum_capital_numeral_height_mm"] for q in quality),
                      minimum_curve_linewidth_pt=min(q["minimum_curve_linewidth_pt"] for q in quality),
                      geometry_issues={key: sum(len(q[key]) if isinstance(q[key], list) else q[key] for q in quality)
                                       for key in LAYOUT_ERRORS},
                      minimum_required_width_inches=max(
                          record["rendering"]["placement_limits"].get("minimum_width_inches", 0) for record in records),
                      grayscale_pairs_pixel_exact=True, source_inputs_verified_unchanged=True,
                      v13_relative_color_paths_identical=True, validation_errors=errors,
                      manifest=input_record(manifest_path, role="validated_plot_bundle"), author_acceptance=None)
        write_manifest(results / "quality_review.json", report)
        figure_path = figures.relative_to(ROOT).as_posix()
        (results / "README.md").write_text(
            "# Full-RS 原始 PNG 审阅包\n\n"
            f"已生成 75 张彩色及 75 张灰度 PNG，600 dpi。图包：{figure_path}\n\n"
            "原始 CSV 的 1819 点、12 个量全部保留；所有绘图顶点逐一与源表核对。"
            "未平滑、插值、删点或替换值，仅在两处已审计的一阶不连续之间断线。\n\n"
            "两张主图、四面板低值细节图及 72 张单图的相对路径与 v13 一一对应。"
            "前两张低值子图的下界由 0.58 改为 0.57 fm，以显示新的 αT=1.2 支。\n\n"
            "图注、单位及精确温度见 [caption_handoff.md](caption_handoff.md)；"
            "机器验收见 [quality_review.json](quality_review.json)。"
            "源码与绘图合同的原字节保存在 source/，计算提交为 "
            f"{CALCULATION_SHA}。\n\n"
            "本包停在用户要求的 PNG 视觉审阅阶段。没有调用 Julia，"
            "没有生成 PDF、修改论文或替换当前图层。作者接受尚待记录。\n",
            encoding="utf-8", newline="\n")
        print(f"Validated 75 raw color/gray pairs: {figures}", flush=True)
        return manifest_path
    except Exception as exc:
        write_manifest(results / "failure.json", {"status": "failed", "error": str(exc), "completed_figures": len(charts)})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--png-review", action="store_true")
    action.add_argument("--preview", type=Path)
    parser.add_argument("--case", default=DEFAULT_CASE)
    parser.add_argument("--reference-root", type=Path)
    args = parser.parse_args()
    if args.preview:
        preview(args.preview)
    else:
        build_review(args.case, reference_root=args.reference_root)


if __name__ == "__main__":
    main()
