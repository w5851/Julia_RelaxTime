"""Render compact transport panels from validated frozen scan rows.

The input adapter preserves source row identity and explicitly declared
first-order gaps. It does not load a historical figure package.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from decimal import Decimal
from itertools import product
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import (
    AutoLocator, FixedLocator, FormatStrFormatter, FuncFormatter, LogLocator,
    MaxNLocator, NullFormatter, ScalarFormatter,
)

from scripts.plotting.plot_quality import measure_figure, visible_legends
from scripts.plotting.plot_series import SeriesStyle, SeriesStyleRegistry
from scripts.plotting.plot_style import configure_axis_ticks
from scripts.plotting.publication_review import LABEL_PT, LEGEND_PT, TICK_PT, WIDTH_IN, figure_context


ROOT = Path(__file__).resolve().parents[2]


MODE_DIRS = {"mode_a": "mode_a_fixed_muB_phase_scaled", "mode_b": "mode_b_fixed_T_sparse_muB"}

PANELS = {"mode_a": ("muB0.0", "muB450.0", "muB900.0"), "mode_b": ("T120.0", "T160.0", "T200.0")}

PARAMETERS = {"mode_a": ("alpha1.0", "alpha1.1", "alpha1.2"), "mode_b": ("muB0.0", "muB450.0", "muB900.0")}

OBSERVABLES = ("tau_u", "tau_d", "tau_s", "tau_ubar", "tau_dbar", "tau_sbar",
               "eta", "zeta", "sigma", "eta_over_s", "zeta_over_s", "sigma_over_T")

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

ENDPOINT_MARKERS = {"quark": "o", "hadron": "s"}

LAYOUT_ERRORS = ("clipped_text", "text_overlap_pairs", "legend_curve_overlap_count",
                 "legend_landmark_overlap_count", "legend_in_axes_overflow_count", "legend_pair_overlap_count")

@dataclass
class TransportInputs:
    groups: dict
    gaps: dict
    sources: dict

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
    return SeriesStyleRegistry(f"transport_review_{mode}_parameters_v1", tuple(
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
    """Keep ordinary 1/2/5-style labels readable in compact physical panels."""
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

def render_panel(axis, mode, panel, observable, data, profile, *, composite=False):
    registry = make_registry(mode, profile)
    curves, endpoints, values = [], [], []
    for parameter, style in zip(PARAMETERS[mode], registry.styles):
        key = (mode, panel, parameter)
        rows, gaps = data.groups[key], data.gaps.get(key, [])
        segments = split_segments(rows, gaps)
        for segment in segments:
            line, = axis.plot([float(r["xi"]) for r in segment], [float(r[observable]) for r in segment],
                              **style.plot_kwargs())
            line.set_gid(style.series_id)
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
    axis.set_yscale(scale)
    axis.set_xlim(-0.52, 0.52)
    axis.set_xticks([-0.5, 0, 0.5], ["-0.5", "0", "0.5"])
    axis.margins(y=0.06)
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
    """Place compact keys using deterministic, measured corner selection."""
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

def parameter_legend(figure, specs, profile, *, host, location, anchor=None, size=11):
    mode = specs[0]["mode_key"]
    registry = make_registry(mode, profile)
    labels = ["1.0", "1.1", "1.2"] if mode == "mode_a" else ["0", "450", "900"]
    title = r"$\alpha_T$" if mode == "mode_a" else r"$\mu_B\;(\mathrm{MeV})$"
    handles = [Line2D([], [], color=style.color, ls=style.linestyle, lw=1.2, label=label)
               for style, label in zip(registry.styles, labels)]
    settings = dict(loc=location, frameon=False, alignment="left", fontsize=size, title_fontsize=13,
                    handlelength=1.35, handletextpad=0.3, labelspacing=0.10, borderpad=0.05, borderaxespad=0.25)
    if anchor is not None:
        settings.update(bbox_to_anchor=anchor, borderaxespad=0, borderpad=0.25, labelspacing=0.20,
                        handlelength=1.6, handletextpad=0.4)
    figure.axes[host].legend(handles=handles, title=title, **settings)
    return dict(role="parameters", placement="in_axes", host_axes_index=host, location=location,
                title=title, labels=labels, scope="all parameter curves", series_keys=list(PARAMETERS[mode]),
                legend_index=0, anchor_axes=anchor, selection_reason="reviewed fixed-temperature panel layout")

def render_fixed_T(data, profile):
    with figure_context(profile) as context:
        figure, axes = plt.subplots(3, 3, figsize=(6.75, 7.4), sharex=True)
        context.bind(figure)
        figure.subplots_adjust(left=0.12, right=0.98, bottom=0.09, top=0.95, hspace=0.26, wspace=0.38)
        specs = []
        for row, observable in enumerate(("eta_over_s", "zeta_over_s", "sigma_over_T")):
            for col, panel in enumerate(PANELS["mode_b"]):
                axis = axes[row, col]
                spec = render_panel(axis, "mode_b", panel, observable, data, profile, composite=True)
                scale = "log" if row < 2 else "linear"
                axis.set_yscale(scale)
                axis.relim()
                axis.autoscale_view(scalex=False)
                axis.margins(y=0.05)
                if scale == "log":
                    low, high = axis.get_ylim()
                    subs = (1.0,) if high / low > 10 ** 1.5 else (1.0, 2.0, 5.0)
                    axis.yaxis.set_major_locator(LogLocator(base=10, subs=subs))
                    axis.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
                    axis.yaxis.set_minor_formatter(NullFormatter())
                    visible = [value for value in axis.get_yticks() if low <= value <= high]
                    if len(visible) > 5 or len(visible) < 2:
                        # LogLocator may fall back to dense linear ticks within
                        # a narrow decade. Keep the logarithmic axis and space
                        # its ordinary-number labels at the existing font size.
                        axis.yaxis.set_major_locator(FixedLocator(readable_log_ticks(low, high)))
                else:
                    axis.yaxis.set_major_locator(AutoLocator())
                    axis.yaxis.set_major_formatter(ScalarFormatter())
                configure_axis_ticks(axis, profile)
                if col:
                    axis.set_ylabel("")
                if row == 2:
                    axis.set_xlabel(r"$\xi$")
                label = f"({chr(97 + row * 3 + col)})"
                if row == 0:
                    label += rf"  $T={float(panel[1:]):.0f}$ MeV"
                axis.text(0, 1.03, label, transform=axis.transAxes, ha="left", va="bottom", fontsize=13)
                for collection in axis.collections:
                    collection.set_sizes([5.5 ** 2])
                    collection.set_linewidths([1.0])
                spec.update(axis_scale=scale, axis_scale_reason="paper_FIG4_fixed_T_row_scale")
                capture_ticks(axis, spec)
                specs.append(spec)
        placements = [parameter_legend(figure, specs, profile, host=6, location="lower left",
                                      anchor=(0.4466411003595683, 0.31293851991526406), size=13)]
        endpoints = [Line2D([], [], ls="None", marker=marker, mfc="white", mec="black", ms=6,
                            label=label) for marker, label in (("o", "restored"), ("s", "broken"))]
        axes[0, 0].legend(handles=endpoints, title="First-order", loc="lower left",
                          bbox_to_anchor=(0.16557086802564835, 0.3187792443606394), borderaxespad=0,
                          borderpad=0.25, labelspacing=0.20, handlelength=1.6, handletextpad=0.4,
                          frameon=False, fontsize=13, title_fontsize=13)
        placements.append(dict(role="endpoints", placement="in_axes", host_axes_index=0,
                               location="lower left", anchor_axes=[0.16557086802564835, 0.3187792443606394],
                               title="First-order", labels=["restored", "broken"],
                               scope="T=120 MeV, muB=900 MeV sampled first-order branch pair only",
                               applies_to=[dict(mode_key="mode_b", plot_panel="T120.0", plot_series="muB900.0")],
                               legend_index=0, selection_reason="reviewed fixed-temperature panel layout"))
        # visible_legends orders by host axes, not by creation order.
        for item in placements:
            item["legend_index"] = 0 if item["host_axes_index"] == 0 else 1
        placements.sort(key=lambda item: item["legend_index"])
        quality = measure_figure(figure, intended_width_inches=6.75)
        if not clear_layout(quality):
            raise ValueError("fixed-temperature layout obstructs data: " + json.dumps(
                {key: quality[key] for key in LAYOUT_ERRORS if quality[key]}))
        return figure, specs, placements
