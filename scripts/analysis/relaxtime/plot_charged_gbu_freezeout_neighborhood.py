#!/usr/bin/env python3
"""Solver-free freeze-out guidance map from immutable screening CSV/manifest.

Only the displayed neighborhood and contour levels change. No fit, solver,
zero-filling, smoothing, gap bridging, or production-default change occurs.
Legacy review mode retains +/-10 MeV parameter probes. Presentation mode shows
only the current freeze-out curve and explicit 0.05 contour spacing.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tomllib


ROOT = Path(__file__).resolve().parents[3]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def freezeout_T(muB_MeV: float, coefficients: dict, shift_MeV: float = 0) -> float:
    mu = muB_MeV / 1000
    return 1000 * (coefficients["a_GeV"] - coefficients["b_GeV_inv1"] * mu**2
                   - coefficients["c_GeV_inv3"] * mu**4) + shift_MeV


def energy_point(energy_GeV: float, coefficients: dict) -> tuple[float, float]:
    muB = 1000 * coefficients["d_GeV"] / (1 + coefficients["e_GeV_inv1"] * energy_GeV)
    return muB, freezeout_T(muB, coefficients)


def load_frozen(csv_path: Path, manifest_path: Path, profile_path: Path) -> tuple[list[dict], dict, dict]:
    source = json.loads(manifest_path.read_text(encoding="utf-8"))
    if digest(csv_path) != source["merged_csv"]["sha256"]:
        raise ValueError("CSV hash differs from frozen source manifest")
    if digest(profile_path) != source["reference_lines"]["freezeout"]["source_sha256"]:
        raise ValueError("freeze-out profile differs from frozen source manifest")
    rows = []
    keys = set()
    with csv_path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            key = (float(raw["T_MeV"]), float(raw["muB_MeV"]))
            if not all(math.isfinite(value) for value in key) or key in keys:
                raise ValueError("nonfinite or duplicate grid key")
            keys.add(key)
            row = {"T": key[0], "muB": key[1], "status": raw["status"], "R": None}
            if raw["status"] == "screened":
                pi, kaon, ratio = (float(raw[name]) for name in
                                  ("pi_plus_density_inv_fm3", "K_plus_density_inv_fm3", "Kplus_over_pi_plus"))
                if not all(math.isfinite(value) for value in (pi, kaon, ratio)) or pi <= 0:
                    raise ValueError("invalid screened density/ratio")
                if not math.isclose(kaon / pi, ratio, rel_tol=1e-12, abs_tol=1e-14):
                    raise ValueError("ratio is inconsistent with the frozen densities")
                if any(raw[name].lower() != "true" for name in ("pi_plus_passed", "K_plus_passed")):
                    raise ValueError("screened row lacks plus-channel support")
                row["R"] = ratio  # Retain signed values; never clip or take abs.
            rows.append(row)
    expected = {(float(T), float(mu)) for T in source["grid"]["T_MeV"]
                for mu in source["grid"]["muB_MeV"]}
    if len(rows) != source["row_count"] or keys != expected:
        raise ValueError("row count or rectangular grid differs from frozen source")
    actual_counts = {status: sum(row["status"] == status for row in rows)
                     for status in {row["status"] for row in rows}}
    if actual_counts != source["status_counts"]:
        raise ValueError("status counts differ from frozen source")
    with profile_path.open("rb") as handle:
        coefficients = tomllib.load(handle)["coefficients"]
    return rows, source, coefficients


def neighborhood(rows: list[dict], coefficients: dict, *, band_MeV: float = 20,
                 muB_max_MeV: float = 750,
                 T_bounds_MeV: tuple[float, float] = (50, 190)) -> dict:
    import numpy as np

    Ts = sorted({row["T"] for row in rows if T_bounds_MeV[0] <= row["T"] <= T_bounds_MeV[1]})
    mus = sorted({row["muB"] for row in rows if 0 <= row["muB"] <= muB_max_MeV})
    lookup = {(row["T"], row["muB"]): row for row in rows}
    values = np.full((len(Ts), len(mus)), np.nan)
    selected = np.zeros(values.shape, dtype=bool)
    failed = np.zeros(values.shape, dtype=bool)
    for i, T in enumerate(Ts):
        for j, mu in enumerate(mus):
            selected[i, j] = abs(T - freezeout_T(mu, coefficients)) <= band_MeV
            row = lookup[(T, mu)]
            if selected[i, j]:
                failed[i, j] = row["R"] is None
                if row["R"] is not None:
                    values[i, j] = row["R"]
    valid = np.isfinite(values)
    cells = valid[:-1, :-1] & valid[1:, :-1] & valid[:-1, 1:] & valid[1:, 1:]
    return {"T": np.array(Ts), "muB": np.array(mus), "R": values,
            "selected": selected, "failed": failed, "valid_cells": cells}


def uniform_levels(minimum: float, maximum: float, step: float = .05) -> list[float]:
    """Exact decimal-friendly contour levels; no normalization of physical values."""
    return [round(index * step, 10) for index in
            range(math.ceil(minimum / step), math.floor(maximum / step) + 1)]


def is_labeled_level(level: float) -> bool:
    return math.isclose(level * 10, round(level * 10), rel_tol=0, abs_tol=1e-8)


def density_route(source: dict) -> str:
    route = source.get("density_route", "direct_finite_q")
    if route not in {"direct_finite_q", "q0_lambda_reference"}:
        raise ValueError(f"unsupported density_route: {route}")
    return route


def plot_presentation(data: dict, coefficients: dict, band_MeV: float, *,
                      route: str = "direct_finite_q", color_limits=None,
                      clean_contours: bool = False, prescription: str | None = None) -> tuple:
    """Large group-meeting canvas, not an APS/final-delivery qualification."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize, ListedColormap
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.ticker import AutoMinorLocator, MultipleLocator
    import matplotlib.patheffects as pe
    import numpy as np

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 17,
                         "mathtext.fontset": "dejavusans", "axes.labelsize": 20,
                         "xtick.labelsize": 16, "ytick.labelsize": 16})
    fig = plt.figure(figsize=(12.2, 8.8), dpi=120)
    ax = fig.add_axes([.09, .19, .75, .70])
    cax = fig.add_axes([.88, .19, .025, .70])
    values = np.ma.masked_invalid(data["R"])
    minimum, maximum = float(values.min()), float(values.max())
    lower = math.floor(minimum / .05) * .05
    upper = math.ceil(maximum / .05) * .05
    if color_limits is not None:
        lower, upper = color_limits
        if not all(math.isfinite(x) for x in color_limits) or lower >= upper:
            raise ValueError("ordered finite comparison color limits required")
    extend = ("both" if minimum < lower and maximum > upper else
              "min" if minimum < lower else "max" if maximum > upper else "neither")
    cmap = plt.get_cmap("YlGnBu").copy()
    cmap.set_bad("#eeeeee")
    cmap.set_under("#c58bbd")
    cmap.set_over("#10234c")
    ax.set_facecolor("#eeeeee")
    mesh = ax.pcolormesh(data["muB"], data["T"], values, shading="nearest",
                         cmap=cmap, norm=Normalize(lower, upper, clip=False), alpha=.75)
    if not clean_contours:
        ax.pcolormesh(data["muB"], data["T"],
                      np.ma.masked_where(~data["failed"], np.ones(values.shape)),
                      shading="nearest", cmap=ListedColormap(["#777777"]), vmin=0, vmax=1)
    colorbar = fig.colorbar(mesh, cax=cax, extend=extend)
    colorbar.set_label(r"$R_+=n_{K^+}/n_{\pi^+}$", fontsize=20, labelpad=12)
    colorbar.set_ticks(uniform_levels(lower, upper, .10))
    colorbar.ax.yaxis.set_minor_locator(MultipleLocator(.05))
    colorbar.ax.yaxis.set_major_formatter(plt.FormatStrFormatter("%.2f"))
    colorbar.ax.tick_params(labelsize=13)

    levels = uniform_levels(minimum, maximum)
    contour = ax.contour(data["muB"], data["T"], values, levels=levels,
                         colors="black", linewidths=1.15, corner_mask=False, zorder=5)
    mus = np.linspace(0, 750, 401)
    default_T = np.array([freezeout_T(mu, coefficients) for mu in mus])
    line, = ax.plot(mus, default_T, color="#c00000", linewidth=3,
                    linestyle="--", label="Current chemical freeze-out", zorder=7)
    for shift in (-band_MeV, band_MeV):
        ax.plot(mus, default_T + shift, color="#888888", linewidth=1,
                 linestyle=":", zorder=4)
    ax.set_xlim(0, 750)
    ax.set_ylim(float(data["T"].min()), float(data["T"].max()))
    ax.set_xticks([0, 150, 300, 450, 600, 750])
    ax.set_yticks([40, 70, 100, 130, 160, 190, 220])
    ax.set_xlabel(r"$\mu_B$ (MeV)")
    ax.set_ylabel(r"$T$ (MeV)")
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    fig.text(.47, .955, r"$K^+/\pi^+$ near chemical freeze-out", ha="center", fontsize=23)
    route_label = "q=0 extrapolated" if route == "q0_lambda_reference" else "finite-$q$"
    if route == "q0_lambda_reference" and prescription == "finite_window":
        route_label += "; finite window"
    fig.text(.47, .910, rf"Quark-only BQS | {route_label} GBU screening | $T_{{\rm fo}}\pm{band_MeV:g}$ MeV",
             ha="center", fontsize=16)

    # Select well-separated on-path positions for one single clabel operation.
    # This avoids the repeated subset-clabel nearest-level mismatch of MPL 3.10.
    fig.canvas.draw()
    positions = []
    expected_labels = []
    accepted_pixels = []
    for index, (level, paths) in enumerate(zip(contour.levels, contour.allsegs)):
        if not is_labeled_level(float(level)):
            continue
        paths = sorted((path for path in paths if len(path) >= 4), key=len, reverse=True)
        candidates = []
        for path in paths:
            for fraction in np.linspace(.1, .9, 33):
                x, y = path[int(fraction * (len(path) - 1))]
                pixel = ax.transData.transform((x, y))
                if (25 < x < 715 and 45 < y < 214
                        and abs(y - freezeout_T(float(x), coefficients)) > 3):
                    distance = min((np.linalg.norm(pixel - other) for other in accepted_pixels), default=1e6)
                    preference = abs(fraction - (.25 if index % 2 == 0 else .70))
                    candidates.append((distance - 35 * preference, (float(x), float(y)), pixel))
        if candidates:
            _, position, pixel = max(candidates, key=lambda item: item[0])
            positions.append(position)
            expected_labels.append(f"{level:.2f}")
            accepted_pixels.append(pixel)
    # Repeat key labels on distant parts of long, nonmonotonic contours so the
    # upper/lower branches need not be traced across the entire slide by eye.
    for level, paths in zip(contour.levels, contour.allsegs):
        if round(float(level), 2) not in {.40, .50, .60}:
            continue
        primary_pixels = [pixel for pixel, value in zip(accepted_pixels, expected_labels)
                          if value == f"{level:.2f}"]
        candidates = []
        for path in paths:
            for fraction in np.linspace(.08, .92, 45):
                if len(path) < 8:
                    continue
                x, y = path[int(fraction * (len(path) - 1))]
                pixel = ax.transData.transform((x, y))
                if (25 < x < 715 and 45 < y < 214
                        and abs(y - freezeout_T(float(x), coefficients)) > 3
                        and all(np.linalg.norm(pixel - other) > 180 for other in primary_pixels)):
                    distance = min((np.linalg.norm(pixel - other) for other in accepted_pixels), default=1e6)
                    if distance > 70:
                        candidates.append((distance, (float(x), float(y)), pixel))
        if candidates:
            _, position, pixel = max(candidates, key=lambda item: item[0])
            positions.append(position)
            expected_labels.append(f"{level:.2f}")
            accepted_pixels.append(pixel)
    labels = ax.clabel(contour, manual=positions, inline=True, inline_spacing=6,
                        fmt="%.2f", fontsize=15, rightside_up=True)
    if [label.get_text() for label in labels] != expected_labels:
        raise ValueError("contour labels do not match their selected source-level paths")
    for label in labels:
        label.set_path_effects([pe.Stroke(linewidth=.7, foreground="white"), pe.Normal()])
        label.set_zorder(10)
    handles = [Line2D([], [], color="#c00000", linestyle="--", linewidth=3,
                     label="Current chemical freeze-out")]
    if not clean_contours:
        handles.append(Patch(facecolor="#777777", label="Failed support"))
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.47, .067),
               ncol=len(handles), fontsize=16)
    fig.text(.47, .043, r"Contour spacing: $\Delta R_+=0.05$; label spacing: $0.10$.",
             ha="center", fontsize=16)
    fig.text(.47, .012, r"Band clipped to $T\geq40$ MeV input; no fit.",
             ha="center", fontsize=14)
    return fig, {"minimum": minimum, "maximum": maximum, "color_limits": [lower, upper],
                 "colorbar_extend": extend, "density_route": route,
                 "density_prescription": prescription, "warning_annotations": False,
                 "failure_overlay": not clean_contours,
                 "under_range_nodes": int((data["R"] < lower).sum()),
                 "over_range_nodes": int((data["R"] > upper).sum()),
                 "contour_levels": levels, "contour_step": .05,
                 "label_step": .10, "label_white_background": False,
                 "label_outline_width_pt": .7,
                 "visible_contour_labels": [label.get_text() for label in labels],
                 "label_policy": "single inline-clabel call; every label verified against selected exact-level path",
                 "candidate_curves": False, "energy_labels": False}


def run_presentation(args, rows: list[dict], source: dict, coefficients: dict) -> int:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    paths = [args.csv, args.source_manifest, args.freezeout_profile, Path(__file__)]
    before = [digest(path) for path in paths]
    data = neighborhood(rows, coefficients, band_MeV=args.band_MeV, T_bounds_MeV=(40, 220))
    route = density_route(source)
    fig, display = plot_presentation(data, coefficients, args.band_MeV,
                                     route=route, color_limits=args.comparison_color_limits,
                                     clean_contours=args.clean_contours,
                                     prescription=source.get("density_prescription"))
    args.output_dir.mkdir(parents=True)
    png = args.output_dir / "freezeout_neighborhood_ratio.png"
    fig.savefig(png, dpi=220, bbox_inches=None)
    plt.close(fig)
    if [digest(path) for path in paths] != before:
        raise ValueError("frozen input changed during presentation rendering")
    manifest = {"schema_version": "charged_gbu_freezeout_guidance_presentation_v1",
                "delivery_stage": "group_meeting_review", "diagnostic_only": True,
                "manuscript_eligible": False, "current_publication_layer": False,
                "formal_plotting_sop_applied": False, "numerical_rescan": False,
                "source_run_id": source["source_run_id"], "calculation_sha": source["source_git_sha"],
                "inputs": [{"path": str(path.resolve()), "sha256": hash_value,
                            "bytes": path.stat().st_size} for path, hash_value in zip(paths, before)],
                "generator_command": subprocess.list2cmdline([sys.executable, *sys.argv]),
                "runtime": {"python": sys.version, "matplotlib": matplotlib.__version__},
                "background_contract": source["background_contract"],
                "density_route": route, "coordinate_contract": source.get("coordinate_contract"),
                "density_prescription": source.get("density_prescription"),
                "warning_display": "none" if args.clean_contours else "legacy failure support",
                "selection_rule": f"0<=muB<=750 MeV; |T-Tfo(muB)|<={args.band_MeV:g} MeV; 40<=T<=220 MeV",
                "interpolation_policy": "heatmap none; display-only contours within four valid corner cells; corner_mask=false",
                "missing_value_policy": "mask; no zero-fill or cross-gap contours",
                "normalization": "linear; explicit comparison range or native min/max rounded outward; overflow colors; no numerical clipping",
                "selected_nodes": int(data["selected"].sum()), "screened_nodes": int((data["selected"] & ~data["failed"]).sum()),
                "failed_nodes": int(data["failed"].sum()), "valid_cells": int(data["valid_cells"].sum()),
                "band_MeV": args.band_MeV, "rendering": display,
                "output": {"path": str(png.resolve()), "sha256": digest(png), "bytes": png.stat().st_size,
                           "dpi": 220, "figure_size_inches": [12.2, 8.8]}}
    (args.output_dir / "plot_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"png": str(png), "screened": manifest["screened_nodes"], "failed": manifest["failed_nodes"],
                      "range": [display["minimum"], display["maximum"]], "labels": display["visible_contour_labels"]}))
    return 0


def plot_map(data: dict, coefficients: dict, style) -> tuple:
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.ticker import FixedLocator, FormatStrFormatter
    import matplotlib.patheffects as pe
    import numpy as np
    from scripts.plotting.plot_style import configure_axis_ticks

    fig = plt.figure(figsize=(6.75, 6.2), dpi=150)
    ax = fig.add_axes([0.12, 0.32, 0.705, 0.575])
    cax = fig.add_axes([0.86, 0.32, 0.03, 0.575])
    ax.set_facecolor("#f2f2f2")
    values = np.ma.masked_invalid(data["R"])
    minimum, maximum = float(values.min()), float(values.max())
    lower, upper = math.floor(minimum / .05) * .05, math.ceil(maximum / .05) * .05
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#f2f2f2")
    mesh = ax.pcolormesh(data["muB"], data["T"], values, shading="nearest",
                         cmap=cmap, norm=Normalize(lower, upper), rasterized=False)
    # A separate mask layer distinguishes failed support from the unselected domain.
    from matplotlib.colors import ListedColormap
    ax.pcolormesh(data["muB"], data["T"], np.ma.masked_where(~data["failed"], np.ones(values.shape)),
                  shading="nearest", cmap=ListedColormap(["#555555"]), vmin=0, vmax=1)
    colorbar = fig.colorbar(mesh, cax=cax)
    colorbar.set_label(r"$R_+=n_{K^+}/n_{\pi^+}$", labelpad=8)
    ticks = np.arange(math.ceil(lower / .1) * .1, upper + .001, .1)
    cax.yaxis.set_major_locator(FixedLocator(ticks))
    cax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    cax.set_xticks([0, 1])
    cax.tick_params(axis="x", labelbottom=False, labeltop=False)
    configure_axis_ticks(cax, style)

    levels = np.arange(math.ceil(minimum / .05) * .05, maximum, .05)
    contour = ax.contour(data["muB"], data["T"], values, levels=levels,
                         corner_mask=False, colors="#eeeeee", linewidths=.75)
    mus = np.linspace(0, 750, 401)
    default_T = np.array([freezeout_T(mu, coefficients) for mu in mus])
    handles = []
    for shift, color, linestyle, label in (
        (0, "#111111", "-", "Default freeze-out"),
        (-10, "#d55e00", "--", r"$a-10$ MeV (probe)"),
        (10, "#0072b2", "-.", r"$a+10$ MeV (probe)"),
    ):
        line, = ax.plot(mus, default_T + shift, color=color, linestyle=linestyle,
                        linewidth=1.65, label=label, zorder=6)
        line.set_path_effects([pe.Stroke(linewidth=3.1, foreground="white"), pe.Normal()])
        handles.append(Line2D([], [], color=color, linestyle=linestyle, linewidth=1.65, label=label))
    for shift in (-20, 20):
        ax.plot(mus, default_T + shift, color="#777777", linewidth=.7, linestyle=":", zorder=4)

    # Labels identify energy mapping on the default curve, not fitted data.
    for energy, offset in ((200, (6, 12)), (20, (2, 12)), (7.7, (-12, -25)),
                           (5, (5, 12)), (3, (-8, -24))):
        mu, T = energy_point(energy, coefficients)
        ax.plot([mu], [T], "o", color="black", markerfacecolor="white", markersize=4.5, zorder=8)
        text = ax.annotate(f"{energy:g}", (mu, T), xytext=offset, textcoords="offset points",
                           ha="left", va="center", fontsize=13, color="black", zorder=9)
        text.set_path_effects([pe.Stroke(linewidth=2.5, foreground="white"), pe.Normal()])

    ax.set_xlim(0, 750)
    ax.set_ylim(50, 190)
    ax.set_xticks([0, 150, 300, 450, 600, 750])
    ax.set_yticks([60, 90, 120, 150, 180])
    ax.set_xlabel(r"$\mu_B\;(\mathrm{MeV})$")
    ax.set_ylabel(r"$T\;(\mathrm{MeV})$")
    configure_axis_ticks(ax, style)
    fig.text(.5, .965, r"Freeze-out neighborhood: $K^+/\pi^+$", ha="center", fontsize=15)
    fig.text(.5, .925, r"Quark-only BQS; finite-$q$ GBU screening", ha="center", fontsize=13)
    handles.append(Patch(facecolor="#555555", label="Failed support"))
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .127), ncol=2,
               columnspacing=1.5, handlelength=2.4, handletextpad=.6, labelspacing=.5)
    fig.text(.5, .088, r"Probe band: $T_{\rm fo}\pm20$ MeV (not an uncertainty band).",
             ha="center", fontsize=13)
    fig.text(.5, .041, r"Default-curve marker labels: $\sqrt{s_{NN}}$ (GeV). No fit or rescan.",
             ha="center", fontsize=13)

    # Add at most one contour label per level, away from the candidate curves.
    fig.canvas.draw()
    contour_labels = []
    for level, paths in zip(contour.levels, contour.allsegs):
        if not paths:
            continue
        path = max(paths, key=len)
        if len(path) < 4:
            continue
        for fraction in (.28, .65, .45, .80, .15):
            x, y = path[int(fraction * (len(path) - 1))]
            difference = y - freezeout_T(float(x), coefficients)
            if 50 < x < 690 and 58 < y < 180 and min(abs(difference - shift) for shift in (-10, 0, 10)) > 4:
                # Use the exact path's level/coordinate; repeated clabel calls can
                # select another level in Matplotlib 3.10's nearest-path search.
                label = ax.text(x, y, f"{level:.2f}", ha="center", va="center",
                                fontsize=13, color="white", zorder=10)
                label.set_path_effects([pe.Stroke(linewidth=2, foreground="#222222"), pe.Normal()])
                contour_labels.append(label)
                break
    # Suppress only overlapping contour text; preserve all paths and values.
    from scripts.plotting.plot_quality import visible_texts
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    label_ids = {id(label) for label in contour_labels}
    occupied = [text.get_window_extent(renderer) for text in visible_texts(fig) if id(text) not in label_ids]
    for label in contour_labels:
        box = label.get_window_extent(renderer)
        if any(box.overlaps(other) for other in occupied):
            label.set_visible(False)
        else:
            occupied.append(box)
    return fig, {"minimum": minimum, "maximum": maximum, "color_limits": [lower, upper],
                 "contour_levels": [float(level) for level in levels],
                 "visible_contour_labels": [label.get_text() for label in contour_labels if label.get_visible()]}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--freezeout-profile", type=Path, default=ROOT / "config/physics/freezeout/default.toml")
    parser.add_argument("--plotting-support-root", type=Path, default=ROOT,
                        help="Read-only checkout containing the current v2 plotting layer")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--presentation", action="store_true", help="Group-meeting PNG; no paper-facing SOP qualification")
    parser.add_argument("--clean-contours", action="store_true",
                        help="Omit warning/failure overlays and labels; retain missing cells")
    parser.add_argument("--band-MeV", dest="band_MeV", type=float, default=20)
    parser.add_argument("--comparison-color-limits", nargs=2, type=float,
                        help="Fixed color range for a like-for-like presentation; overflow is shown explicitly")
    args = parser.parse_args()
    if args.output_dir.exists():
        raise FileExistsError(f"refusing to overwrite existing case: {args.output_dir}")
    if not math.isfinite(args.band_MeV) or args.band_MeV <= 0:
        raise ValueError("band-MeV must be finite and positive")
    rows, source, coefficients = load_frozen(args.csv, args.source_manifest, args.freezeout_profile)
    if args.presentation:
        return run_presentation(args, rows, source, coefficients)
    if args.band_MeV != 20:
        raise ValueError("non-presentation legacy review mode currently requires band-MeV=20")
    support_root = args.plotting_support_root.resolve()
    sys.path.insert(0, str(support_root))
    from scripts.plotting.plot_style import load_profile, configure_matplotlib
    from scripts.plotting.plot_manifest import (input_record, generator_record, runtime_record,
                                               build_manifest, write_manifest, git_commit)
    from scripts.plotting.plot_quality import export_figure
    from scripts.plotting.validate_plot_artifact import validate_manifest
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    paths = [args.csv, args.source_manifest, args.freezeout_profile,
             support_root / "config/plotting/candidate_aps_v2.toml",
             *(support_root / "scripts/plotting" / name for name in
               ("plot_style.py", "plot_manifest.py", "plot_quality.py", "validate_plot_artifact.py"))]
    inputs = [input_record(path, role="frozen_input_or_readonly_plotting_support", root=ROOT) for path in paths]
    data = neighborhood(rows, coefficients)
    profile = load_profile("candidate_aps_v2")
    font = configure_matplotlib(profile)
    fig, display = plot_map(data, coefficients, profile)
    args.output_dir.mkdir(parents=True)
    outputs, quality = export_figure(fig, args.output_dir / "freezeout_neighborhood_ratio", profile, formats=["png"])
    plt.close(fig)
    manifest = build_manifest(
        asset_id="charged_gbu_freezeout_neighborhood_20261002", figure_family="charged_gbu_freezeout_guidance",
        case_slug=args.output_dir.name, figure_mode="audit", semantic_status="diagnostic_screening_not_a_fit",
        style_profile=profile.profile_id, publication_scope="internal_review", root=ROOT,
        generator=generator_record(Path(__file__), root=ROOT, command=subprocess.list2cmdline([sys.executable, *sys.argv]),
                                   runtime=runtime_record({"matplotlib": matplotlib.__version__, "font": font})),
        inputs=inputs,
        axes=[{"field": "muB_MeV", "source_unit": "MeV", "display_unit": "MeV", "transform": "identity", "label": r"$\mu_B\;(\mathrm{MeV})$"},
              {"field": "T_MeV", "source_unit": "MeV", "display_unit": "MeV", "transform": "identity", "label": r"$T\;(\mathrm{MeV})$"}],
        series=[{"series_id": "R_plus", "state": "nonconverged", "support_rule": "screened native nodes in declared probe band",
                 "mask_rule": "failed nodes masked; no cross-mask contours; outside-band domain unselected"},
                {"series_id": "candidate_curves", "state": "parameter_probe_not_fit", "support_rule": "analytic freeze-out profile; a shifts only",
                 "mask_rule": "curves are coordinate guides, not computed yields"}],
        outputs=outputs, selection_rule="0<=muB<=750 MeV; |T-Tfo(muB)|<=20 MeV; 50<=T<=190 MeV",
        interpolation_policy="heatmap none; display-only contour extraction within four valid corner cells; corner_mask=false",
        connector_policy="forbidden across missing support", missing_value_policy="mask; no zero-fill or signed-value clipping",
        validation={"finite": True, "duplicate_keys": True, "support": True, "strict_gate": False},
        rendering={"column": "double_column", "figure_size_inches": quality["figure_size_inches"],
                   "size_override_reason": "single guidance map with external legend and probe-band/energy keys",
                   "color_route": "undecided_review", "legend_policy": "shared_external_then_case_review",
                   "legend_outside": True, "quality": quality, "delivery_stage": "png_review", "vector_delivery_pending": True,
                   "linear_color_scale": True, "normalization": "no clipping or log; native min/max rounded outward", **display},
        calculation_sha=source["source_git_sha"], postprocess_sha=git_commit(ROOT), source_run_id=source["source_run_id"],
    )
    manifest.update({"manuscript_eligible": False, "current_publication_layer": False, "diagnostic_only": True,
                     "numerical_rescan": False, "production_default_changed": False,
                     "plotting_support_git_sha": git_commit(support_root), "plotting_support_exact_hashes": inputs[3:],
                     "plotting_support_worktree_dirty": bool(subprocess.check_output(
                         ["git", "status", "--porcelain"], cwd=support_root, text=True,
                         encoding="utf-8", errors="replace").strip()),
                     "background_contract": source["background_contract"], "full_grid_rows": len(rows),
                     "neighborhood_selected_nodes": int(data["selected"].sum()),
                     "neighborhood_screened_nodes": int((data["selected"] & ~data["failed"]).sum()),
                     "neighborhood_failed_nodes": int(data["failed"].sum()),
                     "four_corner_valid_cells": int(data["valid_cells"].sum()),
                     "probe_definition": {"band_MeV": 20, "a_shift_MeV": [-10, 0, 10],
                                          "other_coefficients_fixed": coefficients, "not_a_confidence_band": True}})
    # Exact support-file hashes, not the dirty support checkout's HEAD, define its content.
    if any(digest(path) != record["sha256"] for path, record in zip(paths, inputs)):
        raise ValueError("read-only input/support changed during rendering")
    manifest_path = args.output_dir / "plot_manifest.json"
    write_manifest(manifest_path, manifest)
    errors = validate_manifest(manifest_path, repo_root=ROOT)
    print(json.dumps({"png": str(args.output_dir / "freezeout_neighborhood_ratio.png"),
                      "screened": manifest["neighborhood_screened_nodes"], "failed": manifest["neighborhood_failed_nodes"],
                      "range": [display["minimum"], display["maximum"]], "quality_errors": errors}, ensure_ascii=False))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
