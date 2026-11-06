#!/usr/bin/env python3
"""Render frozen eta-refinement observations and secant slopes; no numerical solve."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from PIL import Image

from scripts.plotting.plot_accessibility import accessibility_record
from scripts.plotting.plot_bundle import build_bundle
from scripts.plotting.plot_manifest import (build_manifest, generator_record, input_record,
    output_record, runtime_record, sha256_file, write_manifest)
from scripts.plotting.plot_quality import export_figure, inspect_export, placement_limits
from scripts.plotting.plot_series import SeriesStyle, SeriesStyleRegistry
from scripts.plotting.plot_style import load_profile, style_context, configure_axis_ticks
from scripts.plotting.validate_plot_artifact import validate_manifest

STYLES = SeriesStyleRegistry("eta_xi_grid_comparison_v1", (
    SeriesStyle("fine", "#0072B2", "-", "o"),
    SeriesStyle("coarse", "#D55E00", "--", "s"),
    SeriesStyle("tighter", "#000000", "None", "x"),
    SeriesStyle("xi_halves", "#000000", "None", "+"),
    SeriesStyle("eta_mott", "#666666", ":", ""),
))
LABELS = {"fine": r"Refined $\xi$ grid", "coarse": r"$\Delta\xi=0.01$", "tighter": r"Tighter $s$ sampling",
          "xi_halves": r"Halved $\xi$ intervals", "eta_mott": r"$\eta$ Mott"}
OBSERVABLE_LABELS = {
    "tau_ubar": (r"$\tau_{\bar{u}}\;(\mathrm{fm})$", r"$\Delta\ln\tau_{\bar{u}}/\Delta\xi$", "fm"),
    "sigma_over_T": (r"$\sigma_{\rm el}/T$", r"$\Delta\ln(\sigma_{\rm el}/T)/\Delta\xi$", "dimensionless"),
}


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def figure_record(curve, observable, inputs, review_dir, output_dir, result_dir, profile):
    curve_id = curve["curve_id"]
    raw = sorted((r for r in inputs["observables"] if r["curve_id"] == curve_id), key=lambda r: float(r["xi"]))
    coarse = [r for r in raw if r["original_point_id"]]
    tight = sorted((r for r in inputs.get("tighter_observables", []) if r["curve_id"] == curve_id), key=lambda r: float(r["xi"]))
    fs = [r for r in inputs["fine_slopes"] if r["curve_id"] == curve_id]
    cs = [r for r in inputs["coarse_slopes"] if r["curve_id"] == curve_id]
    ts = [r for r in inputs.get("tighter_slopes", []) if r["curve_id"] == curve_id]
    native_intervals = {(r["point_left"], r["point_right"]): r for r in fs}
    matched = []
    for row in inputs.get("slope_convergence", []):
        pair = (row["point_left"], row["point_right"])
        if row["curve_id"] == curve_id and row["observable"] == observable and pair in native_intervals:
            matched.append(dict(delta_xi_midpoint=native_intervals[pair]["delta_xi_midpoint"],
                                **{observable+"_log_slope": row["tighter_log_slope"]}))
    matched.sort(key=lambda r: float(r["delta_xi_midpoint"]))
    if tight and not matched:
        raise ValueError("Tighter reference points require matching native intervals in slope_convergence.csv")
    new_ids = {r["point_id"] for r in tight
               if r["point_kind"] in ("peak_interval_midpoint", "threshold_interval_midpoint")}
    halves = [r for r in ts if r["point_left"] in new_ids or r["point_right"] in new_ids]
    selected = ["fine", "coarse"] + (["tighter"] if tight else []) + (["xi_halves"] if halves else []) + ["eta_mott"]
    styles = {s.series_id: s for s in STYLES.select(selected)}
    ylabel, slope_label, unit = OBSERVABLE_LABELS[observable]
    size = (6.75, 6.0)
    with style_context(profile) as context:
        figure, axes = plt.subplots(2, 1, figsize=size, sharex=True)
        context.bind(figure)
        figure.subplots_adjust(left=.16, right=.965, bottom=.11, top=.74 if halves else .79, hspace=.13)
        for identity, points, slopes in (("fine", raw, fs), ("coarse", coarse, cs), ("tighter", tight, matched)):
            if not points:
                continue
            style = styles[identity]
            axes[0].plot([float(r["delta_xi"]) for r in points], [float(r[observable]) for r in points],
                **style.plot_kwargs(), markersize=4.8, markerfacecolor="white", markeredgewidth=.85,
                linewidth=1.2, zorder=3 if identity == "fine" else 4)
            # These are the secants tabulated remotely, displayed at interval
            # midpoints. The segments are visual guides, not a derivative fit.
            axes[1].plot([float(r["delta_xi_midpoint"]) for r in slopes], [float(r[observable+"_log_slope"]) for r in slopes],
                **style.plot_kwargs(), markersize=4.8, markerfacecolor="white", markeredgewidth=.85,
                linewidth=1.2, zorder=3 if identity == "fine" else 4)
        if halves:
            axes[1].plot([float(r["delta_xi_midpoint"]) for r in halves],
                [float(r[observable+"_log_slope"]) for r in halves],
                **styles["xi_halves"].plot_kwargs(), markersize=5.2, markeredgewidth=.95, linewidth=1.2, zorder=5)
        for ax in axes:
            ax.axvline(0.0, **styles["eta_mott"].plot_kwargs(), linewidth=.8, zorder=1)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
            ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
            ax.ticklabel_format(axis="both", style="plain", useOffset=False)
            configure_axis_ticks(ax, profile)
        axes[0].set_ylabel(ylabel, labelpad=9)
        axes[1].set_ylabel(slope_label, labelpad=9)
        axes[1].set_xlabel(r"$\xi-\xi_\eta$")
        mode = "A" if curve_id.startswith("a_") else "B"
        temperature = float(raw[0]["T_MeV"])
        mu = float(raw[0]["muB_MeV"])
        title = f"Mode {mode}: " + rf"$T={temperature:.4f}\;\mathrm{{MeV}},\quad\mu_B={mu:g}\;\mathrm{{MeV}}$"
        figure.text(.54, .965, title, ha="center", va="top", fontsize=13)
        handles = [Line2D([], [], **styles[k].plot_kwargs(), markersize=4.8, markerfacecolor="white") for k in selected]
        labels = [LABELS[k] for k in selected]
        figure.legend(handles, labels, loc="upper center", bbox_to_anchor=(.54, .91), ncol=2,
            handlelength=2.2, columnspacing=1.4, labelspacing=.35, borderaxespad=0, frameon=False)
        stem = output_dir/f"{curve_id}_{observable}"
        outputs, quality = export_figure(figure, stem, profile, formats=["png"], intended_width_inches=size[0])
    color_path = stem.with_suffix(".png")
    gray_path = stem.with_name(stem.name+"_gray").with_suffix(".png")
    with Image.open(color_path) as im:
        im.convert("L").convert("RGB").save(gray_path, dpi=(profile.dpi, profile.dpi))
    gray = output_record(gray_path, fmt="png", dpi=profile.dpi, vector=False)
    gray["inspection"] = inspect_export(gray_path)
    outputs.append(gray)
    snapshot = result_dir/"source"/Path(__file__).name
    series = [dict(series_id=k, state="numerical_diagnostic" if k != "eta_mott" else "threshold_locator",
        support_rule="frozen numerical observations and remote secants" if k != "eta_mott" else "endpoint-subtracted eta inverse zero",
        mask_rule="none", point_count=len({"fine":raw, "coarse":coarse, "tighter":tight, "xi_halves":halves}.get(k, [curve]))) for k in selected]
    code_names = ("plot_style.py", "plot_series.py", "plot_manifest.py", "plot_quality.py", "plot_accessibility.py",
                  "plot_bundle.py", "plot_provenance.py", "plot_delivery.py", "validate_plot_artifact.py")
    manifest_inputs = [input_record(review_dir/(name+".csv"), role="raw_transport" if name in ("observables", "tighter_observables") else "remote_derived_slopes",
        units={"xi":"dimensionless", "tau":"fm", "log_slope":"dimensionless"}) for name in inputs]
    manifest_inputs.append(input_record(review_dir/"review.json", role="numerical_review_and_provenance"))
    manifest_inputs.extend([input_record(ROOT/"config/plotting/candidate_aps_v2.toml",role="plot_style_contract"),
                           *[input_record(ROOT/"scripts/plotting"/name,role="shared_plotting_code") for name in code_names]])
    for record in manifest_inputs:
        if record["role"] in ("plot_style_contract", "shared_plotting_code"):
            source = ROOT/record["path"]
            target = result_dir/"source"/record["path"]
            target.parent.mkdir(parents=True,exist_ok=True)
            if not target.exists():
                shutil.copy2(source,target)
            record["source_snapshot"] = input_record(target,role="plotting_source_snapshot")
    manifest = build_manifest(asset_id=f"eta_refinement_{curve_id}_{observable}",figure_family="transport_eta_xi_refinement",
        case_slug=output_dir.name,figure_mode="audit",semantic_status="raw_xi_refinement_comparison",
        style_profile=profile.profile_id,publication_scope="internal_review",
        generator=generator_record(__file__,command=" ".join(sys.argv),source_snapshot=snapshot,
            runtime=runtime_record({"matplotlib":matplotlib.__version__})),inputs=manifest_inputs,
        axes=[dict(field="delta_xi",source_unit="dimensionless",display_unit="dimensionless",scale="linear",label=r"$\xi-\xi_\eta$"),
              dict(field=observable,source_unit=unit,display_unit=unit,scale="linear",label=ylabel),
              dict(field=observable+"_log_slope",source_unit="dimensionless",display_unit="dimensionless",scale="linear",label=slope_label)],
        series=series,outputs=outputs,selection_rule=f"all supplied transport points for {curve_id}; original controls by source ID; lower-panel tighter crosses use native matching intervals, plus signs use intervals adjoining added midpoints",
        interpolation_policy="none",connector_policy="forbidden",missing_value_policy="missing required points fail; no replacement",
        validation=dict(finite=True,duplicate_keys=True,support=True,solver_called=False,canonical_data_modified=False,author_review="pending"),
        calculation_sha="2e5c411c716d5cbc524c3fc6cfb335fe05b39444",source_run_id="mott_eta_fine_20261029_v1",
        rendering=dict(delivery_stage="png_review",vector_delivery_pending=True,review_contract="prd_review_v1",
            figure_size_inches=list(size),column="double_column",font=dict(context.font),
            size_override_reason="two vertically aligned panels with readable parameter and comparison legends",
            series_styles=STYLES.record(selected),legend_policy="declared_geometry_checked",legend_outside=True,
            legend_placements=[dict(legend_index=0,host_axes_index=None,title="",labels=labels,
                placement="outside_axes",location="upper center",scope="both panels; halved xi intervals appear only in the lower panel")],
            legend_layout_reason="external shared key leaves the small transport modulation unobscured",
            color_route="undecided_review",quality=quality,placement_limits=placement_limits(quality,profile,outputs),
            accessibility=accessibility_record([dict(value=k,color=styles[k].color,
                noncolor_encoding=f"{styles[k].linestyle}; {styles[k].marker or 'vertical threshold'}") for k in selected],
                backgrounds=["#FFFFFF"],background_scope="white axes without filled layers"),
            connections="straight segments join fine/coarse observations and remote secant midpoints; tighter references and xi halves are marker-only; no fit",
            panel_point_counts=dict(upper=dict(fine=len(raw),coarse=len(coarse),tighter=len(tight)),
                lower=dict(fine=len(fs),coarse=len(cs),tighter=len(matched),xi_halves=len(halves))),
            interval_comparison="tighter crosses use exactly the same endpoints as a native fine-grid secant; plus signs distinguish newly halved intervals; larger reference-only gaps are omitted from the slope panel",
            coordinate_transform="xi-xi_eta is tabulated remotely; xi_eta is the endpoint-subtracted eta threshold zero",
            grayscale=dict(method="Pillow L then RGB",color_sha256=sha256_file(color_path),gray_sha256=sha256_file(gray_path))))
    manifest.update(delivery_stage="png_review",manuscript_eligible=False,current_publication_layer=False,
        vector_delivery_pending=True,caption=input_record(result_dir/"caption.md",role="figure_caption"))
    plt.close(figure)
    return manifest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("review",type=Path)
    parser.add_argument("output",type=Path)
    parser.add_argument("results",type=Path)
    args=parser.parse_args()
    review=json.loads((args.review/"review.json").read_text(encoding="utf-8-sig"))
    assert review["status"] == "completed_diagnostic_review"
    for filename, sha in review["output_sha256"].items():
        assert sha256_file(args.review/filename)==sha,filename
    args.output.mkdir(parents=True,exist_ok=False)
    args.results.mkdir(parents=True,exist_ok=False)
    (args.results/"source").mkdir()
    shutil.copyfile(__file__,args.results/"source"/Path(__file__).name)
    (args.results/"caption.md").write_text(
        "# Eta reverse-Mott xi refinement\n\n"
        "The horizontal coordinate is xi-xi_eta, tabulated remotely relative to each curve's eta threshold. "
        "Upper panels show the directly computed transport values. Lower panels show "
        "ln(y_right/y_left)/(xi_right-xi_left) at each interval midpoint; connecting segments "
        "are guides between these tabulated secants, not a fitted derivative. Orange squares "
        "are the original xi nodes recomputed with the same s-sampling tolerance as the blue refined grid. "
        "The s-table relative tolerance is 1e-5 for both grids and 1e-6 for the black reference markers; "
        "the outer momentum and angular quadrature, inner t quadrature, and production B0 kernel are unchanged. "
        "Upper-panel black crosses, where present, include the located threshold and a peak-interval midpoint. "
        "In the lower panel, black crosses compare exactly the same interval endpoints as a native fine-grid secant; "
        "black plus signs show the newly halved intervals at the tighter tolerance. Slopes over larger gaps between "
        "the sparse reference points are not displayed. All shown slopes and their midpoint coordinates were tabulated remotely. The gray vertical line "
        "locates the zero of the endpoint-subtracted eta light-mode inverse. This diagnostic "
        "locator does not replace the production scattering kernel. No smoothing or replacement "
        "of the formal data is performed. The two mu_B=0 curves have slightly different temperatures "
        "and are not independent physical repetitions. All panels require the full 6.75-inch insertion width. "
        "PNG review only; manuscript eligibility and vector delivery remain pending.\n",encoding="utf-8")
    names=("observables","fine_slopes","coarse_slopes","tighter_observables","tighter_slopes","slope_convergence")
    inputs={name:rows(args.review/(name+".csv")) for name in names if (args.review/(name+".csv")).exists()}
    profile=load_profile("candidate_aps_v2")
    records=[figure_record(curve,observable,inputs,args.review,args.output,args.results,profile)
             for curve in review["curves"] for observable in OBSERVABLE_LABELS]
    bundle=build_bundle(dict(case_slug=args.output.name,delivery_stage="png_review",manuscript_eligible=False,
        current_publication_layer=False,charts=[dict(asset_id=r["asset_id"]) for r in records]),records)
    write_manifest(args.output/"plot_manifest.json",bundle)
    errors=validate_manifest(args.output/"plot_manifest.json",repo_root=ROOT)
    (args.results/"validation.json").write_text(json.dumps(dict(passed=not errors,errors=errors),indent=2)+"\n",encoding="utf-8")
    if errors:
        raise ValueError("; ".join(errors))
    print(json.dumps(dict(figures=len(records),passed=True,manifest=str(args.output/"plot_manifest.json"))))


if __name__ == "__main__":
    main()
