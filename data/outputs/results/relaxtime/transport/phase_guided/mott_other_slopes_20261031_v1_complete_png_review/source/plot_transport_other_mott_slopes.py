#!/usr/bin/env python3
"""Plot frozen other-Mott values and remotely tabulated slopes, without solving."""
from pathlib import Path
import argparse
import csv
import json
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

STYLES = SeriesStyleRegistry("other_mott_slopes_v1", (
    SeriesStyle("fine", "#0072B2", "-", "o"),
    SeriesStyle("original", "#D55E00", "None", "s"),
    SeriesStyle("tighter", "#000000", "None", "x"),
    SeriesStyle("xi_halves", "#008565", "None", "+"),
    SeriesStyle("mott", "#666666", ":", ""),
))
LABELS = dict(fine=r"Refined $\xi$ grid", original=r"Original $\xi$ nodes",
              tighter=r"Tighter $s$ sampling", xi_halves=r"Halved $\xi$ intervals")
SYMBOLS = dict(pi=r"\pi", K="K", sigma=r"\sigma", eta=r"\eta")
QUANTITIES = {
    **{"tau_"+s: (rf"\tau_{{{symbol}}}", "fm") for s, symbol in
       (("u", "u"), ("d", "d"), ("s", "s"), ("ubar", r"\bar{u}"),
        ("dbar", r"\bar{d}"), ("sbar", r"\bar{s}"))},
    "eta_over_s": (r"\eta/s", "dimensionless"),
    "zeta_over_s": (r"\zeta/s", "dimensionless"),
    "sigma_over_T": (r"\sigma_{\rm el}/T", "dimensionless"),
}
SHARED_CODE = ("plot_style.py", "plot_series.py", "plot_manifest.py", "plot_quality.py",
    "plot_accessibility.py", "plot_bundle.py", "plot_provenance.py", "plot_delivery.py", "validate_plot_artifact.py")


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def draw(curve, observable, data, args, profile, manifest_inputs):
    raw = sorted((r for r in data["observations"] if r["curve_id"] == curve), key=lambda r: float(r["xi"]))
    original = [r for r in raw if r["original_point_id"]]
    references = sorted((r for r in data["references"] if r["curve_id"] == curve), key=lambda r: float(r["xi"]))
    slopes = [r for r in data["slopes"] if r["curve_id"] == curve and r["observable"] == observable]
    native = {(r["point_left"], r["point_right"]): r for r in slopes}
    matched = [(native[r["point_left"], r["point_right"]]["delta_xi_midpoint"], r["tighter_log_slope"])
               for r in data["slope_convergence"] if r["curve_id"] == curve and r["observable"] == observable
               and (r["point_left"], r["point_right"]) in native]
    extra = {r["point_id"] for r in references if r["point_kind"] == "xi_midpoint_refinement"}
    halves = [r for r in data["reference_slopes"] if r["curve_id"] == curve and r["observable"] == observable
              and (r["point_left"] in extra or r["point_right"] in extra)]
    selected = ["fine"] + (["original"] if original else []) + (["tighter"] if references else [])
    selected += (["xi_halves"] if halves else []) + ["mott"]
    styles = {s.series_id: s for s in STYLES.select(selected)}
    labels = [LABELS.get(k, "$"+SYMBOLS[raw[0]["meson"]]+"$ Mott") for k in selected]
    symbol, unit = QUANTITIES[observable]
    ylabel = "$"+symbol+(r"\;(\mathrm{fm})" if unit == "fm" else "")+"$"
    slope_label = r"$\Delta\ln("+symbol+r")/\Delta\xi$"
    size = (6.75, 6.)
    with style_context(profile) as context:
        figure, axes = plt.subplots(2, 1, figsize=size, sharex=True)
        context.bind(figure)
        figure.subplots_adjust(left=.16, right=.965, bottom=.11, top=.74 if len(selected)>4 else .79, hspace=.13)
        for key, points in (("fine", raw), ("original", original), ("tighter", references)):
            if points:
                axes[0].plot([float(r["delta_xi"]) for r in points], [float(r[observable]) for r in points],
                    **styles[key].plot_kwargs(), markersize=4.8, markerfacecolor="white", markeredgewidth=.85,
                    linewidth=1.2, zorder=3 if key == "fine" else 4)
        for key, values in (("fine", [(r["delta_xi_midpoint"], r["log_slope"]) for r in slopes]),
                            ("tighter", matched), ("xi_halves", [(r["delta_xi_midpoint"], r["log_slope"]) for r in halves])):
            if values:
                values = sorted(values, key=lambda pair: float(pair[0]))
                axes[1].plot([float(x) for x, _ in values], [float(y) for _, y in values],
                    **styles[key].plot_kwargs(), markersize=4.8, markerfacecolor="white", markeredgewidth=.85,
                    linewidth=1.2, zorder=3 if key == "fine" else 4)
        for ax in axes:
            ax.axvline(0, **styles["mott"].plot_kwargs(), linewidth=.8, zorder=1)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
            ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
            ax.ticklabel_format(axis="both", style="plain", useOffset=False)
            configure_axis_ticks(ax, profile)
        axes[0].set_ylabel(ylabel, labelpad=9)
        axes[1].set_ylabel(slope_label, labelpad=9)
        xlabel = r"$\xi-\xi_{\mathrm{Mott}}$"
        axes[1].set_xlabel(xlabel)
        mode = "A" if raw[0]["window"].startswith("mode_a:") else "B"
        title = f"Mode {mode}: " + rf"$T={float(raw[0]['T_MeV']):.4f}\;\mathrm{{MeV}},\quad\mu_B={float(raw[0]['muB_MeV']):g}\;\mathrm{{MeV}}$"
        figure.text(.54, .965, title, ha="center", va="top", fontsize=13)
        handles = [Line2D([], [], **styles[k].plot_kwargs(), markersize=4.8, markerfacecolor="white") for k in selected]
        figure.legend(handles, labels, loc="upper center", bbox_to_anchor=(.54,.91), ncol=2,
            handlelength=2.2, columnspacing=1.4, labelspacing=.35, borderaxespad=0, frameon=False)
        stem = args.output/f"{curve}_{observable}"
        outputs, quality = export_figure(figure, stem, profile, formats=["png"], intended_width_inches=size[0])
    color, gray = stem.with_suffix(".png"), stem.with_name(stem.name+"_gray").with_suffix(".png")
    with Image.open(color) as im:
        im.convert("L").convert("RGB").save(gray, dpi=(profile.dpi, profile.dpi))
    gray_record = output_record(gray, fmt="png", dpi=profile.dpi, vector=False)
    gray_record["inspection"] = inspect_export(gray)
    outputs.append(gray_record)
    record = build_manifest(asset_id=f"other_mott_{curve}_{observable}", figure_family="transport_other_mott_slopes",
        case_slug=args.output.name, figure_mode="audit", semantic_status="raw_local_slope_diagnostic",
        style_profile=profile.profile_id, publication_scope="internal_review",
        generator=generator_record(__file__, command=" ".join(sys.argv),
            source_snapshot=args.results/"source"/Path(__file__).name,
            runtime=runtime_record({"matplotlib":matplotlib.__version__})), inputs=manifest_inputs,
        axes=[dict(field="delta_xi", source_unit="dimensionless", display_unit="dimensionless", scale="linear", label=xlabel),
              dict(field=observable, source_unit=unit, display_unit=unit, scale="linear", label=ylabel),
              dict(field="log_slope", source_unit="dimensionless", display_unit="dimensionless", scale="linear", label=slope_label)],
        series=[dict(series_id=k, state="threshold_locator" if k=="mott" else "numerical_diagnostic",
            support_rule="located closed threshold" if k=="mott" else "frozen observations or remotely tabulated secants",
            mask_rule="none", point_count=len(dict(fine=raw, original=original, tighter=references, xi_halves=halves).get(k,[raw[0]]))) for k in selected],
        outputs=outputs, selection_rule=f"all supplied completed points for {curve}; sparse precision slopes shown only when endpoints match native intervals",
        interpolation_policy="none", connector_policy="forbidden", missing_value_policy="required sources fail; no imputation",
        validation=dict(finite=True,duplicate_keys=True,support=True,solver_called=False,canonical_data_modified=False,author_review="pending"),
        calculation_sha="2e5c411c716d5cbc524c3fc6cfb335fe05b39444", source_run_id="mott_other_slopes_20261029_v1",
        rendering=dict(delivery_stage="png_review", vector_delivery_pending=True, review_contract="prd_review_v1",
            figure_size_inches=list(size), column="double_column", font=dict(context.font),
            size_override_reason="paired raw-value and slope panels with a readable shared parameter key",
            series_styles=STYLES.record(selected), legend_policy="declared_geometry_checked", legend_outside=True,
            legend_placements=[dict(legend_index=0,host_axes_index=None,title="",labels=labels,placement="outside_axes",
                location="upper center",scope="both panels; original nodes only in upper, halved intervals only in lower")],
            legend_layout_reason="shared external key preserves visibility of small raw transport changes",
            color_route="undecided_review", quality=quality, placement_limits=placement_limits(quality,profile,outputs),
            accessibility=accessibility_record([dict(value=k,color=styles[k].color,
                noncolor_encoding=f"{styles[k].linestyle}; {styles[k].marker or 'vertical threshold'}") for k in selected],
                backgrounds=["#FFFFFF"],background_scope="white axes without fills"),
            connections="straight guides between observed values and tabulated secants; controls are marker-only; no derivative fit",
            coordinate_transform="delta_xi and delta_xi_midpoint were computed on quote-remote relative to the located threshold",
            panel_point_counts=dict(upper=dict(fine=len(raw),original=len(original),tighter=len(references)),
                lower=dict(fine=len(slopes),tighter=len(matched),xi_halves=len(halves))),
            grayscale=dict(method="Pillow L then RGB",color_sha256=sha256_file(color),gray_sha256=sha256_file(gray))))
    record.update(delivery_stage="png_review", manuscript_eligible=False,current_publication_layer=False,
                  vector_delivery_pending=True,caption=input_record(args.results/"caption.md",role="figure_caption"))
    plt.close(figure)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path, help="Verified returned package root")
    parser.add_argument("review", help="Completed review directory in the returned package")
    parser.add_argument("output", type=Path)
    parser.add_argument("results", type=Path)
    parser.add_argument("--runs", required=True)
    parser.add_argument("--comparison-runs", default="")
    parser.add_argument("--observables", default="tau_ubar,tau_sbar,eta_over_s,zeta_over_s,sigma_over_T")
    args = parser.parse_args()
    review_dir = args.case/args.review
    summary = json.loads((review_dir/"summary.json").read_text(encoding="utf-8-sig"))
    assert summary["status"] == "completed_diagnostic_review"
    package_path = args.case/"package_manifest.json"
    entries = {r["path"]:r for r in json.loads(package_path.read_text(encoding="utf-8-sig"))["files"]}
    files = [review_dir/"summary.json", package_path]
    data = dict(observations=[],references=[])
    for key, names in (("observations",args.runs), ("references",args.comparison_runs)):
        for name in filter(None,names.split(",")):
            path = args.case/name/"observables.csv"
            run_summary = args.case/name/"summary.json"
            assert sha256_file(run_summary) == summary["provenance"][name]
            files.extend((path,run_summary))
            data[key].extend(read_csv(path))
        assert len(data[key]) == len({r["point_id"] for r in data[key]})
    for name in ("slopes","slope_convergence","reference_slopes"):
        path = review_dir/(name+".csv")
        data[name] = read_csv(path) if path.exists() else []
        if path.exists():
            files.append(path)
    assert data["slopes"] and len(data["observations"]) == summary["points"]
    for path in files:
        if path == package_path:
            continue
        entry = entries[path.relative_to(args.case).as_posix()]
        assert path.stat().st_size == entry["bytes"] and sha256_file(path) == entry["sha256"],path
    observables = args.observables.split(",")
    assert all(k in QUANTITIES for k in observables) and len(set(observables))==len(observables)
    args.output.mkdir(parents=True,exist_ok=False)
    args.results.mkdir(parents=True,exist_ok=False)
    (args.results/"source").mkdir()
    shutil.copyfile(__file__,args.results/"source"/Path(__file__).name)
    (args.results/"caption.md").write_text(
        "# Other-Mott local slope diagnostics\n\n"
        "Each figure uses the temperature, baryon chemical potential, meson and located threshold of one selected window. "
        "Upper panels show directly computed transport values; lower panels show ln(y_R/y_L)/(xi_R-xi_L) at interval midpoints. "
        "Coordinates and slopes were tabulated on quote-remote. Blue circles use the refined xi grid with s-table rtol=1e-5; "
        "orange squares, if present, identify recomputed original xi nodes only in the upper panel. Black crosses use "
        "rtol=1e-6; in the lower panel they compare identical interval endpoints. Green plus signs distinguish split "
        "intervals adjoining additional xi nodes at the tighter tolerance. Slopes over unmatched, wider gaps between "
        "sparse references are omitted. Lines between observations and secants are visual guides, not fitted derivatives. "
        "The vertical gray line is the diagnostic closed-threshold inverse zero; production B0 and mixed-mode governance "
        "are unchanged. s-table convergence alone is not an error bound for all integration layers. No smoothing, formal "
        "data replacement or publication promotion is performed. Only completed points of the selected diagnostic stages "
        "are shown; a pilot figure does not establish completion of the full meson survey. All figures require the full "
        "6.75-inch insertion width. PNG review only; vector delivery and author review remain pending.\n",encoding="utf-8")
    inputs = [input_record(p,role="frozen_transport_diagnostic",units={"xi":"dimensionless","tau":"fm","log_slope":"dimensionless"}) for p in files]
    for path in [ROOT/"config/plotting/candidate_aps_v2.toml",*[ROOT/"scripts/plotting"/n for n in SHARED_CODE]]:
        target=args.results/"source"/path.relative_to(ROOT)
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(path,target)
        entry=input_record(path,role="plot_style_contract" if path.suffix==".toml" else "shared_plotting_code")
        entry["source_snapshot"]=input_record(target,role="plotting_source_snapshot")
        inputs.append(entry)
    profile=load_profile("candidate_aps_v2")
    records=[draw(curve,observable,data,args,profile,inputs)
             for curve in sorted({r["curve_id"] for r in data["observations"]}) for observable in observables]
    bundle=build_bundle(dict(case_slug=args.output.name,delivery_stage="png_review",manuscript_eligible=False,
        current_publication_layer=False,charts=[dict(asset_id=r["asset_id"]) for r in records]),records)
    manifest=args.output/"plot_manifest.json"
    write_manifest(manifest,bundle)
    errors=validate_manifest(manifest,repo_root=ROOT)
    (args.results/"validation.json").write_text(json.dumps(dict(passed=not errors,errors=errors),indent=2)+"\n",encoding="utf-8")
    if errors:
        raise ValueError("; ".join(errors))
    print(json.dumps(dict(figures=len(records),passed=True,manifest=str(manifest))))


if __name__ == "__main__":
    main()
