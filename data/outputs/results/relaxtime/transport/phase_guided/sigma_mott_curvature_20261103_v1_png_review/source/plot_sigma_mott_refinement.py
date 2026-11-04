#!/usr/bin/env python3
"""Plot hash-verified sigma diagnostics. No numerical physics is evaluated here."""
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
from matplotlib.ticker import MaxNLocator, LogLocator, NullFormatter
from PIL import Image
from scripts.plotting.plot_accessibility import accessibility_record
from scripts.plotting.plot_bundle import build_bundle
from scripts.plotting.plot_manifest import (build_manifest, generator_record, input_record,
    output_record, runtime_record, sha256_file, write_manifest)
from scripts.plotting.plot_quality import export_figure, inspect_export, placement_limits
from scripts.plotting.plot_series import SeriesStyle, SeriesStyleRegistry
from scripts.plotting.plot_style import load_profile, style_context, configure_axis_ticks
from scripts.plotting.validate_plot_artifact import validate_manifest

SHARED = ("plot_style.py", "plot_series.py", "plot_manifest.py", "plot_quality.py",
    "plot_accessibility.py", "plot_bundle.py", "plot_provenance.py", "plot_delivery.py", "validate_plot_artifact.py")
STYLES = SeriesStyleRegistry("sigma_mott_checks_v1", (
    SeriesStyle("sig_m200", "#0072B2", "-", "o"),
    SeriesStyle("sig_m100", "#D55E00", "--", "s"),
    SeriesStyle("sig_m050", "#008565", "-.", "D"),
    SeriesStyle("sig_m020", "#AD5E8D", ":", "^"),
    SeriesStyle("sig_c000", "#000000", "-", ""),
    SeriesStyle("sig_p020", "#666666", "--", ""),
    SeriesStyle("s6", "#D55E00", "-.", "s"),
    SeriesStyle("hybrid_s6_fine_x0.001", "#0072B2", "-", "D"),
    SeriesStyle("combined_s6_x0.001", "#008565", "--", "D"),
    SeriesStyle("inner_s7_controls", "#000000", "None", "x"),
    SeriesStyle("peak_area", "#AD5E8D", "-", "D"),
    SeriesStyle("continuous", "#0072B2", "--", "o"),
    SeriesStyle("discrete", "#D55E00", ":", "s"),
))
LABELS = {
    "sig_m200": r"$\Delta\xi=-2\times10^{-4}$", "sig_m100": r"$\Delta\xi=-10^{-4}$",
    "sig_m050": r"$\Delta\xi=-5\times10^{-5}$", "sig_m020": r"$\Delta\xi=-2\times10^{-5}$",
    "sig_c000": r"$\Delta\xi=0$", "sig_p020": r"$\Delta\xi=+2\times10^{-5}$",
    "s6": r"Original quadrature ($10^{-6}$)",
    "s7": r"Same grid, $s$ tol. $10^{-7}$",
    "hybrid_s6_fine_x0.001": "Continuous near-threshold control",
    "combined_s6_x0.001": "Above + fixed cutoff membership",
    "inner_s7_controls": r"Inner controls: $s$ tol. $10^{-7}$",
    "peak_area": "Unweighted peak area", "continuous": "Continuum thermal weight",
    "discrete": "Original discrete thermal weight",
}


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def verify_folder(folder):
    summary = json.loads((folder/"summary.json").read_text(encoding="utf-8-sig"))
    for name, digest in summary["files"].items():
        if sha256_file(folder/name) != digest:
            raise ValueError("Frozen input changed: " + str(folder/name))
    return summary


def lines(ax, rows, xfield, yfield, style, *, markers=True, markevery=None):
    rows = sorted(rows, key=lambda r: float(r[xfield]))
    kwargs = style.plot_kwargs()
    if not markers:
        kwargs["markevery"] = []
    elif markevery is not None:
        kwargs["markevery"] = markevery
    ax.plot([float(r[xfield]) for r in rows], [float(r[yfield]) for r in rows],
            **kwargs, linewidth=1.2, markersize=5.2, markerfacecolor="white", markeredgewidth=1.0)


def draw(kind, data, args, profile, inputs):
    size = (6.75, 6.7)
    if kind == "spectrum":
        selected = ["sig_m200", "sig_m100", "sig_m050", "sig_m020", "sig_c000", "sig_p020"]
        scope = "Both scattering-channel panels; markers identify sampled mixed-scalar peak positions."
        xlabel = r"$x=s-4M_l^2\;(\mathrm{fm}^{-2})$"
        ylabels = [r"$\Sigma\;(\mathrm{fm}^{2})$"]*2
        title = r"$\sigma$ threshold: full scattering cross sections"
        count = len(data["spectrum"])
        units = ("fm^-2", "fm^2")
    elif kind == "slopes":
        selected = ["s6", "hybrid_s6_fine_x0.001", "combined_s6_x0.001", "inner_s7_controls"]
        scope = "Lines use uniform s6 endpoints at all h. Black crosses show the separate uniform inner s7 triplets for each control. Combined and near-threshold controls coincide in tau_ubar; the combined line is drawn only in conductivity."
        xlabel = r"$h$"
        ylabels = [r"$J_{\tau_{\bar u}}(h)$", r"$J_{\sigma_{\rm el}/T}(h)$"]
        title = r"$\sigma$ threshold: local log-slope difference"
        count = len(data["slopes"])
        units = ("dimensionless", "dimensionless")
    elif kind == "curvature":
        selected = ["s6", "hybrid_s6_fine_x0.001", "inner_s7_controls"]
        scope = "Both panels use tau_ubar. Lines use uniform s6 endpoints; black crosses independently show raw and near-threshold inner s7 triplets. The lower panel is the stored centered second difference of ln(tau_ubar), not a pointwise second derivative."
        xlabel = r"$h$"
        ylabels = [r"$J(h)$", r"$C(h)=J(h)/h$"]
        title = r"$\sigma$ threshold: $\ln\tau_{\bar u}$ curvature diagnostic"
        count = len(data["slopes"])
        units = ("dimensionless", "dimensionless")
    elif kind == "areas":
        selected = ["peak_area", "continuous", "discrete"]
        scope = "Unweighted area applies to upper panel; thermal weights apply to lower panel."
        xlabel = r"$x_r\;(\mathrm{fm}^{-2})$"
        ylabels = [r"$A_{\rm peak}$", r"$R_{\rm peak}\;(\mathrm{fm}^{2})$"]
        title = r"$\sigma$ peak: area and thermal contribution"
        count = len(data["peaks"])+len(data["thermal_peaks"])
        units = ("fm^-2", "dimensionless; fm^2")
    else:
        raise ValueError(kind)
    styles = {s.series_id:s for s in STYLES.select(selected)}
    labels = [LABELS[key] for key in selected]
    if kind=="spectrum":
        counts={k:sum(r["point_id"]==k and r["channel"] in ("uubar_to_uubar","uubar_to_ddbar")
            and 1e-6<=float(r["x_inv_fm2"])<=1e-3 for r in data["spectrum"]) for k in selected}
    elif kind=="slopes":
        counts={k:sum(r["variant"]==k and (r["observable"]=="sigma_over_T" or
                      (r["observable"]=="tau_ubar" and k!="combined_s6_x0.001"))
                      for r in data["slopes"]) for k in selected if k!="inner_s7_controls"}
        counts["inner_s7_controls"] = sum(r["variant"] in ("s7","hybrid_s7_fine_x0.001","combined_s7_x0.001")
            and (r["observable"]=="sigma_over_T" or (r["observable"]=="tau_ubar" and r["variant"]!="combined_s7_x0.001"))
            for r in data["slopes"])
    elif kind=="curvature":
        counts={k:2*sum(r["observable"]=="tau_ubar" and r["variant"] in
                       (("s7","hybrid_s7_fine_x0.001") if k=="inner_s7_controls" else (k,))
                       for r in data["slopes"]) for k in selected}
    else:
        counts={k:sum(r["channel"]=="uubar_to_uubar" for r in
                      (data["peaks"] if k=="peak_area" else data["thermal_peaks"])) for k in selected}
    with style_context(profile) as context:
        figure, axes = plt.subplots(2, 1, figsize=size, sharex=True)
        context.bind(figure)
        figure.subplots_adjust(left=.17, right=.965, bottom=.11, top=.67 if kind=="slopes" else .705, hspace=.18)
        if kind == "spectrum":
            for ax, channel, text in zip(axes,("uubar_to_uubar","uubar_to_ddbar"),
                                        (r"$u\bar u\to u\bar u$",r"$u\bar u\to d\bar d$")):
                for key in selected:
                    rows = [r for r in data["spectrum"] if r["point_id"] == key and r["channel"] == channel
                            and 1e-6 <= float(r["x_inv_fm2"]) <= 1e-3]
                    landmarks = [r for r in data["peaks"] if r["point_id"] == key and r["channel"] == channel]
                    rows=sorted(rows,key=lambda r:float(r["x_inv_fm2"]))
                    landmark_indices=[]
                    if landmarks:
                        landmark_indices=[i for i,r in enumerate(rows)
                                          if any(r["native_sigma_fm2"]==p["full_peak_height_fm2"] for p in landmarks)]
                    lines(ax,rows,"x_inv_fm2","native_sigma_fm2",styles[key],markevery=landmark_indices)
                ax.set_xscale("log");ax.set_yscale("log")
                ax.set_xlim(1e-6,1e-3)
                ax.text(.035,.90,text,transform=ax.transAxes,ha="left",va="top")
        elif kind == "slopes":
            for ax, observable in zip(axes,("tau_ubar","sigma_over_T")):
                for key in selected:
                    if key=="combined_s6_x0.001" and observable=="tau_ubar":
                        continue
                    variants = ("s7","hybrid_s7_fine_x0.001","combined_s7_x0.001") if key=="inner_s7_controls" else (key,)
                    rows = [r for r in data["slopes"] if r["variant"] in variants and r["observable"] == observable
                            and not (observable=="tau_ubar" and r["variant"]=="combined_s7_x0.001")]
                    lines(ax,rows,"h","jump",styles[key])
                ax.set_xscale("log")
                ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
                ax.ticklabel_format(axis="y",style="plain",useOffset=False)
        elif kind == "curvature":
            for ax, field in zip(axes,("jump","jump_over_h")):
                for key in selected:
                    variants=("s7","hybrid_s7_fine_x0.001") if key=="inner_s7_controls" else (key,)
                    rows=[r for r in data["slopes"] if r["observable"]=="tau_ubar" and r["variant"] in variants]
                    lines(ax,rows,"h",field,styles[key])
                ax.set_xscale("log")
                ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
                ax.ticklabel_format(axis="y",style="plain",useOffset=False)
        else:
            raw = [r for r in data["peaks"] if r["channel"] == "uubar_to_uubar"]
            thermal = [r for r in data["thermal_peaks"] if r["channel"] == "uubar_to_uubar"]
            lines(axes[0],raw,"peak_x_inv_fm2","mixed_area_dimensionless",styles["peak_area"])
            lines(axes[1],thermal,"x_peak_inv_fm2","continuous_mixed_square_rate_fm2",styles["continuous"])
            lines(axes[1],thermal,"x_peak_inv_fm2","discrete_exact_mixed_square_rate_fm2",styles["discrete"])
            for ax in axes:
                ax.set_xscale("log");ax.set_yscale("log")
        for ax, label in zip(axes,ylabels):
            ax.set_ylabel(label,labelpad=9)
            configure_axis_ticks(ax,profile)
            if ax.get_yscale()=="log":
                ax.yaxis.set_major_locator(LogLocator(base=10,numticks=12))
                ax.yaxis.set_minor_locator(LogLocator(base=10,subs=range(2,10),numticks=100))
                ax.yaxis.set_minor_formatter(NullFormatter())
            if ax.get_xscale()=="log":
                ax.xaxis.set_minor_locator(LogLocator(base=10,subs=range(2,10),numticks=100))
        axes[-1].set_xlabel(xlabel)
        figure.text(.55,.98,title,ha="center",va="top",fontsize=13)
        figure.text(.55,.925,r"$T=125.736\;\mathrm{MeV},\quad\mu_B=900\;\mathrm{MeV}$",ha="center",va="top",fontsize=13)
        handles = [Line2D([],[],**styles[key].plot_kwargs(),markersize=5.2,markerfacecolor="white") for key in selected]
        figure.legend(handles,labels,loc="upper center",bbox_to_anchor=(.55,.87),ncol=2 if kind=="spectrum" else 1,
                      handlelength=2.1,columnspacing=1.1,labelspacing=.30,borderaxespad=0,frameon=False)
        stem=args.output/("sigma_"+kind)
        outputs,quality=export_figure(figure,stem,profile,formats=["png"],intended_width_inches=size[0])
    color=stem.with_suffix(".png");gray=stem.with_name(stem.name+"_gray").with_suffix(".png")
    with Image.open(color) as im:
        im.convert("L").convert("RGB").save(gray,dpi=(profile.dpi,profile.dpi))
    gray_record=output_record(gray,fmt="png",dpi=profile.dpi,vector=False)
    gray_record["inspection"]=inspect_export(gray);outputs.append(gray_record)
    record=build_manifest(asset_id="sigma_mott_"+kind,figure_family="sigma_mott_refinement",
        case_slug=args.output.name,figure_mode="audit",semantic_status="sigma_mott_numerical_diagnostic",
        style_profile=profile.profile_id,publication_scope="internal_review",
        generator=generator_record(__file__,command=" ".join(sys.argv),
            source_snapshot=args.results/"source"/Path(__file__).name,
            runtime=runtime_record({"matplotlib":matplotlib.__version__})),inputs=inputs,
        axes=[dict(field="h" if kind in ("slopes","curvature") else "x_inv_fm2",source_unit=units[0],display_unit=units[0],scale="log",label=xlabel),
              *[dict(field=kind+"_panel_"+str(i),source_unit=u,display_unit=u,
                      scale="linear" if kind in ("slopes","curvature") else "log",label=label) for i,(label,u) in
                enumerate(zip(ylabels,("dimensionless","fm^2") if kind=="areas" else (units[1],units[1])))]],
        series=[dict(series_id=k,state="numerical_diagnostic",support_rule="frozen completed spectrum or explicitly qualified finite-window analysis",
                     mask_rule="none",point_count=counts[k]) for k in selected],outputs=outputs,
        selection_rule="Declared two-channel near-threshold spectra" if kind=="spectrum" else "Declared observables and completed same-window diagnostics",
        interpolation_policy="none",connector_policy="forbidden",missing_value_policy="required sources fail; no imputation",
        validation=dict(finite=True,duplicate_keys=True,support=True,solver_called=False,
                        canonical_data_modified=False,author_review="pending"),
        calculation_sha="2e5c411c716d5cbc524c3fc6cfb335fe05b39444",source_run_id="trsigma-mott-20261101-v1",
        rendering=dict(delivery_stage="png_review",vector_delivery_pending=True,review_contract="prd_review_v1",
            figure_size_inches=list(size),column="double_column",font=dict(context.font),
            size_override_reason="Paired scientific diagnostics with shared parameter and curve identity keys",
            typography_exception="dense_composite_review_compact_typography",
            typography_exception_reason="Internal PNG review of two panels and nested physics subscripts; minimum 1.5 mm remains enforced",
            series_styles=STYLES.record(selected),legend_policy="declared_geometry_checked",legend_outside=True,
            legend_placements=[dict(legend_index=0,host_axes_index=None,title="",labels=labels,
                placement="outside_axes",location="upper center",scope=scope)],
            legend_layout_reason="Shared key avoids obstructing threshold peaks and small slope differences",
            color_route="undecided_review",quality=quality,placement_limits=placement_limits(quality,profile,outputs),
            accessibility=accessibility_record([dict(value=k,color=styles[k].color,
                noncolor_encoding=f"{styles[k].linestyle}; {styles[k].marker or 'line'}") for k in selected],
                backgrounds=["#FFFFFF"],background_scope="white axes without fills"),
            connections="Straight guides between supplied samples; precision controls are marker-only; no fitted derivative",
            coordinate_transform="All physical coordinates, areas, thermal rates and log-secant differences are precomputed on quote-remote",
            grayscale=dict(method="Pillow L then RGB",color_sha256=sha256_file(color),gray_sha256=sha256_file(gray))))
    record.update(delivery_stage="png_review",manuscript_eligible=False,current_publication_layer=False,
                  vector_delivery_pending=True,caption=input_record(args.results/"caption.md",role="figure_caption"))
    plt.close(figure)
    return record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case",type=Path)
    parser.add_argument("review")
    parser.add_argument("output",type=Path)
    parser.add_argument("results",type=Path)
    parser.add_argument("--figures",default="spectrum,slopes,areas")
    args=parser.parse_args()
    chosen=args.figures.split(",")
    assert set(chosen)<={"spectrum","slopes","areas","curvature"} and len(chosen)==len(set(chosen))
    review_dir=args.case/args.review
    summary=verify_folder(review_dir)
    if set(chosen)-{"spectrum"}:
        assert summary["status"]=="completed_diagnostic_review"
    verify_folder(args.case/"spectrum_v2")
    data=dict(spectrum=read_csv(args.case/"spectrum_v2/cross_section_spectrum.csv"),
              peaks=read_csv(review_dir/"peak_scaling.csv"),slopes=[],thermal_peaks=[])
    files=[args.case/"contract.json",args.case/"spectrum_v2/summary.json",
           args.case/"spectrum_v2/cross_section_spectrum.csv",review_dir/"summary.json",review_dir/"peak_scaling.csv"]
    if {"slopes","curvature"}.intersection(chosen):
        assert "precision_policy" in summary, "Slopes require the uniform-precision review"
        data["slopes"]=read_csv(review_dir/"centered_slopes.csv");files.append(review_dir/"centered_slopes.csv")
        assert {"s6","s7","hybrid_s6_fine_x0.001","hybrid_s7_fine_x0.001","combined_s6_x0.001","combined_s7_x0.001"} <= {r["variant"] for r in data["slopes"]}
    if "areas" in chosen:
        verify_folder(args.case/"thermal_audit_v1")
        data["thermal_peaks"]=read_csv(args.case/"thermal_audit_v1/thermal_peak_areas.csv")
        files.extend((args.case/"thermal_audit_v1/summary.json",args.case/"thermal_audit_v1/thermal_peak_areas.csv"))
    args.output.mkdir(parents=True,exist_ok=False);args.results.mkdir(parents=True,exist_ok=False)
    (args.results/"source").mkdir()
    shutil.copyfile(__file__,args.results/"source"/Path(__file__).name)
    (args.results/"caption.md").write_text(
        "# Sigma Mott neighborhood diagnostics\n\n"
        "Mode A, mu_B=900 MeV, T=125.7360002645551 MeV, xi_M=0.015327819995582102. "
        "The original frozen kernel is unchanged. Delta xi is xi-xi_M; x=s-4M_l^2. The spectrum panels "
        "show the full blocked cross section in the two light quark-antiquark channels; markers identify "
        "sampled mixed-scalar resonance peaks. The mixed scalar contains both sigma and sigma-prime. "
        "Peak areas integrate the positive mixed-scalar square over x_r plus/minus 32 half-widths. "
        "A_peak=integral Sigma_mixed dx is dimensionless; R_peak includes the recorded continuum or "
        "original discrete incoming thermal weight and has units fm^2. The slope diagnostic is "
        "J_y(h)=ln[y(xi_M+h)/y(xi_M)]/h-ln[y(xi_M)/y(xi_M-h)]/h. The continuum near-threshold "
        "control replaces only x<0.001 fm^-2 in two light annihilation channels; it is a diagnostic "
        "and does not promote a new production quadrature. Every line uses uniform s-tolerance 10^-6 "
        "at all three secant endpoints. Black crosses separately show the same inner triplet at 10^-7 "
        "for the original, continuous near-threshold and combined controls; they coincide with the "
        "line markers at this scale. Their numerical differences remain in precision_comparison.csv. "
        "The combined counterfactual also holds the uu and ud hard-upper-cutoff event membership "
        "fixed to the Mott center. It attributes discrete boundary effects, not a new continuum solution. "
        "These two qq processes do not contribute to tau_ubar, so combined and near-threshold values "
        "are identical there and only one line is drawn. The mixed-scalar peak area includes both "
        "mixed poles within the squared flavor propagator, but not its cross terms with other amplitudes; "
        "the full scattering spectra retain those interferences. Lines are sample guides, not interpolation "
        "for numerical analysis. Finite-window slopes are not proof of differentiability. Full insertion "
        "width 6.75 inches, PNG review only. Numerical conclusions and figure author acceptance are separate.\n",
        encoding="utf-8")
    if "curvature" in chosen:
        with (args.results/"caption.md").open("a",encoding="utf-8") as stream:
            stream.write(
                "\nThe curvature figure uses y=tau_ubar and L=ln(y/y_ref), with any fixed positive y_ref. "
                "Its upper panel shows the stored J(h) for left/right log-secant differences. Its lower "
                "panel reads the stored jump_over_h field, C(h)=[L(xi_M+h)-2 L(xi_M)+L(xi_M-h)]/h^2. "
                "C(h) is a finite-window estimator of L'', not of tau_ubar''; when derivatives exist, "
                "L''=tau_ubar''/tau_ubar-(tau_ubar'/tau_ubar)^2. Positive C^1 tau_ubar shares any "
                "uncancelled divergent curvature with L. These are four window sizes plus one inner "
                "s7 precision comparison, without a fit or extrapolation to h=0. The continuum control "
                "replaces only the two declared near-threshold annihilation bands; it is not a full "
                "outer-quadrature convergence certificate. The present curve does not establish a "
                "stable h^(-1/2) divergence or classify full transport differentiability.\n")
    inputs=[input_record(p,role="frozen_sigma_diagnostic") for p in files]
    for path in [ROOT/"config/plotting/candidate_aps_v2.toml",*[ROOT/"scripts/plotting"/name for name in SHARED]]:
        target=args.results/"source"/path.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(path,target)
        entry=input_record(path,role="plot_style_contract" if path.suffix==".toml" else "shared_plotting_code")
        entry["source_snapshot"]=input_record(target,role="plotting_source_snapshot");inputs.append(entry)
    profile=load_profile("candidate_aps_v2")
    records=[draw(kind,data,args,profile,inputs) for kind in chosen]
    bundle=build_bundle(dict(case_slug=args.output.name,delivery_stage="png_review",manuscript_eligible=False,
        current_publication_layer=False,charts=[dict(asset_id=r["asset_id"]) for r in records]),records)
    manifest=args.output/"plot_manifest.json";write_manifest(manifest,bundle)
    errors=validate_manifest(manifest,repo_root=ROOT)
    (args.results/"validation.json").write_text(json.dumps(dict(passed=not errors,errors=errors),indent=2)+"\n",encoding="utf-8")
    if errors:
        raise ValueError("; ".join(errors))
    print(json.dumps(dict(figures=len(records),passed=True,manifest=str(manifest))))


if __name__=="__main__":
    main()
