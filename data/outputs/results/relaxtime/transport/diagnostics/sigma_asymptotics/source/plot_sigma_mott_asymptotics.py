#!/usr/bin/env python3
"""Plot completed remote side-limit diagnostics without recalculating numbers."""
from pathlib import Path
import argparse
import csv
import json
import shutil
import sys

ROOT=Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator,MaxNLocator,NullFormatter
from PIL import Image
from scripts.plotting.plot_accessibility import accessibility_record
from scripts.plotting.plot_bundle import build_bundle
from scripts.plotting.plot_manifest import build_manifest,generator_record,input_record,output_record,runtime_record,sha256_file,write_manifest
from scripts.plotting.plot_quality import export_figure,inspect_export,placement_limits
from scripts.plotting.plot_series import SeriesStyle,SeriesStyleRegistry
from scripts.plotting.plot_style import load_profile,style_context,configure_axis_ticks
from scripts.plotting.validate_plot_artifact import validate_manifest

SHARED=("plot_style.py","plot_series.py","plot_manifest.py","plot_quality.py","plot_accessibility.py",
    "plot_bundle.py","plot_provenance.py","plot_delivery.py","validate_plot_artifact.py")
STYLES=SeriesStyleRegistry("sigma_mott_side_limits_v1",(
    SeriesStyle("left","#0072B2","-","o"),
    SeriesStyle("right","#D55E00","--","s"),
    SeriesStyle("half_power","#666666",":",""),
    SeriesStyle("base_coefficient","#0072B2","-","D"),
    SeriesStyle("smaller_step","#008565","None","s"),
    SeriesStyle("fine_background","#000000","None","x"),
    SeriesStyle("threshold_coefficient","#666666","--",""),
))
LABELS={"left":r"Resonant side: $\xi=\xi_M-h$","right":r"Bound side: $\xi=\xi_M+h$",
    "half_power":r"Reference exponent $1/2$","base_coefficient":"Continuous integral",
    "smaller_step":"Smaller derivative steps","fine_background":"Finer background + smaller steps",
    "threshold_coefficient":"Independent threshold coefficient"}


def read_csv(path):
    with path.open(encoding="utf-8-sig",newline="") as stream:
        return list(csv.DictReader(stream))


def draw(kind,data,args,profile,inputs):
    size=(6.75,7.0 if kind=="coefficients" else 6.8)
    selected=(["base_coefficient","smaller_step","fine_background","threshold_coefficient"] if kind=="coefficients" else
        ["left","right"]+(["half_power"] if kind=="scaling" else []))
    styles={s.series_id:s for s in STYLES.select(selected)}
    labels=[LABELS[s] for s in selected]
    if kind=="coefficients":
        fields=["scaled_curvature","log_coefficient_estimate"]
        ylabels=[r"$\sqrt{h}\,I_S^{\prime\prime}(\xi_M-h)\;(\mathrm{fm}^{2})$",
                 r"$Q_{\log}(h)\;(\mathrm{fm}^{2})$"]
        title=r"$\sigma$: coefficient limits on both sides"
    elif kind=="curvature":
        fields=["curvature","curvature"]
        ylabels=[r"$I_S^{\prime\prime}(\xi_M-h)\;(\mathrm{fm}^{2})$",
                 r"$I_S^{\prime\prime}(\xi_M+h)\;(\mathrm{fm}^{2})$"]
        title=r"$\sigma$: separate side limits of the threshold integral"
    else:
        fields=["scaled_curvature","effective_exponent"]
        ylabels=[r"$\sqrt{h}\,I_S^{\prime\prime}\;(\mathrm{fm}^{2})$",r"$p_{\rm eff}$"]
        title=r"$\sigma$: test of the left-side $h^{-1/2}$ law"
    counts={s:0 for s in selected}
    with style_context(profile) as context:
        fig,axes=plt.subplots(2,1,figsize=size,sharex=True)
        context.bind(fig)
        fig.subplots_adjust(left=.19,right=.965,bottom=.105,top=.64 if kind=="coefficients" else .69,hspace=.23)
        for index,(ax,field) in enumerate(zip(axes,fields)):
            keys=selected if kind=="coefficients" else ["left" if index==0 else "right"] if kind=="curvature" else ["left","right"]
            for key in keys:
                if kind=="coefficients":
                    variants={"base_coefficient":"original","smaller_step":"same_background_small_step","fine_background":"fine_background"}
                    if key=="threshold_coefficient":
                        domain=[min(float(r["h"]) for r in data),max(float(r["h"]) for r in data)]
                        value=args.coefficient_meta["independent_left_coefficient" if index==0 else "independent_right_log_coefficient"]
                        rows=[dict(h=h,**{field:value}) for h in domain]
                    elif index==0:
                        source=data if key=="base_coefficient" else args.reference_curvatures
                        rows=[r for r in source if r["side"]=="left" and (key=="base_coefficient" or r["variant"]==variants[key])]
                    else:
                        rows=[r for r in args.log_coefficients if r["variant"]==variants[key]]
                    rows=sorted(rows,key=lambda r:float(r["h"]))
                else:
                    rows=sorted([r for r in data if r["side"]==key and r.get(field,"")!=""],key=lambda r:float(r["h"]))
                style=styles[key]
                marker_size=7.4 if key in ("smaller_step","fine_background") else 5.2
                ax.plot([float(r["h"]) for r in rows],[float(r[field]) for r in rows],
                    **style.plot_kwargs(),linewidth=1.2,markersize=marker_size,markerfacecolor="white",
                    markeredgewidth=1.4 if key in ("smaller_step","fine_background") else 1.0)
                counts[key]+=len(rows)
            ax.set_xscale("log")
            if kind=="curvature" and index==0:
                ax.set_yscale("log")
            else:
                ax.yaxis.set_major_locator(MaxNLocator(nbins=4))
                ax.ticklabel_format(axis="y",style="plain",useOffset=False)
            configure_axis_ticks(ax,profile)
            ax.xaxis.set_minor_locator(LogLocator(base=10,subs=range(2,10),numticks=100))
            if ax.get_yscale()=="log":
                ax.yaxis.set_minor_locator(LogLocator(base=10,subs=range(2,10),numticks=100))
                ax.yaxis.set_minor_formatter(NullFormatter())
            ax.set_ylabel(ylabels[index],labelpad=10)
            if kind=="coefficients":
                ax.margins(x=.09,y=.20)
        if kind=="scaling":
            domain=[min(float(r["h"]) for r in data),max(float(r["h"]) for r in data)]
            axes[1].plot(domain,[.5,.5],**styles["half_power"].plot_kwargs(),linewidth=1.2)
            counts["half_power"]=2
            axes[1].set_ylim(0.,.62)
        axes[-1].set_xlabel(r"$h=|\xi-\xi_M|$ (Mott limit toward left)")
        fig.text(.55,.98,title,ha="center",va="top",fontsize=13)
        fig.text(.55,.925,r"$T=125.736\;\mathrm{MeV},\quad\mu_B=900\;\mathrm{MeV}$",ha="center",va="top",fontsize=13)
        handles=[Line2D([],[],**styles[s].plot_kwargs(),
            markersize=7.4 if s in ("smaller_step","fine_background") else 5.2,markerfacecolor="white",
            markeredgewidth=1.4 if s in ("smaller_step","fine_background") else 1.0) for s in selected]
        fig.legend(handles,labels,loc="upper center",bbox_to_anchor=(.55,.87),ncol=1,
            handlelength=2.1,labelspacing=.30,borderaxespad=0,frameon=False)
        stem=args.output/("sigma_side_"+kind)
        outputs,quality=export_figure(fig,stem,profile,formats=["png"],intended_width_inches=size[0])
    color=stem.with_suffix(".png");gray=stem.with_name(stem.name+"_gray").with_suffix(".png")
    with Image.open(color) as im:
        im.convert("L").convert("RGB").save(gray,dpi=(profile.dpi,profile.dpi))
    gr=output_record(gray,fmt="png",dpi=profile.dpi,vector=False);gr["inspection"]=inspect_export(gray);outputs.append(gr)
    manifest=build_manifest(asset_id="sigma_side_"+kind,figure_family="sigma_mott_asymptotics",case_slug=args.output.name,
        figure_mode="audit",semantic_status="conditional_continuum_threshold_diagnostic",style_profile=profile.profile_id,
        publication_scope="internal_review",generator=generator_record(__file__,command=" ".join(sys.argv),
            source_snapshot=args.results/"source"/Path(__file__).name,runtime=runtime_record({"matplotlib":matplotlib.__version__})),
        inputs=inputs,axes=[dict(field="h",source_unit="dimensionless",display_unit="dimensionless",scale="log",label=r"$h=|\xi-\xi_M|$"),
            *[dict(field=field,source_unit="dimensionless" if field=="effective_exponent" else "fm^2",
                display_unit="dimensionless" if field=="effective_exponent" else "fm^2",
                scale="log" if kind=="curvature" and i==0 else "linear",label=label)
              for i,(field,label) in enumerate(zip(fields,ylabels))]],
        series=[dict(series_id=s,state="analytical_reference" if s in ("half_power","threshold_coefficient") else "numerical_diagnostic",
            support_rule="Independent threshold coefficient, no fit" if s=="threshold_coefficient" else "Constant one-half reference, no fit" if s=="half_power" else "Precomputed remote side-local curvature",
            mask_rule="first effective exponent absent by definition" if kind=="scaling" else "none",point_count=counts[s]) for s in selected],
        outputs=outputs,selection_rule="Both sides of the same nondegenerate sigma threshold in the continuum local mixed-square integral",
        interpolation_policy="none",connector_policy="forbidden",missing_value_policy="no imputation",
        validation=dict(finite=True,duplicate_keys=True,support=True,solver_called=False,canonical_data_modified=False,author_review="pending"),
        calculation_sha="2e5c411c716d5cbc524c3fc6cfb335fe05b39444",source_run_id="trsigma-mott-20261101-v1/asym_full_v2",
        rendering=dict(delivery_stage="png_review",vector_delivery_pending=True,review_contract="prd_review_v1",figure_size_inches=list(size),
            column="double_column",font=dict(context.font),size_override_reason="Two separate side-limit diagnostics with readable mathematical labels",
            typography_exception="dense_composite_review_compact_typography",typography_exception_reason="Internal review with nested derivatives; enforce minimum 1.5 mm",
            series_styles=STYLES.record(selected),legend_policy="declared_geometry_checked",legend_outside=True,
            legend_placements=[dict(legend_index=0,host_axes_index=None,title="",labels=labels,placement="outside_axes",location="upper center",
                scope="Coefficient figure: left scaled curvature above, right logarithmic increment below; sparse references retain their actual previous h. Curvature: left above, right below; scaling: both sides in both panels")],
            legend_layout_reason="Shared key avoids obscuring the asymptotic platform",color_route="undecided_review",quality=quality,
            placement_limits=placement_limits(quality,profile,outputs),
            accessibility=accessibility_record([dict(value=s,color=styles[s].color,
                noncolor_encoding=f"{styles[s].linestyle}; {styles[s].marker or 'reference line'}") for s in selected],
                backgrounds=["#FFFFFF"],background_scope="white axes without fills"),
            coordinate_transform="Read-only display of remote curvature, sqrt(h) times curvature and successive ratio exponents",
            connections="Straight guides through computed samples; no fit or inferred point at h=0",
            grayscale=dict(method="Pillow L then RGB",color_sha256=sha256_file(color),gray_sha256=sha256_file(gray))))
    manifest.update(delivery_stage="png_review",manuscript_eligible=False,current_publication_layer=False,vector_delivery_pending=True,
        caption=input_record(args.results/"caption.md",role="figure_caption"))
    plt.close(fig)
    return manifest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("review",type=Path);parser.add_argument("output",type=Path);parser.add_argument("results",type=Path)
    parser.add_argument("--figures",default="curvature,scaling")
    parser.add_argument("--background-review",type=Path)
    args=parser.parse_args()
    kinds=args.figures.split(",")
    if not set(kinds)<={"curvature","scaling","coefficients"} or len(kinds)!=len(set(kinds)):
        raise ValueError("Invalid or duplicate figure selection")
    summary=json.loads((args.review/"summary.json").read_text(encoding="utf-8-sig"))
    for name,expected in summary["files"].items():
        if sha256_file(args.review/name)!=expected:
            raise ValueError("Frozen input changed: "+name)
    data=read_csv(args.review/"side_curvatures.csv")
    extra_inputs=[]
    if "coefficients" in kinds:
        if args.background_review is None:raise ValueError("Coefficient plot requires completed background review")
        args.coefficient_meta=json.loads((args.background_review/"summary.json").read_text(encoding="utf-8-sig"))
        for name,expected in args.coefficient_meta["files"].items():
            if sha256_file(args.background_review/name)!=expected:raise ValueError("Frozen reference changed: "+name)
        args.reference_curvatures=read_csv(args.background_review/"reference_curvatures.csv")
        args.log_coefficients=read_csv(args.background_review/"right_log_coefficient.csv")
        extra_inputs=[input_record(args.background_review/name,role="independent_numerical_reference") for name in
            ("summary.json","reference_curvatures.csv","right_log_coefficient.csv")]
    args.output.mkdir(parents=True,exist_ok=False);args.results.mkdir(parents=True,exist_ok=False)
    (args.results/"source").mkdir()
    shutil.copyfile(__file__,args.results/"source"/Path(__file__).name)
    caption=("Sigma side-limit diagnostic. Mode A, T=125.7360002645551 MeV, mu_B=900 MeV. "
        f"The refined Mott coordinate is xi_M={summary['xi_mott']}. I_S is the continuous thermal integral "
        "of the positive full sigma/sigma-prime flavor-propagator square over 0<x=s-4M_l^2<0.01 fm^-2 "
        "in uubar to uubar (the same light square occurs in uubar to ddbar). It includes both mixed poles, "
        "but excludes interference with other amplitudes and the rest of the collision rate. I_S is in fm^2. "
        "The numerical determinant safeguard is omitted in this isolated continuous-formula diagnostic; "
        "production kernels are unchanged. The actual integral is evaluated at 110 solved backgrounds, "
        "using a 192-bit local-kernel representation and a resonance-centered quadrature. No h exponent "
        "is built into the rates. At each u=+/-h, central differences with steps h/4 and h/8 stay entirely "
        "on that side; the displayed derivative is (4 D_h8-D_h4)/3. The normalized plot displays "
        "sqrt(h) I_S'', and p_eff=ln(|I_S''(h)|/|I_S''(2h)|)/ln(2). A half-power law gives a nonzero "
        "platform and p_eff approaching 1/2. First exponents are absent by definition. "
        "Numerical error and independent reference checks are retained alongside the remote tables. "
        "These are local-integral derivatives, not complete tau_ubar or terminal transport derivatives. "
        "Lines only guide the supplied points; no fit, smoothing, or extrapolation to h=0. "
        "Minimum insertion width 6.75 inches. PNG diagnostic review only.\n")
    if "coefficients" in kinds:
        caption+=("The coefficient figure has left-side sqrt(h) I_S'' in the upper panel and right-side "
            "Q_log(h)=[I_S''(h)-I_S''(h_previous)]/ln(h_previous/h) in the lower panel. "
            "The original curve uses h_previous=2h; sparse checks use the stored actual previous h. "
            "Smaller-step checks use steps h/8 and h/16; finer backgrounds use 128/32 background nodes, "
            "64/32 A-field nodes and their own relocated Mott center. Reference lines come from the "
            "independently evaluated threshold kernel and weight, not a fit to curvature. "
            "The finest right-side interval is retained despite its greater numerical sensitivity. "
            "Overlapping reference markers represent close actual values and are not displaced.\n")
    (args.results/"caption.md").write_text(caption,encoding="utf-8")
    inputs=[input_record(args.review/name,role="frozen_remote_diagnostic") for name in ("summary.json","side_curvatures.csv")]+extra_inputs
    for path in [ROOT/"config/plotting/candidate_aps_v2.toml",*[ROOT/"scripts/plotting"/name for name in SHARED]]:
        target=args.results/"source"/path.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(path,target)
        rec=input_record(path,role="plot_style_contract" if path.suffix==".toml" else "shared_plotting_code")
        rec["source_snapshot"]=input_record(target,role="plotting_source_snapshot");inputs.append(rec)
    profile=load_profile("candidate_aps_v2")
    figures=[draw(kind,data,args,profile,inputs) for kind in kinds]
    bundle=build_bundle(dict(case_slug=args.output.name,delivery_stage="png_review",manuscript_eligible=False,
        current_publication_layer=False,charts=[dict(asset_id=r["asset_id"]) for r in figures]),figures)
    path=args.output/"plot_manifest.json";write_manifest(path,bundle)
    errors=validate_manifest(path,repo_root=ROOT)
    (args.results/"validation.json").write_text(json.dumps(dict(passed=not errors,errors=errors),indent=2)+"\n",encoding="utf-8")
    if errors:raise ValueError("; ".join(errors))
    print(json.dumps(dict(passed=True,figures=len(figures),manifest=str(path))))


if __name__=="__main__":main()
