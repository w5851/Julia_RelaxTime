"""Render frozen sampled parameter envelopes and selected direct q0 curves."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from scripts.analysis.relaxtime.plot_freezeout_comparison import load_experiment
from scripts.plotting.plot_manifest import (build_manifest, generator_record, input_record,
    output_record, runtime_record, sha256_file, write_manifest)
from scripts.plotting.plot_quality import export_figure, inspect_export, placement_limits
from scripts.plotting.plot_style import load_profile, configure_matplotlib, configure_axis_ticks, figure_size_inches
from scripts.plotting.plot_accessibility import accessibility_record
from scripts.plotting.validate_plot_artifact import validate_manifest


def table(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def number(row, key):
    return float(row[key]) if row[key] else math.nan


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    source, output = args.summary_dir.resolve(), args.output_dir.resolve()
    results = ROOT / "data/outputs/results/relaxtime/meson_density" / output.name
    if output.exists() or results.exists():
        raise ValueError("plot case exists; choose a new sibling case")
    for name, expected in json.loads((source / "checksums.json").read_text(encoding="utf-8")).items():
        if sha256_file(source / name) != expected:
            raise ValueError(f"summary input changed: {name}")
    summary = json.loads((source / "summary.json").read_text(encoding="utf-8"))
    experiment_path = Path(summary["experiment"]["path"])
    if sha256_file(experiment_path) != summary["experiment"]["sha256"]:
        raise ValueError("experimental input changed")
    screen_path = Path(summary["screen_assessment"]["path"])
    if sha256_file(screen_path) != summary["screen_assessment"]["sha256"]:
        raise ValueError("screening membership input changed")
    screen = json.loads(screen_path.read_text(encoding="utf-8"))
    bands = table(source / "screen_band.csv")
    fine = [r for r in table(source / "fine_curves.csv") if r["profile"] == summary["best_selected"]]
    baseline = [r for r in table(source / "baseline_curves.csv") if r["profile"] == "andronic2018"]
    experiments = [r for r in load_experiment(experiment_path)
        if 4.87 <= r["energy"] <= 200 and r["method"] != "digitized"]
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FuncFormatter, MultipleLocator
    from PIL import Image

    profile = load_profile("candidate_aps_v2")
    font = configure_matplotlib(profile)
    figure, ax = plt.subplots(figsize=figure_size_inches(profile))
    figure.subplots_adjust(left=.115, right=.975, bottom=.16, top=.96)
    x = [number(r, "energy_GeV") for r in bands]
    ax.fill_between(x, [number(r,"minimum") for r in bands], [number(r,"maximum") for r in bands],
        facecolor="#E8E8E8", edgecolor="#595959", linewidth=.65, linestyle="--", label="Sampled envelope", zorder=1)
    near_count = len(screen["near_best_band_complete_profiles"])
    if near_count > 1:
        ax.fill_between(x, [number(r,"near_minimum") for r in bands], [number(r,"near_maximum") for r in bands],
            facecolor="#C6DEF0", edgecolor="#0072B2", linewidth=.65, hatch="/", label="Near-best envelope", zorder=2)
    for rows, color, style, label in ((fine,"#0072B2","-","Best sampled (refined)"),
        (baseline,"#A04800","--","Andronic (2018)")):
        ax.plot([number(r,"energy_GeV") for r in rows], [number(r,"ratio") for r in rows],
            color=color, linestyle=style, linewidth=1.7, label=label, zorder=3)
    for system, marker, color, label in (("Au+Au","o","#B2282E","Au+Au: AGS, STAR"),
        ("Pb+Pb","s","#222222","Pb+Pb: NA49")):
        rows = [r for r in experiments if r["system"] == system and r["method"] != "digitized"]
        ax.errorbar([r["energy"] for r in rows], [r["ratio"] for r in rows],
            yerr=[[r["error_minus"] for r in rows],[r["error_plus"] for r in rows]],
            fmt=marker, color=color, markerfacecolor=color if system == "Au+Au" else "white",
            markersize=5, markeredgewidth=.8, elinewidth=.7, capsize=2, capthick=.7,
            linestyle="none", label=label, zorder=5)
    valid = [number(r,"maximum") for r in bands if r["maximum"]]
    valid += [number(r,"ratio") for r in fine+baseline if r["ratio"]]
    ax.set_xscale("log")
    ax.set_xlim(4.3, 230)
    ax.set_ylim(0, max(.65, math.ceil((max(valid)+.10)/.05)*.05))
    ax.xaxis.set_major_locator(FixedLocator([5,10,20,50,100,200]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, pos: f"{value:g}"))
    ax.yaxis.set_major_locator(MultipleLocator(.1))
    configure_axis_ticks(ax, profile)
    xlabel, ylabel = r"$\sqrt{s_{NN}}\;(\mathrm{GeV})$", r"$K^+/\pi^+$"
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel, labelpad=8)
    ax.text(.026,.962,"PNJL + GBU\nChemical freezeout",transform=ax.transAxes,va="top",fontsize=12)
    ax.legend(loc="upper right",fontsize=12,handlelength=2.7,labelspacing=.4,borderaxespad=.5,frameon=False)
    labels = ax.get_legend_handles_labels()[1]
    output.mkdir(parents=True)
    results.mkdir(parents=True)
    outputs, quality = export_figure(figure, output / "kplus_piplus_parameter_envelope", profile, formats=["png"])
    gray = output / "kplus_piplus_parameter_envelope_gray.png"
    with Image.open(output / "kplus_piplus_parameter_envelope.png") as im:
        im.convert("L").convert("RGB").save(gray,dpi=(profile.dpi,profile.dpi))
    gray_record = output_record(gray,fmt="png",dpi=profile.dpi,vector=False)
    gray_record["inspection"] = inspect_export(gray)
    outputs.append(gray_record)
    plt.close(figure)
    snapshot = results / Path(__file__).name
    shutil.copy2(__file__, snapshot)
    caption = ("Cross-literature chemical-freezeout parameter exploration with fixed quark-only BQS and q0 Lambda extrapolation. "
        f"Gray: min/max of {len(screen['band_complete_profiles'])} complete screening-resolution curves on the 17-point grid. "
        f"Hatched blue, when present: {near_count} complete curves with the predeclared nine-point Au+Au RMSE within 10% of the sampled minimum. "
        f"Solid blue: {summary['best_selected']}, evaluated on the refined energy grid at fine numerical resolution. "
        "Orange dashed: frozen Andronic 2018 default at fine resolution. Filled circles: Au+Au AGS/STAR; open squares: auxiliary NA49 Pb+Pb. "
        "Only nine Au+Au points at 7.7-200 GeV enter the primary RMSE; the 4.87 GeV point and Pb+Pb points are shape diagnostics. "
        "Source centralities and rapidity windows are retained in the versioned experimental input. "
        "The envelopes describe finite sampled parameter sets, not probability intervals or a global optimum. "
        "Straight boundary segments join existing computed support; no added interpolation points or smoothing. "
        "Failed candidates excluded from whole-curve envelopes are identified in the assessment; failures in the best curve remain gaps. "
        "Fine convergence was checked only on the selected curves, not uniformly over the full screening ensemble. "
        "Finite Bose window, no decay or rapidity-acceptance model. PNG review only.\n")
    (results / "caption.md").write_text(caption, encoding="utf-8")
    encodings = [dict(value="sampled envelope",color="#595959",noncolor_encoding="dashed boundary"),
        dict(value="near-best envelope",color="#0072B2",noncolor_encoding="hatched fill"),
        dict(value="best sampled",color="#0072B2",noncolor_encoding="solid"),
        dict(value="Andronic",color="#A04800",noncolor_encoding="dashed"),
        dict(value="AuAu",color="#B2282E",noncolor_encoding="filled circles"),
        dict(value="NA49",color="#222222",noncolor_encoding="open squares")]
    manifest = build_manifest(asset_id="meson_freezeout_parameter_envelope_20261018",figure_family="charged_meson_freezeout",
        case_slug=output.name,figure_mode="audit",semantic_status="cross_literature_parameter_exploration",
        style_profile=profile.profile_id,publication_scope="internal_review",
        generator=generator_record(__file__,command=" ".join(sys.argv),source_snapshot=snapshot,
            runtime=runtime_record({"matplotlib":matplotlib.__version__})),
        inputs=[*[input_record(source / n,role="frozen_exploration") for n in
            ("summary.json","checksums.json","screen_band.csv","fine_curves.csv","baseline_curves.csv","selected_parameters.csv")],
            input_record(Path(summary["screen_assessment"]["path"]),role="screening_membership_and_selection"),
            input_record(experiment_path,role="experimental_observations"),
            input_record(experiment_path.parent / "sources.json",role="experimental_provenance"),
            input_record(ROOT / "config/plotting/candidate_aps_v2.toml",role="plot_style_contract"),
            *[input_record(ROOT / "scripts/plotting" / n,role="shared_plotting_code") for n in
                ("plot_style.py","plot_manifest.py","plot_quality.py","plot_accessibility.py","validate_plot_artifact.py")]],
        axes=[dict(field="sqrt_s_NN_GeV",source_unit="GeV",display_unit="GeV",scale="log",label=xlabel),
            dict(field="Kplus_over_pi_plus",source_unit="dimensionless",display_unit="dimensionless",scale="linear",label=ylabel)],
        series=[dict(series_id="screening_envelope",state="parameter_exploration",point_count=len(bands),
            member_count=len(screen["band_complete_profiles"]),support_rule="direct q0 support; fixed complete curves",
            mask_rule="incomplete trajectories excluded with named records; no zero filling"),
            *([dict(series_id="near_best_envelope",state="parameter_exploration",point_count=len(bands),
                member_count=near_count, members=screen["near_best_band_complete_profiles"],
                selection_rule="fixed nine-point Au+Au RMSE <= 1.10 times sampled minimum",
                support_rule="direct screening q0 support; fixed complete curves",
                mask_rule="incomplete trajectories excluded as whole curves")] if near_count > 1 else []),
            dict(series_id=summary["best_selected"],state="q0_reference_diagnostic",point_count=len(fine),
                support_rule="direct fine solves",mask_rule="NaN at failed support"),
            dict(series_id="andronic2018",state="published_default_reference",point_count=len(baseline)),
            *[dict(series_id=system,state="experimental",point_count=sum(r["system"] == system for r in experiments),
                ids=[r["id"] for r in experiments if r["system"] == system])
                for system in ("Au+Au", "Pb+Pb")]],
        outputs=outputs,selection_rule="4.87-200 GeV; frozen fixed sample envelope and selected refined curve",
        interpolation_policy="none",connector_policy="forbidden",missing_value_policy="gaps retained; no zero fill",
        validation=dict(finite=True,duplicate_keys=True,support=True,solver_called=False,canonical_data_modified=False,
            checks="frozen hashes, fixed membership, source semantics, measured rendering",author_review="pending"),
        rendering=dict(delivery_stage="png_review",vector_delivery_pending=True,review_contract="prd_review_v1",
            figure_size_inches=list(figure_size_inches(profile)),column="double_column",font=font,
            legend_policy="declared_geometry_checked",legend_outside=False,
            legend_placements=[dict(legend_index=0,host_axes_index=0,title="",labels=labels,placement="in_axes",location="upper right",scope="bands, curves and experiments")],
            legend_layout_reason="author requests legends inside one panel; upper-right whitespace retained",
            color_route="undecided_review",quality=quality,placement_limits=placement_limits(quality,profile,outputs),
            accessibility=accessibility_record(encodings,backgrounds=["#FFFFFF","#E8E8E8","#C6DEF0"],background_scope="white axes and light band fills; dark boundaries and hatching carry band semantics"),
            connections="straight segments at existing direct support; no additional points"),
        calculation_sha=summary["calculation_sha"],source_run_id=summary["fine_inputs"][0]["scan_identity"])
    manifest.update(delivery_stage="png_review",manuscript_eligible=False,current_publication_layer=False,
        vector_delivery_pending=True,caption=input_record(results / "caption.md",role="figure_caption"))
    write_manifest(output / "plot_manifest.json",manifest)
    issues = validate_manifest(output / "plot_manifest.json",repo_root=ROOT)
    if issues:
        raise ValueError("plot contract failed: " + "; ".join(issues))
    print(json.dumps(dict(figure=str(output / "kplus_piplus_parameter_envelope.png"),validation="passed")))


if __name__ == "__main__":
    main()
