"""Render frozen sampled parameter envelopes and selected direct q0 curves."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from scripts.analysis.relaxtime.plot_freezeout_comparison import (
    load_experiment, experiment_provenance, experiment_groups, prepare_review_directories,
    export_review_png, write_review_manifest, experimental_caption,
)
from scripts.plotting.plot_manifest import (build_manifest, generator_record, input_record,
    output_record, runtime_record, sha256_file)
from scripts.plotting.plot_quality import inspect_export, placement_limits
from scripts.plotting.plot_style import load_profile, configure_matplotlib, configure_axis_ticks, figure_size_inches
from scripts.plotting.plot_accessibility import accessibility_record
from scripts.plotting.validate_plot_artifact import validate_manifest


def table(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def number(row, key):
    return float(row[key]) if row[key] else math.nan


def retained_data_path(path):
    """Locate a retained data file after its originating worktree is retired."""
    original = Path(path)
    if original.is_file():
        return original
    parts = str(path).replace("\\", "/").split("/")
    if "data" not in parts:
        raise FileNotFoundError(f"no retained data path: {path}")
    relative = Path(*parts[parts.index("data"):])
    current = (ROOT / relative).resolve()
    if ".." in relative.parts or not current.is_relative_to(ROOT / "data"):
        raise ValueError(f"invalid retained data path: {path}")
    return current


def display_experiments(summary, override=None):
    """Select a display table while keeping the frozen scoring input intact."""
    scoring_path = retained_data_path(summary["experiment"]["path"])
    if sha256_file(scoring_path) != summary["experiment"]["sha256"]:
        raise ValueError("experimental scoring input changed")
    display_path = Path(override).resolve() if override is not None else scoring_path
    rows = [r for r in load_experiment(display_path)
        if 2.35 <= r["energy"] <= 200 and r["method"] != "digitized"
        and r["system"] in {"Au+Au", "Pb+Pb"}
        and r["experiment"] != "ALICE"]
    return display_path, scoring_path, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--experiment", type=Path,
        help="optional display-only experimental table; frozen scores and bands are unchanged")
    parser.add_argument("--replace-draft", action="store_true",
        help="replace the unaccepted PNG draft in place; no directories or unrelated files are deleted")
    args = parser.parse_args()
    source, output = args.summary_dir.resolve(), args.output_dir.resolve()
    results = ROOT / "data/outputs/results/relaxtime/meson_density" / output.name
    for name, expected in json.loads((source / "checksums.json").read_text(encoding="utf-8")).items():
        if sha256_file(source / name) != expected:
            raise ValueError(f"summary input changed: {name}")
    summary = json.loads((source / "summary.json").read_text(encoding="utf-8"))
    is_andronic = summary.get("family") == "andronic_logistic_rational_v1"
    experiment_path, scoring_path, experiments = display_experiments(summary, args.experiment)
    screen_path = retained_data_path(summary["screen_assessment"]["path"])
    if sha256_file(screen_path) != summary["screen_assessment"]["sha256"]:
        raise ValueError("screening membership input changed")
    screen = json.loads(screen_path.read_text(encoding="utf-8"))
    bands = table(source / "screen_band.csv")
    fine = [r for r in table(source / "fine_curves.csv") if r["profile"] == summary["best_selected"]]
    baseline = [r for r in table(source / "baseline_curves.csv") if r["profile"] == "andronic2018"]
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FuncFormatter, MultipleLocator
    from PIL import Image

    profile = load_profile("candidate_aps_v2")
    font = configure_matplotlib(profile)
    figure, ax = plt.subplots(figsize=(figure_size_inches(profile)[0], 6.0))
    figure.subplots_adjust(left=.115, right=.975, bottom=.16, top=.96)
    x = [number(r, "energy_GeV") for r in bands]
    ax.fill_between(x, [number(r,"minimum") for r in bands], [number(r,"maximum") for r in bands],
        facecolor="#E8E8E8", edgecolor="#595959", linewidth=.65, linestyle="--", label="Sampled envelope", zorder=1)
    near_count = len(screen["near_best_band_complete_profiles"])
    if near_count > 1:
        ax.fill_between(x, [number(r,"near_minimum") for r in bands], [number(r,"near_maximum") for r in bands],
            facecolor="#FFFFFF", edgecolor="#0072B2", linewidth=.8, hatch="///", label="Near-best envelope", zorder=2)
    for rows, color, style, label in ((fine,"#0072B2","-","Best sampled (refined)"),
        (baseline,"#A04800","--","Andronic (2018)")):
        ax.plot([number(r,"energy_GeV") for r in rows], [number(r,"ratio") for r in rows],
            color=color, linestyle=style, linewidth=1.7, label=label, zorder=3)
    groups = experiment_groups(experiments)
    for group in groups:
        rows, marker, color, label, filled = (group[key] for key in
            ("rows", "marker", "color", "label", "filled"))
        ax.errorbar([r["energy"] for r in rows], [r["ratio"] for r in rows],
            yerr=[[r["error_minus"] for r in rows],[r["error_plus"] for r in rows]],
            fmt=marker, color=color, markerfacecolor=color if filled else "white",
            markersize=6 if group["name"] == "NA61" else 5, markeredgewidth=.8, elinewidth=.7, capsize=2, capthick=.7,
            linestyle="none", label=label, zorder=5)
    valid = [number(r,"maximum") for r in bands if r["maximum"]]
    valid += [number(r,"ratio") for r in fine+baseline if r["ratio"]]
    ax.set_xscale("log")
    ax.set_xlim(2.1, 230)
    ax.set_ylim(0, max(.88, math.ceil((max(valid)+.10)/.05)*.05))
    ax.xaxis.set_major_locator(FixedLocator([3,5,10,20,50,100,200]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, pos: f"{value:g}"))
    ax.yaxis.set_major_locator(MultipleLocator(.1))
    configure_axis_ticks(ax, profile)
    xlabel, ylabel = r"$\sqrt{s_{NN}}\;(\mathrm{GeV})$", r"$K^+/\pi^+$"
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel, labelpad=8)
    title = "PNJL + GBU\nAndronic freezeout" if is_andronic else "PNJL + GBU\nChemical freezeout"
    ax.text(.026,.962,title,transform=ax.transAxes,va="top",fontsize=12)
    ax.legend(loc="upper right",fontsize=12,handlelength=2.3,labelspacing=.32,borderaxespad=.5,frameon=False)
    labels = ax.get_legend_handles_labels()[1]
    prepare_review_directories(output, results, replace_draft=args.replace_draft)
    outputs, quality = export_review_png(figure, output / "kplus_piplus_parameter_envelope", profile)
    gray = output / "kplus_piplus_parameter_envelope_gray.png"
    with Image.open(output / "kplus_piplus_parameter_envelope.png") as im:
        im.convert("L").convert("RGB").save(gray,dpi=(profile.dpi,profile.dpi))
    gray_record = output_record(gray,fmt="png",dpi=profile.dpi,vector=False)
    gray_record["inspection"] = inspect_export(gray)
    outputs.append(gray_record)
    plt.close(figure)
    caption = ("Cross-literature chemical-freezeout parameter exploration with fixed quark-only BQS and q0 Lambda extrapolation. "
        f"Gray: min/max of {len(screen['band_complete_profiles'])} complete screening-resolution curves on the 17-point grid. "
        f"Hatched blue, when present: {near_count} complete curves with the predeclared nine-point Au+Au RMSE within 10% of the sampled minimum. "
        f"Solid blue: {summary['best_selected']}, evaluated on the refined energy grid at fine numerical resolution. "
        "Orange dashed: frozen Andronic 2018 default at fine resolution. "
        + experimental_caption() +
        "Only nine Au+Au points at 7.7-200 GeV enter the primary RMSE; the 4.87 GeV point and Pb+Pb points are shape diagnostics. "
        "Source centralities and rapidity windows are retained in the versioned experimental input. "
        "New observations in a display-only table do not enter the frozen scores or change band membership. "
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
        *[dict(value=g["name"],color=g["color"],
            noncolor_encoding=f"{g['marker']}; {'filled' if g['filled'] else 'open'}") for g in groups]]
    manifest = build_manifest(asset_id="meson_freezeout_"+("andronic_" if is_andronic else "")+"parameter_envelope_20261018",figure_family="charged_meson_freezeout",
        case_slug=output.name,figure_mode="audit",semantic_status="cross_literature_parameter_exploration",
        style_profile=profile.profile_id,publication_scope="internal_review",
        generator=generator_record(__file__,command=" ".join(sys.argv),
            runtime=runtime_record({"matplotlib":matplotlib.__version__})),
        inputs=[*[input_record(source / n,role="frozen_exploration") for n in
            ("summary.json","checksums.json","screen_band.csv","fine_curves.csv","baseline_curves.csv","selected_parameters.csv")],
            input_record(screen_path,role="screening_membership_and_selection"),
            input_record(experiment_path,role="experimental_observations"),
            input_record(scoring_path,role="frozen_scoring_observations"),
            experiment_provenance(experiment_path),
            input_record(ROOT / "scripts/analysis/relaxtime/plot_freezeout_comparison.py",role="shared_plotting_code"),
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
            dict(series_id="andronic2018",state="published_default_reference",point_count=len(baseline),
                support_rule="frozen original fine-grid direct q0 calculations",
                mask_rule="NaN at failed support; no bridging",
                failure_count=sum(not r["ratio"] for r in baseline)),
            *[dict(series_id=g["name"],state=g["state"],point_count=len(g["rows"]),
                support_rule="retained source observations; separate systems/charges/acceptances",
                mask_rule="none; retained source observations only",
                rapidity=sorted({r["rapidity"] for r in g["rows"]}),
                error_methods=sorted({r["error_method"] for r in g["rows"]}),
                ids=[r["id"] for r in g["rows"]]) for g in groups]],
        outputs=outputs,selection_rule="2.35-200 GeV observations; frozen fixed sample envelope and selected refined curve on original support",
        interpolation_policy="none",connector_policy="forbidden",missing_value_policy="gaps retained; no zero fill",
        validation=dict(finite=True,duplicate_keys=True,support=True,solver_called=False,canonical_data_modified=False,
            checks="frozen hashes, fixed membership, source semantics, measured rendering",author_review="pending"),
        rendering=dict(delivery_stage="png_review",vector_delivery_pending=True,review_contract="prd_review_v1",
            figure_size_inches=list(figure.get_size_inches()),size_override_reason="Additional acceptance and preliminary-data groups need an in-panel legend above the high-energy plateau.",column="double_column",font=font,
            legend_policy="declared_geometry_checked",legend_outside=False,
            legend_placements=[dict(legend_index=0,host_axes_index=0,title="",labels=labels,placement="in_axes",location="upper right",scope="bands, curves and experiments")],
            legend_layout_reason="author requests legends inside one panel; upper-right whitespace retained",
            color_route="undecided_review",quality=quality,placement_limits=placement_limits(quality,profile,outputs),
            accessibility=accessibility_record(encodings,backgrounds=["#FFFFFF","#E8E8E8"],background_scope="white axes and gray full envelope; dense blue hatching on white identifies the near-best envelope"),
            connections="straight segments at existing direct support; no additional points"),
        calculation_sha=summary["calculation_sha"],source_run_id=summary["fine_inputs"][0]["scan_identity"])
    manifest.update(delivery_stage="png_review",manuscript_eligible=False,current_publication_layer=False,
        vector_delivery_pending=True,caption=input_record(results / "caption.md",role="figure_caption"))
    write_review_manifest(output / "plot_manifest.json",manifest)
    issues = validate_manifest(output / "plot_manifest.json",repo_root=ROOT)
    if issues:
        raise ValueError("plot contract failed: " + "; ".join(issues))
    print(json.dumps(dict(figure=str(output / "kplus_piplus_parameter_envelope.png"),validation="passed")))


if __name__ == "__main__":
    main()
