"""Plot frozen q0 chemical-freezeout calculations and published experiments.

No solver calls, interpolation, parameter adjustment, or replacement of failures.
The primary comparison uses equal weight per published Au+Au energy, not per
theoretical support point. Charge-summed ALICE and digitized HADES stay auxiliary.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.plotting.plot_manifest import (
    build_manifest, generator_record, input_record, output_record,
    runtime_record, sha256_file, write_manifest,
)
from scripts.plotting.plot_quality import export_figure, inspect_export, placement_limits
from scripts.plotting.plot_style import load_profile, configure_matplotlib, configure_axis_ticks, figure_size_inches
from scripts.plotting.plot_accessibility import accessibility_record
from scripts.plotting.validate_plot_artifact import validate_manifest

PROFILE_LABELS = {"cleymans2006": "Cleymans (2006)", "andronic2018": "Andronic (2018)"}


def truth(value):
    return str(value).lower() == "true"


def checked_json(path):
    checksum = path.with_name(path.name + ".sha256")
    if checksum.read_text(encoding="utf-8").strip() != sha256_file(path):
        raise ValueError(f"checkpoint checksum mismatch: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def load_calculation(run_dir):
    manifest = checked_json(run_dir / "manifest.json")
    state = checked_json(run_dir / "run.json")
    if manifest["scan_identity"] != state["scan_identity"] or manifest["run_sha256"] != sha256_file(run_dir / "run.json"):
        raise ValueError("calculation run and manifest identities disagree")
    if manifest["density_route"] != "q0_lambda_reference":
        raise ValueError("this comparison requires the q0 reference route")
    for path, digest in manifest["outputs"].items():
        if sha256_file(run_dir / path) != digest:
            raise ValueError(f"frozen calculation changed: {path}")
    with (run_dir / "freezeout_ratios.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    result = {}
    for row in rows:
        profile = row["profile_name"]
        energy = float(row["sqrt_s_NN_GeV"])
        if energy in result.setdefault(profile, {}):
            raise ValueError(f"duplicate energy for {profile}: {energy}")
        passed = truth(row["plus_passed"])
        value = float(row["Kplus_over_pi_plus"] or "nan")
        if passed and not math.isfinite(value):
            raise ValueError("passed calculation has nonfinite ratio")
        if not passed and math.isfinite(value):
            raise ValueError("failed calculation must not supply a plotting value")
        result[profile][energy] = {**row, "energy": energy, "ratio": value, "passed": passed,
            "T_MeV": float(row["T_MeV"]), "muB_MeV": float(row["muB_MeV"])}
    if len(rows) != manifest["point_count"]:
        raise ValueError("calculation point count mismatch")
    return result, manifest, state


def load_experiment(path):
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len({r["id"] for r in rows}) != len(rows):
        raise ValueError("duplicate experimental observation")
    return [{**r, "energy": float(r["sqrt_s_NN_GeV"]), "ratio": float(r["ratio"]),
        "error_minus": float(r["error_minus"]), "error_plus": float(r["error_plus"]),
        "primary_auau": truth(r["primary_auau"]), "auxiliary_horn": truth(r["auxiliary_horn"])} for r in rows]


def residual_summary(pairs):
    residuals = [p[0] - p[1] for p in pairs]
    if not residuals:
        return {"count": 0, "rmse": None, "mae": None, "mean_signed_bias": None}
    return {"count": len(residuals), "rmse": math.sqrt(sum(x*x for x in residuals)/len(residuals)),
        "mae": sum(abs(x) for x in residuals)/len(residuals),
        "mean_signed_bias": sum(residuals)/len(residuals)}


def compare(models, experiments):
    profiles = list(models)
    if len(profiles) < 2:
        raise ValueError("at least two profiles are required for comparison")
    requested = [r for r in experiments if r["primary_auau"] or r["auxiliary_horn"]]
    common = []
    excluded = []
    for observation in requested:
        energy = observation["energy"]
        if any(energy not in models[p] for p in profiles):
            raise ValueError(f"calculate the experimental energy directly; no interpolation: {energy}")
        if all(models[p][energy]["passed"] for p in profiles):
            common.append(observation)
        else:
            excluded.append({"id": observation["id"], "energy": energy,
                "failed_profiles": [p for p in profiles if not models[p][energy]["passed"]]})
    primary = [r for r in common if r["primary_auau"]]
    if not primary:
        raise ValueError("no common valid Au+Au comparisons")
    summaries = {}
    residual_rows = []
    for profile, points in models.items():
        def score(predicate):
            return residual_summary([(points[r["energy"]]["ratio"], r["ratio"]) for r in common if predicate(r)])
        window = sorted((p for p in points.values() if p["passed"] and 4 <= p["energy"] <= 20), key=lambda p: p["energy"])
        if not window:
            raise ValueError(f"no valid horn-window support for {profile}")
        index = max(range(len(window)), key=lambda i: window[i]["ratio"])
        peak = window[index]
        bracket = [window[max(0,index-1)]["energy"], window[min(len(window)-1,index+1)]["energy"]]
        at = lambda energy: points[energy]["ratio"] if energy in points and points[energy]["passed"] else None
        summaries[profile] = {
            "AuAu_primary": score(lambda r: r["primary_auau"]),
            "AuAu_horn_4_to_20_GeV": score(lambda r: r["primary_auau"] and 4 <= r["energy"] <= 20),
            "AuAu_19p6_to_200_GeV": score(lambda r: r["primary_auau"] and r["energy"] >= 19.6),
            "NA49_auxiliary": score(lambda r: r["auxiliary_horn"]),
            "sampled_peak_4_to_20_GeV": {"energy_GeV": peak["energy"], "ratio": peak["ratio"],
                "neighbor_bracket_GeV": bracket, "at_window_edge": index in (0,len(window)-1),
                "interpretation": "maximum of directly evaluated support; not a fitted continuous peak"},
            "ratio_at_7p7": at(7.7), "ratio_at_200": at(200.),
            "NA49_7p62_minus_17p27": None if at(7.62) is None or at(17.27) is None else at(7.62)-at(17.27),
            "negative_passed_points": sum(p["passed"] and p["ratio"] < 0 for p in points.values()),
        }
        for r in common:
            point = points[r["energy"]]
            residual_rows.append({"profile_name": profile, "experiment_id": r["id"], "system": r["system"],
                "sqrt_s_NN_GeV": r["energy"], "T_MeV": point["T_MeV"], "muB_MeV": point["muB_MeV"],
                "model_ratio": point["ratio"], "experimental_ratio": r["ratio"],
                "residual": point["ratio"]-r["ratio"], "primary_auau": r["primary_auau"]})
    best = min(profiles, key=lambda p: summaries[p]["AuAu_primary"]["rmse"])
    metric = {"schema": "chemical_freezeout_default_comparison_v1", "primary_metric": "equal-weight AuAu RMSE",
        "selection_scope": "descriptive comparison/model selection on these observations; not independent validation or a fit",
        "parameter_adjustment": False, "common_primary_points": len(primary), "excluded_failed_observations": excluded,
        "best_default_on_primary_metric": best, "profiles": summaries,
        "experimental_features": {"STAR_7p7": next(r["ratio"] for r in experiments if r["id"] == "STAR_BES_7.7"),
            "STAR_200": next(r["ratio"] for r in experiments if r["id"] == "STAR_200"),
            "NA49_7p62_minus_17p27": next(r["ratio"] for r in experiments if r["id"] == "NA49_7.62")-
                next(r["ratio"] for r in experiments if r["id"] == "NA49_17.27")}}
    return metric, residual_rows


def render(models, experiments, output_dir, calculation, state, run_dir, experiment_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FuncFormatter, MultipleLocator
    from PIL import Image

    profile = load_profile("candidate_aps_v2")
    font = configure_matplotlib(profile)
    figure, ax = plt.subplots(figsize=figure_size_inches(profile))
    figure.subplots_adjust(left=.115, right=.975, bottom=.16, top=.96)
    colors = ["#0072B2", "#D55E00"]
    linestyles = ["-", "--"]
    encodings, series = [], []
    profiles = list(dict.fromkeys(item["profile_name"] for item in state["points"]))
    for index, name in enumerate(profiles):
        points = sorted(models[name].values(), key=lambda p: p["energy"])
        color, style = colors[index % len(colors)], linestyles[index % len(linestyles)]
        ax.plot([p["energy"] for p in points], [p["ratio"] if p["passed"] else math.nan for p in points],
            color=color, linestyle=style, linewidth=1.6, label=PROFILE_LABELS.get(name,name), zorder=2)
        encodings.append({"value": name, "color": color, "noncolor_encoding": style})
        series.append({"series_id": name, "state": "q0_reference_diagnostic", "support_rule": "direct energy solves",
            "mask_rule": "failed=NaN; disconnect at every failed support point", "point_count": len(points),
            "failure_count": sum(not p["passed"] for p in points),
            "endpoint_warning_points": sum(int(p["pi_plus_endpoint_warning_shells"])+int(p["K_plus_endpoint_warning_shells"])>0 for p in points)})
    groups = [
        ("AuAu", lambda r: r["system"] == "Au+Au" and r["method"] != "digitized", "o", "#B2282E", "Au+Au: AGS, STAR", True),
        ("NA49", lambda r: r["experiment"] == "NA49", "s", "#222222", "Pb+Pb: NA49", True),
        ("ALICE", lambda r: r["experiment"] == "ALICE", "v", "#222222", "Pb+Pb: ALICE (charge sum)", False),
        ("HADES", lambda r: r["experiment"] == "HADES", "+", "#B2282E", "Au+Au: HADES (digitized)", True),
    ]
    for name, select, marker, color, label, filled in groups:
        rows = [r for r in experiments if select(r)]
        if not rows:
            continue
        kwargs = dict(fmt=marker, color=color, markerfacecolor=color if filled else "white",
            markeredgecolor=color, markersize=5.2, markeredgewidth=.85, elinewidth=.65,
            capsize=2, capthick=.65, linestyle="none", label=label, zorder=4)
        errors = None if name == "HADES" else [[r["error_minus"] for r in rows], [r["error_plus"] for r in rows]]
        ax.errorbar([r["energy"] for r in rows], [r["ratio"] for r in rows], yerr=errors, **kwargs)
        series.append({"series_id": name, "state": "experimental" if name != "HADES" else "digitized_reference",
            "support_rule": "retained source observations; separate systems/charges", "mask_rule": "none",
            "point_count": len(rows), "ids": [r["id"] for r in rows]})
    valid_values = [p["ratio"] for points in models.values() for p in points.values() if p["passed"]]
    ax.set_xscale("log")
    ax.set_xlim(2.1, 7000.)
    # Empty space above the high-energy plateau carries the user-preferred
    # in-panel legend. No computed extrema are clipped to create this space.
    ax.set_ylim(min(0., min(valid_values)-.02), max(.56, math.ceil((max(valid_values)+.08)/.05)*.05))
    ax.xaxis.set_major_locator(FixedLocator([3, 10, 100, 1000]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:g}"))
    ax.yaxis.set_major_locator(MultipleLocator(.1))
    configure_axis_ticks(ax, profile)
    xlabel, ylabel = r"$\sqrt{s_{NN}}\;(\mathrm{GeV})$", r"$K^+/\pi^+$"
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel, labelpad=8)
    ax.text(.025, .962, "PNJL + GBU\nChemical freezeout", transform=ax.transAxes,
        va="top", fontsize=12)
    ax.legend(loc="upper right", fontsize=12, handlelength=2.8, labelspacing=.45,
        borderaxespad=.55, frameon=False)
    outputs, quality = export_figure(figure, output_dir / "kplus_piplus_freezeout_defaults", profile, formats=["png"])
    color_path = output_dir / "kplus_piplus_freezeout_defaults.png"
    gray_path = output_dir / "kplus_piplus_freezeout_defaults_gray.png"
    with Image.open(color_path) as im:
        im.convert("L").convert("RGB").save(gray_path, dpi=(profile.dpi,profile.dpi))
    gray = output_record(gray_path, fmt="png", dpi=profile.dpi, vector=False)
    gray["inspection"] = inspect_export(gray_path); outputs.append(gray)
    plt.close(figure)
    snapshot = output_dir / "source" / Path(__file__).name
    snapshot.parent.mkdir(exist_ok=True); shutil.copyfile(__file__, snapshot)
    manifest = build_manifest(asset_id="meson_freezeout_defaults_20261017", figure_family="charged_meson_freezeout",
        case_slug=output_dir.name, figure_mode="audit", semantic_status="default_parameter_comparison",
        style_profile=profile.profile_id, publication_scope="internal_review",
        generator=generator_record(__file__, command=" ".join(sys.argv), source_snapshot=snapshot,
            runtime=runtime_record({"matplotlib": matplotlib.__version__})),
        inputs=[input_record(run_dir / "freezeout_ratios.csv", role="direct_q0_calculations", units={"T":"MeV","muB":"MeV","sqrt_s_NN":"GeV","density":"fm^-3","ratio":"dimensionless"}),
            input_record(run_dir / "run.json", role="calculation_identity"),
            input_record(run_dir / "manifest.json", role="calculation_manifest"),
            input_record(experiment_path, role="experimental_observations"),
            input_record(experiment_path.parent / "sources.json", role="experimental_provenance"),
            input_record(ROOT / "config/plotting/candidate_aps_v2.toml", role="plot_style_contract"),
            *[input_record(ROOT / "scripts/plotting" / name, role="shared_plotting_code") for name in
                ("plot_style.py", "plot_manifest.py", "plot_quality.py", "plot_accessibility.py")]],
        axes=[{"field":"sqrt_s_NN_GeV","source_unit":"GeV","display_unit":"GeV","scale":"log","label":xlabel},
            {"field":"Kplus_over_pi_plus","source_unit":"dimensionless","display_unit":"dimensionless","scale":"linear","label":ylabel}],
        series=series, outputs=outputs, selection_rule="all directly calculated energies and the retained experimental table",
        interpolation_policy="none", connector_policy="forbidden",
        missing_value_policy="NaN breaks; no zero-fill, smoothing, or gap bridge",
        validation={"finite":True,"duplicate_keys":True,"support":True,
            "solver_called":False,"canonical_data_modified":False,"checks":"hashes, ratio finiteness, experiment energy coverage, measured rendering","author_review":"pending"},
        rendering={"delivery_stage":"png_review","vector_delivery_pending":True,
            "figure_size_inches":list(figure_size_inches(profile)),"column":"double_column","font":font,
            "legend_policy":"shared_in_panel_reviewed_geometry_checked",
            "legend_layout_reason":"author preference for legends inside the same panel; geometry and curves checked",
            "color_route":"undecided_review","quality":quality,
            "placement_limits":placement_limits(quality,profile,outputs),
            "accessibility":accessibility_record(encodings,backgrounds=["#FFFFFF"],background_scope="white axes; no fill"),
            "connections":"straight segments between adjacent valid direct support only",
            "HADES_error_display":"marker only; pixel uncertainty is not an experimental error bar"},
        calculation_sha=state["git_head"], source_run_id=calculation["scan_identity"])
    manifest.update({"delivery_stage":"png_review", "manuscript_eligible":False,
        "current_publication_layer":False,"vector_delivery_pending":True,
        "observable_note":"Local quark-only meson densities compared with final midrapidity hadron yields; no decay or acceptance model. ALICE is a charge sum.",
        "analysis_outputs":[input_record(output_dir/name,role="derived_comparison") for name in ("comparison_metrics.json","comparison_residuals.csv")]})
    write_manifest(output_dir / "plot_manifest.json", manifest)
    issues = validate_manifest(output_dir / "plot_manifest.json", repo_root=ROOT)
    if issues:
        raise ValueError("plot contract failed: " + "; ".join(issues))
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--experiment", type=Path, default=ROOT / "data/inputs/meson_freezeout_experiment_20261017/experimental_ratios.csv")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite plot case: {output}")
    models, calculation, state = load_calculation(args.run_dir.resolve())
    experiments = load_experiment(args.experiment.resolve())
    metrics, residuals = compare(models, experiments)
    output.mkdir(parents=True)
    (output / "comparison_metrics.json").write_text(json.dumps(metrics,indent=2)+"\n",encoding="utf-8")
    with (output / "comparison_residuals.csv").open("w",encoding="utf-8",newline="") as handle:
        writer = csv.DictWriter(handle,fieldnames=list(residuals[0])); writer.writeheader(); writer.writerows(residuals)
    render(models,experiments,output,calculation,state,args.run_dir.resolve(),args.experiment.resolve())
    print(json.dumps({"figure":str(output / "kplus_piplus_freezeout_defaults.png"),
        "best_default":metrics["best_default_on_primary_metric"],"comparison_points":metrics["common_primary_points"]}))


if __name__ == "__main__":
    main()
