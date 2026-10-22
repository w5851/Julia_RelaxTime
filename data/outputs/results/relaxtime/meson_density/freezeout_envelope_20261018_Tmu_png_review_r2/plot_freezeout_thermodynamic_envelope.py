"""Prepare analytic T(muB) envelopes from frozen parameters, then render frozen tables."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import shutil
import sys
import tomllib

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from scripts.analysis.relaxtime.freezeout_parameter_envelope import load_designs, trajectory, write_csv, write_json
from scripts.plotting.plot_manifest import (build_manifest, generator_record, input_record,
    output_record, runtime_record, sha256_file, write_manifest)
from scripts.plotting.plot_quality import export_figure, inspect_export, placement_limits
from scripts.plotting.plot_style import load_profile, configure_matplotlib, configure_axis_ticks, figure_size_inches
from scripts.plotting.plot_accessibility import accessibility_record
from scripts.plotting.validate_plot_artifact import validate_manifest


def table(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def checked_json(path, expected=None):
    if expected is not None and sha256_file(path) != expected:
        raise ValueError(f"frozen input changed: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def polynomial_temperature(parameters, muB_GeV):
    if len(parameters) != 5 or not all(map(math.isfinite, parameters)):
        raise ValueError("five finite polynomial parameters required")
    if not math.isfinite(muB_GeV) or muB_GeV < 0:
        raise ValueError("muB_GeV must be nonnegative and finite")
    a, b, c, _, _ = parameters
    return a - b*muB_GeV**2 - c*muB_GeV**4


def andronic_temperature(coefficients, muB_GeV):
    t0, offset, width, d, e = [coefficients[k] for k in
        ("T_lim_GeV", "log_offset", "log_width", "d_GeV", "e_GeV_inv1")]
    if not all(map(math.isfinite, (t0, offset, width, d, e))) or min(t0, width, d, e) <= 0:
        raise ValueError("invalid Andronic coefficients")
    if not math.isfinite(muB_GeV) or not 0 <= muB_GeV < d:
        raise ValueError("Andronic muB_GeV must satisfy 0 <= muB < d")
    if muB_GeV == 0:
        return t0
    # Invert muB=d/(1+e*sqrt(s)); logarithms avoid overflow near muB=0.
    log_energy = math.log(d-muB_GeV)-math.log(muB_GeV)-math.log(e)
    z = offset-log_energy/width
    if z >= 0:
        exp_minus = math.exp(-z)
        return t0*exp_minus/(1+exp_minus)
    return t0/(1+math.exp(z))


def fixed_temperature_band(candidates, members, mu_grid_GeV):
    if not members or len(set(members)) != len(members):
        raise ValueError("band requires a fixed, nonempty, unique member set")
    if not mu_grid_GeV or any(not math.isfinite(x) or x < 0 for x in mu_grid_GeV):
        raise ValueError("mu grid must be finite and nonnegative")
    if any(b <= a for a, b in zip(mu_grid_GeV, mu_grid_GeV[1:])):
        raise ValueError("mu grid must be strictly increasing")
    rows = []
    for mu in mu_grid_GeV:
        values = [(polynomial_temperature(candidates[n]["values"], mu), n) for n in members]
        if any(not math.isfinite(t) or t <= 0 for t, _ in values):
            raise ValueError("display domain contains nonpositive T; narrow it without dropping band members")
        lo, hi = min(values), max(values)
        rows.append(dict(muB_MeV=1000*mu, minimum_T_MeV=1000*lo[0], maximum_T_MeV=1000*hi[0],
            minimum_profile=lo[1], maximum_profile=hi[1], member_count=len(members)))
    return rows


def prepare(args):
    source, output = args.summary_dir.resolve(), args.output.resolve()
    if output.exists():
        raise ValueError("prepared case exists; choose a new sibling directory")
    if not math.isfinite(args.mu_max_MeV) or args.mu_max_MeV <= 0 or args.points < 2:
        raise ValueError("positive finite mu-max-MeV and at least two points required")
    for name, digest in checked_json(source / "checksums.json").items():
        if sha256_file(source / name) != digest:
            raise ValueError(f"summary input changed: {name}")
    summary = checked_json(source / "summary.json")
    screen_path = Path(summary["screen_assessment"]["path"])
    screen = checked_json(screen_path, summary["screen_assessment"]["sha256"])
    for record in screen["design_inputs"]:
        checked_json(Path(record["path"]), record["sha256"])
    designs, candidates = load_designs([Path(r["path"]) for r in screen["design_inputs"]])
    members, near = screen["band_complete_profiles"], screen["near_best_band_complete_profiles"]
    best = summary["best_selected"]
    if best not in near or not set(near).issubset(members):
        raise ValueError("inconsistent frozen envelope membership")
    andronic_path = None
    for record in summary["baseline_inputs"]:
        run_dir = Path(record["directory"])
        manifest = checked_json(run_dir / "manifest.json", record["manifest_sha256"])
        if "andronic2018" in manifest["profiles"]:
            state = checked_json(run_dir / "run.json", manifest["run_sha256"])
            relative = "config/physics/freezeout/andronic2018.toml"
            andronic_path = run_dir / "source_snapshot" / relative
            if sha256_file(andronic_path) != state["source_hashes"][relative]:
                raise ValueError("frozen Andronic profile changed")
            break
    if andronic_path is None:
        raise ValueError("frozen Andronic default not found")
    andronic_profile = tomllib.loads(andronic_path.read_text(encoding="utf-8"))
    if andronic_profile["family"] != "andronic_logistic_rational_v1":
        raise ValueError("unexpected Andronic temperature family")
    andronic = andronic_profile["coefficients"]
    energy_min, energy_max = [designs[0]["design"][k] for k in ("minimum_energy_GeV", "maximum_energy_GeV")]
    anchors = []
    for name in (best, "andronic2018"):
        for energy in (energy_min, 7.7, 19.6, energy_max):
            if name == best:
                t, mu = trajectory(candidates[name]["values"], energy)
            else:
                mu = andronic["d_GeV"]/(1+andronic["e_GeV_inv1"]*energy)
                t = andronic_temperature(andronic, mu)
            anchors.append(dict(profile=name, energy_GeV=energy, muB_MeV=mu*1000, T_MeV=t*1000,
                inside_display_domain=mu*1000 <= args.mu_max_MeV))
    if any(a["energy_GeV"] == 7.7 and not a["inside_display_domain"] for a in anchors):
        raise ValueError("mu-max-MeV must include both 7.7 GeV marker positions")
    grid = sorted(set([args.mu_max_MeV/1000*i/(args.points-1) for i in range(args.points)] +
        [a["muB_MeV"]/1000 for a in anchors if a["inside_display_domain"]]))
    band = fixed_temperature_band(candidates, members, grid)
    near_band = fixed_temperature_band(candidates, near, grid)
    for row, near_row, mu in zip(band, near_band, grid):
        row.update(near_minimum_T_MeV=near_row["minimum_T_MeV"], near_maximum_T_MeV=near_row["maximum_T_MeV"],
            near_minimum_profile=near_row["minimum_profile"], near_maximum_profile=near_row["maximum_profile"],
            near_member_count=len(near), best_T_MeV=1000*polynomial_temperature(candidates[best]["values"], mu),
            andronic_T_MeV=1000*andronic_temperature(andronic, mu))
    errors = []
    coordinate_count = 0
    for row in table(source / "fine_curves.csv") + table(source / "baseline_curves.csv"):
        name, energy = row["profile"], float(row["energy_GeV"])
        if name in candidates and name in summary["selected_fine_ranking"]:
            t, mu = trajectory(candidates[name]["values"], energy)
            t_at_mu = polynomial_temperature(candidates[name]["values"], float(row["muB_MeV"])/1000)
        elif name == "andronic2018":
            mu = andronic["d_GeV"]/(1+andronic["e_GeV_inv1"]*energy)
            t = andronic_temperature(andronic, mu)
            t_at_mu = andronic_temperature(andronic, float(row["muB_MeV"])/1000)
        else:
            continue
        errors += [abs(t*1000-float(row["T_MeV"])), abs(mu*1000-float(row["muB_MeV"])),
            abs(t_at_mu*1000-float(row["T_MeV"]))]
        coordinate_count += 1
    if not errors or max(errors) > 1e-8:
        raise ValueError("analytic coordinates disagree with frozen calculation coordinates")
    metadata = dict(schema="chemical_freezeout_Tmu_envelope_v1", best_selected=best,
        band_members=members, near_best_members=near,
        near_best_distinct_Tmu_curves=len({tuple(candidates[n]["values"][:3]) for n in near}),
        selected_parameters=summary["selected_parameters"], andronic_default=andronic,
        optimized_family="cleymans_polynomial_rational_v1", andronic_optimized=False,
        muB_range_MeV=[0, args.mu_max_MeV], original_energy_range_GeV=[energy_min, energy_max],
        base_grid_points=args.points, grid_points=len(grid),
        definition="At fixed muB take T minima/maxima over the same frozen full and near-best sets as the ratio figure; no new selection",
        coordinate_validation=dict(points=coordinate_count, maximum_error_MeV=max(errors)),
        support_semantics="Analytic parameterizations on the declared muB domain, including continuations outside each member's original finite energy window; no new density solutions",
        projection_semantics="Only a,b,c determine polynomial T(muB); d,e set the energy positions on that curve",
        failure_semantics="Earlier ratio failures remain unchanged; analytic chemical-freezeout geometry itself has no K+ integration gate",
        calculation_sha=summary["calculation_sha"], source_run_id=summary["fine_inputs"][0]["scan_identity"],
        inputs=[input_record(source / n,role="frozen_parameter_selection") for n in
            ("summary.json", "checksums.json", "fine_curves.csv", "baseline_curves.csv")] +
            [input_record(screen_path,role="fixed_band_membership"), input_record(andronic_path,role="frozen_Andronic_default")] +
            [input_record(Path(r["path"]),role="frozen_sampled_parameters") for r in screen["design_inputs"]])
    output.mkdir(parents=True)
    write_csv(output / "thermodynamic_band.csv", band)
    write_csv(output / "energy_anchors.csv", anchors)
    write_csv(output / "parameters.csv", [dict(profile=n, **dict(zip(("a_GeV","b_GeV_inv1","c_GeV_inv3","d_GeV","e_GeV_inv1"), candidates[n]["values"]))) for n in members])
    write_json(output / "preparation.json", metadata)
    for path in (Path(__file__), ROOT / "scripts/analysis/relaxtime/freezeout_parameter_envelope.py"):
        shutil.copy2(path, output / path.name)
    write_json(output / "checksums.json", {p.name:sha256_file(p) for p in sorted(output.iterdir()) if p.is_file()})
    print(json.dumps(dict(output=str(output), members=len(members), near_members=len(near),
        distinct_near_curves=metadata["near_best_distinct_Tmu_curves"], coordinate_validation=metadata["coordinate_validation"])))


def render(args):
    source, output = args.prepared_dir.resolve(), args.output_dir.resolve()
    results = ROOT / "data/outputs/results/relaxtime/meson_density" / output.name
    if output.exists() or results.exists():
        raise ValueError("figure case exists; choose a new sibling case")
    for name, digest in checked_json(source / "checksums.json").items():
        if sha256_file(source / name) != digest:
            raise ValueError(f"prepared input changed: {name}")
    meta = checked_json(source / "preparation.json")
    bands, anchors = table(source / "thermodynamic_band.csv"), table(source / "energy_anchors.csv")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator
    from PIL import Image
    profile = load_profile("candidate_aps_v2")
    font = configure_matplotlib(profile)
    figure, ax = plt.subplots(figsize=figure_size_inches(profile))
    figure.subplots_adjust(left=.115,right=.975,bottom=.16,top=.96)
    x = [float(r["muB_MeV"]) for r in bands]
    ax.fill_between(x,[float(r["minimum_T_MeV"]) for r in bands],[float(r["maximum_T_MeV"]) for r in bands],
        facecolor="#E8E8E8",edgecolor="#595959",linewidth=.65,linestyle="--",label=f"Sampled band ({len(meta['band_members'])} sets)",zorder=1)
    ax.fill_between(x,[float(r["near_minimum_T_MeV"]) for r in bands],[float(r["near_maximum_T_MeV"]) for r in bands],
        facecolor="#C6DEF0",edgecolor="#0072B2",linewidth=.65,hatch="/",label=f"Near-best band ({len(meta['near_best_members'])} sets)",zorder=2)
    for name,key,color,style,marker,label in (
        (meta["best_selected"],"best_T_MeV","#0072B2","-","o","Best polynomial"),
        ("andronic2018","andronic_T_MeV","#A04800","--","s","Andronic 2018 (default)")):
        anchor = next(r for r in anchors if r["profile"] == name and float(r["energy_GeV"]) == 7.7)
        index = min(range(len(x)),key=lambda i:abs(x[i]-float(anchor["muB_MeV"])))
        if abs(x[index]-float(anchor["muB_MeV"])) > 1e-9:
            raise ValueError("7.7 GeV anchor outside frozen display grid")
        ax.plot(x,[float(r[key]) for r in bands],color=color,linestyle=style,linewidth=1.7,
            marker=marker,markevery=[index],markersize=6,markerfacecolor="white",markeredgewidth=1.1,
            label=label,zorder=3)
    ax.set_xlim(0,meta["muB_range_MeV"][1]+10)
    lowest = min(float(r["minimum_T_MeV"]) for r in bands)
    highest = max(float(r["maximum_T_MeV"]) for r in bands)
    ax.set_ylim(min(50,math.floor((lowest-5)/20)*20),max(180,math.ceil((highest+10)/20)*20))
    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.yaxis.set_major_locator(MultipleLocator(20))
    configure_axis_ticks(ax,profile)
    xlabel,ylabel = r"$\mu_B\;(\mathrm{MeV})$",r"$T\;(\mathrm{MeV})$"
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel,labelpad=8)
    ax.text(.03,.975,"Chemical freezeout (analytic)",transform=ax.transAxes,va="top",fontsize=12)
    legend_title = r"Markers: $\sqrt{s_{NN}}=7.7\;\mathrm{GeV}$"
    ax.legend(loc="lower left",title=legend_title,title_fontsize=13,fontsize=12,handlelength=2.7,
        labelspacing=.4,borderaxespad=.6,frameon=False)
    labels = ax.get_legend_handles_labels()[1]
    output.mkdir(parents=True)
    results.mkdir(parents=True)
    outputs,quality = export_figure(figure,output / "chemical_freezeout_T_muB_band",profile,formats=["png"])
    gray = output / "chemical_freezeout_T_muB_band_gray.png"
    with Image.open(output / "chemical_freezeout_T_muB_band.png") as im:
        im.convert("L").convert("RGB").save(gray,dpi=(profile.dpi,profile.dpi))
    gray_record = output_record(gray,fmt="png",dpi=profile.dpi,vector=False)
    gray_record["inspection"] = inspect_export(gray)
    outputs.append(gray_record)
    plt.close(figure)
    snapshot = results / Path(__file__).name
    shutil.copy2(__file__,snapshot)
    caption = ("Chemical-freezeout T(muB) parameterizations. Gray and hatched blue use exactly the same "
        f"{len(meta['band_members'])} full and {len(meta['near_best_members'])} near-best parameter sets as the frozen K+/pi+ comparison. "
        "At each displayed muB, bounds are the minimum/maximum analytic temperature of those fixed sets. "
        f"The {len(meta['near_best_members'])} near-best parameter sets contain {meta['near_best_distinct_Tmu_curves']} distinct T(muB) polynomials: "
        "d,e only change the mapping to collision energy. "
        "Solid blue is the best sampled polynomial; orange dashed is the unoptimized Andronic 2018 default, obtained by analytically eliminating sqrt(s). "
        "Symbols locate sqrt(sNN)=7.7 GeV on each curve; they are parameterized coordinates, not experimental freezeout measurements. "
        f"The displayed 0-{meta['muB_range_MeV'][1]:g} MeV muB domain includes analytic continuations outside each curve's original "
        f"{meta['original_energy_range_GeV'][0]:g}-{meta['original_energy_range_GeV'][1]:g} GeV ratio-scan interval. "
        "These continuations and muB=0 limits are not additional q0 density results. No crossover curve is shown. "
        "The old Andronic K+ gate failure does not invalidate the analytic temperature parameterization; the original ratio gap is unchanged. "
        "Bands describe finite sampled parameters, not confidence intervals. PNG review only.\n")
    (results / "caption.md").write_text(caption,encoding="utf-8")
    encodings = [dict(value="sampled band",color="#595959",noncolor_encoding="dashed boundary"),
        dict(value="near-best band",color="#0072B2",noncolor_encoding="hatched fill"),
        dict(value="best polynomial",color="#0072B2",noncolor_encoding="solid and circle"),
        dict(value="Andronic default",color="#A04800",noncolor_encoding="dashed and square")]
    manifest = build_manifest(asset_id="meson_freezeout_Tmu_envelope_20261018",figure_family="charged_meson_freezeout",
        case_slug=output.name,figure_mode="audit",semantic_status="analytic_chemical_freezeout_parameter_envelope",
        style_profile=profile.profile_id,publication_scope="internal_review",
        generator=generator_record(__file__,command=" ".join(sys.argv),source_snapshot=snapshot,
            runtime=runtime_record({"matplotlib":matplotlib.__version__})),
        inputs=[input_record(source / n,role="frozen_analytic_parameter_table") for n in
            ("thermodynamic_band.csv","energy_anchors.csv","parameters.csv","preparation.json","checksums.json")] +
            meta["inputs"] + [input_record(ROOT / "config/plotting/candidate_aps_v2.toml",role="plot_style_contract")] +
            [input_record(ROOT / "scripts/plotting" / n,role="shared_plotting_code") for n in
                ("plot_style.py","plot_manifest.py","plot_quality.py","plot_accessibility.py","validate_plot_artifact.py")],
        axes=[dict(field="muB_MeV",source_unit="MeV",display_unit="MeV",scale="linear",label=xlabel),
            dict(field="T_MeV",source_unit="MeV",display_unit="MeV",scale="linear",label=ylabel)],
        series=[dict(series_id=name,state="analytic_parameterization",point_count=len(bands),
            support_rule="frozen direct analytic evaluations, including explicitly declared continuation",
            mask_rule="all fixed members have positive finite T throughout the display domain; no member removal",
            members=members) for name,members in (("sampled_band",meta["band_members"]),
                ("near_best_band",meta["near_best_members"]),(meta["best_selected"],[meta["best_selected"]]),("andronic2018",["andronic2018"]))],
        outputs=outputs,selection_rule=f"fixed membership inherited from ratio study; analytic muB display domain 0-{meta['muB_range_MeV'][1]:g} MeV",
        interpolation_policy="none",connector_policy="forbidden",missing_value_policy="no density values inferred; retain original failures externally",
        validation=dict(finite=True,duplicate_keys=True,support=True,solver_called=False,canonical_data_modified=False,
            checks="frozen hashes, analytical identities and coordinate comparison, positive T, rendering geometry",author_review="pending"),
        rendering=dict(delivery_stage="png_review",vector_delivery_pending=True,review_contract="prd_review_v1",
            figure_size_inches=list(figure_size_inches(profile)),column="double_column",font=font,
            legend_policy="declared_geometry_checked",legend_outside=False,
            legend_placements=[dict(legend_index=0,host_axes_index=0,title=legend_title,labels=labels,
                placement="in_axes",location="lower left",scope="parameter bands, curve families and collision-energy positions")],
            legend_layout_reason="author requests in-panel legends; lower-left region is empty of curves",
            color_route="undecided_review",quality=quality,placement_limits=placement_limits(quality,profile,outputs),
            accessibility=accessibility_record(encodings,backgrounds=["#FFFFFF","#E8E8E8","#C6DEF0"],
                background_scope="white axes and light band fills; dark boundaries and hatching carry band semantics"),
            connections="straight display segments between frozen direct analytic evaluations; no interpolation of ratio data"),
        calculation_sha=meta["calculation_sha"],source_run_id=meta["source_run_id"])
    manifest.update(delivery_stage="png_review",manuscript_eligible=False,current_publication_layer=False,
        vector_delivery_pending=True,caption=input_record(results / "caption.md",role="figure_caption"))
    write_manifest(output / "plot_manifest.json",manifest)
    issues = validate_manifest(output / "plot_manifest.json",repo_root=ROOT)
    if issues:
        raise ValueError("plot contract failed: " + "; ".join(issues))
    print(json.dumps(dict(figure=str(output / "chemical_freezeout_T_muB_band.png"),validation="passed")))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command",required=True)
    prep = sub.add_parser("prepare")
    prep.add_argument("--summary-dir",type=Path,required=True)
    prep.add_argument("--output",type=Path,required=True)
    prep.add_argument("--mu-max-MeV",dest="mu_max_MeV",type=float,default=500.)
    prep.add_argument("--points",type=int,default=501)
    plot = sub.add_parser("render")
    plot.add_argument("--prepared-dir",type=Path,required=True)
    plot.add_argument("--output-dir",type=Path,required=True)
    args = parser.parse_args()
    (prepare if args.command == "prepare" else render)(args)


if __name__ == "__main__":
    main()
