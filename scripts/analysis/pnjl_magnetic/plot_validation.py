"""PNG review figures for the isolated magnetic single-anchor PALC study.

No equation is solved here. Native branch order is preserved, including folds;
paper points and the visibly different MFIR/smooth regulators stay distinct.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scripts.plotting.plot_style import load_profile, style_context, configure_axis_ticks
from scripts.plotting.plot_series import SeriesStyle, SeriesStyleRegistry
from scripts.plotting.plot_manifest import input_record, generator_record, build_manifest, runtime_record, write_manifest
from scripts.plotting.plot_quality import export_figure
from scripts.plotting.plot_bundle import build_bundle
from scripts.plotting.plot_accessibility import accessibility_record


def load_cases(path: Path) -> dict[float, dict]:
    data = json.loads(path.read_text(encoding="utf-8"))
    cases = {}
    for case in data["cases"]:
        if "error" in case:
            raise ValueError(f"failed input case: {case}")
        key = case["eB_GeV2"]
        if key in cases:
            raise ValueError("duplicate magnetic case")
        cases[key] = case
    return cases


def values(case: dict, field: str) -> np.ndarray:
    return np.asarray([p[field] for p in case["points"]], dtype=float)


def draw_branch(ax, case, style, label, field="rhoB_over_rho0"):
    x, y = values(case, "muB_MeV"), values(case, field)
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("nonfinite native branch point")
    ax.plot(x, y, label=label, **style.plot_kwargs())


def make_axes(ctx, nrows=1, ncols=3, size=(6.75, 3.25)):
    fig, axes = plt.subplots(nrows, ncols, figsize=size, squeeze=False)
    ctx.bind(fig)
    fig.subplots_adjust(left=0.115, right=0.97, bottom=0.21 if nrows == 1 else 0.11,
                        top=0.75 if nrows == 1 else 0.92, wspace=0.28, hspace=0.12)
    return fig, axes


def legend_above(fig, ax):
    handles, labels = ax.get_legend_handles_labels()
    return fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.54, 0.99),
                      ncol=len(labels), handlelength=2.0, columnspacing=1.2)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results_dir", type=Path)
    parser.add_argument("figures_dir", type=Path)
    args = parser.parse_args()
    results_dir, figures_dir = args.results_dir.resolve(), args.figures_dir.resolve()
    figures_dir.mkdir(parents=True, exist_ok=True)
    # The explicitly coarse production-kernel probe is p32; p128 is kept as a
    # separate later run and is not silently substituted here.
    current_path = results_dir / "current_mfir_p32_high.json"
    mfir_path = results_dir / "resolved_mfir_o48_fine_v2.json"
    smooth_path = results_dir / "resolved_smooth_o48_fine_v2.json"
    points_path = results_dir / "thesis_fig5_3_points.csv"
    current, mfir, smooth = (load_cases(p) for p in (current_path, mfir_path, smooth_path))
    with points_path.open(encoding="utf-8", newline="") as f:
        paper_points = list(csv.DictReader(f))
    profile = load_profile("candidate_aps_v2")
    registry = SeriesStyleRegistry("magnetic_palc_validation_v1", (
        SeriesStyle("current_p32", "#AD5E8D", ":"),
        SeriesStyle("mfir_reference", "#0072B2", "-"),
        SeriesStyle("smooth_reference", "#D55E00", "--"),
        SeriesStyle("thesis_points", "#000000", "None", "o"),
        SeriesStyle("u_density", "#0072B2", "-."),
        SeriesStyle("d_density", "#D55E00", ":"),
    ))
    styles = {s.series_id: s for s in registry.styles}
    source_paths = (current_path, mfir_path, smooth_path, points_path,
                    ROOT / "config/plotting/candidate_aps_v2.toml")
    inputs = [input_record(p, role="native_validation_results_or_paper_readout", root=ROOT) for p in source_paths]
    records, charts = [], []

    def finish(fig, stem, selected, axes, note):
        outputs, quality = export_figure(fig, figures_dir / stem, profile, formats=["png"])
        plt.close(fig)
        access = accessibility_record([
            {"value": s.series_id, "color": s.color, "noncolor_encoding": f"{s.linestyle}/{s.marker}/{s.series_id}"}
            for s in registry.select(selected)], backgrounds=["#ffffff"],
            background_scope="white axes; no fills or overlapping shaded bands")
        manifest = build_manifest(asset_id=stem, figure_family="magnetic_palc_validation",
            case_slug=figures_dir.name, figure_mode="audit", semantic_status="isolated_low_temperature_validation",
            style_profile=profile.profile_id, publication_scope="internal_review", root=ROOT,
            generator=generator_record(Path(__file__), root=ROOT,
                command=subprocess.list2cmdline([sys.executable, *sys.argv]),
                runtime=runtime_record({"matplotlib": matplotlib.__version__})), inputs=inputs,
            axes=axes, series=[{"series_id": s.series_id, "state": "diagnostic",
                "support_rule": "saved native branch vertices and separately calibrated paper features",
                "mask_rule": "no synthetic point, smoothing, root rejection or gap fill"} for s in registry.select(selected)],
            outputs=outputs, selection_rule=note,
            interpolation_policy="none; connected native continuation vertices in arclength order",
            connector_policy="native adjacent vertices only", missing_value_policy="error; no fill",
            validation={"finite": True, "duplicate_keys": True, "support": True, "strict_gate": False},
            rendering={"column": "double_column", "figure_size_inches": quality["figure_size_inches"],
                "size_override_reason": "three magnetic panels with one shared key; heights vary with quantities",
                "color_route": "undecided_review", "legend_policy": "shared_external_then_case_review",
                "quality": quality, "delivery_stage": "png_review", "vector_delivery_pending": True,
                "review_contract": "prd_review_v1", "accessibility": access,
                "series_styles": registry.record(selected),
                "typography_exception": "dense_composite_review_compact_typography"},
        )
        # This isolated review bundle is generated from the current dirty
        # worktree; validate the exact working-tree hashes rather than asking
        # the historical commit to contain this untracked analysis script.
        manifest["git_commit"] = None
        manifest.update({"manuscript_eligible": False, "current_publication_layer": False,
                         "production_default_changed": False, "scientific_note": note})
        records.append(manifest)
        charts.append({"asset_id": stem, "png": outputs[0]["path"]})

    rho_axes = [{"field": "muB_MeV", "source_unit": "MeV", "display_unit": "MeV",
                 "label": r"$\mu_B\;(\mathrm{MeV})$"},
                {"field": "rhoB_over_rho0", "source_unit": "dimensionless", "display_unit": "dimensionless",
                 "label": r"$\rho_B/\rho_0$"}]
    with style_context(profile) as ctx:
        fig, aa = make_axes(ctx)
        for j, eb in enumerate((0.05, 0.155, 0.4)):
            ax = aa[0,j]
            draw_branch(ax, current[eb], styles["current_p32"], "Current p32")
            draw_branch(ax, mfir[eb], styles["mfir_reference"], "Resolved MFIR")
            draw_branch(ax, smooth[eb], styles["smooth_reference"], "Thesis regulator")
            ax.set_xlim((950,1550) if eb == 0.05 else ((820,1560) if eb == 0.155 else (550,1250)))
            ax.set_ylim(0,8 if eb != .4 else 5)
            ax.set_xlabel(r"$\mu_B\;(\mathrm{MeV})$")
            ax.set_title(rf"$eB={eb:g}\;\mathrm{{GeV}}^2$", pad=7)
            configure_axis_ticks(ax, profile)
        aa[0,0].set_ylabel(r"$\rho_B/\rho_0$")
        legend_above(fig,aa[0,0])
        finish(fig,"rho_mu_three_fields",["current_p32","mfir_reference","smooth_reference"],rho_axes,
            "T=5 MeV; a=0.0108805 for all curves. Current diagnostic uses the unsplit p32 kernel; references use independent Fermi-surface-resolved quadrature. Regulators are different.")

    with style_context(profile) as ctx:
        fig, axx = plt.subplots(figsize=(6.75,4.0));ctx.bind(fig)
        fig.subplots_adjust(left=.115,right=.97,bottom=.17,top=.78)
        draw_branch(axx,smooth[.05],styles["smooth_reference"],"Thesis regulator")
        xx=np.array([float(r['x']) for r in paper_points]); yy=np.array([float(r['y']) for r in paper_points])
        kw=styles["thesis_points"].plot_kwargs()
        axx.errorbar(xx,yy,xerr=np.array([[float(r['x'])-float(r['x_lower']) for r in paper_points],
                      [float(r['x_upper'])-float(r['x']) for r in paper_points]]),
                     yerr=np.array([[float(r['y'])-float(r['y_lower']) for r in paper_points],
                      [float(r['y_upper'])-float(r['y']) for r in paper_points]]),
                     label="Thesis Fig. 5-3",markersize=4.5,capsize=2,
                     color=kw["color"],linestyle="None",marker=kw["marker"])
        # Error-bar helper artists are deliberately unregistered. The single
        # marker line below is the scientific series identity measured by the
        # plotting quality contract.
        axx.plot(xx, yy, label="_nolegend_", color=kw["color"],
                 linestyle=kw["linestyle"], marker=kw["marker"], gid="thesis_points")
        axx.set_xlim(1030,1220);axx.set_ylim(0,2.5)
        axx.set_xlabel(r"$\mu_B\;(\mathrm{MeV})$");axx.set_ylabel(r"$\rho_B/\rho_0$")
        configure_axis_ticks(axx,profile);legend_above(fig,axx)
        finish(fig,"weak_field_thesis_features",["smooth_reference","thesis_points"],rho_axes,
            "T=5 MeV, eB=0.05 GeV2. Five preselected paper fold features; readout ranges are pixel-selection and calibration ranges, not statistical errors. Native curve is unshifted.")

    with style_context(profile) as ctx:
        fig, aa=make_axes(ctx,nrows=2,size=(6.75,6.0))
        fig.subplots_adjust(top=.82, hspace=.28, bottom=.15)
        for j,eb in enumerate((.05,.155,.4)):
            for family,sty in ((mfir,"mfir_reference"),(smooth,"smooth_reference")):
                pts=family[eb]["points"]
                mu=values(family[eb],"muB_MeV")
                aa[0,j].plot(mu,[p['masses_MeV'][0] for p in pts],label="MFIR" if family is mfir else "Thesis regulator",**styles[sty].plot_kwargs())
                aa[1,j].plot(mu,[p['masses_MeV'][1] for p in pts],**styles[sty].plot_kwargs())
            for i in (0,1):
                aa[i,j].set_xlim(500,1400);aa[i,j].set_ylim(0,470);aa[i,j].set_yticks([0,200,400]);configure_axis_ticks(aa[i,j],profile)
                if i == 1:
                    aa[i,j].set_xticks([800, 1200])
            aa[0,j].tick_params(labelbottom=False)
            aa[0,j].set_title(rf"$eB={eb:g}\;\mathrm{{GeV}}^2$",pad=7)
            aa[1,j].set_xlabel(r"$\mu_B\;(\mathrm{MeV})$")
        aa[0,0].set_ylabel(r"$M_u\;(\mathrm{MeV})$");aa[1,0].set_ylabel(r"$M_d\;(\mathrm{MeV})$")
        legend_above(fig, aa[0,0])
        finish(fig,"constituent_masses",["mfir_reference","smooth_reference"],
            [rho_axes[0],{"field":"masses_MeV","source_unit":"MeV","display_unit":"MeV","label":r"$M_{u,d}\;(\mathrm{MeV})$"}],
            "T=5 MeV, flavors in separate rows. Strong-field smooth vacuum masses reproduce the chapter-5 431/400 MeV quoted values; no regulator retuning.")

    with style_context(profile) as ctx:
        fig,aa=make_axes(ctx,ncols=2,size=(6.75,3.55))
        pts=smooth[.4]["points"]
        rho=values(smooth[.4],"rhoB_over_rho0")
        u=np.asarray([pt['rho'][0]/.16 for pt in pts]); d=np.asarray([pt['rho'][1]/.16 for pt in pts])
        aa[0,0].plot(rho,u,label=r"$\rho_u/\rho_0$",**styles["u_density"].plot_kwargs())
        aa[0,0].plot(rho,d,label=r"$\rho_d/\rho_0$",**styles["d_density"].plot_kwargs())
        aa[0,0].set_xlim(0,5);aa[0,0].set_ylim(0,12)
        aa[0,0].set_xlabel(r"$\rho_B/\rho_0$");aa[0,0].set_ylabel(r"$\rho_f/\rho_0$")
        # This zoom includes the anti-S turning; no thermodynamic transition
        # line is inferred from its pressure crossing.
        draw_branch(aa[0,1],smooth[.4],styles["smooth_reference"],"Anti-S")
        aa[0,1].set_xlim(1020,1100);aa[0,1].set_ylim(.3,1.8)
        aa[0,1].set_xlabel(r"$\mu_B\;(\mathrm{MeV})$");aa[0,1].set_ylabel(r"$\rho_B/\rho_0$")
        for ax in aa.ravel():configure_axis_ticks(ax,profile)
        legend_above(fig,aa[0,0])
        finish(fig,"strong_field_flavor_and_antiS",["u_density","d_density","smooth_reference"],
            [rho_axes[0],rho_axes[1],{"field":"rho_f/rho0","source_unit":"dimensionless","display_unit":"dimensionless","label":r"$\rho_f/\rho_0$"}],
            "eB=0.4 GeV2, T=5 MeV, thesis regulator; density flavors split as in Fig5-4 and anti-S as in Fig5-5. No phase/Maxwell calculation.")

    bundle=build_bundle({"case_slug":figures_dir.name,"charts":charts},records)
    write_manifest(figures_dir/"plot_manifest.json",bundle)
    print(json.dumps({"figure_paths":[c['png'] for c in charts]},indent=2))


if __name__ == "__main__":
    main()

