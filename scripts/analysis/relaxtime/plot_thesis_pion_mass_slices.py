"""Render the thesis pion pole-mass slices at three fixed temperatures.

The thesis panel `fig:thesis-pion-slices` previously used a one-off reviewer
discussion figure drawn from the retired `mott_phase/reference_100_300_fine`
scan.  This generator rebuilds it from the promoted numerical layer
(`ward_moller_sphere_20261121_v1/spectra.csv`, 2010 frozen rows) on the
`candidate_aps_v2` profile, so the panel shares its prescription with the
neighbouring `pion_kaon_spectra` delivery and can be traced byte-for-byte.

Panel layout follows the review ruling: one row of three mass panels, dropping
the former Mott-temperature panel because its five points duplicate
`tab:mott-temperature` in the thesis text.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.ticker import FixedLocator

from scripts.plotting.plot_accessibility import accessibility_record
from scripts.plotting.plot_manifest import (build_manifest, generator_record, input_record,
                                            runtime_record, write_manifest)
from scripts.plotting.plot_quality import export_figure, placement_limits
from scripts.plotting.plot_series import SeriesStyle, SeriesStyleRegistry
from scripts.plotting.plot_style import configure_axis_ticks, load_profile
from scripts.plotting.publication_review import (REVIEW_COLORS, figure_context,
                                                 grayscale_review_record,
                                                 register_shared_series_styles,
                                                 write_sop_snapshot)
from scripts.plotting.review_vector_delivery import (COLOR_ROUTE, add_vector_arguments,
                                                    add_vector_pdf, read_accepted_review,
                                                    require_fresh_case)
from scripts.plotting.validate_plot_artifact import validate_manifest, validate_manifest_record

COLORS = REVIEW_COLORS
STEM = "pion_mass_slices"
TEMPERATURES = (180.0, 200.0, 220.0)
XI_POINTS = (-0.3, -0.15, 0.0, 0.15, 0.3)

REGIMES = SeriesStyleRegistry("thesis_pion_mass_regimes_v1", (
    SeriesStyle("pion_mass.slice_line", COLORS[0], "-", "None"),
    SeriesStyle("pion_mass.pole_point", COLORS[0], "None", "o"),
))
REGIME_LINES = {style.series_id: style for style in REGIMES.styles}


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def n(row: dict[str, str], field: str) -> float:
    return float(row[field])


def slice_rows(source: Path) -> dict[float, list[dict[str, float]]]:
    """Return the pi pole mass on the declared xi grid, grouped by temperature."""
    data = [row for row in rows(source) if row["meson"] == "pi"]
    selected: dict[float, list[dict[str, float]]] = {}
    for temperature in TEMPERATURES:
        picked = [row for row in data if abs(n(row, "T_MeV") - temperature) < 1e-9]
        if not picked:
            raise ValueError(f"source has no pion rows at T={temperature} MeV")
        ordered = sorted(picked, key=lambda row: n(row, "xi"))
        if [round(n(row, "xi"), 6) for row in ordered] != [round(x, 6) for x in XI_POINTS]:
            raise ValueError(f"xi grid at T={temperature} MeV is {[n(r, 'xi') for r in ordered]}")
        for row in ordered:
            if row["status"] not in ("bound", "above_threshold", "resonance", "scattering"):
                raise ValueError(f"unexpected pole status {row['status']!r} at T={temperature}")
        selected[temperature] = [
            dict(xi=n(row, "xi"), mass=n(row, "mass_MeV"), threshold=n(row, "threshold_MeV"))
            for row in ordered
        ]
    return selected


def render(slices: dict[float, list[dict[str, float]]], profile):
    with figure_context(profile, compact=True) as context:
        fig, axs = plt.subplots(1, 3, figsize=(6.75, 2.55), sharex=True)
        context.bind(fig)
        fig.subplots_adjust(left=0.115, right=0.985, bottom=0.34, top=0.90, wspace=0.30)
        limits = {180.0: (130.0, 140.0), 200.0: (130.0, 205.0), 220.0: (160.0, 340.0)}
        ticks = {180.0: (130, 135, 140), 200.0: (140, 170, 200), 220.0: (180, 240, 300)}
        for index, (temperature, ax) in enumerate(zip(TEMPERATURES, axs)):
            series = slices[temperature]
            xs = [point["xi"] for point in series]
            ax.plot(xs, [point["mass"] for point in series], lw=1.2, zorder=2,
                    **REGIME_LINES["pion_mass.slice_line"].plot_kwargs())
            for point in series:
                above = point["mass"] > point["threshold"]
                kw = REGIME_LINES["pion_mass.pole_point"].plot_kwargs()
                ax.plot([point["xi"]], [point["mass"]], ms=5.2, zorder=4,
                        markeredgecolor=COLORS[0], markeredgewidth=0.9,
                        markerfacecolor=COLORS[0] if above else "white", **kw)
            low, high = limits[temperature]
            ax.set_ylim(low, high)
            ax.yaxis.set_major_locator(FixedLocator(list(ticks[temperature])))
            ax.set_xlim(-0.34, 0.34)
            ax.set_xticks([-0.3, 0.0, 0.3], labels=[r"$-0.3$", r"$0$", r"$0.3$"])
            ax.xaxis.set_minor_locator(FixedLocator([-0.15, 0.15]))
            ax.text(0.96, 0.95, f"({chr(97 + index)}) " + rf"$T={int(temperature)}$ MeV",
                    transform=ax.transAxes, ha="right", va="top", fontsize=11, fontweight="normal")
            configure_axis_ticks(ax, profile)
        axs[0].set_ylabel(r"$M_\pi$ (MeV)", labelpad=6)
        fig.supxlabel(r"$\xi$", x=0.52, y=0.17, fontsize=13)
        handles = [
            Line2D([], [], color=COLORS[0], ls="none", marker="o", mfc="white", mec=COLORS[0],
                   mew=0.9, ms=5, label=r"$M_\pi<M_u+M_d$"),
            Line2D([], [], color=COLORS[0], ls="none", marker="o", mfc=COLORS[0], mec=COLORS[0],
                   mew=0.9, ms=5, label=r"$M_\pi>M_u+M_d$"),
        ]
        fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.52, 0.035), ncol=2,
                   frameon=False, fontsize=10.5, handlelength=0.9, handletextpad=0.4,
                   columnspacing=1.6, borderpad=0.1)
        return fig


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True,
                        help="ward_moller_sphere case holding spectra.csv")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--review-results-dir", type=Path,
                        help="Case results directory for the frozen SOP snapshot")
    add_vector_arguments(parser, multiple=False)
    args = parser.parse_args()

    out, meta = args.out.resolve(), args.results.resolve()
    results = args.review_results_dir.resolve() if args.review_results_dir else None
    require_fresh_case(out, *([results] if results else []))
    acceptance = None
    if args.deliver_vector:
        if results is None:
            parser.error("vector delivery requires --review-results-dir for the new SOP snapshot")
        acceptance = read_accepted_review(args.acceptance, root=ROOT, generator=Path(__file__))
        if acceptance.record["asset_id"] != STEM:
            parser.error(f"--acceptance is for {acceptance.record['asset_id']}, expected {STEM}")

    profile = load_profile("candidate_aps_v2")
    palette = {**profile.data["palette"], "colors": [*COLORS, *profile.colors[3:]]}
    profile = replace(profile, data={**profile.data, "palette": palette})
    out.mkdir(parents=True, exist_ok=True)

    source = meta / "spectra.csv"
    numerical = json.loads((meta / "numerical_manifest.json").read_text(encoding="utf-8"))
    slices = slice_rows(source)

    inputs = [(meta / "numerical_manifest.json", "frozen numerical provenance"),
              (source, "computed spectral series"),
              (profile.path, "plot profile")]
    inputs.extend((ROOT / f"scripts/plotting/{helper}.py", "plotting helper") for helper in (
        "publication_review", "plot_style", "plot_series", "plot_manifest", "plot_quality",
        "plot_provenance", "plot_bundle", "plot_accessibility", "plot_delivery",
        "review_vector_delivery", "validate_plot_artifact"))
    if results is not None:
        inputs.append((write_sop_snapshot(results), "production figure SOP snapshot, LF as in Git"))

    fig = render(slices, profile)
    for artist in fig.findobj(match=Text):
        artist.set_fontfamily("Times New Roman")
        artist.set_math_fontfamily("stix")
    fig.set_dpi(600)
    stem = out / STEM
    outputs, quality = export_figure(fig, stem, profile, formats=["png"], intended_width_inches=6.75)
    outputs[0]["role"] = "color_review"
    outputs.append(grayscale_review_record(stem.with_suffix(".png")))

    series, specs = [], []
    for index, temperature in enumerate(TEMPERATURES):
        panel = chr(97 + index)
        picked = [row for row in rows(source)
                  if row["meson"] == "pi" and abs(n(row, "T_MeV") - temperature) < 1e-9]
        positions = [1 + i + 2 for i, _ in enumerate(picked)]
        display = "data/outputs/results/" + source.relative_to(ROOT / "data/outputs/results").as_posix()
        series.append(dict(series_id=f"{panel}.pi.mass.slice", style_ids=["pion_mass.slice_line"],
                           state="computed_native_series", source_csv=display, source_column="mass_MeV",
                           source_rows=positions, point_count=len(positions), panel=panel,
                           support_rule="five frozen pole-search rows at the declared temperature",
                           mask_rule="none"))
        for kind, above in (("bound", False), ("above_threshold", True)):
            members = [i for i, point in enumerate(slices[temperature])
                       if (point["mass"] > point["threshold"]) == above]
            if not members:
                continue
            series.append(dict(series_id=f"{panel}.pi.mass.{kind}_points",
                               style_ids=["pion_mass.pole_point"], state="computed_native_series",
                               source_csv=display, source_column="mass_MeV",
                               source_rows=[positions[i] for i in members], point_count=len(members),
                               panel=panel,
                               support_rule="fill encodes the threshold regime; the sampled value is unchanged",
                               mask_rule="none"))
        specs.append(dict(panel=panel, meson="pi", observables=["mass_MeV"],
                          x_limits=[-0.34, 0.34],
                          y_limits=list({180.0: (130.0, 140.0), 200.0: (130.0, 205.0),
                                         220.0: (160.0, 340.0)}[temperature]),
                          axis_scale="linear", labels=["pole mass"],
                          temperature_MeV=temperature))

    rendering = dict(
        column="double_column", figure_size_inches=quality["figure_size_inches"],
        size_override_reason="three fixed-temperature pion mass slices on one row; Mott-temperature panel dropped",
        quality=quality, placement_limits=placement_limits(quality, profile, outputs),
        legend_policy="declared_geometry_checked",
        legend_placements=[dict(legend_index=0, host_axes_index=None, title="", labels=[r"$M_\pi<M_u+M_d$", r"$M_\pi>M_u+M_d$"],
                                placement="outside_axes", location="figure bottom centre", scope="regime encoding")],
        legend_outside=True,
        color_route=COLOR_ROUTE, delivery_stage="png_review", vector_delivery_pending=True,
        review_contract="prd_review_v1",
        accessibility=accessibility_record(
            [dict(value="bound", color=COLORS[0], noncolor_encoding="open circle"),
             dict(value="above_threshold", color=COLORS[0], noncolor_encoding="filled circle")],
            backgrounds=["#ffffff"], background_scope="opaque white panel backgrounds"),
        path_simplification=False, panel_specs=specs,
        series_styles=register_shared_series_styles(series, REGIMES),
        series_identity_scope="shared style identities refer to the threshold-regime encoding of every panel point",
    )
    record = build_manifest(
        asset_id=STEM, figure_family="xjtu_thesis_chapter3", case_slug=out.name, figure_mode="audit",
        semantic_status="new_prescription_computed_samples", style_profile=profile.profile_id,
        publication_scope="internal_review",
        generator=generator_record(__file__, command=" ".join(sys.argv),
                                   runtime=runtime_record({"matplotlib": matplotlib.__version__})),
        inputs=[input_record(path, role=role) for path, role in inputs], outputs=outputs,
        axes=[dict(field="xi", source_unit="dimensionless", display_unit="dimensionless", label=r"$\xi$"),
              dict(field="mass_MeV", source_unit="MeV", display_unit="MeV", label=r"$M_\pi$ (MeV)")],
        series=series,
        selection_rule="pi pole mass at T=180/200/220 MeV on the five frozen anisotropy points; Mott roots live in tab:mott-temperature",
        interpolation_policy="none", connector_policy="straight guide line between sampled points only",
        missing_value_policy="reject missing/nonfinite/duplicate rows",
        validation=dict(finite=True, duplicate_keys=True, support=True, strict_gate=False),
        rendering=rendering, calculation_sha=numerical["calculation_sha"])
    record.update(manuscript_eligible=False, current_publication_layer=False, delivery_stage="png_review",
                  vector_delivery_pending=True, solver_called=False, canonical_data_modified=False)
    if acceptance is not None:
        add_vector_pdf(record, acceptance, figure=fig, stem=stem, profile=profile, root=ROOT)
    errors = validate_manifest_record(record)
    if errors:
        raise ValueError(STEM + ": " + "; ".join(errors))
    write_manifest(out / "plot_manifest.json", record, newline="\n")
    violations = validate_manifest(out / "plot_manifest.json")
    if violations:
        raise ValueError("; ".join(violations))
    print(json.dumps({"figure": STEM, "glyph_mm": quality["minimum_capital_numeral_height_mm"],
                      "text_overlaps": quality["text_overlap_pairs"],
                      "legend_curves": quality["legend_curve_overlap_count"]}), flush=True)
    plt.close(fig)


if __name__ == "__main__":
    main()
