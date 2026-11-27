#!/usr/bin/env python3
"""Render the full-grid single-panel diagnostic set of the Ward/Møller transport layer.

The reviewed manuscript figures are four composites; this generator produces the
matching **single-panel** breakdown (6 slices × 12 observables = 72 figures) from
the same frozen layer, with no solver call and no numeric recomputation.  It reads
`mode_a_grid.csv` / `mode_b_grid.csv` exactly as the manuscript generator does and
reuses `scripts/plotting/transport_panels.py` for curves, endpoints, axis scale and
tick capture, so a panel and its composite parent cannot drift apart.

Semantics: the inputs are the **medium-rest-frame Møller** layer
(`relative_velocity = "medium_moller"`), not the retired centre-of-mass
`publication_clean_v15` tables.  Nothing here is derived from v15.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

# Run from any working directory: the plotting helpers live under the repository
# root, which is not on sys.path when this file is invoked by path.
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import matplotlib.pyplot as plt
from matplotlib.text import Text
from matplotlib.ticker import FixedLocator, LogLocator

from scripts.plotting.plot_style import load_profile
from scripts.plotting.plot_manifest import (build_manifest, input_record, generator_record, runtime_record,
                                            write_manifest)
from scripts.plotting.plot_quality import export_figure, placement_limits
from scripts.plotting.plot_accessibility import accessibility_record
from scripts.plotting.plot_bundle import build_bundle
from scripts.plotting.validate_plot_artifact import validate_manifest, validate_manifest_record
from scripts.plotting import transport_panels as panels
from scripts.plotting.publication_review import (figure_context, grayscale_review_record,
                                                 register_shared_series_styles, write_sop_snapshot)
from scripts.plotting.review_vector_delivery import COLOR_ROUTE, require_fresh_case

ROOT = Path(__file__).resolve().parents[3]
COLUMN = "single_column"
HELPERS = ("publication_review", "transport_panels", "plot_style", "plot_series", "plot_manifest",
           "plot_quality", "plot_provenance", "plot_bundle", "plot_accessibility", "plot_delivery",
           "review_vector_delivery", "validate_plot_artifact")


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def transport_inputs(mode_a: Path, mode_b: Path) -> panels.TransportInputs:
    groups = {**panels.validate_rows(rows(mode_a), "mode_a"), **panels.validate_rows(rows(mode_b), "mode_b")}
    gaps = {("mode_a", "muB900.0", "alpha1.0"): [(-0.003, 0.003)],
            ("mode_b", "T120.0", "muB900.0"): [(-0.13, -0.12)]}
    return panels.TransportInputs(groups, gaps, {"mode_a": mode_a, "mode_b": mode_b})


def stem_for(mode: str, panel: str, observable: str) -> str:
    return f"{mode}_{panel}_{observable}".replace(".", "")


def declared_legends(quality: dict, scope: str) -> list[dict]:
    return [dict(legend_index=item["legend_index"], host_axes_index=item["host_axes_index"],
                 title=item["title"], labels=item["labels"],
                 placement="in_axes" if item["host_axes_index"] is not None else "outside_axes",
                 location=item.get("location", "declared clear panel corner"), scope=scope)
            for item in quality["legend_layout"]]


def export_record(figure, stem: Path, profile, input_paths, calculation_sha, *, description, axes, series,
                  registry, placements, intended_width: float) -> dict:
    # Legends added after a bound style context inherit Matplotlib defaults; bind
    # their ordinary and mathematical faces to the same manuscript fonts.
    for artist in figure.findobj(match=Text):
        artist.set_fontfamily("Times New Roman")
        artist.set_math_fontfamily("stix")
    figure.set_dpi(600)
    outputs, quality = export_figure(figure, stem, profile, formats=["png"],
                                     intended_width_inches=intended_width)
    outputs[0]["role"] = "color_review"
    outputs.append(grayscale_review_record(stem.with_suffix(".png")))
    rendering = dict(column=COLUMN, figure_size_inches=quality["figure_size_inches"],
                     size_override_reason=description, quality=quality,
                     placement_limits=placement_limits(quality, profile, outputs),
                     legend_policy="declared_geometry_checked", legend_placements=placements,
                     legend_outside=all(not item["overlapped_axes"] for item in quality["legend_layout"]),
                     color_route=COLOR_ROUTE, delivery_stage="png_review", vector_delivery_pending=True,
                     review_contract="prd_review_v1",
                     accessibility=accessibility_record(
                         [dict(value=str(index), color=style.color, noncolor_encoding=style.linestyle)
                          for index, style in enumerate(registry.styles)],
                         backgrounds=["#ffffff"], background_scope="opaque white panel backgrounds"),
                     path_simplification=False, panel_specs=[], series_styles=register_shared_series_styles(series, registry),
                     series_identity_scope="shared style identities refer to the panel series and selected frozen rows")
    record = build_manifest(
        asset_id="ward_moller_" + stem.name, figure_family="pnjl_transport_full_grid_panels",
        case_slug=stem.parent.name, figure_mode="audit", semantic_status="new_prescription_computed_samples",
        style_profile=profile.profile_id, publication_scope="internal_review",
        generator=generator_record(__file__, command=" ".join(sys.argv),
                                   runtime=runtime_record({"matplotlib": matplotlib.__version__})),
        inputs=[input_record(path, role=role) for path, role in input_paths], outputs=outputs,
        axes=axes, series=series, selection_rule=description, interpolation_policy="none",
        connector_policy="forbidden", missing_value_policy="reject missing/nonfinite/duplicate rows",
        validation=dict(finite=True, duplicate_keys=True, support=True, strict_gate=False),
        rendering=rendering, calculation_sha=calculation_sha)
    record.update(manuscript_eligible=False, current_publication_layer=False, delivery_stage="png_review",
                  vector_delivery_pending=True, solver_called=False, canonical_data_modified=False)
    errors = validate_manifest_record(record)
    if errors:
        raise ValueError(stem.name + ": " + "; ".join(errors))
    print(json.dumps({"figure": stem.name, "glyph_mm": quality["minimum_capital_numeral_height_mm"],
                      "text_overlaps": quality["text_overlap_pairs"],
                      "legend_curves": quality["legend_curve_overlap_count"]}), flush=True)
    plt.close(figure)
    return record


def series_from_spec(spec: dict) -> list[dict]:
    return [dict(curve, series_id=curve["curve_id"], style_ids=[curve["series_id"]],
                 state="computed_native_series",
                 support_rule="computed source rows, split only at the declared first-order pair",
                 mask_rule="none")
            for curve in spec["series"]]


def render_single(mode: str, panel: str, observable: str, data, profile, *, tick_count: int | None = None):
    size = tuple(float(value) for value in profile.data["figure_size_in"][COLUMN])
    figure, axis = plt.subplots(figsize=size)
    # Reserve the margins a single-column panel needs for a rotated y label, the
    # $\xi$ label and the tick furniture before legends are measured.
    figure.subplots_adjust(left=0.235, right=0.975, bottom=0.215, top=0.965)
    spec = panels.render_panel(axis, mode, panel, observable, data, profile)
    axis.set_xlabel(r"$\xi$", labelpad=2)

    base_ticks = list(axis.get_yticks())
    # A single-column panel is only 2.5 in tall.  On a wide decade span the
    # composite's readable log tick set crowds the labels once a legend is placed,
    # so thin the ticks until the declared layout is genuinely clear.  The
    # admission gate is never relaxed — only the tick density adapts.
    counts = (tick_count,) if tick_count is not None else ((5, 4, 3) if spec["axis_scale"] == "log" else (None,))
    last: ValueError | None = None
    for count in counts:
        if count is not None and len(base_ticks) > count:
            picked = sorted({base_ticks[round(index * (len(base_ticks) - 1) / (count - 1))]
                             for index in range(count)})
            axis.yaxis.set_major_locator(FixedLocator(picked))
            # The readable tick set sits on non-decade values (0.5, 5, 50, ...),
            # which suppresses matplotlib's automatic decade minors.  Restore an
            # explicit minor locator so all four sides keep inward minor ticks.
            axis.yaxis.set_minor_locator(LogLocator(base=10, subs=tuple(range(2, 10)), numticks=100))
            panels.capture_ticks(axis, spec)
        try:
            placements = panels.place_legends(figure, [spec], profile, hosts=(0, 0), single=True)
            return figure, spec, placements
        except ValueError as exc:
            last = exc
            for legend in panels.visible_legends(figure):
                legend.remove()
    raise ValueError(f"{mode}/{panel}/{observable}: {last}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True, help="figure case directory (must not exist)")
    parser.add_argument("--results", type=Path, required=True, help="results case directory (must not exist)")
    parser.add_argument("--mode", choices=("mode_a", "mode_b"), help="restrict to one mode")
    parser.add_argument("--panel", help="restrict to one panel, for example muB450.0")
    parser.add_argument("--observable", help="restrict to one observable, for example tau_s")
    parser.add_argument("--manifest", type=Path, help="read the frozen layer from this directory instead")
    args = parser.parse_args()

    out, results = args.out.resolve(), args.results.resolve()
    require_fresh_case(out, results)
    profile = load_profile("candidate_aps_v2")
    layer = (args.manifest or (ROOT / "data" / "outputs" / "results" / "relaxtime" / "transport" /
                               "ward_moller" / "ward_moller_sphere_20261121_v1")).resolve()
    numerical = json.loads((layer / "numerical_manifest.json").read_text(encoding="utf-8"))
    data = transport_inputs(layer / "mode_a_grid.csv", layer / "mode_b_grid.csv")

    out.mkdir(parents=True, exist_ok=True)
    inputs = [(layer / "numerical_manifest.json", "frozen numerical provenance"),
              (profile.path, "plot profile")]
    inputs.extend((ROOT / f"scripts/plotting/{helper}.py", "plotting helper") for helper in HELPERS)
    inputs.append((write_sop_snapshot(results), "production figure SOP snapshot, LF as in Git"))

    records, charts = [], []
    modes = (args.mode,) if args.mode else tuple(panels.PANELS)
    for mode in modes:
        for panel in panels.PANELS[mode]:
            if args.panel and panel != args.panel:
                continue
            for observable in panels.OBSERVABLES:
                if args.observable and observable != args.observable:
                    continue
                stem = out / stem_for(mode, panel, observable)
                sources = list(inputs)
                sources.extend([(data.sources[mode], "recomputed phase-referenced transport"
                                 if mode == "mode_a" else "recomputed fixed-temperature transport")])
                axis_record = [dict(field="xi", source_unit="dimensionless", display_unit="dimensionless",
                                    label=r"$\xi$"),
                               dict(field=observable,
                                    source_unit="fm" if observable.startswith("tau_") else
                                    ("fm-3" if observable in ("eta", "zeta") else
                                     ("fm-1" if observable == "sigma" else "dimensionless")),
                                    display_unit="fm" if observable.startswith("tau_") else
                                    ("fm-3" if observable in ("eta", "zeta") else
                                     ("fm-1" if observable == "sigma" else "dimensionless")),
                                    label=panels.LABELS[observable])]
                registry = panels.make_registry(mode, profile)
                # Wide tick labels (for example the small eta values at muB=0) can
                # push a single-column panel below the v2 glyph/placement gate.  The
                # SOP allows a documented size override, so retry at a larger size
                # rather than loosening the gate.
                base_width = float(profile.data["figure_size_in"][COLUMN][0])
                # η, ζ and σ carry superscript unit glyphs (fm^-3, fm^-1).  At the
                # nominal single-column width their superscripts measure below the
                # v2 2 mm floor, so declare a wider intended width (the SOP's
                # documented size override) instead of loosening the gate.
                # Two independent admission failures need different remedies:
                # crowded tick labels need fewer ticks, small script glyphs need a
                # wider declared width.  Walk a small explicit plan.
                plan = ((None, base_width), (4, base_width), (3, base_width), (3, 3.9))
                for tick_count, intended in plan:
                    # Draw and export inside the profile's style context: without it
                    # matplotlib falls back to 7/10 pt defaults, which drops
                    # superscript unit glyphs below the 2 mm admission floor.
                    with figure_context(profile, compact=True):
                        figure, spec, placements = render_single(mode, panel, observable, data, profile,
                                                                 tick_count=tick_count)
                        try:
                            record = export_record(
                                figure, stem, profile, sources, numerical["calculation_sha"],
                                description=(f"{panels.panel_title(mode, panel)}; {observable}; "
                                             "medium-rest-frame Moller layer, no recomputation"),
                                axes=axis_record, series=series_from_spec(spec), registry=registry,
                                placements=placements, intended_width=intended)
                            break
                        except ValueError as exc:
                            plt.close(figure)
                            # A rejected attempt may already have written its PNG; the
                            # retry must start from a clean stem or the exporter
                            # refuses to overwrite it.
                            for leftover in (stem.with_suffix(".png"),
                                             stem.with_name(stem.name + "_gray.png")):
                                leftover.unlink(missing_ok=True)
                            message = str(exc)
                            retryable = ("glyph height is below" in message
                                         or "adjacent tick labels closer" in message
                                         or "placement lies outside its usable width" in message)
                            if not retryable or (tick_count, intended) == plan[-1]:
                                raise
                            print(json.dumps({"figure": stem.name, "retry": [tick_count, intended],
                                              "reason": message[:150]}), flush=True)
                records.append(record)
                charts.append(dict(name=stem.name, asset_id=record["asset_id"]))

    if not records:
        parser.error("no panel matched the requested filters")
    manifest = out / "plot_manifest.json"
    write_manifest(manifest, build_bundle(dict(case=out.name, charts=charts), records), newline="\n")
    violations = validate_manifest(manifest)
    if violations:
        raise ValueError("; ".join(violations))
    print(json.dumps({"panels": len(records), "manifest": str(manifest)}), flush=True)


if __name__ == "__main__":
    main()
