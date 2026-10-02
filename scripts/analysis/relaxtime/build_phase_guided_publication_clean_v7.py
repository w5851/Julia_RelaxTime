#!/usr/bin/env python3
"""Render v5 display values under the final-size APS v2 plotting contract.

Outputs are review-only v7 figures: 72 single charts and two native composites.
The default build writes the complete candidate APS v2 PDF+PNG set; the
``--png-review`` stage writes a separate PNG-only sibling for author review.
No solver, numerical replacement, or manuscript promotion runs.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import datetime as dt
import importlib.util
import json
from pathlib import Path
import shutil
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

from scripts.plotting.plot_manifest import (
    build_manifest, generator_record, input_record,
    runtime_record, sha256_file, write_manifest,
)
from scripts.plotting.plot_quality import export_figure
from scripts.plotting.plot_style import configure_axis_ticks, configure_matplotlib, load_profile
from scripts.plotting.validate_plot_artifact import validate_manifest

PARENT_SCRIPT = ROOT / "scripts/analysis/relaxtime/build_phase_guided_publication_clean_v6.py"
SPEC = importlib.util.spec_from_file_location("publication_clean_v6_inputs_for_v7", PARENT_SCRIPT)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load v5 input reader: {PARENT_SCRIPT}")
PARENT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PARENT)
V2 = PARENT.V2

FULL_ANALYSIS_ROOT = PARENT.TRANSPORT_ANALYSIS_ROOT / "phase_guided_transport_publication_clean_v7"
FULL_FIGURE_ROOT = ROOT / "data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v7"
PNG_REVIEW_ANALYSIS_ROOT = PARENT.TRANSPORT_ANALYSIS_ROOT / "phase_guided_transport_publication_clean_v7_png_review"
PNG_REVIEW_FIGURE_ROOT = ROOT / "data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v7_png_review"
ANALYSIS_ROOT = FULL_ANALYSIS_ROOT
FIGURE_ROOT = FULL_FIGURE_ROOT
OBSERVABLES = list(PARENT.OBSERVABLES)
UNITS = {
    **{field: "fm" for field in OBSERVABLES if field.startswith("tau_")},
    "eta": "fm^-3", "zeta": "fm^-3", "sigma": "fm^-1",
    "eta_over_s": "dimensionless", "zeta_over_s": "dimensionless", "sigma_over_T": "dimensionless",
}
LABELS = dict(V2.OBSERVABLE_LABELS)
for _field, _unit in UNITS.items():
    if _unit != "dimensionless":
        _math_unit = {"fm": r"\mathrm{fm}", "fm^-3": r"\mathrm{fm}^{-3}", "fm^-1": r"\mathrm{fm}^{-1}"}[_unit]
        LABELS[_field] = LABELS[_field].removesuffix("$") + rf"\;({_math_unit})$"

ENDPOINT_LABELS = {"quark": "1st-order (restored)", "hadron": "1st-order (broken)"}
ENDPOINT_MARKERS = {"quark": "o", "hadron": "s"}
COMPOSITES = {
    "figure1_relaxation_times_comparison": ["tau_u", "tau_s", "tau_ubar", "tau_sbar"],
    "figure2_transport_coefficients_comparison": ["eta_over_s", "zeta_over_s", "sigma_over_T"],
}
PANELS = ["muB0.0", "muB450.0", "muB900.0"]


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def series_label(mode: str, row: dict[str, str]) -> str:
    if mode == "mode_a":
        alpha = float(row["plot_series"].removeprefix("alpha"))
        return rf"$\alpha_T = {alpha:.1f}$"
    return rf"$\mu_B = {float(row['muB_MeV']):.0f}\,\mathrm{{MeV}}$"


def panel_title(mode: str, panel: str) -> str:
    if mode == "mode_a":
        return rf"$\mu_B = {float(panel.removeprefix('muB')):.0f}\,\mathrm{{MeV}}$"
    return rf"$T = {float(panel.removeprefix('T')):.0f}\,\mathrm{{MeV}}$"


def group_inputs(points: list[dict[str, str]], gaps: list[dict[str, Any]]) -> tuple[dict, dict]:
    grouped = defaultdict(list)
    gap_map = defaultdict(list)
    for row in points:
        grouped[(row["mode_key"], row["plot_panel"], row["observable"], row["plot_series"])].append(row)
    for rows in grouped.values():
        rows.sort(key=lambda row: float(row["xi"]))
    for gap in gaps:
        gap_map[(gap["mode_key"], gap["plot_panel"], gap["plot_series"])].append(gap)
    return grouped, gap_map


def render_panel(ax: Any, mode: str, panel: str, observable: str, grouped: dict, gap_map: dict, profile: Any) -> dict[str, Any]:
    names = sorted(key[3] for key in grouped if key[:3] == (mode, panel, observable))
    rows_all = []
    gaps_all = []
    series_records = []
    endpoint_records = []
    for index, name in enumerate(names):
        rows = grouped[(mode, panel, observable, name)]
        gaps = gap_map.get((mode, panel, name), [])
        segments = V2.split_curve_segments(rows, gaps)
        color = profile.colors[index]
        linestyle = profile.data["palette"]["parameter_linestyles"][index]
        for segment in segments:
            ax.plot([float(row["xi"]) for row in segment], [float(row["clean_value"]) for row in segment], color=color, linestyle=linestyle)
        point_index = {V2.canonical_xi(row["xi"]): row for row in rows}
        for gap in gaps:
            for role in ("left", "right"):
                endpoint = gap["endpoint_rows"][role]
                phase = endpoint["phase_curr"]
                if phase not in ENDPOINT_MARKERS:
                    raise ValueError(f"unresolved first-order endpoint: {gap['boundary_id']}")
                point = point_index[V2.canonical_xi(endpoint["xi"])]
                ax.scatter([float(point["xi"])], [float(point["clean_value"])],
                           marker=ENDPOINT_MARKERS[phase], s=profile.data["marker_size_pt"] ** 2,
                           facecolor="white", edgecolor=color, linewidth=0.8, zorder=6)
                endpoint_records.append({"boundary_id": gap["boundary_id"], "role": role, "phase": phase,
                                         "xi": point["xi"], "clean_value": point["clean_value"]})
        series_records.append({
            "series_id": f"{mode}.{panel}.{observable}.{name}",
            "state": "author_accepted_display_derivative", "label": series_label(mode, rows[0]),
            "support_rule": "v5 frozen clean_value rows; no new values or support points",
            "mask_rule": "preserve audited v2 first-order gap; split each gap into separate segments",
            "row_count": len(rows), "segment_count": len(segments),
            "color": color, "linestyle": linestyle,
        })
        rows_all.extend(rows)
        gaps_all.extend(gaps)
    axis_spec = V2.figure_axis_spec(rows_all, gaps_all)
    ax.set_yscale(axis_spec["axis_scale"])
    ax.set_xlim(-0.52, 0.52)
    ax.set_xticks([-0.5, 0, 0.5])
    if axis_spec["axis_scale"] == "linear":
        ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
    configure_axis_ticks(ax, profile)
    ax.set_ylabel(LABELS[observable])
    ax.yaxis.labelpad = 3
    return {"mode_key": mode, "plot_panel": panel, "observable": observable,
            "series": series_records, "endpoints": endpoint_records, **axis_spec}


def legend_handles(spec: dict[str, Any], profile: Any, *, endpoints: bool) -> list[Any]:
    handles = [Line2D([], [], color=item["color"], linestyle=item["linestyle"], label=item["label"]) for item in spec["series"]]
    if endpoints:
        handles.extend(Line2D([], [], linestyle="None", marker=ENDPOINT_MARKERS[phase], markerfacecolor="white",
                              markeredgecolor="black", markersize=profile.data["marker_size_pt"], label=label)
                       for phase, label in ENDPOINT_LABELS.items())
    return handles


def write_chart(figure: Any, stem: Path, specs: list[dict[str, Any]], profile: Any, inputs: list[dict], font: dict,
                command: str, *, composite: bool, delivery_stage: str, output_formats: tuple[str, ...]) -> tuple[dict, Path]:
    outputs, quality = export_figure(figure, stem, profile, formats=output_formats)
    axes_records = []
    for spec in specs:
        axes_records.extend([
            {"field": "xi", "source_unit": "dimensionless", "display_unit": "dimensionless", "label": r"$\xi$", "transform": "identity", "panel": spec["plot_panel"]},
            {"field": spec["observable"], "source_unit": UNITS[spec["observable"]], "display_unit": UNITS[spec["observable"]],
             "label": LABELS[spec["observable"]], "transform": spec["axis_scale"], "panel": spec["plot_panel"]},
        ])
    manifest = build_manifest(
        asset_id=f"relaxtime.phase_guided.v7.{delivery_stage}.{stem.relative_to(FIGURE_ROOT).as_posix()}",
        figure_family="phase_guided_transport", case_slug=stem.name,
        figure_mode="audit", semantic_status="author_review_display_derivative",
        style_profile=profile.profile_id, publication_scope="internal_review",
        generator=generator_record(Path(__file__), command=command, runtime=runtime_record({"matplotlib": matplotlib.__version__, "font": font})),
        inputs=inputs, axes=axes_records, series=[item for spec in specs for item in spec["series"]], outputs=outputs,
        selection_rule="render frozen v5 clean_value table in original xi order; rows/channels/phase gaps unchanged",
        interpolation_policy="inherited v5 display adjustments; no new interpolation or replacement",
        connector_policy="forbidden", missing_value_policy="preserve first-order and missing-support gaps",
        validation={"finite": True, "duplicate_keys": True, "support": True, "strict_gate": False},
        rendering={"column": "double_column", "figure_size_inches": quality["figure_size_inches"],
                   "size_override_reason": "native manuscript composite at final physical size" if composite else None,
                   "bbox_inches": None, "quality": quality, "legend_policy": "shared_external", "legend_outside": True,
                   "color_route": "undecided_review", "grayscale_review": "required before author acceptance",
                   "y_axis_policy": "independent panel ranges; inherit v2 linear/log rule", "panel_specs": specs,
                   "delivery_stage": delivery_stage, "vector_delivery_pending": delivery_stage == "png_review",
                   "output_formats": list(output_formats)},
        calculation_sha=PARENT.V2.V1.CALCULATION_SHA,
    )
    manifest.update({"manuscript_eligible": False, "current_publication_layer": False, "solver_called": False,
                     "canonical_data_modified": False, "new_display_values": False,
                     "delivery_stage": delivery_stage, "vector_delivery_pending": delivery_stage == "png_review",
                     "workflow_head_sha": PARENT.V2.V1.WORKFLOW_HEAD_SHA,
                     "source_run_id_note": "multiple upstream runs; see frozen v5 per-point run_id and raw manifests",
                     "numerical_status": "inherited_author_accepted_display_only", "raw_manuscript_eligible": False})
    manifest_path = stem.with_suffix(".plot_manifest.json")
    write_manifest(manifest_path, manifest)
    violations = validate_manifest(manifest_path)
    if violations:
        raise ValueError(f"{stem.name}: " + "; ".join(violations))
    plt.close(figure)
    return {"stem": relative(stem), "manifest": relative(manifest_path), "manifest_sha256": sha256_file(manifest_path),
            "kind": "composite" if composite else "single", "mode_key": specs[0]["mode_key"], "outputs": outputs,
            "minimum_capital_numeral_height_mm": quality["minimum_capital_numeral_height_mm"]}, manifest_path


def render_single(mode: str, panel: str, observable: str, grouped: dict, gap_map: dict, profile: Any) -> tuple[Any, list[dict]]:
    figure, ax = plt.subplots(figsize=(6.75, 4.6))
    figure.subplots_adjust(left=0.16, right=0.97, bottom=0.14, top=0.73)
    spec = render_panel(ax, mode, panel, observable, grouped, gap_map, profile)
    ax.set_xlabel(r"$\xi$")
    figure.suptitle(panel_title(mode, panel), y=0.97)
    handles = legend_handles(spec, profile, endpoints=bool(spec["endpoints"]))
    figure.legend(handles=handles[:3], loc="upper center", bbox_to_anchor=(0.5, 0.915), ncol=3, handlelength=2.1, columnspacing=1)
    if len(handles) > 3:
        figure.legend(handles=handles[3:], loc="upper center", bbox_to_anchor=(0.5, 0.835), ncol=2, handlelength=1, columnspacing=1)
    return figure, [spec]


def render_composite(observables: list[str], grouped: dict, gap_map: dict, profile: Any) -> tuple[Any, list[dict]]:
    nrows = len(observables)
    figure, axes = plt.subplots(nrows, 3, figsize=(6.75, 8.2 if nrows == 4 else 6.8), sharex=True)
    figure.subplots_adjust(left=0.12, right=0.98, bottom=0.095, top=0.83, hspace=0.26, wspace=0.48)
    specs = []
    for row, observable in enumerate(observables):
        for col, panel in enumerate(PANELS):
            ax = axes[row, col]
            spec = render_panel(ax, "mode_a", panel, observable, grouped, gap_map, profile)
            specs.append(spec)
            if col:
                ax.set_ylabel("")
            if row == 0:
                ax.set_title(panel_title("mode_a", panel), pad=19)
            if row == nrows - 1:
                ax.set_xlabel(r"$\xi$")
            ax.text(0.01, 1.015, f"({chr(97 + row * 3 + col)})", transform=ax.transAxes, ha="left", va="bottom")
    handles = legend_handles(specs[0], profile, endpoints=True)
    figure.legend(handles=handles[:3], loc="upper center", bbox_to_anchor=(0.52, 0.995), ncol=3, handlelength=2.1)
    figure.legend(handles=handles[3:], loc="upper center", bbox_to_anchor=(0.52, 0.951), ncol=2, handlelength=1)
    return figure, specs


def parent_records(package: dict, profile: Any) -> list[dict]:
    paths = {PARENT.V5_MANIFEST: "accepted_display_parent", PARENT.V5_PLOT_MANIFEST: "parent_plot_manifest",
             PARENT.V5_FIGURE_LAYER_MANIFEST: "parent_figure_layer", PARENT.CURRENT_POINTER: "retained_current_pointer",
             profile.path: "style_profile", PARENT_SCRIPT: "input_reader",
             ROOT / "scripts/plotting/plot_style.py": "shared_style_code",
             ROOT / "scripts/plotting/plot_manifest.py": "shared_manifest_code",
             ROOT / "scripts/plotting/plot_quality.py": "shared_quality_code",
             ROOT / "scripts/plotting/validate_plot_artifact.py": "shared_validator_code"}
    paths.update({path: "inherited_display_table" for path in PARENT.V5_TABLE_ROOT.glob("*.csv")})
    for module in (PARENT.V5, PARENT.V4, PARENT.V3, PARENT.V2, PARENT.V2.V1):
        paths[Path(module.__file__)] = "inherited_input_or_gap_code"
    paths[PARENT.V4.PARENT_MANIFEST] = "gap_provenance_parent_manifest"
    paths[PARENT.V4.PARENT_POINTS] = "gap_provenance_parent_points"
    for source in package["source_inputs"]:
        for role, path in PARENT.V2.V1.case_paths(source["mode_key"]).items():
            if role != "result_dir":
                paths[path] = f"raw_provenance_{role}"
    return [input_record(path, role=role) for path, role in paths.items()]


def write_caption_notes(points: list[dict]) -> Path:
    mapping = {}
    for row in points:
        if row["mode_key"] != "mode_a":
            continue
        key = (row["plot_panel"], row["plot_series"])
        value = float(row["T_MeV"])
        if key in mapping and mapping[key] != value:
            raise ValueError(f"temperature varies along xi for {key}; cannot move T to caption")
        mapping[key] = value
    table = "\n".join(f"| {panel} | " + " | ".join(f"{mapping[(panel, f'alpha{alpha:.1f}')]:.6f}" for alpha in (1, 1.1, 1.2)) + " |" for panel in PANELS)
    path = ANALYSIS_ROOT / "caption_handoff.md"
    path.write_text("""# v7 display-only caption handoff

No paper project files have been changed. Keep Figure 1's rows as
tau_u, tau_s, tau_ubar, tau_sbar (including antiquark bars); Figure 2's rows
are eta/s, zeta/s, sigma/T. Columns are mu_B = 0, 450, 900 MeV.

Suggested additions to the existing scientific captions:

> The curves correspond to alpha_T = 1.0, 1.1, and 1.2, respectively.
> Open circles and squares denote the chirally restored and chirally broken
> endpoints of the first-order transition, respectively. Disconnected
> segments display the two phase branches separately. Different vertical
> scales are used among panels. The relaxation times are given in fm.

Use mathematical notation in LaTeX. For alpha_T = (1.0, 1.1, 1.2), the
temperatures in MeV are approximately (200, 220, 240), (183, 201, 219), and
(126, 138, 151) in the left, middle, and right columns, respectively.
Do not imply a temperature shared across columns.

Exact frozen temperature mapping (MeV):

| panel | alpha 1.0 | alpha 1.1 | alpha 1.2 |
| --- | ---: | ---: | ---: |
""" + table + "\n\nRetain the existing disclosure of display-only local replacements and the raw/adjustment tables. v7 changes no numerical values and does not certify raw convergence.\n", encoding="utf-8")
    return path


def main() -> None:
    global ANALYSIS_ROOT, FIGURE_ROOT

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--composites-only", action="store_true", help="temporary layout probe; the full review build includes 72 singles")
    parser.add_argument(
        "--png-review",
        action="store_true",
        help="generate the PNG-only v7 author-review stage; vector delivery remains pending",
    )
    parser.add_argument("--output-root", type=Path, default=FIGURE_ROOT)
    args = parser.parse_args()
    if args.png_review:
        if args.composites_only:
            raise ValueError("--png-review cannot be combined with --composites-only")
        ANALYSIS_ROOT = PNG_REVIEW_ANALYSIS_ROOT
        FIGURE_ROOT = PNG_REVIEW_FIGURE_ROOT
        delivery_stage = "png_review"
        output_formats = ("png",)
    else:
        delivery_stage = "vector_delivery"
        output_formats = ("pdf", "png")
    if args.composites_only:
        if args.output_root == FIGURE_ROOT:
            raise ValueError("development preview requires a separate --output-root")
    elif args.output_root.resolve() != FIGURE_ROOT.resolve():
        if not args.png_review:
            raise ValueError("full v7 build uses the versioned review root")
    output_root = args.output_root if args.composites_only else FIGURE_ROOT
    if output_root.exists() or (not args.composites_only and ANALYSIS_ROOT.exists()):
        raise FileExistsError("refusing to overwrite an existing v7/preview review package")
    points, package, _, _, gaps, table_hashes = PARENT.load_v5_inputs()
    keys = [(row["mode_key"], row["plot_panel"], row["plot_series"], row["observable"], row["xi"]) for row in points]
    if len(keys) != len(set(keys)):
        raise ValueError("duplicate frozen v5 point keys")
    for item in package["outputs"]:
        path = ROOT / item["path"]
        if sha256_file(path) != item["sha256"]:
            raise ValueError(f"parent package hash mismatch: {path}")
    profile = load_profile("candidate_aps_v2")
    font = configure_matplotlib(profile)
    grouped, gap_map = group_inputs(points, gaps)
    if args.composites_only:
        from scripts.plotting.plot_quality import measure_figure
        args.output_root.mkdir(parents=True)
        for name, observables in COMPOSITES.items():
            figure, _ = render_composite(observables, grouped, gap_map, profile)
            quality = measure_figure(figure, intended_width_inches=6.75)
            figure.savefig(args.output_root / f"{name}.png", dpi=150, bbox_inches=None)
            print(json.dumps({"figure": name, "quality": quality}, ensure_ascii=True))
            plt.close(figure)
        return
    ANALYSIS_ROOT.mkdir(parents=True)
    shutil.copytree(PARENT.V5_TABLE_ROOT, ANALYSIS_ROOT / "tables")
    inputs = parent_records(package, profile)
    command = "python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v7.py"
    if args.png_review:
        command += " --png-review"
    charts = []
    for mode in ("mode_a", "mode_b"):
        panels = sorted({key[1] for key in grouped if key[0] == mode})
        for panel in panels:
            for observable in OBSERVABLES:
                figure, specs = render_single(mode, panel, observable, grouped, gap_map, profile)
                mode_dir = "mode_a_fixed_muB_phase_scaled" if mode == "mode_a" else "mode_b_fixed_T_sparse_muB"
                stem = FIGURE_ROOT / mode_dir / f"plot_panel={panel}" / f"{observable}_vs_xi"
                chart, _ = write_chart(
                    figure, stem, specs, profile, inputs, font, command,
                    composite=False, delivery_stage=delivery_stage, output_formats=output_formats,
                )
                charts.append(chart)
        print(f"[v7] {mode}: 36 single charts validated", flush=True)
    for name, observables in COMPOSITES.items():
        figure, specs = render_composite(observables, grouped, gap_map, profile)
        chart, _ = write_chart(
            figure, FIGURE_ROOT / "composites" / name, specs, profile, inputs, font, command,
            composite=True, delivery_stage=delivery_stage, output_formats=output_formats,
        )
        charts.append(chart)
    for item in inputs:
        if sha256_file(ROOT / item["path"]) != item["sha256"]:
            raise ValueError(f"input changed during v7 render: {item['path']}")
    for name, expected in table_hashes.items():
        if sha256_file(ANALYSIS_ROOT / "tables" / name) != expected:
            raise ValueError(f"v7 changed inherited table: {name}")
    notes = write_caption_notes(points)
    style_path = ANALYSIS_ROOT / "v7_display_style_manifest.json"
    style_schema = "publication_clean_v7_png_review_display_style_v1" if args.png_review else "publication_clean_v7_display_style_v1"
    style = {"schema": style_schema, "style_profile": profile.profile_id,
             "profile_sha256": sha256_file(profile.path), "font": font, "font_size_pt": profile.data["font_size_pt"],
             "formats": list(output_formats), "profile_formats": list(profile.formats), "unit_style": "parentheses", "endpoint_labels": ENDPOINT_LABELS,
             "series_legend": "alpha_T only; temperatures in caption_handoff.md",
             "composite_strategy": "native vector axes from frozen points; no PNG resizing or raster-wrapped PDF",
             "ticks": profile.data["ticks"], "manuscript_eligible": False,
             "delivery_stage": delivery_stage, "vector_delivery_pending": delivery_stage == "png_review",
             "parent_manifest_sha256": sha256_file(PARENT.V5_MANIFEST)}
    write_manifest(style_path, style)
    index_schema = "publication_clean_v7_png_review_figure_index_v1" if args.png_review else "publication_clean_v7_figure_index_v1"
    index = {"schema": index_schema, "status": "author_review_required",
             "manuscript_eligible": False, "current_publication_layer": False, "solver_called": False,
             "delivery_stage": delivery_stage, "vector_delivery_pending": delivery_stage == "png_review",
             "single_figure_count": 72, "mode_counts": {"mode_a": 36, "mode_b": 36}, "composite_figure_count": 2,
             "style_manifest": relative(style_path), "style_manifest_sha256": sha256_file(style_path), "charts": charts}
    write_manifest(FIGURE_ROOT / "plot_manifest.json", index)
    readme = ANALYSIS_ROOT / "README.md"
    readme.write_text(f"""# publication_clean_v7 {'PNG-only author review' if args.png_review else 'final-size plotting review'}

v7 consumes the accepted v5 clean_value table without adding display values.
All inherited CSV tables are byte-identical. The 72 single charts and two
mode-A manuscript composites use candidate_aps_v2, {', '.join(output_formats)},
parentheses units, coherent typography, four-sided inward major/minor ticks,
and color plus line styles. Composite legends are shared outside the data
axes. Temperatures and full endpoint explanations are in caption_handoff.md.
Y-axis ranges remain independent and the audited v2 linear/log choice and
phase gaps are preserved.

Every chart has a shared-contract .plot_manifest.json with measured physical
dimensions, actual glyph heights, and hashes. PDF-specific embedded-font and
vector checks apply only after the vector delivery stage.
The figure-root plot_manifest.json is the index, not a single-chart contract.
Use the per-chart manifests with validate_plot_artifact.py.

Status: author_review_required, delivery_stage={delivery_stage},
manuscript_eligible=false, current_publication_layer=false,
vector_delivery_pending={delivery_stage == 'png_review'}. The accepted
v5/current pointer, raw data, solver, numerical gates, production registry,
and paper project are unchanged.
Color-print versus online-color/grayscale-print remains an author choice;
the latter requires additional PS/EPS production delivery. PNG review is not
the final vector delivery.

Reproduce in a new/absent review directory:
python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v7.py{' --png-review' if args.png_review else ''}
""", encoding="utf-8")
    package_schema = "publication_clean_v7_png_review_package_v1" if args.png_review else "publication_clean_v7_package_v1"
    manifest = {"schema": package_schema, "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
                "status": "author_review_required", "manuscript_eligible": False, "current_publication_layer": False,
                "delivery_stage": delivery_stage, "vector_delivery_pending": delivery_stage == "png_review",
                "author_acceptance": None, "solver_called": False, "canonical_data_modified": False,
                "numerical_status": "inherited_author_accepted_display_only", "raw_manuscript_eligible": False,
                "calculation_sha": package.get("calculation_sha"), "workflow_head_sha": package.get("workflow_head_sha"),
                "parent_manifest": relative(PARENT.V5_MANIFEST), "parent_manifest_sha256": sha256_file(PARENT.V5_MANIFEST),
                "inputs": inputs, "point_row_count": len(points), "inherited_table_hashes": table_hashes,
                "figure_index": relative(FIGURE_ROOT / "plot_manifest.json"), "figure_index_sha256": sha256_file(FIGURE_ROOT / "plot_manifest.json"),
                "outputs": [input_record(path, role="v7_display_document_or_table") for path in [readme, notes, style_path, *sorted((ANALYSIS_ROOT / "tables").glob("*.csv"))]]}
    write_manifest(ANALYSIS_ROOT / "manifest.json", manifest)
    print(json.dumps({"figures": relative(FIGURE_ROOT), "manifest": relative(ANALYSIS_ROOT / "manifest.json"),
                      "single_charts": 72, "composites": 2, "delivery_stage": delivery_stage,
                      "vector_delivery_pending": delivery_stage == "png_review", "manuscript_eligible": False,
                      "minimum_capital_numeral_height_mm": min(chart["minimum_capital_numeral_height_mm"] for chart in charts)}))


if __name__ == "__main__":
    main()
