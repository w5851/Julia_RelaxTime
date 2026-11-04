#!/usr/bin/env python3
"""Compare full-RS publication layout with the explicitly requested legacy forms."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.plotting.plot_bundle import expand_bundle
from scripts.plotting.plot_manifest import input_record, write_manifest
from scripts.analysis.relaxtime import build_phase_guided_publication_clean_v14 as delivery


def indexed(path):
    case = path.parent.name
    pairs = expand_bundle(json.loads(path.read_text(encoding="utf-8")))
    return {chart.get("relative_stem", chart["stem"].split(case + "/")[-1]): record for chart, record in pairs}


def panels(record):
    return [(p["plot_panel"], p["observable"], p["axis_scale"]) for p in record["rendering"]["panel_specs"]]


def legend_form(record):
    return [(p["role"], p["placement"], p["host_axes_index"], p["location"], p["title"], p["labels"])
            for p in record["rendering"]["legend_placements"]]


def palette_form(record):
    return [[(curve["color"].lower(), curve["linestyle"]) for curve in panel["series"]]
            for panel in record["rendering"]["panel_specs"]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    old_path = delivery.LEGACY_ROOT / "publication_clean_v13_png_review/plot_manifest.json"
    fixed_path = delivery.LEGACY_ROOT / "transport_fixed_T_panel_legends_20261021__png_review/plot_manifest.json"
    old, new = indexed(old_path), indexed(args.manifest)
    issues = []
    checks = []
    for name, record in old.items():
        actual = new[name]
        same = dict(size=record["rendering"]["figure_size_inches"] == actual["rendering"]["figure_size_inches"],
                    panel_order_and_scales=panels(record) == panels(actual),
                    legend_positions_and_labels=legend_form(record) == legend_form(actual),
                    colors_and_line_styles=palette_form(record) == palette_form(actual))
        checks.append(dict(figure=name, **same))
        if not all(same.values()):
            issues.append(dict(figure=name, failed=[key for key, value in same.items() if not value]))
    fixed = json.loads(fixed_path.read_text(encoding="utf-8"))
    new_fixed = new[delivery.FIG4.as_posix()]
    expected_legend_hosts = {p["title"]: (ord(p["panel"]) - 97, p["anchor_axes"], p["labels"])
                             for p in fixed["rendering"]["legend_layout"]}
    actual_legend_hosts = {p["title"]: (p["host_axes_index"], list(p["anchor_axes"]), p["labels"])
                           for p in new_fixed["rendering"]["legend_placements"]}
    fixed_checks = dict(size=fixed["rendering"]["figure_size_inches"] == new_fixed["rendering"]["figure_size_inches"],
                        legend_hosts_anchors_and_labels=expected_legend_hosts == actual_legend_hosts,
                        panel_order_and_scales=panels(new_fixed) == [(panel, observable, scale)
                            for observable, scale in (("eta_over_s", "log"), ("zeta_over_s", "log"), ("sigma_over_T", "linear"))
                            for panel in ("T120.0", "T160.0", "T200.0")],
                        raw_vertices_verified=new_fixed["validation"]["original_vertex_count"] == 2727)
    if not all(fixed_checks.values()):
        issues.append(dict(figure=str(delivery.FIG4), failed=[key for key, value in fixed_checks.items() if not value]))
    paper_figure = Path("D:/Desktop/paper/My_Paper/pnjl_aniso_transport/writing_package/outputs/journal_en/"
                        "review_revision_2026-10-12/figure_fixed_T/figure4_transport_fixed_T.pdf")
    result = dict(schema="phase_guided_publication_form_comparison_v1", status="passed" if not issues else "failed",
        v13_figures_checked=len(checks), v13_form_checks=checks, fixed_T_form_checks=fixed_checks, issues=issues,
        permitted_data_dependent_differences=[
            "Axis bounds and numerical tick labels follow the new original data.",
            "First-order marker coordinates follow this grid's audited sampled branch pair.",
            "The first two detail viewports extend from 0.58 to 0.57 fm, preserving all original path vertices.",
            "The narrow FIG4 zeta/s log panel uses spaced ordinary-number tick labels to prevent overlap."],
        raw_review_alignment=["Restore legacy green #009e73.",
            "Restore the two mode-A muB=0 sigma/sigma-over-T single-figure legend corners.",
            "Restore the mode-B T=160 eta/s single figure to a linear axis."],
        manuscript_FIG4_context="main.tex still references the older 6.75x6.9 inch top-legend PDF; this task explicitly requests the 20261021 6.75x7.4 inch in-panel-legend form.",
        references=[input_record(old_path, role="v13_layout"), input_record(fixed_path, role="requested_fixed_T_layout"),
                    input_record(paper_figure, role="current_manuscript_FIG4_read_only"),
                    input_record(args.manifest, role="new_review_bundle")])
    write_manifest(args.output, result)
    print(json.dumps({key: value for key, value in result.items() if key not in ("v13_form_checks", "references")}, indent=2))
    if issues:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
