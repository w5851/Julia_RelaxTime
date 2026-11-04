#!/usr/bin/env python3
"""Deliver full-RS raw transport figures in the established publication layout."""

from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import replace
import hashlib
import io
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoLocator, FixedLocator, FuncFormatter, LogLocator, MaxNLocator, NullFormatter, FormatStrFormatter, ScalarFormatter
from PIL import Image

from scripts.analysis.relaxtime import build_phase_guided_full_rs_raw as raw
from scripts.plotting.plot_bundle import build_bundle, expand_bundle
from scripts.plotting.plot_manifest import generator_record, input_record, output_record, runtime_record, sha256_file, write_manifest
from scripts.plotting.plot_quality import export_figure, inspect_export, measure_figure, placement_limits, visible_legends
from scripts.plotting.plot_style import configure_axis_ticks, figure_style_context, load_profile
from scripts.plotting.validate_plot_artifact import validate_manifest, validate_manifest_record

REVIEW_CASE = "publication_clean_v14_png_review"
FINAL_CASE = "publication_clean_v14_pdf"
LEGACY_ROOT = Path("D:/Desktop/Julia_RelaxTime/data/outputs/figures/relaxtime/transport/phase_guided")
FIG4 = Path("composites/figure4_transport_fixed_T")
LEGACY_COLORS = ("#0072b2", "#d55e00", "#009e73")
SPECIAL_SINGLES = {
    ("mode_a", "muB0.0", "sigma"): "lower left",
    ("mode_a", "muB0.0", "sigma_over_T"): "upper left",
    ("mode_b", "T160.0", "eta_over_s"): None,
}


def publication_profile():
    profile = load_profile("strict_aps_v2")
    data = deepcopy(dict(profile.data))
    data["palette"]["colors"] = [*LEGACY_COLORS, *profile.colors[3:]]
    return replace(profile, data=data)


def parameter_legend(figure, specs, profile, *, host, location, anchor=None, size=11):
    mode = specs[0]["mode_key"]
    registry = raw.make_registry(mode, profile)
    labels = ["1.0", "1.1", "1.2"] if mode == "mode_a" else ["0", "450", "900"]
    title = r"$\alpha_T$" if mode == "mode_a" else r"$\mu_B\;(\mathrm{MeV})$"
    handles = [Line2D([], [], color=style.color, ls=style.linestyle, lw=1.2, label=label)
               for style, label in zip(registry.styles, labels)]
    settings = dict(loc=location, frameon=False, alignment="left", fontsize=size, title_fontsize=13,
                    handlelength=1.35, handletextpad=0.3, labelspacing=0.10, borderpad=0.05, borderaxespad=0.25)
    if anchor is not None:
        settings.update(bbox_to_anchor=anchor, borderaxespad=0, borderpad=0.25, labelspacing=0.20,
                        handlelength=1.6, handletextpad=0.4)
    figure.axes[host].legend(handles=handles, title=title, **settings)
    return dict(role="parameters", placement="in_axes", host_axes_index=host, location=location,
                title=title, labels=labels, scope="all parameter curves", series_keys=list(raw.PARAMETERS[mode]),
                legend_index=0, anchor_axes=anchor, selection_reason="preserve the requested legacy panel layout")


def render_aligned_single(mode, panel, observable, data, profile):
    figure, specs, placements = raw.render_single(mode, panel, observable, data, profile)
    key = (mode, panel, observable)
    if key not in SPECIAL_SINGLES:
        return figure, specs, placements
    axis, spec = figure.axes[0], specs[0]
    if key == ("mode_b", "T160.0", "eta_over_s"):
        axis.set_yscale("linear")
        axis.relim()
        axis.autoscale_view(scalex=False)
        axis.yaxis.set_major_locator(MaxNLocator(nbins=5, min_n_ticks=3, steps=[1, 2, 5, 10]))
        axis.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
        spec.update(axis_scale="linear", axis_scale_reason="preserve_v13_single_panel_axis_type")
        configure_axis_ticks(axis, profile)
        raw.capture_ticks(axis, spec)
        return figure, specs, placements
    for legend in visible_legends(figure):
        legend.remove()
    location = SPECIAL_SINGLES[key]
    placements = [parameter_legend(figure, specs, profile, host=0, location=location)]
    # Keep the old legend corner; only expand the numerical viewport if needed.
    base_low, base_high = axis.get_ylim()
    for margin in (0.0, 0.15, 0.30, 0.50, 0.75):
        span = base_high - base_low
        if location.startswith("lower"):
            axis.set_ylim(max(0.0, base_low - margin * span), base_high)
        else:
            axis.set_ylim(base_low, base_high + margin * span)
        if raw.clear_layout(measure_figure(figure, intended_width_inches=6.75)):
            break
    else:
        raise ValueError(f"legacy legend corner remains obstructed: {key}")
    raw.capture_ticks(axis, spec)
    spec["viewport_reason"] = "data-derived extent expanded only to retain the legacy legend corner"
    return figure, specs, placements


def render_fixed_T(data, profile):
    with raw.figure_context(profile) as context:
        figure, axes = plt.subplots(3, 3, figsize=(6.75, 7.4), sharex=True)
        context.bind(figure)
        figure.subplots_adjust(left=0.12, right=0.98, bottom=0.09, top=0.95, hspace=0.26, wspace=0.38)
        specs = []
        for row, observable in enumerate(("eta_over_s", "zeta_over_s", "sigma_over_T")):
            for col, panel in enumerate(raw.PANELS["mode_b"]):
                axis = axes[row, col]
                spec = raw.render_panel(axis, "mode_b", panel, observable, data, profile, composite=True)
                scale = "log" if row < 2 else "linear"
                axis.set_yscale(scale)
                axis.relim()
                axis.autoscale_view(scalex=False)
                axis.margins(y=0.05)
                if scale == "log":
                    low, high = axis.get_ylim()
                    subs = (1.0,) if high / low > 10 ** 1.5 else (1.0, 2.0, 5.0)
                    axis.yaxis.set_major_locator(LogLocator(base=10, subs=subs))
                    axis.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
                    axis.yaxis.set_minor_formatter(NullFormatter())
                    visible = [value for value in axis.get_yticks() if low <= value <= high]
                    if len(visible) > 5 or len(visible) < 2:
                        # LogLocator may fall back to dense linear ticks within
                        # a narrow decade. Keep the logarithmic axis and space
                        # its ordinary-number labels at the existing font size.
                        axis.yaxis.set_major_locator(FixedLocator(raw.readable_log_ticks(low, high)))
                else:
                    axis.yaxis.set_major_locator(AutoLocator())
                    axis.yaxis.set_major_formatter(ScalarFormatter())
                configure_axis_ticks(axis, profile)
                if col:
                    axis.set_ylabel("")
                if row == 2:
                    axis.set_xlabel(r"$\xi$")
                label = f"({chr(97 + row * 3 + col)})"
                if row == 0:
                    label += rf"  $T={float(panel[1:]):.0f}$ MeV"
                axis.text(0, 1.03, label, transform=axis.transAxes, ha="left", va="bottom", fontsize=13)
                for collection in axis.collections:
                    collection.set_sizes([5.5 ** 2])
                    collection.set_linewidths([1.0])
                spec.update(axis_scale=scale, axis_scale_reason="paper_FIG4_fixed_T_row_scale")
                raw.capture_ticks(axis, spec)
                specs.append(spec)
        placements = [parameter_legend(figure, specs, profile, host=6, location="lower left",
                                      anchor=(0.4466411003595683, 0.31293851991526406), size=13)]
        endpoints = [Line2D([], [], ls="None", marker=marker, mfc="white", mec="black", ms=6,
                            label=label) for marker, label in (("o", "restored"), ("s", "broken"))]
        axes[0, 0].legend(handles=endpoints, title="First-order", loc="lower left",
                          bbox_to_anchor=(0.16557086802564835, 0.3187792443606394), borderaxespad=0,
                          borderpad=0.25, labelspacing=0.20, handlelength=1.6, handletextpad=0.4,
                          frameon=False, fontsize=13, title_fontsize=13)
        placements.append(dict(role="endpoints", placement="in_axes", host_axes_index=0,
                               location="lower left", anchor_axes=[0.16557086802564835, 0.3187792443606394],
                               title="First-order", labels=["restored", "broken"],
                               scope="T=120 MeV, muB=900 MeV sampled first-order branch pair only",
                               applies_to=[dict(mode_key="mode_b", plot_panel="T120.0", plot_series="muB900.0")],
                               legend_index=0, selection_reason="preserve the requested legacy panel layout"))
        # visible_legends orders by host axes, not by creation order.
        for item in placements:
            item["legend_index"] = 0 if item["host_axes_index"] == 0 else 1
        placements.sort(key=lambda item: item["legend_index"])
        quality = measure_figure(figure, intended_width_inches=6.75)
        if not raw.clear_layout(quality):
            raise ValueError("FIG4 legacy layout obstructs raw data: " + json.dumps(
                {key: quality[key] for key in raw.LAYOUT_ERRORS if quality[key]}))
        return figure, specs, placements


def chart_requests():
    yield from raw.chart_requests()
    yield FIG4, "composite", ("figure4_transport_fixed_T",)


def render_request(kind, arguments, data, profile):
    if arguments == ("figure4_transport_fixed_T",):
        return render_fixed_T(data, profile)
    if kind == "single":
        return render_aligned_single(*arguments, data, profile)
    return raw.render_request(kind, arguments, data, profile)


def render_signature(figure, specs, placements):
    payload = dict(panels=specs, legends=placements, axes=[dict(position=list(axis.get_position().bounds),
        xlabel=axis.get_xlabel(), ylabel=axis.get_ylabel(), title=axis.get_title(),
        texts=[text.get_text() for text in axis.texts]) for axis in figure.axes])
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return dict(schema="transport_raw_render_signature_v1", sha256=digest,
                basis="raw source-row identities, panel scales/ranges/ticks, endpoints, labels, axes geometry and legends")


def preview(destination):
    destination.mkdir(parents=True, exist_ok=False)
    data, profile = raw.load_inputs(), publication_profile()
    requests = [(FIG4, "composite", ("figure4_transport_fixed_T",))]
    requests.extend((Path("single") / args[0] / args[1] / args[2], "single", args) for args in SPECIAL_SINGLES)
    for relative, kind, arguments in requests:
        figure, specs, placements = render_request(kind, arguments, data, profile)
        try:
            vertices = raw.verify_raw_vertices(figure, specs, data)
            figure.set_dpi(profile.dpi)
            _, quality = export_figure(figure, destination / relative, profile, formats=["png"])
            print(json.dumps(dict(figure=str(relative), raw_vertices=vertices,
                                  minimum_glyph_mm=quality["minimum_capital_numeral_height_mm"],
                                  issues={key: quality[key] for key in raw.LAYOUT_ERRORS})), flush=True)
        finally:
            plt.close(figure)


def freeze_sources(results, profile):
    records = raw.freeze_sources(results, profile)
    paths = [(Path(__file__), "publication_delivery_renderer"),
             (ROOT / "tests/unit/python/test_phase_guided_publication_clean_v14.py", "publication_renderer_tests"),
             (LEGACY_ROOT / "publication_clean_v13_png_review/plot_manifest.json", "v13_form_reference_only"),
             (LEGACY_ROOT / "transport_fixed_T_panel_legends_20261021__png_review/plot_manifest.json", "fixed_T_form_reference_only")]
    for path, role in paths:
        target = (results / "source" / path.relative_to(ROOT) if path.is_relative_to(ROOT)
                  else results / "source/references" / role / path.name)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        records.append(input_record(target, role=role))
    return records


def write_caption(results):
    path = results / "caption_handoff.md"
    path.write_text("""# Full-RS publication_clean_v14

本包直接读取全网格验收通过的 1819 个原始点和 12 个同名列。
无平滑、插值、删点或局部替换，关闭绘制路径简化；470 个有限 τ 比值提示保留。
仅断开背景审计确认的一阶相变相邻原始点，空心圆/方块表示恢复/破缺侧最近采样点。
它们不是外推的共存端点。所有数值检查在 quote-remote，纯 Python 绘图在本机。

figure1：列为 μB=0/450/900 MeV，行为 τu/τs/τū/τs̄，纵轴统一 log。
figure2：同样三列，行为 η/s、ζ/s、σ/T，纵轴 linear。
两图实线/虚线/点划线依次为 αT=1.0/1.1/1.2，T=αT Tphase(μB,ξ=0)。
mode-A 一阶断线只适用于 μB=900 MeV、αT=1.0 的 ξ=−0.003 到 +0.003。
温度原始精确映射、每条曲线输入行和每幅图轴范围均见 plot_manifest.json。

figure4：列为 T=120/160/200 MeV，行为 η/s、ζ/s、σ/T；前两行 log，第三行 linear。
实线/虚线/点划线依次为 μB=0/450/900 MeV。图例保留用户指定的面板内布局。
一阶断线只适用于 T=120 MeV、μB=900 MeV 的 ξ=−0.13 到 −0.12。
建议 FIG4 图注：Transport ratios at fixed temperatures, shown as functions of ξ.
Columns correspond to T=120, 160 and 200 MeV. Line styles identify μB=0, 450 and
900 MeV. Open circles and squares mark the nearest sampled restored/broken points
on either side of the reviewed first-order discontinuity. All raw samples are retained.

低值细节图保留线性窗口：0.57–0.66、0.57–0.66、0.375–0.43、0.27–0.41 fm。
前两窗比历史 v13 下延 0.01 fm，以容纳本次原始数据；离开窗口不表示物理端点。
其余 72 张单图沿用 v13 的面板、轴类型、图例角落、线型与蓝/红/绿色。
坐标范围随原始值适配。FIG4 尺寸 6.75×7.4 in，其他图保留 v13 尺寸。

τ 单位 fm；η/ζ 单位 fm⁻³；σ 单位 fm⁻¹；η/s、ζ/s、σ/T 与 ξ 无量纲。
μB=3μq；hadron/quark 标签指当前 PNJL 模型的手征破缺/恢复分支。
图的最小插入宽度由各自 placement_limits 指定，不能默认缩成单栏。
PNG 接受、数值晋升及 PDF 验收记录保存在同名 results case；论文文件并未自动替换。
""", encoding="utf-8")
    return path


def build_review(case=REVIEW_CASE):
    figures, results = raw.FIGURE_PARENT / case, raw.RESULT_PARENT / case
    if Path(case).name != case or any(path.exists() for path in (figures, results)):
        raise FileExistsError("review must use a new single-name case")
    data, profile = raw.load_inputs(), publication_profile()
    figures.mkdir(parents=True)
    results.mkdir(parents=True)
    inputs = [input_record(path, role="accepted_raw_grid" if path in data.sources.values() else "numerical_acceptance_evidence")
              for path in data.evidence_paths]
    inputs.extend(freeze_sources(results, profile))
    inputs.append(input_record(write_caption(results), role="raw_figure_caption_and_parameter_mapping"))
    generator = generator_record(Path(__file__), command=f"python {Path(__file__).relative_to(ROOT).as_posix()} --png-review --case {case}",
        source_snapshot=results / "source" / Path(__file__).relative_to(ROOT),
        runtime=runtime_record({"matplotlib": matplotlib.__version__, "Pillow": Image.__version__}))
    charts, cache = [], {}
    for relative, kind, arguments in chart_requests():
        figure, specs, placements = render_request(kind, arguments, data, profile)
        signature = render_signature(figure, specs, placements)
        chart = raw.write_chart(figure, relative, specs, placements, kind, data, profile, figures, results, inputs, generator, cache)
        record = chart["manifest_record"]
        record.update(derived_display_geometry=False, caption_parameters=dict(mode=specs[0]["mode_key"],
            parameter_encoding=list(raw.PARAMETERS[specs[0]["mode_key"]])), frozen_render_signature=signature)
        record["rendering"].update(color_route="color_print_and_online", coexistence_fill="none",
            parameter_encoding=[style.record() for style in raw.make_registry(specs[0]["mode_key"], profile).styles],
            palette_override=dict(colors=list(LEGACY_COLORS), reason="retain established manuscript palette"))
        chart.update(relative_stem=relative.as_posix(), arguments=list(arguments))
        if relative == FIG4:
            record["rendering"]["size_override_reason"] = "user-requested 20261021 fixed-T nine-panel layout"
        accepted_record = results / "accepted_manifests" / relative.with_suffix(".json")
        write_manifest(accepted_record, record)
        chart["accepted_record"] = raw.display_path(accepted_record)
        charts.append(chart)
        print(f"[{len(charts):02d}/76] {relative.as_posix()}", flush=True)
    package = dict(case_slug=case, single_figure_count=72, composite_figure_count=3, detail_figure_count=1,
                   color_figure_count=76, grayscale_figure_count=76, raw_input_point_count=1819,
                   numerical_data_postprocessing="none", delivery_stage="png_review", manuscript_eligible=False,
                   current_publication_layer=False, vector_delivery_pending=True, charts=charts)
    path = figures / "plot_manifest.json"
    write_manifest(path, build_bundle(package, [chart["manifest_record"] for chart in charts]))
    errors = validate_manifest(path)
    if errors:
        raise ValueError("; ".join(errors))
    write_manifest(results / "quality_review.json", dict(status="passed", figures=76, png_count=152,
        original_vertices_verified=sum(chart["manifest_record"]["validation"]["original_vertex_count"] for chart in charts),
        no_postprocessing=True, validation_errors=errors, manifest=input_record(path, role="validated_review_bundle")))


def export_delivery(review_case, case, authorization_path):
    figures, results = raw.FIGURE_PARENT / case, raw.RESULT_PARENT / case
    if Path(case).name != case or any(path.exists() for path in (figures, results)):
        raise FileExistsError("delivery must use a new single-name case")
    source_manifest = raw.FIGURE_PARENT / review_case / "plot_manifest.json"
    authorization = raw.read_json(authorization_path)
    if (authorization.get("schema") != "transport_figure_conditional_author_acceptance_v1"
            or authorization.get("condition_verified") is not True
            or authorization.get("pdf_export_authorized") is not True
            or not authorization.get("author_instruction")):
        raise ValueError("actual author instruction and verified acceptance condition required")
    if authorization["review_manifest_sha256"] != sha256_file(source_manifest):
        raise ValueError("authorization binds a different PNG bundle")
    errors = validate_manifest(source_manifest)
    if errors:
        raise ValueError("; ".join(errors))
    pairs = expand_bundle(raw.read_json(source_manifest))
    if any(record["generator"]["sha256"] != sha256_file(Path(__file__)) for _, record in pairs):
        raise ValueError("renderer changed after PNG acceptance")
    data, profile = raw.load_inputs(), publication_profile()
    figures.mkdir(parents=True)
    results.mkdir(parents=True)
    shutil.copyfile(authorization_path, results / "author_authorization.json")
    charts, cache = [], {}
    for source_chart, accepted in pairs:
        relative = Path(source_chart["relative_stem"])
        arguments = tuple(source_chart["arguments"])
        figure, specs, placements = render_request(source_chart["kind"], arguments, data, profile)
        try:
            if render_signature(figure, specs, placements) != accepted["frozen_render_signature"]:
                raise ValueError(f"rendered geometry/semantics changed: {relative}")
            raw.verify_raw_vertices(figure, specs, data)
            figure.set_dpi(profile.dpi)
            buffer = io.BytesIO()
            with figure_style_context(figure):
                figure.savefig(buffer, format="png", dpi=profile.dpi, bbox_inches=None)
            if hashlib.sha256(buffer.getvalue()).hexdigest() != accepted["outputs"][0]["sha256"]:
                raise ValueError(f"re-rendered PNG is not byte-identical: {relative}")
            pdfs, quality = export_figure(figure, figures / relative, profile, formats=["pdf"])
        finally:
            plt.close(figure)
        outputs = []
        for old in accepted["outputs"]:
            source = ROOT / old["path"]
            target = figures / source.relative_to(raw.FIGURE_PARENT / review_case)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
            item = deepcopy(old)
            item.update(path=raw.display_path(target), inspection=inspect_export(target))
            outputs.append(item)
        outputs.extend(pdfs)
        acceptance_path = results / "png_acceptance" / relative.with_suffix(".json")
        acceptance = dict(schema="plot_png_acceptance_v1", status="author_accepted", pdf_export_authorized=True,
            author_instruction=authorization["author_instruction"], acceptance_basis="explicit conditional author instruction",
            author_personally_viewed_this_new_png=False, condition_verification=authorization["condition_verification"],
            accepted_manifest=input_record(ROOT / source_chart["accepted_record"], role="accepted_png_manifest"),
            accepted_png=accepted["outputs"][0])
        write_manifest(acceptance_path, acceptance)
        record = deepcopy(accepted)
        record.update(asset_id=f"relaxtime.phase_guided.{case}.{relative.as_posix()}", case_slug=case,
            figure_mode="strict", semantic_status="formal_full_rs_raw_transport", publication_scope="main_text_candidate",
            outputs=outputs, manuscript_eligible=True, current_publication_layer=True, formal_production_promoted=True,
            delivery_stage="vector_delivery", vector_delivery_pending=False,
            author_review=dict(status="author_accepted_png", record=input_record(acceptance_path, role="author_png_acceptance")))
        record["validation"].update(strict_gate=True, png_rerender_byte_identical=True)
        record["rendering"].update(quality=quality, delivery_stage="vector_delivery", vector_delivery_pending=False,
            output_formats=["png", "pdf"], placement_limits=placement_limits(quality, profile, outputs))
        record["rendering"]["grayscale_review"].update(status="agent_checked_under_conditional_author_acceptance",
            output=outputs[1]["path"])
        errors = validate_manifest_record(record, hash_cache=cache)
        if errors:
            raise ValueError(f"{relative}: " + "; ".join(errors))
        charts.append(dict(stem=raw.display_path(figures / relative), kind=source_chart["kind"], mode_key=source_chart["mode_key"],
                           relative_stem=relative.as_posix(), outputs=outputs, manifest_record=record))
        print(f"[{len(charts):02d}/76] PDF verified: {relative.as_posix()}", flush=True)
    package = dict(case_slug=case, single_figure_count=72, composite_figure_count=3, detail_figure_count=1,
                   color_figure_count=76, grayscale_figure_count=76, vector_figure_count=76, raw_input_point_count=1819,
                   numerical_data_postprocessing="none", manuscript_eligible=True, current_publication_layer=True,
                   vector_delivery_pending=False, delivery_stage="vector_delivery", charts=charts)
    path = figures / "plot_manifest.json"
    write_manifest(path, build_bundle(package, [chart["manifest_record"] for chart in charts]))
    errors = validate_manifest(path)
    if errors:
        raise ValueError("; ".join(errors))
    write_manifest(results / "delivery_validation.json", dict(status="passed", figures=76, png_count=152, pdf_count=76,
        all_png_rerenders_byte_identical=True, all_pdfs_single_page_embedded_non_type3_fonts_and_vector=True,
        no_postprocessing=True, validation_errors=errors, manifest=input_record(path, role="validated_delivery_bundle")))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    stage = parser.add_mutually_exclusive_group(required=True)
    stage.add_argument("--preview", type=Path)
    stage.add_argument("--png-review", action="store_true")
    stage.add_argument("--vector-delivery", action="store_true")
    parser.add_argument("--case")
    parser.add_argument("--review-case", default=REVIEW_CASE)
    parser.add_argument("--authorization", type=Path)
    args = parser.parse_args()
    if args.preview:
        preview(args.preview)
    elif args.png_review:
        build_review(args.case or REVIEW_CASE)
    elif args.authorization:
        export_delivery(args.review_case, args.case or FINAL_CASE, args.authorization)
    else:
        parser.error("--vector-delivery requires --authorization")


if __name__ == "__main__":
    main()
