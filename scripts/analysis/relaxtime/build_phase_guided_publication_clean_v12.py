#!/usr/bin/env python3
"""Render final-size v12 review figures from the unchanged v5/v11 display data."""

from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.analysis.relaxtime import build_phase_guided_publication_clean_v11 as V11
from scripts.plotting.plot_bundle import build_bundle
from scripts.plotting.plot_manifest import (
    build_manifest, generator_record, input_record, output_record, runtime_record,
    sha256_file, write_manifest,
)
from scripts.plotting.plot_quality import export_figure, inspect_export, measure_figure, placement_limits
from scripts.plotting.plot_style import configure_matplotlib, load_profile
from scripts.plotting.validate_plot_artifact import validate_manifest_record

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MaxNLocator, FormatStrFormatter
from PIL import Image, ImageOps

ANALYSIS_ROOT = V11.V10.TRANSPORT_ANALYSIS_ROOT / "phase_guided_transport_publication_clean_v12_png_review"
FIGURE_ROOT = ROOT / "data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v12_png_review"
WIDTH_IN = 6.75
LABEL_PT = 13.0
TICK_PT = 11.0
LEGEND_PT = 13.0
COMPOSITE_SIZES = {
    "figure1_relaxation_times_comparison": (WIDTH_IN, 7.1),
    "figure2_transport_coefficients_comparison": (WIDTH_IN, 5.9),
}
COMPOSITE_LAYOUT = dict(left=0.115, right=0.98, bottom=0.085, top=0.875, hspace=0.20, wspace=0.31)
DETAILS = (
    ("muB0.0", "tau_u", "(a)", (0.58, 0.66)),
    ("muB0.0", "tau_ubar", "(g)", (0.58, 0.66)),
    ("muB450.0", "tau_ubar", "(h)", (0.375, 0.43)),
    ("muB900.0", "tau_ubar", "(i)", (0.27, 0.41)),
)


def common_legends(figure, specs, profile, *, endpoints=True, top_inches=0):
    handles = V11.legend_handles(specs, profile, endpoints=endpoints)
    height = figure.get_size_inches()[1]
    figure.legend(handles=handles[:3], loc="upper center", bbox_to_anchor=(0.55, 0.999 - top_inches / height),
                  ncol=3, frameon=False, fontsize=LEGEND_PT,
                  handlelength=1.65, columnspacing=1.0, handletextpad=0.4, borderpad=0.1)
    if len(handles) > 3:
        figure.legend(handles=handles[3:], loc="upper center", bbox_to_anchor=(0.55, 1 - (0.30 + top_inches) / height),
                      ncol=2, frameon=False, fontsize=LEGEND_PT,
                      handlelength=1.0, columnspacing=1.0, handletextpad=0.35, borderpad=0.1)


def render_single(mode, panel, observable, grouped, gap_map, profile):
    figure, axis = plt.subplots(figsize=(WIDTH_IN, 4.6))
    spec = V11.render_panel(axis, mode, panel, observable, grouped, gap_map, profile)
    endpoints = bool(spec["endpoints"])
    figure.subplots_adjust(left=0.16, right=0.97, bottom=0.14,
                           top=1 - (1.0 if endpoints else 0.66) / 4.6)
    axis.set_xlabel(r"$\xi$")
    figure.suptitle(V11.V10.panel_title(mode, panel), y=0.99)
    common_legends(figure, [spec], profile, endpoints=endpoints, top_inches=0.26)
    return figure, [spec]


def render_composite(name, grouped, gap_map, profile):
    observables = V11.COMPOSITES[name]
    with matplotlib.rc_context({"font.size": LABEL_PT, "axes.labelsize": LABEL_PT,
                                "axes.titlesize": LABEL_PT, "xtick.labelsize": TICK_PT,
                                "ytick.labelsize": TICK_PT}):
        figure, axes = plt.subplots(len(observables), 3, figsize=COMPOSITE_SIZES[name], sharex=True)
        figure.subplots_adjust(**{**COMPOSITE_LAYOUT, "top": 1 - 0.90 / figure.get_size_inches()[1]})
        specs = []
        for row, observable in enumerate(observables):
            row_specs = []
            for col, panel in enumerate(V11.PANELS):
                axis = axes[row, col]
                spec = V11.render_panel(axis, "mode_a", panel, observable, grouped, gap_map, profile, composite=True)
                row_specs.append(spec)
                if col:
                    axis.set_ylabel("")
                if row == 0:
                    axis.set_title(V11.V10.panel_title("mode_a", panel), pad=5)
                if row == len(observables) - 1:
                    axis.set_xlabel(r"$\xi$")
                axis.text(0, 1.02, f"({chr(97 + row * 3 + col)})", transform=axis.transAxes,
                          ha="left", va="bottom", fontsize=11)
            V11.apply_row_precision(axes[row], row_specs, observable)
            specs.extend(row_specs)
        common_legends(figure, specs, profile)
        return figure, specs


def render_detail(grouped, gap_map, profile):
    with matplotlib.rc_context({"axes.labelsize": LABEL_PT, "axes.titlesize": LABEL_PT,
                                "xtick.labelsize": TICK_PT, "ytick.labelsize": TICK_PT}):
        figure, axes = plt.subplots(2, 2, figsize=(WIDTH_IN, 5.5), sharex=True)
        figure.subplots_adjust(left=0.12, right=0.97, bottom=0.12, top=0.86, hspace=0.28, wspace=0.32)
        specs = []
        for axis, (panel, observable, source_panel, limits) in zip(axes.flat, DETAILS):
            spec = V11.render_panel(axis, "mode_a", panel, observable, grouped, gap_map, profile, composite=True)
            # This is a disclosed viewport on the original curves, not new data.
            axis.set_yscale("linear")
            axis.set_ylim(*limits)
            axis.yaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 5, 10]))
            axis.yaxis.set_minor_locator(AutoMinorLocator(2))
            axis.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
            axis.set_title(V11.V10.panel_title("mode_a", panel), pad=5)
            axis.text(0, 1.02, source_panel, transform=axis.transAxes, va="bottom", fontsize=11)
            spec.update(axis_scale="linear", axis_scale_reason="explicit_low_value_detail_view",
                        source_figure="figure1_relaxation_times_comparison", source_panel=source_panel,
                        viewport_clipping="values outside the declared y view are intentionally not visible; all curve vertices retained",
                        y_tick_formatter_policy="fixed_decimal_from_detail_major_tick_spacing", y_tick_decimal_places=2)
            V11.capture_tick_spec(axis, spec)
            specs.append(spec)
        for axis in axes[-1]:
            axis.set_xlabel(r"$\xi$")
        common_legends(figure, specs, profile, endpoints=False)
        return figure, specs


def preview(output_root):
    output_root.mkdir(parents=True, exist_ok=False)
    points, _, _, _, gaps, _ = V11.PARENT.load_v5_inputs()
    grouped, gap_map = V11.V10.group_inputs(points, gaps)
    profile = load_profile("candidate_aps_v2")
    configure_matplotlib(profile)
    for name in [*COMPOSITE_SIZES, "figure1_low_value_details"]:
        figure, specs = (render_detail(grouped, gap_map, profile) if name == "figure1_low_value_details"
                         else render_composite(name, grouped, gap_map, profile))
        quality = measure_figure(figure, intended_width_inches=WIDTH_IN)
        figure.savefig(output_root / f"{name}.png", dpi=160)
        plt.close(figure)
        print(json.dumps({"figure": name, **{k: quality[k] for k in (
            "minimum_capital_numeral_height_mm", "clipped_text", "text_overlap_pairs",
            "legend_axes_overlap_count", "legend_curve_overlap_count", "legend_landmark_overlap_count")}}))


def retained_v11_paths():
    roots = [V11.FIGURE_ROOT, V11.ANALYSIS_ROOT,
             V11.FIGURE_ROOT.with_name("publication_clean_v11_pdf_review"),
             V11.ANALYSIS_ROOT.with_name("phase_guided_transport_publication_clean_v11_pdf_review")]
    return sorted(path for root in roots for path in root.rglob("*") if path.is_file())


def parent_records(package, profile):
    records = V11.parent_records(package, profile)
    records.extend(input_record(path, role="retained_v11_artifact") for path in retained_v11_paths())
    code = [*Path(__file__).parent.glob("build_phase_guided_publication_clean_v*.py"),
            ROOT / "scripts/plotting/plot_quality.py", ROOT / "scripts/plotting/plot_manifest.py",
            ROOT / "scripts/plotting/plot_style.py", ROOT / "scripts/plotting/plot_bundle.py",
            ROOT / "scripts/plotting/validate_plot_artifact.py",
            ROOT / "docs/analysis/relaxtime/phase_guided_transport/plotting_case_contract.md"]
    records.extend(input_record(path, role="v12_renderer_or_contract") for path in code)
    return list({record["path"]: record for record in records}.values())


def grayscale_record(color, target, profile):
    if target.exists():
        raise FileExistsError(f"refusing to overwrite grayscale review: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    with Image.open(ROOT / color["path"]) as image:
        ImageOps.grayscale(image.convert("RGB")).save(target, dpi=(profile.dpi, profile.dpi))
    record = output_record(target, fmt="png", dpi=profile.dpi, vector=False)
    record.update(inspection=inspect_export(target), role="grayscale_review",
                  source_color_sha256=color["sha256"], conversion="Pillow RGB to L; 0.299 R + 0.587 G + 0.114 B")
    return record


def write_chart(figure, stem, specs, profile, inputs, generator, kind, hash_cache):
    try:
        outputs, quality = export_figure(figure, stem, profile, formats=("png",))
    finally:
        plt.close(figure)
    outputs[0]["role"] = "color_review"
    outputs.append(grayscale_record(outputs[0], FIGURE_ROOT / "grayscale_review" / stem.relative_to(FIGURE_ROOT).with_suffix(".png"), profile))
    axes = []
    for spec in specs:
        axes.extend([
            {"field": "xi", "source_unit": "dimensionless", "display_unit": "dimensionless",
             "label": r"$\xi$", "transform": "identity", "panel": spec["plot_panel"]},
            {"field": spec["observable"], "source_unit": V11.UNITS[spec["observable"]],
             "display_unit": V11.UNITS[spec["observable"]], "label": V11.LABELS[spec["observable"]],
             "transform": spec["axis_scale"], "panel": spec["plot_panel"]},
        ])
    limits = placement_limits(quality, profile, outputs)
    manifest = build_manifest(
        asset_id=f"relaxtime.phase_guided.v12.png_review.{stem.relative_to(FIGURE_ROOT).as_posix()}",
        figure_family="phase_guided_transport", case_slug=stem.name, figure_mode="audit",
        semantic_status="author_review_display_derivative", style_profile=profile.profile_id,
        publication_scope="internal_review", generator=generator, inputs=inputs, axes=axes,
        series=[item for spec in specs for item in spec["series"]], outputs=outputs,
        selection_rule="all frozen v5 clean_value vertices in original xi order; v11 phase gaps preserved",
        interpolation_policy="inherited v5 display adjustments; no new interpolation or replacement",
        connector_policy="forbidden", missing_value_policy="preserve first-order and missing-support gaps",
        validation={"finite": True, "duplicate_keys": True, "support": True, "strict_gate": False},
        rendering={
            "column": "double_column", "figure_size_inches": quality["figure_size_inches"],
            "size_override_reason": "retained v11 main canvas with external legends" if kind == "composite" else
                                    "explicit four-panel low-value detail" if kind == "detail" else None,
            "typography_exception": None, "quality": quality, "placement_limits": limits,
            "font_overrides_pt": {"labels_titles": LABEL_PT, "ticks": TICK_PT, "legends": LEGEND_PT,
                                  "panel_labels": 11} if kind != "single" else None,
            "legend_policy": "shared_external", "legend_outside": True, "legend_alignment": "center",
            "legend_scope": {"parameters": "all panels in this figure", "endpoints":
                "muB900.0 alpha1.0 first-order curves only" if kind == "composite" else
                "endpoints outside detail view; no endpoint key" if kind == "detail" else
                "actual endpoints in panel_specs only"},
            "panel_specs": specs, "color_route": "undecided_review",
            "grayscale_review": {"status": "author_review_required", "output": outputs[1]["path"],
                                  "machine_checks": "same geometry and pixel dimensions; grayscale conversion only",
                                  "visual_checks": ["trace each line style", "identify branch endpoints",
                                                    "inspect low-value detail", "read every label at insertion width"]},
            "caption_handoff": (ANALYSIS_ROOT / "caption_handoff.md").relative_to(ROOT).as_posix(),
            "delivery_stage": "png_review", "vector_delivery_pending": True, "output_formats": ["png"],
        }, calculation_sha=V11.PARENT.V5.V4.V3.V2.V1.CALCULATION_SHA,
    )
    manifest.update(manuscript_eligible=False, current_publication_layer=False, solver_called=False,
                    canonical_data_modified=False, new_display_values=False, delivery_stage="png_review",
                    vector_delivery_pending=True, numerical_status="inherited_author_accepted_display_only",
                    raw_manuscript_eligible=False, workflow_head_sha=V11.PARENT.V5.V4.V3.V2.V1.WORKFLOW_HEAD_SHA)
    violations = validate_manifest_record(manifest, hash_cache=hash_cache)
    # In-axes legends on singles must also pass the explicit visibility checks.
    if quality["legend_curve_overlap_count"] or quality["legend_landmark_overlap_count"]:
        violations.append("v12 legend intersects a curve or landmark")
    if violations:
        raise ValueError(f"{stem}: " + "; ".join(violations))
    return {"stem": stem.relative_to(ROOT).as_posix(), "kind": kind, "mode_key": specs[0]["mode_key"],
            "outputs": outputs, "manifest_record": manifest}


def write_handoff(points):
    temperatures = {(r["plot_panel"], r["plot_series"]): float(r["T_MeV"]) for r in points if r["mode_key"] == "mode_a"}
    mapping = "\n".join("| " + panel + " | " + " | ".join(
        f"{temperatures[(panel, f'alpha{alpha:.1f}')]:.6f}" for alpha in (1, 1.1, 1.2)) + " |" for panel in V11.PANELS)
    text = r"""# v12 图注与插入尺寸交接

两张主图均在 171.45 mm（6.75 in）宽度直接绘制。不要默认缩为单栏；
逐图 `rendering.placement_limits` 给出由字形、线宽、端点和 PNG 有效 DPI
共同约束的插入区间。推荐保持原生宽度，最终整篇稿件另行测量。

Figure 1（弛豫时间）建议图注：

> Relaxation times as functions of the anisotropy parameter $\xi$. Columns
> correspond to $\mu_B=0$, 450, and 900 MeV; rows show $\tau_u$, $\tau_s$,
> $\tau_{\bar u}$, and $\tau_{\bar s}$ in fm. All vertical axes are logarithmic,
> with independent ranges. The common color and line-style key above the
> panels applies throughout: solid, dashed and dash-dotted lines correspond
> to $\alpha_T=1.0$, 1.1 and 1.2. For $\mu_B=900$ MeV and $\alpha_T=1.0$,
> open circles and squares mark the endpoints of the chirally restored and
> chirally broken branches at the first-order transition, respectively.
> Disconnected segments retain the branch gap. Low-value details of panels
> (a), (g), (h) and (i) are provided separately with linear vertical axes.

Figure 2（输运系数）建议图注：

> Transport ratios as functions of $\xi$. Columns correspond to $\mu_B=0$,
> 450 and 900 MeV; rows show $\eta/s$, $\zeta/s$ and $\sigma/T$. Vertical axes
> are linear with independent ranges. The common line-style key and the
> first-order endpoint convention are the same as in Figure 1. The endpoint
> key applies only to $\mu_B=900$ MeV, $\alpha_T=1.0$ curves.

局部视图建议图注：

> Low-value views of Figure 1 panels (a), (g), (h) and (i), using the same
> display values and line styles. The vertical axes are linear and restricted
> to the indicated ranges; curves leaving the view continue in Figure 1.
> First-order endpoints are outside these windows. Panel labels refer to the
> corresponding main-figure panels; these views introduce no new data.

注意：这里的 chirally broken quark branch 不表示另行计算了强子输运。
显示曲线继承 v5 已披露的 display adjustments；没有新增平滑或数值收敛证据。
论文讨论不得把不同纵轴范围下的视觉斜率当成可直接比较的增长率。

精确冻结温度映射（MeV；表格精度不代表不确定度）：

| panel | alpha_T=1.0 | alpha_T=1.1 | alpha_T=1.2 |
| --- | ---: | ---: | ---: |
""" + mapping + "\n\nPNG 审查层：作者接受待定，PDF 待后续交付，manuscript_eligible=false，current=v5。\n"
    (ANALYSIS_ROOT / "caption_handoff.md").write_text(text, encoding="utf-8", newline="\n")


def build_review():
    if FIGURE_ROOT.exists() or ANALYSIS_ROOT.exists():
        raise FileExistsError("refusing to overwrite an existing v12 case")
    points, package, _, _, gaps, table_hashes = V11.PARENT.load_v5_inputs()
    keys = [(row["mode_key"], row["plot_panel"], row["plot_series"], row["observable"], row["xi"]) for row in points]
    if len(keys) != len(set(keys)):
        raise ValueError("duplicate frozen display coordinates")
    for record in package["outputs"]:
        if sha256_file(ROOT / record["path"]) != record["sha256"]:
            raise ValueError(f"v5 input hash mismatch: {record['path']}")
    profile = load_profile("candidate_aps_v2")
    font = configure_matplotlib(profile)
    grouped, gap_map = V11.V10.group_inputs(points, gaps)
    inputs = parent_records(package, profile)
    generator = generator_record(Path(__file__), command="python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v12.py --png-review",
                                 runtime=runtime_record({"matplotlib": matplotlib.__version__, "font": font}))
    ANALYSIS_ROOT.mkdir(parents=True, exist_ok=False)
    charts, hash_cache = [], {}
    for mode in ("mode_a", "mode_b"):
        for panel in sorted({key[1] for key in grouped if key[0] == mode}):
            for observable in V11.OBSERVABLES:
                figure, specs = render_single(mode, panel, observable, grouped, gap_map, profile)
                mode_dir = "mode_a_fixed_muB_phase_scaled" if mode == "mode_a" else "mode_b_fixed_T_sparse_muB"
                stem = FIGURE_ROOT / mode_dir / f"plot_panel={panel}" / f"{observable}_vs_xi"
                charts.append(write_chart(figure, stem, specs, profile, inputs, generator, "single", hash_cache))
        print(f"[v12] {mode}: 36 color/gray chart pairs validated", flush=True)
    for name in COMPOSITE_SIZES:
        figure, specs = render_composite(name, grouped, gap_map, profile)
        charts.append(write_chart(figure, FIGURE_ROOT / "composites" / name, specs, profile, inputs, generator, "composite", hash_cache))
    figure, specs = render_detail(grouped, gap_map, profile)
    charts.append(write_chart(figure, FIGURE_ROOT / "details/figure1_low_value_details", specs, profile, inputs, generator, "detail", hash_cache))
    for record in inputs:
        if sha256_file(ROOT / record["path"]) != record["sha256"]:
            raise ValueError(f"protected input changed: {record['path']}")
    write_handoff(points)
    index_path = FIGURE_ROOT / "plot_manifest.json"
    index = {
        "schema": "publication_clean_v12_png_review_figure_index_v1", "status": "author_review_required",
        "single_figure_count": 72, "composite_figure_count": 2, "detail_figure_count": 1,
        "color_png_count": 75, "grayscale_png_count": 75, "mode_counts": {"mode_a": 36, "mode_b": 36},
        "manuscript_eligible": False, "current_publication_layer": False, "delivery_stage": "png_review",
        "vector_delivery_pending": True, "charts": charts,
    }
    write_manifest(index_path, build_bundle(index, [chart["manifest_record"] for chart in charts]))
    report = {
        "schema": "publication_clean_v12_placement_report_v1", "recommended_width_mm": WIDTH_IN * 25.4,
        "source": "native-size renderer measurements; not a compiled manuscript or physical-printer test",
        "figures": [{"stem": chart["stem"], "kind": chart["kind"],
                     "physical_size_inches": chart["manifest_record"]["rendering"]["figure_size_inches"],
                     "placement_limits": chart["manifest_record"]["rendering"]["placement_limits"]} for chart in charts],
    }
    write_manifest(ANALYSIS_ROOT / "placement_report.json", report)
    (ANALYSIS_ROOT / "README.md").write_text("""# publication_clean_v12 PNG 审查层

v12 从相同冻结点表生成 72 张单图、2 张主复合图和 1 张局部视图；
每张都交付彩色及灰度 PNG，共 150 张。所有图按 171.45 mm 宽、600 dpi
导出，保持 v11 两张主图的原生画布尺寸，不使用紧裁切。

所有图改用图外公共图例。复合图的图例和带数学下标的列标题为 13 pt，刻度为
11 pt；最小大写/数字字形必须达到 2 mm，不再使用紧凑字号例外。
主图坐标语义和全部曲线顶点沿用 v11。新增局部图仅给出明确的线性纵轴
窗口，帮助辨认 Figure 1 的低值细节；它不替代完整主图，也不新增数据。

查看 `placement_report.json` 的逐图插入区间和 `caption_handoff.md` 的图注。
推荐按 171.45 mm 宽插入；单图也不能直接缩成单栏。超过推荐宽度时，PNG
有效 DPI 会降低；低于最小宽度时，小字或标记不合格。PDF 后续交付时需要
重新计算其尺寸约束，不用 PNG 分辨率判定矢量曲线质量。

`plot_manifest.json` 位于对应 figure case，集中记录所有输入、生成器、
逐图布局、尺寸限制及彩色/灰度输出 hash。灰度图片已经生成；作者仍需
在最终尺寸逐项检查线型追踪、端点、局部细节和文字。

这是独立绘图任务：manuscript_eligible=false，current=v5，PDF pending。
v5/v10/v11、原始计算和分支断线均保留；不调用求解器或数值收敛门禁。
本包没有验证当前论文整页排版，也不等同于投稿资格或数值生产晋升。

```powershell
python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v12.py --png-review
```
""", encoding="utf-8", newline="\n")
    auxiliary = [ANALYSIS_ROOT / name for name in ("README.md", "caption_handoff.md", "placement_report.json")]
    write_manifest(ANALYSIS_ROOT / "manifest.json", {
        "schema": "publication_clean_v12_png_review_package_v1", "task_classification": "independent",
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(), "status": "author_review_required",
        "author_acceptance": None, "manuscript_eligible": False, "current_publication_layer": False,
        "delivery_stage": "png_review", "vector_delivery_pending": True, "solver_called": False,
        "canonical_data_modified": False, "new_display_values": False,
        "raw_manuscript_eligible": False, "numerical_status": "inherited_author_accepted_display_only",
        "generator": generator, "inputs": inputs, "point_row_count": len(points),
        "inherited_table_hashes": table_hashes,
        "outputs": [input_record(path, role="v12_review_contract") for path in [index_path, *auxiliary]],
        "figure_index": index_path.relative_to(ROOT).as_posix(),
        "figure_index_sha256": sha256_file(index_path), "protected_inputs_verified_unchanged": True,
    })
    print(f"[v12] 75 color/gray pairs validated: {FIGURE_ROOT}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preview", type=Path, help="write a disposable preview to a new directory")
    parser.add_argument("--png-review", action="store_true", help="generate the new frozen review case")
    args = parser.parse_args()
    if args.preview and args.png_review:
        parser.error("--preview and --png-review are mutually exclusive")
    if args.preview:
        preview(args.preview)
    elif args.png_review:
        build_review()
    else:
        parser.error("choose --preview or --png-review")


if __name__ == "__main__":
    main()
