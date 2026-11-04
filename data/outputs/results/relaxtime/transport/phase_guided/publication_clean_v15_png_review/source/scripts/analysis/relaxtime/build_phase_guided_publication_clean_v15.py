#!/usr/bin/env python3
"""Render the accepted eta-augmented grid with the unchanged v14 figure layout."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg")
from PIL import Image

from scripts.analysis.relaxtime import build_phase_guided_full_rs_raw as raw
from scripts.analysis.relaxtime import build_phase_guided_publication_clean_v14 as layout
from scripts.plotting.plot_bundle import build_bundle
from scripts.plotting.plot_manifest import generator_record, input_record, runtime_record, sha256_file, write_manifest
from scripts.plotting.validate_plot_artifact import validate_manifest, validate_manifest_record

CASE = "publication_clean_v15_png_review"
PRODUCTION = raw.PRODUCTION_ROOT.parent / "full_rs_eta_augmented_20261029_v1"
REFERENCE = raw.FIGURE_PARENT / "publication_clean_v14_pdf/plot_manifest.json"
REPLAY_IDS = {"a0087", "a0088", "a0089", "a0379", "a0380", "a0381", "b0693", "b0694", "b0695"}
COUNTS = {"v14_retained": 1810, "recomputed_original": 9, "xi_inserted": 48}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def load_inputs(production=PRODUCTION):
    returned, promotion_path = production / "returned", production / "formal_promotion_v1/promotion.json"
    grid = returned / "grid"
    review_path = grid / "review.json"
    review, promotion = raw.read_json(review_path), raw.read_json(promotion_path)
    require(review["calculation_sha"] == raw.CALCULATION_SHA and review["points"] == 1867
            and review["full_grid_numerical_acceptance_passed"] is True and review["artifact_contract_passed"] is True,
            "augmented numerical grid is not accepted")
    require(promotion["formal_production_promoted"] is True
            and promotion["grid_review_sha256"] == sha256_file(review_path), "promotion does not bind this grid")
    package_path = returned / "package_manifest.json"
    package = raw.read_json(package_path)
    receipt_path = production / "local_return_verification.json"
    receipt = raw.read_json(receipt_path)
    require(receipt["status"] == "verified_return" and receipt["package_manifest_sha256"] == sha256_file(package_path),
            "returned package has not passed integrity verification")
    for item in package["files"]:
        path = returned / item["path"]
        require(path.stat().st_size == item["bytes"] and sha256_file(path) == item["sha256"],
                f"returned file changed: {item['path']}")
    evidence = [review_path, promotion_path, package_path, receipt_path]
    for name, digest in review["output_sha256"].items():
        path = grid / name
        require(sha256_file(path) == digest, f"accepted grid changed: {name}")
        evidence.append(path)
    for name in ("input_manifest.json", "selection/selection.json", "selection/points.csv", "selection/channels.csv",
                 "thermo/summary.json", "thermo/replay_checks.csv", "assembly_verification.json"):
        evidence.append(returned / name)
    evidence.extend(returned / item["path"] for item in raw.read_json(returned / "input_manifest.json")["files"])
    groups, sources, all_rows = {}, {}, []
    for mode, expected in (("mode_a", 942), ("mode_b", 925)):
        path = grid / mode / "phase_guided_transport_grid.csv"
        rows = read_csv(path)
        require(len(rows) == expected, "wrong augmented mode count")
        groups.update(raw.validate_rows(rows, mode))
        sources[mode] = path
        all_rows.extend(rows)
    require(len({r["point_id"] for r in all_rows}) == 1867, "duplicate augmented point identity")
    require(dict(Counter(r["source_kind"] for r in all_rows)) == COUNTS, "wrong source selection")
    require({r["point_id"] for r in all_rows if r["source_kind"] == "recomputed_original"} == REPLAY_IDS,
            "unexpected original-coordinate replacement")
    expected_keys = {(m, p, s) for m in raw.MODE_DIRS for p in raw.PANELS[m] for s in raw.PARAMETERS[m]}
    require(set(groups) == expected_keys, "incomplete panel/parameter coverage")
    require(sum(r["quality_flag"] == "true" for r in all_rows) == 470, "quality flags changed")
    gaps = raw.boundary_map(raw.read_json(grid / "base_background_audit.json"), groups)
    require(sum(map(len, gaps.values())) == 2, "reviewed first-order gaps changed")
    acceptance = dict(review, formal_production_promoted=True)
    return raw.RawInputs(groups, gaps, sources, list(dict.fromkeys(evidence)), acceptance)


def freeze_sources(results, profile):
    records = raw.freeze_sources(results, profile)
    for source, role in ((Path(__file__), "eta_augmented_renderer"),
                         (Path(layout.__file__), "unchanged_v14_renderer"),
                         (ROOT / "tests/unit/python/test_phase_guided_publication_clean_v15.py", "augmented_renderer_tests"),
                         (REFERENCE, "v14_form_reference_only")):
        target = results / "source" / source.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        records.append(input_record(target, role=role))
    return records


def write_caption(results):
    path = results / "caption_handoff.md"
    path.write_text(r"""# publication_clean_v15：η 阈值局部加密后的原始网格

本包直接绘制已正式验收的 1867 个数值点（mode A 942、mode B 925），读取源表的
六种 τ、η、ζ、σ、η/s、ζ/s、σ/T 共 12 个原始列，不在绘图时重算物理量。
1810 个 v14 原值原样保留；9 个原坐标采用本轮重算值，另加入 48 个新 ξ 点。
三条加密曲线是 mode A 的 μB=0/450 MeV、αT=1.0，以及 mode B 的 T=200 MeV、μB=0。
这57点含39个 s 表 rtol=10⁻⁵、18个 rtol=10⁻⁶ 的已验证结果；其余生产参数与 v14 相同。
精确参数、残差要求的实际覆盖、逐点选源和接缝审查见数值包的 manifest 与 review。
图中没有平滑、插值、拟合、删点或值混合；关闭路径简化，相邻实算点以直线相连。

新增点不伪造历史求根迭代数或相态追踪信息。48个新点的 phase 标签为 unknown；
9个原坐标保留旧相态元数据。两处已审计一阶断线及470个有限τ比值提示点保持不变。
断线只在 mode A μB=900、αT=1.0 的 ξ=−0.003/+0.003，以及 mode B T=120、μB=900
的 ξ=−0.13/−0.12。空心圆／方块为两侧最近的手征恢复／破缺采样点，不是外推共存端点。

全部76张彩色图与76张同尺寸灰度图复用 v14 renderer、面板、轴类型、单位、图例和线型：
72张单图、Figure 1/2/4 三张主组合图，以及一张四面板低值细节图。
Figure 1：列 μB=0/450/900 MeV，行 τu/τs/τū/τs̄，纵轴 log。
Figure 2：同样三列，行 η/s、ζ/s、σ/T，纵轴 linear。
这两图实线／虚线／点划线对应 αT=1.0/1.1/1.2，T=αT Tphase(μB,ξ=0)。
Figure 4：列 T=120/160/200 MeV，行 η/s、ζ/s、σ/T；前两行 log、第三行 linear。
其参数线型对应 μB=0/450/900 MeV，保留20261021版面板内图例。
低值细节窗口仍为0.57–0.66、0.57–0.66、0.375–0.43、0.27–0.41 fm。
τ单位 fm；η/ζ 为 fm⁻³；σ为 fm⁻¹；ξ、η/s、ζ/s、σ/T 无量纲；μB=3μq。
Figure 4 为6.75×7.4 in；其余尺寸沿用v14。最小插入宽度以各图 placement_limits 为准。

建议补充图注：The accepted full-RS grid includes additional computed samples around three
reverse eta-Mott thresholds. Nine original-coordinate values are replaced by the verified
recomputations at the refined integration precision; all other original values are retained.
The plotted curves connect the accepted samples without smoothing or fitted interpolation.
The open markers retain their limited meaning as the nearest sampled points across the two
reviewed first-order discontinuities.

全部数值回放与派生审计在 quote-remote 完成；本机仅进行文件核验及 Python 绘图。
本包为 PNG 审阅阶段：manuscript_eligible=false、current_publication_layer=false，
作者对新PNG的接受与后续PDF交付另行记录；数值产物的正式资格不代替视觉接受。
""", encoding="utf-8")
    return path


def build_review(case=CASE, production=PRODUCTION):
    figures, results = raw.FIGURE_PARENT / case, raw.RESULT_PARENT / case
    if not case or Path(case).name != case or any(c in case for c in ("/", "\\", ":")) or any(p.exists() for p in (figures, results)):
        raise FileExistsError("review must use a new single-name case")
    data, profile = load_inputs(production), layout.publication_profile()
    figures.mkdir(parents=True)
    results.mkdir(parents=True)
    inputs = [input_record(path, role="accepted_augmented_grid" if path in data.sources.values() else "numerical_acceptance_evidence")
              for path in data.evidence_paths]
    inputs.extend(freeze_sources(results, profile))
    inputs.append(input_record(write_caption(results), role="caption_and_parameter_mapping"))
    generator = generator_record(Path(__file__), command=f"python {Path(__file__).relative_to(ROOT).as_posix()} --png-review --case {case}",
        source_snapshot=results / "source" / Path(__file__).relative_to(ROOT),
        runtime=runtime_record({"matplotlib": matplotlib.__version__, "Pillow": Image.__version__}))
    charts, cache = [], {}
    for relative, kind, arguments in layout.chart_requests():
        figure, specs, placements = layout.render_request(kind, arguments, data, profile)
        signature = layout.render_signature(figure, specs, placements)
        chart = raw.write_chart(figure, relative, specs, placements, kind, data, profile, figures, results, inputs, generator, cache)
        record = chart["manifest_record"]
        record.update(semantic_status="eta_augmented_full_rs_grid_png_review", derived_display_geometry=False,
            selection_rule="all accepted augmented rows of each listed curve, ascending xi; 1810 unchanged, 9 recomputed and 48 added grid points",
            upstream_numerical_refinement=dict(point_count=1867, selection=COUNTS, source="accepted grid and point_provenance.csv"),
            caption_parameters=dict(mode=specs[0]["mode_key"], parameter_encoding=list(raw.PARAMETERS[specs[0]["mode_key"]])),
            frozen_render_signature=signature)
        for series in record["series"]:
            series.update(state="accepted_eta_augmented_raw", support_rule="every accepted source row of each listed curve; computed xi refinement included")
        record["rendering"].update(color_route="color_print_and_online", coexistence_fill="none",
            parameter_encoding=[style.record() for style in raw.make_registry(specs[0]["mode_key"], profile).styles],
            palette_override=dict(colors=list(layout.LEGACY_COLORS), reason="retain established manuscript palette"),
            drawing_rule="straight segments between consecutive accepted computed samples; no generated support vertices")
        if relative == layout.FIG4:
            record["rendering"]["size_override_reason"] = "user-requested 20261021 fixed-T nine-panel layout, unchanged from v14"
        require(not validate_manifest_record(record, hash_cache=cache), "invalid augmented figure record")
        chart.update(relative_stem=relative.as_posix(), arguments=list(arguments))
        charts.append(chart)
        print(f"[{len(charts):02d}/76] {relative.as_posix()}", flush=True)
    require(len(charts) == 76, "unexpected figure inventory")
    package = dict(case_slug=case, single_figure_count=72, composite_figure_count=3, detail_figure_count=1,
        color_figure_count=76, grayscale_figure_count=76, raw_input_point_count=1867, source_selection=COUNTS,
        numerical_data_postprocessing="none; upstream selection of computed refinements is explicit",
        delivery_stage="png_review", manuscript_eligible=False, current_publication_layer=False,
        vector_delivery_pending=True, layout_reference="publication_clean_v14_pdf", charts=charts)
    path = figures / "plot_manifest.json"
    write_manifest(path, build_bundle(package, [chart["manifest_record"] for chart in charts]))
    errors = validate_manifest(path)
    require(not errors, "; ".join(errors))
    write_manifest(results / "quality_review.json", dict(status="passed", figures=76, png_count=152,
        original_vertices_verified=sum(c["manifest_record"]["validation"]["original_vertex_count"] for c in charts),
        source_selection=COUNTS, no_smoothing_or_interpolation=True, validation_errors=errors,
        manifest=input_record(path, role="validated_png_review_bundle")))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--png-review", action="store_true", required=True)
    parser.add_argument("--case", default=CASE)
    parser.add_argument("--production", type=Path, default=PRODUCTION)
    args = parser.parse_args()
    build_review(args.case, args.production)


if __name__ == "__main__":
    main()
