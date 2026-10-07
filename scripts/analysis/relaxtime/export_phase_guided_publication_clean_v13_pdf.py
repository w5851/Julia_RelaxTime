#!/usr/bin/env python3
"""Export the accepted v13 PNG case as vector PDFs using its frozen renderer.

The PNG bundle stays immutable. Every new PDF records the original PNGs,
version-level author acceptance, exact renderer reproduction, and PDF inspection.
"""

from __future__ import annotations

import argparse
import copy
import datetime as dt
import hashlib
import io
import json
import os
from pathlib import Path
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from PIL import Image

from scripts.analysis.relaxtime import build_phase_guided_publication_clean_v13 as V13
from scripts.plotting.plot_bundle import build_bundle, load_chart_records
from scripts.plotting.plot_manifest import (
    generator_record, input_record, runtime_record, sha256_file, write_manifest,
)
from scripts.plotting.plot_provenance import validate_hash_record
from scripts.plotting.plot_quality import export_figure, placement_limits
from scripts.plotting.plot_style import configure_matplotlib, load_profile
from scripts.plotting.validate_plot_artifact import validate_manifest, validate_manifest_record

FIGURE_ROOT = V13.FIGURE_ROOT.with_name("publication_clean_v13_pdf")
ANALYSIS_ROOT = V13.ANALYSIS_ROOT.with_name("phase_guided_transport_publication_clean_v13_pdf")
PNG_PACKAGE = V13.ANALYSIS_ROOT / "manifest.json"
PNG_INDEX = V13.FIGURE_ROOT / "plot_manifest.json"
ACCEPTANCE = V13.ANALYSIS_ROOT.parent / "publication_clean_v13_png_acceptance_v1.json"
SOP = ROOT / "docs/guides/sop/workflows/figure_production.md"
COUNTS = {"single": 72, "composite": 2, "detail": 1}

# Everything else, including numerical status and the renderer generator,
# must remain exactly the same as the accepted PNG record.
CHANGED_FIELDS = {"generated_at_utc", "asset_id", "outputs", "rendering",
                  "delivery_stage", "vector_delivery_pending"}
CHANGED_RENDERING = {"quality", "placement_limits", "pdf_placement_limits", "caption_handoff",
                     "delivery_stage", "vector_delivery_pending", "output_formats"}


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def require_current(record, *, root=ROOT):
    """Require live bytes; archived source cannot authorize executing drifted code."""
    path = root / record["path"]
    if (not path.is_file() or sha256_file(path) != record["sha256"]
            or ("bytes" in record and path.stat().st_size != record["bytes"])):
        raise ValueError(f"frozen live file changed: {record['path']}")


def require_acceptance(acceptance):
    if (acceptance.get("schema") != "publication_clean_v13_png_acceptance_v1"
            or acceptance.get("author_png_accepted") is not True
            or acceptance.get("pdf_export_authorized") is not True
            or acceptance.get("manuscript_eligible") is not False
            or acceptance.get("current_publication_layer") is not False):
        raise ValueError("v13 PDF export requires explicit PNG acceptance and unchanged eligibility")


def verify_source_snapshot(acceptance, package):
    snapshot = acceptance["source_snapshot"]
    require_current(snapshot)
    with zipfile.ZipFile(ROOT / snapshot["path"]) as archive:
        metadata = json.loads(archive.read("snapshot_manifest.json"))
        if metadata["source_package_sha256"] != acceptance["png_package"]["sha256"]:
            raise ValueError("source archive is not bound to the accepted PNG package")
        archived = {item["path"]: item for item in metadata["files"]}
        for path, item in archived.items():
            payload = archive.read(path)
            if len(payload) != item["bytes"] or hashlib.sha256(payload).hexdigest() != item["sha256"]:
                raise ValueError(f"invalid archived source: {path}")
        executable = []
        for record in [*package["inputs"], package["generator"]]:
            path = record["path"]
            runs = ((path.startswith("scripts/") and Path(path).suffix == ".py")
                    or path == "config/plotting/candidate_aps_v2.toml")
            if runs:
                if path not in archived or archived[path]["sha256"] != record["sha256"]:
                    raise ValueError(f"renderer dependency missing from accepted snapshot: {path}")
                require_current(record)
                executable.append(path)
    return sorted(set(executable))


def load_accepted_case():
    acceptance = read_json(ACCEPTANCE)
    require_acceptance(acceptance)
    for key, expected in (("png_package", PNG_PACKAGE), ("png_bundle", PNG_INDEX)):
        if (ROOT / acceptance[key]["path"]).resolve() != expected.resolve():
            raise ValueError(f"acceptance references a different {key}")
        require_current(acceptance[key])
    package = read_json(PNG_PACKAGE)
    live_sources = verify_source_snapshot(acceptance, package)
    for record in package["inputs"]:
        if errors := validate_hash_record(record, root=ROOT, label="accepted inputs"):
            raise ValueError(f"accepted input changed: {errors}")
    for record in package["outputs"]:
        require_current(record)
    if errors := validate_manifest(PNG_INDEX, repo_root=ROOT):
        raise ValueError(f"invalid accepted PNG bundle: {errors}")
    index, pairs = load_chart_records(PNG_INDEX, root=ROOT)
    if {kind: sum(chart["kind"] == kind for chart, _ in pairs) for kind in COUNTS} != COUNTS or len(pairs) != 75:
        raise ValueError("v13 requires 72 singles, two main composites, and one detail figure")
    for _, record in pairs:
        if (record["generator"] != package["generator"] or record["inputs"] != package["inputs"]
                or record["rendering"]["typography_exception"] is not None):
            raise ValueError("accepted PNG record has inconsistent source or typography")
    return acceptance, package, index, pairs, live_sources


def rgba_digest(image):
    converted = image.convert("RGBA")
    return {"size_pixels": list(converted.size),
            "rgba_pixel_sha256": hashlib.sha256(converted.tobytes()).hexdigest()}


def confirm_png_pixels(figure, reference, dpi):
    """Compare actual RGBA pixels, ignoring non-visual PNG metadata."""
    with io.BytesIO() as buffer:
        figure.savefig(buffer, format="png", dpi=dpi, bbox_inches=None)
        buffer.seek(0)
        with Image.open(buffer) as generated, Image.open(reference) as retained:
            actual, expected = rgba_digest(generated), rgba_digest(retained)
    if actual != expected:
        raise ValueError(f"renderer no longer reproduces the accepted PNG: {reference}")
    return {"pixel_identical": True, "method": "SHA-256 of decoded RGBA pixels", **actual}


def make_vector_record(source, outputs, quality, profile, export_provenance, pixels):
    record = copy.deepcopy(source)
    record.update(
        generated_at_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
        asset_id=source["asset_id"].replace(".png_review.", ".vector_delivery."),
        outputs=[*copy.deepcopy(source["outputs"]), *outputs],
        delivery_stage="vector_delivery", vector_delivery_pending=False,
        vector_delivery_complete=True, author_png_accepted=True,
        vector_export=copy.deepcopy(export_provenance),
        source_png_manifest={**input_record(PNG_INDEX, role="accepted_png_bundle"),
                             "figure_id": source["asset_id"]},
        source_png_pixel_comparison=pixels,
    )
    record["rendering"].update(
        delivery_stage="vector_delivery", vector_delivery_pending=False,
        output_formats=["png", "pdf"], quality=quality,
        placement_limits=placement_limits(quality, profile, record["outputs"]),
        pdf_placement_limits=placement_limits(quality, profile, outputs),
        caption_handoff=(ANALYSIS_ROOT / "caption_handoff.md").relative_to(ROOT).as_posix(),
    )
    return record


def export_chart(chart, source, grouped, gap_map, profile, export_provenance):
    if chart["kind"] == "single":
        spec = source["rendering"]["panel_specs"][0]
        figure, specs, placements = V13.render_single(
            spec["mode_key"], spec["plot_panel"], spec["observable"], grouped, gap_map, profile)
    elif chart["kind"] == "composite":
        figure, specs, placements = V13.render_composite(source["case_slug"], grouped, gap_map, profile)
    elif chart["kind"] == "detail":
        figure, specs, placements = V13.render_detail(grouped, gap_map, profile)
    else:
        raise ValueError(f"unknown chart kind: {chart['kind']}")
    try:
        if specs != source["rendering"]["panel_specs"] or placements != source["rendering"]["legend_placements"]:
            raise ValueError("panel specifications or legend placements changed after PNG acceptance")
        figure.set_dpi(profile.dpi)
        pixels = confirm_png_pixels(figure, ROOT / source["outputs"][0]["path"], profile.dpi)
        stem = FIGURE_ROOT / (ROOT / chart["stem"]).relative_to(V13.FIGURE_ROOT)
        outputs, quality = export_figure(figure, stem, profile, formats=("pdf",))
        outputs[0]["role"] = "vector_delivery"
        record = make_vector_record(source, outputs, quality, profile, export_provenance, pixels)
        return {"stem": stem.relative_to(ROOT).as_posix(), "kind": chart["kind"],
                "mode_key": chart["mode_key"], "outputs": record["outputs"]}, record
    finally:
        V13.plt.close(figure)


def verify_derivative(source, record, profile):
    for key in source.keys() - CHANGED_FIELDS:
        if record.get(key) != source[key]:
            raise ValueError(f"accepted source field changed: {key}")
    for key in source["rendering"].keys() - CHANGED_RENDERING:
        if record["rendering"].get(key) != source["rendering"][key]:
            raise ValueError(f"accepted rendering field changed: {key}")
    if record["outputs"][:2] != source["outputs"] or len(record["outputs"]) != 3:
        raise ValueError("vector case must retain both accepted PNG references and add one PDF")
    for container in (record, record["rendering"]):
        if container.get("delivery_stage") != "vector_delivery" or container.get("vector_delivery_pending") is not False:
            raise ValueError("vector delivery stage is inconsistent")
    if record.get("author_png_accepted") is not True or record.get("vector_delivery_complete") is not True:
        raise ValueError("missing accepted vector delivery state")
    pdf = record["outputs"][2]
    if pdf["format"] != "pdf" or pdf["vector"] is not True:
        raise ValueError("the new output must be a vector PDF")
    limits = placement_limits(record["rendering"]["quality"], profile, [pdf])
    if record["rendering"].get("pdf_placement_limits") != limits or not limits["measured_width_qualified"]:
        raise ValueError("PDF placement limits do not match measured quality")
    if record["source_png_manifest"]["figure_id"] != source["asset_id"]:
        raise ValueError("wrong accepted PNG figure reference")
    require_current(record["source_png_manifest"])
    with Image.open(ROOT / source["outputs"][0]["path"]) as image:
        expected = rgba_digest(image)
    pixels = record["source_png_pixel_comparison"]
    if pixels.get("pixel_identical") is not True or any(pixels.get(key) != value for key, value in expected.items()):
        raise ValueError("accepted PNG pixel evidence does not match")


def export_provenance(acceptance, live_sources, font):
    return {
        "exporter": generator_record(Path(__file__),
            command="python scripts/analysis/relaxtime/export_phase_guided_publication_clean_v13_pdf.py",
            runtime=runtime_record({"matplotlib": V13.matplotlib.__version__, "font": font})),
        "author_acceptance": input_record(ACCEPTANCE, role="version_level_png_acceptance"),
        "source_package": acceptance["png_package"], "source_bundle": acceptance["png_bundle"],
        "source_snapshot": acceptance["source_snapshot"],
        "export_contract": input_record(SOP, role="active_vector_export_contract"),
        "live_renderer_dependencies_verified": live_sources,
        "acceptance_scope": "version-level acceptance; per-file author inspection is not asserted",
        "source_authority": "exact PNG-bound source archive plus matching live renderer and profile bytes",
    }


def write_handoff(pairs):
    caption = (V13.ANALYSIS_ROOT / "caption_handoff.md").read_text(encoding="utf-8")
    caption = caption.replace("# v13 图注与插入尺寸交接", "# v13 PDF 图注与插入尺寸交接", 1)
    caption = caption.replace("逐图 `rendering.placement_limits` 给出由字形、线宽、端点和 PNG 有效 DPI\n共同约束的插入区间。",
                              "逐图 `rendering.pdf_placement_limits` 给出矢量 PDF 由字形、线宽、端点\n和项目最大宽度共同约束的插入区间；PNG 有效 DPI 不限制矢量 PDF。")
    caption = caption.replace("PNG 审查层：作者接受待定，PDF 待后续交付，manuscript_eligible=false，current=v5。",
                              "v13 PNG 已获版本级人工接受；对应矢量 PDF 已交付。manuscript_eligible=false，current=v5。")
    (ANALYSIS_ROOT / "caption_handoff.md").write_text(caption, encoding="utf-8", newline="\n")
    report = {"schema": "publication_clean_v13_pdf_placement_report_v1",
              "recommended_width_mm": V13.WIDTH_IN * 25.4,
              "scope": "native figure dimensions; manuscript placement and physical printout not assessed",
              "figures": [{"pdf": record["outputs"][2], "kind": chart["kind"],
                           "pdf_placement_limits": record["rendering"]["pdf_placement_limits"],
                           "png_and_pdf_placement_limits": record["rendering"]["placement_limits"],
                           "accepted_png_pixel_comparison": record["source_png_pixel_comparison"]}
                          for chart, record in pairs]}
    write_manifest(ANALYSIS_ROOT / "placement_report.json", report)
    rows = ["# v13 PDF 图件索引", "", "全部图件保留已接受 PNG 的排版；推荐插入宽度 171.45 mm。", "",
            "| 图件 | 类型 | PDF |", "| --- | --- | --- |"]
    for chart, record in sorted(pairs, key=lambda item: (item[0]["kind"] == "single", item[0]["stem"])):
        pdf = ROOT / record["outputs"][2]["path"]
        link = Path(os.path.relpath(pdf, ANALYSIS_ROOT)).as_posix()
        label = pdf.relative_to(FIGURE_ROOT).with_suffix("").as_posix()
        rows.append(f"| {label} | {chart['kind']} | [PDF]({link}) |")
    (ANALYSIS_ROOT / "pdf_index.md").write_text("\n".join(rows) + "\n", encoding="utf-8", newline="\n")
    (ANALYSIS_ROOT / "README.md").write_text("""# publication_clean_v13 矢量 PDF

v13 PNG 已获作者版本级接受；本包交付 72 张单图、2 张主复合图、1 张局部图，
共 75 份单页彩色矢量 PDF。[图件索引](pdf_index.md)列出全部文件；
[图注交接](caption_handoff.md)保留科学说明与冻结温度映射。

所有图保持 171.45 mm（6.75 in）原生宽度。每张图均通过同一冻结 renderer
重绘与已接受 PNG 的 RGBA 像素逐一比对；坐标、曲线顶点、分支 gap 和图例不变。
PDF 逐文件检查物理尺寸、单页、嵌入字体、非 Type 3 字体和无嵌入栅格图像；
完整逐图测量及文件 hash 存于对应图像目录的唯一 `plot_manifest.json`。

[尺寸报告](placement_report.json)分别给出 PDF 与 PNG 的插入限制。
推荐维持原生宽度；最小实测字形超过 2 mm，未使用紧凑字形例外。
矢量 PDF 不受 PNG 分辨率上限约束；双栏图仍不可直接缩成单栏图。

原 PNG/灰度图及其 manifest 保持冻结，新 manifest 引用原文件。
接受记录仅表明版本级人工接受，不声称每份灰度图均被逐项人工审查。
本包完成 `vector_delivery`，保持 `manuscript_eligible=false`、current=v5；
没有调用求解器、改变数值数据、新增平滑或取得数值生产／论文装配资格。

```powershell
python scripts/analysis/relaxtime/export_phase_guided_publication_clean_v13_pdf.py --check
```

首次生成不带 `--check`；已存在的输出目录拒绝覆盖。
后续独立复核记录见 `verification/`（生成器本身不写入人工审阅结论）。
""", encoding="utf-8", newline="\n")


def build_delivery():
    if FIGURE_ROOT.exists() or ANALYSIS_ROOT.exists():
        raise FileExistsError("refusing to overwrite an existing v13 PDF case")
    acceptance, package, _, accepted, live_sources = load_accepted_case()
    print("[v13-pdf] accepted PNG bundle, frozen inputs and live renderer verified", flush=True)
    points, _, _, _, gaps, hashes = V13.V11.PARENT.load_v5_inputs()
    if len(points) != package["point_row_count"] or hashes != package["inherited_table_hashes"]:
        raise ValueError("frozen v5 point tables changed")
    grouped, gap_map = V13.V11.V10.group_inputs(points, gaps)
    profile = load_profile("candidate_aps_v2")
    font = configure_matplotlib(profile)
    provenance = export_provenance(acceptance, live_sources, font)
    ANALYSIS_ROOT.mkdir(parents=True, exist_ok=False)
    delivered, cache = [], {}
    for number, (chart, source) in enumerate(accepted, 1):
        new_chart, record = export_chart(chart, source, grouped, gap_map, profile, provenance)
        verify_derivative(source, record, profile)
        if errors := validate_manifest_record(record, repo_root=ROOT, hash_cache=cache):
            raise ValueError(f"invalid vector delivery: {new_chart['stem']}: {errors}")
        delivered.append((new_chart, record))
        if number % 6 == 0 or number == len(accepted):
            print(f"[v13-pdf] {number}/{len(accepted)} PDFs inspected; accepted PNG pixels identical", flush=True)
    verify_source_snapshot(acceptance, package)
    for record in (acceptance["png_package"], acceptance["png_bundle"], *package["outputs"]):
        require_current(record)
    write_handoff(delivered)
    figure_index = FIGURE_ROOT / "plot_manifest.json"
    metadata = {"schema": "publication_clean_v13_pdf_figure_index_v1", "status": "vector_delivery_complete",
                "single_figure_count": 72, "composite_figure_count": 2, "detail_figure_count": 1,
                "pdf_count": 75, "referenced_color_png_count": 75, "referenced_grayscale_png_count": 75,
                "mode_counts": {"mode_a": 36, "mode_b": 36}, "delivery_stage": "vector_delivery",
                "vector_delivery_pending": False, "manuscript_eligible": False, "current_publication_layer": False,
                "charts": [chart for chart, _ in delivered]}
    write_manifest(figure_index, build_bundle(metadata, [record for _, record in delivered]))
    write_manifest(ANALYSIS_ROOT / "manifest.json", {
        "schema": "publication_clean_v13_pdf_package_v1", "task_classification": "independent",
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "status": "vector_delivery_complete", "delivery_stage": "vector_delivery", "vector_delivery_pending": False,
        "manuscript_eligible": False, "current_publication_layer": False, "solver_called": False,
        "canonical_data_modified": False, "new_display_values": False,
        "vector_export": provenance, "point_row_count": len(points), "inherited_table_hashes": hashes,
        "figure_index": figure_index.relative_to(ROOT).as_posix(), "figure_index_sha256": sha256_file(figure_index),
        "outputs": [input_record(path, role="v13_pdf_delivery") for path in [figure_index,
                    *(ANALYSIS_ROOT / name for name in ("README.md", "caption_handoff.md", "pdf_index.md", "placement_report.json"))]],
        "validation_summary": {"pdf_count": len(delivered), "pixel_identical_count": len(delivered),
            "minimum_glyph_height_mm": min(record["rendering"]["quality"]["minimum_capital_numeral_height_mm"]
                                            for _, record in delivered),
            "all_single_page": True, "all_fonts_embedded_non_type3": True, "embedded_raster_image_count": 0,
            "scope": "measured figures and inspected PDFs; separate from numerical and manuscript qualification"},
    })
    print(f"[v13-pdf] completed: {FIGURE_ROOT}", flush=True)


def check_delivery():
    _, _, _, accepted, _ = load_accepted_case()
    package = read_json(ANALYSIS_ROOT / "manifest.json")
    for record in package["outputs"]:
        require_current(record)
    provenance = package["vector_export"]
    for key in ("exporter", "author_acceptance", "source_package", "source_bundle", "source_snapshot", "export_contract"):
        require_current(provenance[key])
    index_path = FIGURE_ROOT / "plot_manifest.json"
    if package["figure_index_sha256"] != sha256_file(index_path):
        raise ValueError("PDF bundle hash changed")
    index, delivered = load_chart_records(index_path, root=ROOT)
    if len(delivered) != 75 or {kind: sum(chart["kind"] == kind for chart, _ in delivered) for kind in COUNTS} != COUNTS:
        raise ValueError("PDF bundle figure counts differ from accepted PNG")
    if errors := validate_manifest(index_path, repo_root=ROOT):
        raise ValueError(f"PDF contract validation failed: {errors}")
    source_by_id = {record["asset_id"]: record for _, record in accepted}
    seen = set()
    profile = load_profile("candidate_aps_v2")
    for chart, record in delivered:
        source_id = record["source_png_manifest"]["figure_id"]
        if source_id in seen or source_id not in source_by_id:
            raise ValueError("PDF references a duplicate or unknown accepted PNG")
        seen.add(source_id)
        verify_derivative(source_by_id[source_id], record, profile)
        if record["vector_export"] != provenance or chart["outputs"] != record["outputs"]:
            raise ValueError("PDF chart or export provenance differs from package")
    actual = {path.resolve() for path in FIGURE_ROOT.rglob("*.pdf")}
    expected = {(ROOT / record["outputs"][2]["path"]).resolve() for _, record in delivered}
    if actual != expected or list(FIGURE_ROOT.rglob("*.png")):
        raise ValueError("PDF case inventory differs from manifest; original PNGs must only be referenced")
    if list(FIGURE_ROOT.rglob("*.json")) != [index_path]:
        raise ValueError("PDF case must contain exactly one plot manifest")
    print("[v13-pdf] OK: all 75 vector PDFs, accepted PNG references, source bytes and geometry validated", flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="verify the retained delivery without regenerating PDFs")
    args = parser.parse_args(argv)
    check_delivery() if args.check else build_delivery()


if __name__ == "__main__":
    main()
