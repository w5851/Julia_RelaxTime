#!/usr/bin/env python3
"""Measure a disposable manuscript insertion and prepare A4 print previews.

Run with the bundled document Python (pdfplumber and pypdf). The manuscript
project and accepted PNG/PDF figure packages are read-only inputs.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys

import pdfplumber
from pypdf import PdfReader, PdfWriter, Transformation

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.plotting.plot_bundle import load_chart_records

PNG_ROOT = ROOT / "data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v11_png_review/composites"
PNG_MANIFEST = PNG_ROOT.parent / "plot_manifest.json"
DEFAULT_REPORT_ROOT = ROOT / "docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v11_paper_size_review"
FIGURES = ("figure1_relaxation_times_comparison", "figure2_transport_coefficients_comparison")
MM_PER_PDF_POINT = 25.4 / 72


def file_record(path: Path) -> dict:
    return {"path": str(path.resolve()), "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def figure_manifest(figure: str) -> dict:
    _, records = load_chart_records(PNG_MANIFEST, root=ROOT)
    stem = (PNG_ROOT / figure).relative_to(ROOT).as_posix()
    matches = [record for chart, record in records if chart.get("stem") == stem]
    if len(matches) != 1:
        raise ValueError(f"expected exactly one composite record for {figure}, got {len(matches)}")
    return matches[0]


def a4_scale(page_width_pt: float, page_height_pt: float, margin_mm: float) -> float:
    return min((210 - 2 * margin_mm) / (page_width_pt * MM_PER_PDF_POINT),
               (297 - 2 * margin_mm) / (page_height_pt * MM_PER_PDF_POINT))


def a4_preview(source: Path, destination: Path, pages: list[int], margin_mm: float) -> None:
    if destination.exists():
        raise FileExistsError(f"refusing to overwrite print preview: {destination}")
    reader, writer = PdfReader(source), PdfWriter()
    width, height = 210 / MM_PER_PDF_POINT, 297 / MM_PER_PDF_POINT
    for index in pages:
        page = reader.pages[index]
        source_width, source_height = float(page.mediabox.width), float(page.mediabox.height)
        scale = a4_scale(source_width, source_height, margin_mm)
        target = writer.add_blank_page(width=width, height=height)
        transform = Transformation().scale(scale).translate(
            (width - scale * source_width) / 2, (height - scale * source_height) / 2)
        target.merge_transformed_page(copy.deepcopy(page), transform, expand=False)
    writer.write(destination)


def measure_placement(pdf: Path) -> list[dict]:
    results = []
    with pdfplumber.open(pdf) as document:
        for figure in FIGURES:
            manifest = figure_manifest(figure)
            quality = manifest["rendering"]["quality"]
            pixels = tuple(manifest["outputs"][0]["inspection"]["size_pixels"])
            matches = [(index, page, image) for index, page in enumerate(document.pages)
                       for image in page.images if tuple(image["srcsize"]) == pixels]
            if len(matches) != 1:
                raise ValueError(f"expected exactly one inserted PNG for {figure}, got {len(matches)}")
            index, page, image = matches[0]
            native_width_in, native_height_in = quality["figure_size_inches"]
            placement_scale = image["width"] / (72 * native_width_in)
            caption = [word for word in page.extract_words() if word["top"] > image["bottom"]
                       and word["text"].startswith("FIG.")]
            if len(caption) != 1:
                raise ValueError(f"missing caption below {figure}")
            caption_top = caption[0]["top"]
            scenarios = []
            for name, margin in (("A4_fit_page", 0), ("A4_fit_5mm_printable_margin", 5)):
                scale = a4_scale(page.width, page.height, margin)
                scenarios.append({
                    "name": name, "margin_mm": margin, "print_scale": scale,
                    "figure_width_mm": image["width"] * MM_PER_PDF_POINT * scale,
                    "figure_height_mm": image["height"] * MM_PER_PDF_POINT * scale,
                    "minimum_primary_glyph_height_mm": quality["minimum_primary_capital_numeral_height_mm"] * placement_scale * scale,
                    "minimum_math_script_height_mm": quality["minimum_script_capital_numeral_height_mm"] * placement_scale * scale,
                    "effective_png_dpi": pixels[0] / (image["width"] / 72 * scale),
                })
            results.append({
                "figure": figure, "page": index + 1,
                "page_size_pdf_pt": [page.width, page.height], "page_size_mm": [page.width * MM_PER_PDF_POINT, page.height * MM_PER_PDF_POINT],
                "native_size_inches": [native_width_in, native_height_in],
                "placement_scale": placement_scale,
                "placed_width_pdf_pt": image["width"], "placed_width_mm": image["width"] * MM_PER_PDF_POINT,
                "placed_height_mm": image["height"] * MM_PER_PDF_POINT,
                "minimum_primary_glyph_height_mm": quality["minimum_primary_capital_numeral_height_mm"] * placement_scale,
                "minimum_math_script_height_mm": quality["minimum_script_capital_numeral_height_mm"] * placement_scale,
                "effective_png_dpi": pixels[0] / (image["width"] / 72),
                "image_bbox_pdf_pt": [image["x0"], image["top"], image["x1"], image["bottom"]],
                "caption_top_pdf_pt": caption_top,
                "caption_gap_mm": (caption_top - image["bottom"]) * MM_PER_PDF_POINT,
                "print_scenarios": scenarios,
                "native_png_manifest": {**file_record(PNG_MANIFEST), "figure_id": manifest["asset_id"]},
            })
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--png-placement-pdf", type=Path, required=True)
    parser.add_argument("--vector-paper-preview", type=Path, required=True)
    parser.add_argument("--preview-tex", type=Path, required=True)
    parser.add_argument("--preview-root", type=Path, required=True)
    parser.add_argument("--report-root", type=Path, default=DEFAULT_REPORT_ROOT)
    args = parser.parse_args()
    if args.report_root.exists():
        raise FileExistsError("refusing to overwrite a paper-size review")
    results = measure_placement(args.png_placement_pdf)
    pages = [record["page"] - 1 for record in results]
    outputs = []
    for name, margin in (("a4_fit_page", 0), ("a4_fit_5mm_margin", 5)):
        destination = args.preview_root / f"{name}.pdf"
        a4_preview(args.vector_paper_preview, destination, pages, margin)
        outputs.append(file_record(destination))
    original = ROOT.parent.parent / "Desktop/paper/My_Paper/pnjl_aniso_transport"
    source_paths = [args.png_placement_pdf, args.vector_paper_preview, args.preview_tex,
                    original / "main.tex", original / "main.pdf", Path(__file__)]
    report = {
        "schema": "publication_clean_v11_manuscript_placement_review_v1",
        "task_classification": "independent", "manuscript_eligible": False, "paper_project_modified": False,
        "method": "PNG image bounding boxes in disposable manuscript; PDF companions keep the identical physical canvas",
        "units": {"PDF_point": "1/72 inch", "TeX_point": "1/72.27 inch"},
        "placement": results, "inputs": [file_record(path) for path in source_paths], "outputs": outputs,
        "float_options": {"original": "[t]", "temporary_preview": "[!t]",
                          "reason": "large figure plus expanded caption exceeded top-float fraction; !t removes this restriction without shrinking"},
        "unresolved_formal_gates": [
            "math scripts below the repository 2 mm final-size threshold",
            "paper textwidth is approximately 7.057 in, above the candidate profile 7 in limit",
            "600 dpi native PNG falls below 600 effective dpi when enlarged at textwidth; PDF is vector and has no raster DPI limit",
        ],
        "scope": "visual size and export review only; frozen data, eligibility and current=v5 unchanged",
    }
    args.report_root.mkdir(parents=True, exist_ok=False)
    (args.report_root / "placement_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    first = results[0]
    rows = []
    for record in results:
        rows.append(f"| {record['figure']} | {record['page']} | {record['placed_width_mm']:.2f} | {record['placed_height_mm']:.2f} |")
    printing = []
    for scenario in first["print_scenarios"]:
        printing.append(f"| {scenario['name']} | {100 * scenario['print_scale']:.2f}% | {scenario['figure_width_mm']:.2f} | {scenario['minimum_primary_glyph_height_mm']:.2f} | {scenario['minimum_math_script_height_mm']:.2f} |")
    text = f"""# v11 manuscript placement and A4 print review

The current manuscript is Letter (215.9 x 279.4 mm), not A4. In the disposable
RevTeX preview, textwidth inserts the 6.75-in native figures at
{first['placed_width_mm']:.2f} mm, or {100 * first['placement_scale']:.2f}% of native size.
The main glyphs are at least {first['minimum_primary_glyph_height_mm']:.2f} mm;
math scripts are approximately {first['minimum_math_script_height_mm']:.2f} mm.
PNG effective resolution at that placement is {first['effective_png_dpi']:.1f} dpi.

| Figure | Page | Width (mm) | Height (mm) |
| --- | ---: | ---: | ---: |
""" + "\n".join(rows) + """

| A4 printing assumption | Scale | Figure width (mm) | Main glyph (mm) | Script (mm) |
| --- | ---: | ---: | ---: | ---: |
""" + "\n".join(printing) + """

Fit-page is a geometric page-box calculation. A real printer may impose a
different printable area; the 5-mm case is an explicit example, not a claim
about the user's printer. Actual-size printing keeps the Letter PDF glyph
sizes, subject to the printer's printable-area limits.

The temporary source updates only paths, captions and figure float options.
The original [t] setting produced float-placement warnings for the taller
figure/caption; [!t] in the preview allows the available page area and keeps
10 pages without those warnings. The formal paper files are not changed.

Visual placement is usable at the tested width. This is not formal submission
compliance: math scripts remain below 2 mm, textwidth slightly exceeds the
repository 7-in profile cap, and an enlarged PNG falls below 600 effective
dpi. The review PDFs remove only the raster-resolution limitation, not the
font-size or width gates. Do not silently promote manuscript eligibility.

The two A4 print-preview PDFs are disposable, two-page excerpts of the vector
manuscript preview. See placement_report.json for absolute paths and hashes.
"""
    (args.report_root / "README.md").write_text(text, encoding="utf-8")
    print(json.dumps({"report_root": str(args.report_root), "placement": results}, indent=2))


if __name__ == "__main__":
    main()
