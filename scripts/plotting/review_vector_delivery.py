"""Prepare and validate a PNG review before exporting its matching vector figure."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from scripts.plotting.plot_bundle import BUNDLE_SCHEMA, expand_bundle
from scripts.plotting.plot_delivery import validate_png_delivery
from scripts.plotting.plot_manifest import input_record, sha256_file
from scripts.plotting.plot_provenance import validate_hash_record
from scripts.plotting.plot_quality import export_figure, placement_limits


COLOR_ROUTE = "color_print_and_online"


@dataclass(frozen=True)
class AcceptedReview:
    record: dict
    author_review: dict


def add_vector_arguments(parser, *, multiple=False):
    parser.add_argument("--deliver-vector", action="store_true",
                        help="Export a vector PDF after verifying the accepted PNG review")
    parser.add_argument("--acceptance", type=Path, action="append" if multiple else "store",
                        help="Author's hash-bound PNG acceptance" + (" (repeat for each figure)" if multiple else ""))


def require_paired_vector_arguments(parser, args):
    if args.deliver_vector != bool(args.acceptance):
        parser.error("--deliver-vector and --acceptance must be used together")


def require_fresh_case(*paths: Path):
    for path in paths:
        if path.exists() and (not path.is_dir() or any(path.iterdir())):
            raise FileExistsError(f"Choose a new case directory: {path}")


def read_accepted_review(path: Path, *, root: Path, generator: Path) -> AcceptedReview:
    """Check authorization and the accepted inputs before drawing or writing anything."""
    author_review = dict(status="author_accepted_png",
                         record=input_record(path, role="author_png_acceptance", root=root))
    evidence = json.loads(path.read_text(encoding="utf-8-sig"))
    if (not isinstance(evidence, dict) or evidence.get("schema") != "plot_png_acceptance_v1"
            or evidence.get("status") != "author_accepted"
            or evidence.get("pdf_export_authorized") is not True
            or not isinstance(evidence.get("author_instruction"), str)
            or not evidence["author_instruction"].strip()):
        raise ValueError("PNG acceptance must record an actual author instruction and PDF authorization")
    manifest_pin = evidence.get("accepted_manifest")
    if not isinstance(manifest_pin, dict):
        raise ValueError("acceptance must bind the accepted_manifest")
    errors = validate_hash_record(manifest_pin, root=root, label="accepted_manifest", allow_historical=False)
    if errors:
        raise ValueError("; ".join(errors))
    record = json.loads((root / manifest_pin["path"]).read_text(encoding="utf-8-sig"))
    if record.get("schema_version") == BUNDLE_SCHEMA:
        accepted_id = evidence.get("accepted_figure_id")
        matches = [item for _, item in expand_bundle(record) if item["asset_id"] == accepted_id]
        if len(matches) != 1:
            raise ValueError("accepted_figure_id must select exactly one bundled figure")
        record = matches[0]

    generator_pin = record.get("generator", {})
    if ((root / generator_pin.get("path", "")).resolve() != generator.resolve()
            or generator_pin.get("sha256") != sha256_file(generator)):
        raise ValueError("PNG/vector generator hash or path mismatch")
    if record.get("rendering", {}).get("color_route") != COLOR_ROUTE:
        raise ValueError("accepted manifest must declare color_print_and_online")

    # Verify the declared acceptance against the original review, then verify
    # the exact bytes in its dependencies and both review outputs on disk.
    candidate = dict(record, rendering=dict(record["rendering"], vector_delivery_pending=False),
                     author_review=author_review)
    errors = validate_png_delivery(candidate, root=root)
    for collection, label in (("inputs", "frozen input"), ("outputs", "accepted output")):
        for pin in record[collection]:
            errors.extend(validate_hash_record(pin, root=root, label=f"{label} {pin['path']}",
                                               allow_historical=False))
    if errors:
        raise ValueError("; ".join(errors))
    return AcceptedReview(record, author_review)


def add_vector_pdf(record: dict, acceptance: AcceptedReview, *, figure, stem: Path, profile, root: Path):
    """Verify the re-rendered review against acceptance before writing the PDF."""
    if record["asset_id"] != acceptance.record["asset_id"]:
        raise ValueError("acceptance selects a different figure")

    def pngs(value):
        return sorted((item["role"], item["sha256"], item["bytes"])
                      for item in value["outputs"] if item["format"] == "png")

    if pngs(record) != pngs(acceptance.record):
        raise ValueError("color or grayscale PNG differs from the frozen review")
    record.update(delivery_stage="vector_delivery", vector_delivery_pending=False,
                  manuscript_eligible=False, current_publication_layer=False,
                  author_review=acceptance.author_review)
    rendering = record["rendering"]
    rendering.update(delivery_stage="vector_delivery", vector_delivery_pending=False,
                     color_route=COLOR_ROUTE, output_formats=["png", "pdf"])
    errors = validate_png_delivery(json.loads(json.dumps(record)), root=root)
    if errors:
        raise ValueError("; ".join(errors))
    pdfs, quality = export_figure(figure, stem, profile, formats=["pdf"],  # also inspects the actual PDF
                                   intended_width_inches=rendering["figure_size_inches"][0])
    inspection = pdfs[0]["inspection"]
    if (inspection["page_count"] != 1 or not inspection["fonts_embedded"]
            or inspection["type3_font_count"] or inspection["raster_image_count"]):
        raise ValueError(f"invalid single-page vector PDF: {stem}")
    if quality["actual_fonts"] != rendering["quality"]["actual_fonts"]:
        raise ValueError("vector export changed the measured font sources")
    record["outputs"].extend(pdfs)
    rendering["placement_limits"] = placement_limits(rendering["quality"], profile, record["outputs"])
    record["validation"]["png_rerender_byte_identical"] = True
