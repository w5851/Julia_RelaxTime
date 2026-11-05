#!/usr/bin/env python3
"""Export the accepted v15 PNG bundle with identical rendering and frozen inputs."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from copy import deepcopy
from dataclasses import replace
import hashlib
import io
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from matplotlib import font_manager
import matplotlib.pyplot as plt

from scripts.analysis.relaxtime import build_phase_guided_publication_clean_v15 as frozen
from scripts.plotting.plot_bundle import build_bundle, expand_bundle
from scripts.plotting.plot_manifest import input_record, generator_record, sha256_file, write_manifest
from scripts.plotting.plot_quality import export_figure, placement_limits
from scripts.plotting.plot_style import figure_style_context
from scripts.plotting.validate_plot_artifact import validate_manifest, validate_manifest_record

CASE = "publication_clean_v15_pdf"


def checked_authorization(path, source_manifest):
    value = frozen.raw.read_json(path)
    frozen.require(value.get("schema") == "transport_v15_author_acceptance_v1"
                   and all(value.get(key) is True for key in
                           ("png_accepted", "pdf_export_authorized", "figure_promotion_authorized",
                            "current_update_authorized", "manuscript_use_authorized"))
                   and bool(value.get("author_instruction", "").strip()),
                   "explicit v15 author acceptance and delivery authorization required")
    frozen.require(value.get("review_manifest_sha256") == sha256_file(source_manifest),
                   "authorization binds a different PNG bundle")
    return value


def normalized_path(path):
    return str(Path(path).resolve()).casefold()


@contextmanager
def exact_fonts(pairs, sources):
    """Temporarily resolve retained font bytes without changing system fonts."""
    expected = {}
    for _, record in pairs:
        rendering = record["rendering"]
        for item in [rendering["font"], *rendering["quality"]["actual_fonts"]]:
            previous = expected.setdefault(item["path"], item)
            frozen.require((previous["sha256"], previous["bytes"]) == (item["sha256"], item["bytes"]),
                           "inconsistent accepted font identity")
    frozen.require(set(sources) <= set(expected), "font map contains a source outside the accepted bundle")
    for original, item in expected.items():
        path = Path(sources.get(original, original))
        frozen.require(path.is_file() and path.stat().st_size == item["bytes"]
                       and sha256_file(path) == item["sha256"], f"accepted font bytes unavailable: {original}")
    mapping = {normalized_path(key): str(Path(value).resolve()) for key, value in sources.items()}
    previous = font_manager.fontManager.ttflist
    # FontEntry is immutable. Replace entries only for this process/context.
    font_manager.fontManager.ttflist = [replace(item, fname=mapping[normalized_path(item.fname)])
        if normalized_path(item.fname) in mapping else item for item in previous]
    font_manager.fontManager._findfont_cached.cache_clear()
    try:
        yield
    finally:
        font_manager.fontManager.ttflist = previous
        font_manager.fontManager._findfont_cached.cache_clear()


def canonical_font_evidence(quality, accepted, sources):
    """Keep original font identities after verifying the relocated measurements."""
    reverse = {normalized_path(value): key for key, value in sources.items()}
    measured = deepcopy(quality["actual_fonts"])
    for item in measured:
        item["path"] = reverse.get(normalized_path(item["path"]), item["path"])
    measured.sort(key=lambda item: item["path"])
    expected = sorted(accepted["rendering"]["quality"]["actual_fonts"], key=lambda item: item["path"])
    frozen.require(measured == expected, "measured font contents or roles differ from accepted PNG")
    quality["actual_fonts"] = deepcopy(accepted["rendering"]["quality"]["actual_fonts"])
    if "font" in quality:
        actual = deepcopy(quality["font"])
        actual["path"] = reverse.get(normalized_path(actual["path"]), actual["path"])
        frozen.require(actual == accepted["rendering"]["font"], "resolved font differs from accepted PNG")
        quality["font"] = actual


def export(case, authorization_path, font_sources=None):
    if (not case or case in {".", ".."} or Path(case).name != case
            or any(char in case for char in ("/", "\\", ":"))):
        raise ValueError("delivery requires a single new case name")
    figures, results = frozen.raw.FIGURE_PARENT / case, frozen.raw.RESULT_PARENT / case
    if figures.exists() or results.exists():
        raise FileExistsError("refusing to overwrite an existing delivery")
    source_manifest = frozen.raw.FIGURE_PARENT / frozen.CASE / "plot_manifest.json"
    authorization = checked_authorization(authorization_path, source_manifest)
    sources = font_sources or {}
    errors = validate_manifest(source_manifest, font_sources=sources)
    frozen.require(not errors, "; ".join(errors))
    bundle = frozen.raw.read_json(source_manifest)
    pairs = expand_bundle(bundle)
    frozen.require(len(pairs) == 76 and bundle["package"]["raw_input_point_count"] == 1867,
                   "wrong v15 figure inventory")
    # The exporter is new; every rendering module remains the accepted version.
    frozen.require(all(record["generator"]["sha256"] == sha256_file(Path(frozen.__file__))
                       for _, record in pairs), "PNG generator changed after acceptance")
    snapshot_root = frozen.raw.RESULT_PARENT / frozen.CASE / "source"
    for relative in ("scripts/analysis/relaxtime/build_phase_guided_publication_clean_v14.py",
                     "scripts/analysis/relaxtime/build_phase_guided_full_rs_raw.py",
                     "scripts/plotting/plot_style.py", "scripts/plotting/plot_quality.py",
                     "scripts/plotting/plot_series.py", "scripts/plotting/plot_accessibility.py"):
        frozen.require(sha256_file(ROOT / relative) == sha256_file(snapshot_root / relative),
                       f"accepted rendering dependency changed: {relative}")
    data, profile = frozen.load_inputs(), frozen.layout.publication_profile()
    figures.mkdir(parents=True)
    results.mkdir(parents=True)
    shutil.copyfile(authorization_path, results / "author_authorization.json")
    write_manifest(results / "font_sources.json", sources)
    charts, cache = [], {}
    with exact_fonts(pairs, sources):
        for chart, accepted in pairs:
            relative = Path(chart["relative_stem"])
            figure, specs, placements = frozen.layout.render_request(
                chart["kind"], tuple(chart["arguments"]), data, profile)
            try:
                frozen.require(frozen.layout.render_signature(figure, specs, placements)
                               == accepted["frozen_render_signature"], f"display semantics changed: {relative}")
                frozen.raw.verify_raw_vertices(figure, specs, data)
                figure.set_dpi(profile.dpi)
                buffer = io.BytesIO()
                with figure_style_context(figure):
                    figure.savefig(buffer, format="png", dpi=profile.dpi, bbox_inches=None)
                frozen.require(hashlib.sha256(buffer.getvalue()).hexdigest() == accepted["outputs"][0]["sha256"],
                               f"PNG rerender differs from accepted bytes: {relative}")
                pdfs, quality = export_figure(figure, figures / relative, profile, formats=["pdf"])
                canonical_font_evidence(quality, accepted, sources)
            finally:
                plt.close(figure)
            # Reuse the accepted PNGs. A vector companion needs no duplicate PNG/source tree.
            outputs = deepcopy(accepted["outputs"]) + pdfs
            acceptance_path = results / "png_acceptance" / relative.with_suffix(".json")
            write_manifest(acceptance_path, dict(schema="plot_png_acceptance_v1", status="author_accepted",
                pdf_export_authorized=True, author_instruction=authorization["author_instruction"],
                acceptance_basis=authorization["acceptance_basis"],
                accepted_manifest=input_record(source_manifest, role="accepted_png_manifest"),
                accepted_figure_id=accepted["asset_id"], accepted_png=accepted["outputs"][0]))
            record = deepcopy(accepted)
            record.update(asset_id=f"relaxtime.phase_guided.{case}.{relative.as_posix()}", case_slug=case,
                figure_mode="strict", semantic_status="formal_full_rs_eta_augmented_transport",
                publication_scope="main_text_candidate", outputs=outputs, manuscript_eligible=True,
                current_publication_layer=True, formal_production_promoted=True,
                delivery_stage="vector_delivery", vector_delivery_pending=False,
                author_review=dict(status="author_accepted_png",
                    record=input_record(acceptance_path, role="author_png_acceptance")))
            record["validation"].update(strict_gate=True, png_rerender_byte_identical=True)
            record["rendering"].update(quality=quality, delivery_stage="vector_delivery",
                vector_delivery_pending=False, output_formats=["png", "pdf"],
                placement_limits=placement_limits(quality, profile, outputs))
            errors = validate_manifest_record(record, hash_cache=cache, font_sources=sources)
            frozen.require(not errors, f"{relative}: " + "; ".join(errors))
            charts.append(dict(chart, stem=frozen.raw.display_path(figures / relative),
                               outputs=outputs, manifest_record=record))
            print(f"[{len(charts):02d}/76] PDF and identical PNG verified: {relative.as_posix()}", flush=True)
    package = dict(bundle["package"], case_slug=case, vector_figure_count=76, charts=charts,
        manuscript_eligible=True, current_publication_layer=True, vector_delivery_pending=False,
        delivery_stage="vector_delivery", png_storage="reference_accepted_png_case",
        vector_exporter=generator_record(Path(__file__), command=(
            "python scripts/analysis/relaxtime/export_phase_guided_publication_clean_v15_pdf.py "
            f"--case {case} --authorization {frozen.raw.display_path(authorization_path)}")),
        font_sources_record=input_record(results / "font_sources.json", role="identical_font_locations"))
    path = figures / "plot_manifest.json"
    write_manifest(path, build_bundle(package, [chart["manifest_record"] for chart in charts]))
    errors = validate_manifest(path, font_sources=sources)
    frozen.require(not errors, "; ".join(errors))
    write_manifest(results / "delivery_validation.json", dict(status="passed", figures=76,
        pdf_count=76, referenced_png_count=152, duplicated_png_count=0,
        all_png_rerenders_byte_identical=True, all_pdfs_single_page_embedded_non_type3_fonts_and_vector=True,
        no_smoothing_or_interpolation=True, source_selection=frozen.COUNTS,
        manifest=input_record(path, role="validated_delivery_bundle"), validation_errors=errors))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", default=CASE)
    parser.add_argument("--authorization", required=True, type=Path)
    parser.add_argument("--font-sources", type=Path)
    args = parser.parse_args()
    sources = frozen.raw.read_json(args.font_sources) if args.font_sources else None
    export(args.case, args.authorization, sources)


if __name__ == "__main__":
    main()
