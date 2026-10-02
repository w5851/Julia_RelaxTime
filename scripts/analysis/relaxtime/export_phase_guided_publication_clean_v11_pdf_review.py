#!/usr/bin/env python3
"""Export PDF review companions without changing the hash-bound v11 PNG case.

The retained renderer is reused directly. This is a technical PDF preflight,
not vector_delivery: unresolved final-size violations remain explicit.
"""

from __future__ import annotations

import copy
import argparse
import datetime as dt
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from PIL import Image

from scripts.plotting.plot_manifest import (
    generator_record, input_record, runtime_record, sha256_file, write_manifest,
)
from scripts.plotting.plot_quality import export_figure
from scripts.plotting.plot_style import configure_matplotlib, load_profile
from scripts.plotting.validate_plot_artifact import _check_final_size_contract, validate_manifest_record, validate_review_companion
from scripts.plotting.plot_bundle import build_bundle, load_chart_records
from scripts.plotting.plot_provenance import code_ref_for_manifest, validate_hash_record

V11_SCRIPT = ROOT / "scripts/analysis/relaxtime/build_phase_guided_publication_clean_v11.py"
SPEC = importlib.util.spec_from_file_location("retained_publication_v11_pdf_renderer", V11_SCRIPT)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load frozen renderer: {V11_SCRIPT}")
V11 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V11)

FIGURE_ROOT = ROOT / "data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v11_pdf_review"
ANALYSIS_ROOT = V11.V10.TRANSPORT_ANALYSIS_ROOT / "phase_guided_transport_publication_clean_v11_pdf_review"
PNG_PACKAGE = V11.ANALYSIS_ROOT / "manifest.json"
PNG_INDEX = V11.FIGURE_ROOT / "plot_manifest.json"
GLYPH_BLOCKER = "v2 final capital/numeral glyph height is below 2 mm"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def verify_records(records: list[dict[str, Any]]) -> None:
    for record in records:
        if errors := validate_hash_record(record, root=ROOT, label="frozen",
                                         code_ref=code_ref_for_manifest(PNG_INDEX, ROOT)):
            raise ValueError(f"frozen input/output changed: {errors}")


def load_frozen_case() -> tuple[dict, dict, list[dict]]:
    package = read_json(PNG_PACKAGE)
    index, pairs = load_chart_records(PNG_INDEX, root=ROOT)
    index["charts"] = [{**chart, "manifest_record": record} for chart, record in pairs]
    verify_records([*package["inputs"], *package["outputs"]])
    if len(index["charts"]) != 74:
        raise ValueError("v11 must contain 72 singles and two composites")
    retained = [input_record(PNG_PACKAGE, role="frozen_v11_png_package")]
    retained.extend(package["inputs"])
    retained.extend(package["outputs"])
    for chart, record in pairs:
        errors = validate_manifest_record(record, code_ref=code_ref_for_manifest(PNG_INDEX, ROOT))
        if errors:
            raise ValueError(f"invalid retained PNG: {chart}: {errors}")
    retained.append(input_record(PNG_INDEX, role="frozen_v11_png_bundle"))
    return package, index, retained


def confirm_composite_pixels(figure: Any, reference: Path, dpi: int) -> dict[str, Any]:
    """Require the reused renderer to reproduce the accepted PNG pixel-for-pixel."""
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=dpi, bbox_inches=None)
    buffer.seek(0)
    with Image.open(buffer) as generated, Image.open(reference) as retained:
        left, right = generated.convert("RGBA"), retained.convert("RGBA")
        expected = hashlib.sha256(right.tobytes()).hexdigest()
        actual = hashlib.sha256(left.tobytes()).hexdigest()
        if left.size != right.size or actual != expected:
            raise ValueError(f"renderer no longer reproduces accepted composite: {reference}")
        return {"pixel_identical": True, "size_pixels": list(left.size), "rgba_pixel_sha256": actual}


def submission_preflight(source: dict, output: dict, quality: dict, profile: Any) -> list[str]:
    candidate = copy.deepcopy(source)
    candidate["outputs"] = [*source["outputs"], output]
    candidate["rendering"].update({
        "delivery_stage": "vector_delivery", "output_formats": ["png", "pdf"],
        "typography_exception": None, "quality": quality,
    })
    errors: list[str] = []
    _check_final_size_contract(candidate, profile, ROOT, errors)
    return errors


def export_chart(chart: dict, grouped: dict, gap_map: dict, profile: Any, font: dict) -> dict:
    source_path = PNG_INDEX
    source = chart["manifest_record"]
    old_specs = source["rendering"]["panel_specs"]
    composite = chart["kind"] == "composite"
    if composite:
        figure, specs = V11.render_composite(V11.COMPOSITES[source["case_slug"]], grouped, gap_map, profile)
    else:
        first = old_specs[0]
        figure, specs = V11.render_single(first["mode_key"], first["plot_panel"], first["observable"], grouped, gap_map, profile)
    try:
        if json.dumps(specs, sort_keys=True) != json.dumps(old_specs, sort_keys=True):
            raise ValueError(f"plot specifications differ from the accepted PNG: {source_path}")
        reference = source["outputs"][0]
        pixels = confirm_composite_pixels(figure, ROOT / reference["path"], profile.dpi) if composite else None
        relative_stem = (ROOT / chart["stem"]).relative_to(V11.FIGURE_ROOT)
        stem = FIGURE_ROOT / relative_stem
        outputs, quality = export_figure(figure, stem, profile, formats=("pdf",))
        errors = submission_preflight(source, outputs[0], quality, profile)
        unexpected = [error for error in errors if error != GLYPH_BLOCKER]
        if unexpected:
            raise ValueError(f"unexpected PDF preflight failures: {source_path}: {unexpected}")
        # The original PNG case and its SOP remain frozen. A companion records
        # failures rather than extending the PNG-only exception to delivery.
        manifest = {
            "schema": "publication_clean_v11_pdf_review_companion_v1",
            "delivery_stage": "pdf_review_preflight", "status": "technical_pdf_review_only",
            "manuscript_eligible": False, "current_publication_layer": False,
            "formal_vector_delivery_complete": False, "solver_called": False,
            "canonical_data_modified": False, "new_display_values": False,
            "author_png_visual_acceptance": "user accepted the current PNG version before requesting paper-size checks and PDFs",
            "source_png_manifest": {**input_record(source_path, role="accepted_png_review_manifest"),
                                    "figure_id": source["asset_id"]},
            "source_png_output": reference, "source_png_pixel_comparison": pixels,
            "generator": generator_record(Path(__file__),
                command="python scripts/analysis/relaxtime/export_phase_guided_publication_clean_v11_pdf_review.py",
                runtime=runtime_record({"matplotlib": V11.matplotlib.__version__, "font": font})),
            "renderer": input_record(V11_SCRIPT, role="frozen_v11_renderer"),
            "kind": chart["kind"], "mode_key": chart["mode_key"], "panel_specs": specs,
            "axes": source["axes"], "series": source["series"],
            "selection_rule": source["selection_rule"], "interpolation_policy": source["interpolation_policy"],
            "connector_policy": source["connector_policy"], "missing_value_policy": source["missing_value_policy"],
            "quality": quality, "outputs": outputs,
            "submission_preflight": {"contract": "final_size_v2", "assessed_delivery_stage": "vector_delivery",
                "typography_exception_applied": False, "passed": not errors, "violations": errors,
                "scope": "native-size visual/export contract only; no numerical eligibility promotion"},
            "calculation_sha": source["calculation_sha"], "workflow_head_sha": source["workflow_head_sha"],
        }
        return {"manifest_record": manifest,
                "kind": chart["kind"], "mode_key": chart["mode_key"], "outputs": outputs,
                "submission_preflight_violations": errors}
    finally:
        V11.plt.close(figure)


def verify_companion(path: Path) -> list[str]:
    return verify_companion_record(read_json(path))


def verify_companion_record(manifest: dict) -> list[str]:
    errors = validate_review_companion(manifest, code_ref=code_ref_for_manifest(PNG_INDEX, ROOT))
    if any(error != GLYPH_BLOCKER for error in manifest["submission_preflight"]["violations"]):
        errors.append("unresolved non-typography export failures")
    return errors


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="validate existing companions without rendering")
    args = parser.parse_args(argv)
    if args.check:
        _, pairs = load_chart_records(FIGURE_ROOT / "plot_manifest.json", root=ROOT)
        for _, record in pairs:
            if errors := verify_companion_record(record):
                raise ValueError(f"invalid PDF companion: {errors}")
        print(f"[v11-pdf-review] verified {len(pairs)} retained PDFs; no rendering")
        return
    if FIGURE_ROOT.exists() or ANALYSIS_ROOT.exists():
        raise FileExistsError("refusing to overwrite an existing v11 PDF review companion")
    package, index, retained = load_frozen_case()
    points, _, _, _, gaps, _ = V11.PARENT.load_v5_inputs()
    grouped, gap_map = V11.V10.group_inputs(points, gaps)
    profile = load_profile("candidate_aps_v2")
    font = configure_matplotlib(profile)
    ANALYSIS_ROOT.mkdir(parents=True, exist_ok=False)
    charts = []
    for number, chart in enumerate(index["charts"], 1):
        result = export_chart(chart, grouped, gap_map, profile, font)
        if errors := verify_companion_record(result["manifest_record"]):
            raise ValueError(f"invalid PDF companion: {errors}")
        charts.append(result)
        if number % 12 == 0 or number == len(index["charts"]):
            print(f"[v11-pdf-review] {number}/{len(index['charts'])} vector charts checked", flush=True)
    verify_records(retained)
    figure_index = FIGURE_ROOT / "plot_manifest.json"
    write_manifest(figure_index, build_bundle({
        "schema": "publication_clean_v11_pdf_review_index_v1", "delivery_stage": "pdf_review_preflight",
        "manuscript_eligible": False, "current_publication_layer": False, "formal_vector_delivery_complete": False,
        "single_figure_count": 72, "composite_figure_count": 2, "mode_counts": {"mode_a": 36, "mode_b": 36},
        "charts": charts,
    }, [chart["manifest_record"] for chart in charts]))
    readme = ANALYSIS_ROOT / "README.md"
    readme.write_text("""# publication_clean_v11 PDF review companion

The user accepted the v11 PNG layout. This sibling supplements its 72 single
charts and two composites with genuine vector PDFs from the unchanged v11
renderer and frozen v5 values. The two composites reproduce the accepted PNG
pixels exactly when rendered again at 600 dpi. All PDFs are single-page,
contain embedded non-Type-3 fonts, and have no embedded raster images.

This is a technical review/preflight companion, not SOP vector_delivery.
The two composite final-size preflights retain the sub-2-mm math-script
violation. No typography exception is applied to those preflights and no
SOP/profile/validator threshold is changed. The accepted PNG artifacts,
their hash-bound generators/contracts, v5/current, numerical values, raw
provenance, and paper-project files remain unchanged. The package keeps
manuscript_eligible=false and formal_vector_delivery_complete=false.

One plot_manifest.json retains shared provenance and all 74 per-chart records:
the accepted PNG reference, PDF inspection, panel specifications, and all
native-size preflight failures. Native dimensions are unchanged;
paper insertion/printing dimensions require a separate placement assessment.

```powershell
python scripts/analysis/relaxtime/export_phase_guided_publication_clean_v11_pdf_review.py
```

Do not use this companion as an automatic numerical or manuscript promotion.
""", encoding="utf-8")
    write_manifest(ANALYSIS_ROOT / "manifest.json", {
        "schema": "publication_clean_v11_pdf_review_package_v1",
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(), "task_classification": "independent",
        "delivery_stage": "pdf_review_preflight", "status": "technical_pdf_review_only",
        "manuscript_eligible": False, "current_publication_layer": False, "formal_vector_delivery_complete": False,
        "solver_called": False, "canonical_data_modified": False, "new_display_values": False,
        "source_package": V11.relative(PNG_PACKAGE), "source_package_sha256": sha256_file(PNG_PACKAGE),
        "parent_manifest_sha256": package["parent_manifest_sha256"],
        "inputs": retained,
        "outputs": [input_record(readme, role="pdf_review_readme"), input_record(figure_index, role="pdf_review_index")],
        "figure_index": V11.relative(figure_index), "figure_index_sha256": sha256_file(figure_index),
        "single_figure_count": 72, "composite_figure_count": 2,
        "native_size_preflight_failure_count": sum(bool(chart["submission_preflight_violations"]) for chart in charts),
    })
    print(f"[v11-pdf-review] completed: {FIGURE_ROOT}", flush=True)


if __name__ == "__main__":
    main()
