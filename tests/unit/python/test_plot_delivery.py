from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest
from PIL import Image

from scripts.plotting.plot_accessibility import accessibility_record, contrast_ratio, validate_accessibility
from scripts.plotting.plot_delivery import validate_png_delivery
from scripts.plotting.plot_manifest import input_record, output_record
from scripts.plotting.plot_style import load_profile


def test_known_wcag_contrasts_and_opaque_color_requirement():
    assert contrast_ratio("black", "white") == 21
    assert contrast_ratio("#E69F00", "white") == pytest.approx(2.2522876008)
    assert contrast_ratio("#E69F00", "#D9D7CA") == pytest.approx(1.5576576321)
    with pytest.raises(ValueError, match="opaque"):
        contrast_ratio("#00000080", "white")
    for name in ("candidate_aps_v2", "strict_aps_v2"):
        assert all(contrast_ratio(color, "white") >= 3 for color in load_profile(name).colors)


def test_accessibility_rejects_color_only_and_fabricated_contrast():
    record = accessibility_record([
        {"value": -0.5, "color": "#0072B2", "noncolor_encoding": "circle"},
        {"value": 0.5, "color": "#A86F00", "noncolor_encoding": "triangle"}],
        backgrounds=["white"], background_scope="opaque lines; no fill")
    assert validate_accessibility(record) == []
    record["parameters"][1]["noncolor_encoding"] = "circle"
    assert any("distinct noncolor" in issue for issue in validate_accessibility(record))
    record["parameters"][1].update(color="#E69F00", contrasts=[10])
    errors = validate_accessibility(record)
    assert any("disagree" in issue for issue in errors)
    assert any("below 3:1" in issue for issue in errors)


@pytest.fixture
def delivery(tmp_path):
    review_png = tmp_path / "review.png"
    delivered_png = tmp_path / "delivered.png"
    Image.new("RGB", (12, 12), "white").save(review_png)
    delivered_png.write_bytes(review_png.read_bytes())
    source = tmp_path / "input.csv"
    source.write_text("x,y\n1,2\n", encoding="utf-8")
    accepted = {"style_profile": "candidate_aps_v2", "axes": [{"field": "x"}],
        "series": [{"state": "estimated_density", "value": 2.0}], "selection_rule": "frozen",
        "interpolation_policy": "none", "connector_policy": "explicit display closure", "missing_value_policy": "preserve gaps",
        "derived_display_geometry": [{"start": [1, 2], "end": [2, 3]}], "caption_parameters": {"estimate": "density only"},
        "frozen_render_signature": {"input": "fixture"}, "generator": {"sha256": "fixture-generator"},
        "inputs": [input_record(source, role="calculation_result", root=tmp_path)],
        "outputs": [output_record(review_png, fmt="png", dpi=600, vector=False, root=tmp_path)],
        "manuscript_eligible": False, "current_publication_layer": False,
        "rendering": {"delivery_stage": "png_review", "vector_delivery_pending": True,
            "figure_size_inches": [6.75, 4.6], "column": "double_column", "legend_policy": "declared_geometry_checked",
            "legend_placements": [], "parameter_encoding": "shape", "coexistence_fill": "none", "accessibility": {"fixture": True}}}
    png_manifest = tmp_path / "accepted_manifest.json"
    png_manifest.write_text(json.dumps(accepted), encoding="utf-8")
    receipt = tmp_path / "acceptance.json"
    response = {"schema": "plot_png_acceptance_v1", "status": "author_accepted", "pdf_export_authorized": True,
        "author_instruction": "Synthetic test fixture; not a real author approval.",
        "accepted_manifest": input_record(png_manifest, role="accepted_png_manifest", root=tmp_path),
        "accepted_png": input_record(review_png, role="accepted_png", root=tmp_path)}
    receipt.write_text(json.dumps(response), encoding="utf-8")
    result = copy.deepcopy(accepted)
    result["outputs"] = [output_record(delivered_png, fmt="png", dpi=600, vector=False, root=tmp_path)]
    result["rendering"].update(delivery_stage="vector_delivery", vector_delivery_pending=False)
    result["author_review"] = {"status": "author_accepted_png", "record": input_record(receipt, role="author_png_acceptance", root=tmp_path)}
    return result, receipt, tmp_path


def test_hash_bound_acceptance_and_relocated_identical_png_pass(delivery):
    manifest, _, root = delivery
    assert validate_png_delivery(manifest, root=root) == []


def test_missing_acceptance_is_rejected_even_with_old_pass_flags(delivery):
    manifest, _, root = delivery
    manifest.pop("author_review")
    manifest["vector_delivery_verification"] = {"status": "passed", "png_byte_identity": True}
    assert any("author_accepted_png" in error for error in validate_png_delivery(manifest, root=root))


def test_pending_author_cannot_export_vectors(delivery):
    manifest, receipt, root = delivery
    response = json.loads(receipt.read_text())
    response["status"] = "pending"
    receipt.write_text(json.dumps(response))
    manifest["author_review"]["record"] = input_record(receipt, role="author_png_acceptance", root=root)
    assert any("actual author instruction" in error for error in validate_png_delivery(manifest, root=root))


@pytest.mark.parametrize("field", ["axes", "series", "derived_display_geometry", "caption_parameters", "frozen_render_signature"])
def test_accepted_display_semantics_cannot_drift(delivery, field):
    manifest, _, root = delivery
    manifest[field] = {"changed": True}
    assert any(f"semantic mismatch: {field}" in error for error in validate_png_delivery(manifest, root=root))


def test_palette_or_legend_changes_need_a_new_png_review(delivery):
    manifest, _, root = delivery
    manifest["rendering"]["parameter_encoding"] = "another shape mapping"
    assert any("rendering mismatch" in error for error in validate_png_delivery(manifest, root=root))


def test_input_removal_and_png_substitution_are_rejected(delivery):
    manifest, _, root = delivery
    manifest["inputs"] = []
    manifest["outputs"][0]["sha256"] = "0" * 64
    errors = validate_png_delivery(manifest, root=root)
    assert any("frozen input" in error for error in errors)
    assert any("differs from" in error for error in errors)
