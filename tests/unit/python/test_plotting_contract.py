from __future__ import annotations

import json
from pathlib import Path

from PIL import Image

from scripts.plotting.plot_manifest import (
    PROJECT_ROOT,
    build_manifest,
    generator_record,
    input_record,
    output_record,
    write_manifest,
)
from scripts.plotting.plot_style import figure_size_inches, load_profile
from scripts.plotting.validate_plot_artifact import validate_manifest


def _make_manifest(tmp_path: Path, *, mode: str = "strict", scope: str = "main_text_candidate") -> Path:
    input_path = tmp_path / "input.csv"
    input_path.write_text("x,y\n1,2\n", encoding="utf-8")
    png_path = tmp_path / "figure.png"
    Image.new("RGB", (20, 20), "white").save(png_path, dpi=(600, 600))
    svg_path = tmp_path / "figure.svg"
    svg_path.write_text(
        '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20"></svg>\n',
        encoding="utf-8",
    )
    manifest = build_manifest(
        asset_id="test.plotting.contract",
        figure_family="test",
        case_slug="test_case",
        figure_mode=mode,
        semantic_status="estimated_midpoint" if mode == "estimated_midpoint" else "confirmed",
        style_profile="candidate_origin_like_v1" if mode == "estimated_midpoint" else "strict_origin_like_v1",
        publication_scope=scope,
        generator=generator_record(
            PROJECT_ROOT / "scripts" / "plotting" / "plot_manifest.py",
            command="pytest",
        ),
        inputs=[input_record(input_path, role="fixture", schema="fixture_v1", units={"x": "a.u.", "y": "a.u."})],
        axes=[{"field": "x", "source_unit": "a.u.", "display_unit": "a.u.", "transform": "identity"}],
        series=[
            {
                "series_id": "fixture_series",
                "state": "model_support" if mode != "estimated_midpoint" else "estimated_midpoint",
                "support_rule": "fixture rows",
                "mask_rule": "none",
            }
        ],
        outputs=[
            output_record(png_path, fmt="png", dpi=600, vector=False),
            output_record(svg_path, fmt="svg", dpi=600, vector=True),
        ],
        selection_rule="fixture rows",
        interpolation_policy="none",
        connector_policy="forbidden",
        missing_value_policy="reject nonfinite",
        validation={"finite": True, "duplicate_keys": True, "support": True, "strict_gate": mode == "strict"},
        rendering={"column": "single_column", "figure_size_inches": [3.375, 2.5], "bbox_inches": "tight", "pad_inches": 0.08},
    )
    manifest_path = tmp_path / "plot_manifest.json"
    write_manifest(manifest_path, manifest)
    return manifest_path


def test_profiles_freeze_current_aps_like_baseline():
    audit = load_profile("audit_v1")
    candidate = load_profile("candidate_origin_like_v1")
    strict = load_profile("strict_origin_like_v1")
    assert figure_size_inches(audit, "single_column") == (3.375, 2.5)
    assert figure_size_inches(audit, "double_column") == (6.75, 4.6)
    assert candidate.formats == ("png", "svg")
    assert strict.formats == ("png", "svg")
    assert strict.optional_formats == ("pdf",)
    assert strict.dpi == 600
    assert audit.draw_support_points is True
    assert candidate.draw_support_points is False
    assert strict.draw_support_points is False
    assert strict.support_visual_policy == "line_first_landmark_only"
    assert strict.legend_policy == "dense_aware_best_then_external"
    assert strict.max_in_axes_legend_entries == 4
    assert strict.allow_external_legend is True


def test_strict_png_svg_manifest_passes(tmp_path):
    manifest_path = _make_manifest(tmp_path)
    assert validate_manifest(manifest_path) == []


def test_estimated_midpoint_requires_internal_scope(tmp_path):
    manifest_path = _make_manifest(tmp_path, mode="estimated_midpoint", scope="main_text_candidate")
    errors = validate_manifest(manifest_path)
    assert any("supplement_or_internal_review" in error for error in errors)


def test_strict_rejects_estimated_or_connector_series(tmp_path):
    manifest_path = _make_manifest(tmp_path)
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload["series"][0]["state"] = "plotting_connector"
    manifest_path.write_text(json.dumps(payload), encoding="utf-8")
    errors = validate_manifest(manifest_path)
    assert any("forbidden state" in error for error in errors)


def test_input_hash_mismatch_is_rejected(tmp_path):
    manifest_path = _make_manifest(tmp_path)
    input_path = tmp_path / "input.csv"
    input_path.write_text("x,y\n1,999\n", encoding="utf-8")
    errors = validate_manifest(manifest_path)
    assert any("sha256 mismatch" in error for error in errors)
