from __future__ import annotations

import copy
import json
from pathlib import Path
import shutil

import matplotlib
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scripts.plotting.plot_manifest import (
    build_manifest, generator_record, input_record, output_record, write_manifest,
)
from scripts.plotting.plot_quality import export_figure, inspect_export, measure_figure, visible_texts
from scripts.plotting.plot_style import configure_axis_ticks, configure_matplotlib, load_profile
from scripts.plotting.validate_plot_artifact import validate_manifest


@pytest.fixture(scope="module")
def exported_fixture(tmp_path_factory):
    if not all(shutil.which(tool) for tool in ("pdfinfo", "pdffonts", "pdfimages")):
        pytest.skip("Poppler is required for actual PDF artifact checks")
    folder = tmp_path_factory.mktemp("plotting_quality")
    source = folder / "input.csv"
    source.write_text("x,y\n-0.5,1\n0,2\n0.5,3\n", encoding="utf-8")
    profile = load_profile("strict_aps_v2")
    with matplotlib.rc_context():
        configure_matplotlib(profile)
        figure, axis = plt.subplots(figsize=(6.75, 4.6))
        figure.subplots_adjust(left=0.16, right=0.95, top=0.73, bottom=0.18)
        axis.plot([-0.5, 0, 0.5], [1, 2, 3], label=r"$\alpha_T = 1.0$")
        axis.set_xlabel(r"$\xi$")
        axis.set_ylabel(r"$\tau_u\;(\mathrm{fm})$")
        configure_axis_ticks(axis, profile)
        handles, labels = axis.get_legend_handles_labels()
        figure.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.97))
        outputs, quality = export_figure(figure, folder / "figure", profile)
        plt.close(figure)
    manifest = build_manifest(
        asset_id="quality.test", figure_family="test", case_slug="test", figure_mode="strict",
        semantic_status="confirmed", style_profile=profile.profile_id, publication_scope="main_text_candidate",
        generator=generator_record(Path(__file__), command="pytest"), inputs=[input_record(source, role="fixture")],
        axes=[{"field": "tau_u", "source_unit": "fm", "display_unit": "fm", "label": r"$\tau_u\;(\mathrm{fm})$"}],
        series=[{"series_id": "fixture", "state": "model_support", "support_rule": "frozen rows", "mask_rule": "none"}],
        outputs=outputs, selection_rule="fixture", interpolation_policy="none", connector_policy="forbidden",
        missing_value_policy="reject", validation={"finite": True, "duplicate_keys": True, "support": True, "strict_gate": True},
        rendering={"column": "double_column", "figure_size_inches": quality["figure_size_inches"], "quality": quality,
                   "color_route": "color_print_and_online"},
    )
    return manifest


def test_pdf_png_actual_export_passes_final_size_gate(exported_fixture, tmp_path):
    path = tmp_path / "plot_manifest.json"
    write_manifest(path, exported_fixture)
    assert validate_manifest(path) == []
    pdf = next(item for item in exported_fixture["outputs"] if item["format"] == "pdf")
    assert pdf["inspection"]["raster_image_count"] == 0
    assert pdf["inspection"]["fonts_embedded"] is True


def test_measure_figure_detects_curve_under_in_axes_legend():
    profile = load_profile("candidate_aps_v2")
    with matplotlib.rc_context():
        configure_matplotlib(profile)
        figure, axis = plt.subplots(figsize=(6.75, 4.6))
        axis.plot([-0.5, 0.0, 0.5], [1.75, 1.75, 1.75], label="curve")
        axis.set_xlim(-0.5, 0.5)
        axis.set_ylim(0.0, 3.5)
        axis.legend(loc="center")
        quality = measure_figure(figure, intended_width_inches=6.75)
        plt.close(figure)
    assert quality["legend_curve_overlap_count"] == 1
    assert quality["legend_curve_overlaps"][0]["label"] == "curve"


def test_legend_geometry_preserves_nan_gaps_and_detects_endpoint_markers():
    with matplotlib.rc_context():
        figure, axis = plt.subplots(figsize=(6.75, 4.6))
        axis.plot([-0.5, -0.25, float("nan"), 0.25, 0.5], [1, 1, float("nan"), 1, 1], label="split")
        axis.set_xlim(-0.5, 0.5)
        axis.set_ylim(0, 2)
        legend = axis.legend(loc="center", title="1st-order")
        quality = measure_figure(figure, intended_width_inches=6.75)
        assert quality["legend_curve_overlap_count"] == 0
        assert legend.get_title() in visible_texts(figure)
        axis.scatter([0.0], [1.0], s=25, facecolor="white", edgecolor="blue")
        quality = measure_figure(figure, intended_width_inches=6.75)
        assert quality["legend_landmark_overlap_count"] == 1
        plt.close(figure)


@pytest.mark.parametrize("violation", ["curve", "landmark", "missing_evidence", "overflow", "eligible", "vector_delivery"])
def test_geometry_and_typography_review_exceptions_cannot_bypass_gates(exported_fixture, tmp_path, violation):
    payload = copy.deepcopy(exported_fixture)
    payload.update(figure_mode="audit", publication_scope="internal_review",
                   manuscript_eligible=False, current_publication_layer=False)
    payload["outputs"] = [record for record in payload["outputs"] if record["format"] == "png"]
    rendering = payload["rendering"]
    rendering.update(delivery_stage="png_review", vector_delivery_pending=True,
                     legend_policy="shared_in_panel_reviewed_geometry_checked",
                     typography_exception="dense_composite_review_compact_typography")
    quality = rendering["quality"]
    quality["legend_axes_overlap_count"] = 1
    if violation in {"curve", "landmark"}:
        quality[f"legend_{violation}_overlap_count"] = 1
        quality[f"legend_{violation}_overlaps"] = [{"axis_index": 0}]
    elif violation == "missing_evidence":
        del quality["legend_curve_overlap_count"]
    elif violation == "overflow":
        quality["legend_in_axes_overflow_count"] = 1
    elif violation == "eligible":
        payload["manuscript_eligible"] = True
    else:
        rendering["delivery_stage"] = "vector_delivery"
    path = tmp_path / "bad_review.json"
    write_manifest(path, payload)
    assert validate_manifest(path)


def test_png_review_stage_allows_only_png_and_marks_vector_pending(exported_fixture, tmp_path):
    payload = copy.deepcopy(exported_fixture)
    payload["figure_mode"] = "audit"
    payload["publication_scope"] = "internal_review"
    payload["outputs"] = [item for item in payload["outputs"] if item["format"] == "png"]
    payload["manuscript_eligible"] = False
    payload["current_publication_layer"] = False
    payload["rendering"]["delivery_stage"] = "png_review"
    payload["rendering"]["vector_delivery_pending"] = True
    payload["rendering"]["output_formats"] = ["png"]
    path = tmp_path / "png_review_manifest.json"
    write_manifest(path, payload)
    assert validate_manifest(path) == []


@pytest.mark.parametrize("violation,expected", [
    ("no_pdf", "missing required formats"),
    ("square_unit", "parentheses"),
    ("strict_connector", "connector_policy"),
    ("strict_interpolation", "interpolation_policy"),
    ("small_height", "below 2 mm"),
    ("missing_geometry", "measurement evidence"),
    ("outward_ticks", "inward major/minor"),
    ("legend_overlap", "legend must be outside"),
    ("size_mismatch", "physical size differs"),
    ("print_route", "requires PS/EPS"),
    ("undecided_strict", "select its color/print route"),
])
def test_v2_rejects_invalid_delivery(exported_fixture, tmp_path, violation, expected):
    payload = copy.deepcopy(exported_fixture)
    quality = payload["rendering"]["quality"]
    if violation == "no_pdf":
        payload["outputs"] = [item for item in payload["outputs"] if item["format"] != "pdf"]
    elif violation == "square_unit":
        payload["axes"][0]["label"] = r"$\tau_u\;[\mathrm{fm}]$"
    elif violation == "strict_connector":
        payload["connector_policy"] = "allowed"
    elif violation == "strict_interpolation":
        payload["interpolation_policy"] = "linear"
    elif violation == "small_height":
        quality["minimum_capital_numeral_height_mm"] = 1.9
        quality["smallest_glyphs"][0]["final_height_mm"] = 1.9
    elif violation == "missing_geometry":
        del quality["clipped_text"]
    elif violation == "outward_ticks":
        quality["tick_axes"][0]["inward"] = False
    elif violation == "legend_overlap":
        quality["legend_axes_overlap_count"] = 1
    elif violation == "size_mismatch":
        quality["figure_size_inches"][1] += 0.1
        payload["rendering"]["figure_size_inches"][1] += 0.1
    elif violation == "print_route":
        payload["rendering"]["color_route"] = "color_online_grayscale_print"
    elif violation == "undecided_strict":
        payload["rendering"]["color_route"] = "undecided_review"
    path = tmp_path / "plot_manifest.json"
    write_manifest(path, payload)
    assert any(expected in error for error in validate_manifest(path))


def test_raster_wrapped_pdf_is_rejected_even_when_declared_vector(exported_fixture, tmp_path):
    payload = copy.deepcopy(exported_fixture)
    with matplotlib.rc_context():
        configure_matplotlib(load_profile("strict_aps_v2"))
        figure = plt.figure(figsize=(6.75, 4.6))
        axis = figure.add_axes([0, 0, 1, 1])
        axis.imshow([[0, 1], [1, 0]], interpolation="nearest")
        axis.set_axis_off()
        pdf_path = tmp_path / "raster.pdf"
        figure.savefig(pdf_path, bbox_inches=None)
        plt.close(figure)
    pdf = output_record(pdf_path, fmt="pdf", dpi=None, vector=True)
    pdf["inspection"] = inspect_export(pdf_path)
    payload["outputs"] = [pdf, next(item for item in payload["outputs"] if item["format"] == "png")]
    manifest = tmp_path / "plot_manifest.json"
    write_manifest(manifest, payload)
    assert any("not embedded raster images" in error for error in validate_manifest(manifest))


def test_nominal_font_size_does_not_certify_scaled_or_scripted_glyphs():
    with matplotlib.rc_context():
        configure_matplotlib(load_profile("strict_aps_v2"))
        figure, axis = plt.subplots(figsize=(6.75, 4.6))
        figure.subplots_adjust(top=0.75)
        axis.plot([0, 1], [0, 1])
        axis.set_xlabel(r"$\alpha_T$")
        configure_axis_ticks(axis, load_profile("strict_aps_v2"))
        native = measure_figure(figure, intended_width_inches=6.75)
        reduced = measure_figure(figure, intended_width_inches=6.75 * 0.8)
        plt.close(figure)
    assert native["minimum_capital_numeral_height_mm"] >= 2
    assert reduced["minimum_capital_numeral_height_mm"] < 2
    assert native["smallest_glyphs"][0]["glyph"] == "T"
    assert reduced["minimum_capital_numeral_height_mm"] == pytest.approx(native["minimum_capital_numeral_height_mm"] * 0.8)
