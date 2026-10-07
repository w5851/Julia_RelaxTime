from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pytest
from PIL import Image, ImageOps

from scripts.analysis.relaxtime import build_phase_guided_publication_clean_v12 as v12
from scripts.plotting.plot_bundle import load_chart_records
from scripts.plotting.plot_provenance import validate_hash_record
from scripts.plotting.plot_quality import measure_figure, placement_limits
from scripts.plotting.plot_style import configure_matplotlib, load_profile
from scripts.plotting.validate_plot_artifact import validate_manifest


@pytest.fixture(scope="module")
def plot_inputs():
    points, _, _, _, gaps, _ = v12.V11.PARENT.load_v5_inputs()
    grouped, gap_map = v12.V11.V10.group_inputs(points, gaps)
    profile = load_profile("candidate_aps_v2")
    configure_matplotlib(profile)
    return grouped, gap_map, profile


def assert_frozen_vertices(figure, specs, grouped, gap_map):
    for axis, spec in zip(figure.axes, specs):
        expected = []
        mode, panel, observable = spec["mode_key"], spec["plot_panel"], spec["observable"]
        for series in sorted(key[3] for key in grouped if key[:3] == (mode, panel, observable)):
            rows = grouped[(mode, panel, observable, series)]
            expected.extend(v12.V11.V2.split_curve_segments(rows, gap_map.get((mode, panel, series), [])))
        assert len(axis.lines) == len(expected)
        for line, segment in zip(axis.lines, expected):
            assert list(line.get_xdata()) == [float(row["xi"]) for row in segment]
            assert list(line.get_ydata()) == [float(row["clean_value"]) for row in segment]


@pytest.mark.parametrize("name", [*v12.COMPOSITE_SIZES, "figure1_low_value_details"])
def test_final_size_geometry_and_all_original_vertices(name, plot_inputs):
    grouped, gap_map, profile = plot_inputs
    figure, specs = (v12.render_detail(*plot_inputs) if name == "figure1_low_value_details"
                     else v12.render_composite(name, *plot_inputs))
    try:
        quality = measure_figure(figure, intended_width_inches=6.75)
        assert quality["minimum_capital_numeral_height_mm"] >= 2
        assert quality["clipped_text"] == []
        assert quality["text_overlap_pairs"] == []
        assert quality["legend_axes_overlap_count"] == 0
        assert quality["legend_curve_overlap_count"] == 0
        assert quality["legend_landmark_overlap_count"] == 0
        assert all(tick["inward"] and tick["both_sides"] and tick["minor_count"] for tick in quality["tick_axes"])
        assert_frozen_vertices(figure, specs, grouped, gap_map)
        if name == "figure1_low_value_details":
            assert all(axis.get_yscale() == "linear" for axis in figure.axes)
            for spec, detail in zip(specs, v12.DETAILS):
                assert spec["y_limits"] == list(detail[-1])
                assert spec["source_panel"] == detail[2]
                assert all(float(label) == pytest.approx(value) for value, label in
                           zip(spec["y_major_ticks"], spec["y_major_tick_labels"]))
        else:
            assert list(figure.get_size_inches()) == list(v12.V11.COMPOSITE_SIZE_IN[name])
            expected_scale = "log" if name.startswith("figure1") else "linear"
            assert all(axis.get_yscale() == expected_scale for axis in figure.axes)
    finally:
        plt.close(figure)


def test_all_single_views_keep_values_and_expose_single_column_failure(plot_inputs):
    grouped, gap_map, profile = plot_inputs
    for mode, panel, observable in sorted({key[:3] for key in grouped}):
        figure, specs = v12.render_single(mode, panel, observable, *plot_inputs)
        try:
            quality = measure_figure(figure, intended_width_inches=6.75)
            assert quality["minimum_capital_numeral_height_mm"] >= 2
            assert quality["clipped_text"] == []
            assert quality["text_overlap_pairs"] == []
            assert quality["legend_curve_overlap_count"] == 0
            assert quality["legend_landmark_overlap_count"] == 0
            assert quality["legend_axes_overlap_count"] == 0
            assert placement_limits(quality, profile, [])["single_column_reuse_qualified"] is False
            assert_frozen_vertices(figure, specs, grouped, gap_map)
        finally:
            plt.close(figure)


def test_artifact_bundle_has_verified_color_gray_pairs_and_unchanged_sources():
    root = v12.ROOT
    path = v12.FIGURE_ROOT / "plot_manifest.json"
    assert path.is_file(), "generate the v12 case before artifact validation"
    index, pairs = load_chart_records(path, root=root)
    assert (index["single_figure_count"], index["composite_figure_count"], index["detail_figure_count"]) == (72, 2, 1)
    assert len(pairs) == 75
    assert len(list(v12.FIGURE_ROOT.rglob("*.png"))) == 150
    assert list(v12.FIGURE_ROOT.rglob("*.json")) == [path]
    assert not list(v12.FIGURE_ROOT.rglob("*.pdf"))
    assert validate_manifest(path, repo_root=root) == []
    for chart, record in pairs:
        assert record["manuscript_eligible"] is False
        assert record["rendering"]["typography_exception"] is None
        assert record["rendering"]["placement_limits"]["measured_width_qualified"]
        assert record["rendering"]["grayscale_review"]["status"] == "author_review_required"
        color, gray = record["outputs"]
        assert gray["source_color_sha256"] == color["sha256"]
        with Image.open(root / color["path"]) as a, Image.open(root / gray["path"]) as b:
            assert a.size == b.size
            assert b.mode == "L"
            assert ImageOps.grayscale(a.convert("RGB")).tobytes() == b.tobytes()
    package = json.loads((v12.ANALYSIS_ROOT / "manifest.json").read_text(encoding="utf-8"))
    for record in package["inputs"]:
        assert validate_hash_record(record, root=root, label="inputs") == []
    for record in package["outputs"]:
        assert validate_hash_record(record, root=root, label="outputs", allow_historical=False) == []
    assert package["author_acceptance"] is None
    assert package["canonical_data_modified"] is False
    assert package["new_display_values"] is False
    assert package["solver_called"] is False


def test_rebuild_refuses_to_overwrite_artifact(monkeypatch, tmp_path):
    monkeypatch.setattr(v12, "FIGURE_ROOT", tmp_path)
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        v12.build_review()
