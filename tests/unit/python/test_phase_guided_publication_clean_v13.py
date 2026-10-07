from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pytest
from PIL import Image, ImageOps

from scripts.analysis.relaxtime import build_phase_guided_publication_clean_v13 as v13
from scripts.plotting.plot_bundle import load_chart_records
from scripts.plotting.plot_provenance import validate_hash_record
from scripts.plotting.plot_quality import measure_figure, placement_limits
from scripts.plotting.plot_style import configure_matplotlib, load_profile
from scripts.plotting.validate_plot_artifact import _check_declared_legends, validate_manifest


@pytest.fixture(scope="module")
def plot_inputs():
    points, _, _, _, gaps, _ = v13.V11.PARENT.load_v5_inputs()
    grouped, gap_map = v13.V11.V10.group_inputs(points, gaps)
    profile = load_profile("candidate_aps_v2")
    configure_matplotlib(profile)
    return grouped, gap_map, profile


def assert_frozen_vertices(figure, specs, grouped, gap_map):
    for axis, spec in zip(figure.axes, specs):
        expected = []
        mode, panel, observable = spec["mode_key"], spec["plot_panel"], spec["observable"]
        for series in sorted(key[3] for key in grouped if key[:3] == (mode, panel, observable)):
            rows = grouped[(mode, panel, observable, series)]
            expected.extend(v13.V11.V2.split_curve_segments(rows, gap_map.get((mode, panel, series), [])))
        assert len(axis.lines) == len(expected)
        for line, segment in zip(axis.lines, expected):
            assert list(line.get_xdata()) == [float(row["xi"]) for row in segment]
            assert list(line.get_ydata()) == [float(row["clean_value"]) for row in segment]


def assert_declared_layout(quality, placements):
    errors = []
    _check_declared_legends({"legend_placements": placements,
                            "legend_outside": all(p["placement"] == "outside_axes" for p in placements)},
                           quality, errors)
    assert errors == []


@pytest.mark.parametrize("name", [*v13.COMPOSITE_SIZES, "figure1_low_value_details"])
def test_final_size_geometry_and_all_original_vertices(name, plot_inputs):
    grouped, gap_map, profile = plot_inputs
    figure, specs, placements = (v13.render_detail(*plot_inputs) if name == "figure1_low_value_details"
                     else v13.render_composite(name, *plot_inputs))
    try:
        figure.set_dpi(profile.dpi)
        quality = measure_figure(figure, intended_width_inches=6.75)
        assert quality["minimum_capital_numeral_height_mm"] >= 2
        assert quality["clipped_text"] == []
        assert quality["text_overlap_pairs"] == []
        assert_declared_layout(quality, placements)
        assert [item["host_axes_index"] for item in placements] == ([0] if name == "figure1_low_value_details" else [0, 2])
        assert all(item["location"] == "upper left" and item["placement"] == "in_axes" for item in placements)
        assert quality["legend_curve_overlap_count"] == 0
        assert quality["legend_landmark_overlap_count"] == 0
        assert all(tick["inward"] and tick["both_sides"] and tick["minor_count"] for tick in quality["tick_axes"])
        assert_frozen_vertices(figure, specs, grouped, gap_map)
        if name == "figure1_low_value_details":
            assert all(axis.get_yscale() == "linear" for axis in figure.axes)
            for spec, detail in zip(specs, v13.DETAILS):
                assert spec["y_limits"] == list(detail[-1])
                assert spec["source_panel"] == detail[2]
                assert all(float(label) == pytest.approx(value) for value, label in
                           zip(spec["y_major_ticks"], spec["y_major_tick_labels"]))
        else:
            assert list(figure.get_size_inches()) == list(v13.V11.COMPOSITE_SIZE_IN[name])
            expected_scale = "log" if name.startswith("figure1") else "linear"
            assert all(axis.get_yscale() == expected_scale for axis in figure.axes)
    finally:
        plt.close(figure)


def test_all_single_views_keep_values_and_expose_single_column_failure(plot_inputs):
    grouped, gap_map, profile = plot_inputs
    for mode, panel, observable in sorted({key[:3] for key in grouped}):
        figure, specs, placements = v13.render_single(mode, panel, observable, *plot_inputs)
        try:
            figure.set_dpi(profile.dpi)
            quality = measure_figure(figure, intended_width_inches=6.75)
            assert quality["minimum_capital_numeral_height_mm"] >= 2
            assert quality["clipped_text"] == []
            assert quality["text_overlap_pairs"] == []
            assert quality["legend_curve_overlap_count"] == 0
            assert quality["legend_landmark_overlap_count"] == 0
            assert_declared_layout(quality, placements)
            parameter = next(item for item in placements if item["role"] == "parameters")
            if mode == "mode_a":
                assert parameter["title"] == r"$\alpha_T$"
                assert parameter["labels"] == ["1.0", "1.1", "1.2"]
                assert parameter["series_keys"] == ["alpha1.0", "alpha1.1", "alpha1.2"]
            else:
                assert parameter["title"] == r"$\mu_B\;(\mathrm{MeV})$"
                assert parameter["labels"] == ["0", "450", "900"]
                assert parameter["series_keys"] == ["muB0.0", "muB450.0", "muB900.0"]
            for key in (item for item in placements if item["role"] == "endpoints"):
                assert key["labels"] == ["restored", "broken"]
                assert key["applies_to"] == [{"mode_key": mode, "plot_panel": panel, "plot_series": series}
                                              for series in sorted({p["series"] for p in specs[0]["endpoints"]})]
            assert placement_limits(quality, profile, [])["single_column_reuse_qualified"] is False
            assert_frozen_vertices(figure, specs, grouped, gap_map)
        finally:
            plt.close(figure)


def test_artifact_bundle_has_verified_color_gray_pairs_and_unchanged_sources():
    root = v13.ROOT
    path = v13.FIGURE_ROOT / "plot_manifest.json"
    assert path.is_file(), "generate the v13 case before artifact validation"
    index, pairs = load_chart_records(path, root=root)
    assert (index["single_figure_count"], index["composite_figure_count"], index["detail_figure_count"]) == (72, 2, 1)
    assert len(pairs) == 75
    assert len(list(v13.FIGURE_ROOT.rglob("*.png"))) == 150
    assert list(v13.FIGURE_ROOT.rglob("*.json")) == [path]
    assert not list(v13.FIGURE_ROOT.rglob("*.pdf"))
    assert validate_manifest(path, repo_root=root) == []
    for chart, record in pairs:
        assert record["manuscript_eligible"] is False
        assert record["rendering"]["typography_exception"] is None
        assert record["rendering"]["legend_policy"] == "declared_geometry_checked"
        assert_declared_layout(record["rendering"]["quality"], record["rendering"]["legend_placements"])
        assert record["rendering"]["placement_limits"]["measured_width_qualified"]
        assert record["rendering"]["grayscale_review"]["status"] == "author_review_required"
        color, gray = record["outputs"]
        assert gray["source_color_sha256"] == color["sha256"]
        with Image.open(root / color["path"]) as a, Image.open(root / gray["path"]) as b:
            assert a.size == b.size
            assert b.mode == "L"
            assert ImageOps.grayscale(a.convert("RGB")).tobytes() == b.tobytes()
    package = json.loads((v13.ANALYSIS_ROOT / "manifest.json").read_text(encoding="utf-8"))
    for record in package["inputs"]:
        assert validate_hash_record(record, root=root, label="inputs") == []
    for record in package["outputs"]:
        assert validate_hash_record(record, root=root, label="outputs", allow_historical=False) == []
    assert package["author_acceptance"] is None
    assert package["canonical_data_modified"] is False
    assert package["new_display_values"] is False
    assert package["solver_called"] is False
    paths = {record["path"] for record in package["inputs"]}
    assert "docs/guides/sop/workflows/figure_production.md" in paths
    assert "config/plotting/candidate_aps_v2.toml" in paths
    assert "scripts/plotting/plot_provenance.py" in paths
    assert "docs/analysis/relaxtime/phase_guided_transport/literature_figure_style_audit_20261005/artifact_manifest.json" in paths
    assert "docs/analysis/relaxtime/phase_guided_transport/publication_clean_v12_code_snapshot_v1.zip" in paths


def test_rebuild_refuses_to_overwrite_artifact(monkeypatch, tmp_path):
    monkeypatch.setattr(v13, "FIGURE_ROOT", tmp_path)
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        v13.build_review()
