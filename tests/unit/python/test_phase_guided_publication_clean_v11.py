from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

import pytest

from scripts.plotting.plot_manifest import sha256_file
from scripts.plotting.plot_quality import measure_figure
from scripts.plotting.plot_style import configure_matplotlib, load_profile
from scripts.plotting.validate_plot_artifact import validate_manifest


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts/analysis/relaxtime/build_phase_guided_publication_clean_v11.py"


@pytest.fixture(scope="module")
def module():
    spec = importlib.util.spec_from_file_location("publication_clean_v11_test_module", SCRIPT)
    assert spec is not None and spec.loader is not None
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


@pytest.fixture(scope="module")
def plot_inputs(module):
    points, _, _, _, gaps, _ = module.PARENT.load_v5_inputs()
    grouped, gap_map = module.V10.group_inputs(points, gaps)
    profile = load_profile("candidate_aps_v2")
    configure_matplotlib(profile)
    return grouped, gap_map, profile


def test_v11_style_contract_keeps_the_retained_v10_renderer_unchanged(module):
    assert module.X_MAJOR_TICKS == (-0.5, 0, 0.5)
    assert module.LINEAR_MINOR_SUBDIVISIONS == 2
    assert module.ENDPOINT_TITLE == "First-order"
    assert module.ENDPOINT_LABELS == {"quark": "First-order (restored)", "hadron": "First-order (broken)"}
    assert module.LOG_MAJOR_SUBS == (1, 2, 5)
    assert module.LOG_MINOR_SUBS == tuple(range(2, 10))
    assert module.LABELS["tau_u"] == r"$\tau_u\;(\mathrm{fm})$"
    assert module.V10.ENDPOINT_LABELS["quark"] == "1st-order (restored)"
    assert module.V10.FIGURE_ROOT.name == "publication_clean_v10_png_review"
    assert module.COMPOSITE_LAYOUT == module.V10.COMPOSITE_LAYOUT
    assert module.COMPOSITE_SIZE_IN == module.V10.COMPOSITE_SIZE_IN


@pytest.mark.parametrize("limits", [(0.55, 2.85), (2.4, 45), (0.3, 16), (0.001, 2000), (1.01, 1.1)])
def test_v11_log_ticks_are_readable_exact_and_within_the_view(module, limits):
    ticks = module.readable_log_ticks(*limits)
    assert 3 <= len(ticks) <= 6
    assert ticks == sorted(set(ticks))
    assert all(limits[0] <= value <= limits[1] for value in ticks)
    assert all(math.log(right / left) / math.log(limits[1] / limits[0]) >= module.LOG_MIN_LABEL_FRACTION - 1e-12
               for left, right in zip(ticks, ticks[1:]))
    labels = [module.format_log_tick(value) for value in ticks]
    assert len(set(labels)) == len(labels)
    for value, label in zip(ticks, labels):
        assert "e" not in label.lower()
        assert math.isclose(float(label), value, rel_tol=1e-9)
        assert not ("." in label and label.endswith("0"))
    assert module.format_log_tick(0.5) == "0.5"
    assert module.format_log_tick(1) == "1"
    assert module.format_log_tick(20) == "20"


@pytest.mark.parametrize("limits", [(0, 1), (-1, 2), (2, 1), (1, math.inf)])
def test_v11_log_ticks_reject_invalid_limits(module, limits):
    with pytest.raises(ValueError, match="finite, positive, and increasing"):
        module.readable_log_ticks(*limits)


def test_v11_linear_row_uses_finest_tick_spacing_and_log_row_keeps_plain_labels(module):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator

    figure, axes = plt.subplots(1, 3)
    try:
        specs = [{}, {}, {}]
        for axis, ticks in zip(axes, ([0.1, 0.15, 0.2], [0.2, 0.3, 0.4], [1, 2, 3])):
            axis.set_ylim(ticks[0], ticks[-1])
            axis.yaxis.set_major_locator(FixedLocator(ticks))
        assert module.apply_row_precision(axes, specs, "eta_over_s") == 2
        assert [item["y_tick_decimal_places"] for item in specs] == [2, 2, 2]
        assert specs[1]["y_major_tick_labels"] == ["0.20", "0.30", "0.40"]
        axes[0].set_yscale("log")
        with pytest.raises(ValueError, match="cannot mix"):
            module.apply_row_precision(axes, specs, "tau_u")
        for axis in axes:
            axis.set_yscale("log")
        assert module.apply_row_precision(axes, specs, "tau_u") is None
    finally:
        plt.close(figure)


def test_v11_composites_use_uniform_scales_row_precision_and_public_legend_alignment(module, plot_inputs):
    import matplotlib.pyplot as plt

    grouped, gap_map, profile = plot_inputs
    for name, observables in module.COMPOSITES.items():
        figure, specs = module.render_composite(observables, grouped, gap_map, profile)
        try:
            figure.canvas.draw()
            expected_scale = "log" if name.startswith("figure1_") else "linear"
            assert all(axis.get_yscale() == expected_scale for axis in figure.axes)
            assert all(item["axis_scale"] == expected_scale for item in specs)
            if expected_scale == "log":
                assert len(specs) == 12
                assert all(3 <= len(item["y_major_ticks"]) <= 6 for item in specs)
            else:
                assert [item["y_tick_decimal_places"] for item in specs] == [2] * 6 + [3] * 3
            for index in (0, 2):
                legend = figure.axes[index].get_legend()
                assert legend.get_alignment() == "left"
                assert all(text.get_fontsize() == 10.5 for text in legend.get_texts())
            legend = figure.axes[2].get_legend()
            assert legend.get_title().get_text() == "First-order"
            assert legend.get_title().get_fontsize() == 10.5
            assert legend.get_title().get_fontweight() == "normal"
            assert [text.get_text() for text in legend.get_texts()] == ["restored", "broken"]
            for axis in figure.axes:
                if axis.get_yscale() == "log":
                    assert all(axis.yaxis.get_minor_formatter()(tick) == "" for tick in axis.yaxis.get_minorticklocs())
            quality = measure_figure(figure, intended_width_inches=6.75)
            assert quality["clipped_text"] == []
            assert quality["text_overlap_pairs"] == []
            assert quality["legend_curve_overlap_count"] == 0
            assert quality["legend_landmark_overlap_count"] == 0
            assert quality["legend_in_axes_overflow_count"] == 0
        finally:
            plt.close(figure)


def test_v11_preserves_all_frozen_plot_coordinates_gaps_and_mode_b_scales(module, plot_inputs):
    import matplotlib.pyplot as plt

    grouped, gap_map, profile = plot_inputs
    figure, axis = plt.subplots()
    try:
        for mode, panel, observable in sorted({key[:3] for key in grouped}):
            axis.clear()
            spec = module.render_panel(axis, mode, panel, observable, grouped, gap_map, profile)
            expected = []
            rows_all = []
            gaps_all = []
            for series in sorted(key[3] for key in grouped if key[:3] == (mode, panel, observable)):
                rows = grouped[(mode, panel, observable, series)]
                gaps = gap_map.get((mode, panel, series), [])
                rows_all.extend(rows)
                gaps_all.extend(gaps)
                expected.extend(module.V2.split_curve_segments(rows, gaps))
            assert len(axis.lines) == len(expected)
            for line, segment in zip(axis.lines, expected):
                assert list(line.get_xdata()) == [float(row["xi"]) for row in segment]
                assert list(line.get_ydata()) == [float(row["clean_value"]) for row in segment]
            if mode == "mode_b":
                assert spec["axis_scale"] == module.V10.axis_spec_for_rows(rows_all, gaps_all)["axis_scale"]
            elif observable in module.TAU_OBSERVABLES:
                assert spec["axis_scale"] == "log"
    finally:
        plt.close(figure)


def test_v11_review_artifacts_counts_hashes_and_geometry(module):
    path = module.FIGURE_ROOT / "plot_manifest.json"
    assert path.is_file(), "build v11 with --png-review before checking the artifact package"
    index = json.loads(path.read_text(encoding="utf-8"))
    assert index["schema"] == "publication_clean_v11_png_review_figure_index_v1"
    assert index["single_figure_count"] == 72
    assert index["composite_figure_count"] == 2
    assert index["mode_counts"] == {"mode_a": 36, "mode_b": 36}
    assert len(index["charts"]) == len(list(module.FIGURE_ROOT.rglob("*.png"))) == 74
    assert not list(module.FIGURE_ROOT.rglob("*.pdf"))
    for chart in index["charts"]:
        manifest_path = ROOT / chart["manifest"]
        assert sha256_file(manifest_path) == chart["manifest_sha256"]
        assert validate_manifest(manifest_path) == []
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        assert manifest["manuscript_eligible"] is False
        assert manifest["current_publication_layer"] is False
        assert manifest["solver_called"] is False
        assert manifest["new_display_values"] is False
        assert manifest["rendering"]["legend_alignment"] == "left"
        if chart["kind"] == "single":
            assert manifest["rendering"]["legend_scope"]["endpoints"] == "actual endpoints listed in panel_specs only"
        quality = manifest["rendering"]["quality"]
        assert quality["clipped_text"] == []
        assert quality["text_overlap_pairs"] == []
        if chart["kind"] == "composite":
            assert quality["legend_curve_overlap_count"] == 0
            assert quality["legend_landmark_overlap_count"] == 0
            assert quality["legend_in_axes_overflow_count"] == 0
            assert quality["minimum_primary_capital_numeral_height_mm"] >= 2
            assert quality["minimum_capital_numeral_height_mm"] >= 1.5
        for axis in quality["tick_axes"]:
            assert axis["both_sides"] and axis["inward"]
            assert axis["minor_count"] > 0
            assert axis["major_count"] >= 3
        for output in manifest["outputs"]:
            assert output["format"] == "png"
            assert output["dpi"] == 600
            assert sha256_file(ROOT / output["path"]) == output["sha256"]


def test_v11_package_preserves_v5_v10_raw_and_current_and_records_caption_scope(module):
    package = json.loads((module.ANALYSIS_ROOT / "manifest.json").read_text(encoding="utf-8"))
    assert package["schema"] == "publication_clean_v11_png_review_package_v1"
    assert package["task_classification"] == "independent"
    assert package["author_acceptance"] is None
    assert package["manuscript_eligible"] is False
    assert package["new_display_values"] is False
    for item in [*package["inputs"], *package["outputs"]]:
        assert sha256_file(ROOT / item["path"]) == item["sha256"]
    retained = {item["path"]: item["sha256"] for item in package["inputs"] if item["role"] == "retained_v10_review_artifact"}
    assert retained == {item["path"]: item["sha256"] for item in module.retained_v10_records()}
    for name, expected in package["inherited_table_hashes"].items():
        assert sha256_file(module.V5_TABLE_ROOT / name) == expected
        assert sha256_file(module.TABLE_ROOT / name) == expected
    pointer = json.loads(module.CURRENT_POINTER.read_text(encoding="utf-8"))
    assert pointer["current_analysis_package"].endswith("publication_clean_v5")
    text = module.CAPTION_HANDOFF.read_text(encoding="utf-8")
    assert "Do not say both legends apply to all panels" in text
    assert "all on logarithmic y axes" in text
    assert "Vertical ranges differ between panels" in text
    assert "chirally restored and chirally broken branches" in text
