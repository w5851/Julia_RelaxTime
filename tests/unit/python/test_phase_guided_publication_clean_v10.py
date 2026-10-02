from __future__ import annotations

import importlib.util
from functools import lru_cache
import json
from pathlib import Path

from scripts.plotting.plot_manifest import sha256_file
from scripts.plotting.plot_style import configure_matplotlib, load_profile
from scripts.plotting.validate_plot_artifact import validate_manifest
from scripts.plotting.plot_provenance import validate_hash_record


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "relaxtime" / "build_phase_guided_publication_clean_v10.py"


def _module():
    spec = importlib.util.spec_from_file_location("publication_clean_v10_test_module", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@lru_cache(maxsize=1)
def _stage_module():
    path = ROOT / "scripts/analysis/relaxtime/formalize_phase_guided_publication_clean_v11_stage.py"
    spec = importlib.util.spec_from_file_location("v10_snapshot_stage_contract", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _assert_only_recorded_historical_drift(manifest_path):
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    stage = _stage_module()
    drift_paths = stage.HISTORICAL_CONTRACT_PATHS
    expected = []
    for index, record in enumerate(manifest["inputs"]):
        if record["path"] in drift_paths:
            assert sha256_file(ROOT / record["path"]) != record["sha256"]
            expected.extend([
                f"inputs[{index}].bytes mismatch for {record['path']}",
                f"inputs[{index}].sha256 mismatch for {record['path']}",
            ])
    assert len(expected) == 4
    assert set(validate_manifest(manifest_path, code_ref=stage.FROZEN_CODE_REF)) == set(expected)


def test_v10_style_contract_is_png_review_only():
    module = _module()
    assert module.X_MAJOR_TICKS == (-0.5, 0.0, 0.5)
    assert module.LINEAR_MINOR_SUBDIVISIONS == 2
    assert module.FORCED_LOG_PANEL == "muB900.0"
    assert module.ENDPOINT_LABELS == {
        "quark": "1st-order (restored)",
        "hadron": "1st-order (broken)",
    }
    assert module.LABELS["tau_u"] == r"$\tau_u\;(\mathrm{fm})$"


def test_v10_axis_rule_can_force_log_without_changing_values():
    module = _module()
    rows = [{"clean_value": "1.0"}, {"clean_value": "2.0"}]
    spec = module.axis_spec_for_rows(rows, [], force_log=True)
    assert spec["axis_scale"] == "log"
    assert spec["axis_scale_reason"] == "forced_log_for_complete_muB900_relaxation_column"
    assert spec["data_min"] == 1.0
    assert spec["data_max"] == 2.0


def test_v10_review_package_counts_hashes_and_layout_contract():
    module = _module()
    index_path = module.FIGURE_ROOT / "plot_manifest.json"
    assert index_path.is_file(), "build v10 with --png-review before running focused artifact checks"
    index = json.loads(index_path.read_text(encoding="utf-8"))
    assert index["schema"] == "publication_clean_v10_png_review_figure_index_v1"
    assert index["single_figure_count"] == 72
    assert index["mode_counts"] == {"mode_a": 36, "mode_b": 36}
    assert index["composite_figure_count"] == 2
    assert index["manuscript_eligible"] is False
    assert index["current_publication_layer"] is False
    assert len(index["charts"]) == 74
    assert len(list(module.FIGURE_ROOT.rglob("*.png"))) == 74
    assert len(list(module.FIGURE_ROOT.rglob("*.pdf"))) == 0

    for chart in index["charts"]:
        manifest_path = ROOT / chart["manifest"]
        assert sha256_file(manifest_path) == chart["manifest_sha256"]
        _assert_only_recorded_historical_drift(manifest_path)
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        assert manifest["manuscript_eligible"] is False
        assert manifest["current_publication_layer"] is False
        assert manifest["canonical_data_modified"] is False
        assert manifest["solver_called"] is False
        assert manifest["rendering"]["legend_outside"] is False
        assert {output["format"] for output in manifest["outputs"]} == {"png"}
        assert manifest["rendering"]["tick_policy"]["linear_minor_subdivisions"] == 2
        assert manifest["rendering"]["x_tick_label_policy"] == "three labels at -0.5, 0, and 0.5"
        for output in manifest["outputs"]:
            assert output["dpi"] == 600
            assert sha256_file(ROOT / output["path"]) == output["sha256"]
        quality = manifest["rendering"]["quality"]
        assert quality["clipped_text"] == []
        assert quality["text_overlap_pairs"] == []
        assert quality["legend_axes_overlap_count"] > 0
        if chart["kind"] == "composite":
            assert quality["legend_curve_overlap_count"] == 0
            assert quality["legend_landmark_overlap_count"] == 0
            assert quality["legend_in_axes_overflow_count"] == 0
            assert quality["minimum_primary_capital_numeral_height_mm"] >= 2
            assert quality["minimum_capital_numeral_height_mm"] >= 1.5
        for tick in quality["tick_axes"]:
            assert tick["inward"] is True
            assert tick["both_sides"] is True
            assert tick["minor_count"] > 0
            assert tick["major_count"] > 0
            if tick["axis"] == "x":
                assert tick["major_count"] == 3
            elif tick["scale"] == "linear":
                assert tick["major_count"] >= (3 if chart["kind"] == "composite" else 4)

    composite_manifests = [
        json.loads((ROOT / chart["manifest"]).read_text(encoding="utf-8"))
        for chart in index["charts"]
        if chart["kind"] == "composite"
    ]
    figure1 = next(item for item in composite_manifests if item["case_slug"] == "figure1_relaxation_times_comparison")
    assert all(spec["axis_scale"] == "log" for spec in figure1["rendering"]["panel_specs"] if spec["plot_panel"] == "muB900.0")
    assert figure1["rendering"]["legend_host_panel"] == {"parameters": "row0_col0", "endpoints": "row0_col2"}
    assert "endpoint marker key" in figure1["rendering"]["legend_contents"]
    figure2 = next(item for item in composite_manifests if item["case_slug"] == "figure2_transport_coefficients_comparison")
    assert figure2["rendering"]["legend_host_panel"] == {"parameters": "row0_col0", "endpoints": "row0_col2"}
    assert figure2["rendering"]["legend_font_size_pt"] == 10.5


def test_v10_preserves_v5_tables_and_current_pointer():
    module = _module()
    package_path = module.ANALYSIS_ROOT / "manifest.json"
    package = json.loads(package_path.read_text(encoding="utf-8"))
    assert package["schema"] == "publication_clean_v10_png_review_package_v1"
    assert package["author_acceptance"] is None
    assert package["manuscript_eligible"] is False
    _stage_module().audit_v10_snapshot()
    for record in [*package["inputs"], *package["outputs"]]:
        if record["path"] not in _stage_module().HISTORICAL_CONTRACT_PATHS:
            assert validate_hash_record(record, root=ROOT, label="inputs", code_ref=_stage_module().FROZEN_CODE_REF) == []
    for name, expected in package["inherited_table_hashes"].items():
        assert sha256_file(module.V5_TABLE_ROOT / name) == expected
        assert sha256_file(module.TABLE_ROOT / name) == expected
    pointer = json.loads(module.CURRENT_POINTER.read_text(encoding="utf-8"))
    assert pointer["current_analysis_package"].endswith("publication_clean_v5")
    assert pointer["manuscript_eligible"] is True


def test_v10_sigma_ticks_are_exact_at_three_decimal_places():
    module = _module()
    ticks = module.sigma_ticks(0.0135, 0.0271)
    assert 4 <= len(ticks) <= 6
    assert all(abs(value - float(f"{value:.3f}")) < 1e-12 for value in ticks)
    assert len(set(f"{value:.3f}" for value in ticks)) == len(ticks)
    assert len(module.sigma_ticks(0.0249, 0.099)) >= 4


def test_v10_plot_coordinates_preserve_every_frozen_point_and_phase_gap():
    import matplotlib.pyplot as plt

    module = _module()
    points, _, _, _, gaps, _ = module.PARENT.load_v5_inputs()
    grouped, gap_map = module.group_inputs(points, gaps)
    profile = load_profile("candidate_aps_v2")
    configure_matplotlib(profile)
    figure, axis = plt.subplots()
    try:
        for mode, panel, observable in sorted({key[:3] for key in grouped}):
            axis.clear()
            module.render_panel(axis, mode, panel, observable, grouped, gap_map, profile)
            expected = []
            for name in sorted(key[3] for key in grouped if key[:3] == (mode, panel, observable)):
                expected.extend(module.V2.split_curve_segments(
                    grouped[(mode, panel, observable, name)], gap_map.get((mode, panel, name), [])))
            assert len(axis.lines) == len(expected)
            for line, segment in zip(axis.lines, expected):
                assert list(line.get_xdata()) == [float(row["xi"]) for row in segment]
                assert list(line.get_ydata()) == [float(row["clean_value"]) for row in segment]
    finally:
        plt.close(figure)
