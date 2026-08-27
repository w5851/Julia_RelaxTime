from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "relaxtime" / "build_phase_guided_publication_clean_v1.py"
SPEC = importlib.util.spec_from_file_location("phase_guided_publication_clean_v1", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _row(*, panel: str = "muB450.0", series: str = "alpha1.0", xi: str = "0.0", value: str = "1.0") -> dict[str, str]:
    return {
        "T_MeV": "200.0",
        "muB_MeV": "450.0",
        "xi": xi,
        "mode": "mode_a_fixed_muB_phase_scaled",
        "phase_reference_kind": "crossover",
        "phase_structure": "no_transition",
        "plot_panel": panel,
        "plot_series": series,
        "plot_series_label": series,
        "converged": "true",
        "quality_flag": "false",
        "quality_reason": "ok",
        "eta_over_s": value,
        "zeta_over_s": str(float(value) + 1.0),
        "sigma_over_T": str(float(value) + 2.0),
        "run_id": "fixture",
    }


def _loaded(rows: list[dict[str, str]], mode_key: str = "mode_a") -> dict[str, dict[str, object]]:
    index = {(row["plot_panel"], row["plot_series"], MODULE.canonical_xi(row["xi"])): row for row in rows}
    curves: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in rows:
        curves.setdefault((row["plot_panel"], row["plot_series"]), []).append(row)
    for curve in curves.values():
        curve.sort(key=lambda row: float(row["xi"]))
    return {mode_key: {"rows": rows, "index": index, "curves": curves, "manifest": {}}}


def test_replacement_recomputes_from_current_neighbours() -> None:
    rows = [_row(xi="0.0", value="10.0"), _row(xi="0.1", value="20.0"), _row(xi="0.2", value="40.0")]
    loaded = _loaded(rows)
    recipe = [{
        "window_id": "w",
        "mode_key": "mode_a",
        "plot_panel": "muB450.0",
        "plot_series": "alpha1.0",
        "observable": "eta_over_s",
        "xi": "0.1",
        "left_xi": "0.0",
        "right_xi": "0.2",
        "raw_production_value": "999.0",
        "paper_display_value": "998.0",
    }]
    derived = MODULE.build_replacement_map(loaded, recipe)
    assert len(derived) == 1
    assert derived[0]["raw_production_value_current"] == 20.0
    assert derived[0]["derived_display_value"] == 25.0
    assert derived[0]["recipe_raw_production_value"] == "999.0"


def test_replacement_requires_all_three_current_points() -> None:
    loaded = _loaded([_row(xi="0.0"), _row(xi="0.1")])
    recipe = [{
        "window_id": "w",
        "mode_key": "mode_a",
        "plot_panel": "muB450.0",
        "plot_series": "alpha1.0",
        "observable": "eta_over_s",
        "xi": "0.1",
        "left_xi": "0.0",
        "right_xi": "0.2",
        "raw_production_value": "1.0",
        "paper_display_value": "1.0",
    }]
    with pytest.raises(ValueError, match="recipe key missing"):
        MODULE.build_replacement_map(loaded, recipe)


def test_direct_coexistence_marker_reconciles_to_two_raw_side_points() -> None:
    rows = [_row(xi="-0.003", value="2.0"), _row(xi="0.003", value="4.0")]
    loaded = _loaded(rows)
    recipe = [{
        "window_id": "first_order",
        "mode_key": "mode_a",
        "plot_panel": "muB450.0",
        "plot_series": "alpha1.0",
        "observable": "eta_over_s",
        "xi": "0.0",
        "raw_production_value": "3.0",
        "marker": "star",
    }]
    markers = MODULE.build_marker_map(loaded, recipe)
    assert [row["render_xi"] for row in markers] == ["-0.0030000000", "0.0030000000"]
    assert all(row["marker_status"] == "reconciled_direct_coexistence_side_point" for row in markers)


def test_missing_non_mode_a_marker_is_an_error() -> None:
    loaded = _loaded([_row(xi="-0.003"), _row(xi="0.003")], mode_key="mode_b")
    recipe = [{
        "window_id": "first_order",
        "mode_key": "mode_b",
        "plot_panel": "muB450.0",
        "plot_series": "alpha1.0",
        "observable": "eta_over_s",
        "xi": "0.0",
        "raw_production_value": "3.0",
        "marker": "star",
    }]
    with pytest.raises(ValueError, match="marker key missing"):
        MODULE.build_marker_map(loaded, recipe)


def test_clean_points_separate_raw_interpolated_and_marker_status() -> None:
    rows = [_row(xi="-0.1", value="1.0"), _row(xi="0.0", value="2.0"), _row(xi="0.1", value="3.0")]
    loaded = _loaded(rows)
    replacement = [{
        "window_id": "w",
        "mode_key": "mode_a",
        "plot_panel": "muB450.0",
        "plot_series": "alpha1.0",
        "observable": "eta_over_s",
        "xi": "0.0000000000",
        "derived_display_value": 7.0,
    }]
    marker = [{
        "window_id": "m",
        "mode_key": "mode_a",
        "plot_panel": "muB450.0",
        "plot_series": "alpha1.0",
        "render_xi": "0.1000000000",
        "observable": "eta_over_s",
        "raw_production_value_current": 3.0,
    }]
    points = MODULE.build_clean_points(loaded, replacement, marker)
    by_xi = {row["xi"]: row for row in points if row["observable"] == "eta_over_s"}
    assert by_xi["0.0000000000"]["display_status"] == "interpolated_noncertified"
    assert by_xi["0.0000000000"]["clean_value"] == 7.0
    assert by_xi["0.1000000000"]["display_status"] == "first_order_raw_marker"
    assert by_xi["-0.1000000000"]["display_status"] == "raw"
    assert all(row["canonical_data_modified"] is False for row in points)


def test_curve_index_counts_derived_rows() -> None:
    points = [
        {"mode_key": "mode_a", "plot_panel": "p", "plot_series": "s", "observable": "eta_over_s", "display_status": "raw", "xi": "0.0"},
        {"mode_key": "mode_a", "plot_panel": "p", "plot_series": "s", "observable": "eta_over_s", "display_status": "interpolated_noncertified", "xi": "0.1"},
        {"mode_key": "mode_a", "plot_panel": "p", "plot_series": "s", "observable": "eta_over_s", "display_status": "first_order_raw_marker", "xi": "0.2"},
    ]
    rows = MODULE.build_curve_index(points)
    assert rows == [{
        "mode_key": "mode_a", "plot_panel": "p", "plot_series": "s", "observable": "eta_over_s",
        "point_count": 3, "replacement_count": 1, "marker_count": 1, "xi_min": 0.0, "xi_max": 0.2,
        "canonical_data_modified": False,
    }]
