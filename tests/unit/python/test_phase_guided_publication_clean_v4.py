import importlib.util
import math
from pathlib import Path


ROOT = Path(__file__).parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "relaxtime" / "build_phase_guided_publication_clean_v4.py"
SPEC = importlib.util.spec_from_file_location("phase_guided_publication_clean_v4", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _point(
    *,
    mode_key: str,
    panel: str,
    series: str,
    observable: str,
    xi: str,
    value: str,
    phase_kind: str,
    phase_structure: str,
) -> dict[str, str]:
    return {
        "mode_key": mode_key,
        "mode": "mode_b_fixed_T_sparse_muB" if mode_key == "mode_b" else "mode_a_fixed_muB_phase_scaled",
        "plot_panel": panel,
        "plot_series": series,
        "plot_series_label": series,
        "T_MeV": "200.0",
        "muB_MeV": "0.0",
        "xi": xi,
        "observable": observable,
        "raw_value": value,
        "v2_clean_value": value,
        "clean_value": value,
        "display_status": "raw",
        "value_source": "fixture",
        "phase_structure": phase_structure,
        "phase_reference_kind": phase_kind,
        "quality_flag": "false",
        "quality_reason": "ok",
        "run_id": "fixture",
        "canonical_data_modified": "False",
        "smoothing_window": "",
        "smoothing_anchor_points": "",
        "smoothing_fit_method": "",
        "protected_by_phase_gate": "",
    }


def _interior_adjustment() -> dict[str, object]:
    return {
        "adjustment_id": "fixture_v4_interior",
        "mode_key": "mode_b",
        "plot_panel": "T200.0",
        "plot_series": "muB0.0",
        "observable": "zeta",
        "target_xi": 0.36,
        "anchor_left_xi": 0.35,
        "anchor_right_xi": 0.37,
        "phase_policy": "crossover_only",
        "method": "linear_interpolation",
        "fit_space": "linear",
        "mechanism_status": "inherited_parent_scope",
        "reason": "fixture",
    }


def _endpoint_adjustment() -> dict[str, object]:
    return {
        "adjustment_id": "fixture_v4_endpoint",
        "mode_key": "mode_a",
        "plot_panel": "muB900.0",
        "plot_series": "alpha1.0",
        "observable": "tau_sbar",
        "target_xi": -0.003,
        "anchor_left_xi": -0.02,
        "anchor_right_xi": -0.01,
        "phase_policy": "left_branch_endpoint_only",
        "method": "log_linear_extrapolation",
        "fit_space": "log",
        "mechanism_status": "small_denominator_supported_with_upstream_endpoint_sensitivity",
        "reason": "fixture",
    }


def test_v4_interpolates_composite_point_and_preserves_parent_inputs() -> None:
    points = [
        _point(
            mode_key="mode_b",
            panel="T200.0",
            series="muB0.0",
            observable="zeta",
            xi=xi,
            value=value,
            phase_kind="crossover",
            phase_structure="crossover_possible",
        )
        for xi, value in (("0.35", "0.68"), ("0.36", "0.70"), ("0.37", "0.74"))
    ]
    final, audit = MODULE.apply_adjustments(points, [_interior_adjustment()], [])
    assert math.isclose(float(final[1]["clean_value"]), 0.71)
    assert final[1]["raw_value"] == "0.70"
    assert final[1]["v2_clean_value"] == "0.70"
    assert final[1]["display_status"] == "v4_residual_smoothed"
    assert audit[0]["phase_gate"] == "passed_crossover_only"
    assert audit[0]["adjustment_type"] == "interior_crossover_display_interpolation"


def test_v4_allows_only_a_left_branch_endpoint_extrapolation() -> None:
    points = [
        _point(
            mode_key="mode_a",
            panel="muB900.0",
            series="alpha1.0",
            observable="tau_sbar",
            xi=xi,
            value=value,
            phase_kind="first_order",
            phase_structure="first_order_possible",
        )
        for xi, value in (("-0.02", "2.2993789650830054"), ("-0.01", "2.426591623490907"), ("-0.003", "2.547560459778322"))
    ]
    gap = [{
        "mode_key": "mode_a",
        "plot_panel": "muB900.0",
        "plot_series": "alpha1.0",
        "gap_xi_low": "-0.003",
        "gap_xi_high": "0.003",
    }]
    final, audit = MODULE.apply_adjustments(points, [_endpoint_adjustment()], gap)
    expected = math.exp(
        math.log(2.426591623490907)
        + 0.7 * (math.log(2.426591623490907) - math.log(2.2993789650830054))
    )
    assert math.isclose(float(final[-1]["clean_value"]), expected)
    assert final[-1]["raw_value"] == "2.547560459778322"
    assert audit[0]["phase_gate"] == "passed_left_branch_endpoint_only"
    assert audit[0]["endpoint_log_slope_before"] != audit[0]["endpoint_log_slope_after"]
    assert float(audit[0]["relative_change_from_parent"]) < 0.0


def test_v4_endpoint_gate_rejects_non_audited_or_cross_gap_target() -> None:
    points = [
        _point(
            mode_key="mode_a",
            panel="muB900.0",
            series="alpha1.0",
            observable="tau_sbar",
            xi=xi,
            value=value,
            phase_kind="first_order",
            phase_structure="first_order_possible",
        )
        for xi, value in (("-0.02", "2.3"), ("-0.01", "2.4"), ("-0.003", "2.5"))
    ]
    wrong_gap = [{
        "mode_key": "mode_a",
        "plot_panel": "muB900.0",
        "plot_series": "alpha1.0",
        "gap_xi_low": "-0.004",
        "gap_xi_high": "0.003",
    }]
    try:
        MODULE.apply_adjustments(points, [_endpoint_adjustment()], wrong_gap)
    except ValueError as error:
        assert "phase gap" in str(error) or "audited left gap endpoint" in str(error)
    else:  # pragma: no cover
        raise AssertionError("endpoint adjustment must require the audited gap")


def test_v4_crossover_gate_rejects_first_order_window() -> None:
    points = [
        _point(
            mode_key="mode_b",
            panel="T200.0",
            series="muB0.0",
            observable="zeta",
            xi=xi,
            value=value,
            phase_kind="first_order",
            phase_structure="first_order_possible",
        )
        for xi, value in (("0.35", "0.68"), ("0.36", "0.70"), ("0.37", "0.74"))
    ]
    try:
        MODULE.apply_adjustments(points, [_interior_adjustment()], [])
    except ValueError as error:
        assert "crossover gate failed" in str(error)
    else:  # pragma: no cover
        raise AssertionError("crossover-only adjustment must reject first-order rows")
