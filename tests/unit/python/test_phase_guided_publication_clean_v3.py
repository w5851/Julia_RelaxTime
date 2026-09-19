import importlib.util
from pathlib import Path


ROOT = Path(__file__).parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "relaxtime" / "build_phase_guided_publication_clean_v3.py"
SPEC = importlib.util.spec_from_file_location("phase_guided_publication_clean_v3", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _point(xi: str, value: str, *, phase_structure: str = "crossover_possible") -> dict[str, str]:
    return {
        "mode_key": "mode_b",
        "mode": "mode_b_fixed_T_sparse_muB",
        "plot_panel": "T200.0",
        "plot_series": "muB0.0",
        "plot_series_label": "muB=0.0 MeV",
        "T_MeV": "200.0",
        "muB_MeV": "0.0",
        "xi": xi,
        "observable": "tau_u",
        "raw_value": value,
        "v2_clean_value": value,
        "clean_value": value,
        "display_status": "raw",
        "value_source": "fixture",
        "phase_structure": phase_structure,
        "phase_reference_kind": "crossover",
        "quality_flag": "false",
        "quality_reason": "ok",
        "run_id": "fixture",
        "canonical_data_modified": "False",
        "smoothing_window": "",
        "smoothing_anchor_points": "",
        "smoothing_fit_method": "",
        "protected_by_phase_gate": "",
    }


def _adjustment() -> dict[str, object]:
    return {
        "adjustment_id": "fixture_v3",
        "mode_key": "mode_b",
        "plot_panel": "T200.0",
        "plot_series": "muB0.0",
        "observable": "tau_u",
        "target_xi": 0.36,
        "anchor_left_xi": 0.35,
        "anchor_right_xi": 0.37,
        "phase_policy": "crossover_only",
        "method": "linear_interpolation",
        "reason": "fixture",
    }


def test_v3_interpolates_only_target_and_preserves_raw_and_v2() -> None:
    points = [_point(xi, value) for xi, value in (("0.35", "2.0"), ("0.36", "2.1"), ("0.37", "2.4"))]
    final, audit = MODULE.apply_adjustments(points, [_adjustment()], [])

    target = final[1]
    assert float(target["clean_value"]) == 2.2
    assert target["raw_value"] == "2.1"
    assert target["v2_clean_value"] == "2.1"
    assert target["display_status"] == "v3_residual_smoothed"
    assert target["protected_by_phase_gate"] == "True"
    assert final[0]["clean_value"] == "2.0"
    assert final[2]["clean_value"] == "2.4"
    assert audit[0]["phase_gate"] == "passed"
    assert audit[0]["method"] == "linear_interpolation"


def test_v3_rejects_first_order_labelled_window() -> None:
    points = [_point(xi, value, phase_structure="first_order_possible")
              for xi, value in (("0.35", "2.0"), ("0.36", "2.1"), ("0.37", "2.4"))]
    try:
        MODULE.apply_adjustments(points, [_adjustment()], [])
    except ValueError as error:
        assert "crossover gate failed" in str(error)
    else:  # pragma: no cover
        raise AssertionError("first-order labelled window must be rejected")
