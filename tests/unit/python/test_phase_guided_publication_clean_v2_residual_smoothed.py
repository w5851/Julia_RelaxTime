import importlib.util
from pathlib import Path


ROOT = Path(__file__).parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "relaxtime" / "build_phase_guided_publication_clean_v2_residual_smoothed.py"
SPEC = importlib.util.spec_from_file_location("phase_guided_publication_clean_v2_residual_smoothed", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _point(*, xi: str, value: str, phase_kind: str = "crossover", phase_structure: str = "no_transition") -> dict[str, str]:
    return {
        "mode_key": "mode_b",
        "mode": "mode_b_fixed_T_sparse_muB",
        "plot_panel": "T200.0",
        "plot_series": "muB0.0",
        "plot_series_label": "muB=0.0 MeV",
        "T_MeV": "200.0",
        "muB_MeV": "0.0",
        "xi": xi,
        "observable": "tau_dbar",
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
        "canonical_data_modified": False,
    }


def _window(*, policy: str = "crossover_only") -> dict[str, object]:
    return {
        "window_id": "fixture_window",
        "scope": "test",
        "mode_key": "mode_b",
        "plot_panel": "T200.0",
        "plot_series": "muB0.0",
        "anchor_left_xi": -0.2,
        "anchor_right_xi": 0.2,
        "strength": 1.0,
        "value_space": "log",
        "phase_policy": policy,
        "observables_tuple": ("tau_dbar",),
    }


def _phase_index(rows: list[dict[str, str]], phase_curr: str = "quark") -> dict[tuple[str, str, str, str], dict[str, str]]:
    return {
        ("mode_b", "T200.0", "muB0.0", MODULE.V2.canonical_xi(row["xi"])): {
            "phase_curr": phase_curr,
        }
        for row in rows
    }


def test_monotone_log_trend_preserves_anchors_and_is_positive() -> None:
    rows = [_point(xi=str(x), value=str(y)) for x, y in ((-0.2, 1.0), (-0.1, 1.2), (0.0, 2.0), (0.1, 2.1), (0.2, 2.2))]
    trend = MODULE._trend_values(rows, 0, 4, "clean_value", [], "log")
    assert trend["-0.2000000000"] == 1.0
    assert trend["0.2000000000"] == 2.2
    assert all(value > 0.0 for value in trend.values())
    ordered = [trend[MODULE.V2.canonical_xi(row["xi"])] for row in rows]
    assert ordered == sorted(ordered)


def test_residual_smoothing_changes_interior_only() -> None:
    rows = [_point(xi=str(x), value=str(y)) for x, y in ((-0.2, 1.0), (-0.1, 1.01), (0.0, 0.4), (0.1, 1.3), (0.2, 1.4))]
    out, point_audit, window_audit, curve_audit = MODULE.apply_residual_smoothing(
        rows,
        [_window()],
        _phase_index(rows),
        [],
    )
    assert len(point_audit) == 3
    assert window_audit[0]["action"] == "applied"
    assert float(out[0]["clean_value"]) == 1.0
    assert float(out[-1]["clean_value"]) == 1.4
    assert all(row["display_status"] == "residual_smoothed" for row in out[1:-1])
    assert all(row["smoothing_window"] == "fixture_window" for row in out[1:-1])
    assert all(row["protected_by_phase_gate"] is True for row in out[1:-1])
    assert all("local_residual_before" in item and "local_residual_after" in item for item in point_audit)
    assert curve_audit[0]["action"] == "applied"


def test_roughness_guard_retains_a_curve_when_candidate_is_worse(monkeypatch) -> None:
    xi_values = ("-0.2", "-0.1", "0.0", "0.1", "0.2")
    values = ("1.0", "1.1", "1.2", "1.3", "1.4")
    rows = [_point(xi=xi, value=value) for xi, value in zip(xi_values, values)]
    window = _window()
    def deliberately_rough_trend(rows, left_index, right_index, value_field, gaps, value_space):
        return {
            MODULE.V2.canonical_xi(row["xi"]): value
            for row, value in zip(rows[left_index : right_index + 1], (1.0, 1.8, 0.5, 1.8, 1.4))
        }

    monkeypatch.setattr(MODULE, "_trend_values", deliberately_rough_trend)
    out, point_audit, window_audit, curve_audit = MODULE.apply_residual_smoothing(
        rows,
        [window],
        _phase_index(rows),
        [],
    )
    assert not point_audit
    assert window_audit[0]["retained_curve_count"] == 1
    assert curve_audit[0]["action"] == "retained_v2"
    assert curve_audit[0]["candidate_roughness"] > curve_audit[0]["source_roughness"]
    assert [row["clean_value"] for row in out] == list(values)


def test_branch_window_may_use_one_labelled_branch_but_gap_is_hard_barrier() -> None:
    rows = [_point(xi=str(x), value=str(y), phase_kind="first_order", phase_structure="first_order_possible")
            for x, y in ((-0.2, 1.0), (-0.1, 1.1), (0.0, 1.2), (0.1, 1.3), (0.2, 1.4))]
    phase = _phase_index(rows, phase_curr="quark")
    branch_window = _window(policy="branch_local")
    out, audit, _, _ = MODULE.apply_residual_smoothing(rows, [branch_window], phase, [])
    assert len(audit) == 3
    assert out[2]["display_status"] == "residual_smoothed"

    gap = [{"gap_xi_low": -0.05, "gap_xi_high": 0.05}]
    out, audit, window_audit, _ = MODULE.apply_residual_smoothing(
        [_point(xi=str(x), value=str(y)) for x, y in ((-0.2, 1.0), (-0.1, 1.1), (0.0, 1.2), (0.1, 1.3), (0.2, 1.4))],
        [branch_window],
        phase,
        gap,
    )
    assert not audit
    assert window_audit[0]["action"] == "retained_v2"
    assert "gap" in window_audit[0]["gate_errors"]


def test_curve_scoped_gap_does_not_block_another_curve() -> None:
    rows = [_point(xi=str(x), value=str(y)) for x, y in ((-0.2, 1.0), (-0.1, 1.01), (0.0, 0.4), (0.1, 1.3), (0.2, 1.4))]
    phase = _phase_index(rows)
    unrelated_gap = [{
        "mode_key": "mode_a",
        "plot_panel": "muB900.0",
        "plot_series": "alpha1.0",
        "gap_xi_low": -0.05,
        "gap_xi_high": 0.05,
    }]
    out, audit, window_audit, _ = MODULE.apply_residual_smoothing(
        rows,
        [_window()],
        phase,
        unrelated_gap,
    )
    assert len(audit) == 3
    assert window_audit[0]["action"] == "applied"
    assert out[2]["display_status"] == "residual_smoothed"
