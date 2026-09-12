from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "relaxtime" / "build_phase_guided_publication_clean_v2.py"
SPEC = importlib.util.spec_from_file_location("phase_guided_publication_clean_v2", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _row(*, panel: str, series: str, xi: str, phase_curr: str) -> dict[str, str]:
    return {
        "plot_panel": panel,
        "plot_series": series,
        "xi": xi,
        "T_MeV": "120.0" if panel == "T120.0" else "125.73725802686799",
        "muB_MeV": "900.0",
        "phase_curr": phase_curr,
        "phase_reference_kind": "first_order",
        "phase_structure": "first_order_possible",
        "quality_flag": "false",
        "quality_reason": "ok",
        "run_id": "fixture",
    }


def _loaded(rows: list[dict[str, str]], mode_key: str) -> dict[str, dict[str, object]]:
    index = {(row["plot_panel"], row["plot_series"], MODULE.canonical_xi(row["xi"])): row for row in rows}
    curves: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in rows:
        curves.setdefault((row["plot_panel"], row["plot_series"]), []).append(row)
    for curve in curves.values():
        curve.sort(key=lambda row: float(row["xi"]))
    return {mode_key: {"rows": rows, "index": index, "curves": curves, "manifest": {}}}


def test_boundary_gap_map_uses_phase_curr_endpoints_and_preserves_source_bracket() -> None:
    mode_a_rows = [
        _row(panel="muB900.0", series="alpha1.0", xi="-0.003", phase_curr="quark"),
        _row(panel="muB900.0", series="alpha1.0", xi="0.003", phase_curr="hadron"),
    ]
    mode_b_rows = [
        _row(panel="T120.0", series="muB900.0", xi="-0.14", phase_curr="quark"),
        _row(panel="T120.0", series="muB900.0", xi="-0.13", phase_curr="quark"),
        _row(panel="T120.0", series="muB900.0", xi="-0.12", phase_curr="hadron"),
    ]
    loaded = {
        **_loaded(mode_a_rows, "mode_a"),
        **_loaded(mode_b_rows, "mode_b"),
    }
    records, table_rows = MODULE.build_boundary_gap_map(loaded)
    assert [(r["mode_key"], r["gap_xi_low"], r["gap_xi_high"]) for r in records] == [
        ("mode_a", -0.003, 0.003),
        ("mode_b", -0.13, -0.12),
    ]
    mode_b_rows_out = [r for r in table_rows if r["mode_key"] == "mode_b"]
    assert {r["source_interval_xi_low"] for r in mode_b_rows_out} == {"-0.14"}
    assert {r["source_interval_xi_high"] for r in mode_b_rows_out} == {"-0.13"}
    assert {r["phase_label"] for r in mode_b_rows_out} == {
        "手征恢复相 (quark)",
        "手征破缺相 (hadron)",
    }


def test_split_curve_segments_never_bridge_gap() -> None:
    rows = [{"xi": str(x)} for x in (-0.2, -0.003, 0.003, 0.1)]
    gaps = [{"gap_xi_low": -0.003, "gap_xi_high": 0.003}]
    segments = MODULE.split_curve_segments(rows, gaps)
    assert [[float(row["xi"]) for row in segment] for segment in segments] == [
        [-0.2, -0.003],
        [0.003, 0.1],
    ]


def test_midpoint_markers_are_audit_only_in_v2() -> None:
    rows = [{"marker_status": "confirmed_interval_midpoint", "render_marker": True, "reason": "old"}]
    out = MODULE.suppress_midpoint_markers(rows)
    assert out[0]["render_marker"] is False
    assert out[0]["marker_status"] == "suppressed_phase_switch_midpoint"
    assert "no CEP" in out[0]["marker_semantics"]


def test_phase_switch_inventory_distinguishes_crossover_bookkeeping() -> None:
    rows = [
        _row(panel="muB450.0", series="alpha1.0", xi="0.05", phase_curr="quark"),
        _row(panel="muB450.0", series="alpha1.0", xi="0.06", phase_curr="hadron"),
    ]
    rows[0]["phase_reference_kind"] = rows[1]["phase_reference_kind"] = "crossover"
    rows[0]["phase_structure"] = rows[1]["phase_structure"] = "crossover_possible"
    loaded = _loaded(rows, "mode_a")
    inventory = MODULE.build_phase_switch_inventory(loaded, [])
    assert len(inventory) == 1
    assert inventory[0]["physical_first_order_candidate"] is False
    assert inventory[0]["render_gap"] is False
