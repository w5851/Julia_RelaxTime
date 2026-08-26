from __future__ import annotations

import csv
import importlib.util
from pathlib import Path


ROOT = Path(__file__).parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "relaxtime" / "render_issue130_rs_candidate_legacy_figures.py"
SPEC = importlib.util.spec_from_file_location("issue130_rs_candidate_legacy_figures", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _write_scan(path: Path, *, mode: str, duplicate: bool = False) -> None:
    fields = list(MODULE.COMMON_COLUMNS) + ["alpha_T"]
    row = {field: "1.0" for field in fields}
    row.update(
        {
            "T_MeV": "120.0",
            "muB_MeV": "450.0",
            "xi": "-0.5",
            "mode": mode,
            "plot_panel": "T120.0" if mode.startswith("mode_b") else "muB450.0",
            "plot_series_label": "muB=450.0 MeV" if mode.startswith("mode_b") else "alpha_T=1.0",
        }
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("# diagnostic input\n")
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerow(row)
        if duplicate:
            writer.writerow(row)


def test_validate_scan_accepts_comment_header_and_checks_plot_contract(tmp_path: Path) -> None:
    path = tmp_path / "mode_a.csv"
    _write_scan(path, mode="mode_a_fixed_muB_phase_scaled")
    result = MODULE.validate_scan(path, "mode_a_fixed_muB_phase_scaled")
    assert result["rows"] == 1
    assert result["unique_keys"] == 1
    assert result["solver_called"] is False


def test_validate_scan_rejects_duplicate_transport_keys(tmp_path: Path) -> None:
    path = tmp_path / "mode_b.csv"
    _write_scan(path, mode="mode_b_fixed_T_sparse_muB", duplicate=True)
    try:
        MODULE.validate_scan(path, "mode_b_fixed_T_sparse_muB")
    except ValueError as exc:
        assert "duplicate numerical keys" in str(exc)
    else:
        raise AssertionError("duplicate keys must be rejected")


def test_validate_scan_rejects_wrong_mode(tmp_path: Path) -> None:
    path = tmp_path / "mode_a.csv"
    _write_scan(path, mode="mode_a_fixed_muB_phase_scaled")
    try:
        MODULE.validate_scan(path, "mode_b_fixed_T_sparse_muB")
    except ValueError as exc:
        assert "mode values" in str(exc)
    else:
        raise AssertionError("wrong mode must be rejected")
