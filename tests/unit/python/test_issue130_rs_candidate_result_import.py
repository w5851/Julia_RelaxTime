from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "relaxtime" / "import_issue130_rs_candidate_results.py"
SPEC = importlib.util.spec_from_file_location("issue130_rs_candidate_result_import", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _scan_row(*, mode: str, xi: str = "0.1", muB: str = "900.0", alpha: str = "1.0") -> dict[str, str]:
    row = {field: "1.0" for field in (
        "T_MeV", "muq_MeV", "muB_MeV", "xi", "alpha_T", "rho_baryon",
        "tau_u", "tau_d", "tau_s", "tau_ubar", "tau_dbar", "tau_sbar",
        "eta", "sigma", "zeta", "eta_over_s", "zeta_over_s", "sigma_over_T",
    )}
    row.update({
        "mode": mode,
        "phase_reference_kind": "crossover",
        "converged": "true",
        "phase_structure": "no_transition",
        "quality_flag": "false",
        "quality_reason": "ok",
        "muB_MeV": muB,
        "xi": xi,
        "alpha_T": alpha,
    })
    return row


def test_scan_validation_accepts_nan_alpha_only_for_fixed_temperature() -> None:
    row = _scan_row(mode=MODULE.MODES[1])
    row["alpha_T"] = "NaN"
    summary = MODULE.validate_scan_rows(
        MODULE.MODES[1], list(row), [row], expected_rows=1
    )
    assert summary["rows"] == 1


def test_scan_validation_rejects_duplicate_keys() -> None:
    first = _scan_row(mode=MODULE.MODES[0], xi="0.1")
    second = _scan_row(mode=MODULE.MODES[0], xi="0.1")
    with pytest.raises(ValueError, match="duplicate keys"):
        MODULE.validate_scan_rows(MODULE.MODES[0], list(first), [first, second])


def test_scan_validation_rejects_nonfinite_transport_value() -> None:
    row = _scan_row(mode=MODULE.MODES[0])
    row["tau_u"] = "NaN"
    with pytest.raises(ValueError, match="non-finite"):
        MODULE.validate_scan_rows(MODULE.MODES[0], list(row), [row])


def test_direct_coexistence_requires_plus_minus_replacement() -> None:
    rows = [
        _scan_row(mode=MODULE.MODES[0], xi="-0.003"),
        _scan_row(mode=MODULE.MODES[0], xi="0.003"),
    ]
    gate = MODULE.validate_direct_coexistence(rows)
    assert gate["has_minus_0003"] is True
    assert gate["has_plus_0003"] is True
    assert gate["has_zero"] is False

    rows.append(_scan_row(mode=MODULE.MODES[0], xi="0.0"))
    with pytest.raises(ValueError, match="replacement contract"):
        MODULE.validate_direct_coexistence(rows)


def test_diagnostic_validation_rejects_negative_contribution() -> None:
    row = {
        "T_MeV": "120.0",
        "muB_MeV": "0.0",
        "xi": "0.1",
        "species": "u",
        "channel": "elastic",
        "density_key": "u",
        "multiplicity": "1.0",
        "density": "1.0",
        "rate": "1.0",
        "contribution": "-1e-3",
        "total": "1.0",
        "tau_inv_species": "1.0",
    }
    with pytest.raises(ValueError, match="negative"):
        MODULE.validate_diag_rows(MODULE.MODES[0], list(row), [row])


def test_tree_hash_is_order_independent(tmp_path: Path) -> None:
    root = tmp_path / "case"
    (root / "b").mkdir(parents=True)
    (root / "a.txt").write_text("a", encoding="ascii")
    (root / "b" / "c.txt").write_text("c", encoding="ascii")
    first = MODULE.tree_hash(root)
    (root / "b" / "d.txt").write_text("d", encoding="ascii")
    assert MODULE.tree_hash(root) != first
