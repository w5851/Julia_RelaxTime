from __future__ import annotations

import importlib.util
import math
from pathlib import Path


ROOT = Path(__file__).parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "relaxtime" / "audit_issue130_rs_pilot_quality.py"
SPEC = importlib.util.spec_from_file_location("issue130_rs_pilot_quality", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _row(*, temperature: str, xi: str = "-0.5", quality: bool = True) -> dict[str, str]:
    row = {
        "muB_MeV": "900.0",
        "xi": xi,
        "mode": "mode_a_fixed_muB_phase_scaled",
        "alpha_T": "1.0",
        "T_MeV": temperature,
        "quality_flag": "true" if quality else "false",
        "quality_reason": "tau_u_ubar_ratio_high" if quality else "ok",
        "quality_metric": "7.0" if quality else "1.5",
        "n_u": "10.0",
        "n_ubar": "0.1",
        "phase_reference_kind": "first_order",
        "phase_structure": "first_order_possible",
    }
    for field in MODULE.TRANSPORT_FIELDS:
        row[field] = "1.0"
    return row


def test_quality_audit_pairs_common_warning_and_reports_temperature_shift() -> None:
    rows = {
        "candidate_runtime": [_row(temperature="125.8")],
        "legacy": [_row(temperature="125.7")],
    }
    result = MODULE.build_quality_rows(rows)
    assert len(result) == 1
    assert math.isclose(result[0]["delta_T_MeV"], 0.1)
    assert result[0]["phase_kind_candidate"] == result[0]["phase_kind_legacy"]
    assert result[0]["phase_structure_candidate"] == result[0]["phase_structure_legacy"]


def test_quality_audit_ignores_nonwarning_rows() -> None:
    rows = {
        "candidate_runtime": [_row(temperature="125.8", quality=False)],
        "legacy": [_row(temperature="125.7", quality=False)],
    }
    assert MODULE.build_quality_rows(rows) == []
