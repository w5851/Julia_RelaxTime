"""Contract tests for the frozen Ward/Møller Mott-response assembler."""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts/analysis/relaxtime/assemble_ward_moller_mott_response.py"
SPEC = importlib.util.spec_from_file_location("ward_moller_mott_assembler", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_scan_steps_match_driver_contract() -> None:
    assert len(MODULE.STEPS["fine"]) == 33
    assert len(MODULE.STEPS["reference"]) == 11
    assert MODULE.STEPS["reference"][0] == -64
    assert MODULE.STEPS["reference"][-1] == 64
    assert 0 in MODULE.STEPS["fine"]


def test_secant_is_logarithmic_and_requires_positive_rates() -> None:
    left = {"eta_over_s": 1.0, "delta_xi": -0.1}
    right = {"eta_over_s": 2.0, "delta_xi": 0.1}
    # log(2) over a width 0.2 is the stored logarithmic slope.
    assert MODULE.secant(left, right) == pytest.approx(5.0 * math.log(2.0))
    with pytest.raises(ValueError, match="positive eta/s"):
        MODULE.secant({"eta_over_s": 0.0, "delta_xi": 0.0}, right)


def test_metadata_rejects_a_nonproduction_flux() -> None:
    metadata = {
        "core_commit": MODULE.CORE,
        "julia_version": "1.12.5",
        "runtime_sampling": {"relative_velocity": "cm"},
        "profile": {
            "source_sha256": MODULE.PROFILE,
            "a_polarization": "ward_matched_RS",
            "a_coupling": "ward_matched_RS",
            "b0_direction": "sphere_average",
            "background": {"p_nodes": 96, "cos_nodes": 24, "thermal_p_max_inv_fm": 25.0},
        },
    }
    with pytest.raises(ValueError, match="pinned Ward/Møller"):
        MODULE.check_meta(metadata)
