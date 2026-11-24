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


def test_display_secants_use_native_endpoints_and_resolve_the_pseudoscalar_crossing() -> None:
    def grid(tier):
        rows = []
        for step in MODULE.STEPS[tier]:
            x = 0.000125 * step
            rows.append({"delta_xi": x, "eta_over_s": math.exp(2.0*x+3.0*x*x)})
        return rows

    rows, checks = MODULE.local_slopes({(meson, tier): grid(tier)
                                       for meson in MODULE.MESONS for tier in MODULE.TIERS})
    sigma = [row for row in rows if row["meson"] == "sigma"]
    assert len(sigma) == 16
    for row in sigma:
        assert row["u_right"] - row["u_left"] == pytest.approx(0.001)
        assert row["u_left"] * row["u_right"] >= 0
        assert row["slope"] == pytest.approx(2.0 + 6.0*row["u_midpoint"], abs=2e-12)
    kaon = [row for row in rows if row["meson"] == "K"]
    assert all(row["u_right"] - row["u_left"] >= 0.0005
               for row in kaon if -.002 <= row["u_left"] < -.0005)
    assert any(row["u_left"] == 0 and row["u_right"] == 0.000125 for row in kaon)
    assert sum(row["meson"] == "pi" for row in rows) == 32
    assert checks["sigma"]["max_shared_endpoint_slope_defect"] == 0.0


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
