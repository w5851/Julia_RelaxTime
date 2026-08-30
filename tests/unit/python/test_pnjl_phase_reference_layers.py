from scripts.analysis.pnjl.build_issue130_phase_reference_layers import (
    bracket_xi,
    interpolate_surface,
    uniform_xi,
)


def test_uniform_xi_contract_is_deterministic():
    values = uniform_xi()
    assert len(values) == 161
    assert values[0] == -0.5
    assert values[-1] == 0.5
    assert all(round(right - left, 10) == 0.00625 for left, right in zip(values, values[1:]))


def test_surface_interpolation_uses_common_axis_support_only():
    rows = [
        {"xi": -0.5, "T_MeV": 10.0, "mu_MeV": 100.0, "rho": 1.0},
        {"xi": -0.5, "T_MeV": 20.0, "mu_MeV": 110.0, "rho": 2.0},
        {"xi": 0.0, "T_MeV": 15.0, "mu_MeV": 105.0, "rho": 1.5},
        {"xi": 0.0, "T_MeV": 20.0, "mu_MeV": 110.0, "rho": 2.0},
    ]
    derived, coverage = interpolate_surface(
        rows,
        surface="maxwell",
        axis="T_MeV",
        fields=("mu_MeV", "rho"),
        source_layer="synthetic",
    )
    assert derived
    assert {row["layer"] for row in derived} == {"strict_reference_v1", "interpolated_noncertified"}
    assert any(row["layer"] == "interpolated_noncertified" for row in derived)
    assert all(10.0 <= row["T_MeV"] <= 20.0 for row in derived)
    assert any(item["coverage_status"] == "interpolated_common_support" for item in coverage)


def test_xi_bracket_does_not_extrapolate():
    source = {"-0.5": [1], "0.0": [2]}
    assert bracket_xi(source, -0.25) == (-0.5, 0.0)
    assert bracket_xi(source, 0.5) is None
