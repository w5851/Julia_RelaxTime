import pytest

from scripts.analysis.pnjl.render_issue130_phase_surface_v9 import (
    extend_maxwell_groups_to_temperature_floor,
)


def _row(xi, temperature, mu):
    return {
        "surface": "maxwell",
        "xi": str(xi),
        "T_MeV": str(temperature),
        "mu_MeV": str(mu),
    }


def test_low_temperature_extension_is_display_only_and_linear_from_two_native_rows():
    groups, trace = extend_maxwell_groups_to_temperature_floor(
        {0.0: [_row(0.0, 5.0, 10.0), _row(0.0, 7.0, 9.0), _row(0.0, 9.0, 8.0)]},
        floor_T_MeV=0.0,
    )

    assert len(groups[0.0]) == 4
    synthetic = groups[0.0][0]
    assert synthetic["T_MeV"] == 0.0
    assert synthetic["mu_MeV"] == 12.5
    assert synthetic["layer"] == "display_only_extrapolated_noncertified"
    assert trace[0]["synthetic_row"] is True
    assert trace[0]["extension_gap_MeV"] == 5.0


def test_native_support_is_not_duplicated():
    groups, trace = extend_maxwell_groups_to_temperature_floor(
        {0.0: [_row(0.0, 0.0, 10.0), _row(0.0, 2.0, 9.0)]},
        floor_T_MeV=0.0,
    )

    assert len(groups[0.0]) == 2
    assert trace[0]["synthetic_row"] is False
    assert trace[0]["extension_method"] == "native_support"


def test_extension_requires_two_native_temperature_rows():
    with pytest.raises(ValueError, match="lacks two rows"):
        extend_maxwell_groups_to_temperature_floor(
            {0.0: [_row(0.0, 5.0, 10.0)]},
            floor_T_MeV=0.0,
        )
