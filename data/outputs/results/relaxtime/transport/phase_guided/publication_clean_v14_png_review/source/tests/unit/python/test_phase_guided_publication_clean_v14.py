from __future__ import annotations

import json

import matplotlib.pyplot as plt
import pytest

from scripts.analysis.relaxtime import build_phase_guided_publication_clean_v14 as delivery
from scripts.analysis.relaxtime import build_phase_guided_full_rs_raw as raw
from scripts.plotting.validate_plot_artifact import _check_declared_legends


def test_chart_inventory_preserves_v13_and_adds_fixed_T():
    requests = list(delivery.chart_requests())
    assert requests[:-1] == list(raw.chart_requests())
    assert requests[-1] == (delivery.FIG4, "composite", ("figure4_transport_fixed_T",))
    assert len(requests) == 76


def test_publication_palette_is_explicit_and_does_not_change_shared_profile():
    original = raw.load_profile("strict_aps_v2")
    profile = delivery.publication_profile()
    assert profile.colors[:3] == delivery.LEGACY_COLORS
    assert profile.data["palette"]["parameter_linestyles"][:3] == ["-", "--", "-."]
    assert original.colors[2] == "#008565"
    assert profile.data["line_width_pt"] == original.data["line_width_pt"] == 1.2


@pytest.mark.parametrize("kind,arguments", [
    ("composite", ("figure4_transport_fixed_T",)),
    *[("single", args) for args in delivery.SPECIAL_SINGLES],
])
def test_aligned_artists_preserve_original_values_and_have_clear_legends(kind, arguments):
    data, profile = raw.load_inputs(), delivery.publication_profile()
    figure, specs, placements = delivery.render_request(kind, arguments, data, profile)
    try:
        vertices = raw.verify_raw_vertices(figure, specs, data)
        quality = raw.measure_figure(figure, intended_width_inches=6.75)
        assert raw.clear_layout(quality), {key: quality[key] for key in raw.LAYOUT_ERRORS}
        assert quality["minimum_capital_numeral_height_mm"] >= 2
        issues = []
        _check_declared_legends(dict(legend_placements=placements, legend_outside=False), quality, issues)
        assert not issues
        if kind == "composite":
            assert vertices == 909 * 3
            assert [axis.get_yscale() for axis in figure.axes] == ["log"] * 6 + ["linear"] * 3
            assert [spec["plot_panel"] for spec in specs] == ["T120.0", "T160.0", "T200.0"] * 3
            assert [spec["observable"] for spec in specs] == ["eta_over_s"] * 3 + ["zeta_over_s"] * 3 + ["sigma_over_T"] * 3
            assert [(p["role"], p["host_axes_index"]) for p in placements] == [("endpoints", 0), ("parameters", 6)]
            endpoints = [e for spec in specs for e in spec["endpoints"]]
            assert len(endpoints) == 6
            assert {e["xi"] for e in endpoints} == {-0.13, -0.12}
        elif arguments == ("mode_b", "T160.0", "eta_over_s"):
            assert figure.axes[0].get_yscale() == "linear"
        else:
            assert placements[0]["location"] == delivery.SPECIAL_SINGLES[arguments]
        original = delivery.render_signature(figure, specs, placements)
        figure.axes[0].set_ylabel("incorrect units")
        assert delivery.render_signature(figure, specs, placements) != original
    finally:
        plt.close(figure)


def test_vector_export_requires_actual_bound_authorization(tmp_path, monkeypatch):
    monkeypatch.setattr(raw, "FIGURE_PARENT", tmp_path / "figures")
    monkeypatch.setattr(raw, "RESULT_PARENT", tmp_path / "results")
    authorization = tmp_path / "authorization.json"
    authorization.write_text(json.dumps({"condition_verified": True}), encoding="utf-8")
    with pytest.raises(ValueError, match="actual author instruction"):
        delivery.export_delivery("review", "delivery", authorization)
    assert not (tmp_path / "figures/delivery").exists()
