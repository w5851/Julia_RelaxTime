from __future__ import annotations

from copy import deepcopy

import matplotlib.pyplot as plt
import pytest

from scripts.analysis.relaxtime import build_phase_guided_full_rs_raw as raw
from scripts.plotting.plot_series import validate_series_styles
from scripts.plotting.validate_plot_artifact import _check_declared_legends


def fixture_rows():
    rows = []
    for xi, value in ((0.5, 1), (-0.5, 1), (0, 100)):
        rows.append({
            **{field: str(value) for field in raw.OBSERVABLES},
            "xi": str(xi), "T_MeV": "200", "muq_MeV": "0", "muB_MeV": "0",
            "plot_panel": "muB0.0", "plot_series": "alpha1.0", "mode": raw.MODE_DIRS["mode_a"],
            "converged": "true", "quality_flag": "true" if xi == 0 else "false",
            "quality_reason": "tau_u_ubar_ratio_high" if xi == 0 else "ok",
            "phase_curr": "hadron" if xi > 0 else "quark", "phase_structure": "crossover_possible",
            "alpha_T": "1.0", "T_phase_base_MeV": "200",
        })
    return rows


def test_sorting_retains_flagged_spike_and_does_not_mutate_source():
    rows = fixture_rows()
    before = deepcopy(rows)
    groups = raw.validate_rows(rows, "mode_a")
    curve = groups[("mode_a", "muB0.0", "alpha1.0")]
    assert rows == before
    assert [row["tau_u"] for row in curve] == ["1", "100", "1"]
    assert curve[1]["quality_flag"] == "true"
    assert [row["csv_row"] for row in curve] == [3, 4, 2]
    # A crossover/large spike never creates a gap without independent evidence.
    assert raw.split_segments(curve, []) == [curve]


@pytest.mark.parametrize("field,value,message", [
    ("tau_u", "NaN", "non-finite"), ("eta", "-1", "non-positive"),
    ("converged", "false", "unconverged"), ("T_MeV", "201", "temperature identity"),
    ("quality_reason", "unreviewed_error", "unreviewed quality flag"),
])
def test_bad_source_is_rejected_instead_of_hidden(field, value, message):
    rows = fixture_rows()
    rows[2][field] = value
    with pytest.raises(ValueError, match=message):
        raw.validate_rows(rows, "mode_a")


def test_duplicate_source_point_is_rejected():
    rows = fixture_rows()
    rows.append(dict(rows[0]))
    with pytest.raises(ValueError, match="duplicate source point"):
        raw.validate_rows(rows, "mode_a")


def test_physical_split_retains_both_original_sides_and_rejects_missing_pair():
    groups = raw.validate_rows(fixture_rows(), "mode_a")
    rows = groups[("mode_a", "muB0.0", "alpha1.0")]
    segments = raw.split_segments(rows, [(0, 0.5)])
    assert [point for segment in segments for point in segment] == rows
    assert [len(segment) for segment in segments] == [2, 1]
    assert segments[0][-1]["tau_u"] == "100"
    with pytest.raises(ValueError, match="adjacent original points"):
        raw.split_segments(rows, [(-0.25, 0.25)])


@pytest.fixture(scope="module")
def accepted_inputs():
    return raw.load_inputs(), raw.load_profile("candidate_aps_v2")


def test_new_grid_identity_and_reviewed_boundaries(accepted_inputs):
    data, _ = accepted_inputs
    assert sum(map(len, data.groups.values())) == 1819
    assert sum(row["quality_flag"] == "true" for rows in data.groups.values() for row in rows) == 470
    assert data.gaps == {
        ("mode_a", "muB900.0", "alpha1.0"): [(-0.003, 0.003)],
        ("mode_b", "T120.0", "muB900.0"): [(-0.13, -0.12)],
    }
    assert 0.0 not in {float(row["xi"]) for row in data.groups[("mode_a", "muB900.0", "alpha1.0")]}


@pytest.mark.parametrize("mode", ["mode_a", "mode_b"])
def test_parameter_encoding_is_independent_of_selection_order(mode):
    profile = raw.load_profile("candidate_aps_v2")
    registry = raw.make_registry(mode, profile)
    identities = [style.series_id for style in registry.styles]
    assert registry.select(reversed(identities)) == tuple(reversed(registry.styles))
    assert registry.select(identities[1:]) == registry.styles[1:]


@pytest.mark.parametrize("kind,args", [
    ("composite", ("figure1_relaxation_times_comparison",)),
    ("composite", ("figure2_transport_coefficients_comparison",)),
    ("detail", ()),
    ("single", ("mode_b", "T120.0", "tau_u")),
])
def test_raw_artists_and_final_size_layout(kind, args, accepted_inputs):
    data, profile = accepted_inputs
    figure, specs, placements = raw.render_request(kind, args, data, profile)
    try:
        assert raw.verify_raw_vertices(figure, specs, data) > 0
        figure.set_dpi(profile.dpi)
        quality = raw.measure_figure(figure, intended_width_inches=6.75)
        assert raw.clear_layout(quality), {key: quality[key] for key in raw.LAYOUT_ERRORS}
        assert quality["minimum_capital_numeral_height_mm"] >= 2
        assert quality["minimum_curve_linewidth_pt"] >= 0.5
        assert all(tick["inward"] and tick["both_sides"] and tick["minor_count"] for tick in quality["tick_axes"])
        errors = []
        _check_declared_legends(dict(legend_placements=placements,
                                    legend_outside=all(p["placement"] == "outside_axes" for p in placements)),
                               quality, errors)
        assert errors == []
        registry = raw.make_registry(specs[0]["mode_key"], profile)
        assert validate_series_styles(registry.record([s.series_id for s in registry.styles]),
                                      [s.record() for s in registry.styles], quality["series_styles"]) == []
        if kind == "detail":
            assert [axis.get_ylim() for axis in figure.axes] == [detail[-1] for detail in raw.DETAILS]
        elif kind == "composite":
            assert all(axis.get_yscale() == ("log" if args[0].startswith("figure1") else "linear") for axis in figure.axes)
        # Deliberately corrupt one artist: the audit must actually detect it.
        line = figure.axes[0].lines[0]
        corrupted = list(line.get_ydata())
        corrupted[0] *= 1.01
        line.set_ydata(corrupted)
        with pytest.raises(ValueError, match="raw vertices changed"):
            raw.verify_raw_vertices(figure, specs, data)
    finally:
        plt.close(figure)


def test_build_refuses_existing_case_before_reading_sources(monkeypatch, tmp_path):
    monkeypatch.setattr(raw, "FIGURE_PARENT", tmp_path)
    (tmp_path / "already_exists").mkdir()
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        raw.build_review("already_exists")
