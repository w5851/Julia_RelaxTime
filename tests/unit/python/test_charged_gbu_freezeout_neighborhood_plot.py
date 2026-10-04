from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts/analysis/relaxtime/plot_charged_gbu_freezeout_neighborhood.py"
SPEC = importlib.util.spec_from_file_location("freezeout_neighborhood_plot", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)
COEFFICIENTS = {"a_GeV": .100, "b_GeV_inv1": 0, "c_GeV_inv3": 0,
                "d_GeV": 1.308, "e_GeV_inv1": .273}


def test_probe_changes_only_a_not_energy_mapping():
    mu, T = MODULE.energy_point(3, COEFFICIENTS)
    assert mu == pytest.approx(1308 / (1 + .273 * 3))
    assert MODULE.freezeout_T(mu, COEFFICIENTS, -10) == pytest.approx(T - 10)
    assert MODULE.freezeout_T(mu, COEFFICIENTS, 10) == pytest.approx(T + 10)


def test_masked_corner_blocks_contours_and_negative_ratio_is_preserved():
    import numpy as np

    rows = [{"T": T, "muB": mu, "status": "screened", "R": -.2 if mu == 25 else .5}
            for T in (90, 100, 110) for mu in (0, 25, 50)]
    rows[4]["R"] = None
    data = MODULE.neighborhood(rows, COEFFICIENTS)
    assert data["failed"].sum() == 1
    assert np.isnan(data["R"][1, 1])
    assert data["R"][0, 1] == -.2
    assert not data["valid_cells"].any()


def test_outside_band_is_not_zero_filled():
    import numpy as np

    rows = [{"T": T, "muB": mu, "status": "screened", "R": .5}
            for T in (70, 100, 130) for mu in (0, 25)]
    data = MODULE.neighborhood(rows, COEFFICIENTS)
    assert data["selected"].sum() == 2
    assert np.isnan(data["R"][[0, 2]]).all()
    assert not data["failed"].any()


def test_entrypoint_is_solver_free_and_no_overwrite():
    source = SCRIPT.read_text(encoding="utf-8")
    assert "Solver-free" in source
    assert "refusing to overwrite existing case" in source
    assert "corner_mask=False" in source
    assert "numerical_rescan\": False" in source
    assert "No fit or rescan" in source


def test_loader_rejects_changed_csv_hash(tmp_path):
    csv = tmp_path / "points.csv"
    csv.write_text("changed", encoding="utf-8")
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"merged_csv": {"sha256": "0" * 64}}), encoding="utf-8")
    with pytest.raises(ValueError, match="CSV hash"):
        MODULE.load_frozen(csv, manifest, tmp_path / "profile.toml")


def test_loader_rejects_changed_freezeout_profile(tmp_path):
    csv = tmp_path / "points.csv"
    csv.write_text("fixture", encoding="utf-8")
    profile = tmp_path / "profile.toml"
    profile.write_text("changed", encoding="utf-8")
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"merged_csv": {"sha256": MODULE.digest(csv)},
                                   "reference_lines": {"freezeout": {"source_sha256": "0" * 64}}}), encoding="utf-8")
    with pytest.raises(ValueError, match="freeze-out profile"):
        MODULE.load_frozen(csv, manifest, profile)


def test_fifty_MeV_neighborhood_uses_all_existing_temperature_support():
    rows = [{"T": T, "muB": mu, "status": "screened", "R": .5}
            for T in (40, 50, 100, 150, 200, 220) for mu in (0, 25)]
    data = MODULE.neighborhood(rows, COEFFICIENTS, band_MeV=50, T_bounds_MeV=(40, 220))
    assert list(data["T"]) == [40, 50, 100, 150, 200, 220]
    assert data["selected"].sum() == 6
    assert not data["failed"].any()


def test_contour_levels_are_uniform_point_zero_five():
    levels = MODULE.uniform_levels(.279, .744)
    assert levels == [.30, .35, .40, .45, .50, .55, .60, .65, .70]
    assert all(b - a == pytest.approx(.05) for a, b in zip(levels, levels[1:]))


def test_presentation_has_no_probe_curves_or_energy_labels():
    source = SCRIPT.read_text(encoding="utf-8")
    presentation = source.split("def plot_presentation", 1)[1].split("def run_presentation", 1)[0]
    assert "energy_point(" not in presentation
    assert '"candidate_curves": False' in presentation
    assert '"energy_labels": False' in presentation
    assert "inline=True" in presentation
    assert "expected_labels" in presentation


def test_labels_are_every_point_one_while_contours_stay_point_zero_five():
    levels = MODULE.uniform_levels(.0475, 1.129)
    labeled = [value for value in levels if MODULE.is_labeled_level(value)]
    assert labeled == [.1, .2, .3, .4, .5, .6, .7, .8, .9, 1., 1.1]
    assert .35 in levels and .35 not in labeled
    assert MODULE.is_labeled_level(.3 + 1e-15)


def test_presentation_has_no_white_label_boxes():
    source = SCRIPT.read_text(encoding="utf-8")
    presentation = source.split("def plot_presentation", 1)[1].split("def run_presentation", 1)[0]
    assert "set_bbox" not in presentation
    assert '"label_white_background": False' in presentation
    assert '"label_step": .10' in presentation


def test_route_title_is_source_bound_and_legacy_defaults_to_direct():
    assert MODULE.density_route({}) == "direct_finite_q"
    assert MODULE.density_route({"density_route": "q0_lambda_reference"}) == "q0_lambda_reference"
    with pytest.raises(ValueError, match="density_route"):
        MODULE.density_route({"density_route": "folded"})


def test_fixed_comparison_colors_show_overflow_without_changing_values():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    rows = [{"T": T, "muB": mu, "status": "screened", "R": value}
            for (T, mu), value in zip(((90, 0), (90, 25), (100, 0), (100, 25)), (-.1, .3, .5, 1.2))]
    data = MODULE.neighborhood(rows, COEFFICIENTS)
    before = data["R"].copy()
    fig, display = MODULE.plot_presentation(data, COEFFICIENTS, 50,
                                             route="q0_lambda_reference", color_limits=(0, 1.15))
    assert display["color_limits"] == [0, 1.15]
    assert display["colorbar_extend"] == "both"
    assert display["under_range_nodes"] == display["over_range_nodes"] == 1
    assert any("q=0 extrapolated" in text.get_text() for text in fig.texts)
    np.testing.assert_array_equal(data["R"], before)
    plt.close(fig)


def test_clean_contours_keep_values_and_gaps_without_warning_artists():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    rows = [{"T": T, "muB": mu, "status": "screened", "R": .2 + .002*T + .001*mu}
            for T in (80, 90, 100, 110) for mu in (0, 25, 50, 75)]
    rows[0]["R"] = None
    data = MODULE.neighborhood(rows, COEFFICIENTS, band_MeV=50)
    before = data["R"].copy()
    fig, display = MODULE.plot_presentation(data, COEFFICIENTS, 50,
        route="q0_lambda_reference", clean_contours=True, prescription="finite_window")
    np.testing.assert_array_equal(data["R"], before)
    assert data["failed"].sum() == 1
    assert not display["failure_overlay"] and not display["warning_annotations"]
    legend_text = [text.get_text() for legend in fig.legends for text in legend.get_texts()]
    assert legend_text == ["Current chemical freeze-out"]
    assert any("finite window" in text.get_text() for text in fig.texts)
    plt.close(fig)
