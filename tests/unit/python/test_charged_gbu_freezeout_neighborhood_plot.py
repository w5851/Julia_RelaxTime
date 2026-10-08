from __future__ import annotations

import importlib.util
import csv
import json
from pathlib import Path
import sys

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


def test_low_ratio_levels_select_only_requested_supported_values():
    assert MODULE.select_low_ratio_levels(.002, .035, [.04, .005, .001, .02]) == [.005, .02]
    assert MODULE.select_low_ratio_levels(.001, .04, None) == []
    for invalid in ([0], [-.01], [.05], [.1], [float("nan")], [float("inf")], [.01, .01]):
        with pytest.raises(ValueError, match="low-ratio-levels"):
            MODULE.select_low_ratio_levels(0, 1, invalid)


def test_low_ratio_contours_keep_regular_paths_values_and_exact_small_labels():
    import math
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    rows = []
    for T in range(40, 221, 5):
        for mu in range(0, 751, 75):
            ratio = (.0006 * math.exp((T-40)/13) if T <= 90 else .028 + .006*(T-90)) * (1+mu/750)
            failed = (T, mu) == (70, 300)
            rows.append({"T": T, "muB": mu, "R": None if failed else ratio,
                         "status": "gate_failed" if failed else "screened"})
    data = MODULE.neighborhood(rows, COEFFICIENTS, band_MeV=50, T_bounds_MeV=(40, 220),
                               include_rectangle_MeV=(40, 130, 0, 750))
    before = data["R"].copy()
    options = {"route": "q0_lambda_reference", "color_limits": (0, 1.15),
               "clean_contours": True, "prescription": "finite_window"}
    base_fig, base = MODULE.plot_presentation(data, COEFFICIENTS, 50, **options)
    levels = [.002, .005, .01, .02, .03, .04]
    fig, display = MODULE.plot_presentation(data, COEFFICIENTS, 50, low_ratio_levels=levels, **options)
    assert display["low_ratio_levels"] == levels
    assert display["visible_low_ratio_labels"] == ["0.002", "0.005", "0.01", "0.02", "0.03", "0.04"]
    assert display["regular_contour_levels"] == base["contour_levels"]
    assert display["contour_levels"] == sorted(base["contour_levels"] + levels)
    assert display["contour_step"] is None and display["label_step"] is None
    assert display["regular_contour_step"] == .05 and display["regular_label_step"] == .10
    assert display["color_limits"] == base["color_limits"]
    base_contour, regular_contour, low_contour = base_fig.axes[0].collections[1], fig.axes[0].collections[1], fig.axes[0].collections[2]
    for old_path, path in zip(base_contour.get_paths(), regular_contour.get_paths()):
        np.testing.assert_array_equal(old_path.vertices, path.vertices)
        np.testing.assert_array_equal(old_path.codes, path.codes)
    assert [(t.get_text(), t.get_position()) for t in regular_contour.labelTexts] == [
        (t.get_text(), t.get_position()) for t in base_contour.labelTexts]
    assert low_contour._corner_mask is False
    np.testing.assert_array_equal(low_contour.levels, levels)
    np.testing.assert_array_equal(data["R"], before)
    assert data["failed"].sum() == 1
    assert not display["failure_overlay"] and not display["warning_annotations"]
    assert any("finer below" in text.get_text() for text in fig.texts)
    plt.close(base_fig)
    plt.close(fig)


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


def test_rectangle_adds_only_requested_native_nodes_and_keeps_failed_corners():
    import numpy as np

    rows = [{"T": T, "muB": mu, "status": "screened", "R": -.2 if T == 70 else .5}
            for T in (40, 70, 100, 130, 160, 190, 220) for mu in (0, 25, 50, 75)]
    rows[1].update(status="gate_failed", R=None)
    original = MODULE.neighborhood(rows, COEFFICIENTS, T_bounds_MeV=(40, 220))
    data = MODULE.neighborhood(rows, COEFFICIENTS, T_bounds_MeV=(40, 220),
                               include_rectangle_MeV=(40, 100, 0, 50))
    assert data["selected"].sum() == 10
    assert data["rectangle_selected"].sum() == 9
    assert data["added"].sum() == 6
    assert data["failed"].sum() == 1
    assert np.isnan(data["R"][0, 1])
    assert data["R"][1, 2] == -.2
    assert np.isnan(data["R"][:2, 3]).all()
    assert np.isnan(data["R"][3:]).all()
    assert not data["valid_cells"][0].any()
    assert data["valid_cells"].sum() == 2
    np.testing.assert_array_equal(data["R"][original["selected"]], original["R"][original["selected"]])
    np.testing.assert_array_equal(data["band_selected"], original["selected"])


@pytest.mark.parametrize("rectangle", [(130, 40, 0, 750), (40, 130, 750, 0),
                                       (39, 130, 0, 750), (40, 130, 0, 775),
                                       (40, float("nan"), 0, 750), (40, 130, 0)])
def test_rectangle_rejects_invalid_or_unavailable_bounds(rectangle):
    rows = [{"T": T, "muB": mu, "status": "screened", "R": .5}
            for T in (40, 100, 220) for mu in (0, 250, 500, 750)]
    with pytest.raises(ValueError, match="include-rectangle-MeV"):
        MODULE.neighborhood(rows, COEFFICIENTS, T_bounds_MeV=(40, 220),
                             include_rectangle_MeV=rectangle)


def test_extended_presentation_legend_does_not_overlap_xlabel_or_footer():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = [{"T": T, "muB": mu, "status": "screened", "R": .002*T + .001*mu}
            for T in (40, 70, 100, 130, 160, 190, 220) for mu in (0, 250, 500, 750)]
    data = MODULE.neighborhood(rows, COEFFICIENTS, band_MeV=50, T_bounds_MeV=(40, 220),
                               include_rectangle_MeV=(40, 130, 0, 750))
    fig, _ = MODULE.plot_presentation(data, COEFFICIENTS, 50, route="q0_lambda_reference",
                                      color_limits=(0, 1.15), clean_contours=True,
                                      prescription="finite_window")
    fig.set_dpi(220)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    legend = fig.legends[0].get_window_extent(renderer)
    assert not legend.overlaps(fig.axes[0].xaxis.label.get_window_extent(renderer))
    assert all(not legend.overlaps(text.get_window_extent(renderer)) for text in fig.texts)
    plt.close(fig)


def test_extended_presentation_exports_original_fields_and_reproducible_inputs(tmp_path, monkeypatch):
    csv_path = tmp_path / "points.csv"
    profile_path = tmp_path / "profile.toml"
    source_path = tmp_path / "manifest.json"
    output_dir = tmp_path / "extended"
    Ts, mus = [40, 70, 100, 130, 160, 190, 220], [0, 25, 50, 75]
    rows = []
    for T in Ts:
        for mu in mus:
            failed = (T, mu) == (40, 25)
            ratio = "" if failed else f"{.002*T + .001*mu:.12f}"
            rows.append({"T_MeV": str(T), "muB_MeV": str(mu),
                         "status": "gate_failed" if failed else "screened",
                         "pi_plus_density_inv_fm3": "1.0000000000000000",
                         "K_plus_density_inv_fm3": ratio, "Kplus_over_pi_plus": ratio,
                         "pi_plus_passed": "True", "K_plus_passed": "False" if failed else "True",
                         "K_plus_endpoint_warning": "false" if failed else "true",
                         "K_plus_failure_reason": "q0 reference onset is not positive real" if failed else "",
                         "preserved_diagnostic": "0007"})
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    profile_path.write_text("[coefficients]\n" + "\n".join(f"{key} = {value}" for key, value in COEFFICIENTS.items()),
                            encoding="utf-8")
    source_path.write_text(json.dumps({
        "merged_csv": {"sha256": MODULE.digest(csv_path)},
        "reference_lines": {"freezeout": {"source_sha256": MODULE.digest(profile_path)}},
        "grid": {"T_MeV": Ts, "muB_MeV": mus}, "row_count": len(rows),
        "status_counts": {"screened": 27, "gate_failed": 1},
        "source_run_id": "fixture", "source_git_sha": "fixture", "background_contract": {},
        "density_route": "q0_lambda_reference", "density_prescription": "finite_window"}), encoding="utf-8")
    before = [path.read_bytes() for path in (csv_path, source_path, profile_path)]
    monkeypatch.setattr(sys, "argv", [str(SCRIPT), "--csv", str(csv_path),
        "--source-manifest", str(source_path), "--freezeout-profile", str(profile_path),
        "--presentation", "--band-MeV", "20", "--include-rectangle-MeV", "40", "100", "0", "50",
        "--clean-contours", "--comparison-color-limits", "0", "1.15", "--output-dir", str(output_dir)])
    assert MODULE.main() == 0
    manifest = json.loads((output_dir / "plot_manifest.json").read_text(encoding="utf-8"))
    extension = manifest["coverage_extension"]
    assert extension["requested_nodes"] == 9
    assert extension["requested_screened_nodes"] == 8
    assert extension["requested_failed_nodes"] == 1
    assert extension["added_display_nodes"] == 6
    assert extension["added_screened_nodes"] == 5
    assert extension["newly_computed_nodes"] == extension["missing_native_grid_nodes"] == 0
    assert manifest["selected_nodes"] == 10 and manifest["screened_nodes"] == 9
    assert manifest["failed_nodes"] == 1
    assert " OR " in manifest["selection_rule"]
    assert manifest["warning_display"] == "none"
    assert not manifest["rendering"]["failure_overlay"]
    for record, expected_rows in zip(manifest["data_tables"],
            ([rows[i] for i in (0, 1, 2, 4, 5, 6, 8, 9, 10, 11)],
             [rows[i] for i in (0, 1, 2, 4, 5, 6)],
             [rows[i] for i in (0, 1, 2, 4, 5, 6, 8, 9, 10)])):
        path = Path(record["path"])
        with path.open(newline="", encoding="utf-8") as handle:
            assert list(csv.DictReader(handle)) == expected_rows
        assert record["row_count"] == len(expected_rows)
        assert MODULE.digest(path) == record["sha256"]
    for record in manifest["source_snapshots"]:
        assert Path(record["path"]).read_bytes() == Path(record["source_path"]).read_bytes()
        assert MODULE.digest(Path(record["path"])) == record["sha256"]
    assert [path.read_bytes() for path in (csv_path, source_path, profile_path)] == before
    assert MODULE.digest(Path(manifest["output"]["path"])) == manifest["output"]["sha256"]
    with pytest.raises(FileExistsError, match="refusing to overwrite existing case"):
        MODULE.main()
