from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
PLOTTER = ROOT / "scripts" / "analysis" / "relaxtime" / "plot_charged_gbu_contour_screening.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("charged_gbu_contour_screening_plot", PLOTTER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_merged_numerical_product_preserves_endpoint_diagnostics(tmp_path):
    module = _load_module()
    path = tmp_path / "merged.csv"
    row = {"T_MeV": 140, "muB_MeV": 425, "status": "screened",
           "K_plus_endpoint_warning": True, "K_plus_endpoint_warning_shells": 2,
           "K_plus_static_imaginary_max": 1.2e-6,
           "K_plus_density_prescription": "finite_window", "K_plus_failure_reason": ""}
    module._write_merged_csv({"rows": [row]}, path)
    with path.open(newline="", encoding="utf-8") as handle:
        saved = list(csv.DictReader(handle))[0]
    assert saved["status"] == "screened"
    assert saved["K_plus_endpoint_warning"] == "True"
    assert saved["K_plus_endpoint_warning_shells"] == "2"
    assert float(saved["K_plus_static_imaginary_max"]) == 1.2e-6
    assert saved["K_plus_density_prescription"] == "finite_window"


def _write_shard(root: Path, shard: int, rows: list[dict[str, str]]) -> None:
    shard_root = root / f"shard-{shard}"
    shard_root.mkdir(parents=True)
    fields = list(rows[0])
    with (shard_root / "contour_points.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "schema": "charged_gbu_contour_scan_v2",
        "point_count": len(rows),
        "T_grid": [40.0, 50.0],
        "muB_grid": [0.0, 50.0],
        "channels": ["pi_plus", "K_plus", "pi_minus", "K_minus"],
        "settings": {"mesh": 64},
        "config": "config/models/pnjl/charged_gbu_infinite_v1.toml",
        "git_head": "fixture",
        "source_hashes": {"fixture": "abc"},
    }
    (shard_root / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")


def _row(T: str, muB: str, *, status: str = "screened") -> dict[str, str]:
    values = {
        "T_MeV": T,
        "muB_MeV": muB,
        "status": status,
        "pi_plus_density_inv_fm3": "1.0" if status == "screened" else "",
        "K_plus_density_inv_fm3": "0.5" if status == "screened" else "",
        "Kplus_over_pi_plus": "0.5" if status == "screened" else "",
        "pi_minus_density_inv_fm3": "1.2" if status == "screened" else "",
        "K_minus_density_inv_fm3": "0.6" if status == "screened" else "",
        "Kminus_over_pi_minus": "0.5" if status == "screened" else "",
        "pi_plus_passed": "true",
        "K_plus_passed": "true",
        "pi_minus_passed": "true",
        "K_minus_passed": "true",
    }
    return values


def test_loader_preserves_failed_points_as_mask(tmp_path):
    module = _load_module()
    rows = [_row(str(T), str(muB)) for T in (40, 50) for muB in (0, 50)]
    rows[-1] = _row("50", "50", status="gate_failed")
    _write_shard(tmp_path, 0, rows)
    dataset = module.load_dataset(tmp_path)
    assert len(dataset["rows"]) == 4
    assert dataset["status_counts"] == {"screened": 3, "gate_failed": 1}
    values = module.matrix(dataset, "Kplus_over_pi_plus")
    assert values[1][1] is None
    assert module.mask_matrix(dataset)[1][1] == 2


def test_loader_rejects_duplicate_grid_keys(tmp_path):
    module = _load_module()
    rows = [_row("40", "0"), _row("40", "0"), _row("40", "50"), _row("50", "50")]
    _write_shard(tmp_path, 0, rows)
    with pytest.raises(ValueError, match="duplicate screening key"):
        module.load_dataset(tmp_path)


def test_plotter_is_solver_free_and_explicitly_no_interpolation():
    source = PLOTTER.read_text(encoding="utf-8")
    assert "solver-free" in source
    assert '"interpolation_policy": "none"' in source
    assert "failed rows remain masked" in source
    assert "zero-fill" in source


def test_reference_lines_keep_bqs_and_equal_flavor_coordinate_contract():
    module = _load_module()
    lines = module.load_reference_lines()
    assert set(lines) == {"freezeout", "crossover", "first_order_maxwell"}
    assert lines["freezeout"]["points"]
    assert lines["crossover"]["coordinate_convention"].startswith("input mu_MeV")
    assert "equal-flavor" in lines["first_order_maxwell"]["background_scope"]
    dataset = {"T_grid": [40.0, 50.0], "muB_grid": [0.0, 50.0]}
    x_limits, y_limits = module._axis_limits(dataset, {"freezeout": {"points": [{"T_MeV": 100.0, "muB_MeV": 200.0}]}})
    assert x_limits[1] == 200.0
    assert y_limits[1] == 100.0


def test_contours_block_every_cell_touching_a_failed_corner(tmp_path):
    module = _load_module()
    rows = [_row(str(T), str(muB)) for T in (40, 50) for muB in (0, 50)]
    rows[-1] = _row("50", "50", status="gate_failed")
    _write_shard(tmp_path, 0, rows)
    dataset = module.load_dataset(tmp_path)
    values, info = module.contour_payload(dataset, "Kplus_over_pi_plus")
    assert values.mask.tolist() == [[False, False], [False, True]]
    assert info["valid_cells"] == 0
    assert info["blocked_cells"] == 1
    assert info["corner_mask"] is False
    assert info["levels"] == []


def test_contour_levels_preserve_signed_values_and_use_no_smoothing(tmp_path):
    module = _load_module()
    rows = [_row(str(T), str(muB)) for T in (40, 50) for muB in (0, 50)]
    for row, value in zip(rows, (-1, 0, 1, 2)):
        row["Kplus_over_pi_plus"] = str(value)
    _write_shard(tmp_path, 0, rows)
    dataset = module.load_dataset(tmp_path)
    values, info = module.contour_payload(dataset, "Kplus_over_pi_plus")
    assert values.tolist() == [[-1.0, 0.0], [1.0, 2.0]]
    assert min(info["levels"]) < 0 < max(info["levels"])
    assert info["valid_cells"] == 1
    assert info["smoothing"] == "none"
    assert info["extrapolation"] == "none"


def test_wide_positive_density_levels_are_log_spaced_without_changing_data(tmp_path):
    module = _load_module()
    rows = [_row(str(T), str(muB)) for T in (40, 50) for muB in (0, 50)]
    for row, value in zip(rows, (1e-6, 1e-4, 1e-2, 1.0)):
        row["pi_plus_density_inv_fm3"] = str(value)
    _write_shard(tmp_path, 0, rows)
    dataset = module.load_dataset(tmp_path)
    values, info = module.contour_payload(dataset, "pi_plus_density_inv_fm3")
    assert values.tolist() == [[1e-6, 1e-4], [1e-2, 1.0]]
    assert len(info["levels"]) == 7
    assert "logarithmically spaced" in info["level_selection"]


def test_png_render_smoke_retains_failed_mask_and_emits_contour_figures(tmp_path):
    module = _load_module()
    rows = [_row(str(T), str(muB)) for T in (40, 50) for muB in (0, 50)]
    for row, value in zip(rows, (-0.1, 0.3, 0.5, 0.7)):
        row["Kplus_over_pi_plus"] = str(value)
    _write_shard(tmp_path / "input", 0, rows)
    output = tmp_path / "figures"
    assert module.main(["--input-root", str(tmp_path / "input"), "--output-dir", str(output), "--no-reference-lines"]) == 0
    manifest = json.loads((output / "plot_manifest.json").read_text(encoding="utf-8"))
    assert len(manifest["figures"]) == 7
    assert len(list(output.glob("*_contours.png"))) == 3
    assert manifest["contour_fields"]["Kplus_over_pi_plus"]["corner_mask"] is False
    assert manifest["manuscript_eligible"] is False
    assert manifest["value_summaries"]["Kplus_over_pi_plus"]["negative_count"] == 1
    assert manifest["source_git_sha"] == "fixture"
    assert "postprocess_git_sha" in manifest
    assert manifest["generator"]["sha256"] == module.sha256_file(PLOTTER)
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        module.main(["--input-root", str(tmp_path / "input"), "--output-dir", str(output)])


def test_contour_labels_skip_tiny_paths_and_stagger_long_ones():
    module = _load_module()
    segments = [[[[0, 0], [1, 0]]], [[[0, 0], [0, 1]]], [[[0, 0], [0.001, 0]]]]
    positions = module.contour_label_positions(segments, 1.0, 1.0)
    assert len(positions) == 2
    assert positions[0][0] == 0
    assert positions[1][0] == 1
    assert positions[0][1] == pytest.approx((0.15, 0.0))
    assert positions[1][1] == pytest.approx((0.0, 0.5))


def test_replot_workflow_skips_scan_and_keeps_source_and_render_runs_separate():
    source = (ROOT / ".github" / "workflows" / "relaxtime-charged-gbu-contour-scan.yml").read_text(encoding="utf-8")
    assert "if: ${{ inputs.replot_run_id == '' && !inputs.endpoint_audit && !inputs.q0_window_audit }}" in source
    assert "run-id: ${{ inputs.replot_run_id || github.run_id }}" in source
    assert "--postprocess-run-id" in source
    assert "always() && !inputs.endpoint_audit && !inputs.q0_window_audit && needs.prepare.result == 'success'" in source


def test_label_collision_suppression_changes_only_text_visibility():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    module = _load_module()
    fig, ax = plt.subplots()
    labels = [ax.text(0.2, 0.2, "0.8"), ax.text(0.2, 0.2, "1.6"), ax.text(0.8, 0.8, "2.4")]
    fig.canvas.draw()
    counts = module.suppress_overlapping_contour_labels(labels, fig.canvas.get_renderer())
    assert counts == {"visible": 2, "hidden_overlapping_text": 1}
    assert [label.get_visible() for label in labels] == [True, False, True]
    assert [label.get_text() for label in labels] == ["0.8", "1.6", "2.4"]
    plt.close(fig)


def test_loader_rejects_mixed_density_routes(tmp_path):
    module = _load_module()
    _write_shard(tmp_path, 0, [_row("40", "0"), _row("40", "50")])
    _write_shard(tmp_path, 1, [_row("50", "0"), _row("50", "50")])
    path = tmp_path / "shard-1/manifest.json"
    manifest = json.loads(path.read_text())
    manifest["density_route"] = "q0_lambda_reference"
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="density_route"):
        module.load_dataset(tmp_path)


def test_q0_workflow_uses_same_infinite_kernel_grid_and_new_scan():
    source = (ROOT / ".github/workflows/relaxtime-charged-gbu-contour-scan.yml").read_text()
    assert 'default: "direct_finite_q"' in source
    assert '--density-route "$DENSITY_ROUTE"' in source
    assert "--comparison-color-limits 0 1.15" in source
    assert "--presentation --band-MeV 50" in source
