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
