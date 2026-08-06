import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "plot_pnjl_cep_hybrid_production_shadow.py"


def _module():
    spec = importlib.util.spec_from_file_location("cep_hybrid_plotter", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_local_y_bounds_preserve_weak_s_shape_scale():
    module = _module()
    low, high = module._local_y_bounds([(1.8, 291.2516), (1.9, 291.2622), (2.0, 291.2735)])
    assert high - low < 0.05
    assert low < 291.2516 < 291.2735 < high


def test_local_y_bounds_has_finite_fallback_for_missing_points():
    module = _module()
    assert module._local_y_bounds([]) == (0.0, 1.0)


def test_smooth_window_uses_longest_low_slope_run_without_phase_labels():
    module = _module()
    values = [
        0.0,
        10.0,
        10.5,
        11.0,
        11.5,
        12.0,
        22.0,
    ]
    rows = [
        {"xi": "0.0", "method": "production_hybrid", "T_MeV": "10.0", "rho": str(index), "muq_MeV": str(mu)}
        for index, mu in enumerate(values)
    ]
    low, high, metadata = module._smooth_window(rows, "0.0", 10.0)
    assert metadata["smooth_window_status"] == "selected"
    assert low == 1.0
    assert high == 5.0


def test_slice_row_matches_serialized_temperature():
    module = _module()
    rows = [{"xi": "0.0", "method": "production_hybrid", "T_MeV": "130.9619140625"}]
    row = module._slice_row(rows, "0.0", 130.9619140625)
    assert row is rows[0]
