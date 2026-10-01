from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "relaxtime" / "select_charged_gbu_contour_candidates.py"


def _load():
    spec = importlib.util.spec_from_file_location("charged_gbu_contour_candidates", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _rows():
    rows = []
    for T in (40.0, 50.0, 60.0):
        for muB in (0.0, 50.0, 100.0):
            rows.append({
                "T_MeV": T,
                "muB_MeV": muB,
                "status": "screened",
                "Kplus_over_pi_plus": 0.2 + 0.01 * T / 10 + 0.001 * muB,
            })
    rows[-1]["status"] = "gate_failed"
    rows[-1]["Kplus_over_pi_plus"] = None
    return rows


def _refs():
    return {
        "freezeout": {"points": [{"sqrt_s_NN_GeV": 7.7, "T_MeV": 50.0, "muB_MeV": 50.0}]},
        "crossover": {"points": [{"T_MeV": 50.0, "muB_MeV": 100.0}], "source": "fixture"},
        "first_order_maxwell": {"points": [{"T_MeV": 200.0, "muB_MeV": 1000.0}], "source": "fixture"},
    }


def test_selector_preserves_failed_points_as_non_dispatch_and_records_categories():
    module = _load()
    candidates, summary = module.select_candidates(_rows(), [40.0, 50.0, 60.0], [0.0, 50.0, 100.0], _refs(), max_gradient=2, max_freezeout=1, max_phase=2)
    assert candidates
    assert all(row["status"] == "screened" for row in candidates)
    assert all(":" in row["production_point"] for row in candidates)
    assert any(any("freezeout" in origin for origin in row["provenance"]) for row in candidates)
    assert any(row["category"] == "mask_boundary" for row in candidates)
    assert summary["failed_count"] == 1
    assert summary["phase_reference_covered_anchor_counts"]["first_order_maxwell"] == 0


def test_selector_does_not_cross_mask_for_gradient():
    module = _load()
    rows = _rows()
    # Mask the entire T=50 row; no gradient may be estimated through it.
    for row in rows:
        if row["T_MeV"] == 50.0:
            row["status"] = "gate_failed"
            row["Kplus_over_pi_plus"] = None
    candidates, _ = module.select_candidates(rows, [40.0, 50.0, 60.0], [0.0, 50.0, 100.0], {}, max_gradient=3, max_freezeout=1, max_phase=1)
    assert not any(row["category"] == "gradient_extremum" for row in candidates)


def test_candidate_script_has_diagnostic_boundary():
    source = SCRIPT.read_text(encoding="utf-8")
    assert "solver-free" in source
    assert "never dispatched" in source
    assert '"production_default": False' in source


def test_actions_chain_candidate_selection_to_sparse_full_gate():
    scan = (ROOT / ".github" / "workflows" / "relaxtime-charged-gbu-contour-scan.yml").read_text(encoding="utf-8")
    full = (ROOT / ".github" / "workflows" / "relaxtime-charged-gbu-contour-full-gate.yml").read_text(encoding="utf-8")
    assert "select_charged_gbu_contour_candidates.py" in scan
    assert "charged-gbu-contour-${{ inputs.output_tag || format('run-{0}', github.run_id) }}-candidates" in scan
    assert "BENCH_CGBU_WORKLOAD: production" in full
    assert "full_smoke_production_gates" in full or "full production-gate point" in full
    assert "production_default=false" in full
