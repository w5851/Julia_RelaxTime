import csv
import hashlib
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "pnjl_stagec_density_certificate_feasibility.py"
WORKFLOW = ROOT / ".github" / "workflows" / "pnjl-stagec-density-certificate-feasibility.yml"


def load_module():
    spec = importlib.util.spec_from_file_location("stagec_density_feasibility", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_anchor_matrix_and_caps_are_fixed():
    module = load_module()
    assert len(module.DENSITY_XIS) == 10
    assert len(module.METHODS) == 3
    assert module.EXPECTED_JOB_COUNT == 30
    assert module.CAPS == (12, 16, 24)
    assert "STAGE_C_STEP = 0.003125" in (ROOT / "scripts" / "analysis" / "pnjl_stagec_density_certificate_job.jl").read_text(encoding="utf-8")
    assert module.DENSITY_ANCHORS + module.FIRST_ORDER_CONTROLS + module.MONOTONE_CONTROLS == module.ALL_ANCHORS


def test_stage_b_s_shape_candidate_and_deterministic_feature_selection():
    module = load_module()
    points = [(0.0, 300.0), (1.0, 301.0), (2.0, 303.0), (3.0, 302.0), (4.0, 301.0), (5.0, 302.0), (6.0, 304.0)]
    candidates = module._candidate_intervals(points)
    assert len(candidates) == 1
    assert candidates[0]["rho_low"] == 2.0
    assert candidates[0]["rho_high"] == 4.0

    targets = module._feature_targets(points)
    pool = [(rho + 0.1, mu + 0.01) for rho, mu in points]
    first = module._select_local_points(pool, targets, 4)
    second = module._select_local_points(pool, targets, 4)
    assert first == second
    assert len(first) == 4
    assert [row[0] for row in first] == sorted(row[0] for row in first)


def test_oracle_ambiguous_is_recorded_as_unsupported_confirmation():
    module = load_module()
    slices = []
    curves = []
    costs = [{"method": "memoized_dense", "unique_solves": "100000"}]
    for index, (xi, temperature) in enumerate(module.ALL_ANCHORS):
        is_monotone = (xi, temperature) in module.MONOTONE_CONTROLS
        values = [300.0 + rho if is_monotone else value for rho, value in (
            (0.0, 300.0), (1.0, 301.0), (2.0, 303.0), (3.0, 302.0), (4.0, 301.0), (5.0, 302.0), (6.0, 304.0)
        )]
        for rho, mu in zip((0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0), values):
            curves.append({
                "method": "independent_oracle", "xi": str(xi), "T_MeV": str(temperature),
                "rho_level": "0", "rho": str(rho), "muq_MeV": str(mu),
                "residual_norm": "0.0", "converged": "true", "finite": "true",
                "calculation_sha": module.SOURCE_CALCULATION_SHA,
            })
        curves.append({
            "method": "independent_oracle", "xi": str(xi), "T_MeV": str(temperature),
            "rho_level": "1", "rho": "2.5", "muq_MeV": "302.5",
            "residual_norm": "0.0", "converged": "true", "finite": "true",
            "calculation_sha": module.SOURCE_CALCULATION_SHA,
        })
        status = "confirmed_monotone" if is_monotone else "confirmed_first_order"
        slices.append({
            "method": "production_hybrid", "xi": str(xi), "T_MeV": str(temperature),
            "stage_b_status": status, "geometry_converged": "true",
            "position_error_MeV": "0.001", "density_error": "0.001",
            "maxwell_area_gate": "0.000001", "area_residual": "0.000001",
        })
        slices.append({
            "method": "independent_oracle", "xi": str(xi), "T_MeV": str(temperature),
            "result_status": "ambiguous_near_critical" if index == 0 else status,
        })
    bundle = {"curves": curves, "slices": slices, "costs": costs}
    frontier, _rows, _selected = module._evaluate_route(bundle, "balanced_density_features_v2", 12, 100000)
    assert frontier["unsupported_confirmations"] >= 1
    assert not frontier["state_gate"]


def test_validation_rejects_duplicate_nan_and_wrong_provenance():
    module = load_module()
    manifests = [
        {"method": method, "xi": xi, "calculation_sha": module.SOURCE_CALCULATION_SHA, "postprocess_sha": "b" * 40}
        for xi in module.DENSITY_XIS
        for method in module.METHODS
    ]
    curves = [
        {"method": "production_hybrid", "xi": "-0.5", "T_MeV": "60", "rho_level": "0", "rho": "0", "muq_MeV": "300", "residual_norm": "0", "converged": "true", "finite": "true", "calculation_sha": module.SOURCE_CALCULATION_SHA},
    ]
    valid = module._validate_bundle(curves, [{"method": "production_hybrid"}], [{"method": "memoized_dense", "unique_solves": "1"}], manifests, module.SOURCE_CALCULATION_SHA, "b" * 40, job_count=30)
    assert valid == []

    duplicate = module._validate_bundle(curves + [dict(curves[0])], [], [], manifests, module.SOURCE_CALCULATION_SHA, "b" * 40, job_count=30)
    assert any("duplicate curve key" in error for error in duplicate)

    invalid = dict(curves[0])
    invalid["muq_MeV"] = "nan"
    invalid_errors = module._validate_bundle([invalid], [], [], manifests, module.SOURCE_CALCULATION_SHA, "c" * 40, job_count=30)
    assert any("postprocess SHA mismatch" in error for error in invalid_errors)
    assert any("invalid curve coordinate" in error for error in invalid_errors)


def test_bundle_loader_checks_declared_file_hashes(tmp_path):
    module = load_module()
    for index, xi in enumerate(module.DENSITY_XIS):
        for method_index, method in enumerate(module.METHODS):
            job = tmp_path / f"job-{index}-{method}"
            job.mkdir()
            rows = {
                "curve_points.csv": [{
                    "method": method, "xi": xi, "T_MeV": 1.0, "rho_level": 0, "rho": index * 10 + method_index,
                    "muq_MeV": 300.0, "residual_norm": 0.0, "converged": "true", "finite": "true",
                    "calculation_sha": module.SOURCE_CALCULATION_SHA,
                }],
                "slice_metrics.csv": [{"method": method, "xi": xi, "T_MeV": 1.0}],
                "method_costs.csv": [{"method": method, "unique_solves": 1}],
                "cep_accuracy.csv": [{"method": method, "xi": xi, "anchor_T_MeV": 1.0}],
            }
            for name, values in rows.items():
                with (job / name).open("w", newline="", encoding="utf-8") as handle:
                    writer = csv.DictWriter(handle, fieldnames=list(values[0]))
                    writer.writeheader()
                    writer.writerows(values)
            hashes = {name: hashlib.sha256((job / name).read_bytes()).hexdigest() for name in rows}
            (job / "manifest.json").write_text(json.dumps({
                "method": method, "xi": xi, "calculation_sha": module.SOURCE_CALCULATION_SHA,
                "postprocess_sha": "b" * 40, "files": hashes,
            }), encoding="utf-8")

    broken = next(tmp_path.rglob("curve_points.csv"))
    broken.write_text(broken.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    try:
        module._load_bundle(tmp_path, module.SOURCE_CALCULATION_SHA, "b" * 40)
    except ValueError as error:
        assert "manifest hash mismatch" in str(error)
    else:
        raise AssertionError("tampered artifact was accepted")


def test_workflow_has_dual_sha_matrix_and_replay_contract():
    module = load_module()
    text = WORKFLOW.read_text(encoding="utf-8")
    assert "calculation_ref" in text
    assert "source_run_id" in text and "source_calculation_sha" in text
    assert "expected-source-postprocess-sha" in text
    assert "git worktree add --detach calculation" in text
    assert "matrix.xi" in text and "matrix.method" in text
    assert "max-parallel: 10" in text
    assert "test \"$numeric_jobs\" = 30" in text
    assert "aggregate_replay" in text
    assert "gh run download" in text
    assert "matplotlib==3.9.2" in text
    assert "solver_called" in SCRIPT.read_text(encoding="utf-8")
    assert "run_production_phase_pipeline" not in SCRIPT.read_text(encoding="utf-8")
