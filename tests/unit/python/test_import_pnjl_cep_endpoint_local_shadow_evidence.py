import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "import_pnjl_cep_endpoint_local_shadow_evidence.py"


def _module():
    spec = importlib.util.spec_from_file_location("endpoint_local_evidence_import", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_csv(path: Path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _source(tmp_path: Path) -> Path:
    source = tmp_path / "source"
    figures = source / "figures"
    figures.mkdir(parents=True)
    methods = ("production_hybrid", "memoized_dense", "independent_oracle")
    xis = (-0.5, 0.0, 0.5)
    temperatures = tuple(float(value) for value in range(8))
    slices = []
    curves = []
    endpoint_index = 0
    for method in methods:
        for xi in xis:
            for temperature in temperatures:
                certificate = "none"
                if method == "production_hybrid" and endpoint_index < 3:
                    certificate = (
                        "endpoint_limited_first_order"
                        if endpoint_index == 0
                        else "endpoint_local_geometry_first_order"
                    )
                    endpoint_index += 1
                slices.append({
                    "xi": xi, "method": method, "T_MeV": temperature,
                    "result_status": "confirmed_first_order",
                    "certificate_type": certificate,
                    "support_low": 0.0 if certificate != "none" else "",
                    "support_high": 0.00625 if certificate != "none" else "",
                    "endpoint_anchor_rho": 0.003125 if certificate != "none" else "",
                    "endpoint_left_bracket_low": 0.0 if certificate != "none" else "",
                    "endpoint_left_bracket_high": 0.00625 if certificate != "none" else "",
                    "endpoint_refinement_count": 12 if certificate != "none" else 0,
                })
                curves.append({
                    "xi": xi, "method": method, "T_MeV": temperature,
                    "rho": 1.0, "muq_MeV": 300.0,
                    "converged": "true", "finite": "true",
                })
    _write_csv(source / "slice_metrics.csv", slices)
    _write_csv(source / "geometry_accuracy.csv", slices)
    _write_csv(source / "curve_points.csv", curves)
    for name in ("method_costs.csv", "actions_costs.csv", "cep_accuracy.csv", "curve_index.csv"):
        _write_csv(source / name, [{"name": name, "value": 1}])
    (figures / "example.png").write_bytes(b"png")
    (figures / "plot_manifest.json").write_text(json.dumps({
        "schema_version": "cep_maxwell_endpoint_local_production_shadow_v4",
        "figures": [{"file": "example.png", "xi": "0.0", "T_MeV": 5.0, "reason": "test"}],
    }), encoding="utf-8")
    hashes = {
        name: hashlib.sha256((source / name).read_bytes()).hexdigest()
        for name in (
            "slice_metrics.csv", "geometry_accuracy.csv", "curve_points.csv",
            "method_costs.csv", "actions_costs.csv", "cep_accuracy.csv", "curve_index.csv",
        )
    }
    manifest = {
        "schema_version": "cep_maxwell_endpoint_local_production_shadow_v4",
        "scope": "full", "evidence_state": "final",
        "expected_calculation_sha": "a" * 40, "postprocess_sha": "b" * 40,
        "source_run_id": "123", "run_id": "123",
        "oracle_overlay": {"deep_run_id": "456"},
        "gate": {
            "verdict": "full_hybrid_candidate",
            "workflow_contract_errors": [], "oracle_errors": [],
            "classification_errors": [], "endpoint_errors": [],
            "coverage_errors": [], "performance_errors": [],
        },
        "file_sha256": hashes,
    }
    (source / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return source


def test_import_keeps_raw_curves_external_and_records_dual_sha(tmp_path):
    module = _module()
    source = _source(tmp_path)
    output = tmp_path / "evidence"
    manifest = module.import_evidence(source, output, "789")
    assert manifest["verdict"] == "full_hybrid_candidate"
    assert manifest["aggregate_replay_run_id"] == "789"
    assert manifest["calculation_sha"] == "a" * 40
    assert manifest["postprocess_sha"] == "b" * 40
    assert manifest["raw_curve_points"]["raw_curve_copy_in_repository"] is False
    assert not (output / "curve_points.csv").exists()
    assert not (output / "tables" / "curve_points.csv").exists()
    assert (output / "tables" / "slice_metrics.csv").is_file()
    assert (output / "figures" / "example.png").is_file()


def test_import_rejects_unexpected_nan_sentinel(tmp_path):
    module = _module()
    source = _source(tmp_path)
    (source / "method_costs.csv").write_text("name,value\nmethod_costs.csv,NaN\n", encoding="utf-8")
    with pytest.raises(ValueError, match="unexpected NaN sentinel"):
        module._validate_derived_tables(source)
