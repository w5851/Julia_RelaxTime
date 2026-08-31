import csv
import hashlib
import importlib.util
import io
import json
from pathlib import Path
from email.message import Message
from urllib.error import HTTPError

import pytest


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "pnjl_c2_targeted_closure.py"
WORKFLOW = ROOT / "docs" / "analysis" / "governance" / "diagnostic_workflow_retirement_wave2_v1" / "definitions" / "pnjl-c2-targeted-closure-v1.yml"


def load_module():
    spec = importlib.util.spec_from_file_location("c2_targeted_closure", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_csv(path: Path, fieldnames, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_method(root: Path, module, target_id: str, method: str, *, status="confirmed_first_order"):
    target_dir = root / target_id / method
    target_dir.mkdir(parents=True)
    xi, temperature = module.REGRESSION_TARGETS[target_id]
    curve = target_dir / "curve_points.csv"
    _write_csv(
        curve,
        ["xi", "T_MeV", "rho", "muq_MeV", "residual_norm", "converged", "finite"],
        [
            {"xi": xi, "T_MeV": temperature, "rho": 0.0, "muq_MeV": 300.0, "residual_norm": 0.0, "converged": "true", "finite": "true"},
            {"xi": xi, "T_MeV": temperature, "rho": 0.003125, "muq_MeV": 301.0, "residual_norm": 0.0, "converged": "true", "finite": "true"},
        ],
    )
    diagnostics = target_dir / "slice_diagnostics.csv"
    _write_csv(
        diagnostics,
        ["result_status", "maxwell_candidate_count", "maxwell_crossing_count", "geometry_converged", "position_error_MeV", "density_error", "area_residual", "finite_and_converged"],
        [{"result_status": status, "maxwell_candidate_count": 1, "maxwell_crossing_count": 3, "geometry_converged": "true", "position_error_MeV": 0.0, "density_error": 0.0, "area_residual": 0.0, "finite_and_converged": "true"}],
    )
    accuracy = target_dir / "accuracy.csv"
    _write_csv(accuracy, ["result_status"], [{"result_status": status}])
    costs = target_dir / "method_costs.csv"
    _write_csv(costs, ["unique_solves", "point_requests", "cache_hits", "failed_points", "runner_seconds"], [{"unique_solves": 2, "point_requests": 2, "cache_hits": 0, "failed_points": 0, "runner_seconds": 1.0}])
    files = {name: module.sha256(target_dir / name) for name in ("curve_points.csv", "slice_diagnostics.csv", "accuracy.csv", "method_costs.csv")}
    manifest = {
        "schema_version": module.JOB_SCHEMA,
        "scope": "regression_curves",
        "target_id": target_id,
        "target": {"id": target_id, "xi": xi, "T_MeV": temperature},
        "method": method,
        "source_run_id": "123",
        "calculation_sha": module.CALCULATION_SHA,
        "postprocess_sha": "a" * 40,
        "solver_called": True,
        "files": files,
        "provenance": {"reference_write": False, "oracle_labels_used_for_routing": False},
    }
    (target_dir / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")


def test_workflow_freezes_scope_sha_and_exact_target_matrix():
    module = load_module()
    text = WORKFLOW.read_text(encoding="utf-8")
    assert module.SCHEMA_VERSION in text
    assert module.CALCULATION_SHA in text
    assert "regression_curves" in text and "cep_brackets" in text and "cross_axis_audit" in text
    assert len(module.artifact_names("regression_curves")) == 9
    assert len(module.artifact_names("cep_brackets")) == 3
    assert "oracle_labels_used_for_routing" in SCRIPT.read_text(encoding="utf-8")


def test_source_validation_reconciles_two_methods_without_solver_call(tmp_path):
    module = load_module()
    target_id = next(iter(module.REGRESSION_TARGETS))
    _write_method(tmp_path, module, target_id, "production_hybrid")
    _write_method(tmp_path, module, target_id, "independent_oracle")
    result = module.validate_source(tmp_path, module.CALCULATION_SHA, "a" * 40, "123")
    assert result["manifest_count"] == 2
    assert result["failed_points"] == 0
    assert result["solver_called"] is True


def test_source_validation_rejects_duplicate_raw_curve_key(tmp_path):
    module = load_module()
    target_id = next(iter(module.REGRESSION_TARGETS))
    _write_method(tmp_path, module, target_id, "production_hybrid")
    _write_method(tmp_path, module, target_id, "independent_oracle")
    curve = tmp_path / target_id / "production_hybrid" / "curve_points.csv"
    text = curve.read_text(encoding="utf-8")
    curve.write_text(text + text.splitlines()[-1] + "\n", encoding="utf-8")
    manifest = json.loads((tmp_path / target_id / "production_hybrid" / "manifest.json").read_text(encoding="utf-8"))
    manifest["files"]["curve_points.csv"] = module.sha256(curve)
    (tmp_path / target_id / "production_hybrid" / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate raw curve key"):
        module.validate_source(tmp_path, module.CALCULATION_SHA, "a" * 40, "123")


def test_aggregate_keeps_classification_mismatch_diagnostic(tmp_path):
    module = load_module()
    target_id = next(iter(module.REGRESSION_TARGETS))
    for current_id in module.REGRESSION_TARGETS:
        for method in ("production_hybrid", "independent_oracle"):
            status = "confirmed_first_order"
            if current_id == target_id and method == "production_hybrid":
                status = "ambiguous_near_critical"
            _write_method(tmp_path / "source", module, current_id, method, status=status)
    manifest = module.aggregate_validate(tmp_path / "source", tmp_path / "aggregate", "regression_curves", module.CALCULATION_SHA, "a" * 40, "123")
    assert manifest["verdict"] == "targeted_classification_inconclusive"
    overlay = (tmp_path / "aggregate" / "classification_overlay.csv").read_text(encoding="utf-8")
    assert "classification_match" in overlay
    assert manifest["solver_called"] is False


def test_safe_extract_rejects_path_traversal(tmp_path):
    module = load_module()
    import zipfile

    archive = tmp_path / "bad.zip"
    with zipfile.ZipFile(archive, "w") as handle:
        handle.writestr("../escape.txt", "bad")
    with pytest.raises(ValueError, match="unsafe artifact path"):
        module._safe_extract(archive, tmp_path / "out")


def test_artifact_zip_redirect_drops_api_authorization(monkeypatch):
    module = load_module()
    api_url = "https://api.github.com/repos/example/repo/actions/artifacts/1/zip"
    signed_url = "https://artifact-storage.example/signed.zip"
    headers = Message()
    headers["Location"] = signed_url

    class FakeOpener:
        def open(self, request, timeout):
            assert request.headers["Authorization"] == "Bearer secret"
            raise HTTPError(api_url, 302, "Found", headers, io.BytesIO())

    class FakeResponse:
        def __enter__(self):
            return self

        def __exit__(self, *_):
            return False

        def read(self):
            return b"zip-bytes"

    monkeypatch.setattr(module.urllib.request, "build_opener", lambda *_: FakeOpener())

    def fake_urlopen(request, timeout):
        assert request.full_url == signed_url
        assert "Authorization" not in request.headers
        return FakeResponse()

    monkeypatch.setattr(module.urllib.request, "urlopen", fake_urlopen)
    assert module._download_artifact_zip(api_url, "secret") == b"zip-bytes"
