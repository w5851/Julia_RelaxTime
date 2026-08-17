import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[3]
MODULE_PATH = REPO_ROOT / "scripts" / "analysis" / "pnjl_raw_curve_archive.py"
SPEC = importlib.util.spec_from_file_location("pnjl_raw_curve_archive", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def _write_curve(path: Path, xi: float, temperature: float) -> None:
    rows = []
    for rho in reversed(MODULE._rho_grid()):
        rows.append({
            "T_MeV": f"{temperature:.8f}",
            "rho": f"{rho:.8f}",
            "xi": f"{xi:.8f}",
            "mu_u_MeV": "100.0",
            "mu_d_MeV": "100.0",
            "mu_s_MeV": "120.0",
            "mu_avg_MeV": "106.6666666667",
            "mu_B_MeV": "300.0",
            "mu_Q_MeV": "0.0",
            "mu_S_MeV": "-20.0",
            "pressure_fm4": "1.0",
            "entropy_fm3": "0.5",
            "energy_fm4": "2.0",
            "rho_u_fm3": "0.1",
            "rho_d_fm3": "0.1",
            "rho_s_fm3": "0.1",
            "phi_u": "-0.1",
            "phi_d": "-0.1",
            "phi_s": "-0.2",
            "Phi1": "0.2",
            "Phi2": "0.3",
            "M_u_MeV": "300.0",
            "M_d_MeV": "300.0",
            "M_s_MeV": "500.0",
            "iterations": "4",
            "residual_norm": "1e-10",
            "converged": "true",
            "message": "oracle",
        })
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MODULE.RAW_CURVE_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _write_grid(path: Path) -> None:
    path.write_text(
        "axis,xi,T_MeV,level\n"
        "rho,-0.5,1.0,0\n"
        "rho,-0.5,147.0625,1\n"
        "rho,0.0,1.0,0\n"
        "temperature,-0.5,,0\n"
        "xi,,,-1\n",
        encoding="utf-8",
    )


def _build_sources(tmp_path: Path):
    source = tmp_path / "sources"
    source_postprocess = "a" * 40
    workflow = "b" * 40
    for shard_id, xi, temperatures in (
        ("xi_m0p5", -0.5, [1.0, 147.0625]),
        ("xi_0", 0.0, [1.0]),
    ):
        process_root = source / shard_id
        for temperature in temperatures:
            _write_curve(
                process_root
                / "production_eval"
                / f"prod_eval_xi_{MODULE._token(xi)}_T_{MODULE._token(temperature)}.csv",
                xi,
                temperature,
            )
        MODULE.build_source_manifest(
            process_root,
            process_root / MODULE.DENSE_SOURCE_MANIFEST_NAME,
            source_run_id="raw-run",
            source_artifact_name=f"raw-artifact-{shard_id}",
            shard_id=shard_id,
            xi=xi,
            calculation_sha="c" * 40,
            source_postprocess_sha=source_postprocess,
            source_workflow_sha=workflow,
            audit_workflow_sha=workflow,
            expected_temperatures=temperatures,
            source_grid_run_id="31862752226",
            source_grid_artifact_name="c2-aggregate",
            source_grid_manifest_path="reference/manifest.json",
            source_grid_manifest_sha256="d" * 64,
            audit_run_id="31941614867",
        )
    grid = tmp_path / "grid.csv"
    _write_grid(grid)
    coverage = tmp_path / "expected_coverage.json"
    MODULE.build_expected_coverage(
        grid,
        coverage,
        source_grid_run_id="31862752226",
        source_grid_artifact_name="c2-aggregate",
        source_grid_manifest_path="reference/manifest.json",
        source_grid_manifest_sha256="d" * 64,
        calculation_sha="c" * 40,
        source_postprocess_sha=source_postprocess,
        audit_run_id="31941614867",
    )
    return source, coverage, source_postprocess, workflow


def test_coverage_is_exact_irregular_set_and_not_cartesian(tmp_path):
    source, coverage_path, _, _ = _build_sources(tmp_path)
    coverage = json.loads(coverage_path.read_text(encoding="utf-8"))
    assert coverage["expected_curve_count"] == 3
    assert coverage["coordinate_policy"].startswith("exact set")
    assert coverage["temperatures_by_xi"]["-0.5"] == [1.0, 147.0625]
    matrix_path = tmp_path / "matrix.json"
    MODULE.build_expected_coverage(
        tmp_path / "grid.csv",
        tmp_path / "coverage2.json",
        source_grid_run_id="31862752226",
        source_grid_artifact_name="c2-aggregate",
        source_grid_manifest_path="reference/manifest.json",
        source_grid_manifest_sha256="d" * 64,
        calculation_sha="c" * 40,
        source_postprocess_sha="a" * 40,
        matrix_output=matrix_path,
    )
    matrix = json.loads(matrix_path.read_text(encoding="utf-8"))["include"]
    assert len(matrix) == 2
    assert {entry["xi"] for entry in matrix} == {"-0.5", "0"}
    assert source.is_dir()


def test_build_validate_and_restore_are_byte_preserving(tmp_path):
    source, coverage_path, source_postprocess, workflow = _build_sources(tmp_path)
    archive_dir = tmp_path / "archive"
    manifest = MODULE.build_archive(
        source,
        archive_dir,
        expected_coverage=coverage_path,
        require_full_domain=True,
        calculation_sha="c" * 40,
        source_postprocess_sha=source_postprocess,
        audit_workflow_sha=workflow,
        expected_source_run_id="raw-run",
        audit_run_id="31941614867",
    )
    assert manifest["status"] == "full_domain_verified"
    assert manifest["method"] == "independent_oracle"
    result = MODULE.validate_archive(
        archive_dir,
        expected_coverage=coverage_path,
        require_full_domain=True,
        calculation_sha="c" * 40,
        source_postprocess_sha=source_postprocess,
        audit_workflow_sha=workflow,
    )
    assert result["curve_count"] == 3
    source_path = source / "xi_m0p5" / "production_eval" / "prod_eval_xi_m0p5_T_1.csv"
    restored = MODULE.restore_curve_bytes(archive_dir, xi=-0.5, temperature_MeV=1.0)
    assert restored == source_path.read_bytes()
    assert hashlib.sha256(restored).hexdigest() == next(
        row["raw_curve_sha256"]
        for row in MODULE._read_index(archive_dir)
        if row["xi"] == "-0.5" and row["T_MeV"] == "1.0"
    )
    representative = (archive_dir / "representative_index.csv").read_text(encoding="utf-8")
    assert "full_archive_reference" in representative
    assert not (archive_dir / "representative_curves").exists()


def test_strict_archive_rejects_missing_coordinate(tmp_path):
    source, coverage_path, source_postprocess, workflow = _build_sources(tmp_path)
    source_manifest = source / "xi_0" / MODULE.DENSE_SOURCE_MANIFEST_NAME
    source_manifest.unlink()
    (source / "xi_0" / MODULE.DENSE_SOURCE_MANIFEST_HASH_NAME).unlink()
    with pytest.raises(ValueError, match="not full-domain verified|no raw_curve_source_manifest"):
        MODULE.build_archive(
            source,
            tmp_path / "archive",
            expected_coverage=coverage_path,
            require_full_domain=True,
            calculation_sha="c" * 40,
            source_postprocess_sha=source_postprocess,
            audit_workflow_sha=workflow,
        )


def test_archive_validation_rejects_curve_tampering(tmp_path):
    source, coverage_path, source_postprocess, workflow = _build_sources(tmp_path)
    archive_dir = tmp_path / "archive"
    MODULE.build_archive(
        source,
        archive_dir,
        expected_coverage=coverage_path,
        require_full_domain=True,
        calculation_sha="c" * 40,
        source_postprocess_sha=source_postprocess,
        audit_workflow_sha=workflow,
    )
    curve = archive_dir / "curves" / "xi_m0p5" / "T_1.csv"
    curve.write_bytes(curve.read_bytes() + b"#tampered\n")
    with pytest.raises(ValueError, match="hash mismatch"):
        MODULE.validate_archive(archive_dir, expected_coverage=coverage_path, require_full_domain=True)


def test_current_c2_artifact_has_no_recoverable_raw_source():
    source_root = Path(r"D:\Temp\issue130_raw_only_c2_source_31862752226")
    if not source_root.is_dir():
        pytest.skip("local C2 audit artifact is not present")
    audit = MODULE.assess_source_recovery(source_root)
    assert audit["status"] == "not_recoverable"
    assert audit["source_manifest_count"] == 0
    assert "production_eval" in audit["reason"]
