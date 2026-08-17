"""Contract and recovery tools for the PNJL raw rho-mu archive.

This module is deliberately solver-free.  The numerical producer writes the
source CSV files; this module checks their schema and provenance, copies the
source bytes without CSV reserialization, builds a full-domain archive, and
verifies recovery by hash.

The contract is Oracle-only.  The old nine-point targeted bundle is not a
source for this archive and is intentionally not accepted here.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import shutil
from collections import Counter
from pathlib import Path
from typing import Any, Iterable


ARCHIVE_SCHEMA_VERSION = "raw_curve_archive_v1"
SOURCE_SCHEMA_VERSION = "pnjl_c2_raw_only_oracle_job_v1"
DENSE_SOURCE_SCHEMA_VERSION = SOURCE_SCHEMA_VERSION
DENSE_SOURCE_MANIFEST_NAME = "raw_curve_source_manifest.json"
DENSE_SOURCE_MANIFEST_HASH_NAME = "raw_curve_source_manifest.sha256"

SOURCE_REPOSITORY = "w5851/Julia_RelaxTime"
AUDIT_RUN_ID = "31941614867"
SOURCE_GRID_RUN_ID = "31862752226"
CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
SOURCE_POSTPROCESS_SHA = "67d73f871578e35759c08b3c75200c51646cf6cd"
AUDIT_WORKFLOW_SHA = SOURCE_POSTPROCESS_SHA

METHOD = "independent_oracle"
METHODS = (METHOD,)
DENSE_DEFAULT_METHOD = METHOD

RHO_MIN = 0.0
RHO_MAX = 4.0
RHO_STEP = 0.003125
RHO_COUNT = 1281

RAW_CURVE_COLUMNS = (
    "T_MeV",
    "rho",
    "xi",
    "mu_u_MeV",
    "mu_d_MeV",
    "mu_s_MeV",
    "mu_avg_MeV",
    "mu_B_MeV",
    "mu_Q_MeV",
    "mu_S_MeV",
    "pressure_fm4",
    "entropy_fm3",
    "energy_fm4",
    "rho_u_fm3",
    "rho_d_fm3",
    "rho_s_fm3",
    "phi_u",
    "phi_d",
    "phi_s",
    "Phi1",
    "Phi2",
    "M_u_MeV",
    "M_d_MeV",
    "M_s_MeV",
    "iterations",
    "residual_norm",
    "converged",
    "message",
)

# Compatibility names for callers that used the earlier dense-only draft.
DENSE_CURVE_COLUMNS = RAW_CURVE_COLUMNS

RAW_NUMERIC_COLUMNS = (
    "T_MeV",
    "rho",
    "xi",
    "mu_u_MeV",
    "mu_d_MeV",
    "mu_s_MeV",
    "mu_avg_MeV",
    "mu_B_MeV",
    "mu_Q_MeV",
    "mu_S_MeV",
    "pressure_fm4",
    "entropy_fm3",
    "energy_fm4",
    "rho_u_fm3",
    "rho_d_fm3",
    "rho_s_fm3",
    "phi_u",
    "phi_d",
    "phi_s",
    "Phi1",
    "Phi2",
    "M_u_MeV",
    "M_d_MeV",
    "M_s_MeV",
    "residual_norm",
)

RAW_CURVE_SCHEMA = {
    "schema_version": ARCHIVE_SCHEMA_VERSION,
    "source_schema_version": SOURCE_SCHEMA_VERSION,
    "file_format": "UTF-8 CSV; archive curve bytes are copied without reserialization",
    "method": {
        "allowed_values": list(METHODS),
        "value": METHOD,
        "semantics": "independent_oracle fixed-rho equilibrium solve; no production phase pipeline",
    },
    "coordinate_key": ["xi", "T_MeV", "method", "rho"],
    "curve_key": ["xi", "T_MeV", "method"],
    "rho_grid": {
        "min": RHO_MIN,
        "max": RHO_MAX,
        "step": RHO_STEP,
        "count": RHO_COUNT,
        "unit": "rho/rho0 (dimensionless), rho0=0.16 fm^-3",
        "ordering": "source ordering may be ascending or descending; the complete set is required",
    },
    "coordinate_domain": {
        "xi": "dimensionless RS anisotropy parameter",
        "T_MeV": "MeV",
        "rho": "rho/rho0 (dimensionless)",
    },
    "retention": {
        "full": "retain every exact resolved (xi,T_MeV) source CSV in curves/",
        "representative": "representative_index.csv contains references into curve_index.csv only",
        "no_independent_representative_copy": True,
    },
    "columns": [
        {"name": "T_MeV", "type": "float64", "unit": "MeV", "description": "fixed temperature coordinate"},
        {"name": "rho", "type": "float64", "unit": "rho/rho0 (dimensionless)", "description": "normalized density coordinate"},
        {"name": "xi", "type": "float64", "unit": "dimensionless", "description": "RS anisotropy parameter"},
        {"name": "mu_u_MeV", "type": "float64", "unit": "MeV", "description": "up-quark chemical potential"},
        {"name": "mu_d_MeV", "type": "float64", "unit": "MeV", "description": "down-quark chemical potential"},
        {"name": "mu_s_MeV", "type": "float64", "unit": "MeV", "description": "strange-quark chemical potential"},
        {"name": "mu_avg_MeV", "type": "float64", "unit": "MeV", "description": "arithmetic mean of the three quark chemical potentials"},
        {"name": "mu_B_MeV", "type": "float64", "unit": "MeV", "description": "source baryon chemical potential"},
        {"name": "mu_Q_MeV", "type": "float64", "unit": "MeV", "description": "source electric-charge chemical potential"},
        {"name": "mu_S_MeV", "type": "float64", "unit": "MeV", "description": "source strangeness chemical potential"},
        {"name": "pressure_fm4", "type": "float64", "unit": "fm^-4", "description": "natural-unit pressure"},
        {"name": "entropy_fm3", "type": "float64", "unit": "fm^-3", "description": "natural-unit entropy density"},
        {"name": "energy_fm4", "type": "float64", "unit": "fm^-4", "description": "natural-unit energy density"},
        {"name": "rho_u_fm3", "type": "float64", "unit": "fm^-3", "description": "up-quark number density"},
        {"name": "rho_d_fm3", "type": "float64", "unit": "fm^-3", "description": "down-quark number density"},
        {"name": "rho_s_fm3", "type": "float64", "unit": "fm^-3", "description": "strange-quark number density"},
        {"name": "phi_u", "type": "float64", "unit": "model output", "description": "up mean-field output as emitted by the source"},
        {"name": "phi_d", "type": "float64", "unit": "model output", "description": "down mean-field output as emitted by the source"},
        {"name": "phi_s", "type": "float64", "unit": "model output", "description": "strange mean-field output as emitted by the source"},
        {"name": "Phi1", "type": "float64", "unit": "dimensionless", "description": "first Polyakov-loop output"},
        {"name": "Phi2", "type": "float64", "unit": "dimensionless", "description": "second Polyakov-loop output"},
        {"name": "M_u_MeV", "type": "float64", "unit": "MeV", "description": "up constituent mass"},
        {"name": "M_d_MeV", "type": "float64", "unit": "MeV", "description": "down constituent mass"},
        {"name": "M_s_MeV", "type": "float64", "unit": "MeV", "description": "strange constituent mass"},
        {"name": "iterations", "type": "int", "unit": "count", "description": "source nonlinear iteration count"},
        {"name": "residual_norm", "type": "float64", "unit": "solver-defined", "description": "source solver residual norm"},
        {"name": "converged", "type": "bool", "unit": "none", "description": "source solver convergence flag"},
        {"name": "message", "type": "string", "unit": "none", "description": "source solver message"},
    ],
}

INDEX_COLUMNS = (
    "curve_key",
    "xi",
    "T_MeV",
    "method",
    "raw_curve_file",
    "raw_curve_sha256",
    "raw_curve_rows",
    "rho_min",
    "rho_max",
    "source_artifact_name",
    "source_artifact_run_id",
    "source_artifact_url",
    "source_file",
    "source_file_sha256",
    "source_manifest_file",
    "source_manifest_sha256",
    "calculation_sha",
    "source_postprocess_sha",
    "source_workflow_sha",
    "audit_workflow_sha",
    "representative",
    "representative_reason",
)

REPRESENTATIVE_COLUMNS = (
    "curve_key",
    "xi",
    "T_MeV",
    "method",
    "raw_curve_file",
    "representative_reason",
    "full_archive_reference",
)

SOURCE_ARTIFACT_COLUMNS = (
    "source_manifest_file",
    "source_manifest_sha256",
    "source_artifact_name",
    "source_run_id",
    "source_artifact_url",
    "xi",
    "curve_count",
    "source_grid_run_id",
    "source_grid_artifact_name",
    "source_grid_manifest_path",
    "source_grid_manifest_sha256",
)


def _source_run_url(run_id: str) -> str:
    return f"https://github.com/{SOURCE_REPOSITORY}/actions/runs/{run_id}"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != RAW_CURVE_COLUMNS:
            raise ValueError(f"raw curve schema mismatch: {path}")
        rows = list(reader)
    return rows


def _write_csv(path: Path, rows: Iterable[dict[str, Any]], fieldnames: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _float(value: Any, field: str, path: Path) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"invalid {field} in {path}: {value!r}") from error
    if not math.isfinite(result):
        raise ValueError(f"non-finite {field} in {path}: {value!r}")
    return result


def _integer(value: Any, field: str, path: Path) -> int:
    result = _float(value, field, path)
    if not result.is_integer():
        raise ValueError(f"invalid integer {field} in {path}: {value!r}")
    return int(result)


def _true(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def _sha(value: Any, field: str) -> str:
    text = str(value or "").strip().lower()
    if not re.fullmatch(r"[0-9a-f]{40}", text):
        raise ValueError(f"{field} must be a 40-character git SHA: {value!r}")
    return text


def _digest(value: Any, field: str) -> str:
    text = str(value or "").strip().lower()
    if not re.fullmatch(r"[0-9a-f]{64}", text):
        raise ValueError(f"{field} must be a 64-character SHA-256 digest: {value!r}")
    return text


def _canonical_float(value: float) -> str:
    if value == 0.0:
        value = 0.0
    return format(value, ".12g")


def curve_key(xi: float, temperature_MeV: float, method: str = METHOD) -> str:
    return f"{_canonical_float(float(xi))}|{_canonical_float(float(temperature_MeV))}|{method}"


def _token(value: float) -> str:
    text = format(float(value), ".8f").rstrip("0").rstrip(".")
    if text in {"", "-0"}:
        text = "0"
    return text.replace("-", "m").replace(".", "p")


def _rho_grid() -> list[float]:
    return [round(RHO_MIN + index * RHO_STEP, 12) for index in range(RHO_COUNT)]


def _coordinate_key(xi: Any, temperature: Any, method: str = METHOD) -> tuple[float, float, str]:
    return round(float(xi), 10), round(float(temperature), 8), method


def _source_relative(root: Path, path: Path) -> str:
    return path.resolve().relative_to(root.resolve()).as_posix()


def _safe_relative(root: Path, relative: str) -> Path:
    root = root.resolve()
    candidate = (root / relative).resolve()
    if candidate != root and root not in candidate.parents:
        raise ValueError(f"path escapes root: {relative}")
    return candidate


def _hash_sidecar(path: Path) -> str:
    if not path.is_file():
        raise ValueError(f"missing hash sidecar: {path}")
    tokens = path.read_text(encoding="utf-8").strip().split()
    if not tokens:
        raise ValueError(f"empty hash sidecar: {path}")
    return _digest(tokens[0], "manifest hash")


def _write_hash_sidecar(path: Path, digest: str, target_name: str) -> None:
    path.write_text(f"{digest}  {target_name}\n", encoding="utf-8")


def _validate_curve(
    path: Path,
    *,
    xi_expected: float,
    temperature_expected: float,
    method_expected: str = METHOD,
) -> dict[str, Any]:
    rows = _read_csv(path)
    if len(rows) != RHO_COUNT:
        raise ValueError(f"raw curve must contain {RHO_COUNT} rows: {path} ({len(rows)})")
    expected_rhos = _rho_grid()
    observed_rhos: list[float] = []
    seen: set[float] = set()
    for row in rows:
        xi = _float(row["xi"], "xi", path)
        temperature = _float(row["T_MeV"], "T_MeV", path)
        rho = _float(row["rho"], "rho", path)
        if not math.isclose(xi, xi_expected, abs_tol=1e-9, rel_tol=0.0):
            raise ValueError(f"xi mismatch in {path}: {xi} != {xi_expected}")
        if not math.isclose(temperature, temperature_expected, abs_tol=1e-7, rel_tol=0.0):
            raise ValueError(f"temperature mismatch in {path}: {temperature} != {temperature_expected}")
        if not (RHO_MIN - 1e-9 <= rho <= RHO_MAX + 1e-9):
            raise ValueError(f"rho outside [0,4] in {path}: {rho}")
        rounded_rho = round(rho, 10)
        if rounded_rho in seen:
            raise ValueError(f"duplicate rho coordinate in {path}: {rho}")
        seen.add(rounded_rho)
        observed_rhos.append(rho)
        for field in RAW_NUMERIC_COLUMNS:
            _float(row[field], field, path)
        iterations = _integer(row["iterations"], "iterations", path)
        if iterations < 0:
            raise ValueError(f"negative iterations in {path}")
        if row["method"] if "method" in row else False:
            raise ValueError(f"method must be in the archive index, not the 28-column source CSV: {path}")
        if not _true(row["converged"]):
            raise ValueError(f"raw curve contains a non-converged row: {path}")
    expected_set = {round(value, 10) for value in expected_rhos}
    if set(observed_rhos) != expected_set:
        missing = sorted(expected_set.difference(seen))
        extra = sorted(seen.difference(expected_set))
        raise ValueError(f"rho grid mismatch in {path}: missing={missing[:3]}, extra={extra[:3]}")
    return {
        "rows": len(rows),
        "rho_min": min(observed_rhos),
        "rho_max": max(observed_rhos),
        "ordering": "ascending" if observed_rhos[0] < observed_rhos[-1] else "descending",
        "method": method_expected,
    }


def _manifest_files(manifest_path: Path, manifest: dict[str, Any]) -> list[dict[str, Any]]:
    files = manifest.get("files")
    if not isinstance(files, list) or not files:
        raise ValueError(f"source manifest has no files: {manifest_path}")
    result: list[dict[str, Any]] = []
    for entry in files:
        if not isinstance(entry, dict):
            raise ValueError(f"invalid source file entry: {manifest_path}")
        relative = str(entry.get("path", ""))
        if not relative or "production_eval" not in Path(relative).parts:
            raise ValueError(f"source file is outside production_eval: {manifest_path}: {relative}")
        path = _safe_relative(manifest_path.parent, relative)
        if not path.is_file():
            raise ValueError(f"missing source curve: {path}")
        expected_hash = _digest(entry.get("sha256"), "source file sha256")
        actual_hash = sha256(path)
        if actual_hash != expected_hash:
            raise ValueError(f"source file hash mismatch: {path}")
        result.append({**entry, "_path": path, "_sha256": actual_hash})
    return result


def _load_source_manifest(path: Path) -> tuple[dict[str, Any], list[dict[str, Any]], str]:
    manifest = _read_json(path)
    if manifest.get("schema_version") != SOURCE_SCHEMA_VERSION:
        raise ValueError(f"unexpected source schema: {path}")
    if manifest.get("scope") != "c2_full_grid_raw_rho_mu":
        raise ValueError(f"unexpected source scope: {path}")
    if manifest.get("method") != METHOD:
        raise ValueError(f"source method must be {METHOD}: {path}")
    calculation = _sha(manifest.get("calculation_sha"), "calculation_sha")
    if manifest.get("solver_called") is not True:
        raise ValueError(f"source manifest does not record solver_called=true: {path}")
    provenance = manifest.get("provenance", {})
    if provenance.get("phase_pipeline_called") is not False:
        raise ValueError(f"source manifest does not prove phase pipeline was skipped: {path}")
    if provenance.get("higher_order_postprocess_called") is not False:
        raise ValueError(f"source manifest does not prove higher-order postprocess was skipped: {path}")
    files = _manifest_files(path, manifest)
    sidecar = path.with_name(DENSE_SOURCE_MANIFEST_HASH_NAME)
    digest = sha256(path)
    if _hash_sidecar(sidecar) != digest:
        raise ValueError(f"source manifest hash mismatch: {path}")
    return manifest, files, digest


def build_source_manifest(
    process_root: Path,
    output: Path,
    *,
    source_run_id: str,
    source_artifact_name: str,
    shard_id: str | None = None,
    stage: str | None = None,
    xi: float | None = None,
    method: str = METHOD,
    calculation_sha: str = CALCULATION_SHA,
    source_postprocess_sha: str | None = None,
    postprocess_sha: str | None = None,
    source_workflow_sha: str | None = None,
    workflow_sha: str | None = None,
    audit_workflow_sha: str = AUDIT_WORKFLOW_SHA,
    expected_temperatures: Iterable[float] | None = None,
    source_grid_run_id: str = SOURCE_GRID_RUN_ID,
    source_grid_artifact_name: str | None = None,
    source_grid_manifest_path: str | None = None,
    source_grid_manifest_sha256: str | None = None,
    audit_run_id: str = AUDIT_RUN_ID,
) -> dict[str, Any]:
    process_root = process_root.resolve()
    output = output.resolve()
    method == METHOD or (_ for _ in ()).throw(ValueError(f"only method {METHOD} is accepted"))
    calculation_sha = _sha(calculation_sha, "calculation_sha")
    source_postprocess_sha = _sha(source_postprocess_sha or postprocess_sha or SOURCE_POSTPROCESS_SHA, "source_postprocess_sha")
    source_workflow_sha = _sha(source_workflow_sha or workflow_sha or audit_workflow_sha, "source_workflow_sha")
    audit_workflow_sha = _sha(audit_workflow_sha, "audit_workflow_sha")
    source_run_id = str(source_run_id)
    source_artifact_name = str(source_artifact_name)
    if not source_run_id or not source_artifact_name:
        raise ValueError("source_run_id and source_artifact_name are required")
    if source_grid_manifest_path and not source_grid_manifest_sha256:
        raise ValueError("source_grid_manifest_sha256 is required when source_grid_manifest_path is recorded")
    if source_grid_manifest_sha256:
        source_grid_manifest_sha256 = _digest(source_grid_manifest_sha256, "source_grid_manifest_sha256")

    eval_dirs = sorted(path for path in process_root.rglob("production_eval") if path.is_dir())
    if not eval_dirs:
        raise ValueError(f"missing production_eval below {process_root}")
    csv_paths = sorted(path for directory in eval_dirs for path in directory.glob("*.csv"))
    if not csv_paths:
        raise ValueError(f"missing production_eval CSV files below {process_root}")

    entries: list[dict[str, Any]] = []
    seen: set[tuple[float, float, str]] = set()
    expected = None if expected_temperatures is None else {round(float(value), 8) for value in expected_temperatures}
    for path in csv_paths:
        rows = _read_csv(path)
        if not rows:
            raise ValueError(f"empty source curve: {path}")
        xi_value = _float(rows[0]["xi"], "xi", path)
        temperature_value = _float(rows[0]["T_MeV"], "T_MeV", path)
        if xi is not None and not math.isclose(xi_value, float(xi), abs_tol=1e-9, rel_tol=0.0):
            raise ValueError(f"source xi mismatch: {path}")
        summary = _validate_curve(
            path,
            xi_expected=xi_value,
            temperature_expected=temperature_value,
            method_expected=method,
        )
        key = _coordinate_key(xi_value, temperature_value, method)
        if key in seen:
            raise ValueError(f"duplicate source curve coordinate: {key}")
        seen.add(key)
        entries.append({
            "path": _source_relative(process_root, path),
            "sha256": sha256(path),
            "xi": xi_value,
            "T_MeV": temperature_value,
            "rows": summary["rows"],
            "rho_min": summary["rho_min"],
            "rho_max": summary["rho_max"],
            "columns": list(RAW_CURVE_COLUMNS),
        })

    actual_temperatures = {round(float(entry["T_MeV"]), 8) for entry in entries}
    if expected is not None and actual_temperatures != expected:
        raise ValueError(
            f"source temperature coverage mismatch: missing={sorted(expected - actual_temperatures)}, "
            f"extra={sorted(actual_temperatures - expected)}"
        )
    xis = sorted({float(entry["xi"]) for entry in entries})
    if len(xis) != 1:
        raise ValueError(f"one raw-only shard must contain exactly one xi: {xis}")

    payload: dict[str, Any] = {
        "schema_version": SOURCE_SCHEMA_VERSION,
        "scope": "c2_full_grid_raw_rho_mu",
        "method": METHOD,
        "shard_id": str(shard_id or stage or f"xi_{_token(xis[0])}"),
        "xi": xis[0],
        "source_run_id": source_run_id,
        "source_artifact_name": source_artifact_name,
        "source_artifact_url": _source_run_url(source_run_id),
        "calculation_sha": calculation_sha,
        "source_postprocess_sha": source_postprocess_sha,
        "source_workflow_sha": source_workflow_sha,
        "audit_workflow_sha": audit_workflow_sha,
        "audit_run_id": str(audit_run_id),
        "source_grid": {
            "run_id": str(source_grid_run_id),
            "artifact_name": source_grid_artifact_name,
            "manifest_path": source_grid_manifest_path,
            "manifest_sha256": source_grid_manifest_sha256,
        },
        "rho_grid": {
            "min": RHO_MIN,
            "max": RHO_MAX,
            "step": RHO_STEP,
            "count": RHO_COUNT,
            "unit": "rho/rho0 (dimensionless)",
        },
        "solver": {
            "method": METHOD,
            "model_kind": "PNJL",
            "constraint_mode": "fixed_rho",
            "reverse_rho": True,
            "seed_policy": "candidates",
            "p_num": 24,
            "t_num": 8,
            "thermo_quadrature_policy": "rs_reduced_adaptive",
            "thermo_quadrature_rtol": 1e-8,
            "thermo_quadrature_atol": 1e-10,
            "thermo_quadrature_maxevals": 10**7,
            "iterations": 80,
        },
        "solver_called": True,
        "files": sorted(entries, key=lambda entry: (float(entry["T_MeV"]), str(entry["path"]))),
        "provenance": {
            "phase_pipeline_called": False,
            "higher_order_postprocess_called": False,
            "maxwell_called": False,
            "reference_write": False,
            "transport_called": False,
            "source_bytes_are_raw_production_eval": True,
        },
    }
    _write_json(output, payload)
    digest = sha256(output)
    _write_hash_sidecar(output.with_name(DENSE_SOURCE_MANIFEST_HASH_NAME), digest, output.name)
    return payload


def _find_source_manifests(source_root: Path) -> list[Path]:
    return sorted(path for path in source_root.rglob(DENSE_SOURCE_MANIFEST_NAME) if path.is_file())


def assess_source_recovery(source_root: Path) -> dict[str, Any]:
    source_root = source_root.resolve()
    manifests = _find_source_manifests(source_root) if source_root.is_dir() else []
    if not manifests:
        production_dirs = list(source_root.rglob("production_eval")) if source_root.is_dir() else []
        return {
            "status": "not_recoverable",
            "reason": "no raw_curve_source_manifest.json was found; current artifacts do not contain raw production_eval source",
            "production_eval_directories": len(production_dirs),
            "source_manifest_count": 0,
            "coordinate_count": 0,
        }
    records: list[dict[str, Any]] = []
    issues: list[str] = []
    coordinates: set[tuple[float, float, str]] = set()
    for path in manifests:
        try:
            manifest, files, digest = _load_source_manifest(path)
            for entry in files:
                coordinates.add(_coordinate_key(entry["xi"], entry["T_MeV"], manifest["method"]))
            records.append({
                "path": path.relative_to(source_root).as_posix(),
                "sha256": digest,
                "curve_count": len(files),
                "xi": manifest.get("xi"),
            })
        except (OSError, ValueError, json.JSONDecodeError) as error:
            issues.append(f"{path}: {error}")
    return {
        "status": "raw_source_available" if not issues else "raw_source_invalid",
        "source_manifest_count": len(manifests),
        "coordinate_count": len(coordinates),
        "curve_count": sum(record["curve_count"] for record in records),
        "manifests": records,
        "issues": issues,
    }


assess_dense_recovery = assess_source_recovery


def build_expected_coverage(
    grid_csv: Path,
    output: Path,
    *,
    source_grid_run_id: str,
    source_grid_artifact_name: str,
    source_grid_manifest_path: str,
    source_grid_manifest_sha256: str,
    calculation_sha: str = CALCULATION_SHA,
    source_postprocess_sha: str = SOURCE_POSTPROCESS_SHA,
    audit_run_id: str = AUDIT_RUN_ID,
    method: str = METHOD,
    matrix_output: Path | None = None,
) -> dict[str, Any]:
    calculation_sha = _sha(calculation_sha, "calculation_sha")
    source_postprocess_sha = _sha(source_postprocess_sha, "source_postprocess_sha")
    method == METHOD or (_ for _ in ()).throw(ValueError(f"only method {METHOD} is accepted"))
    coordinates: set[tuple[float, float]] = set()
    with grid_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"axis", "xi", "T_MeV"}
        if not required.issubset(set(reader.fieldnames or ())):
            raise ValueError(f"phase_grid_convergence CSV lacks {sorted(required)}: {grid_csv}")
        for row in reader:
            if row.get("axis") != "rho":
                continue
            if str(row.get("xi", "")).strip() == "" or str(row.get("T_MeV", "")).strip() == "":
                raise ValueError(f"rho coverage row lacks xi/T_MeV: {grid_csv}")
            xi = _float(row["xi"], "xi", grid_csv)
            temperature = _float(row["T_MeV"], "T_MeV", grid_csv)
            coordinates.add((round(xi, 10), round(temperature, 8)))
    if not coordinates:
        raise ValueError(f"no axis=rho coordinates in {grid_csv}")
    ordered = sorted(coordinates)
    by_xi: dict[str, list[float]] = {}
    for xi, temperature in ordered:
        by_xi.setdefault(_canonical_float(xi), []).append(temperature)
    xi_values = sorted({xi for xi, _ in ordered})
    coverage = {
        "schema_version": "pnjl_c2_raw_only_expected_coverage_v1",
        "scope": "c2_full_grid_exact_resolved_rho_slices",
        "method": METHOD,
        "calculation_sha": calculation_sha,
        "source_postprocess_sha": source_postprocess_sha,
        "audit_run_id": str(audit_run_id),
        "source_grid": {
            "run_id": str(source_grid_run_id),
            "artifact_name": str(source_grid_artifact_name),
            "manifest_path": str(source_grid_manifest_path),
            "manifest_sha256": _digest(source_grid_manifest_sha256, "source_grid_manifest_sha256"),
            "grid_csv": grid_csv.name,
            "grid_csv_sha256": sha256(grid_csv),
        },
        "coordinate_policy": "exact set of finite axis=rho (xi,T_MeV) rows; no Cartesian completion",
        "coordinates": [{"xi": xi, "T_MeV": temperature, "method": METHOD} for xi, temperature in ordered],
        "xi_values": xi_values,
        "temperatures_by_xi": {key: values for key, values in sorted(by_xi.items())},
        "rho_grid": {
            "min": RHO_MIN,
            "max": RHO_MAX,
            "step": RHO_STEP,
            "count": RHO_COUNT,
            "unit": "rho/rho0 (dimensionless)",
        },
        "expected_curve_count": len(ordered),
    }
    _write_json(output, coverage)
    if matrix_output is not None:
        matrix = [
            {
                "shard_id": f"xi_{_token(xi)}",
                "xi": _canonical_float(xi),
                "temperatures_json": json.dumps(by_xi[_canonical_float(xi)], separators=(",", ":")),
            }
            for xi in xi_values
        ]
        _write_json(matrix_output, {"include": matrix})
    return coverage


def _load_expected_coverage(value: Path | dict[str, Any] | None) -> dict[str, Any] | None:
    if value is None:
        return None
    payload = _read_json(value) if isinstance(value, Path) else value
    if payload.get("schema_version") != "pnjl_c2_raw_only_expected_coverage_v1":
        raise ValueError("unexpected expected coverage schema")
    if payload.get("method") != METHOD:
        raise ValueError("expected coverage method must be independent_oracle")
    coordinates = payload.get("coordinates")
    if not isinstance(coordinates, list) or not coordinates:
        raise ValueError("expected coverage has no coordinates")
    normalized: set[tuple[float, float, str]] = set()
    for entry in coordinates:
        if not isinstance(entry, dict) or entry.get("method", METHOD) != METHOD:
            raise ValueError("invalid expected coverage coordinate")
        normalized.add(_coordinate_key(entry["xi"], entry["T_MeV"], METHOD))
    payload = dict(payload)
    payload["_coordinate_keys"] = normalized
    return payload


def _representative_selection(keys: Iterable[tuple[float, float, str]]) -> dict[tuple[float, float, str], str]:
    grouped: dict[float, list[tuple[float, float, str]]] = {}
    for key in keys:
        grouped.setdefault(key[0], []).append(key)
    selected: dict[tuple[float, float, str], str] = {}
    for xi, rows in grouped.items():
        ordered = sorted(rows, key=lambda key: key[1])
        selected[ordered[0]] = "lowest available T for each xi"
        selected[ordered[-1]] = "highest available T for each xi"
        middle = min(ordered, key=lambda key: (abs(key[1] - 147.0), key[1]))
        selected[middle] = "nearest available T to C2 CEP-region anchor 147 MeV for each xi"
    for xi_target, label in (
        (min(grouped), "global xi lower endpoint"),
        (0.0, "central xi anchor") if 0.0 in grouped else (min(grouped), "global xi lower endpoint"),
        (max(grouped), "global xi upper endpoint"),
    ):
        rows = sorted(grouped[xi_target], key=lambda key: key[1])
        for target_T in (1.0, 147.0, 240.0):
            selected[min(rows, key=lambda key: (abs(key[1] - target_T), key[1]))] = label + f"; nearest T={target_T:g} MeV"
    return selected


def _source_artifact_row(manifest_path: Path, manifest: dict[str, Any], digest: str, source_root: Path) -> dict[str, Any]:
    grid = manifest.get("source_grid", {})
    return {
        "source_manifest_file": _source_relative(source_root, manifest_path),
        "source_manifest_sha256": digest,
        "source_artifact_name": manifest.get("source_artifact_name", ""),
        "source_run_id": manifest.get("source_run_id", ""),
        "source_artifact_url": manifest.get("source_artifact_url", ""),
        "xi": manifest.get("xi", ""),
        "curve_count": len(manifest.get("files", [])),
        "source_grid_run_id": grid.get("run_id", ""),
        "source_grid_artifact_name": grid.get("artifact_name", ""),
        "source_grid_manifest_path": grid.get("manifest_path", ""),
        "source_grid_manifest_sha256": grid.get("manifest_sha256", ""),
    }


def build_archive(
    source_root: Path,
    output_dir: Path,
    *,
    expected_coverage: Path | dict[str, Any] | None = None,
    require_full_domain: bool = False,
    calculation_sha: str = CALCULATION_SHA,
    source_postprocess_sha: str | None = None,
    expected_source_postprocess_sha: str | None = None,
    audit_workflow_sha: str = AUDIT_WORKFLOW_SHA,
    expected_source_run_id: str | None = None,
    audit_run_id: str = AUDIT_RUN_ID,
) -> dict[str, Any]:
    source_root = source_root.resolve()
    output_dir = output_dir.resolve()
    calculation_sha = _sha(calculation_sha, "calculation_sha")
    source_postprocess_sha = _sha(
        source_postprocess_sha or expected_source_postprocess_sha or SOURCE_POSTPROCESS_SHA,
        "source_postprocess_sha",
    )
    audit_workflow_sha = _sha(audit_workflow_sha, "audit_workflow_sha")
    expected = _load_expected_coverage(expected_coverage)
    manifests = _find_source_manifests(source_root)
    if not manifests:
        raise ValueError(f"no {DENSE_SOURCE_MANIFEST_NAME} below {source_root}")

    records: list[tuple[Path, dict[str, Any], list[dict[str, Any]], str]] = []
    actual: dict[tuple[float, float, str], tuple[Path, dict[str, Any], dict[str, Any], str]] = {}
    source_artifacts: dict[str, dict[str, Any]] = {}
    for manifest_path in manifests:
        manifest, files, digest = _load_source_manifest(manifest_path)
        if _sha(manifest["calculation_sha"], "calculation_sha") != calculation_sha:
            raise ValueError(f"source calculation SHA mismatch: {manifest_path}")
        manifest_source_postprocess = _sha(
            manifest.get("source_postprocess_sha"), "source_postprocess_sha"
        )
        if manifest_source_postprocess != source_postprocess_sha:
            raise ValueError(f"source postprocess SHA mismatch: {manifest_path}")
        if expected_source_run_id is not None and str(manifest.get("source_run_id")) != str(expected_source_run_id):
            raise ValueError(f"source run ID mismatch: {manifest_path}")
        records.append((manifest_path, manifest, files, digest))
        source_artifacts[digest] = _source_artifact_row(manifest_path, manifest, digest, source_root)
        for entry in files:
            key = _coordinate_key(entry["xi"], entry["T_MeV"], manifest["method"])
            if key in actual:
                raise ValueError(f"duplicate archive coordinate: {key}")
            actual[key] = (manifest_path, manifest, entry, digest)

    actual_keys = set(actual)
    expected_keys = set(expected["_coordinate_keys"]) if expected is not None else set()
    missing = sorted(expected_keys - actual_keys)
    extra = sorted(actual_keys - expected_keys) if expected is not None else []
    status = "full_domain_verified" if expected is not None and not missing and not extra else "partial_only"
    if require_full_domain and status != "full_domain_verified":
        raise ValueError(f"raw archive is not full-domain verified: missing={len(missing)}, extra={len(extra)}")

    if output_dir.exists():
        shutil.rmtree(output_dir)
    (output_dir / "curves").mkdir(parents=True)
    (output_dir / "source_manifests").mkdir(parents=True)

    index_rows: list[dict[str, Any]] = []
    for key in sorted(actual_keys):
        manifest_path, manifest, entry, digest = actual[key]
        xi, temperature, method = key
        raw_name = f"curves/xi_{_token(xi)}/T_{_token(temperature)}.csv"
        destination = output_dir / raw_name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(entry["_path"], destination)
        if destination.read_bytes() != entry["_path"].read_bytes():
            raise ValueError(f"source byte copy mismatch: {entry['_path']}")
        source_manifest_name = f"source_manifests/{digest}.json"
        source_manifest_destination = output_dir / source_manifest_name
        if not source_manifest_destination.exists():
            shutil.copyfile(manifest_path, source_manifest_destination)
            source_sidecar = manifest_path.with_name(DENSE_SOURCE_MANIFEST_HASH_NAME)
            shutil.copyfile(source_sidecar, output_dir / f"source_manifests/{digest}.sha256")
        index_rows.append({
            "curve_key": curve_key(xi, temperature, method),
            "xi": xi,
            "T_MeV": temperature,
            "method": method,
            "raw_curve_file": raw_name,
            "raw_curve_sha256": sha256(destination),
            "raw_curve_rows": entry["rows"],
            "rho_min": entry["rho_min"],
            "rho_max": entry["rho_max"],
            "source_artifact_name": manifest.get("source_artifact_name", ""),
            "source_artifact_run_id": manifest.get("source_run_id", ""),
            "source_artifact_url": manifest.get("source_artifact_url", ""),
            "source_file": _source_relative(source_root, entry["_path"]),
            "source_file_sha256": entry["_sha256"],
            "source_manifest_file": source_manifest_name,
            "source_manifest_sha256": digest,
            "calculation_sha": manifest["calculation_sha"],
            "source_postprocess_sha": manifest["source_postprocess_sha"],
            "source_workflow_sha": manifest["source_workflow_sha"],
            "audit_workflow_sha": audit_workflow_sha,
            "representative": "false",
            "representative_reason": "",
        })

    selections = _representative_selection(actual_keys)
    by_key = {row["curve_key"]: row for row in index_rows}
    representative_rows: list[dict[str, Any]] = []
    for key, reason in sorted(selections.items()):
        row = by_key[curve_key(*key)]
        row["representative"] = "true"
        row["representative_reason"] = reason
        representative_rows.append({
            "curve_key": row["curve_key"],
            "xi": row["xi"],
            "T_MeV": row["T_MeV"],
            "method": row["method"],
            "raw_curve_file": row["raw_curve_file"],
            "representative_reason": reason,
            "full_archive_reference": "curve_index.csv",
        })

    _write_csv(output_dir / "curve_index.csv", index_rows, INDEX_COLUMNS)
    _write_csv(output_dir / "representative_index.csv", representative_rows, REPRESENTATIVE_COLUMNS)
    _write_csv(
        output_dir / "source_artifacts.csv",
        [source_artifacts[key] for key in sorted(source_artifacts)],
        SOURCE_ARTIFACT_COLUMNS,
    )
    _write_json(output_dir / "schema.json", RAW_CURVE_SCHEMA)

    file_hashes: dict[str, str] = {}
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and path.name not in {"archive_manifest.json", "archive_manifest.sha256"}:
            file_hashes[_source_relative(output_dir, path)] = sha256(path)

    source_postprocess_values = sorted({_sha(manifest["source_postprocess_sha"], "source_postprocess_sha") for _, manifest, _, _ in records})
    source_workflow_values = sorted({_sha(manifest["source_workflow_sha"], "source_workflow_sha") for _, manifest, _, _ in records})
    source_run_values = sorted({str(manifest.get("source_run_id", "")) for _, manifest, _, _ in records})
    source_grid_values = sorted({
        (
            str(manifest.get("source_grid", {}).get("run_id", "")),
            str(manifest.get("source_grid", {}).get("artifact_name", "")),
            str(manifest.get("source_grid", {}).get("manifest_path", "")),
            str(manifest.get("source_grid", {}).get("manifest_sha256", "")),
        )
        for _, manifest, _, _ in records
    })
    archive_manifest: dict[str, Any] = {
        "schema_version": ARCHIVE_SCHEMA_VERSION,
        "archive_id": "issue130-rho-mu-curve-retention",
        "scope": "c2_full_grid_raw_rho_mu",
        "status": status,
        "method": METHOD,
        "calculation_sha": calculation_sha,
        "source_postprocess_sha": source_postprocess_values[0] if len(source_postprocess_values) == 1 else source_postprocess_values,
        "source_workflow_sha": source_workflow_values[0] if len(source_workflow_values) == 1 else source_workflow_values,
        "audit_workflow_sha": audit_workflow_sha,
        "audit_run_id": str(audit_run_id),
        "source_run_ids": source_run_values,
        "source_grid": [
            {
                "run_id": run_id,
                "artifact_name": artifact_name,
                "manifest_path": manifest_path,
                "manifest_sha256": manifest_sha,
            }
            for run_id, artifact_name, manifest_path, manifest_sha in source_grid_values
        ],
        "schema": RAW_CURVE_SCHEMA,
        "rho_grid": RAW_CURVE_SCHEMA["rho_grid"],
        "coverage": {
            "expected_coordinate_count": len(expected_keys) if expected is not None else None,
            "actual_coordinate_count": len(actual_keys),
            "missing_curve_keys": [curve_key(*key) for key in missing],
            "extra_curve_keys": [curve_key(*key) for key in extra],
            "coordinate_policy": expected.get("coordinate_policy") if expected is not None else "source manifests only",
            "expected_coverage_file_sha256": sha256(expected_coverage) if isinstance(expected_coverage, Path) else None,
        },
        "retention": {
            "full_curve_count": len(index_rows),
            "representative_curve_count": len(representative_rows),
            "full_archive_root": "curves/",
            "full_index": "curve_index.csv",
            "representative_index": "representative_index.csv",
            "representatives_are_references_only": True,
        },
        "provenance": {
            "source_artifact_type": "GitHub Actions raw-only shard artifacts",
            "source_manifest_type": DENSE_SOURCE_MANIFEST_NAME,
            "source_manifest_hash_type": "sha256",
            "source_curve_bytes_copied_without_reserialization": True,
            "phase_pipeline_called": False,
            "higher_order_postprocess_called": False,
            "maxwell_called": False,
            "reference_write": False,
            "transport_called": False,
        },
        "files": file_hashes,
    }
    _write_json(output_dir / "archive_manifest.json", archive_manifest)
    _write_hash_sidecar(
        output_dir / "archive_manifest.sha256",
        sha256(output_dir / "archive_manifest.json"),
        "archive_manifest.json",
    )
    return archive_manifest


build_dense_source_manifest = build_source_manifest
build_dense_archive = build_archive


def _load_archive_manifest(archive_dir: Path) -> dict[str, Any]:
    path = archive_dir / "archive_manifest.json"
    manifest = _read_json(path)
    if manifest.get("schema_version") != ARCHIVE_SCHEMA_VERSION:
        raise ValueError(f"unexpected archive schema: {path}")
    if _hash_sidecar(archive_dir / "archive_manifest.sha256") != sha256(path):
        raise ValueError(f"archive manifest hash mismatch: {path}")
    return manifest


def _read_index(archive_dir: Path) -> list[dict[str, str]]:
    path = archive_dir / "curve_index.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != INDEX_COLUMNS:
            raise ValueError(f"archive curve index schema mismatch: {path}")
        return list(reader)


def validate_archive(
    archive_dir: Path,
    *,
    expected_coverage: Path | dict[str, Any] | None = None,
    require_full_domain: bool = False,
    calculation_sha: str | None = None,
    source_postprocess_sha: str | None = None,
    expected_source_postprocess_sha: str | None = None,
    audit_workflow_sha: str | None = None,
) -> dict[str, Any]:
    archive_dir = archive_dir.resolve()
    manifest = _load_archive_manifest(archive_dir)
    if manifest.get("method") != METHOD:
        raise ValueError("archive method is not independent_oracle")
    if calculation_sha is not None and _sha(manifest["calculation_sha"], "calculation_sha") != _sha(calculation_sha, "calculation_sha"):
        raise ValueError("archive calculation SHA mismatch")
    expected_postprocess = source_postprocess_sha or expected_source_postprocess_sha
    if expected_postprocess is not None:
        actual_postprocess = manifest["source_postprocess_sha"]
        if isinstance(actual_postprocess, list) or _sha(actual_postprocess, "source_postprocess_sha") != _sha(expected_postprocess, "source_postprocess_sha"):
            raise ValueError("archive source postprocess SHA mismatch")
    if audit_workflow_sha is not None and _sha(manifest["audit_workflow_sha"], "audit_workflow_sha") != _sha(audit_workflow_sha, "audit_workflow_sha"):
        raise ValueError("archive audit workflow SHA mismatch")

    expected = _load_expected_coverage(expected_coverage)
    index_rows = _read_index(archive_dir)
    keys: set[tuple[float, float, str]] = set()
    indexed_files: set[str] = set()
    for row in index_rows:
        key = _coordinate_key(row["xi"], row["T_MeV"], row["method"])
        if key in keys:
            raise ValueError(f"duplicate archive index coordinate: {key}")
        keys.add(key)
        if row["method"] != METHOD:
            raise ValueError("archive index contains a non-Oracle method")
        relative = row["raw_curve_file"]
        path = _safe_relative(archive_dir, relative)
        if not path.is_file():
            raise ValueError(f"missing archive curve: {path}")
        indexed_files.add(relative)
        actual_hash = sha256(path)
        if actual_hash != row["raw_curve_sha256"]:
            raise ValueError(f"archive curve hash mismatch: {path}")
        summary = _validate_curve(path, xi_expected=float(row["xi"]), temperature_expected=float(row["T_MeV"]))
        if int(row["raw_curve_rows"]) != summary["rows"]:
            raise ValueError(f"archive curve row count mismatch: {path}")
        source_manifest_path = _safe_relative(archive_dir, row["source_manifest_file"])
        if sha256(source_manifest_path) != row["source_manifest_sha256"]:
            raise ValueError(f"archive source manifest hash mismatch: {source_manifest_path}")
        if row["source_file_sha256"] != row["raw_curve_sha256"]:
            raise ValueError(f"archive source/raw byte hash mismatch: {path}")

    actual_curve_files = {
        _source_relative(archive_dir, path)
        for path in (archive_dir / "curves").rglob("*.csv")
    }
    if actual_curve_files != indexed_files:
        raise ValueError("archive curve files and curve_index.csv differ")

    representative_path = archive_dir / "representative_index.csv"
    with representative_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != REPRESENTATIVE_COLUMNS:
            raise ValueError(f"representative index schema mismatch: {representative_path}")
        representatives = list(reader)
    index_by_key = {row["curve_key"]: row for row in index_rows}
    for row in representatives:
        full = index_by_key.get(row["curve_key"])
        if full is None or full["raw_curve_file"] != row["raw_curve_file"]:
            raise ValueError("representative index does not reference the full archive")
        if full["representative"].lower() != "true":
            raise ValueError("representative index points to an unmarked full curve")

    expected_keys = set(expected["_coordinate_keys"]) if expected is not None else {
        _coordinate_key(row["xi"], row["T_MeV"], row["method"]) for row in index_rows
    }
    missing = sorted(expected_keys - keys)
    extra = sorted(keys - expected_keys) if expected is not None else []
    status = "full_domain_verified" if not missing and not extra else "partial_only"
    if require_full_domain and status != "full_domain_verified":
        raise ValueError(f"archive is not full-domain verified: missing={len(missing)}, extra={len(extra)}")
    if manifest.get("status") != status:
        raise ValueError(f"archive status mismatch: manifest={manifest.get('status')}, actual={status}")

    for relative, digest in manifest.get("files", {}).items():
        path = _safe_relative(archive_dir, relative)
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"archive file hash mismatch: {relative}")
    return {
        "status": status,
        "curve_count": len(index_rows),
        "representative_count": len(representatives),
        "missing_curve_count": len(missing),
        "extra_curve_count": len(extra),
        "archive_manifest_sha256": sha256(archive_dir / "archive_manifest.json"),
    }


validate_dense_archive = validate_archive


def restore_curve_bytes(
    archive_dir: Path,
    *,
    xi: float,
    temperature_MeV: float,
    method: str = METHOD,
) -> bytes:
    if method != METHOD:
        raise ValueError(f"only method {METHOD} is accepted")
    archive_dir = archive_dir.resolve()
    _load_archive_manifest(archive_dir)
    rows = _read_index(archive_dir)
    matches = [
        row for row in rows
        if _coordinate_key(row["xi"], row["T_MeV"], row["method"]) == _coordinate_key(xi, temperature_MeV, method)
    ]
    if len(matches) != 1:
        raise KeyError(f"archive curve not found uniquely: {curve_key(xi, temperature_MeV, method)}")
    row = matches[0]
    path = _safe_relative(archive_dir, row["raw_curve_file"])
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != row["raw_curve_sha256"]:
        raise ValueError(f"archive curve hash mismatch during restore: {path}")
    return data


restore_dense_curve_bytes = restore_curve_bytes


def restore_curve(archive_dir: Path, *, xi: float, temperature_MeV: float, method: str = METHOD) -> list[dict[str, str]]:
    data = restore_curve_bytes(archive_dir, xi=xi, temperature_MeV=temperature_MeV, method=method)
    text = data.decode("utf-8")
    reader = csv.DictReader(text.splitlines())
    return list(reader)


restore_dense_curve = restore_curve


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    assess = subparsers.add_parser("assess-source-recovery", aliases=["assess-dense-recovery"])
    assess.add_argument("--source-root", type=Path, required=True)

    coverage = subparsers.add_parser("build-coverage")
    coverage.add_argument("--grid-csv", type=Path, required=True)
    coverage.add_argument("--output", type=Path, required=True)
    coverage.add_argument("--matrix-output", type=Path)
    coverage.add_argument("--source-grid-run-id", required=True)
    coverage.add_argument("--source-grid-artifact-name", required=True)
    coverage.add_argument("--source-grid-manifest-path", required=True)
    coverage.add_argument("--source-grid-manifest-sha256", required=True)
    coverage.add_argument("--calculation-sha", default=CALCULATION_SHA)
    coverage.add_argument("--source-postprocess-sha", default=SOURCE_POSTPROCESS_SHA)
    coverage.add_argument("--audit-run-id", default=AUDIT_RUN_ID)

    source = subparsers.add_parser("build-source-manifest", aliases=["build-dense-source-manifest"])
    source.add_argument("--process-root", type=Path, required=True)
    source.add_argument("--output", type=Path, required=True)
    source.add_argument("--source-run-id", required=True)
    source.add_argument("--source-artifact-name", required=True)
    source.add_argument("--shard-id")
    source.add_argument("--stage")
    source.add_argument("--xi", type=float)
    source.add_argument("--method", default=METHOD)
    source.add_argument("--calculation-sha", default=CALCULATION_SHA)
    source.add_argument("--source-postprocess-sha", dest="source_postprocess_sha", default=SOURCE_POSTPROCESS_SHA)
    source.add_argument("--postprocess-sha", dest="source_postprocess_sha_legacy")
    source.add_argument("--source-workflow-sha", dest="source_workflow_sha", default=None)
    source.add_argument("--workflow-sha", dest="source_workflow_sha_legacy")
    source.add_argument("--audit-workflow-sha", default=AUDIT_WORKFLOW_SHA)
    source.add_argument("--expected-temperatures-json", type=Path)
    source.add_argument("--source-grid-run-id", default=SOURCE_GRID_RUN_ID)
    source.add_argument("--source-grid-artifact-name")
    source.add_argument("--source-grid-manifest-path")
    source.add_argument("--source-grid-manifest-sha256")
    source.add_argument("--audit-run-id", default=AUDIT_RUN_ID)

    build = subparsers.add_parser("build-archive", aliases=["build-dense-archive"])
    build.add_argument("--source-root", type=Path, required=True)
    build.add_argument("--output-dir", type=Path, required=True)
    build.add_argument("--expected-coverage", type=Path)
    build.add_argument("--require-full-domain", action="store_true")
    build.add_argument("--calculation-sha", default=CALCULATION_SHA)
    build.add_argument("--source-postprocess-sha", default=None)
    build.add_argument("--expected-source-postprocess-sha", default=None)
    build.add_argument("--audit-workflow-sha", default=AUDIT_WORKFLOW_SHA)
    build.add_argument("--expected-source-run-id")
    build.add_argument("--audit-run-id", default=AUDIT_RUN_ID)

    validate = subparsers.add_parser("validate-archive", aliases=["validate-dense-archive"])
    validate.add_argument("--archive-dir", type=Path, required=True)
    validate.add_argument("--expected-coverage", type=Path)
    validate.add_argument("--require-full-domain", action="store_true")
    validate.add_argument("--calculation-sha")
    validate.add_argument("--source-postprocess-sha")
    validate.add_argument("--expected-source-postprocess-sha")
    validate.add_argument("--audit-workflow-sha")

    restore = subparsers.add_parser("restore-curve", aliases=["restore-dense-curve"])
    restore.add_argument("--archive-dir", type=Path, required=True)
    restore.add_argument("--xi", type=float, required=True)
    restore.add_argument("--temperature-MeV", type=float, required=True)
    restore.add_argument("--method", default=METHOD)
    restore.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.command in {"assess-source-recovery", "assess-dense-recovery"}:
        print(json.dumps(assess_source_recovery(args.source_root), indent=2, sort_keys=True))
        return 0
    if args.command == "build-coverage":
        payload = build_expected_coverage(
            args.grid_csv,
            args.output,
            source_grid_run_id=args.source_grid_run_id,
            source_grid_artifact_name=args.source_grid_artifact_name,
            source_grid_manifest_path=args.source_grid_manifest_path,
            source_grid_manifest_sha256=args.source_grid_manifest_sha256,
            calculation_sha=args.calculation_sha,
            source_postprocess_sha=args.source_postprocess_sha,
            audit_run_id=args.audit_run_id,
            matrix_output=args.matrix_output,
        )
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0
    if args.command in {"build-source-manifest", "build-dense-source-manifest"}:
        expected_temperatures = None
        if args.expected_temperatures_json is not None:
            values = json.loads(args.expected_temperatures_json.read_text(encoding="utf-8"))
            expected_temperatures = [float(value) for value in values]
        payload = build_source_manifest(
            args.process_root,
            args.output,
            source_run_id=args.source_run_id,
            source_artifact_name=args.source_artifact_name,
            shard_id=args.shard_id,
            stage=args.stage,
            xi=args.xi,
            method=args.method,
            calculation_sha=args.calculation_sha,
            source_postprocess_sha=args.source_postprocess_sha or args.source_postprocess_sha_legacy,
            source_workflow_sha=args.source_workflow_sha or args.source_workflow_sha_legacy,
            audit_workflow_sha=args.audit_workflow_sha,
            expected_temperatures=expected_temperatures,
            source_grid_run_id=args.source_grid_run_id,
            source_grid_artifact_name=args.source_grid_artifact_name,
            source_grid_manifest_path=args.source_grid_manifest_path,
            source_grid_manifest_sha256=args.source_grid_manifest_sha256,
            audit_run_id=args.audit_run_id,
        )
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0
    if args.command in {"build-archive", "build-dense-archive"}:
        payload = build_archive(
            args.source_root,
            args.output_dir,
            expected_coverage=args.expected_coverage,
            require_full_domain=args.require_full_domain,
            calculation_sha=args.calculation_sha,
            source_postprocess_sha=args.source_postprocess_sha,
            expected_source_postprocess_sha=args.expected_source_postprocess_sha,
            audit_workflow_sha=args.audit_workflow_sha,
            expected_source_run_id=args.expected_source_run_id,
            audit_run_id=args.audit_run_id,
        )
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0
    if args.command in {"validate-archive", "validate-dense-archive"}:
        payload = validate_archive(
            args.archive_dir,
            expected_coverage=args.expected_coverage,
            require_full_domain=args.require_full_domain,
            calculation_sha=args.calculation_sha,
            source_postprocess_sha=args.source_postprocess_sha,
            expected_source_postprocess_sha=args.expected_source_postprocess_sha,
            audit_workflow_sha=args.audit_workflow_sha,
        )
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0
    if args.command in {"restore-curve", "restore-dense-curve"}:
        data = restore_curve_bytes(
            args.archive_dir,
            xi=args.xi,
            temperature_MeV=args.temperature_MeV,
            method=args.method,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_bytes(data)
        print(json.dumps({"output": str(args.output), "sha256": hashlib.sha256(data).hexdigest()}))
        return 0
    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
