"""Validate and aggregate Issue #130 targeted-closure artifacts.

The numerical runner is Julia and executes only from the frozen calculation
checkout.  This module is intentionally solver-free: it checks hashes,
finite/converged raw rho-mu rows, cache-cost reconciliation, and the overlay
between the production hybrid and independent fine-grid oracle.  It never
changes a calculation or reference artifact.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import re
import shutil
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path
from typing import Any, Iterable

SCHEMA_VERSION = "pnjl_c2_targeted_closure_v1"
JOB_SCHEMA = "pnjl_c2_targeted_closure_job_v1"
CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
METHODS = ("production_hybrid", "independent_oracle")
REGRESSION_TARGETS = {
    "xi_m0p35_T51": (-0.35, 51.0),
    "xi_m0p25_T41": (-0.25, 41.0),
    "xi_m0p20_T41": (-0.20, 41.0),
    "xi_m0p15_T41": (-0.15, 41.0),
    "xi_m0p10_T41": (-0.10, 41.0),
    "xi_0p00_T51": (0.0, 51.0),
    "xi_0p30_T21": (0.30, 21.0),
    "xi_0p35_T51": (0.35, 51.0),
    "xi_0p35_T101": (0.35, 101.0),
}
CEP_TARGETS = {
    "xi_0p125_CEP_midpoint": (0.125, 126.25, 126.1875, 126.3125),
    "xi_0p39375_CEP_midpoint": (0.39375, 113.5, 113.4375, 113.5625),
    "xi_0p50_CEP_midpoint": (0.5, 107.0, 106.9375, 107.0625),
}
TARGETED_METHODS = set(METHODS)
RHO_GRID_STEP = 0.003125
RHO_MAX = 4.0
COST_STOP_MINUTES = 180.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def finite(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def _target_ids(scope: str) -> tuple[str, ...]:
    if scope == "regression_curves":
        return tuple(REGRESSION_TARGETS)
    if scope == "cep_brackets":
        return tuple(CEP_TARGETS)
    return tuple()


def artifact_names(scope: str) -> tuple[str, ...]:
    prefix = "c2-targeted-regression" if scope == "regression_curves" else "c2-targeted-cep"
    return tuple(f"{prefix}-{target_id}" for target_id in _target_ids(scope))


def _safe_extract(zip_path: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    root = destination.resolve()
    with zipfile.ZipFile(zip_path) as archive:
        for member in archive.infolist():
            candidate = (destination / member.filename).resolve()
            if root != candidate and root not in candidate.parents:
                raise ValueError(f"unsafe artifact path: {member.filename}")
        archive.extractall(destination)


def _api_json(url: str, token: str) -> Any:
    request = urllib.request.Request(
        url,
        headers={
            "Accept": "application/vnd.github+json",
            "Authorization": f"Bearer {token}",
            "X-GitHub-Api-Version": "2022-11-28",
        },
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.loads(response.read().decode("utf-8"))


def download_artifacts(repo: str, run_id: str, token: str, scope: str, output_dir: Path) -> dict[str, Any]:
    if not re.fullmatch(r"\d+", str(run_id)):
        raise ValueError("source run id must be numeric")
    if not token:
        raise ValueError("GH_TOKEN is required for artifact download")
    base = f"https://api.github.com/repos/{repo}/actions/runs/{run_id}/artifacts"
    payload = _api_json(base + "?per_page=100", token)
    available = {str(item["name"]): item for item in payload.get("artifacts", [])}
    names = artifact_names(scope)
    missing = [name for name in names if name not in available]
    if missing:
        raise ValueError(f"source artifacts missing: {missing}")
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)
    downloaded: list[dict[str, Any]] = []
    for name in names:
        item = available[name]
        artifact_id = int(item["id"])
        url = f"https://api.github.com/repos/{repo}/actions/artifacts/{artifact_id}/zip"
        request = urllib.request.Request(
            url,
            headers={
                "Accept": "application/vnd.github+json",
                "Authorization": f"Bearer {token}",
                "X-GitHub-Api-Version": "2022-11-28",
            },
        )
        with urllib.request.urlopen(request, timeout=120) as response:
            archive_bytes = response.read()
        archive_path = output_dir / f".{name}.zip"
        archive_path.write_bytes(archive_bytes)
        artifact_dir = output_dir / name
        _safe_extract(archive_path, artifact_dir)
        archive_path.unlink()
        downloaded.append({"name": name, "id": artifact_id, "size_in_bytes": item.get("size_in_bytes")})
    return {"source_run_id": str(run_id), "scope": scope, "artifacts": downloaded}


def _manifest_paths(input_dir: Path) -> list[Path]:
    return sorted(
        path
        for path in input_dir.rglob("manifest.json")
        if path.parent.name in TARGETED_METHODS
    )


def _declared_file(manifest: dict[str, Any], manifest_path: Path, name: str) -> Path:
    declared = manifest.get("files", {}).get(name)
    if not declared:
        raise ValueError(f"manifest missing files.{name}: {manifest_path}")
    path = manifest_path.parent / name
    if not path.is_file():
        raise ValueError(f"missing {name}: {manifest_path}")
    if str(declared) != sha256(path):
        raise ValueError(f"hash mismatch for {name}: {manifest_path}")
    return path


def _validate_curve(path: Path, xi: float, temperature: float) -> int:
    rows = read_csv(path)
    if not rows:
        raise ValueError(f"empty raw curve: {path}")
    seen: set[tuple[float, float, float]] = set()
    required = ("xi", "T_MeV", "rho", "muq_MeV", "residual_norm", "converged", "finite")
    for row in rows:
        if not all(field in row for field in required):
            raise ValueError(f"raw curve schema missing fields: {path}")
        key = (float(row["xi"]), float(row["T_MeV"]), float(row["rho"]))
        if key in seen:
            raise ValueError(f"duplicate raw curve key {key}: {path}")
        seen.add(key)
        if not math.isclose(key[0], xi, abs_tol=1e-8, rel_tol=0.0) or not math.isclose(
            key[1], temperature, abs_tol=1e-8, rel_tol=0.0
        ):
            raise ValueError(f"raw curve target mismatch {key}: {path}")
        if not all(finite(row[field]) for field in ("xi", "T_MeV", "rho", "muq_MeV", "residual_norm")):
            raise ValueError(f"non-finite raw curve row: {path}")
        if not _truthy(row["converged"]) or not _truthy(row["finite"]):
            raise ValueError(f"non-converged raw curve row: {path}")
    return len(rows)


def _validate_cost(path: Path) -> dict[str, Any]:
    rows = read_csv(path)
    if len(rows) != 1:
        raise ValueError(f"expected one method cost row: {path}")
    row = rows[0]
    fields = ("unique_solves", "point_requests", "cache_hits", "failed_points", "runner_seconds")
    if not all(finite(row.get(field)) for field in fields):
        raise ValueError(f"non-finite method cost: {path}")
    unique = int(float(row["unique_solves"]))
    requests = int(float(row["point_requests"]))
    hits = int(float(row["cache_hits"]))
    failed = int(float(row["failed_points"]))
    if min(unique, requests, hits, failed) < 0 or requests != unique + hits:
        raise ValueError(f"point request reconciliation failure: {path}")
    return {**row, "unique_solves": unique, "point_requests": requests, "cache_hits": hits, "failed_points": failed}


def _validate_one_manifest(path: Path, expected_sha: str, expected_postprocess_sha: str, source_run_id: str) -> dict[str, Any]:
    manifest = read_json(path)
    if manifest.get("schema_version") != JOB_SCHEMA:
        raise ValueError(f"unexpected job schema: {path}")
    method = str(manifest.get("method", ""))
    if method not in TARGETED_METHODS:
        raise ValueError(f"unexpected method in manifest: {path}")
    if str(manifest.get("scope", "")) not in {"regression_curves", "cep_brackets"}:
        raise ValueError(f"unexpected targeted scope: {path}")
    if str(manifest.get("calculation_sha", "")).lower() != expected_sha.lower():
        raise ValueError(f"calculation SHA mismatch: {path}")
    if expected_postprocess_sha and str(manifest.get("postprocess_sha", "")) != expected_postprocess_sha:
        raise ValueError(f"postprocess SHA mismatch: {path}")
    if str(manifest.get("source_run_id", source_run_id)) != str(source_run_id):
        raise ValueError(f"source run mismatch: {path}")
    if manifest.get("solver_called") is not True:
        raise ValueError(f"numerical artifact must record solver_called=true: {path}")
    provenance = manifest.get("provenance", {})
    if provenance.get("reference_write") is not False:
        raise ValueError(f"reference write contract violated: {path}")
    if provenance.get("oracle_labels_used_for_routing") is not False:
        raise ValueError(f"oracle routing leakage: {path}")
    target = manifest.get("target", {})
    target_id = str(target.get("id", manifest.get("target_id", "")))
    if not target_id:
        raise ValueError(f"missing target id: {path}")
    if target_id in REGRESSION_TARGETS:
        xi, temperature = REGRESSION_TARGETS[target_id]
    elif target_id in CEP_TARGETS:
        xi, temperature = CEP_TARGETS[target_id][:2]
    else:
        raise ValueError(f"unexpected target id {target_id}: {path}")
    curve = _declared_file(manifest, path, "curve_points.csv")
    diagnostics = _declared_file(manifest, path, "slice_diagnostics.csv")
    accuracy = _declared_file(manifest, path, "accuracy.csv")
    costs = _declared_file(manifest, path, "method_costs.csv")
    curve_rows = _validate_curve(curve, xi, temperature)
    diagnostics_rows = read_csv(diagnostics)
    if len(diagnostics_rows) != 1:
        raise ValueError(f"expected one diagnostics row: {diagnostics}")
    cost = _validate_cost(costs)
    accuracy_rows = read_csv(accuracy)
    if len(accuracy_rows) != 1:
        raise ValueError(f"expected one accuracy row: {accuracy}")
    return {
        "manifest_path": path,
        "manifest": manifest,
        "method": method,
        "target_id": target_id,
        "xi": xi,
        "T_MeV": temperature,
        "curve_path": curve,
        "curve_rows": curve_rows,
        "diagnostics": diagnostics_rows[0],
        "accuracy": accuracy_rows[0],
        "cost": cost,
        "files": [
            {"path": str(item.relative_to(path.parents[2])).replace("\\", "/"), "sha256": sha256(item)}
            for item in (path, curve, diagnostics, accuracy, costs)
        ],
    }


def validate_source(
    input_dir: Path,
    expected_sha: str,
    expected_postprocess_sha: str,
    source_run_id: str,
    scope: str | None = None,
) -> dict[str, Any]:
    manifests = _manifest_paths(input_dir)
    expected_targets = set(REGRESSION_TARGETS) | set(CEP_TARGETS)
    if not manifests:
        raise ValueError(f"no targeted method manifests under {input_dir}")
    records = [_validate_one_manifest(path, expected_sha, expected_postprocess_sha, source_run_id) for path in manifests]
    seen = {(record["target_id"], record["method"]) for record in records}
    # The validator is used per scope directory.  A job directory has exactly
    # two methods; infer the target set from the first manifest and require the
    # complete two-method pair for every target in that directory.
    target_ids = {record["target_id"] for record in records}
    if scope is not None:
        expected_scope_targets = set(_target_ids(scope))
        if target_ids != expected_scope_targets:
            raise ValueError(
                f"incomplete {scope} target matrix; expected={sorted(expected_scope_targets)} "
                f"actual={sorted(target_ids)}"
            )
    if not target_ids <= expected_targets:
        raise ValueError(f"unexpected targeted ids: {sorted(target_ids - expected_targets)}")
    expected_pairs = {(target_id, method) for target_id in target_ids for method in METHODS}
    if seen != expected_pairs:
        raise ValueError(f"incomplete method matrix; missing={sorted(expected_pairs - seen)} extra={sorted(seen - expected_pairs)}")
    input_files = [item for record in records for item in record["files"]]
    total_seconds = sum(float(record["cost"]["runner_seconds"]) for record in records)
    total_unique = sum(int(record["cost"]["unique_solves"]) for record in records)
    failed_points = sum(int(record["cost"]["failed_points"]) for record in records)
    return {
        "schema_version": SCHEMA_VERSION,
        "source_run_id": str(source_run_id),
        "source_calculation_sha": expected_sha,
        "source_postprocess_sha": expected_postprocess_sha,
        "target_count": len(target_ids),
        "method_count": len(METHODS),
        "manifest_count": len(records),
        "solver_called": True,
        "failed_points": failed_points,
        "runner_minutes": total_seconds / 60.0,
        "unique_solves": total_unique,
        "input_files": input_files,
    }


def _as_float(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _plot_curves(output_dir: Path, records: Iterable[dict[str, Any]]) -> list[dict[str, str]]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []
    by_target: dict[str, list[dict[str, Any]]] = {}
    for record in records:
        by_target.setdefault(record["target_id"], []).append(record)
    plots: list[dict[str, str]] = []
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    for target_id, target_records in sorted(by_target.items()):
        figure, axis = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
        for record in sorted(target_records, key=lambda value: value["method"]):
            rows = read_csv(record["curve_path"])
            axis.plot(
                [float(row["rho"]) for row in rows],
                [float(row["muq_MeV"]) for row in rows],
                linewidth=1.1,
                label=record["method"],
            )
        axis.set_xlabel(r"$\rho$")
        axis.set_ylabel(r"$\mu_q$ [MeV]")
        axis.set_title(f"Issue #130 targeted closure: {target_id}")
        axis.grid(True, alpha=0.25)
        axis.legend(frameon=False)
        path = figure_dir / f"rho_mu_{target_id}.png"
        figure.savefig(path, dpi=180)
        plt.close(figure)
        plots.append({"path": str(path.relative_to(output_dir)).replace("\\", "/"), "sha256": sha256(path)})
    return plots


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0]) if rows else ["status"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _aggregate_records_for_scope(
    input_dir: Path,
    expected_sha: str,
    expected_postprocess_sha: str,
    source_run_id: str,
    scope: str,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    validation = validate_source(input_dir, expected_sha, expected_postprocess_sha, source_run_id, scope)
    records = [_validate_one_manifest(path, expected_sha, expected_postprocess_sha, source_run_id) for path in _manifest_paths(input_dir)]
    return validation, records


def aggregate_validate(input_dir: Path, output_dir: Path, scope: str, expected_sha: str, expected_postprocess_sha: str, source_run_id: str) -> dict[str, Any]:
    validation, records = _aggregate_records_for_scope(
        input_dir, expected_sha, expected_postprocess_sha, source_run_id, scope
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    grouped: dict[str, dict[str, dict[str, Any]]] = {}
    for record in records:
        grouped.setdefault(record["target_id"], {})[record["method"]] = record
    overlay_rows: list[dict[str, Any]] = []
    cost_rows: list[dict[str, Any]] = []
    for target_id, methods in sorted(grouped.items()):
        hybrid = methods["production_hybrid"]
        oracle = methods["independent_oracle"]
        hdiag, odiag = hybrid["diagnostics"], oracle["diagnostics"]
        hstatus = hdiag.get("result_status", "")
        ostatus = odiag.get("result_status", "")
        overlay_rows.append({
            "target_id": target_id,
            "xi": hybrid["xi"],
            "T_MeV": hybrid["T_MeV"],
            "hybrid_status": hstatus,
            "oracle_status": ostatus,
            "classification_match": hstatus == ostatus,
            "hybrid_candidate_count": hdiag.get("maxwell_candidate_count", ""),
            "oracle_candidate_count": odiag.get("maxwell_candidate_count", ""),
            "hybrid_crossing_count": hdiag.get("maxwell_crossing_count", ""),
            "oracle_crossing_count": odiag.get("maxwell_crossing_count", ""),
            "hybrid_geometry_converged": hdiag.get("geometry_converged", ""),
            "oracle_geometry_converged": odiag.get("geometry_converged", ""),
            "hybrid_position_error_MeV": hdiag.get("position_error_MeV", ""),
            "hybrid_density_error": hdiag.get("density_error", ""),
            "hybrid_area_residual": hdiag.get("area_residual", ""),
            "oracle_position_error_MeV": odiag.get("position_error_MeV", ""),
            "oracle_density_error": odiag.get("density_error", ""),
            "oracle_area_residual": odiag.get("area_residual", ""),
            "finite_and_converged": _truthy(hdiag.get("finite_and_converged")) and _truthy(odiag.get("finite_and_converged")),
        })
        for record in (hybrid, oracle):
            cost_rows.append({"target_id": target_id, "method": record["method"], **record["cost"]})
    matches = all(row["classification_match"] for row in overlay_rows)
    finite_ok = all(row["finite_and_converged"] for row in overlay_rows)
    total_runner_minutes = sum(float(row["runner_seconds"]) for row in cost_rows) / 60.0
    cost_ok = total_runner_minutes <= COST_STOP_MINUTES
    if scope == "regression_curves":
        verdict = (
            "targeted_cost_inconclusive"
            if not cost_ok
            else "targeted_closure_candidate"
            if matches and finite_ok
            else "targeted_classification_inconclusive"
        )
    else:
        # A single midpoint is evidence for the next CEP refinement step, not
        # a closed bracket certificate.
        verdict = (
            "cep_cost_inconclusive"
            if not cost_ok
            else "cep_midpoint_diagnostic"
            if matches and finite_ok
            else "cep_candidate_inconclusive"
        )
    _write_csv(output_dir / "classification_overlay.csv", overlay_rows)
    _write_csv(output_dir / "cost_frontier.csv", cost_rows)
    plots = _plot_curves(output_dir, records)
    claim_rows = [
        {"claim_id": "raw_curve_finite", "status": str(finite_ok).lower(), "evidence": "classification_overlay.csv"},
        {"claim_id": "hybrid_oracle_classification_match", "status": str(matches).lower(), "evidence": "classification_overlay.csv"},
        {"claim_id": "production_promotion", "status": "not_claimed", "evidence": "manifest.json"},
    ]
    _write_csv(output_dir / "claim_ledger.csv", claim_rows)
    output_files = [output_dir / name for name in ("classification_overlay.csv", "cost_frontier.csv", "claim_ledger.csv")]
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "scope": scope,
        "verdict": verdict,
        "source_run_id": str(source_run_id),
        "source_calculation_sha": expected_sha,
        "source_postprocess_sha": expected_postprocess_sha,
        "solver_called": False,
        "oracle_labels_used_for_routing": False,
        "source_validation": validation,
        "classification_match": matches,
        "finite_and_converged": finite_ok,
        "runner_minutes": total_runner_minutes,
        "cost_stop_minutes": COST_STOP_MINUTES,
        "cost_gate": cost_ok,
        "plots": plots,
        "files": {path.name: sha256(path) for path in output_files},
        "provenance": {
            "reference_write": False,
            "overlay_only": True,
            "route_before_oracle_gate": True,
        },
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def cross_axis_audit(input_dir: Path, output_dir: Path, expected_sha: str) -> dict[str, Any]:
    """Summarize available cross-axis records without solver access.

    The C2 aggregate did not guarantee retention of raw geometry rows.  When
    the caller supplies a directory without a recognized geometry table this
    function emits an explicit input-missing verdict instead of inferring a
    zero-failure result.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    recognized = False
    for path in sorted(input_dir.rglob("*.csv")):
        for row in read_csv(path):
            reason = str(row.get("reason", row.get("failure_reason", row.get("status", ""))))
            axis = str(row.get("axis", row.get("grid_axis", "")))
            if not axis and "rho" in path.name.lower():
                axis = "rho"
            if axis in {"rho", "T", "xi", "temperature", "density"} or "geometry" in path.name.lower():
                recognized = True
                rows.append({"source": str(path.relative_to(input_dir)).replace("\\", "/"), "axis": axis, "reason": reason, **row})
    if recognized:
        verdict = "cross_axis_audit_complete"
        _write_csv(output_dir / "cross_axis_records.csv", rows)
    else:
        verdict = "cross_axis_input_missing"
        _write_csv(output_dir / "cross_axis_records.csv", [{"status": "input_missing"}])
    claim_rows = [{"claim_id": "cross_axis_geometry_audit", "status": verdict, "evidence": "cross_axis_records.csv"}]
    _write_csv(output_dir / "claim_ledger.csv", claim_rows)
    files = {name: sha256(output_dir / name) for name in ("cross_axis_records.csv", "claim_ledger.csv")}
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "scope": "cross_axis_audit",
        "verdict": verdict,
        "source_calculation_sha": expected_sha,
        "solver_called": False,
        "reference_write": False,
        "files": files,
        "row_count": len(rows),
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scope", required=True, choices=["regression_curves", "cep_brackets", "cross_axis_audit"])
    parser.add_argument("--mode", required=True, choices=["download", "source_validate", "aggregate_validate", "cross_axis_audit"])
    parser.add_argument("--input-dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--expected-calculation-sha", default=CALCULATION_SHA)
    parser.add_argument("--expected-source-postprocess-sha", default="")
    parser.add_argument("--source-run-id", default="")
    parser.add_argument("--repo", default="w5851/Julia_RelaxTime")
    parser.add_argument("--run-id", default="")
    args = parser.parse_args()
    if args.expected_calculation_sha.lower() != CALCULATION_SHA:
        raise SystemExit("unexpected calculation SHA")
    if args.mode == "download":
        if args.output_dir is None:
            raise SystemExit("--output-dir is required")
        result = download_artifacts(args.repo, args.run_id, os.environ.get("GH_TOKEN", ""), args.scope, args.output_dir)
        print(json.dumps(result, sort_keys=True))
        return 0
    if args.input_dir is None:
        raise SystemExit("--input-dir is required")
    if args.mode == "source_validate":
        result = validate_source(
            args.input_dir,
            args.expected_calculation_sha,
            args.expected_source_postprocess_sha,
            args.source_run_id,
            args.scope,
        )
        output = args.output_dir or args.input_dir
        output.mkdir(parents=True, exist_ok=True)
        (output / "source_validation.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps(result, sort_keys=True))
        return 0
    if args.mode == "cross_axis_audit":
        if args.output_dir is None:
            raise SystemExit("--output-dir is required")
        result = cross_axis_audit(args.input_dir, args.output_dir, args.expected_calculation_sha)
    else:
        if args.output_dir is None:
            raise SystemExit("--output-dir is required")
        result = aggregate_validate(args.input_dir, args.output_dir, args.scope, args.expected_calculation_sha, args.expected_source_postprocess_sha, args.source_run_id)
    print(json.dumps(result, sort_keys=True))
    return 0 if result.get("verdict") not in {
        "targeted_classification_inconclusive",
        "targeted_cost_inconclusive",
        "cep_candidate_inconclusive",
        "cep_cost_inconclusive",
        "cross_axis_input_missing",
    } else 2


if __name__ == "__main__":
    raise SystemExit(main())
