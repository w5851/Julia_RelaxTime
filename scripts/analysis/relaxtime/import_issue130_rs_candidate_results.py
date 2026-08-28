#!/usr/bin/env python3
"""Import the reviewed Issue #130 RS candidate aggregate into versioned results.

The numerical source is an immutable external Actions aggregate.  This script
validates the aggregate and its selected shard sidecars, copies the two scan
and diagnostic CSVs into versioned ``prod_v2`` result cases, derives only the
portable plan/config/empty-failure sidecars, and records all provenance.  It
never calls a solver, changes runtime defaults, overwrites ``prod_v1``, or
claims numerical production promotion.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any, Iterable


CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
WORKFLOW_HEAD_SHA = "22874505877491754eed27519ad8a7b871c82571"
CASE_NAME = "first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2"
LEGACY_SNAPSHOT_VERSION = "legacy_prod_v1_snapshot_v1"
AGGREGATE_NAME = "aggregate_replay_20260826_v4"
AUDIT_NAME = "post_repair_audit_20260826_v1"
AGGREGATE_MANIFEST_SHA256 = (
    "42e88e028832311a95e45816eaf8ab1bb2dbbe8bb094ed6b59681c51d1cad754"
)
AUDIT_MANIFEST_SHA256 = (
    "587f79db64476f777b25ca5ee234a43b5a31e17bea375864d9e13f0226d1d6f0"
)
AUDIT_SUMMARY_SHA256 = (
    "dc7d5605c5b40e0903659d0ef3c092ec098f41ac9ac98379b59333adf31d2e15"
)
MODES = (
    "mode_a_fixed_muB_phase_scaled",
    "mode_b_fixed_T_sparse_muB",
)
SCAN_FILES = {
    MODES[0]: "mode_a_scan.csv",
    MODES[1]: "mode_b_scan.csv",
}
DIAG_FILES = {
    MODES[0]: "mode_a_diag.csv",
    MODES[1]: "mode_b_diag.csv",
}
EXPECTED_SCAN_ROWS = {MODES[0]: 910, MODES[1]: 909}
EXPECTED_DIAG_ROWS = {MODES[0]: 38220, MODES[1]: 38178}
SCAN_KEYS = {
    MODES[0]: ("muB_MeV", "xi", "alpha_T"),
    MODES[1]: ("T_MeV", "muB_MeV", "xi"),
}
DIAG_KEYS = (
    "T_MeV",
    "muB_MeV",
    "xi",
    "species",
    "channel",
    "density_key",
)
PLAN_FIELDS = ("T_MeV", "muB_MeV", "xi", "mode", "plot_panel", "plot_series")
FAILED_FIELDS = (
    "T_MeV",
    "muB_MeV",
    "xi",
    "seed_source",
    "phase_prev",
    "phase_curr_hint",
    "error_type",
    "error_message",
    "timestamp",
)
DIAG_NUMERIC_FIELDS = (
    "multiplicity",
    "density",
    "rate",
    "contribution",
    "total",
    "tau_inv_species",
)
NONNEGATIVE_DIAG_FIELDS = ("multiplicity", "density", "rate", "contribution", "total", "tau_inv_species")
CONVERGENCE_COPY_FILES = (
    "aggregate_summary.json",
    "candidate_legacy_comparison.csv",
    "candidate_legacy_comparison_summary.json",
    "claim_ledger.csv",
    "cost_summary.csv",
    "per_shard_summary.csv",
    "action_run_summary.csv",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, dict):
        raise ValueError(f"JSON object required: {path}")
    return payload


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def read_csv(path: Path) -> tuple[list[str], list[str], list[dict[str, str]]]:
    text = path.read_text(encoding="utf-8-sig")
    comments = [line for line in text.splitlines() if line.startswith("#")]
    data = "\n".join(line for line in text.splitlines() if not line.startswith("#"))
    reader = csv.DictReader(data.splitlines())
    fieldnames = list(reader.fieldnames or [])
    return comments, fieldnames, [dict(row) for row in reader]


def write_csv(path: Path, comments: Iterable[str], fieldnames: Iterable[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for comment in comments:
            handle.write(f"{comment.rstrip()}\n")
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in writer.fieldnames})


def parse_float(value: Any) -> float | None:
    if value is None or str(value).strip() == "":
        return None
    try:
        return float(str(value).strip())
    except ValueError:
        return None


def normalized_key(row: dict[str, str], fields: Iterable[str]) -> tuple[Any, ...]:
    result: list[Any] = []
    for field in fields:
        value = row.get(field, "")
        parsed = parse_float(value)
        result.append(round(parsed, 12) if parsed is not None else value)
    return tuple(result)


def _require_fields(fieldnames: Iterable[str], required: Iterable[str], label: str) -> None:
    missing = sorted(set(required) - set(fieldnames))
    if missing:
        raise ValueError(f"{label}: missing fields {missing}")


def _numeric_nonfinite_fields(row: dict[str, str], *, mode: str) -> list[str]:
    bad: list[str] = []
    for field, raw in row.items():
        parsed = parse_float(raw)
        if parsed is None:
            continue
        if math.isfinite(parsed):
            continue
        if mode == MODES[1] and field == "alpha_T":
            continue
        bad.append(field)
    return bad


def validate_scan_rows(
    mode: str,
    fieldnames: list[str],
    rows: list[dict[str, str]],
    *,
    expected_rows: int | None = None,
) -> dict[str, Any]:
    required = (
        "T_MeV",
        "muq_MeV",
        "muB_MeV",
        "xi",
        "mode",
        "phase_reference_kind",
        "converged",
        "phase_structure",
        "quality_flag",
        "quality_reason",
        "tau_u",
        "tau_d",
        "tau_s",
        "tau_ubar",
        "tau_dbar",
        "tau_sbar",
        "eta",
        "sigma",
        "zeta",
        "eta_over_s",
        "zeta_over_s",
        "sigma_over_T",
    )
    _require_fields(fieldnames, required, f"{mode} scan")
    if expected_rows is not None and len(rows) != expected_rows:
        raise ValueError(f"{mode} scan rows {len(rows)} != expected {expected_rows}")
    keys = [normalized_key(row, SCAN_KEYS[mode]) for row in rows]
    duplicates = sorted({key for key in keys if keys.count(key) > 1}, key=str)
    if duplicates:
        raise ValueError(f"{mode} scan duplicate keys: {duplicates[:3]}")
    invalid_converged = [index for index, row in enumerate(rows) if row.get("converged", "").lower() != "true"]
    if invalid_converged:
        raise ValueError(f"{mode} scan has non-converged rows: {invalid_converged[:3]}")
    nonfinite = {
        str(index): fields
        for index, row in enumerate(rows)
        if (fields := _numeric_nonfinite_fields(row, mode=mode))
    }
    if nonfinite:
        raise ValueError(f"{mode} scan non-finite fields: {next(iter(nonfinite.items()))}")
    quality_flags = Counter(row.get("quality_flag", "") for row in rows)
    quality_reasons = Counter(row.get("quality_reason", "") for row in rows)
    return {
        "rows": len(rows),
        "duplicate_keys": 0,
        "converged_counts": dict(Counter(row.get("converged", "") for row in rows)),
        "quality_flag_counts": dict(quality_flags),
        "quality_reason_counts": dict(quality_reasons),
        "key_fields": list(SCAN_KEYS[mode]),
    }


def validate_diag_rows(
    mode: str,
    fieldnames: list[str],
    rows: list[dict[str, str]],
    *,
    expected_rows: int | None = None,
) -> dict[str, Any]:
    required = DIAG_KEYS + DIAG_NUMERIC_FIELDS + ("species", "channel")
    _require_fields(fieldnames, required, f"{mode} diagnostic")
    if expected_rows is not None and len(rows) != expected_rows:
        raise ValueError(f"{mode} diagnostic rows {len(rows)} != expected {expected_rows}")
    keys = [normalized_key(row, DIAG_KEYS) for row in rows]
    duplicates = sorted({key for key in keys if keys.count(key) > 1}, key=str)
    if duplicates:
        raise ValueError(f"{mode} diagnostic duplicate keys: {duplicates[:3]}")
    negative: list[tuple[int, str]] = []
    nonfinite: list[tuple[int, str]] = []
    for index, row in enumerate(rows):
        for field in DIAG_NUMERIC_FIELDS:
            value = parse_float(row.get(field))
            if value is None or not math.isfinite(value):
                nonfinite.append((index, field))
            elif field in NONNEGATIVE_DIAG_FIELDS and value < -1e-12:
                negative.append((index, field))
    if nonfinite:
        raise ValueError(f"{mode} diagnostic non-finite fields: {nonfinite[:3]}")
    if negative:
        raise ValueError(f"{mode} diagnostic negative fields: {negative[:3]}")
    return {"rows": len(rows), "duplicate_keys": 0, "nonfinite": 0, "negative": 0}


def validate_direct_coexistence(rows: list[dict[str, str]]) -> dict[str, Any]:
    selected = [
        row
        for row in rows
        if abs(float(row["muB_MeV"]) - 900.0) < 1e-9
        and abs(float(row["alpha_T"]) - 1.0) < 1e-9
    ]
    xis = {round(float(row["xi"]), 12) for row in selected}
    if -0.003 not in xis or 0.003 not in xis or 0.0 in xis:
        raise ValueError(f"direct coexistence replacement contract failed: {sorted(xis)}")
    return {
        "rows": len(selected),
        "has_minus_0003": True,
        "has_plus_0003": True,
        "has_zero": False,
        "xi_min": min(xis),
        "xi_max": max(xis),
    }


def validate_manifest_file(path: Path, expected_sha256: str, schema: str) -> dict[str, Any]:
    actual = sha256_file(path)
    if actual != expected_sha256:
        raise ValueError(f"{path} hash {actual} != expected {expected_sha256}")
    payload = read_json(path)
    if payload.get("schema_version") != schema:
        raise ValueError(f"unexpected schema in {path}: {payload.get('schema_version')}")
    return payload


def validate_external_source(source_root: Path) -> dict[str, Any]:
    aggregate_root = source_root / AGGREGATE_NAME
    audit_root = source_root / AUDIT_NAME
    aggregate_manifest_path = aggregate_root / "manifest.json"
    aggregate_summary_path = aggregate_root / "aggregate_summary.json"
    audit_manifest_path = audit_root / "manifest.json"
    audit_summary_path = audit_root / "audit_summary.json"
    aggregate_manifest = validate_manifest_file(
        aggregate_manifest_path,
        AGGREGATE_MANIFEST_SHA256,
        "issue130_rs_sharded_production_v2_aggregate_manifest_v4",
    )
    aggregate_summary_sha = next(
        item["sha256"]
        for item in aggregate_manifest["files"]
        if item["path"] == "aggregate_summary.json"
    )
    aggregate_summary = validate_manifest_file(
        aggregate_summary_path,
        aggregate_summary_sha,
        "issue130_rs_sharded_production_v2_aggregate_replay_v4",
    )
    audit_manifest = validate_manifest_file(
        audit_manifest_path,
        AUDIT_MANIFEST_SHA256,
        "issue130_rs_sharded_production_v2_post_repair_audit_v1_manifest",
    )
    audit_summary = validate_manifest_file(
        audit_summary_path,
        AUDIT_SUMMARY_SHA256,
        "issue130_rs_sharded_production_v2_post_repair_audit_v1",
    )
    if aggregate_manifest.get("calculation_sha") != CALCULATION_SHA:
        raise ValueError("aggregate calculation SHA mismatch")
    if aggregate_manifest.get("workflow_head_sha") != WORKFLOW_HEAD_SHA:
        raise ValueError("aggregate workflow head SHA mismatch")
    if aggregate_summary.get("selected_unique_run_count") != 30:
        raise ValueError("aggregate selected run count is not 30")
    if not aggregate_summary.get("all_hard_gates_including_timing"):
        raise ValueError("aggregate hard gates are not all passing")
    if audit_summary.get("verdict") != "post_repair_audit_pass_diagnostic_only":
        raise ValueError("post-repair audit is not the accepted diagnostic-only verdict")
    if audit_summary.get("calculation_sha") != CALCULATION_SHA:
        raise ValueError("post-repair audit calculation SHA mismatch")
    for item in aggregate_manifest.get("files", []):
        path = aggregate_root / item["path"]
        if not path.is_file() or sha256_file(path) != item["sha256"]:
            raise ValueError(f"aggregate file hash mismatch: {path}")
    for item in audit_manifest.get("files", []):
        path = audit_root / item["path"]
        if not path.is_file() or sha256_file(path) != item["sha256"]:
            raise ValueError(f"post-repair audit file hash mismatch: {path}")
    for mode in MODES:
        scan_path = aggregate_root / SCAN_FILES[mode]
        diag_path = aggregate_root / DIAG_FILES[mode]
        _, scan_fields, scan_rows = read_csv(scan_path)
        _, diag_fields, diag_rows = read_csv(diag_path)
        scan_summary = validate_scan_rows(
            mode, scan_fields, scan_rows, expected_rows=EXPECTED_SCAN_ROWS[mode]
        )
        diag_summary = validate_diag_rows(
            mode, diag_fields, diag_rows, expected_rows=EXPECTED_DIAG_ROWS[mode]
        )
        if mode == MODES[0]:
            direct = validate_direct_coexistence(scan_rows)
        else:
            direct = None
        aggregate_manifest.setdefault("_validated", {})[mode] = {
            "scan": scan_summary,
            "diagnostic": diag_summary,
            "direct_coexistence": direct,
        }
    if aggregate_summary.get("direct_coexistence_gate", {}).get("has_minus_0003") is not True:
        raise ValueError("aggregate direct-coexistence gate missing xi=-0.003")
    return {
        "source_root": source_root,
        "aggregate_root": aggregate_root,
        "audit_root": audit_root,
        "aggregate_manifest": aggregate_manifest,
        "aggregate_summary": aggregate_summary,
        "audit_manifest": audit_manifest,
        "audit_summary": audit_summary,
        "aggregate_manifest_sha256": AGGREGATE_MANIFEST_SHA256,
        "aggregate_summary_sha256": aggregate_summary_sha,
        "audit_manifest_sha256": AUDIT_MANIFEST_SHA256,
        "audit_summary_sha256": AUDIT_SUMMARY_SHA256,
    }


def _find_run_case(source_root: Path, run_id: str, mode: str) -> Path:
    root = source_root / "collector" / str(run_id) / "results" / "relaxtime" / "transport" / "phase_guided"
    candidates = list(root.glob(f"{mode}/{CASE_NAME}"))
    if len(candidates) != 1:
        raise ValueError(f"expected one collector case for run {run_id} {mode}, found {candidates}")
    return candidates[0]


def build_source_run_index(source: dict[str, Any]) -> list[dict[str, Any]]:
    aggregate_root: Path = source["aggregate_root"]
    _, fields, rows = read_csv(aggregate_root / "per_shard_summary.csv")
    _require_fields(fields, ("run_id", "mode", "label", "scan_rows", "diag_rows", "failed_rows", "config_hash", "calc_sha_in_audit", "workflow_sha_in_audit"), "per_shard_summary")
    if len(rows) != 30:
        raise ValueError(f"per_shard_summary rows {len(rows)} != 30")
    action_rows = read_csv(aggregate_root / "action_run_summary.csv")[2]
    action_by_run = {str(row.get("run_id")): row for row in action_rows}
    result: list[dict[str, Any]] = []
    for row in rows:
        run_id = str(row["run_id"])
        mode = row["mode"]
        action = action_by_run.get(run_id)
        if action is None or action.get("status") != "completed" or action.get("conclusion") != "success":
            raise ValueError(f"run {run_id} is not a selected successful Actions run")
        if row.get("hard_gate_ok", "").lower() != "true" or row.get("scan_all_converged", "").lower() != "true":
            raise ValueError(f"run {run_id} does not satisfy shard hard gates")
        if int(row.get("failed_rows", "0")) != 0:
            raise ValueError(f"run {run_id} has failed rows in shard summary")
        if row.get("calc_sha_in_audit", "").lower() != "true" or row.get("workflow_sha_in_audit", "").lower() != "true":
            raise ValueError(f"run {run_id} source SHA audit failed")
        case = _find_run_case(source["source_root"], run_id, mode)
        required = {
            name: case / name
            for name in (
                "channel_diagnostics.csv",
                "effective_config.json",
                "failed_points.csv",
                "README.md",
                "run_manifest.json",
                "sampling_plan.csv",
            )
        }
        if any(not path.is_file() for path in required.values()):
            raise ValueError(f"run {run_id} is missing a required sidecar")
        run_manifest = read_json(required["run_manifest.json"])
        if run_manifest.get("git_commit") != WORKFLOW_HEAD_SHA:
            raise ValueError(f"run {run_id} workflow head mismatch")
        if run_manifest.get("summary", {}).get("error_count") != 0:
            raise ValueError(f"run {run_id} has solver/error count")
        failed_rows = read_csv(required["failed_points.csv"])[2]
        if failed_rows:
            raise ValueError(f"run {run_id} has failed points")
        result.append(
            {
                "run_id": run_id,
                "mode": mode,
                "label": row.get("label", ""),
                "url": action.get("url", ""),
                "scan_rows": int(row["scan_rows"]),
                "diag_rows": int(row["diag_rows"]),
                "failed_rows": int(row["failed_rows"]),
                "config_hash": row.get("config_hash", ""),
                "source_case_root": str(case),
                "run_manifest": run_manifest,
                "sidecar_hashes": {name: sha256_file(path) for name, path in required.items()},
                "known_stale_manifest_hash_files": row.get("manifest_hash_mismatch_files", ""),
            }
        )
    if Counter(row["mode"] for row in result) != Counter({MODES[0]: 15, MODES[1]: 15}):
        raise ValueError("expected 15 selected shards per mode")
    return result


def validate_figure_reference(repo_root: Path, mode: str, scan_sha256: str) -> dict[str, Any]:
    figure_root = repo_root / "data" / "outputs" / "figures" / "relaxtime" / "transport" / "phase_guided" / mode / CASE_NAME
    manifest_path = figure_root / "plot_manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(manifest_path)
    manifest = read_json(manifest_path)
    if manifest.get("source_csv_sha256") != scan_sha256:
        raise ValueError(f"{mode} figure source CSV hash does not match aggregate scan")
    figure_hashes = manifest.get("figure_hashes", {})
    if manifest.get("figure_count") != 36 or len(figure_hashes) != 36:
        raise ValueError(f"{mode} figure count is not 36")
    for relative, expected in figure_hashes.items():
        path = figure_root / Path(relative)
        if not path.is_file() or sha256_file(path) != expected:
            raise ValueError(f"{mode} figure hash mismatch: {relative}")
    return {
        "path": figure_root.relative_to(repo_root).as_posix(),
        "plot_manifest_sha256": sha256_file(manifest_path),
        "figure_count": len(figure_hashes),
        "figure_status": manifest.get("figure_status"),
        "numerical_status": manifest.get("numerical_status"),
    }


def _unique_numeric(rows: list[dict[str, str]], field: str, *, finite_only: bool = True) -> list[float]:
    values: set[float] = set()
    for row in rows:
        value = parse_float(row.get(field))
        if value is None or (finite_only and not math.isfinite(value)):
            continue
        values.add(round(value, 12))
    return sorted(values)


def build_sampling_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [{field: row.get(field, "") for field in PLAN_FIELDS} for row in rows]


def build_effective_config(
    mode: str,
    rows: list[dict[str, str]],
    source_runs: list[dict[str, Any]],
    source: dict[str, Any],
) -> dict[str, Any]:
    first = source_runs[0]["run_manifest"]
    return {
        "schema_version": "issue130_rs_aggregate_effective_config_v1",
        "script": "scripts/relaxtime/run_phase_guided_transport_scan.jl",
        "case_name": CASE_NAME,
        "mode": mode,
        "options": {
            "muB_values": _unique_numeric(rows, "muB_MeV"),
            "T_values": _unique_numeric(rows, "T_MeV"),
            "alpha_T_values": _unique_numeric(rows, "alpha_T"),
            "xi_values": _unique_numeric(rows, "xi"),
            "propagator_xi_policy": "match_thermo",
            "sigma_cache_policy": "validated_anchored",
            "phase_anchor_policy": "direct_coexistence",
            "channel_diagnostics": True,
            "compute_bulk": True,
            "tau_p_nodes": int(source_runs[0]["run_manifest"]["argv"][source_runs[0]["run_manifest"]["argv"].index("--tau-p-nodes") + 1]),
            "tau_angle_nodes": int(source_runs[0]["run_manifest"]["argv"][source_runs[0]["run_manifest"]["argv"].index("--tau-angle-nodes") + 1]),
            "tau_phi_nodes": int(source_runs[0]["run_manifest"]["argv"][source_runs[0]["run_manifest"]["argv"].index("--tau-phi-nodes") + 1]),
            "tau_n_sigma_points": int(source_runs[0]["run_manifest"]["argv"][source_runs[0]["run_manifest"]["argv"].index("--tau-n-sigma") + 1]),
            "sigma_grid_n": int(source_runs[0]["run_manifest"]["argv"][source_runs[0]["run_manifest"]["argv"].index("--sigma-grid-n") + 1]),
        },
        "source": {
            "calculation_sha": CALCULATION_SHA,
            "workflow_head_sha": WORKFLOW_HEAD_SHA,
            "aggregate_name": AGGREGATE_NAME,
            "aggregate_manifest_sha256": source["aggregate_manifest_sha256"],
            "selected_shards": [row["run_id"] for row in source_runs],
            "source_solver_called": True,
            "aggregate_replay_solver_called": False,
            "import_solver_called": False,
            "derived_fields": ["sampling_plan.csv", "failed_points.csv", "effective_config.json"],
            "run_manifest_schema": first.get("schema_version"),
        },
    }


def build_import_manifest(
    case_dir: Path,
    repo_root: Path,
    mode: str,
    source: dict[str, Any],
    source_runs: list[dict[str, Any]],
    scan_summary: dict[str, Any],
    diag_summary: dict[str, Any],
    figure_summary: dict[str, Any],
    source_files: list[str],
    legacy_prod_v1_tree_hash: str,
) -> dict[str, Any]:
    result_files = sorted(
        [
            "phase_guided_transport_scan.csv",
            "sampling_plan.csv",
            "channel_diagnostics.csv",
            "failed_points.csv",
            "effective_config.json",
            "README.md",
            "PRODUCTION_AUDIT.md",
            *[f"convergence/{name}" for name in source_files],
            "convergence/post_repair_audit_pointer.json",
            "convergence/source_run_manifest_index.json",
            "manifest.json",
        ]
    )
    hashes = {
        relative: sha256_file(case_dir / relative)
        for relative in result_files
        if (case_dir / relative).is_file()
    }
    return {
        "schema_version": "issue130_rs_candidate_runtime_result_import_v1",
        "case_slug": CASE_NAME,
        "mode": mode,
        "artifact_status": "imported_candidate",
        "numerical_status": "diagnostic_only",
        "registry_status": "versioned_candidate_not_default",
        "runtime_default_switch": False,
        "legacy_fallback": True,
        "explicit_rollback": "--phase-reference-mode legacy",
        "calculation_sha": CALCULATION_SHA,
        "workflow_head_sha": WORKFLOW_HEAD_SHA,
        "source_root": str(source["source_root"]),
        "aggregate_root": AGGREGATE_NAME,
        "aggregate_manifest_sha256": source["aggregate_manifest_sha256"],
        "aggregate_summary_sha256": source["aggregate_summary_sha256"],
        "aggregate_source_file_hashes": {
            item["path"]: item["sha256"]
            for item in source["aggregate_manifest"].get("files", [])
        },
        "post_repair_audit": {
            "root": AUDIT_NAME,
            "manifest_sha256": source["audit_manifest_sha256"],
            "summary_sha256": source["audit_summary_sha256"],
            "verdict": source["audit_summary"]["verdict"],
        },
        "source_solver_called": True,
        "aggregate_replay_solver_called": False,
        "import_solver_called": False,
        "production_write": False,
        "legacy_prod_v1_tree_hash": legacy_prod_v1_tree_hash,
        "source_runs": [
            {"run_id": row["run_id"], "mode": row["mode"], "label": row["label"], "url": row["url"]}
            for row in source_runs
        ],
        "scan_summary": scan_summary,
        "diagnostic_summary": diag_summary,
        "figure_reference": figure_summary,
        "direct_coexistence_gate": source["aggregate_summary"]["direct_coexistence_gate"],
        "quality_reason_counts": source["aggregate_summary"]["quality_reason_counts"][mode],
        "known_metadata_defects": source["audit_summary"].get("known_historical_mismatch_counts", {}),
        "non_goals": [
            "no solver rerun",
            "no runtime default switch",
            "no old prod_v1 rewrite",
            "no legacy deletion",
            "no claim of numerical parity or convergence promotion",
        ],
        "result_files": result_files,
        "hashes": hashes,
        "figure_files": [
            f"{figure_summary['path']}/plot_manifest.json",
            f"{figure_summary['path']}/**/*.png",
        ],
    }


def write_readme(case_dir: Path, mode: str, scan_summary: dict[str, Any], source: dict[str, Any]) -> None:
    scope = "fixed mu_B with phase-scaled temperature bands" if mode == MODES[0] else "fixed temperature with sparse mu_B values"
    text = f"""# {CASE_NAME} ({mode})

This is the versioned Issue #130 RS candidate result imported from the
author-reviewed numerical aggregate. It is not the default runtime result.

## Status

- artifact status: `imported_candidate`
- numerical status: `diagnostic_only`
- runtime default switch: `false`
- legacy fallback: retained; explicit rollback is `--phase-reference-mode legacy`
- source calculation SHA: `{CALCULATION_SHA}`
- source workflow head: `{WORKFLOW_HEAD_SHA}`

## Scope

- mode: `{mode}` ({scope})
- aggregate scan rows: `{scan_summary['rows']}`
- selected source shards: `{len([row for row in source['source_runs'] if row['mode'] == mode])}`
- direct-coexistence policy: `xi=+/-0.003` retained for the mode-A `muB=900, alpha_T=1` anchor

## Provenance

- immutable aggregate: `{AGGREGATE_NAME}`
- aggregate manifest SHA256: `{source['aggregate_manifest_sha256']}`
- post-repair audit: `{AUDIT_NAME}`
- post-repair verdict: `post_repair_audit_pass_diagnostic_only`
- source numerical solver was used; this import and aggregate replay were solver-free

## Quality and limitations

Common `tau_u_ubar_ratio_high` warnings are retained in the scan and manifest.
They are diagnostic warnings, not silently filtered failures. The known
historical sidecar hash defect remains recorded; PR #269 fixed the producer
contract but did not rewrite immutable historical artifacts.

This case does not switch runtime consumers, delete `prod_v1`, or claim that
candidate/legacy numerical differences are a convergence result.
"""
    case_dir.joinpath("README.md").write_text(text, encoding="utf-8", newline="\n")


def write_audit(case_dir: Path, mode: str, scan_summary: dict[str, Any], diag_summary: dict[str, Any], source: dict[str, Any]) -> None:
    text = f"""# RS Candidate Result Import Audit ({mode})

## Verdict

`imported_candidate_diagnostic_only`

## Input gates

- calculation SHA: `{CALCULATION_SHA}`
- workflow head: `{WORKFLOW_HEAD_SHA}`
- aggregate manifest: `{source['aggregate_manifest_sha256']}`
- post-repair audit manifest: `{source['audit_manifest_sha256']}`
- post-repair audit verdict: `{source['audit_summary']['verdict']}`
- scan rows: `{scan_summary['rows']}`
- diagnostic rows: `{diag_summary['rows']}`
- scan duplicate keys: `{scan_summary['duplicate_keys']}`
- diagnostic duplicate keys: `{diag_summary['duplicate_keys']}`
- solver/failed rows: `0` in selected source shards

## Provenance semantics

- `source_solver_called=true`: the selected numerical Actions produced the source data.
- `aggregate_replay_solver_called=false`: aggregate replay did not call the solver.
- `import_solver_called=false`: this repository import only copied and derived sidecars.
- `production_write=false`: no runtime/default production switch occurred.

## Preserved boundaries

- old `prod_v1` result trees are not modified;
- legacy phase-reference fallback and explicit rollback remain available;
- candidate/legacy differences remain input-difference diagnostics, not numerical parity;
- quality warnings and historical sidecar hash defects are preserved in `convergence/`.
"""
    case_dir.joinpath("PRODUCTION_AUDIT.md").write_text(text, encoding="utf-8", newline="\n")


def import_case(repo_root: Path, source: dict[str, Any], mode: str, source_runs: list[dict[str, Any]]) -> dict[str, Any]:
    aggregate_root: Path = source["aggregate_root"]
    scan_path = aggregate_root / SCAN_FILES[mode]
    diag_path = aggregate_root / DIAG_FILES[mode]
    _, scan_fields, scan_rows = read_csv(scan_path)
    _, diag_fields, diag_rows = read_csv(diag_path)
    scan_summary = validate_scan_rows(mode, scan_fields, scan_rows, expected_rows=EXPECTED_SCAN_ROWS[mode])
    diag_summary = validate_diag_rows(mode, diag_fields, diag_rows, expected_rows=EXPECTED_DIAG_ROWS[mode])
    figure_summary = validate_figure_reference(repo_root, mode, sha256_file(scan_path))
    official_root = (repo_root / "data" / "outputs" / "results" / "relaxtime" / "transport" / "phase_guided").resolve()
    target = official_root / mode / CASE_NAME
    old_case = (
        official_root
        / LEGACY_SNAPSHOT_VERSION
        / mode
        / CASE_NAME.replace("_prod_v2", "_prod_v1")
    )
    if target.exists():
        raise FileExistsError(f"refusing to overwrite existing result case: {target}")
    if not old_case.is_dir():
        raise FileNotFoundError(f"legacy prod_v1 result missing: {old_case}")
    old_before = tree_hash(old_case)
    source_files = list(CONVERGENCE_COPY_FILES)
    temp_parent = official_root / mode
    temp_parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{CASE_NAME}.", dir=temp_parent))
    try:
        case_dir = staging / CASE_NAME
        case_dir.mkdir(parents=True, exist_ok=False)
        shutil.copy2(scan_path, case_dir / "phase_guided_transport_scan.csv")
        shutil.copy2(diag_path, case_dir / "channel_diagnostics.csv")
        write_csv(
            case_dir / "sampling_plan.csv",
            [
                f"# aggregate-derived sampling plan; source={AGGREGATE_NAME}",
                f"# calculation_sha={CALCULATION_SHA}; workflow_head_sha={WORKFLOW_HEAD_SHA}",
            ],
            PLAN_FIELDS,
            build_sampling_rows(scan_rows),
        )
        write_csv(
            case_dir / "failed_points.csv",
            [
                f"# aggregate-derived empty failure list; source={AGGREGATE_NAME}",
                f"# selected source shards={len(source_runs)}; failed rows=0",
            ],
            FAILED_FIELDS,
            [],
        )
        write_json(case_dir / "effective_config.json", build_effective_config(mode, scan_rows, source_runs, source))
        convergence = case_dir / "convergence"
        convergence.mkdir(parents=True, exist_ok=True)
        for name in source_files:
            shutil.copy2(aggregate_root / name, convergence / name)
        write_json(
            convergence / "post_repair_audit_pointer.json",
            {
                "schema_version": "issue130_rs_post_repair_audit_pointer_v1",
                "source_root": str(source["source_root"]),
                "audit_root": AUDIT_NAME,
                "manifest_sha256": source["audit_manifest_sha256"],
                "summary_sha256": source["audit_summary_sha256"],
                "verdict": source["audit_summary"]["verdict"],
            },
        )
        write_json(
            convergence / "source_run_manifest_index.json",
            {
                "schema_version": "issue130_rs_source_run_manifest_index_v1",
                "calculation_sha": CALCULATION_SHA,
                "workflow_head_sha": WORKFLOW_HEAD_SHA,
                "mode": mode,
                "selected_runs": source_runs,
            },
        )
        write_readme(case_dir, mode, scan_summary, {**source, "source_runs": source_runs})
        write_audit(case_dir, mode, scan_summary, diag_summary, source)
        manifest = build_import_manifest(
            case_dir,
            repo_root,
            mode,
            source,
            source_runs,
            scan_summary,
            diag_summary,
            figure_summary,
            source_files,
            old_before,
        )
        write_json(case_dir / "manifest.json", manifest)
        if tree_hash(old_case) != old_before:
            raise RuntimeError(f"legacy prod_v1 changed while staging {mode}")
        case_dir.rename(target)
        if tree_hash(old_case) != old_before:
            raise RuntimeError(f"legacy prod_v1 changed during import {mode}")
        return {"mode": mode, "target": target.relative_to(repo_root).as_posix(), "manifest_sha256": sha256_file(target / "manifest.json"), "scan_summary": scan_summary, "diagnostic_summary": diag_summary, "figure_summary": figure_summary}
    finally:
        if staging.exists():
            shutil.rmtree(staging, ignore_errors=True)


def tree_hash(root: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(path for path in root.rglob("*") if path.is_file()):
        digest.update(path.relative_to(root).as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(sha256_file(path).encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--validate-only", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    repo_root = args.repo_root.resolve()
    source_root = args.source_root.resolve()
    source = validate_external_source(source_root)
    source_runs = build_source_run_index(source)
    results = []
    for mode in MODES:
        _, fields, rows = read_csv(source["aggregate_root"] / SCAN_FILES[mode])
        validate_scan_rows(mode, fields, rows, expected_rows=EXPECTED_SCAN_ROWS[mode])
        if args.validate_only:
            results.append({"mode": mode, "ready": True, "scan_rows": len(rows), "source_shards": len([r for r in source_runs if r["mode"] == mode])})
        else:
            results.append(import_case(repo_root, source, mode, [r for r in source_runs if r["mode"] == mode]))
    print(json.dumps({"ready": True, "validate_only": args.validate_only, "calculation_sha": CALCULATION_SHA, "workflow_head_sha": WORKFLOW_HEAD_SHA, "results": results}, ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
