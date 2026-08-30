#!/usr/bin/env python3
"""Collect and gate the PNJL CEP production-shadow artifacts.

The collector treats missing files, duplicate matrix keys, SHA/provenance
drift, and non-finite final points as workflow-contract failures.  Physical
or performance-gate failures remain diagnostic verdicts and never promote a
reference.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

XIS = {-0.5, 0.0, 0.5}
METHODS = {"production_cascade", "memoized_dense", "independent_oracle"}
REQUIRED = ("curve_points.csv", "slice_metrics.csv", "method_costs.csv", "cep_accuracy.csv", "job_summary.json")
ANCHORS = {
    -0.5: (5.0, 20.0, 60.0, 100.0, 130.0, 147.0947265625, 147.2197265625, 160.0),
    0.0: (5.0, 20.0, 60.0, 100.0, 120.0, 130.9619140625, 131.0869140625, 145.0),
    0.5: (5.0, 20.0, 60.0, 90.0, 100.0, 106.9599609375, 107.0849609375, 120.0),
}


def _json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("\n", encoding="utf-8")
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _float(value: Any, default: float = math.nan) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _find_jobs(input_dir: Path) -> list[tuple[Path, dict[str, Any]]]:
    jobs = []
    for path in sorted(input_dir.rglob("job_summary.json")):
        summary = _json(path)
        jobs.append((path.parent, summary))
    return jobs


def _validate_jobs(jobs: list[tuple[Path, dict[str, Any]]], expected_sha: str) -> list[str]:
    errors: list[str] = []
    seen: set[tuple[float, str]] = set()
    for directory, summary in jobs:
        key = (_float(summary.get("xi")), str(summary.get("method", "")))
        if key in seen:
            errors.append(f"duplicate matrix key {key}")
        seen.add(key)
        if key[0] not in XIS or key[1] not in METHODS:
            errors.append(f"unexpected matrix key {key}")
        if summary.get("schema_version") != "cep_cascade_production_shadow_v1":
            errors.append(f"unexpected schema at {directory}")
        for name in REQUIRED:
            if not (directory / name).is_file():
                errors.append(f"missing {name} at {directory}")
        sha = str(summary.get("calculation_sha", ""))
        if sha != expected_sha:
            errors.append(f"calculation SHA mismatch at {directory}: {sha}")
        if summary.get("workflow_head_sha") != expected_sha:
            errors.append(f"workflow head SHA mismatch at {directory}: {summary.get('workflow_head_sha', '')}")
        provenance = summary.get("provenance", {})
        if (
            provenance.get("calculation_sha") != expected_sha
            or provenance.get("reference_write") is not False
            or provenance.get("anchor_run_count") != len(ANCHORS.get(key[0], ()))
        ):
            errors.append(f"invalid provenance at {directory}")
        if not bool(summary.get("finite_and_converged_final", False)):
            errors.append(f"non-finite/non-converged final points at {directory}")
        expected_anchors = ANCHORS.get(key[0], ())
        summary_anchors = tuple(_float(value) for value in summary.get("anchors", ()))
        if len(summary_anchors) != len(expected_anchors) or any(
            not math.isclose(actual, expected, rel_tol=0.0, abs_tol=1e-9)
            for actual, expected in zip(summary_anchors, expected_anchors)
        ):
            errors.append(f"anchor declaration mismatch at {directory}")
        slice_rows = _rows(directory / "slice_metrics.csv") if (directory / "slice_metrics.csv").is_file() else []
        accuracy_rows = _rows(directory / "cep_accuracy.csv") if (directory / "cep_accuracy.csv").is_file() else []
        cost_rows = _rows(directory / "method_costs.csv") if (directory / "method_costs.csv").is_file() else []
        curve_rows = _rows(directory / "curve_points.csv") if (directory / "curve_points.csv").is_file() else []
        for table_rows, table_name in ((slice_rows, "slice_metrics.csv"), (accuracy_rows, "cep_accuracy.csv"), (cost_rows, "method_costs.csv")):
            for row in table_rows:
                if row.get("method") not in (None, "", key[1]):
                    errors.append(f"{table_name} method mismatch at {directory}")
                if row.get("xi") not in (None, "") and not math.isclose(_float(row.get("xi")), key[0], rel_tol=0.0, abs_tol=1e-9):
                    errors.append(f"{table_name} xi mismatch at {directory}")
        if len(slice_rows) != len(expected_anchors):
            errors.append(f"slice_metrics anchor count {len(slice_rows)} != {len(expected_anchors)} at {directory}")
        if len(accuracy_rows) != len(expected_anchors):
            errors.append(f"cep_accuracy anchor count {len(accuracy_rows)} != {len(expected_anchors)} at {directory}")
        if len(cost_rows) != 1:
            errors.append(f"method_costs row count {len(cost_rows)} != 1 at {directory}")
        elif cost_rows[0].get("method") not in (None, "", key[1]):
            errors.append(f"method_costs method mismatch at {directory}")
        if not curve_rows:
            errors.append(f"curve_points is empty at {directory}")
        curve_keys: set[tuple[float, float, float, str]] = set()
        for curve in curve_rows:
            curve_key = (_float(curve.get("T_MeV")), _float(curve.get("rho")), _float(curve.get("xi")), str(curve.get("method", "")))
            if curve_key in curve_keys:
                errors.append(f"duplicate curve point key {curve_key} at {directory}")
            curve_keys.add(curve_key)
            if not (str(curve.get("converged", "")).lower() == "true" and str(curve.get("finite", "")).lower() == "true"):
                errors.append(f"non-finite/non-converged curve point {curve_key} at {directory}")
        declared_hashes = summary.get("curve_file_sha256", {})
        for name in ("curve_points.csv", "slice_metrics.csv", "method_costs.csv", "cep_accuracy.csv"):
            if (directory / name).is_file():
                actual_hash = hashlib.sha256((directory / name).read_bytes()).hexdigest()
                if declared_hashes.get(name) != actual_hash:
                    errors.append(f"file hash mismatch for {name} at {directory}")
        actual_slice_anchors = sorted({_float(row.get("T_MeV")) for row in slice_rows})
        if len(actual_slice_anchors) != len(expected_anchors) or any(
            not any(math.isclose(actual, expected, rel_tol=0.0, abs_tol=1e-5) for expected in expected_anchors)
            for actual in actual_slice_anchors
        ):
            errors.append(f"slice_metrics anchors mismatch at {directory}: {actual_slice_anchors}")
        for row in slice_rows:
            if key[1] == "production_cascade" and _float(row.get("targeted_additions", 0.0)) > 12:
                errors.append(f"targeted cap exceeded at {directory}, T={row.get('T_MeV')}")
    expected = {(xi, method) for xi in XIS for method in METHODS}
    if seen != expected:
        errors.append(f"matrix incomplete: missing={sorted(expected - seen)} extra={sorted(seen - expected)}")
    return errors


def _normalized_cost_row(directory: Path, summary: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    """Normalize per-job cache counters without rewriting source artifacts.

    Historical runners wrote cache counters from the first anchor only.  The
    slice table contains all anchor-level counters; fixed-rho telemetry closes
    the possible one-off probe gap.  The derived row is written only to the
    aggregate output and the raw job hash remains validated separately.
    """
    raw_rows = _rows(directory / "method_costs.csv")
    raw = dict(raw_rows[0])
    slice_rows = _rows(directory / "slice_metrics.csv")
    summed_unique = sum(_float(row.get("unique_solves"), 0.0) for row in slice_rows)
    summed_requests = sum(_float(row.get("point_requests"), 0.0) for row in slice_rows)
    summed_cache = sum(_float(row.get("cache_hits"), 0.0) for row in slice_rows)
    summed_targeted = sum(_float(row.get("targeted_additions"), 0.0) for row in slice_rows)
    fixedrho = _float(raw.get("fixedrho_requests"), 0.0)
    retries = _float(raw.get("retry_count"), 0.0)
    fallbacks = _float(raw.get("fallback_count"), 0.0)
    unattributed = max(fixedrho - summed_unique, 0.0) if retries == 0.0 and fallbacks == 0.0 else 0.0
    total_unique = summed_unique + unattributed
    total_requests = summed_requests + unattributed
    normalized = dict(raw)
    normalized.update({
        "unique_solves": int(round(total_unique)),
        "requested_point_calls": int(round(total_requests)),
        "uncached_equivalent_requests": int(round(total_requests)),
        "cache_hits": int(round(summed_cache)),
        "targeted_additions": int(round(summed_targeted)),
        "unattributed_unique_solves": int(round(unattributed)),
        "unattributed_point_requests": int(round(unattributed)),
        "aggregation_scope": "all_anchors",
        "point_request_reconciliation": str(
            math.isclose(total_requests, total_unique + summed_cache, rel_tol=0.0, abs_tol=1e-9)
        ).lower(),
    })
    corrections = {
        "xi": summary.get("xi", ""),
        "method": summary.get("method", ""),
        "raw_unique_solves": raw.get("unique_solves", ""),
        "normalized_unique_solves": normalized["unique_solves"],
        "raw_requested_point_calls": raw.get("requested_point_calls", ""),
        "normalized_requested_point_calls": normalized["requested_point_calls"],
        "raw_cache_hits": raw.get("cache_hits", ""),
        "normalized_cache_hits": normalized["cache_hits"],
        "raw_targeted_additions": raw.get("targeted_additions", ""),
        "normalized_targeted_additions": normalized["targeted_additions"],
        "unattributed_unique_solves": normalized["unattributed_unique_solves"],
    }
    return normalized, corrections


def _collect_tables(jobs: list[tuple[Path, dict[str, Any]]], output_dir: Path) -> list[dict[str, Any]]:
    corrections: list[dict[str, Any]] = []
    for name in ("curve_points.csv", "slice_metrics.csv", "method_costs.csv", "cep_accuracy.csv"):
        rows: list[dict[str, Any]] = []
        for directory, summary in jobs:
            source_rows = _rows(directory / name)
            if name == "method_costs.csv":
                source_rows = [_normalized_cost_row(directory, summary)[0]]
                corrections.append(_normalized_cost_row(directory, summary)[1])
            for row in source_rows:
                row["xi"] = summary.get("xi", "")
                row["method"] = summary.get("method", "")
                row["calculation_sha"] = summary.get("calculation_sha", "")
                rows.append(row)
        _write_csv(output_dir / name, rows)

    slice_rows = _rows(output_dir / "slice_metrics.csv")
    _write_csv(output_dir / "geometry_accuracy.csv", [
        {
            "xi": row.get("xi", ""),
            "method": row.get("method", ""),
            "T_MeV": row.get("T_MeV", ""),
            "coarse_status": row.get("coarse_status", ""),
            "fine_status": row.get("fine_status", ""),
            "result_status": row.get("result_status", ""),
            "geometry_converged": row.get("geometry_converged", ""),
            "position_error_MeV": row.get("position_error_MeV", ""),
            "density_error": row.get("density_error", ""),
            "maxwell_area_gate": row.get("maxwell_area_gate", ""),
            "area_residual": row.get("area_residual", ""),
            "rho_hadron": row.get("rho_hadron", ""),
            "rho_quark": row.get("rho_quark", ""),
            "mu_spinodal_hadron_MeV": row.get("mu_spinodal_hadron_MeV", ""),
            "mu_spinodal_quark_MeV": row.get("mu_spinodal_quark_MeV", ""),
            "rho_spinodal_hadron": row.get("rho_spinodal_hadron", ""),
            "rho_spinodal_quark": row.get("rho_spinodal_quark", ""),
        }
        for row in slice_rows
    ])
    _write_csv(output_dir / "curve_index.csv", [
        {
            "xi": summary.get("xi", ""),
            "method": summary.get("method", ""),
            "calculation_sha": summary.get("calculation_sha", ""),
            "source_job": str(directory),
            "curve_file": "curve_points.csv",
            "curve_sha256": summary.get("curve_file_sha256", {}).get("curve_points.csv", ""),
            "raw_curve_copy_in_repository": False,
        }
        for directory, summary in jobs
    ])
    return corrections


def _actions(run_id: str | None, output_dir: Path) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    metadata: dict[str, Any] = {"run_id": run_id or "", "source": "unavailable"}
    if run_id:
        try:
            payload = json.loads(subprocess.check_output(
                ["gh", "run", "view", run_id, "--json", "jobs,headSha,url,status,conclusion"],
                text=True,
            ))
            metadata = {"run_id": run_id, "source": "gh run view", **{k: payload.get(k, "") for k in ("headSha", "url", "status", "conclusion")}}
            for job in payload.get("jobs", []):
                if not job.get("startedAt") or not job.get("completedAt"):
                    continue
                start = datetime.fromisoformat(job["startedAt"].replace("Z", "+00:00"))
                end = datetime.fromisoformat(job["completedAt"].replace("Z", "+00:00"))
                seconds = max(0.0, (end - start).total_seconds())
                rows.append({"job_name": job.get("name", ""), "job_id": job.get("databaseId", ""), "status": job.get("status", ""), "conclusion": job.get("conclusion", ""), "elapsed_seconds": seconds, "runner_minutes_rounded": math.ceil(seconds / 60.0)})
        except (OSError, subprocess.CalledProcessError, json.JSONDecodeError, ValueError) as exc:
            metadata["error"] = str(exc)
    _write_csv(output_dir / "actions_costs.csv", rows)
    elapsed = [float(row["elapsed_seconds"]) for row in rows]
    metadata.update({"job_count": len(rows), "critical_path_seconds": max(elapsed, default=0.0), "raw_total_seconds": sum(elapsed), "runner_minutes": sum(int(row["runner_minutes_rounded"]) for row in rows)})
    return metadata


def _performance(cost_rows: list[dict[str, str]]) -> tuple[list[str], dict[str, Any]]:
    grouped: dict[str, dict[str, float]] = {}
    for row in cost_rows:
        group = grouped.setdefault(row.get("method", ""), {})
        for key in ("equilibrium_requests", "fixedrho_requests", "unique_solves", "residual_calls", "jacobian_calls", "newton_iterations", "runner_seconds", "fallback_count", "retry_count"):
            group[key] = group.get(key, 0.0) + _float(row.get(key, 0.0), 0.0)
    errors: list[str] = []
    cascade = grouped.get("production_cascade", {})
    dense = grouped.get("memoized_dense", {})
    for key in ("equilibrium_requests", "fixedrho_requests", "unique_solves", "residual_calls", "jacobian_calls", "newton_iterations", "runner_seconds"):
        allowance = 1.10 if key == "runner_seconds" else 1.0
        if cascade.get(key, 0.0) > allowance * dense.get(key, 0.0) + 1e-9:
            errors.append(f"cascade {key} exceeds memoized dense")
    cascade_rate = (cascade.get("fallback_count", 0.0) + cascade.get("retry_count", 0.0)) / max(cascade.get("unique_solves", 0.0), 1.0)
    dense_rate = (dense.get("fallback_count", 0.0) + dense.get("retry_count", 0.0)) / max(dense.get("unique_solves", 0.0), 1.0)
    if (dense_rate == 0.0 and cascade_rate > 0.0) or (dense_rate > 0.0 and cascade_rate > 1.25 * dense_rate):
        errors.append("cascade fallback/retry rate exceeds the 25% risk threshold")
    return errors, {"grouped": grouped, "cascade_fallback_retry_rate": cascade_rate, "dense_fallback_retry_rate": dense_rate}


def _gate(jobs: list[tuple[Path, dict[str, Any]]], output_dir: Path, contract_errors: list[str]) -> dict[str, Any]:
    slice_rows = _rows(output_dir / "slice_metrics.csv")
    costs = _rows(output_dir / "method_costs.csv")
    oracle_rows = [row for row in slice_rows if row.get("method") == "independent_oracle"]
    oracle_errors: list[str] = []
    if len(oracle_rows) != 24:
        oracle_errors.append("oracle anchors are incomplete")
    if any(_float(row.get("solver_failure_count", 0), 0.0) != 0 for row in oracle_rows):
        oracle_errors.append("oracle anchors have solver failures")
    for row in oracle_rows:
        status = row.get("result_status", "")
        if status != "ambiguous_near_critical" and str(row.get("geometry_converged", "")).lower() != "true":
            oracle_errors.append(f"oracle geometry is not converged at xi={row.get('xi')} T={row.get('T_MeV')}")
        if status == "confirmed_first_order" and (
            row.get("coarse_status") != "valid" or row.get("fine_status") != "valid"
            or row.get("coarse_reason") != "ok" or row.get("fine_reason") != "ok"
        ):
            oracle_errors.append(f"oracle first-order two-layer certificate is incomplete at xi={row.get('xi')} T={row.get('T_MeV')}")
        if status == "confirmed_monotone" and (
            row.get("coarse_status") != "invalid" or row.get("fine_status") != "invalid"
            or row.get("coarse_reason") != "no_s_shape" or row.get("fine_reason") != "no_s_shape"
        ):
            oracle_errors.append(f"oracle monotone two-layer certificate is incomplete at xi={row.get('xi')} T={row.get('T_MeV')}")
    for xi in sorted(XIS):
        statuses = {row.get("result_status", "") for row in oracle_rows if math.isclose(_float(row.get("xi")), xi, rel_tol=0.0, abs_tol=1e-9)}
        if "confirmed_first_order" not in statuses or "confirmed_monotone" not in statuses:
            oracle_errors.append(f"oracle xi={xi} lacks confirmed first-order/monotone two-sided evidence")
    cascade_errors = []
    cascade_low_temperature_failure = False
    cascade_non_low_temperature_failure = False
    by_anchor: dict[tuple[str, str], dict[str, str]] = {}
    for row in slice_rows:
        by_anchor.setdefault((row.get("xi", ""), row.get("T_MeV", "")), {})[row.get("method", "")] = row.get("result_status", "")
    for key, methods in by_anchor.items():
        if "production_cascade" in methods and "independent_oracle" in methods:
            if methods["production_cascade"] != "ambiguous_near_critical" and methods["production_cascade"] != methods["independent_oracle"]:
                cascade_errors.append(f"classification mismatch at xi={key[0]}, T={key[1]}")
                if _float(key[1]) < 60.0:
                    cascade_low_temperature_failure = True
                else:
                    cascade_non_low_temperature_failure = True
    for row in slice_rows:
        if row.get("method") != "production_cascade" or str(row.get("result_status", "")) == "ambiguous_near_critical":
            continue
        if str(row.get("geometry_converged", "")).lower() != "true":
            if _float(row.get("T_MeV")) < 60.0:
                cascade_low_temperature_failure = True
            else:
                cascade_non_low_temperature_failure = True
    oracle_by_anchor = {
        (row.get("xi", ""), row.get("T_MeV", "")): row
        for row in slice_rows
        if row.get("method") == "independent_oracle"
    }
    for row in slice_rows:
        if row.get("method") != "production_cascade" or row.get("result_status") != "confirmed_first_order":
            continue
        oracle = oracle_by_anchor.get((row.get("xi", ""), row.get("T_MeV", "")))
        if oracle is None or oracle.get("result_status") != "confirmed_first_order":
            continue
        for field in ("mu_transition_MeV", "mu_spinodal_hadron_MeV", "mu_spinodal_quark_MeV"):
            cascade_value = _float(row.get(field))
            oracle_value = _float(oracle.get(field))
            if math.isfinite(cascade_value) and math.isfinite(oracle_value) and abs(cascade_value - oracle_value) > 0.025:
                cascade_errors.append(f"{field} differs from oracle by more than 0.025 MeV at xi={row.get('xi')} T={row.get('T_MeV')}")
        for field in ("rho_hadron", "rho_quark", "rho_spinodal_hadron", "rho_spinodal_quark"):
            cascade_value = _float(row.get(field))
            oracle_value = _float(oracle.get(field))
            if math.isfinite(cascade_value) and math.isfinite(oracle_value) and abs(cascade_value - oracle_value) > 0.0025:
                cascade_errors.append(f"{field} differs from oracle by more than 0.0025 at xi={row.get('xi')} T={row.get('T_MeV')}")
        cascade_area = _float(row.get("area_residual"))
        oracle_area = _float(oracle.get("area_residual"))
        if math.isfinite(cascade_area) and math.isfinite(oracle_area) and max(cascade_area, oracle_area) > 5e-5:
            cascade_errors.append(f"Maxwell area exceeds 5e-5 at xi={row.get('xi')} T={row.get('T_MeV')}")
    cascade_coverage_errors = []
    for xi in sorted(XIS):
        rows = [row for row in slice_rows if row.get("method") == "production_cascade" and math.isclose(_float(row.get("xi")), xi, rel_tol=0.0, abs_tol=1e-9)]
        statuses = {row.get("result_status", "") for row in rows}
        if "confirmed_first_order" not in statuses or "confirmed_monotone" not in statuses:
            cascade_coverage_errors.append(f"cascade xi={xi} lacks confirmed first-order/monotone two-sided evidence")
    cascade_errors.extend(cascade_coverage_errors)
    performance_errors, performance = _performance(costs)
    all_errors = contract_errors + oracle_errors + cascade_errors + performance_errors
    if contract_errors:
        verdict = "workflow_failure"
    elif oracle_errors:
        verdict = "oracle_inconclusive"
    elif performance_errors or cascade_non_low_temperature_failure:
        verdict = "integration_failed"
    elif cascade_low_temperature_failure or cascade_coverage_errors:
        verdict = "hybrid_required"
    elif not all_errors:
        verdict = "full_cascade_candidate"
    else:
        verdict = "integration_failed"
    return {"verdict": verdict, "workflow_contract_errors": contract_errors, "oracle_errors": oracle_errors, "cascade_errors": cascade_errors, "cascade_low_temperature_failure": cascade_low_temperature_failure, "cascade_non_low_temperature_failure": cascade_non_low_temperature_failure, "performance_errors": performance_errors, "performance": performance, "automatic_gate_is_not_promotion": True}


def _write_docs(output_dir: Path, gate: dict[str, Any], actions: dict[str, Any]) -> None:
    errors = gate["workflow_contract_errors"] + gate["oracle_errors"] + gate["cascade_errors"] + gate["performance_errors"]
    if gate.get("verdict") == "hybrid_required":
        errors = errors + ["仅低温 cascade 未通过；等待作者判断 hybrid 收口"]
    text = (
        f"# PNJL CEP cascade production shadow v1\n\n"
        f"verdict: `{gate['verdict']}`。这是 opt-in production integration 的诊断产物，"
        f"不覆盖 reference，不启动 C0/C1/C2 或 transport。\n\n"
        f"- oracle: {'stable' if not gate['oracle_errors'] else 'inconclusive'}\n"
        f"- Actions critical path: {actions.get('critical_path_seconds', 0.0)} s\n"
        f"- runner-minutes: {actions.get('runner_minutes', 0)}\n"
        f"- errors: {'；'.join(errors) if errors else '无'}\n\n"
        "`curve_points.csv` 保留原始 rho–mu 曲线索引，`slice_metrics.csv` 分开记录"
        " cascade/dense/oracle 的切片状态，`method_costs.csv` 记录 solver telemetry 和缓存。\n"
    )
    (output_dir / "README.md").write_text(text, encoding="utf-8")
    _write_csv(output_dir / "claim_ledger.csv", [
        {"claim_id": "oracle", "claim": "independent oracle anchors are stable", "status": "pass" if not gate["oracle_errors"] else "oracle_inconclusive", "boundary": "diagnostic only"},
        {"claim_id": "cascade", "claim": "production cascade classifications agree on non-ambiguous anchors", "status": "hybrid_required" if gate["verdict"] == "hybrid_required" else ("pass" if not gate["cascade_errors"] else "integration_failed"), "boundary": "author physical review required"},
        {"claim_id": "performance", "claim": "cascade solver work is no higher than memoized dense", "status": "pass" if not gate["performance_errors"] else "performance_risk", "boundary": "runner noise allowance applies"},
    ])
    (output_dir / "AUDIT.md").write_text(
        "# Shadow audit\n\n"
        "本目录只记录 production cascade 的 shadow 证据。原始 rho–mu 曲线留在"
        " Actions/local artifact，仓库内表格通过 job SHA 与 curve index 追溯；"
        " verdict 不自动晋升 reference。\n",
        encoding="utf-8",
    )
    (output_dir / "plot_manifest.json").write_text(
        json.dumps({"schema_version": "cep_cascade_production_shadow_v1", "figures": [], "reason": "plots are generated from the aggregated curve index"}, indent=2) + "\n",
        encoding="utf-8",
    )


def collect(
    input_dir: Path,
    output_dir: Path,
    run_id: str | None,
    expected_sha: str,
    *,
    postprocess_sha: str | None = None,
    source_run_id: str | None = None,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    jobs = _find_jobs(input_dir)
    contract_errors = _validate_jobs(jobs, expected_sha)
    corrections = _collect_tables(jobs, output_dir)
    actions = _actions(run_id, output_dir)
    gate = _gate(jobs, output_dir, contract_errors)
    _write_docs(output_dir, gate, actions)
    hashes = {}
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and path.name != "manifest.json":
            hashes[str(path.relative_to(output_dir)).replace(os.sep, "/")] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = {
        "schema_version": "cep_cascade_production_shadow_v1",
        "run_id": run_id or "",
        "source_run_id": source_run_id or run_id or "",
        "expected_calculation_sha": expected_sha,
        "postprocess_sha": postprocess_sha or "",
        "job_keys": sorted([[summary.get("xi"), summary.get("method")] for _, summary in jobs]),
        "gate": gate,
        "actions": actions,
        "aggregation_corrections": corrections,
        "file_sha256": hashes,
        "generated_at_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return gate


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--expected-calculation-sha", required=True)
    parser.add_argument("--postprocess-sha", default=None)
    parser.add_argument("--source-run-id", default=None)
    args = parser.parse_args()
    gate = collect(
        args.input_dir,
        args.output_dir,
        args.run_id,
        args.expected_calculation_sha,
        postprocess_sha=args.postprocess_sha,
        source_run_id=args.source_run_id,
    )
    print(json.dumps(gate, indent=2, sort_keys=True))
    return 1 if gate["verdict"] == "workflow_failure" else 2 if gate["verdict"] != "full_cascade_candidate" else 0


if __name__ == "__main__":
    raise SystemExit(main())
