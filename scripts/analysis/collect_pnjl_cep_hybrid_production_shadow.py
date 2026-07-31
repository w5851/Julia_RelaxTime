#!/usr/bin/env python3
"""Collect and gate the versioned CEP hybrid-production shadow evidence.

This collector is deliberately independent from the v1 cascade collector.  It
validates the nine-job matrix, provenance and support-bounded Stage-C points,
then emits a diagnostic verdict; no verdict promotes a reference.
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
METHODS = {"production_hybrid", "memoized_dense", "independent_oracle"}
REQUIRED = ("curve_points.csv", "slice_metrics.csv", "method_costs.csv", "cep_accuracy.csv", "job_summary.json")
ANCHORS = {
    -0.5: (5.0, 20.0, 60.0, 100.0, 130.0, 147.0947265625, 147.2197265625, 160.0),
    0.0: (5.0, 20.0, 60.0, 100.0, 120.0, 130.9619140625, 131.0869140625, 145.0),
    0.5: (5.0, 20.0, 60.0, 90.0, 100.0, 106.9599609375, 107.0849609375, 120.0),
}

# These anchors are deliberately stricter than the generic "each xi has both
# sides" check.  They are the strong first-order/high-temperature controls
# declared by the production-shadow contract.  Only the two CEP-neighbour
# anchors per xi may remain ambiguous.
REQUIRED_FIRST_ORDER = {
    -0.5: (5.0, 20.0, 60.0, 100.0, 130.0),
    0.0: (5.0, 20.0, 60.0, 100.0, 120.0),
    0.5: (5.0, 20.0, 60.0, 90.0, 100.0),
}
REQUIRED_MONOTONE = {
    -0.5: (160.0,),
    0.0: (145.0,),
    0.5: (120.0,),
}
CEP_NEIGHBOURS = {
    -0.5: (147.0947265625, 147.2197265625),
    0.0: (130.9619140625, 131.0869140625),
    0.5: (106.9599609375, 107.0849609375),
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
    return [(path.parent, _json(path)) for path in sorted(input_dir.rglob("job_summary.json"))]


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
        if summary.get("schema_version") != "cep_cascade_production_shadow_v2":
            errors.append(f"unexpected schema at {directory}")
        for name in REQUIRED:
            if not (directory / name).is_file():
                errors.append(f"missing {name} at {directory}")
        if summary.get("calculation_sha") != expected_sha or summary.get("workflow_head_sha") != expected_sha:
            errors.append(f"calculation/workflow SHA mismatch at {directory}")
        provenance = summary.get("provenance", {})
        if provenance.get("calculation_sha") != expected_sha or provenance.get("reference_write") is not False:
            errors.append(f"invalid provenance at {directory}")
        expected = ANCHORS.get(key[0], ())
        declared = tuple(_float(value) for value in summary.get("anchors", ()))
        if len(declared) != len(expected) or any(not math.isclose(a, b, abs_tol=1e-9, rel_tol=0.0) for a, b in zip(declared, expected)):
            errors.append(f"anchor declaration mismatch at {directory}")
        if not bool(summary.get("finite_and_converged_final", False)):
            errors.append(f"non-finite/non-converged final points at {directory}")
        for name in ("slice_metrics.csv", "cep_accuracy.csv", "method_costs.csv"):
            rows = _rows(directory / name) if (directory / name).is_file() else []
            if name == "slice_metrics.csv" and len(rows) != len(expected):
                errors.append(f"slice_metrics anchor count mismatch at {directory}")
            if name == "cep_accuracy.csv" and len(rows) != len(expected):
                errors.append(f"cep_accuracy anchor count mismatch at {directory}")
            if name == "method_costs.csv" and len(rows) != 1:
                errors.append(f"method_costs row count mismatch at {directory}")
            for row in rows:
                if row.get("method") not in (None, "", key[1]) or (
                    row.get("xi") not in (None, "") and not math.isclose(_float(row.get("xi")), key[0], abs_tol=1e-9, rel_tol=0.0)
                ):
                    errors.append(f"table provenance mismatch at {directory}/{name}")
        curves = _rows(directory / "curve_points.csv") if (directory / "curve_points.csv").is_file() else []
        curve_keys: set[tuple[float, float, float, str]] = set()
        for row in curves:
            key_curve = (_float(row.get("T_MeV")), _float(row.get("rho")), _float(row.get("xi")), str(row.get("method", "")))
            if key_curve in curve_keys:
                errors.append(f"duplicate curve point key {key_curve} at {directory}")
            curve_keys.add(key_curve)
            if str(row.get("converged", "")).lower() != "true" or str(row.get("finite", "")).lower() != "true":
                errors.append(f"non-finite/non-converged curve point {key_curve} at {directory}")
        if not curves:
            errors.append(f"curve_points is empty at {directory}")
        declared_hashes = summary.get("curve_file_sha256", {})
        for name in ("curve_points.csv", "slice_metrics.csv", "method_costs.csv", "cep_accuracy.csv"):
            path = directory / name
            if path.is_file() and declared_hashes.get(name) != hashlib.sha256(path.read_bytes()).hexdigest():
                errors.append(f"file hash mismatch for {name} at {directory}")
        if key[1] == "production_hybrid":
            for row in _rows(directory / "slice_metrics.csv") if (directory / "slice_metrics.csv").is_file() else []:
                if _float(row.get("targeted_additions"), 0.0) > 12:
                    errors.append(f"Stage-A targeted cap exceeded at {directory}, T={row.get('T_MeV')}")
        for row in _rows(directory / "slice_metrics.csv") if (directory / "slice_metrics.csv").is_file() else []:
            if row.get("stage_c_status") in (None, "", "not_run"):
                continue
            low, high = _float(row.get("support_low")), _float(row.get("support_high"))
            if not (math.isfinite(low) and math.isfinite(high) and low < high):
                errors.append(f"Stage-C support is missing at {directory}, T={row.get('T_MeV')}")
                continue
            for curve in _rows(directory / "curve_points.csv") if (directory / "curve_points.csv").is_file() else []:
                if curve.get("sampling_role") != "stage_c_support" or not math.isclose(_float(curve.get("T_MeV")), _float(row.get("T_MeV")), abs_tol=1e-8, rel_tol=0.0):
                    continue
                rho = _float(curve.get("rho"))
                if not (low - 1e-9 <= rho <= high + 1e-9):
                    errors.append(f"Stage-C point lies outside declared support at {directory}, rho={curve.get('rho')}")
    expected_keys = {(xi, method) for xi in XIS for method in METHODS}
    if seen != expected_keys:
        errors.append(f"matrix incomplete: missing={sorted(expected_keys - seen)} extra={sorted(seen - expected_keys)}")
    return errors


def _collect_tables(jobs: list[tuple[Path, dict[str, Any]]], output_dir: Path) -> list[str]:
    tables: dict[str, list[dict[str, Any]]] = {}
    for name in ("curve_points.csv", "slice_metrics.csv", "method_costs.csv", "cep_accuracy.csv"):
        rows: list[dict[str, Any]] = []
        for directory, summary in jobs:
            for row in _rows(directory / name):
                row.setdefault("xi", summary.get("xi", ""))
                row.setdefault("method", summary.get("method", ""))
                row.setdefault("calculation_sha", summary.get("calculation_sha", ""))
                rows.append(row)
        tables[name] = rows

    corrections: list[str] = []
    accuracy_by_key = {
        (_float(row.get("xi")), row.get("method", ""), _float(row.get("anchor_T_MeV"))): row
        for row in tables["cep_accuracy.csv"]
    }
    # The first shadow runner wrote the resolved Maxwell chemical potential
    # under mu_transition_MeV in sweep_records, while its CSV adapter looked
    # only for the legacy mu_transition field.  Rebuild the historical field
    # from the per-anchor CEP row during replay; future jobs are fixed in the
    # Julia adapter as well.
    for row in tables["slice_metrics.csv"]:
        if row.get("result_status") != "confirmed_first_order" or math.isfinite(_float(row.get("mu_transition_MeV"))):
            continue
        key = (_float(row.get("xi")), row.get("method", ""), _float(row.get("T_MeV")))
        accuracy = accuracy_by_key.get(key)
        repaired = _float(accuracy.get("muq_last_first_order_MeV")) if accuracy else math.nan
        if math.isfinite(repaired):
            row["mu_transition_MeV"] = f"{repaired:.17g}"
            row["mu_transition_source"] = "cep_accuracy_replay"
            corrections.append(f"reconstructed mu_transition_MeV for {key}")

    for name, rows in tables.items():
        _write_csv(output_dir / name, rows)
    slice_rows = tables["slice_metrics.csv"]
    _write_csv(output_dir / "geometry_accuracy.csv", [
        {key: row.get(key, "") for key in (
            "xi", "method", "T_MeV", "stage_a_status", "stage_b_status", "stage_c_status",
            "stage_used", "upgrade_reason", "result_status", "geometry_converged",
            "position_error_MeV", "density_error", "maxwell_area_gate", "area_residual",
            "mu_transition_MeV", "rho_hadron", "rho_quark", "mu_spinodal_hadron_MeV",
            "mu_spinodal_quark_MeV", "rho_spinodal_hadron", "rho_spinodal_quark",
            "support_low", "support_high", "support_point_count",
        )} for row in slice_rows
    ])
    _write_csv(output_dir / "curve_index.csv", [
        {
            "xi": summary.get("xi", ""), "method": summary.get("method", ""),
            "calculation_sha": summary.get("calculation_sha", ""),
            "source_job": str(directory), "curve_file": "curve_points.csv",
            "curve_sha256": summary.get("curve_file_sha256", {}).get("curve_points.csv", ""),
            "raw_curve_copy_in_repository": False,
        } for directory, summary in jobs
    ])
    return corrections


def _actions(run_id: str | None, output_dir: Path) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    metadata: dict[str, Any] = {"run_id": run_id or "", "source": "unavailable"}
    if run_id:
        try:
            payload = json.loads(subprocess.check_output(["gh", "run", "view", run_id, "--json", "jobs,headSha,url,status,conclusion"], text=True))
            metadata = {"run_id": run_id, "source": "gh run view", **{key: payload.get(key, "") for key in ("headSha", "url", "status", "conclusion")}}
            metadata["snapshot_phase"] = "final" if payload.get("status") == "completed" else "provisional"
            metadata["source_run_completed_success"] = payload.get("status") == "completed" and payload.get("conclusion") == "success"
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
    metadata.setdefault("snapshot_phase", "unknown")
    metadata.setdefault("source_run_completed_success", False)
    return metadata


def _performance(cost_rows: list[dict[str, str]]) -> tuple[list[str], dict[str, Any]]:
    grouped: dict[str, dict[str, float]] = {}
    for row in cost_rows:
        group = grouped.setdefault(row.get("method", ""), {})
        for key in ("equilibrium_requests", "fixedrho_requests", "unique_solves", "residual_calls", "jacobian_calls", "newton_iterations", "runner_seconds", "fallback_count", "retry_count"):
            group[key] = group.get(key, 0.0) + _float(row.get(key, 0.0), 0.0)
    errors: list[str] = []
    hybrid, dense = grouped.get("production_hybrid", {}), grouped.get("memoized_dense", {})
    for key in ("equilibrium_requests", "fixedrho_requests", "unique_solves", "residual_calls", "jacobian_calls", "newton_iterations", "runner_seconds"):
        allowance = 1.10 if key == "runner_seconds" else 1.0
        if hybrid.get(key, 0.0) > allowance * dense.get(key, 0.0) + 1e-9:
            errors.append(f"hybrid {key} exceeds memoized dense")
    hybrid_rate = (hybrid.get("fallback_count", 0.0) + hybrid.get("retry_count", 0.0)) / max(hybrid.get("unique_solves", 0.0), 1.0)
    dense_rate = (dense.get("fallback_count", 0.0) + dense.get("retry_count", 0.0)) / max(dense.get("unique_solves", 0.0), 1.0)
    if (dense_rate == 0.0 and hybrid_rate > 0.0) or (dense_rate > 0.0 and hybrid_rate > 1.25 * dense_rate):
        errors.append("hybrid fallback/retry rate exceeds the 25% risk threshold")
    return errors, {"grouped": grouped, "hybrid_fallback_retry_rate": hybrid_rate, "dense_fallback_retry_rate": dense_rate}


def _gate(output_dir: Path, contract_errors: list[str]) -> dict[str, Any]:
    rows = _rows(output_dir / "slice_metrics.csv")
    costs = _rows(output_dir / "method_costs.csv")
    oracle = [row for row in rows if row.get("method") == "independent_oracle"]
    hybrid = [row for row in rows if row.get("method") == "production_hybrid"]
    oracle_errors: list[str] = []
    for row in oracle:
        if _float(row.get("solver_failure_count"), 0.0) != 0:
            oracle_errors.append(f"oracle solver failure at xi={row.get('xi')} T={row.get('T_MeV')}")
        if row.get("result_status") not in {"confirmed_first_order", "confirmed_monotone", "ambiguous_near_critical"}:
            oracle_errors.append(f"oracle has invalid state at xi={row.get('xi')} T={row.get('T_MeV')}")
        if row.get("result_status") != "ambiguous_near_critical" and str(row.get("geometry_converged", "")).lower() != "true":
            oracle_errors.append(f"oracle geometry not converged at xi={row.get('xi')} T={row.get('T_MeV')}")
    for xi in sorted(XIS):
        xi_rows = [row for row in oracle if math.isclose(_float(row.get("xi")), xi, abs_tol=1e-9, rel_tol=0.0)]
        statuses = {row.get("result_status") for row in xi_rows}
        if "confirmed_first_order" not in statuses or "confirmed_monotone" not in statuses:
            oracle_errors.append(f"oracle xi={xi} lacks two-sided evidence")

    by_oracle_anchor = {
        (round(_float(row.get("xi")), 9), round(_float(row.get("T_MeV")), 9)): row
        for row in oracle
    }
    for xi, temperatures in REQUIRED_FIRST_ORDER.items():
        for temperature in temperatures:
            row = by_oracle_anchor.get((round(xi, 9), round(temperature, 9)))
            if row is None or row.get("result_status") != "confirmed_first_order":
                oracle_errors.append(f"oracle required first-order anchor is not confirmed xi={xi} T={temperature}")
    for xi, temperatures in REQUIRED_MONOTONE.items():
        for temperature in temperatures:
            row = by_oracle_anchor.get((round(xi, 9), round(temperature, 9)))
            if row is None or row.get("result_status") != "confirmed_monotone":
                oracle_errors.append(f"oracle required monotone anchor is not confirmed xi={xi} T={temperature}")
    for xi, (low_temperature, high_temperature) in CEP_NEIGHBOURS.items():
        low = by_oracle_anchor.get((round(xi, 9), round(low_temperature, 9)))
        high = by_oracle_anchor.get((round(xi, 9), round(high_temperature, 9)))
        if low is not None and low.get("result_status") == "confirmed_monotone":
            oracle_errors.append(f"oracle CEP low anchor cannot be confirmed monotone xi={xi} T={low_temperature}")
        if high is not None and high.get("result_status") == "confirmed_first_order":
            oracle_errors.append(f"oracle CEP high anchor cannot be confirmed first-order xi={xi} T={high_temperature}")

    by_anchor = {(row.get("xi"), row.get("T_MeV")): row for row in oracle}
    classification_errors: list[str] = []
    for row in hybrid:
        oracle_row = by_anchor.get((row.get("xi"), row.get("T_MeV")))
        if oracle_row is None:
            classification_errors.append(f"missing oracle anchor for hybrid xi={row.get('xi')} T={row.get('T_MeV')}")
            continue
        oracle_status = oracle_row.get("result_status")
        hybrid_status = row.get("result_status")
        if oracle_status == "ambiguous_near_critical" and hybrid_status == "confirmed_first_order":
            classification_errors.append(f"hybrid confirmed an oracle-ambiguous anchor xi={row.get('xi')} T={row.get('T_MeV')}")
        elif oracle_status != "ambiguous_near_critical" and hybrid_status != oracle_status:
            classification_errors.append(f"hybrid/oracle classification mismatch xi={row.get('xi')} T={row.get('T_MeV')}")
        if hybrid_status == "confirmed_first_order" and oracle_status == "confirmed_first_order":
            for field, tol in (("mu_transition_MeV", 0.025), ("mu_spinodal_hadron_MeV", 0.025), ("mu_spinodal_quark_MeV", 0.025), ("rho_hadron", 0.0025), ("rho_quark", 0.0025), ("rho_spinodal_hadron", 0.0025), ("rho_spinodal_quark", 0.0025)):
                left, right = _float(row.get(field)), _float(oracle_row.get(field))
                if not (math.isfinite(left) and math.isfinite(right)):
                    classification_errors.append(f"missing finite {field} for first-order anchor xi={row.get('xi')} T={row.get('T_MeV')}")
                elif abs(left - right) > tol:
                    classification_errors.append(f"{field} exceeds {tol} at xi={row.get('xi')} T={row.get('T_MeV')}")
            area = max(_float(row.get("area_residual")), _float(oracle_row.get("area_residual")))
            if not math.isfinite(area):
                classification_errors.append(f"missing finite Maxwell area at xi={row.get('xi')} T={row.get('T_MeV')}")
            elif area > 5e-5:
                classification_errors.append(f"Maxwell area exceeds 5e-5 at xi={row.get('xi')} T={row.get('T_MeV')}")
        if row.get("stage_c_status") not in (None, "", "not_run"):
            low, high = _float(row.get("support_low")), _float(row.get("support_high"))
            if not (math.isfinite(low) and math.isfinite(high) and low < high):
                classification_errors.append(f"Stage-C support is absent at xi={row.get('xi')} T={row.get('T_MeV')}")
    coverage_errors: list[str] = []
    for xi in sorted(XIS):
        statuses = {row.get("result_status") for row in hybrid if math.isclose(_float(row.get("xi")), xi, abs_tol=1e-9, rel_tol=0.0)}
        if "confirmed_first_order" not in statuses or "confirmed_monotone" not in statuses:
            coverage_errors.append(f"hybrid xi={xi} lacks confirmed first-order/monotone evidence")
    performance_errors, performance = _performance(costs)
    if contract_errors:
        verdict = "workflow_failure"
    elif oracle_errors:
        verdict = "oracle_inconclusive"
    elif classification_errors:
        verdict = "hybrid_integration_failed"
    elif coverage_errors:
        verdict = "hybrid_integration_failed"
    elif performance_errors:
        verdict = "hybrid_performance_risk"
    else:
        verdict = "full_hybrid_candidate"
    return {"verdict": verdict, "workflow_contract_errors": contract_errors, "oracle_errors": oracle_errors, "classification_errors": classification_errors, "coverage_errors": coverage_errors, "performance_errors": performance_errors, "performance": performance, "automatic_gate_is_not_promotion": True}


def _write_docs(output_dir: Path, gate: dict[str, Any], actions: dict[str, Any]) -> None:
    errors = gate["workflow_contract_errors"] + gate["oracle_errors"] + gate["classification_errors"] + gate["coverage_errors"] + gate["performance_errors"]
    (output_dir / "README.md").write_text(
        f"# PNJL CEP hybrid production shadow v2\n\nverdict: `{gate['verdict']}`。这是显式 opt-in 的诊断产物，"
        "不覆盖 reference，不启动 C0/C1/C2、C3/O1、formal production 或 transport。\n\n"
        f"- Actions critical path: {actions.get('critical_path_seconds', 0.0)} s\n"
        f"- runner-minutes: {actions.get('runner_minutes', 0)}\n"
        f"- errors: {'；'.join(errors) if errors else '无'}\n\n"
        "Stage A 为 rho-support cascade，Stage B 为 memoized dense，Stage C 仅在声明的有限 support 内使用局部 oracle 分辨率；"
        "完整 rho–mu 原始曲线保留在 Actions/local artifact。\n",
        encoding="utf-8",
    )
    _write_csv(output_dir / "claim_ledger.csv", [
        {"claim_id": "oracle", "claim": "independent oracle anchors are stable", "status": "pass" if not gate["oracle_errors"] else "oracle_inconclusive", "boundary": "diagnostic only"},
        {"claim_id": "hybrid", "claim": "hybrid classifications and geometry agree with oracle", "status": "pass" if not gate["classification_errors"] else "hybrid_integration_failed", "boundary": "author physical review required"},
        {"claim_id": "performance", "claim": "hybrid solver work is no higher than memoized dense", "status": "pass" if not gate["performance_errors"] else "hybrid_performance_risk", "boundary": "runner noise allowance applies"},
    ])
    (output_dir / "AUDIT.md").write_text(
        "# Hybrid shadow audit\n\n"
        "本目录只记录 CEP cascade→dense→local-oracle hybrid 的 shadow 证据。"
        "原始曲线留在 Actions/local artifact，仓库表格通过 job SHA 和 curve index 追溯；"
        "任何 verdict 都不自动晋升 reference。\n",
        encoding="utf-8",
    )
    (output_dir / "plot_manifest.json").write_text(json.dumps({"schema_version": "cep_cascade_production_shadow_v2", "figures": [], "reason": "plots are generated by the versioned hybrid plotter"}, indent=2) + "\n", encoding="utf-8")


def collect(input_dir: Path, output_dir: Path, run_id: str | None, expected_sha: str, *, postprocess_sha: str | None = None, source_run_id: str | None = None) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    jobs = _find_jobs(input_dir)
    contract_errors = _validate_jobs(jobs, expected_sha)
    corrections = _collect_tables(jobs, output_dir)
    actions = _actions(run_id, output_dir)
    gate = _gate(output_dir, contract_errors)
    _write_docs(output_dir, gate, actions)
    hashes: dict[str, str] = {}
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and path.name != "manifest.json":
            hashes[str(path.relative_to(output_dir)).replace(os.sep, "/")] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = {
        "schema_version": "cep_cascade_production_shadow_v2",
        "run_id": run_id or "", "source_run_id": source_run_id or run_id or "",
        "expected_calculation_sha": expected_sha, "postprocess_sha": postprocess_sha or "",
        "job_keys": sorted([[summary.get("xi"), summary.get("method")] for _, summary in jobs]),
        "evidence_state": "final" if actions.get("snapshot_phase") == "final" else "provisional",
        "aggregation_corrections": corrections,
        "gate": gate, "actions": actions, "file_sha256": hashes,
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
    gate = collect(args.input_dir, args.output_dir, args.run_id, args.expected_calculation_sha, postprocess_sha=args.postprocess_sha, source_run_id=args.source_run_id)
    print(json.dumps(gate, indent=2, sort_keys=True))
    return 1 if gate["verdict"] == "workflow_failure" else 2 if gate["verdict"] != "full_hybrid_candidate" else 0


if __name__ == "__main__":
    raise SystemExit(main())
