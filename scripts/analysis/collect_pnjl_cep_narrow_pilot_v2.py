#!/usr/bin/env python3
"""Collect and gate PNJL CEP narrow-window pilot v2 evidence."""

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

METHODS = {"c2_dense_baseline", "rho_support_cascade", "high_resolution_oracle"}
XIS = {-0.5, 0.0, 0.5}
DISCOVERY_KEYS = {(xi, "rho_support_cascade", "discovery") for xi in XIS}
VALIDATION_KEYS = {
    (xi, method, "validation") for xi in XIS for method in {"c2_dense_baseline", "high_resolution_oracle"}
}


def _read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _rows(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _float(value: Any) -> float:
    return float(value)


def _finite(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _bool(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


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


def _find_job_dirs(input_dir: Path) -> list[Path]:
    return sorted({path.parent for path in input_dir.rglob("job_summary.json")})


def _load_jobs(input_dir: Path) -> list[dict[str, Any]]:
    jobs = []
    for job_dir in _find_job_dirs(input_dir):
        summary = _read_json(job_dir / "job_summary.json")
        xi = float(summary["xi"])
        method = str(summary["method"])
        stage = str(summary.get("stage", "validation"))
        if method not in METHODS:
            raise ValueError(f"unexpected method {method} in {job_dir}")
        jobs.append({"dir": job_dir, "summary": summary, "xi": xi, "method": method, "stage": stage})
    return jobs


def _validate_jobs(
    jobs: list[dict[str, Any]], expected_calculation_sha: str | None = None
) -> tuple[list[str], list[str]]:
    contract_errors: list[str] = []
    diagnostic_errors: list[str] = []
    seen: set[tuple[float, str, str]] = set()
    for job in jobs:
        key = (job["xi"], job["method"], job["stage"])
        if key in seen:
            contract_errors.append(f"duplicate job key: {key}")
        seen.add(key)
        summary = job["summary"]
        calculation_sha = str(summary.get("calculation_sha", ""))
        if not calculation_sha:
            contract_errors.append(f"missing calculation SHA: {key}")
        elif expected_calculation_sha and calculation_sha != expected_calculation_sha:
            contract_errors.append(
                f"calculation SHA mismatch: {key} expected {expected_calculation_sha}, got {calculation_sha}"
            )
        if str(summary.get("schema_version", "")) != "cep_narrow_pilot_v2":
            contract_errors.append(f"unexpected job schema: {key}")
        for table in ("curve_points.csv", "slice_metrics.csv", "method_costs.csv", "cep_accuracy.csv"):
            if not (job["dir"] / table).is_file():
                contract_errors.append(f"missing {table}: {key}")
        if not bool(summary.get("finite_and_converged_final", False)):
            diagnostic_errors.append(f"non-finite or non-converged final curve: {key}")
        if job["method"] == "rho_support_cascade":
            slice_path = job["dir"] / "slice_metrics.csv"
            slice_rows = _rows(slice_path) if slice_path.is_file() else []
            for row in slice_rows:
                try:
                    targeted = int(float(row.get("targeted_additions", 0)))
                except (TypeError, ValueError):
                    contract_errors.append(f"invalid targeted-point count: {key}")
                    continue
                if targeted > 12:
                    contract_errors.append(f"targeted-point cap exceeded for slice: {key}")
    seen_discovery = seen & DISCOVERY_KEYS
    seen_validation = seen & VALIDATION_KEYS
    if seen_discovery != DISCOVERY_KEYS:
        contract_errors.append(f"missing discovery jobs: {sorted(DISCOVERY_KEYS - seen_discovery)}")
    if seen_validation != VALIDATION_KEYS:
        contract_errors.append(f"missing validation jobs: {sorted(VALIDATION_KEYS - seen_validation)}")
    return contract_errors, diagnostic_errors


def _collect_tables(jobs: list[dict[str, Any]], output_dir: Path) -> None:
    tables = {"curve_points.csv": [], "slice_metrics.csv": [], "method_costs.csv": [], "cep_accuracy.csv": []}
    for job in jobs:
        for name in tables:
            path = job["dir"] / name
            if not path.is_file():
                continue
            for row in _rows(path):
                row["xi"] = job["xi"]
                row["method"] = job["method"]
                row["stage"] = job["stage"]
                tables[name].append(row)
    for name, rows in tables.items():
        _write_csv(output_dir / name, rows)


def _actions_costs(run_id: str | None, output_dir: Path) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    metadata: dict[str, Any] = {"run_id": run_id or "", "source": "unavailable"}
    if run_id:
        try:
            raw = subprocess.check_output(
                ["gh", "run", "view", run_id, "--json", "jobs,headSha,url,status,conclusion"],
                text=True,
            )
            payload = json.loads(raw)
            metadata = {
                "run_id": run_id,
                "source": "gh run view",
                "head_sha": payload.get("headSha", ""),
                "url": payload.get("url", ""),
                "status": payload.get("status", ""),
                "conclusion": payload.get("conclusion", ""),
            }
            for job in payload.get("jobs", []):
                name = str(job.get("name", ""))
                started = job.get("startedAt")
                completed = job.get("completedAt")
                if not started or not completed:
                    continue
                start = datetime.fromisoformat(started.replace("Z", "+00:00"))
                end = datetime.fromisoformat(completed.replace("Z", "+00:00"))
                seconds = max(0.0, (end - start).total_seconds())
                rows.append(
                    {
                        "job_name": name,
                        "job_id": job.get("databaseId", job.get("id", "")),
                        "status": job.get("status", ""),
                        "conclusion": job.get("conclusion", ""),
                        "elapsed_seconds": seconds,
                        "runner_minutes_rounded": math.ceil(seconds / 60.0),
                    }
                )
        except (OSError, subprocess.CalledProcessError, json.JSONDecodeError, ValueError) as exc:
            metadata["error"] = str(exc)
    _write_csv(output_dir / "actions_costs.csv", rows)
    elapsed = [float(row["elapsed_seconds"]) for row in rows]
    metadata["job_count"] = len(rows)
    metadata["critical_path_seconds"] = max(elapsed, default=0.0)
    metadata["raw_total_seconds"] = sum(elapsed)
    metadata["runner_minutes"] = sum(int(row["runner_minutes_rounded"]) for row in rows)
    return metadata


def _accuracy_rows(output_dir: Path) -> list[dict[str, Any]]:
    return _rows(output_dir / "cep_accuracy.csv")


def _cost_rows(output_dir: Path) -> list[dict[str, Any]]:
    return _rows(output_dir / "method_costs.csv")


def _by_xi(rows: list[dict[str, Any]], method: str, stage: str = "validation") -> dict[float, dict[str, Any]]:
    return {
        float(row["xi"]): row
        for row in rows
        if row.get("method") == method and row.get("stage") == stage
    }


def _oracle_gate(accuracy: list[dict[str, Any]], slices: list[dict[str, Any]]) -> tuple[list[str], str]:
    errors: list[str] = []
    oracle = _by_xi(accuracy, "high_resolution_oracle")
    if set(oracle) != XIS:
        errors.append("oracle accuracy rows do not cover all xi")
        return errors, "oracle_inconclusive"
    for xi, row in sorted(oracle.items()):
        if not (_finite(row.get("fine_T_last_first_order_MeV")) and _finite(row.get("fine_T_first_monotone_MeV"))):
            errors.append(f"oracle lacks double-end evidence at xi={xi}")
            continue
        if _finite(row.get("coarse_T_last_first_order_MeV")) and _finite(row.get("coarse_T_first_monotone_MeV")):
            if abs(float(row["coarse_T_last_first_order_MeV"]) - float(row["fine_T_last_first_order_MeV"])) > 0.25:
                errors.append(f"oracle low endpoint is not rho-resolution stable at xi={xi}")
            if abs(float(row["coarse_T_first_monotone_MeV"]) - float(row["fine_T_first_monotone_MeV"])) > 0.25:
                errors.append(f"oracle high endpoint is not rho-resolution stable at xi={xi}")
            if _finite(row.get("coarse_muq_last_first_order_MeV")) and _finite(row.get("fine_muq_last_first_order_MeV")):
                if abs(float(row["coarse_muq_last_first_order_MeV"]) - float(row["fine_muq_last_first_order_MeV"])) > 0.5:
                    errors.append(f"oracle low Maxwell mu is not rho-resolution stable at xi={xi}")
        xi_rows = [item for item in slices if item.get("method") == "high_resolution_oracle" and item.get("stage") == "validation" and abs(float(item["xi"]) - xi) <= 1e-9]
        for item in xi_rows:
            coarse = item.get("coarse_status", "")
            fine = item.get("fine_status", "")
            final = item.get("result_status", "")
            if final == "confirmed_first_order" and not (coarse == "valid" and fine == "valid"):
                errors.append(f"oracle first-order state is not stable across rho levels at xi={xi}, T={item.get('T_MeV')}")
            if final == "confirmed_monotone" and not (coarse == "invalid" and fine == "invalid" and item.get("coarse_reason") == "no_s_shape" and item.get("fine_reason") == "no_s_shape"):
                errors.append(f"oracle monotone state is not stable across rho levels at xi={xi}, T={item.get('T_MeV')}")
    return errors, "stable" if not errors else "oracle_inconclusive"


def _cascade_gate(accuracy: list[dict[str, Any]], slices: list[dict[str, Any]]) -> tuple[list[str], str]:
    errors: list[str] = []
    cascade = _by_xi(accuracy, "rho_support_cascade", "discovery")
    oracle = _by_xi(accuracy, "high_resolution_oracle")
    if set(cascade) != XIS:
        return ["cascade discovery rows do not cover all xi"], "inconclusive"
    for xi in sorted(XIS):
        c = cascade[xi]
        o = oracle.get(xi)
        if o is None:
            errors.append(f"missing oracle row at xi={xi}")
            continue
        c_low = _float(c.get("T_last_first_order_MeV"))
        c_high = _float(c.get("T_first_monotone_MeV"))
        o_low = _float(o.get("fine_T_last_first_order_MeV"))
        o_high = _float(o.get("fine_T_first_monotone_MeV"))
        if not all(_finite(value) for value in (c_low, c_high, o_low, o_high)):
            errors.append(f"cascade/oracle interval is incomplete at xi={xi}")
        else:
            if max(c_low, o_low) > min(c_high, o_high):
                errors.append(f"cascade and oracle ambiguity intervals do not intersect at xi={xi}")
            if c_low > o_low + 0.125 + 1e-9:
                errors.append(f"cascade excludes the oracle low endpoint by more than one resolution step at xi={xi}")
            if c_high < o_high - 0.125 - 1e-9:
                errors.append(f"cascade excludes the oracle high endpoint by more than one resolution step at xi={xi}")
        common = {}
        for row in slices:
            if row.get("stage") != "validation" or row.get("method") not in {"rho_support_cascade", "high_resolution_oracle"}:
                continue
            if abs(float(row["xi"]) - xi) > 1e-9:
                continue
            key = float(row["T_MeV"])
            common.setdefault(key, {})[row["method"]] = row.get("result_status", "")
        comparisons = [pair for pair in common.values() if len(pair) == 2 and all(value != "ambiguous_near_critical" for value in pair.values())]
        for pair in comparisons:
            if pair["rho_support_cascade"] != pair["high_resolution_oracle"]:
                errors.append(f"cascade/oracle classification mismatch at xi={xi}")
                break
    return errors, "within_oracle" if not errors else "inconclusive"


def _performance_gate(costs: list[dict[str, Any]]) -> tuple[list[str], dict[str, Any]]:
    errors: list[str] = []
    grouped: dict[str, dict[str, float]] = {}
    for row in costs:
        method = row.get("method", "")
        group = grouped.setdefault(method, {})
        for key in ("equilibrium_requests", "unique_solves", "residual_calls", "jacobian_calls", "newton_iterations", "runner_seconds", "fallback_count", "retry_count"):
            group[key] = group.get(key, 0.0) + float(row.get(key, 0) or 0)
    cascade = grouped.get("rho_support_cascade", {})
    dense = grouped.get("c2_dense_baseline", {})
    for key in ("equilibrium_requests", "unique_solves", "residual_calls", "jacobian_calls", "newton_iterations", "runner_seconds"):
        if cascade.get(key, 0.0) > dense.get(key, 0.0) + 1e-9:
            errors.append(f"cascade {key} exceeds dense baseline")
    cascade_rate = (cascade.get("fallback_count", 0.0) + cascade.get("retry_count", 0.0)) / max(cascade.get("unique_solves", 0.0), 1.0)
    dense_rate = (dense.get("fallback_count", 0.0) + dense.get("retry_count", 0.0)) / max(dense.get("unique_solves", 0.0), 1.0)
    risk = (dense_rate == 0.0 and cascade_rate > 0.0) or (dense_rate > 0.0 and cascade_rate > 1.25 * dense_rate)
    return errors, {"grouped": grouped, "cascade_fallback_retry_rate": cascade_rate, "dense_fallback_retry_rate": dense_rate, "fallback_retry_risk": risk}


def _gate(
    jobs: list[dict[str, Any]], output_dir: Path, expected_calculation_sha: str | None = None
) -> dict[str, Any]:
    contract_errors, diagnostic_errors = _validate_jobs(jobs, expected_calculation_sha)
    slices = _rows(output_dir / "slice_metrics.csv")
    accuracy = _accuracy_rows(output_dir)
    costs = _cost_rows(output_dir)
    oracle_errors, oracle_status = _oracle_gate(accuracy, slices)
    cascade_errors, cascade_status = _cascade_gate(accuracy, slices)
    performance_errors, performance = _performance_gate(costs)
    all_errors = contract_errors + diagnostic_errors + oracle_errors + cascade_errors + performance_errors
    status = "workflow_failure" if contract_errors else ("pass" if not all_errors else "diagnostic_only")
    if contract_errors:
        overall = "workflow_failure"
    elif oracle_status != "stable":
        overall = "oracle_inconclusive"
    else:
        overall = status
    return {
        "status": overall,
        "workflow_contract_errors": contract_errors,
        "diagnostic_errors": diagnostic_errors,
        "oracle_errors": oracle_errors,
        "cascade_errors": cascade_errors,
        "performance_errors": performance_errors,
        "oracle_status": oracle_status,
        "cascade_status": cascade_status,
        "performance": performance,
        "automatic_gate_is_not_promotion": True,
    }


def _write_manifest(
    output_dir: Path,
    jobs: list[dict[str, Any]],
    gate: dict[str, Any],
    actions: dict[str, Any],
    run_id: str | None,
    expected_calculation_sha: str | None = None,
) -> None:
    hashes: dict[str, str] = {}
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and path.name not in {"manifest.json"}:
            hashes[str(path.relative_to(output_dir)).replace(os.sep, "/")] = hashlib.sha256(path.read_bytes()).hexdigest()
    shas = sorted({str(job["summary"].get("calculation_sha", "")) for job in jobs})
    manifest = {
        "schema_version": "cep_narrow_pilot_v2",
        "generated_at_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "run_id": run_id or "",
        "expected_calculation_sha": expected_calculation_sha or "",
        "job_keys": sorted([[job["xi"], job["method"], job["stage"]] for job in jobs]),
        "calculation_shas": shas,
        "workflow_head_sha": actions.get("head_sha", ""),
        "gate": gate,
        "actions": actions,
        "file_sha256": hashes,
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_claim_ledger(output_dir: Path, gate: dict[str, Any]) -> None:
    rows = [
        {"claim_id": "oracle_interval", "claim": "oracle 在两层 rho 分辨率下形成稳定双端证据", "status": gate["oracle_status"], "evidence": "cep_accuracy.csv, slice_metrics.csv", "boundary": "诊断 gate，不自动晋升"},
        {"claim_id": "cascade_interval", "claim": "cascade ambiguity interval 与 oracle 相交且端点偏差受控", "status": gate["cascade_status"], "evidence": "cep_accuracy.csv, slice_metrics.csv", "boundary": "不代表 production 正确性"},
        {"claim_id": "solver_cost", "claim": "solver 工作量由 request-scoped telemetry 和 exact memoization 记录", "status": "observed", "evidence": "method_costs.csv, actions_costs.csv", "boundary": "分类层微秒耗时不是 blocker"},
        {"claim_id": "performance", "claim": "cascade 聚合 solver 工作量不高于 memoized dense baseline", "status": "pass" if not gate["performance_errors"] else "performance_risk", "evidence": "method_costs.csv, actions_costs.csv", "boundary": "需要作者审核"},
    ]
    _write_csv(output_dir / "claim_ledger.csv", rows)


def _write_readme(output_dir: Path, gate: dict[str, Any], actions: dict[str, Any]) -> None:
    errors = "；".join(
        gate["workflow_contract_errors"] + gate["diagnostic_errors"] + gate["oracle_errors"] + gate["cascade_errors"] + gate["performance_errors"]
    ) or "无"
    text = f"""# PNJL CEP 窄窗口 pilot v2

这是三态 CEP 合同的 diagnostic candidate，不覆盖 `data/reference/pnjl`，不晋升
reference，也不启动 phase-reference 或 transport。完整数值只在 GitHub Actions
运行，旧 canonical 只作为历史叠图。

## 自动 gate（不替代作者物理判断）

- 总 verdict：`{gate['status']}`
- oracle：`{gate['oracle_status']}`
- cascade：`{gate['cascade_status']}`
- Actions critical path：{actions.get('critical_path_seconds', 0.0)} s
- Actions runner-minutes：{actions.get('runner_minutes', 0)}
- 诊断信息：{errors}

`slice_metrics.csv` 分开保留 coarse/fine S-shape、Maxwell、geometry 与三态状态；
`curve_points.csv` 保留每层 rho–mu 原始点及 solver 状态；`method_costs.csv` 记录
equilibrium、residual/Jacobian、Newton、fallback/retry、cache 和 runner 工作量。
`ambiguity_width_T_MeV` 是物理 ambiguity interval，不要求小于
`temperature_resolution_target_MeV=0.125 MeV`。
"""
    (output_dir / "README.md").write_text(text, encoding="utf-8")


def collect(
    input_dir: Path,
    output_dir: Path,
    run_id: str | None = None,
    expected_calculation_sha: str | None = None,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    jobs = _load_jobs(input_dir)
    _collect_tables(jobs, output_dir)
    actions = _actions_costs(run_id, output_dir)
    gate = _gate(jobs, output_dir, expected_calculation_sha)
    _write_claim_ledger(output_dir, gate)
    _write_readme(output_dir, gate, actions)
    _write_manifest(output_dir, jobs, gate, actions, run_id, expected_calculation_sha)
    return gate


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--expected-calculation-sha", default=None)
    args = parser.parse_args()
    gate = collect(args.input_dir, args.output_dir, args.run_id, args.expected_calculation_sha)
    print(json.dumps(gate, indent=2, sort_keys=True))
    return 1 if gate["status"] == "workflow_failure" else 2 if gate["status"] in {"diagnostic_only", "oracle_inconclusive"} else 0


if __name__ == "__main__":
    raise SystemExit(main())
