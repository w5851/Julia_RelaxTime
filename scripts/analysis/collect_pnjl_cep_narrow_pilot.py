#!/usr/bin/env python3
"""Collect and gate PNJL CEP narrow pilot job artifacts.

The collector is intentionally standard-library only so it can run in the
workflow summary job and on a downloaded artifact directory.  It never edits
`data/reference/pnjl`; evidence is written to the caller-selected directory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

METHODS = {
    "c2_dense_baseline",
    "rho_support_cascade",
    "high_resolution_oracle",
}
XIS = {-0.5, 0.0, 0.5}
REQUIRED_JOB_KEYS = {(xi, method) for xi in XIS for method in METHODS}


def _read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _rows(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _float(value: Any) -> float:
    return float(value)


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
    candidates = []
    for path in sorted(input_dir.rglob("job_summary.json")):
        candidates.append(path.parent)
    if not candidates and (input_dir / "job_summary.json").is_file():
        candidates.append(input_dir)
    return candidates


def _load_jobs(input_dir: Path) -> list[dict[str, Any]]:
    jobs = []
    for job_dir in _find_job_dirs(input_dir):
        summary = _read_json(job_dir / "job_summary.json")
        xi = float(summary["xi"])
        method = str(summary["method"])
        if method not in METHODS:
            raise ValueError(f"unexpected method {method} in {job_dir}")
        jobs.append({"dir": job_dir, "summary": summary, "xi": xi, "method": method})
    return jobs


def _validate_jobs(jobs: list[dict[str, Any]]) -> list[str]:
    errors: list[str] = []
    seen: set[tuple[float, str]] = set()
    for job in jobs:
        key = (job["xi"], job["method"])
        if key in seen:
            errors.append(f"duplicate job key: {key}")
        seen.add(key)
        summary = job["summary"]
        if not str(summary.get("calculation_sha", "")):
            errors.append(f"missing calculation SHA: {key}")
        if not bool(summary.get("finite_and_converged", False)):
            errors.append(f"non-finite or non-converged solve: {key}")
        for table in ("slice_metrics.csv", "method_costs.csv", "cep_accuracy.csv"):
            if not (job["dir"] / table).is_file():
                errors.append(f"missing {table}: {key}")
        slice_path = job["dir"] / "slice_metrics.csv"
        if slice_path.is_file():
            for row in _rows(slice_path):
                try:
                    targeted = int(float(row.get("targeted_additions", 0)))
                except (TypeError, ValueError):
                    errors.append(f"invalid targeted-point count: {key}")
                    continue
                if targeted > 12:
                    errors.append(f"targeted-point cap exceeded for slice: {key}")
    missing = REQUIRED_JOB_KEYS - seen
    if missing:
        errors.append(f"missing matrix jobs: {sorted(missing)}")
    return errors


def _collect_tables(jobs: list[dict[str, Any]], output_dir: Path) -> None:
    slice_rows: list[dict[str, Any]] = []
    cost_rows: list[dict[str, Any]] = []
    accuracy_rows: list[dict[str, Any]] = []
    for job in jobs:
        for row in _rows(job["dir"] / "slice_metrics.csv"):
            row["xi"] = job["xi"]
            row["method"] = job["method"]
            slice_rows.append(row)
        for row in _rows(job["dir"] / "method_costs.csv"):
            row["xi"] = job["xi"]
            row["method"] = job["method"]
            cost_rows.append(row)
        for row in _rows(job["dir"] / "cep_accuracy.csv"):
            row["xi"] = job["xi"]
            row["method"] = job["method"]
            accuracy_rows.append(row)
    _write_csv(output_dir / "slice_metrics.csv", slice_rows)
    _write_csv(output_dir / "method_costs.csv", cost_rows)
    _write_csv(output_dir / "cep_accuracy.csv", accuracy_rows)


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
                if "pilot" not in name.lower() and "cep" not in name.lower():
                    continue
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


def _gate(jobs: list[dict[str, Any]], accuracy_rows: list[dict[str, Any]]) -> dict[str, Any]:
    errors = _validate_jobs(jobs)
    oracle = [row for row in accuracy_rows if row.get("method") == "high_resolution_oracle"]
    oracle_refine_present = all(str(row.get("oracle_refine", "")).strip() for row in oracle)
    oracle_refine_stable = oracle_refine_present and all(_bool(row.get("oracle_refine_stable", False)) for row in oracle)
    if not oracle_refine_present:
        errors.append("oracle refinement evidence is missing")
    if not oracle_refine_stable:
        errors.append("oracle 0.00625 -> 0.003125 classification is not stable")
    canonical_match = True
    for row in oracle:
        try:
            canonical_match &= abs(_float(row["delta_T_MeV"])) <= 0.5
            canonical_match &= abs(_float(row["delta_muq_MeV"])) <= 2.0
        except (KeyError, ValueError):
            canonical_match = False
    if not canonical_match:
        errors.append("oracle does not satisfy the predeclared canonical CEP tolerance")
    cascade = [row for row in accuracy_rows if row.get("method") == "rho_support_cascade"]
    baseline = [row for row in accuracy_rows if row.get("method") == "c2_dense_baseline"]
    cascade_in_oracle = len(cascade) == len(oracle) and all(
        abs(_float(c["estimated_T_CEP_MeV"]) - _float(o["estimated_T_CEP_MeV"])) <= max(_float(o["bracket_width_MeV"]), 0.125)
        for c, o in zip(sorted(cascade, key=lambda r: float(r["xi"])), sorted(oracle, key=lambda r: float(r["xi"])))
    )
    if not cascade_in_oracle:
        errors.append("cascade CEP is outside the oracle bracket/uncertainty")
    return {
        "status": "pass" if not errors else "diagnostic_only",
        "errors": errors,
        "oracle_status": "stable_and_canonical" if oracle_refine_present and oracle_refine_stable and canonical_match else "oracle_inconclusive",
        "cascade_status": "within_oracle" if cascade_in_oracle else "outside_oracle",
        "automatic_gate_is_not_promotion": True,
    }


def _write_manifest(output_dir: Path, jobs: list[dict[str, Any]], gate: dict[str, Any], actions: dict[str, Any]) -> None:
    hashes: dict[str, str] = {}
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and path.name not in {"manifest.json"}:
            hashes[str(path.relative_to(output_dir)).replace(os.sep, "/")] = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = {
        "schema_version": "cep_narrow_pilot_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "job_keys": sorted([[job["xi"], job["method"]] for job in jobs]),
        "calculation_shas": sorted({job["summary"].get("calculation_sha", "") for job in jobs}),
        "gate": gate,
        "actions": actions,
        "file_sha256": hashes,
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_claim_ledger(output_dir: Path, gate: dict[str, Any]) -> None:
    rows = [
        {
            "claim_id": "oracle_stability",
            "claim": "oracle 0.00625 -> 0.003125 的 CEP bracket/分类稳定",
            "status": gate["oracle_status"],
            "evidence": "cep_accuracy.csv, slice_metrics.csv",
            "boundary": "仅为 diagnostic gate，不自动晋升 reference",
        },
        {
            "claim_id": "oracle_canonical_agreement",
            "claim": "oracle 与旧 canonical 的 |ΔT|<=0.5 MeV 且 |Δmuq|<=2 MeV",
            "status": gate["oracle_status"],
            "evidence": "cep_accuracy.csv",
            "boundary": "若失败则 verdict=oracle_inconclusive",
        },
        {
            "claim_id": "cascade_accuracy",
            "claim": "cascade CEP 位于 oracle bracket/不确定度内",
            "status": gate["cascade_status"],
            "evidence": "cep_accuracy.csv",
            "boundary": "不作物理正确性或 production 晋升结论",
        },
        {
            "claim_id": "solver_cost",
            "claim": "solver 工作量由 request-scoped telemetry 与 exact memoization 记录",
            "status": "observed",
            "evidence": "method_costs.csv, actions_costs.csv",
            "boundary": "分类层微秒耗时不是性能 gate",
        },
    ]
    _write_csv(output_dir / "claim_ledger.csv", rows)


def _write_readme(output_dir: Path, gate: dict[str, Any], actions: dict[str, Any]) -> None:
    status = gate["status"]
    errors = "；".join(gate["errors"]) if gate["errors"] else "无"
    text = f"""# PNJL CEP 窄窗口 pilot v1

这是一个候选诊断产物，不覆盖 `data/reference/pnjl`，也不代表 reference 晋升。
每个矩阵 job 固定同一 calculation SHA、`p_num=24`、`t_num=8`、
`rs_reduced_adaptive` 和 ground-state pressure governance；完整数值计算只在
GitHub Actions 执行。

## 自动 gate（作者物理判断仍然必需）

- 总 verdict：`{status}`
- oracle verdict：`{gate["oracle_status"]}`
- cascade verdict：`{gate["cascade_status"]}`
- 诊断信息：`{errors}`
- Actions critical path：{actions.get("critical_path_seconds", 0.0)} s
- Actions runner-minutes（向上取整）：{actions.get("runner_minutes", 0)}

`method_costs.csv` 区分 equilibrium request、unique solve、memo cache hit、
residual/Jacobian、Newton/trust-region、fallback/rescue/retry；`slice_metrics.csv`
保留每个 `(T, xi, method)` 切片的局部计数。`claim_ledger.csv` 明确观察、gate
与不可推出的物理结论边界。
"""
    (output_dir / "README.md").write_text(text, encoding="utf-8")


def collect(input_dir: Path, output_dir: Path, run_id: str | None = None) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    jobs = _load_jobs(input_dir)
    validation_errors = _validate_jobs(jobs)
    if validation_errors:
        # Still write a manifest and tables so failure artifacts remain useful.
        pass
    _collect_tables(jobs, output_dir)
    accuracy_rows = _rows(output_dir / "cep_accuracy.csv")
    actions = _actions_costs(run_id, output_dir)
    gate = _gate(jobs, accuracy_rows)
    _write_claim_ledger(output_dir, gate)
    _write_readme(output_dir, gate, actions)
    _write_manifest(output_dir, jobs, gate, actions)
    return {"gate": gate, "actions": actions, "validation_errors": validation_errors}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--run-id", default=None)
    args = parser.parse_args()
    result = collect(args.input_dir.resolve(), args.output_dir.resolve(), args.run_id)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["gate"]["status"] == "pass" else 2


if __name__ == "__main__":
    raise SystemExit(main())
