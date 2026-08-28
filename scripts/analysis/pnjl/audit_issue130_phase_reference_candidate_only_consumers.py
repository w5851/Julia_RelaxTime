#!/usr/bin/env python3
"""Audit candidate-only phase-reference consumer routes without a solver.

The audit exercises the Python candidate adapter and the explicit plot loader,
then checks that Julia/transport entrypoints expose the same opt-in
``candidate_only`` route and the explicit ``legacy`` rollback.  It does not run
the PNJL equilibrium solver, transport scan, or write either reference tree.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping


SCHEMA_VERSION = "pnjl_issue130_phase_reference_candidate_only_consumer_audit_v1"
CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
CANDIDATE_RELATIVE = Path("data/reference/pnjl/issue130_phase_reference_v1")
LEGACY_RELATIVE = Path("data/reference/pnjl/legacy_phase_reference_v1")
AUDIT_RELATIVE = Path(
    "docs/analysis/pnjl/phase_reference/issue130_phase_reference_legacy_audit_v1"
)
OUTPUT_RELATIVE = Path(
    "docs/analysis/pnjl/phase_reference/issue130_phase_reference_candidate_only_consumer_audit_v1"
)
CANONICAL_CONSUMER_XIS = tuple(round(-0.5 + 0.05 * index, 12) for index in range(21))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def write_csv(path: Path, fieldnames: Iterable[str], rows: Iterable[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        writer.writerows(rows)


def import_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def read_text(repo_root: Path, relative_path: str) -> str:
    return (repo_root / relative_path).read_text(encoding="utf-8-sig", errors="replace")


def _stage_a_decision(repo_root: Path) -> dict[str, Any]:
    path = repo_root / AUDIT_RELATIVE / "decision.json"
    if not path.is_file():
        raise FileNotFoundError(f"Stage A decision is required: {path}")
    decision = read_json(path)
    if decision.get("calculation_sha") != CALCULATION_SHA:
        raise ValueError("Stage A calculation SHA does not match the fixed candidate SHA")
    return decision


def build_audit(repo_root: Path, output_root: Path) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"refusing to overwrite non-empty audit output: {output_root}")
    output_root.mkdir(parents=True, exist_ok=True)

    candidate_root = repo_root / CANDIDATE_RELATIVE
    legacy_root = repo_root / LEGACY_RELATIVE
    candidate_manifest_path = candidate_root / "manifest.json"
    strict_manifest_path = candidate_root / "strict" / "manifest.json"
    legacy_manifest_path = legacy_root / "RETIREMENT_MANIFEST.json"
    candidate_manifest = read_json(candidate_manifest_path)
    strict_manifest = read_json(strict_manifest_path)
    legacy_manifest = read_json(legacy_manifest_path)
    stage_a = _stage_a_decision(repo_root)

    adapter = import_module(
        repo_root / "scripts/pnjl/phase_reference_adapter.py",
        "issue130_phase_reference_adapter_candidate_only_audit",
    )
    candidate = adapter.load_phase_reference(candidate_root, layer="strict", allow_runtime=False)
    candidate_only = adapter.build_candidate_only_view(candidate)

    all_candidate_rows_certified = all(
        bool(row.get("certified", False))
        for rows in candidate_only.tables.values()
        for row in rows
    )
    candidate_only_fallback_disabled = (
        candidate_only.diagnostics.get("runtime_view") == "certified_candidate_only"
        and candidate_only.diagnostics.get("fallback_enabled") is False
        and not any(candidate_only.diagnostics.get("fallback_row_counts", {}).values())
    )

    plot_loader = import_module(
        repo_root / "scripts/pnjl/plot_phase_diagram.py",
        "issue130_phase_diagram_candidate_only_audit",
    )
    plot_data = plot_loader.load_candidate_phase_data(candidate_root, "strict")
    plot_counts = {
        "boundary": len(plot_data[0]),
        "cep": len(plot_data[1]),
        "spinodals": len(plot_data[2]),
        "crossover": len(plot_data[3]),
    }

    # The transport/phase-guided consumers request a canonical xi list and
    # then use the adapter's nearest-xi behavior.  For a retirement decision we
    # must record whether the candidate has an exact certified xi first; a
    # nearest neighbor is not evidence that a legacy fallback is unreachable.
    dynamic_request_rows: list[dict[str, Any]] = []
    dynamic_incomplete = False
    for consumer in ("transport_phase_guided", "phase_guided_plan"):
        for table in ("boundary", "crossover", "cep", "spinodals"):
            table_rows = candidate_only.tables.get(table, ())
            for xi in CANONICAL_CONSUMER_XIS:
                exact = [
                    row
                    for row in table_rows
                    if abs(float(row.get("xi", float("nan"))) - xi) <= 1e-10
                ]
                complete = bool(exact)
                dynamic_incomplete = dynamic_incomplete or not complete
                dynamic_request_rows.append(
                    {
                        "consumer": consumer,
                        "table": table,
                        "requested_xi": xi,
                        "candidate_exact": complete,
                        "candidate_certified_rows": len(exact),
                        "fallback_needed": not complete,
                        "status": "pass" if complete else "fallback_required",
                    }
                )

    run_gap = read_text(repo_root, "scripts/relaxtime/run_gap_transport_scan.jl")
    gap_cli = read_text(repo_root, "scripts/relaxtime/gap_transport_scan_cli.jl")
    phase_cli = read_text(repo_root, "scripts/relaxtime/phase_guided_transport_scan_cli.jl")
    phase_adapter = read_text(repo_root, "scripts/relaxtime/phase_reference_adapter.jl")
    validator = read_text(repo_root, "scripts/pnjl/validate_phase_data.py")
    rollback_declared = (
        "mode === :legacy" in run_gap
        and "--phase-reference-mode" in gap_cli
        and "source === :legacy" in phase_adapter
        and legacy_manifest_path.is_file()
    )
    candidate_only_declared = (
        "mode === :candidate_only" in run_gap
        and ":candidate_only" in gap_cli
        and ":candidate_only" in phase_cli
        and "load_phase_reference_candidate_only" in phase_adapter
    )

    consumer_rows = [
        {
            "consumer": "python_phase_reference_adapter",
            "route": "candidate_only",
            "status": "pass" if candidate_only_fallback_disabled and all_candidate_rows_certified else "fail",
            "solver_called": False,
            "fallback_enabled": candidate_only.diagnostics.get("fallback_enabled"),
            "detail": "certified candidate rows only; unresolved rows filtered",
        },
        {
            "consumer": "julia_phase_reference_adapter",
            "route": "candidate_only",
            "status": "pass" if candidate_only_declared else "fail",
            "solver_called": False,
            "fallback_enabled": False,
            "detail": "explicit loader and diagnostics route declared",
        },
        {
            "consumer": "phase_guided_transport",
            "route": "candidate_only",
            "status": "pass" if candidate_only_declared else "fail",
            "solver_called": False,
            "fallback_enabled": False,
            "detail": "CLI/parser route is opt-in; production default remains fallback",
        },
        {
            "consumer": "pnjl_plot_phase_diagram",
            "route": "explicit_candidate_root",
            "status": "pass" if all(plot_counts.values()) else "fail",
            "solver_called": False,
            "fallback_enabled": False,
            "detail": json.dumps(plot_counts, sort_keys=True),
        },
        {
            "consumer": "legacy_rollback",
            "route": "explicit_legacy",
            "status": "pass" if rollback_declared else "fail",
            "solver_called": False,
            "fallback_enabled": False,
            "detail": "snapshot and explicit legacy mode remain declared",
        },
        {
            "consumer": "pnjl_validate_phase_data",
            "route": "candidate_only",
            "status": "pass" if "def validate_candidate_reference" in validator else "migration_required",
            "solver_called": False,
            "fallback_enabled": False,
            "detail": "explicit candidate schema validation route; legacy CSV validation remains backwards compatible",
        },
    ]

    hard_failures = [row["consumer"] for row in consumer_rows if row["status"] == "fail"]
    stop_reasons = [
        "legacy_fallback_key_count_nonzero",
        "default_runtime_fallback_preserved",
    ]
    if dynamic_incomplete:
        stop_reasons.append("dynamic_request_key_coverage_incomplete")
    if not any(
        row["consumer"] == "pnjl_validate_phase_data" and row["status"] == "pass"
        for row in consumer_rows
    ):
        stop_reasons.append("validator_candidate_schema_migration_pending")

    decision = {
        "schema_version": SCHEMA_VERSION,
        "verdict": "candidate_only_contract_supported" if not hard_failures else "candidate_only_contract_inconclusive",
        "candidate_only_contract_supported": not hard_failures,
        "legacy_rollback_contract_supported": rollback_declared,
        "dynamic_request_key_coverage_complete": not dynamic_incomplete,
        "candidate_only_consumer_failures": hard_failures,
        "candidate_uncertified_rows": int(stage_a.get("candidate_uncertified_rows", 0)),
        "legacy_fallback_key_count": int(stage_a.get("legacy_fallback_key_count", 0)),
        "physical_deletion_eligible": False,
        "solver_called": False,
        "reference_write": False,
        "runtime_consumption": False,
        "calculation_sha": CALCULATION_SHA,
        "source_run_id": stage_a.get("source_run_id", ""),
        "replay_run_id": stage_a.get("replay_run_id", ""),
        "stop_reasons": stop_reasons,
    }
    candidate_contract_ok = candidate_only_fallback_disabled and all_candidate_rows_certified
    claims = [
        {
            "claim_id": "candidate_only_runtime_view",
            "claim": "Candidate-only view exposes only certified candidate rows and disables legacy fallback.",
            "status": "supported" if candidate_contract_ok else "not_supported",
            "evidence": "consumer_matrix.csv; candidate manifest hashes",
            "boundary": "does not prove all dynamic consumer request keys are covered",
        },
        {
            "claim_id": "explicit_legacy_rollback",
            "claim": "The explicit legacy rollback route and immutable snapshot remain declared.",
            "status": "supported" if rollback_declared else "not_supported",
            "evidence": "consumer_matrix.csv; legacy RETIREMENT_MANIFEST.json",
            "boundary": "rollback preservation is not evidence that legacy can be deleted",
        },
        {
            "claim_id": "physical_deletion_ready",
            "claim": "PNJL legacy snapshot is ready for physical deletion.",
            "status": "not_supported",
            "evidence": "decision.json; Stage A fallback matrix",
            "boundary": "requires fallback=0 and a separate author-authorized allowlist PR",
        },
        {
            "claim_id": "dynamic_request_key_coverage",
            "claim": "All canonical transport and phase-guided request xi keys are covered by certified candidate rows.",
            "status": "supported" if not dynamic_incomplete else "not_supported",
            "evidence": "dynamic_request_matrix.csv",
            "boundary": "an exact xi gap still requires the explicit legacy fallback; nearest-xi behavior is not coverage evidence",
        },
    ]

    write_csv(
        output_root / "consumer_matrix.csv",
        ["consumer", "route", "status", "solver_called", "fallback_enabled", "detail"],
        consumer_rows,
    )
    write_csv(
        output_root / "dynamic_request_matrix.csv",
        [
            "consumer",
            "table",
            "requested_xi",
            "candidate_exact",
            "candidate_certified_rows",
            "fallback_needed",
            "status",
        ],
        dynamic_request_rows,
    )
    write_json(output_root / "decision.json", decision)
    write_json(output_root / "claim_ledger.json", claims)
    write_csv(output_root / "claim_ledger.csv", list(claims[0]), claims)

    readme = f"""# Issue #130 candidate-only consumer audit v1

这是阶段 B 的 solver-free consumer contract 审计；不运行 PNJL/equilibrium/transport，
不写 candidate 或 legacy reference。verdict：`{decision['verdict']}`。

- candidate calculation SHA：`{CALCULATION_SHA}`
- candidate-only rows：`{sum(len(rows) for rows in candidate_only.tables.values())}`（全部 certified）
- candidate unresolved rows（Stage A）：`{decision['candidate_uncertified_rows']}`
- legacy fallback keys（Stage A）：`{decision['legacy_fallback_key_count']}`
- canonical request-key coverage：`{'complete' if decision['dynamic_request_key_coverage_complete'] else 'incomplete'}`
- remaining fallback-required requests：`{sum(1 for row in dynamic_request_rows if row['fallback_needed'])}`
- physical deletion：`False`

`consumer_matrix.csv` 区分 candidate-only、显式 legacy rollback 和 validator 的
candidate-schema 入口；`dynamic_request_matrix.csv` 按 21 个 canonical xi、四张
phase 表和两个实际消费者记录 exact certified coverage。当前默认 `runtime`
fallback 不变；只有后续请求键覆盖和消费者迁移完成后才可评估物理清理。
"""
    (output_root / "README.md").write_text(readme, encoding="utf-8", newline="\n")
    audit = (
        "# Candidate-only consumer audit\n\n"
        "The candidate-only route is an explicit opt-in migration path. "
        "It does not promote unresolved rows and does not authorize legacy deletion.\n"
    )
    (output_root / "AUDIT.md").write_text(audit, encoding="utf-8", newline="\n")

    generated = sorted(path for path in output_root.rglob("*") if path.is_file())
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "repo_head": "" if not (repo_root / ".git").exists() else _git_head(repo_root),
        "calculation_sha": CALCULATION_SHA,
        "candidate_root": str(CANDIDATE_RELATIVE).replace("\\", "/"),
        "legacy_root": str(LEGACY_RELATIVE).replace("\\", "/"),
        "candidate_manifest_sha256": sha256(candidate_manifest_path),
        "candidate_strict_manifest_sha256": sha256(strict_manifest_path),
        "legacy_retirement_manifest_sha256": sha256(legacy_manifest_path),
        "candidate_reference_status": candidate_manifest.get("reference_status"),
        "candidate_strict_reference_write": strict_manifest.get("reference_write"),
        "legacy_snapshot_status": legacy_manifest.get("status"),
        "solver_called": False,
        "reference_write": False,
        "runtime_consumption": False,
        "decision": "decision.json",
        "verdict": decision["verdict"],
        "files": [
            {"path": path.relative_to(output_root).as_posix(), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in generated
        ],
    }
    write_json(output_root / "manifest.json", manifest)
    return manifest


def _git_head(repo_root: Path) -> str:
    result = subprocess.run(
        ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--output-root", type=Path, default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    output_root = (args.output_root or repo_root / OUTPUT_RELATIVE).resolve()
    manifest = build_audit(repo_root, output_root)
    print(json.dumps({"output_root": str(output_root), "verdict": manifest["verdict"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
