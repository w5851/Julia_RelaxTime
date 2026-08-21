#!/usr/bin/env python3
"""Audit Issue #130 phase-reference runtime consumer compatibility.

This is a solver-free, static compatibility audit.  It inventories the
candidate tables, the repository's phase/transport entrypoints, and their
reference-path contracts.  It deliberately does not create a legacy-shaped
adapter, alter defaults, import the candidate at runtime, or call PNJL.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


SCHEMA = "pnjl_issue130_phase_reference_runtime_consumer_audit_v1"
CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
SOURCE_RUN_ID = "32354095831"
REPLAY_RUN_ID = "32451053476"
CANDIDATE_RELATIVE = Path("data/reference/pnjl/issue130_phase_reference_v1")
LEGACY_FILES = (
    Path("data/reference/pnjl/boundary.csv"),
    Path("data/reference/pnjl/cep.csv"),
    Path("data/reference/pnjl/spinodals.csv"),
    Path("data/reference/pnjl/crossover_dense.csv"),
    Path("data/reference/pnjl/phase_reference_dense_manifest.json"),
)

CONSUMERS = (
    {
        "consumer_id": "gap_transport_defaults",
        "entrypoint": "scripts/relaxtime/run_gap_transport_scan.jl",
        "kind": "transport_phase_hint",
        "contract": "boundary_dense.csv|boundary.csv; cep_dense.csv|cep.csv; crossover_dense.csv|crossover.csv",
        "candidate_reachable": "false",
        "compatibility": "requires_adapter",
        "reason": "candidate sibling is not in preferred_phase_reference_path and its table names/schema differ",
    },
    {
        "consumer_id": "phase_guided_transport",
        "entrypoint": "scripts/relaxtime/phase_guided_transport_scan_plan.jl",
        "kind": "transport_phase_guided_plan",
        "contract": "delegates to GapTransportScanPhaseEquilibrium defaults",
        "candidate_reachable": "false",
        "compatibility": "requires_adapter",
        "reason": "plan uses load_phase_boundary_data/load_crossover_reference backed by legacy defaults",
    },
    {
        "consumer_id": "paper_p1_reference_root",
        "entrypoint": "scripts/relaxtime/run_paper_p1_pipeline.jl",
        "kind": "explicit_phase_reference_root",
        "contract": "boundary_<tag>.csv; spinodals_<tag>.csv; crossover_<tag>.csv; cep_<tag>.csv",
        "candidate_reachable": "explicit_path_only",
        "compatibility": "requires_adapter",
        "reason": "CLI accepts an explicit root/tag but candidate has strict/derived/render table layout, not tagged legacy files",
    },
    {
        "consumer_id": "legacy_phase_plot",
        "entrypoint": "scripts/pnjl/plot_phase_diagram.py",
        "kind": "legacy_plot_validation",
        "contract": "data/reference/pnjl/{boundary.csv,cep.csv,spinodals.csv}",
        "candidate_reachable": "false",
        "compatibility": "requires_adapter",
        "reason": "legacy plotting CLI has no candidate-layer selector or schema mapper",
    },
    {
        "consumer_id": "phase_artifacts_api",
        "entrypoint": "src/models/phase/PhaseArtifacts.jl",
        "kind": "generic_phase_artifact_promotion",
        "contract": "data/reference/<model>/phase_diagram",
        "candidate_reachable": "not_applicable",
        "compatibility": "not_a_consumer",
        "reason": "generic processed phase artifacts are a separate contract from the PNJL candidate sibling",
    },
)

TABLE_SPECS = {
    "strict/tables/maxwell_surface_strict_reference_v1.csv": {
        "surface": "maxwell",
        "keys": ("xi", "T_MeV"),
        "numeric": ("xi", "T_MeV", "mu_MeV", "rho_hadron", "rho_quark", "area_residual"),
        "legacy": "boundary.csv",
    },
    "strict/tables/crossover_surface_strict_reference_v1.csv": {
        "surface": "crossover",
        "keys": ("xi", "mu_MeV"),
        "numeric": ("xi", "mu_MeV", "T_MeV", "rho", "mu_CEP_proxy_MeV"),
        "legacy": "crossover_dense.csv",
    },
    "strict/tables/spinodal_surface_strict_reference_v1.csv": {
        "surface": "spinodal",
        "keys": ("xi", "T_MeV"),
        "numeric": ("xi", "T_MeV", "mu_spinodal_hadron_MeV", "mu_spinodal_quark_MeV"),
        "legacy": "spinodals.csv",
    },
    "strict/tables/cep_boundary_strict_reference_v1.csv": {
        "surface": "cep_boundary",
        "keys": ("xi",),
        "numeric": ("xi", "mu_CEP_proxy_MeV", "T_low_MeV", "T_high_MeV", "T_midpoint_MeV"),
        "legacy": "cep.csv",
    },
    "derived/tables/maxwell_surface_derived_reference_v1.csv": {
        "surface": "maxwell",
        "keys": ("xi", "T_MeV"),
        "numeric": ("xi", "T_MeV", "mu_MeV", "rho_hadron", "rho_quark", "area_residual"),
        "legacy": "boundary.csv",
    },
    "derived/tables/crossover_surface_derived_reference_v1.csv": {
        "surface": "crossover",
        "keys": ("xi", "mu_MeV"),
        "numeric": ("xi", "mu_MeV", "T_MeV", "rho", "mu_CEP_proxy_MeV"),
        "legacy": "crossover_dense.csv",
    },
    "derived/tables/spinodal_surface_derived_reference_v1.csv": {
        "surface": "spinodal",
        "keys": ("xi", "T_MeV"),
        "numeric": ("xi", "T_MeV", "mu_spinodal_hadron_MeV", "mu_spinodal_quark_MeV"),
        "legacy": "spinodals.csv",
    },
    "derived/tables/cep_boundary_derived_reference_v1.csv": {
        "surface": "cep_boundary",
        "keys": ("xi",),
        "numeric": ("xi", "mu_CEP_proxy_MeV", "T_low_MeV", "T_high_MeV", "T_midpoint_MeV"),
        "legacy": "cep.csv",
    },
    "render/tables/maxwell_surface_render.csv": {
        "surface": "maxwell",
        "keys": ("xi", "T_MeV"),
        "numeric": ("xi", "T_MeV", "mu_MeV", "rho_hadron", "rho_quark", "area_residual"),
        "legacy": "boundary.csv",
    },
    "render/tables/crossover_surface_render.csv": {
        "surface": "crossover",
        "keys": ("xi", "mu_MeV"),
        "numeric": ("xi", "mu_MeV", "T_MeV", "rho", "mu_CEP_proxy_MeV"),
        "legacy": "crossover_dense.csv",
    },
    "render/tables/cep_boundary_render.csv": {
        "surface": "cep_boundary",
        "keys": ("xi",),
        "numeric": ("xi", "mu_CEP_proxy_MeV", "T_low_MeV", "T_high_MeV", "T_midpoint_MeV"),
        "legacy": "cep.csv",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--reference-root", type=Path, default=None)
    parser.add_argument("--output-root", type=Path, default=None)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"missing CSV header: {path}")
        return list(reader.fieldnames), list(reader)


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def finite(value: str | None) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def duplicate_count(rows: Iterable[dict[str, str]], keys: tuple[str, ...]) -> int:
    values = [tuple(row.get(key, "") for key in keys) for row in rows]
    return len(values) - len(set(values))


def git_value(repo_root: Path, *args: str) -> str:
    try:
        return subprocess.check_output(["git", *args], cwd=repo_root, text=True, stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        return ""


def rel(path: Path, repo_root: Path) -> str:
    return path.resolve().relative_to(repo_root.resolve()).as_posix()


def write_csv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def candidate_inventory(reference_root: Path, repo_root: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    inventory: list[dict[str, Any]] = []
    checks: list[dict[str, Any]] = []
    for relative, spec in TABLE_SPECS.items():
        path = reference_root / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        header, rows = read_csv(path)
        missing = [field for field in spec["numeric"] if field not in header]
        nonfinite = sum(not finite(row.get(field)) for row in rows for field in spec["numeric"])
        duplicates = duplicate_count(rows, spec["keys"])
        inventory.append({
            "path": rel(path, repo_root),
            "surface": spec["surface"],
            "layer": relative.split("/", 1)[0],
            "rows": len(rows),
            "columns": len(header),
            "sha256": sha256(path),
            "keys": ",".join(spec["keys"]),
            "missing_numeric_fields": ",".join(missing),
            "nonfinite_numeric_values": nonfinite,
            "duplicate_keys": duplicates,
            "legacy_contract": spec["legacy"],
        })
        checks.append({
            "table": rel(path, repo_root),
            "status": "pass" if not missing and nonfinite == 0 and duplicates == 0 else "fail",
            "details": f"rows={len(rows)}; missing={','.join(missing)}; nonfinite={nonfinite}; duplicates={duplicates}",
        })
    return inventory, checks


def legacy_inventory(repo_root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for relative in LEGACY_FILES:
        path = repo_root / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        rows.append({
            "path": relative.as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
            "header": path.read_text(encoding="utf-8").splitlines()[0] if path.suffix == ".csv" else "json",
        })
    return rows


def source_contract_checks(repo_root: Path) -> list[dict[str, Any]]:
    checks: list[dict[str, Any]] = []
    source_paths = [
        repo_root / "scripts/relaxtime/run_gap_transport_scan.jl",
        repo_root / "scripts/relaxtime/gap_transport_scan_phase_equilibrium.jl",
        repo_root / "scripts/relaxtime/phase_guided_transport_scan_plan.jl",
        repo_root / "scripts/relaxtime/run_paper_p1_pipeline.jl",
        repo_root / "scripts/pnjl/plot_phase_diagram.py",
    ]
    for path in source_paths:
        text = path.read_text(encoding="utf-8")
        candidate_mentions = bool(re.search(r"issue130_phase_reference_v1", text))
        checks.append({
            "entrypoint": rel(path, repo_root),
            "candidate_literal_present": candidate_mentions,
            "status": "pass" if not candidate_mentions else "review",
            "details": "no implicit candidate path literal" if not candidate_mentions else "candidate path is explicitly referenced",
        })
    return checks


def main() -> int:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    reference_root = (args.reference_root or repo_root / CANDIDATE_RELATIVE).resolve()
    output_root = (args.output_root or repo_root / "docs/analysis/pnjl/phase_reference/issue130_phase_reference_runtime_consumer_audit_v1").resolve()
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"refusing to overwrite non-empty audit output: {output_root}")
    output_root.mkdir(parents=True, exist_ok=True)
    tables_root = output_root / "tables"

    root_manifest = read_json(reference_root / "manifest.json")
    if root_manifest.get("calculation_sha") != CALCULATION_SHA:
        raise ValueError("candidate calculation SHA mismatch")
    if root_manifest.get("runtime_consumption") is not False:
        raise ValueError("candidate must remain runtime_consumption=false")

    inventory, table_checks = candidate_inventory(reference_root, repo_root)
    legacy_rows = legacy_inventory(repo_root)
    contract_checks = source_contract_checks(repo_root)
    write_csv(tables_root / "candidate_table_inventory.csv", list(inventory[0].keys()), inventory)
    write_csv(tables_root / "candidate_table_checks.csv", ["table", "status", "details"], table_checks)
    write_csv(tables_root / "legacy_reference_snapshot.csv", list(legacy_rows[0].keys()), legacy_rows)
    write_csv(tables_root / "consumer_inventory.csv", list(CONSUMERS[0].keys()), CONSUMERS)
    write_csv(tables_root / "source_contract_checks.csv", list(contract_checks[0].keys()), contract_checks)

    schema_rows = []
    for relative, spec in TABLE_SPECS.items():
        path = reference_root / relative
        header, _ = read_csv(path)
        legacy_required = {
            "boundary.csv": ("xi", "T_MeV", "mu_transition_MeV", "rho_hadron", "rho_quark"),
            "cep.csv": ("xi", "T_CEP_MeV", "muq_CEP_MeV", "muB_CEP_MeV"),
            "spinodals.csv": ("xi", "T_MeV", "mu_spinodal_hadron_MeV", "mu_spinodal_quark_MeV", "rho_spinodal_hadron", "rho_spinodal_quark"),
            "crossover_dense.csv": ("xi", "mu_MeV", "T_crossover_MeV", "rho", "method", "converged", "derivative", "variable"),
        }[spec["legacy"]]
        missing_legacy = sorted(set(legacy_required) - set(header))
        schema_rows.append({
            "candidate_table": rel(path, repo_root),
            "legacy_contract": spec["legacy"],
            "candidate_columns": len(header),
            "missing_legacy_fields": ",".join(missing_legacy),
            "direct_compatible": not missing_legacy,
            "verdict": "direct_compatible" if not missing_legacy else "requires_adapter",
        })
    write_csv(tables_root / "schema_compatibility.csv", list(schema_rows[0].keys()), schema_rows)

    candidate_isolated = (
        root_manifest.get("runtime_consumption") is False
        and all(not item["candidate_literal_present"] for item in contract_checks)
        and all(item["compatibility"] != "direct_compatible" for item in CONSUMERS)
    )
    adapter_required = any(row["verdict"] == "requires_adapter" for row in schema_rows)
    verdict = "candidate_isolation_confirmed_requires_explicit_adapter" if candidate_isolated and adapter_required else "audit_inconclusive"
    decision = {
        "schema_version": SCHEMA,
        "verdict": verdict,
        "candidate_auto_selected": False,
        "runtime_consumption": False,
        "solver_called": False,
        "legacy_defaults_unchanged": True,
        "explicit_adapter_required": adapter_required,
        "calculation_sha": CALCULATION_SHA,
        "source_run_id": SOURCE_RUN_ID,
        "replay_run_id": REPLAY_RUN_ID,
        "non_goals": [
            "do not rewrite legacy reference files",
            "do not add a runtime adapter in this audit",
            "do not run PNJL or transport production",
        ],
        "next_action": "design a separate explicit candidate adapter/consumer contract and author-review it before any runtime switch",
    }
    (output_root / "decision.json").write_text(json.dumps(decision, indent=2) + "\n", encoding="utf-8")

    claims = [
        {
            "claim_id": "candidate_integrity",
            "claim": "Candidate strict/derived/render tables are finite and key-unique under the imported package contract.",
            "status": "supported" if all(row["status"] == "pass" for row in table_checks) else "not_supported",
            "evidence": "tables/candidate_table_checks.csv; data/reference/pnjl/issue130_phase_reference_v1/manifest.json",
            "boundary": "integrity does not authorize runtime use",
        },
        {
            "claim_id": "runtime_isolation",
            "claim": "Existing phase/transport entrypoints do not implicitly select the Issue #130 candidate sibling.",
            "status": "supported" if candidate_isolated else "inconclusive",
            "evidence": "tables/consumer_inventory.csv; tables/source_contract_checks.csv",
            "boundary": "static source audit; no runtime execution",
        },
        {
            "claim_id": "schema_adapter_required",
            "claim": "The candidate cannot be passed directly to legacy consumers without a schema/semantic adapter.",
            "status": "supported" if adapter_required else "inconclusive",
            "evidence": "tables/schema_compatibility.csv",
            "boundary": "adapter design is a separate task and must preserve unresolved/non-certified semantics",
        },
        {
            "claim_id": "promotion",
            "claim": "This audit promotes or switches the runtime phase reference.",
            "status": "not_supported",
            "evidence": "decision.json; data/reference/pnjl/issue130_phase_reference_v1/manifest.json",
            "boundary": "runtime_consumption=false; RS transport remains blocked",
        },
    ]
    write_csv(output_root / "claim_ledger.csv", list(claims[0].keys()), claims)

    script_path = Path(__file__).resolve()
    manifest = {
        "schema_version": SCHEMA,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "repo_head": git_value(repo_root, "rev-parse", "HEAD"),
        "calculation_sha": CALCULATION_SHA,
        "source_run_id": SOURCE_RUN_ID,
        "replay_run_id": REPLAY_RUN_ID,
        "candidate_root": rel(reference_root, repo_root),
        "candidate_manifest_sha256": sha256(reference_root / "manifest.json"),
        "script": rel(script_path, repo_root),
        "script_sha256": sha256(script_path),
        "solver_called": False,
        "runtime_consumption": False,
        "candidate_auto_selected": False,
        "decision": "decision.json",
        "tables": sorted(rel(path, repo_root) for path in tables_root.glob("*.csv")),
        "input_files": sorted(rel(path, repo_root) for path in reference_root.rglob("*") if path.is_file()),
    }
    (output_root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    (output_root / "README.md").write_text(
        "# Issue #130 phase-reference runtime consumer audit v1\n\n"
        "这是一个 solver-free 静态兼容性审计。它核验 versioned candidate 的表结构、"
        "旧 reference 默认路径和 phase/transport/paper 入口，确认 candidate 不会被隐式消费。\n\n"
        f"当前 verdict：`{verdict}`。candidate 的 `runtime_consumption=false` 保持不变；"
        "candidate 表与旧 consumer schema 不同，任何 runtime 切换都必须另行设计并审核显式 adapter。\n\n"
        "本包不修改旧 reference、不运行 PNJL、不运行 transport，也不把 strict unresolved 或 derived "
        "non-certified 行静默升级为 production certificate。\n",
        encoding="utf-8",
    )
    (output_root / "AUDIT.md").write_text(
        "# Issue #130 runtime consumer compatibility audit\n\n"
        "## 结论\n\n"
        f"- verdict: `{verdict}`\n"
        "- 现有 consumer 仍通过旧 `data/reference/pnjl` 文件名/路径解析；未发现 candidate sibling 的隐式路径。\n"
        "- strict/derived/render 表不满足旧 boundary/CEP/spinodal/crossover schema，不能直接替换。\n"
        "- 下一步只能是独立的显式 adapter/consumer contract 设计；它必须保留 strict unresolved、derived `interpolated_noncertified`、CEP bracket 和 first-order/crossover 互斥语义。\n\n"
        "## 边界\n\n"
        "本审计为静态 source/schema 检查，不执行 Julia include，不调用 equilibrium solver，不运行 transport。"
        "因此它确认的是路径隔离与字段契约，不是 transport 数值正确性，也不是 reference promotion。\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
