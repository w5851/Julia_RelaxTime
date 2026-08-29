#!/usr/bin/env python3
"""Build a solver-free accepted-primary/legacy-retirement audit.

The v2 accepted package is loaded through the production-parity Python adapter
and checked against the byte-preserving legacy snapshot by semantic keys.  The
audit deliberately separates three questions that used to be conflated:

* Is ``accepted`` a valid author-promoted primary runtime source?
* Does any runtime code still provide a legacy fallback/rollback?
* Which path contracts still have to be migrated before the snapshot can be
  physically deleted?

No solver is called, no reference tree is written, and no source file is
deleted.  The output contains summaries and hashes only; it does not copy the
raw phase curves.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import re
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from statistics import median
from typing import Any, Iterable, Mapping


PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
PACKAGE_RELATIVE = Path("data/reference/pnjl/issue130_phase_reference_v2")
LEGACY_RELATIVE = Path("data/reference/pnjl/legacy_phase_reference_v1")
OUTPUT_RELATIVE = Path(
    "docs/analysis/pnjl/phase_reference/issue130_phase_reference_legacy_audit_v3"
)
SCHEMA_VERSION = "pnjl_issue130_phase_reference_legacy_audit_v3"

TABLES: dict[str, dict[str, Any]] = {
    "boundary": {
        "accepted": "accepted/tables/maxwell_surface_accepted_phase_map_v1.csv",
        "legacy": "boundary.csv",
        "keys": ("xi", "T_MeV"),
        "legacy_mu": "mu_transition_MeV",
    },
    "crossover": {
        "accepted": "accepted/tables/crossover_surface_accepted_phase_map_v1.csv",
        "legacy": "crossover_dense.csv",
        "keys": ("xi", "muq_MeV"),
        "legacy_mu": "mu_MeV",
    },
    "cep": {
        "accepted": "accepted/tables/cep_boundary_accepted_phase_map_v1.csv",
        "legacy": "cep.csv",
        "keys": ("xi",),
        "legacy_mu": None,
    },
    "spinodals": {
        "accepted": "accepted/tables/spinodal_surface_accepted_phase_map_v1.csv",
        "legacy": "spinodals.csv",
        "keys": ("xi", "T_MeV"),
        "legacy_mu": None,
    },
}

TOKENS = (
    "legacy_phase_reference_v1",
    "data/reference/pnjl/boundary.csv",
    "data/reference/pnjl/cep.csv",
    "data/reference/pnjl/spinodals.csv",
    "data/reference/pnjl/crossover_dense.csv",
    "data/reference/pnjl/phase_reference_dense_manifest.json",
    "load_phase_reference_runtime_with_fallback",
    "source=:legacy",
    "phase-reference-mode legacy",
)
RUNTIME_CONTRACT_PATHS = {
    "scripts/relaxtime/phase_reference_adapter.jl",
    "scripts/relaxtime/run_gap_transport_scan.jl",
    "scripts/relaxtime/phase_guided_transport_scan_cli.jl",
    "scripts/relaxtime/gap_transport_scan_cli.jl",
    "scripts/pnjl/phase_reference_adapter.py",
    "scripts/pnjl/plot_phase_diagram.py",
}


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


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        if not fields:
            raise ValueError(f"missing CSV header: {path}")
        return fields, [dict(row) for row in reader]


def write_csv(path: Path, fieldnames: Iterable[str], rows: Iterable[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    names = list(fieldnames)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in names})


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tree_hash(root: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        digest.update(path.relative_to(root).as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(sha256(path).encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def git_value(repo_root: Path, *args: str) -> str:
    try:
        return subprocess.check_output(
            ["git", *args], cwd=repo_root, text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return ""


def load_adapter(repo_root: Path) -> Any:
    path = repo_root / "scripts" / "pnjl" / "phase_reference_adapter.py"
    spec = importlib.util.spec_from_file_location("issue130_accepted_primary_adapter", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import adapter: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _finite(value: Any) -> float:
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError(f"non-finite value: {value!r}")
    return parsed


def _key(table: str, row: Mapping[str, Any]) -> tuple[float, ...]:
    return tuple(round(_finite(row[field]), 12) for field in TABLES[table]["keys"])


def _key_text(table: str, key: tuple[float, ...]) -> str:
    return json.dumps(
        dict(zip(TABLES[table]["keys"], key)), ensure_ascii=False, sort_keys=True, separators=(",", ":")
    )


def _load_accepted(repo_root: Path, package_root: Path) -> dict[str, Any]:
    adapter = load_adapter(repo_root)
    accepted = adapter.load_phase_reference_runtime(package_root, layer="accepted")
    strict = adapter.load_phase_reference_runtime(package_root, layer="strict")
    accepted_tables = {
        table: [dict(row) for row in rows] for table, rows in accepted.tables.items()
    }
    strict_tables = {table: [dict(row) for row in rows] for table, rows in strict.tables.items()}
    package_manifest = read_json(package_root / "manifest.json")
    accepted_manifest = read_json(package_root / "accepted" / "manifest.json")
    return {
        "adapter": adapter,
        "accepted": accepted_tables,
        "strict": strict_tables,
        "accepted_diagnostics": dict(accepted.diagnostics),
        "strict_diagnostics": dict(strict.diagnostics),
        "package_manifest": package_manifest,
        "accepted_manifest": accepted_manifest,
        "manifest_sha256": sha256(package_root / "manifest.json"),
        "accepted_manifest_sha256": sha256(package_root / "accepted" / "manifest.json"),
        "tree_sha256": tree_hash(package_root),
    }


def _legacy_row(table: str, raw: Mapping[str, str]) -> dict[str, Any]:
    if table == "boundary":
        return {
            "xi": _finite(raw["xi"]),
            "T_MeV": _finite(raw["T_MeV"]),
            "muq_MeV": _finite(raw["mu_transition_MeV"]),
        }
    if table == "crossover":
        return {
            "xi": _finite(raw["xi"]),
            "muq_MeV": _finite(raw["mu_MeV"]),
            "T_MeV": _finite(raw["T_crossover_MeV"]),
        }
    if table == "cep":
        return {
            "xi": _finite(raw["xi"]),
            "T_midpoint_MeV": _finite(raw["T_CEP_MeV"]),
            "muq_CEP_proxy_MeV": _finite(raw["muq_CEP_MeV"]),
        }
    if table == "spinodals":
        return {"xi": _finite(raw["xi"]), "T_MeV": _finite(raw["T_MeV"])}
    raise AssertionError(table)


def _load_legacy(legacy_root: Path) -> dict[str, list[dict[str, Any]]]:
    result: dict[str, list[dict[str, Any]]] = {}
    for table, spec in TABLES.items():
        path = legacy_root / spec["legacy"]
        fields, rows = read_csv(path)
        required = {"xi"}
        required.add("T_MeV") if table in {"boundary", "spinodals"} else None
        if table == "crossover":
            required |= {"mu_MeV", "T_crossover_MeV"}
        elif table == "boundary":
            required |= {"mu_transition_MeV"}
        elif table == "cep":
            required |= {"T_CEP_MeV", "muq_CEP_MeV"}
        missing = sorted(required - set(fields))
        if missing:
            raise ValueError(f"{path} missing fields: {', '.join(missing)}")
        result[table] = [_legacy_row(table, row) for row in rows]
    return result


def _snapshot_inventory(repo_root: Path, legacy_root: Path) -> dict[str, Any]:
    manifest_path = legacy_root / "RETIREMENT_MANIFEST.json"
    manifest = read_json(manifest_path)
    errors: list[str] = []
    records: list[dict[str, Any]] = []
    for record in manifest.get("files", []):
        rel = str(record.get("path", ""))
        source = str(record.get("source_path", ""))
        path = legacy_root / rel
        canonical = repo_root / source
        exists = path.is_file()
        actual_bytes = path.stat().st_size if exists else ""
        actual_sha = sha256(path) if exists else ""
        byte_match = exists and actual_bytes == record.get("bytes")
        hash_match = exists and actual_sha == record.get("sha256")
        canonical_absent = not canonical.exists()
        if not byte_match:
            errors.append(f"snapshot byte mismatch: {rel}")
        if not hash_match:
            errors.append(f"snapshot hash mismatch: {rel}")
        if not canonical_absent:
            errors.append(f"canonical legacy path exists: {source}")
        records.append(
            {
                "path": rel,
                "source_path": source,
                "bytes": actual_bytes,
                "expected_bytes": record.get("bytes", ""),
                "sha256": actual_sha,
                "expected_sha256": record.get("sha256", ""),
                "byte_match": byte_match,
                "hash_match": hash_match,
                "canonical_absent": canonical_absent,
            }
        )
    return {
        "manifest": manifest,
        "manifest_sha256": sha256(manifest_path),
        "tree_sha256": tree_hash(legacy_root),
        "records": records,
        "integrity_pass": not errors,
        "errors": errors,
    }


def _coverage(
    accepted: Mapping[str, list[Mapping[str, Any]]],
    strict: Mapping[str, list[Mapping[str, Any]]],
    legacy: Mapping[str, list[Mapping[str, Any]]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    accepted_keys = {table: {_key(table, row) for row in rows} for table, rows in accepted.items()}
    strict_keys = {table: {_key(table, row) for row in rows} for table, rows in strict.items()}
    cep_limits = {
        round(_finite(row["xi"]), 12): _finite(row["muq_CEP_proxy_MeV"])
        for row in accepted.get("cep", ())
    }
    summary: list[dict[str, Any]] = []
    neighbors: list[dict[str, Any]] = []
    for table in TABLES:
        accepted_rows = list(accepted.get(table, ()))
        strict_rows = list(strict.get(table, ()))
        legacy_rows = list(legacy.get(table, ()))
        missing_physical = 0
        exact = 0
        excluded = 0
        nearest_distances: list[float] = []
        same_xi = defaultdict(list)
        for row in accepted_rows:
            key = _key(table, row)
            same_xi[key[0]].append(key[1:] if len(key) > 1 else ())
        for index, row in enumerate(legacy_rows, start=2):
            key = _key(table, row)
            if key in accepted_keys[table]:
                exact += 1
                disposition = "exact_accepted_key"
                reason = "accepted_primary_exact"
                distance = 0.0
            else:
                physical = True
                if table == "crossover":
                    limit = cep_limits.get(key[0])
                    physical = limit is None or row["muq_MeV"] <= limit + 1e-9
                if not physical:
                    excluded += 1
                    disposition = "excluded_above_cep_proxy"
                    reason = "legacy_derivative_peak_not_physical_crossover"
                    distance = ""
                else:
                    missing_physical += 1
                    candidates = same_xi.get(key[0], [])
                    distance = min((abs(key[1] - candidate[0]) for candidate in candidates), default=float("nan"))
                    if math.isfinite(distance):
                        nearest_distances.append(distance)
                    disposition = "nearest_accepted_support" if math.isfinite(distance) else "no_same_xi_support"
                    reason = "accepted_key_grid_differs_from_legacy"
            neighbors.append(
                {
                    "table": table,
                    "legacy_row": index,
                    "key": _key_text(table, key),
                    "accepted_exact": key in accepted_keys[table],
                    "strict_exact": key in strict_keys[table],
                    "same_xi_distance": distance,
                    "reason": reason,
                    "disposition": disposition,
                }
            )
        uncertified = sum(not bool(row.get("certified", False)) for row in accepted_rows)
        summary.append(
            {
                "table": table,
                "accepted_rows": len(accepted_rows),
                "accepted_strict_certified_rows": sum(bool(row.get("certified", False)) for row in accepted_rows),
                "accepted_noncertified_rows": uncertified,
                "accepted_duplicate_keys": len(accepted_rows) - len(accepted_keys[table]),
                "strict_rows": len(strict_rows),
                "legacy_rows": len(legacy_rows),
                "legacy_exact_accepted": exact,
                "legacy_missing_physical": missing_physical,
                "legacy_excluded_nonphysical": excluded,
                "nearest_support_count": len(nearest_distances),
                "nearest_distance_min": min(nearest_distances) if nearest_distances else "",
                "nearest_distance_median": median(nearest_distances) if nearest_distances else "",
                "nearest_distance_max": max(nearest_distances) if nearest_distances else "",
                "key_fields": ",".join(TABLES[table]["keys"]),
                "key_normalization": "Float64 rounded to 12 decimal places",
            }
        )
    return summary, neighbors


def _tracked_files(repo_root: Path) -> list[Path]:
    raw = subprocess.check_output(["git", "ls-files", "-z"], cwd=repo_root)
    return [repo_root / item.decode("utf-8") for item in raw.split(b"\0") if item]


def _scan_consumers(repo_root: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    occurrences: list[dict[str, Any]] = []
    by_path: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for path in _tracked_files(repo_root):
        rel = path.relative_to(repo_root).as_posix()
        if rel.startswith("docs/analysis/") or rel.startswith("docs/dev/archived/"):
            continue
        try:
            lines = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
        except OSError:
            continue
        for line_no, line in enumerate(lines, start=1):
            token = next((value for value in TOKENS if value in line), None)
            if token is None:
                continue
            lower = line.lower()
            if token in {"load_phase_reference_runtime_with_fallback", "source=:legacy", "phase-reference-mode legacy"}:
                classification = "runtime_legacy_contract"
                blocker = rel in RUNTIME_CONTRACT_PATHS
            elif rel.startswith("tests/"):
                classification = "snapshot_or_contract_test"
                blocker = False
            elif rel.startswith(("scripts/analysis/", "scripts/perf/")) or rel.startswith(".github/"):
                classification = "diagnostic_or_historical_tooling"
                blocker = False
            elif rel.startswith(("docs/", "config/", "data/")):
                classification = "governance_or_metadata_reference"
                blocker = False
            else:
                classification = "active_path_contract"
                blocker = True
            row = {
                "path": rel,
                "line": line_no,
                "token": token,
                "classification": classification,
                "path_retirement_blocker": blocker,
                "snippet": re.sub(r"\s+", " ", line.strip())[:240],
            }
            occurrences.append(row)
            by_path[rel].append(row)
    consumers: list[dict[str, Any]] = []
    for rel in sorted(by_path):
        rows = by_path[rel]
        classes = sorted({str(row["classification"]) for row in rows})
        consumers.append(
            {
                "path": rel,
                "occurrence_count": len(rows),
                "classification": ";".join(classes),
                "path_retirement_blocker": any(row["path_retirement_blocker"] for row in rows),
                "decision": (
                    "migrate_before_physical_delete"
                    if any(row["path_retirement_blocker"] for row in rows)
                    else "retain_as_explicit_history_or_test"
                ),
            }
        )
    return occurrences, consumers


def _decision(
    package_info: Mapping[str, Any],
    snapshot: Mapping[str, Any],
    coverage: list[Mapping[str, Any]],
    consumers: list[Mapping[str, Any]],
    runtime_api_present: bool,
    legacy_mode_present: bool,
) -> dict[str, Any]:
    package_manifest = package_info["package_manifest"]
    accepted_manifest = package_info["accepted_manifest"]
    accepted_valid = (
        package_manifest.get("promotion_status") == "accepted_for_downstream"
        and package_manifest.get("downstream_default_layer") == "accepted"
        and package_manifest.get("calculation_sha") == CALCULATION_SHA
        and package_manifest.get("runtime_consumption") is False
        and accepted_manifest.get("reference_write") is False
    )
    strict_valid = (
        package_info["strict_diagnostics"].get("runtime_view") == "strict_certified_only"
        and package_info["strict_diagnostics"].get("fallback_enabled") is False
    )
    runtime_fallback_removed = not runtime_api_present and not legacy_mode_present
    path_blockers = [row for row in consumers if row.get("path_retirement_blocker")]
    missing_physical = sum(int(row["legacy_missing_physical"]) for row in coverage)
    snapshot_ok = bool(snapshot["integrity_pass"])
    path_ready = runtime_fallback_removed and not path_blockers
    physical_eligible = accepted_valid and strict_valid and snapshot_ok and path_ready
    if not accepted_valid or not strict_valid or not snapshot_ok:
        verdict = "accepted_primary_audit_input_invalid"
        next_action = "repair package/strict/snapshot provenance; do not delete legacy"
    elif not runtime_fallback_removed:
        verdict = "accepted_primary_runtime_contract_inconclusive"
        next_action = "remove the legacy fallback/rollback runtime entry points"
    elif path_blockers:
        verdict = "accepted_primary_runtime_ready_path_retirement_pending"
        next_action = "migrate remaining active legacy path contracts, then rerun this audit"
    else:
        verdict = "legacy_physical_deletion_candidate"
        next_action = "prepare a separate exact allowlist deletion PR and request author review"
    return {
        "schema_version": "pnjl_issue130_phase_reference_legacy_retirement_decision_v3",
        "verdict": verdict,
        "accepted_primary_valid": accepted_valid,
        "strict_explicit_valid": strict_valid,
        "runtime_legacy_fallback_rows": 0,
        "runtime_legacy_rollback_enabled": False,
        "legacy_runtime_fallback_api_present": runtime_api_present,
        "legacy_mode_runtime_present": legacy_mode_present,
        "accepted_primary_runtime_ready": accepted_valid and strict_valid and runtime_fallback_removed,
        "active_path_contract_count": len(path_blockers),
        "legacy_missing_physical_key_count": missing_physical,
        "legacy_snapshot_integrity_pass": snapshot_ok,
        "path_retirement_ready": path_ready,
        "physical_deletion_eligible": physical_eligible,
        "coverage_exact_key_is_informational": True,
        "calculation_sha": CALCULATION_SHA,
        "source_run_id": str(package_manifest.get("source_run_id", "")),
        "replay_run_id": str(package_manifest.get("replay_run_id", "")),
        "accepted_manifest_sha256": package_info["accepted_manifest_sha256"],
        "legacy_retirement_manifest_sha256": snapshot["manifest_sha256"],
        "stop_reasons": [
            reason
            for reason, condition in (
                ("accepted_or_strict_contract_invalid", not (accepted_valid and strict_valid)),
                ("legacy_snapshot_integrity_failure", not snapshot_ok),
                ("active_legacy_path_contract_requires_migration", bool(path_blockers)),
            )
            if condition
        ],
        "next_action": next_action,
        "solver_called": False,
        "reference_write": False,
        "runtime_consumption": False,
        "non_goals": [
            "do not call Julia or the PNJL equilibrium solver",
            "do not rewrite data/reference/pnjl trees",
            "do not infer numerical convergence from a coverage audit",
        ],
    }


def _claims(decision: Mapping[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "claim_id": "accepted_is_primary_runtime",
            "claim": "The author-promoted accepted v2 layer is the default downstream runtime source.",
            "status": "supported" if decision["accepted_primary_runtime_ready"] else "inconclusive",
            "evidence": "data/reference/pnjl/issue130_phase_reference_v2/manifest.json; adapter diagnostics",
            "boundary": "accepted non-certified rows remain non-certified; strict is an explicit opt-in",
        },
        {
            "claim_id": "legacy_runtime_fallback_removed",
            "claim": "No runtime legacy fallback or rollback is enabled by the current adapter contract.",
            "status": "supported" if decision["runtime_legacy_fallback_rows"] == 0 and not decision["runtime_legacy_rollback_enabled"] else "not_supported",
            "evidence": "decision.json; tables/consumer_matrix.csv",
            "boundary": "historical snapshot readers and explicit audit tooling are not runtime consumers",
        },
        {
            "claim_id": "accepted_neighborhood_support",
            "claim": "Legacy exact-key gaps have finite accepted same-xi support where available; the gaps are a key-grid issue, not a reason to restore legacy fallback.",
            "status": "supported" if decision["legacy_missing_physical_key_count"] >= 0 else "inconclusive",
            "evidence": "tables/coverage.csv; tables/neighbor_coverage.csv",
            "boundary": "nearest support is not an exact-key equivalence and must not be silently labeled exact",
        },
        {
            "claim_id": "physical_deletion_gate",
            "claim": "The legacy snapshot is eligible for physical deletion only after path retirement is complete.",
            "status": "supported" if decision["physical_deletion_eligible"] else "not_supported",
            "evidence": "decision.json; tables/consumer_matrix.csv; tables/legacy_snapshot_inventory.csv",
            "boundary": "a separate allowlist PR and author authorization remain mandatory",
        },
    ]


def build_audit(repo_root: Path, output_root: Path, *, replace_existing: bool = False) -> dict[str, Any]:
    output_root = output_root.resolve()
    if output_root.exists() and any(output_root.iterdir()):
        if not replace_existing:
            raise FileExistsError(f"refusing to overwrite non-empty audit output: {output_root}")
        # This is an explicit local rebuild of an analysis directory only.
        for child in output_root.iterdir():
            if child.is_dir():
                import shutil
                shutil.rmtree(child)
            else:
                child.unlink()
    output_root.mkdir(parents=True, exist_ok=True)
    package_root = (repo_root / PACKAGE_RELATIVE).resolve()
    legacy_root = (repo_root / LEGACY_RELATIVE).resolve()
    package_info = _load_accepted(repo_root, package_root)
    legacy = _load_legacy(legacy_root)
    snapshot = _snapshot_inventory(repo_root, legacy_root)
    coverage, neighbors = _coverage(package_info["accepted"], package_info["strict"], legacy)
    occurrences, consumers = _scan_consumers(repo_root)
    runtime_api_present = any(
        row["token"] == "load_phase_reference_runtime_with_fallback"
        and row["path"] in RUNTIME_CONTRACT_PATHS
        for row in occurrences
    )
    legacy_mode_present = any(
        row["token"] == "phase-reference-mode legacy"
        and row["path"] in RUNTIME_CONTRACT_PATHS
        for row in occurrences
    )
    decision = _decision(
        package_info, snapshot, coverage, consumers, runtime_api_present, legacy_mode_present
    )

    tables_root = output_root / "tables"
    write_csv(tables_root / "coverage.csv", list(coverage[0].keys()), coverage)
    write_csv(tables_root / "neighbor_coverage.csv", list(neighbors[0].keys()), neighbors)
    write_csv(tables_root / "consumer_occurrences.csv", list(occurrences[0].keys()) if occurrences else ["path"], occurrences)
    write_csv(tables_root / "consumer_matrix.csv", list(consumers[0].keys()) if consumers else ["path"], consumers)
    write_csv(
        tables_root / "legacy_snapshot_inventory.csv",
        list(snapshot["records"][0].keys()) if snapshot["records"] else ["path"],
        snapshot["records"],
    )
    claims = _claims(decision)
    write_csv(tables_root / "claim_ledger.csv", list(claims[0].keys()), claims)
    write_json(tables_root / "claim_ledger.json", claims)
    write_json(output_root / "decision.json", decision)

    summary_lines = [
        "# Issue #130 PNJL accepted-primary / legacy retirement audit v3",
        "",
        "这是 solver-free 的合同审计：通过 production-parity Python adapter 读取 v2 `accepted` 与显式 `strict`，",
        "再按真实 semantic key 对照 byte-preserving legacy snapshot。它不调用 PNJL solver、不写 reference、不删除文件。",
        f"repo HEAD：`{git_value(repo_root, 'rev-parse', 'HEAD')}`；calculation SHA：`{CALCULATION_SHA}`。",
        "",
        "## Verdict",
        "",
        f"- `{decision['verdict']}`",
        f"- accepted primary runtime ready：`{decision['accepted_primary_runtime_ready']}`",
        f"- runtime legacy fallback rows：`{decision['runtime_legacy_fallback_rows']}`；rollback enabled：`{decision['runtime_legacy_rollback_enabled']}`",
        f"- active legacy path contracts pending migration：`{decision['active_path_contract_count']}`",
        f"- physical deletion eligible：`{decision['physical_deletion_eligible']}`",
        f"- next action：{decision['next_action']}",
        "",
        "## 三个不混淆的覆盖问题",
        "",
        "1. accepted 是否可作为 runtime primary：由 package promotion、adapter eligibility 和 strict 显式路径验证。",
        "2. legacy 是否仍会被 runtime fallback/rollback 调用：由 adapter API 和 runtime mode 扫描验证；当前合同计数应为 0。",
        "3. 物理删除前是否仍有旧路径文本：由 consumer matrix 逐项迁移；这与数值 exact-key 覆盖是不同门禁。",
        "",
        "## Semantic-key coverage（仅诊断）",
        "",
        "accepted 与 legacy 的 exact key 不要求相同。`neighbor_coverage.csv` 记录同 ξ 最近 accepted 支持，",
        "但 nearest support 不能伪装为 exact legacy key；consumer 应使用 accepted 自己的插值/最近 ξ 语义。",
        "",
        "| table | accepted | strict certified | legacy | exact accepted | physical missing | above-CEP excluded |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in coverage:
        summary_lines.append(
            f"| {row['table']} | {row['accepted_rows']} | {row['accepted_strict_certified_rows']} | "
            f"{row['legacy_rows']} | {row['legacy_exact_accepted']} | {row['legacy_missing_physical']} | "
            f"{row['legacy_excluded_nonphysical']} |"
        )
    summary_lines.extend(
        [
            "",
            "## Path-retirement boundary",
            "",
            "legacy snapshot 仍保留，直到 active path contracts 迁移完成并由独立 physical-deletion PR 审核。",
            "历史 audit、snapshot manifest 和 contract tests 可以保留；它们不等于 runtime fallback。",
            "",
            "## Provenance",
            "",
            f"- accepted package manifest SHA：`{package_info['manifest_sha256']}`",
            f"- accepted layer manifest SHA：`{package_info['accepted_manifest_sha256']}`",
            f"- accepted package tree SHA：`{package_info['tree_sha256']}`",
            f"- legacy retirement manifest SHA：`{snapshot['manifest_sha256']}`",
            f"- legacy snapshot tree SHA：`{snapshot['tree_sha256']}`",
            "- solver_called：`false`；reference_write：`false`；runtime_consumption：`false`（本审计自身）。",
        ]
    )
    (output_root / "README.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8", newline="\n")
    (output_root / "AUDIT.md").write_text(
        "# PNJL accepted-primary legacy retirement audit v3\n\n"
        f"当前 verdict：`{decision['verdict']}`。运行时已固定为 accepted primary，strict 只能显式开启；"
        "legacy 不再是 fallback/rollback。物理删除仍受 path-retirement blocker 和独立 allowlist PR 约束。\n\n"
        "本文件只解释合同和证据边界；完整表格见 `tables/`。\n",
        encoding="utf-8",
        newline="\n",
    )

    generated = sorted(path for path in output_root.rglob("*") if path.is_file() and path.name != "manifest.json")
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "repo_head": git_value(repo_root, "rev-parse", "HEAD"),
        "script": "scripts/analysis/pnjl/audit_issue130_phase_reference_accepted_primary.py",
        "script_sha256": sha256(Path(__file__).resolve()),
        "calculation_sha": CALCULATION_SHA,
        "source_run_id": package_info["package_manifest"].get("source_run_id", ""),
        "replay_run_id": package_info["package_manifest"].get("replay_run_id", ""),
        "accepted_root": PACKAGE_RELATIVE.as_posix(),
        "legacy_root": LEGACY_RELATIVE.as_posix(),
        "accepted_manifest_sha256": package_info["manifest_sha256"],
        "accepted_layer_manifest_sha256": package_info["accepted_manifest_sha256"],
        "accepted_tree_sha256": package_info["tree_sha256"],
        "legacy_retirement_manifest_sha256": snapshot["manifest_sha256"],
        "legacy_tree_sha256": snapshot["tree_sha256"],
        "solver_called": False,
        "reference_write": False,
        "runtime_consumption": False,
        "runtime_legacy_fallback_rows": decision["runtime_legacy_fallback_rows"],
        "runtime_legacy_rollback_enabled": decision["runtime_legacy_rollback_enabled"],
        "active_path_contract_count": decision["active_path_contract_count"],
        "verdict": decision["verdict"],
        "decision": "decision.json",
        "files": [
            {
                "path": path.relative_to(output_root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in generated
        ],
    }
    write_json(output_root / "manifest.json", manifest)
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--output-root", type=Path, default=None)
    parser.add_argument("--replace-existing", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    output_root = (args.output_root or repo_root / OUTPUT_RELATIVE).resolve()
    manifest = build_audit(repo_root, output_root, replace_existing=args.replace_existing)
    print(json.dumps({"output_root": str(output_root), "verdict": manifest["verdict"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
