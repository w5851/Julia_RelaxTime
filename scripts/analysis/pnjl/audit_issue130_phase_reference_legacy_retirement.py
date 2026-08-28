#!/usr/bin/env python3
"""Build the solver-free Issue #130 PNJL legacy-reference retirement audit.

The candidate and legacy phase references intentionally use different schemas and
different grids.  This audit therefore compares the adapter's semantic keys
instead of CSV row numbers, verifies the byte-preserving legacy snapshot, and
scans tracked source/config/test files for active legacy consumers.  It never
calls Julia or a PNJL solver, never writes either reference tree, and never
deletes a file.
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
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping


SCHEMA_VERSION = "pnjl_issue130_phase_reference_legacy_retirement_audit_v1"
CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
CANDIDATE_RELATIVE = Path("data/reference/pnjl/issue130_phase_reference_v1")
LEGACY_RELATIVE = Path("data/reference/pnjl/legacy_phase_reference_v1")
LEGACY_MANIFEST_NAME = "RETIREMENT_MANIFEST.json"
OUTPUT_RELATIVE = Path(
    "docs/analysis/pnjl/phase_reference/issue130_phase_reference_legacy_audit_v1"
)

TABLES: dict[str, dict[str, Any]] = {
    "boundary": {
        "candidate": "strict/tables/maxwell_surface_strict_reference_v1.csv",
        "legacy": "boundary.csv",
        "keys": ("xi", "T_MeV"),
    },
    "crossover": {
        "candidate": "strict/tables/crossover_surface_strict_reference_v1.csv",
        "legacy": "crossover_dense.csv",
        "keys": ("xi", "mu_MeV"),
    },
    "cep": {
        "candidate": "strict/tables/cep_boundary_strict_reference_v1.csv",
        "legacy": "cep.csv",
        "keys": ("xi",),
    },
    "spinodals": {
        "candidate": "strict/tables/spinodal_surface_strict_reference_v1.csv",
        "legacy": "spinodals.csv",
        "keys": ("xi", "T_MeV"),
    },
}

# These are repository consumers, not a claim that external/generated callers
# have been exhaustively discovered.  Unknown active occurrences are blockers.
KNOWN_CONSUMERS: dict[str, dict[str, Any]] = {
    ".github/workflows/pnjl-phase-diagram.yml": {
        "kind": "workflow_consumer",
        "role": "legacy_input",
        "retirement_blocker": True,
    },
    ".github/workflows/relaxtime-issue130-rs-numerical-pilot-v1.yml": {
        "kind": "diagnostic_workflow",
        "role": "historical_snapshot_input",
        "retirement_blocker": True,
    },
    "scripts/relaxtime/phase_reference_adapter.jl": {
        "kind": "runtime_adapter",
        "role": "fallback_or_rollback",
        "retirement_blocker": True,
    },
    "scripts/pnjl/phase_reference_adapter.py": {
        "kind": "runtime_adapter",
        "role": "fallback_or_rollback",
        "retirement_blocker": True,
    },
    "scripts/relaxtime/run_gap_transport_scan.jl": {
        "kind": "runtime_consumer",
        "role": "candidate_with_legacy_fallback",
        "retirement_blocker": True,
    },
    "scripts/relaxtime/generate_xi_smoothness_params.jl": {
        "kind": "diagnostic_consumer",
        "role": "legacy_input",
        "retirement_blocker": True,
    },
    "scripts/relaxtime/xi_smoothness_sampling_lib.jl": {
        "kind": "diagnostic_consumer",
        "role": "legacy_input",
        "retirement_blocker": True,
    },
    "scripts/pnjl/plot_phase_diagram.py": {
        "kind": "plot_consumer",
        "role": "legacy_input",
        "retirement_blocker": True,
    },
    "scripts/pnjl/validate_phase_data.py": {
        "kind": "validation_consumer",
        "role": "legacy_input",
        "retirement_blocker": True,
    },
    "scripts/analysis/pnjl/audit_issue130_phase_reference_runtime_consumers.py": {
        "kind": "migration_audit_tooling",
        "role": "historical_path_contract",
        "retirement_blocker": True,
    },
    "scripts/analysis/pnjl/audit_issue130_phase_reference_legacy_retirement.py": {
        "kind": "migration_audit_tooling",
        "role": "current_retirement_audit",
        "retirement_blocker": False,
    },
    "scripts/analysis/pnjl/audit_issue130_phase_reference_candidate_only_consumers.py": {
        "kind": "migration_audit_tooling",
        "role": "candidate_only_consumer_audit",
        "retirement_blocker": False,
    },
    "scripts/analysis/pnjl/build_issue130_rs_runtime_consumer_smoke_v2.py": {
        "kind": "diagnostic_tooling",
        "role": "historical_snapshot_input",
        "retirement_blocker": True,
    },
    "scripts/analysis/pnjl/build_issue130_rs_runtime_parity_evidence.py": {
        "kind": "diagnostic_tooling",
        "role": "historical_snapshot_input",
        "retirement_blocker": True,
    },
    "scripts/analysis/pnjl/import_issue130_phase_reference.py": {
        "kind": "migration_tooling",
        "role": "historical_path_contract",
        "retirement_blocker": True,
    },
    "scripts/analysis/relaxtime/audit_issue130_rs_old_reference_retirement.py": {
        "kind": "migration_audit_tooling",
        "role": "historical_path_contract",
        "retirement_blocker": True,
    },
    "scripts/analysis/relaxtime/import_issue130_rs_candidate_results.py": {
        "kind": "migration_tooling",
        "role": "historical_path_contract",
        "retirement_blocker": True,
    },
    "scripts/analysis/relaxtime/build_phase_guided_pole_sensitive_rendering.py": {
        "kind": "analysis_tooling",
        "role": "historical_legacy_input",
        "retirement_blocker": True,
    },
    "scripts/dev/export_pnjl_chi_b_taylordiff_baseline.jl": {
        "kind": "baseline_export",
        "role": "legacy_cep_baseline",
        "retirement_blocker": True,
    },
    "scripts/perf/pnjl_chi_b_taylordiff_probe.jl": {
        "kind": "performance_probe",
        "role": "legacy_cep_baseline",
        "retirement_blocker": True,
    },
    "tests/regression/pnjl/test_pnjl_chi_b_taylordiff_cep_smoke.jl": {
        "kind": "regression_test",
        "role": "legacy_cep_baseline",
        "retirement_blocker": True,
    },
    "scripts/analysis/pnjl_cep_narrow_pilot.jl": {
        "kind": "analysis_tooling",
        "role": "historical_diagnostic_input",
        "retirement_blocker": False,
    },
    "scripts/analysis/pnjl_cep_narrow_pilot_v2.jl": {
        "kind": "analysis_tooling",
        "role": "historical_diagnostic_input",
        "retirement_blocker": False,
    },
    "tests/unit/python/test_issue130_phase_reference_retirement.py": {
        "kind": "snapshot_contract_test",
        "role": "snapshot_integrity_contract",
        "retirement_blocker": True,
    },
    "tests/unit/python/test_issue130_rs_numerical_pilot.py": {
        "kind": "historical_contract_test",
        "role": "historical_snapshot_input",
        "retirement_blocker": True,
    },
    "tests/unit/python/test_issue130_rs_old_reference_physical_deletion.py": {
        "kind": "snapshot_contract_test",
        "role": "snapshot_integrity_contract",
        "retirement_blocker": True,
    },
    "tests/unit/python/test_pnjl_phase_reference_adapter.py": {
        "kind": "adapter_contract_test",
        "role": "fallback_or_rollback",
        "retirement_blocker": True,
    },
    "tests/unit/relaxtime/test_phase_reference_adapter.jl": {
        "kind": "adapter_contract_test",
        "role": "fallback_or_rollback",
        "retirement_blocker": True,
    },
}

REFERENCE_TOKENS = (
    "legacy_phase_reference_v1",
    "data/reference/pnjl/boundary.csv",
    "data/reference/pnjl/cep.csv",
    "data/reference/pnjl/spinodals.csv",
    "data/reference/pnjl/crossover_dense.csv",
    "data/reference/pnjl/phase_reference_dense_manifest.json",
    "load_phase_reference_runtime_with_fallback",
    "load_default_phase_reference_runtime",
    "source=:legacy",
    "--phase-reference-mode legacy",
)

LEGACY_FILE_EXPECTATIONS = {
    "boundary.csv": ("xi", "T_MeV", "mu_transition_MeV", "rho_hadron", "rho_quark"),
    "cep.csv": ("xi", "T_CEP_MeV", "muq_CEP_MeV", "muB_CEP_MeV"),
    "crossover_dense.csv": (
        "xi",
        "mu_MeV",
        "T_crossover_MeV",
        "rho",
        "method",
        "converged",
        "derivative",
        "variable",
    ),
    "spinodals.csv": (
        "xi",
        "T_MeV",
        "mu_spinodal_hadron_MeV",
        "mu_spinodal_quark_MeV",
        "rho_spinodal_hadron",
        "rho_spinodal_quark",
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tree_hash(root: Path) -> str:
    """Hash relative paths and file hashes, independent of filesystem order."""

    digest = hashlib.sha256()
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        digest.update(path.relative_to(root).as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(sha256(path).encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        if not fields:
            raise ValueError(f"missing CSV header: {path}")
        rows = [dict(row) for row in reader]
    return fields, rows


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def write_csv(path: Path, fieldnames: Iterable[str], rows: Iterable[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    names = list(fieldnames)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in names})


def git_value(repo_root: Path, *args: str) -> str:
    try:
        return subprocess.check_output(
            ["git", *args], cwd=repo_root, text=True, stderr=subprocess.DEVNULL
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return ""


def relative(path: Path, repo_root: Path) -> str:
    return path.resolve().relative_to(repo_root.resolve()).as_posix()


def tracked_files(repo_root: Path) -> list[Path]:
    raw = subprocess.check_output(["git", "ls-files", "-z"], cwd=repo_root)
    return [repo_root / item.decode("utf-8") for item in raw.split(b"\0") if item]


def semantic_float(value: Any) -> float:
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError(f"non-finite key value: {value!r}")
    # The grids are decimal grids; twelve places avoids representation-only
    # differences while not conflating distinct physical keys.
    return round(parsed, 12)


def key_tuple(table: str, row: Mapping[str, Any]) -> tuple[float, ...]:
    return tuple(semantic_float(row[field]) for field in TABLES[table]["keys"])


def key_text(table: str, key: tuple[float, ...]) -> str:
    return json.dumps(
        {field: value for field, value in zip(TABLES[table]["keys"], key)},
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )


def load_adapter(repo_root: Path) -> Any:
    adapter_path = repo_root / "scripts" / "pnjl" / "phase_reference_adapter.py"
    spec = importlib.util.spec_from_file_location("issue130_phase_reference_adapter", adapter_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import adapter: {adapter_path}")
    module = importlib.util.module_from_spec(spec)
    # dataclasses resolve forward annotations through sys.modules while the
    # adapter is executed; register the temporary module before exec_module.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def candidate_tables(repo_root: Path, candidate_root: Path) -> tuple[dict[str, list[dict[str, Any]]], dict[str, Any]]:
    """Load candidate rows through the production-parity Python adapter."""

    adapter = load_adapter(repo_root)
    bundle = adapter.load_phase_reference(candidate_root, layer="strict", allow_runtime=False)
    tables = {table: [dict(row) for row in rows] for table, rows in bundle.tables.items()}
    root_manifest = read_json(candidate_root / "manifest.json")
    strict_manifest = read_json(candidate_root / "strict" / "manifest.json")
    return tables, {
        "root_manifest": root_manifest,
        "strict_manifest": strict_manifest,
        "manifest_sha256": sha256(candidate_root / "manifest.json"),
        "strict_manifest_sha256": sha256(candidate_root / "strict" / "manifest.json"),
        "tree_sha256": tree_hash(candidate_root),
        "adapter_schema": adapter.SCHEMA_VERSION,
    }


def legacy_tables(legacy_root: Path) -> dict[str, list[dict[str, str]]]:
    result: dict[str, list[dict[str, str]]] = {}
    for table, spec in TABLES.items():
        path = legacy_root / spec["legacy"]
        fields, rows = read_csv(path)
        missing = sorted(set(LEGACY_FILE_EXPECTATIONS[spec["legacy"]]) - set(fields))
        if missing:
            raise ValueError(f"{path} missing fields: {', '.join(missing)}")
        result[table] = rows
    return result


def verify_snapshot_manifest(repo_root: Path, legacy_root: Path) -> dict[str, Any]:
    manifest_path = legacy_root / LEGACY_MANIFEST_NAME
    manifest = read_json(manifest_path)
    errors: list[str] = []
    if manifest.get("schema_version") != "pnjl_legacy_phase_reference_retirement_v1":
        errors.append("legacy retirement schema mismatch")
    if manifest.get("status") != "retired_canonical_snapshot":
        errors.append("legacy retirement status mismatch")
    if manifest.get("canonical_root_status") != "dense_legacy_paths_absent":
        errors.append("canonical root status is not dense_legacy_paths_absent")
    records = manifest.get("files")
    if not isinstance(records, list) or not records:
        errors.append("legacy manifest has no file records")
        records = []
    files: list[dict[str, Any]] = []
    for record in records:
        if not isinstance(record, dict):
            errors.append("legacy manifest contains a non-object file record")
            continue
        rel_path = str(record.get("path", ""))
        source_path = str(record.get("source_path", ""))
        if not rel_path:
            errors.append("legacy manifest record has empty path")
            continue
        path = legacy_root / rel_path
        canonical_path = repo_root / source_path
        exists = path.is_file()
        actual_bytes = path.stat().st_size if exists else None
        actual_sha = sha256(path) if exists else ""
        byte_ok = exists and actual_bytes == record.get("bytes")
        hash_ok = exists and actual_sha == record.get("sha256")
        canonical_absent = not canonical_path.exists()
        if not byte_ok:
            errors.append(f"legacy snapshot byte mismatch: {rel_path}")
        if not hash_ok:
            errors.append(f"legacy snapshot hash mismatch: {rel_path}")
        if not canonical_absent:
            errors.append(f"canonical legacy path still exists: {source_path}")
        files.append(
            {
                "path": rel_path,
                "source_path": source_path,
                "bytes": actual_bytes if actual_bytes is not None else "",
                "expected_bytes": record.get("bytes", ""),
                "sha256": actual_sha,
                "expected_sha256": record.get("sha256", ""),
                "byte_match": byte_ok,
                "hash_match": hash_ok,
                "canonical_absent": canonical_absent,
            }
        )
    return {
        "manifest": manifest,
        "manifest_sha256": sha256(manifest_path),
        "tree_sha256": tree_hash(legacy_root),
        "files": files,
        "errors": errors,
        "integrity_pass": not errors,
    }


def coverage_matrix(
    candidate: Mapping[str, list[Mapping[str, Any]]],
    legacy: Mapping[str, list[Mapping[str, str]]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    summary: list[dict[str, Any]] = []
    fallback_rows: list[dict[str, Any]] = []
    for table in TABLES:
        candidate_rows = list(candidate.get(table, ()))
        legacy_rows = list(legacy.get(table, ()))
        all_keys = [key_tuple(table, row) for row in candidate_rows]
        certified_keys = {
            key_tuple(table, row)
            for row in candidate_rows
            if bool(row.get("certified", False))
        }
        candidate_duplicate_count = len(all_keys) - len(set(all_keys))
        legacy_keys = [key_tuple(table, row) for row in legacy_rows]
        legacy_duplicate_count = len(legacy_keys) - len(set(legacy_keys))
        overlap = 0
        fallback = 0
        for index, (row, key) in enumerate(zip(legacy_rows, legacy_keys), start=2):
            matched = key in certified_keys
            if matched:
                overlap += 1
                reason = "certified_candidate_key_overlap"
                disposition = "candidate_preferred"
            else:
                fallback += 1
                reason = "candidate_key_absent_or_uncertified"
                disposition = "legacy_fallback_required"
                fallback_rows.append(
                    {
                        "table": table,
                        "legacy_row": index,
                        "key": key_text(table, key),
                        "reason": reason,
                        "disposition": disposition,
                    }
                )
            # Keep the overlap rows out of fallback_matrix.csv; the summary and
            # key counts retain the complete relation without duplicating data.
        uncertified = sum(
            1 for row in candidate_rows if not bool(row.get("certified", False))
        )
        summary.append(
            {
                "table": table,
                "candidate_rows": len(candidate_rows),
                "candidate_certified_rows": len(certified_keys),
                "candidate_uncertified_rows": uncertified,
                "candidate_duplicate_keys": candidate_duplicate_count,
                "legacy_rows": len(legacy_rows),
                "legacy_duplicate_keys": legacy_duplicate_count,
                "legacy_certified_overlap": overlap,
                "legacy_fallback_required": fallback,
                "key_fields": ",".join(TABLES[table]["keys"]),
                "key_normalization": "float rounded to 12 decimal places",
            }
        )
    return summary, fallback_rows


def classify_path(path: str) -> str:
    if path.startswith("docs/dev/archived/") or path.startswith("docs/analysis/"):
        return "historical_evidence"
    if path.startswith("data/reference/pnjl/legacy_phase_reference_v1/"):
        return "snapshot_metadata"
    if path.startswith("config/governance/"):
        return "active_governance"
    if path.startswith(("scripts/", "src/", "tests/", "config/", ".github/")):
        return "active_repository_contract"
    if path.startswith("docs/dev/active/"):
        return "active_governance"
    return "other_reference"


def occurrence_role(path: str, snippet: str) -> str:
    if path in KNOWN_CONSUMERS:
        return str(KNOWN_CONSUMERS[path]["role"])
    lower = snippet.lower()
    if "audit" in path or "import" in path or "migration" in lower:
        return "migration_or_audit_tooling"
    if path.startswith("tests/") or path.startswith("scripts/perf/"):
        return "baseline_or_regression"
    if "plot" in path or "validate" in path:
        return "visualization_or_validation"
    if "rollback" in lower or "fallback" in lower or "source=:legacy" in lower:
        return "fallback_or_rollback"
    return "legacy_path_contract"


def scan_references(repo_root: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    occurrences: list[dict[str, Any]] = []
    by_path: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for path in tracked_files(repo_root):
        path_text = relative(path, repo_root)
        try:
            lines = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
        except OSError:
            continue
        for line_number, line in enumerate(lines, start=1):
            token = next((value for value in REFERENCE_TOKENS if value in line), None)
            if token is None:
                continue
            row = {
                "path": path_text,
                "line": line_number,
                "token": token,
                "classification": classify_path(path_text),
                "role": occurrence_role(path_text, line),
                "snippet": re.sub(r"\s+", " ", line.strip())[:240],
            }
            occurrences.append(row)
            by_path[path_text].append(row)

    consumers: list[dict[str, Any]] = []
    all_paths = set(by_path) | set(KNOWN_CONSUMERS)
    for path_text in sorted(all_paths):
        rows = by_path.get(path_text, [])
        known = KNOWN_CONSUMERS.get(path_text)
        classification = classify_path(path_text)
        active = classification in {"active_repository_contract", "active_governance"}
        if known is not None:
            kind = str(known["kind"])
            role = str(known["role"])
            blocker = bool(known["retirement_blocker"])
            decision = "migration_required_before_delete" if blocker else "retain_as_history_or_reclassify"
        elif classification == "active_governance":
            kind = "governance_reference"
            role = "task_or_governance_record"
            blocker = False
            decision = "retain_as_governance_evidence"
        elif active:
            kind = "unclassified_active_reference"
            role = occurrence_role(path_text, rows[0]["snippet"] if rows else "")
            blocker = True
            decision = "manual_review_required"
        else:
            kind = "historical_or_metadata_reference"
            role = occurrence_role(path_text, rows[0]["snippet"] if rows else "")
            blocker = False
            decision = "not_an_active_runtime_consumer"
        consumers.append(
            {
                "path": path_text,
                "kind": kind,
                "role": role,
                "classification": classification,
                "occurrence_count": len(rows),
                "reference_tokens": ",".join(sorted({str(row["token"]) for row in rows})),
                "active": active,
                "retirement_blocker": blocker,
                "decision": decision,
            }
        )
    return occurrences, consumers


def compute_decision(
    candidate_info: Mapping[str, Any],
    legacy_info: Mapping[str, Any],
    coverage: Iterable[Mapping[str, Any]],
    consumers: Iterable[Mapping[str, Any]],
) -> dict[str, Any]:
    coverage_rows = list(coverage)
    consumer_rows = list(consumers)
    fallback_total = sum(int(row["legacy_fallback_required"]) for row in coverage_rows)
    candidate_uncertified_total = sum(int(row["candidate_uncertified_rows"]) for row in coverage_rows)
    duplicate_total = sum(
        int(row["candidate_duplicate_keys"]) + int(row["legacy_duplicate_keys"])
        for row in coverage_rows
    )
    active_blockers = [row for row in consumer_rows if bool(row["retirement_blocker"])]
    unknown_active = [
        row for row in consumer_rows if row["kind"] == "unclassified_active_reference"
    ]
    candidate_root = candidate_info["root_manifest"]
    strict_manifest = candidate_info["strict_manifest"]
    provenance_ok = (
        candidate_root.get("calculation_sha") == CALCULATION_SHA
        and candidate_root.get("runtime_consumption") is False
        and strict_manifest.get("calculation_sha") == CALCULATION_SHA
        and strict_manifest.get("solver_called") is False
        and strict_manifest.get("reference_write") is False
    )
    integrity_ok = bool(legacy_info["integrity_pass"]) and provenance_ok and duplicate_total == 0
    # Only candidate rows that correspond to a legacy key can require the
    # snapshot.  Uncertified candidate-only rows remain a separate diagnostic
    # quality issue and must not be counted as a phantom legacy dependency.
    coverage_complete = fallback_total == 0
    consumers_migrated = not active_blockers and not unknown_active
    physical_deletion_eligible = integrity_ok and coverage_complete and consumers_migrated
    if not integrity_ok:
        verdict = "legacy_audit_input_invalid"
        next_action = "repair immutable input/hash/provenance errors; retain the legacy snapshot"
    elif not coverage_complete or not consumers_migrated:
        verdict = "retirement_inconclusive"
        next_action = "retain snapshot; migrate consumers and/or close candidate certified-key gaps before path retirement"
    else:
        verdict = "legacy_physical_deletion_candidate"
        next_action = "prepare a separate deletion allowlist and request explicit author authorization"
    return {
        "schema_version": SCHEMA_VERSION,
        "verdict": verdict,
        "solver_called": False,
        "reference_write": False,
        "runtime_consumption": False,
        "candidate_auto_selected": False,
        "candidate_provenance_ok": provenance_ok,
        "legacy_snapshot_integrity_pass": bool(legacy_info["integrity_pass"]),
        "candidate_duplicate_or_legacy_duplicate_keys": duplicate_total,
        "candidate_uncertified_rows": candidate_uncertified_total,
        "legacy_fallback_key_count": fallback_total,
        "active_consumer_blocker_count": len(active_blockers),
        "unknown_active_reference_count": len(unknown_active),
        "coverage_complete": coverage_complete,
        "consumers_migrated": consumers_migrated,
        "physical_deletion_eligible": physical_deletion_eligible,
        "calculation_sha": CALCULATION_SHA,
        "source_run_id": str(candidate_root.get("source_run_id", "")),
        "replay_run_id": str(candidate_root.get("replay_run_id", "")),
        "stop_reasons": [
            reason
            for reason, condition in (
                ("candidate_or_legacy_key_coverage_incomplete", not coverage_complete),
                ("active_consumer_or_rollback_dependency_present", not consumers_migrated),
                ("unknown_active_reference_requires_review", bool(unknown_active)),
            )
            if condition
        ],
        "next_action": next_action,
        "non_goals": [
            "do not call Julia or the PNJL equilibrium solver",
            "do not rewrite candidate or legacy reference files",
            "do not delete the PNJL legacy snapshot",
            "do not treat a solver-free audit as numerical convergence",
        ],
    }


def build_claims(decision: Mapping[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "claim_id": "candidate_strict_integrity",
            "claim": "Candidate strict tables are finite and key-unique under the production-parity adapter.",
            "status": "supported" if decision["candidate_provenance_ok"] and decision["candidate_duplicate_or_legacy_duplicate_keys"] == 0 else "inconclusive",
            "evidence": "tables/coverage_matrix.csv; candidate manifest hashes",
            "boundary": "does not authorize runtime-only or physical deletion",
        },
        {
            "claim_id": "legacy_snapshot_integrity",
            "claim": "The versioned PNJL legacy snapshot matches its byte/hash manifest and canonical paths are absent.",
            "status": "supported" if decision["legacy_snapshot_integrity_pass"] else "not_supported",
            "evidence": "tables/legacy_snapshot_inventory.csv; data/reference/pnjl/legacy_phase_reference_v1/RETIREMENT_MANIFEST.json",
            "boundary": "integrity preserves the fallback; it is not evidence that fallback is unused",
        },
        {
            "claim_id": "candidate_replaces_all_legacy_keys",
            "claim": "Certified candidate keys fully replace all legacy keys.",
            "status": "supported" if decision["coverage_complete"] else "not_supported",
            "evidence": "tables/coverage_matrix.csv; tables/fallback_matrix.csv",
            "boundary": "candidate uncertified rows remain excluded from runtime view",
        },
        {
            "claim_id": "active_consumers_migrated",
            "claim": "All active repository consumers and explicit rollback paths are independent of the legacy snapshot.",
            "status": "supported" if decision["consumers_migrated"] else "not_supported",
            "evidence": "tables/consumer_matrix.csv; tables/reference_occurrences.csv",
            "boundary": "static tracked-source scan; external/generated consumers are not proven absent",
        },
        {
            "claim_id": "physical_deletion_ready",
            "claim": "The PNJL legacy snapshot is ready for physical deletion.",
            "status": "supported" if decision["physical_deletion_eligible"] else "not_supported",
            "evidence": "decision.json; tables/coverage_matrix.csv; tables/consumer_matrix.csv",
            "boundary": "even a candidate verdict still requires a separate allowlist PR and author authorization",
        },
    ]


def write_analysis_text(
    output_root: Path,
    decision: Mapping[str, Any],
    coverage: list[Mapping[str, Any]],
    consumers: list[Mapping[str, Any]],
    repo_root: Path,
) -> None:
    fallback_total = int(decision["legacy_fallback_key_count"])
    active_blockers = int(decision["active_consumer_blocker_count"])
    rows = [
        "# Issue #130 PNJL legacy phase-reference retirement audit v1",
        "",
        "这是 solver-free 的阶段 A 审计，不调用 Julia/PNJL，不修改 reference，不删除文件。",
        f"本次 repo HEAD：`{git_value(repo_root, 'rev-parse', 'HEAD')}`。",
        "",
        "## Verdict",
        "",
        f"- verdict：`{decision['verdict']}`",
        f"- legacy fallback keys：`{fallback_total}`",
        f"- active consumer/rollback blockers：`{active_blockers}`",
        f"- candidate uncertified rows：`{decision['candidate_uncertified_rows']}`",
        f"- physical deletion eligible：`{decision['physical_deletion_eligible']}`",
        f"- next action：{decision['next_action']}",
        "",
        "## 关键语义",
        "",
        "- 覆盖使用 adapter 的真实语义 key：boundary/spinodals 为 `(xi,T_MeV)`，",
        "  crossover 为 `(xi,mu_MeV)`，CEP 为 `(xi)`；浮点 key 统一四舍五入到 12 位。",
        "- candidate 的 unresolved/non-certified 行不会被视为 certified；因此 legacy fallback",
        "  只要仍有缺键或未认证键就不能物理删除。",
        "- 静态 consumer 扫描只覆盖当前 Git tracked source/config/test；不能证明外部生成器、",
        "  本地脚本副本或未来代码不存在引用。",
        "",
        "## 四表起点",
        "",
        "| table | candidate rows | candidate certified | candidate uncertified | legacy rows | fallback keys |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in coverage:
        rows.append(
            f"| {row['table']} | {row['candidate_rows']} | {row['candidate_certified_rows']} | "
            f"{row['candidate_uncertified_rows']} | {row['legacy_rows']} | {row['legacy_fallback_required']} |"
        )
    rows.extend(
        [
            "",
            "## Consumer 处置",
            "",
            "`tables/consumer_matrix.csv` 给出每个 tracked occurrence 的角色和是否阻塞物理清理；",
            "`tables/claim_ledger.json` 保留机器可读的 claim/boundary 对象，CSV 供表格审阅。",
            f"当前共记录 `{len(consumers)}` 个含 legacy token 的路径，其中 `{active_blockers}` 个仍需迁移或保留 rollback。",
            "",
            "## 下一阶段",
            "",
            "1. 先保留 `data/reference/pnjl/legacy_phase_reference_v1/`，不创建删除 PR。",
            "2. 对 blocker consumer 设计 candidate-only/explicit-rollback contract，并增加 focused tests。",
            "3. 迁移完成后重跑本 audit；只有 fallback=0 且未知 active 引用为 0 才能生成 physical-deletion proposal。",
        ]
    )
    (output_root / "README.md").write_text("\n".join(rows) + "\n", encoding="utf-8", newline="\n")

    audit_lines = [
        "# PNJL legacy phase-reference retirement audit",
        "",
        "## 结论",
        "",
        f"本次 verdict 为 `{decision['verdict']}`。当前仍有 `{fallback_total}` 个 legacy key 未被 certified candidate 精确覆盖，",
        f"并记录 `{active_blockers}` 个 active consumer/rollback blocker，因此本次不具备 physical deletion 条件。",
        "",
        "## 保留与回退",
        "",
        "legacy snapshot 继续作为 candidate 未认证/缺键行的逐键 fallback、显式 legacy rollback 和历史复现输入。",
        "这与 RS prod_v1 的物理删除是不同任务；RS allowlist 不得用于 PNJL 文件。",
        "",
        "## 证据边界",
        "",
        "候选数据通过现有 Python phase-reference adapter 读取，审计过程设置 `solver_called=false`、",
        "`reference_write=false`。本包只保存摘要、key 覆盖和引用矩阵，不复制全量曲线。",
        "",
        "## 停止条件",
        "",
        "若 hash/bytes 不一致、发现未知 active consumer、candidate certified coverage 下降，或 rollback 不可用，",
        "停止 path retirement 和物理删除，恢复到当前 immutable snapshot。",
    ]
    (output_root / "AUDIT.md").write_text("\n".join(audit_lines) + "\n", encoding="utf-8", newline="\n")


def build_audit(
    repo_root: Path,
    output_root: Path,
    *,
    candidate_root: Path | None = None,
    legacy_root: Path | None = None,
) -> dict[str, Any]:
    if output_root.exists() and any(output_root.iterdir()):
        raise FileExistsError(f"refusing to overwrite non-empty audit output: {output_root}")
    output_root.mkdir(parents=True, exist_ok=True)
    candidate_root = (candidate_root or repo_root / CANDIDATE_RELATIVE).resolve()
    legacy_root = (legacy_root or repo_root / LEGACY_RELATIVE).resolve()

    candidate, candidate_info = candidate_tables(repo_root, candidate_root)
    legacy = legacy_tables(legacy_root)
    legacy_info = verify_snapshot_manifest(repo_root, legacy_root)
    coverage, fallback_rows = coverage_matrix(candidate, legacy)
    occurrences, consumers = scan_references(repo_root)
    decision = compute_decision(candidate_info, legacy_info, coverage, consumers)

    tables_root = output_root / "tables"
    write_csv(tables_root / "coverage_matrix.csv", list(coverage[0].keys()), coverage)
    write_csv(
        tables_root / "fallback_matrix.csv",
        ["table", "legacy_row", "key", "reason", "disposition"],
        fallback_rows,
    )
    write_csv(
        tables_root / "legacy_snapshot_inventory.csv",
        list(legacy_info["files"][0].keys()) if legacy_info["files"] else ["path"],
        legacy_info["files"],
    )
    write_csv(
        tables_root / "reference_occurrences.csv",
        ["path", "line", "token", "classification", "role", "snippet"],
        occurrences,
    )
    write_csv(
        tables_root / "consumer_matrix.csv",
        list(consumers[0].keys()) if consumers else ["path"],
        consumers,
    )
    claims = build_claims(decision)
    write_csv(
        tables_root / "claim_ledger.csv",
        list(claims[0].keys()),
        claims,
    )
    # Keep a machine-readable JSON form alongside the tabular view.  The CSV
    # is convenient for review; the JSON preserves the claim objects without
    # flattening evidence/boundary fields and matches the task-contract name.
    write_json(tables_root / "claim_ledger.json", claims)
    write_json(output_root / "decision.json", decision)
    write_analysis_text(output_root, decision, coverage, consumers, repo_root)

    generated_files = sorted(
        path for path in output_root.rglob("*") if path.is_file() and path.name != "manifest.json"
    )
    root_manifest = candidate_info["root_manifest"]
    strict_manifest = candidate_info["strict_manifest"]
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "repo_head": git_value(repo_root, "rev-parse", "HEAD"),
        "calculation_sha": CALCULATION_SHA,
        "source_run_id": str(root_manifest.get("source_run_id", "")),
        "replay_run_id": str(root_manifest.get("replay_run_id", "")),
        "candidate_root": relative(candidate_root, repo_root),
        "legacy_root": relative(legacy_root, repo_root),
        "candidate_manifest_sha256": candidate_info["manifest_sha256"],
        "candidate_strict_manifest_sha256": candidate_info["strict_manifest_sha256"],
        "candidate_tree_sha256": candidate_info["tree_sha256"],
        "legacy_retirement_manifest_sha256": legacy_info["manifest_sha256"],
        "legacy_tree_sha256": legacy_info["tree_sha256"],
        "candidate_runtime_consumption": root_manifest.get("runtime_consumption"),
        "candidate_reference_write": root_manifest.get("reference_write"),
        "candidate_strict_reference_write": strict_manifest.get("reference_write"),
        "solver_called": False,
        "reference_write": False,
        "runtime_consumption": False,
        "candidate_auto_selected": False,
        "decision": "decision.json",
        "verdict": decision["verdict"],
        "legacy_fallback_key_count": decision["legacy_fallback_key_count"],
        "active_consumer_blocker_count": decision["active_consumer_blocker_count"],
        "files": [
            {
                "path": relative(path, output_root),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in generated_files
        ],
        "input_manifests": {
            "candidate_root": relative(candidate_root / "manifest.json", repo_root),
            "candidate_strict": relative(candidate_root / "strict" / "manifest.json", repo_root),
            "legacy_retirement": relative(legacy_root / LEGACY_MANIFEST_NAME, repo_root),
        },
    }
    write_json(output_root / "manifest.json", manifest)
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--candidate-root", type=Path, default=None)
    parser.add_argument("--legacy-root", type=Path, default=None)
    parser.add_argument("--output-root", type=Path, default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    candidate_root = args.candidate_root.resolve() if args.candidate_root is not None else None
    legacy_root = args.legacy_root.resolve() if args.legacy_root is not None else None
    output_root = (args.output_root or repo_root / OUTPUT_RELATIVE).resolve()
    manifest = build_audit(
        repo_root,
        output_root,
        candidate_root=candidate_root,
        legacy_root=legacy_root,
    )
    print(json.dumps({"output_root": str(output_root), "verdict": manifest["verdict"]}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
