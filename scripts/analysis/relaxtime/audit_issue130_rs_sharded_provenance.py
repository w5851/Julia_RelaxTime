#!/usr/bin/env python3
"""Audit the immutable Issue #130 RS shard aggregate without running a solver.

The 2026-08-26 shard artifacts were produced before PR #269.  This audit keeps
those files byte-for-byte unchanged and classifies the two known sidecar hash
defects separately from missing or unexpected artifact corruption.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from collections import Counter
from pathlib import Path
from typing import Any, Iterable


SCHEMA = "issue130_rs_sharded_production_v2_post_repair_audit_v1"
CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
WORKFLOW_HEAD_SHA = "22874505877491754eed27519ad8a7b871c82571"
CASE_NAME = "first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2"
REPAIR_MERGE_SHA = "472e3eb0047b3a9380b27f1935880b0e473ca9b5"
EXPECTED_RUNS = 30
ALLOWED_HISTORICAL_MISMATCHES = {
    "effective_config.json": "written after the pre-repair manifest hash snapshot",
    "failed_points.csv": "header stream was flushed after the pre-repair manifest hash snapshot",
}
NUMERIC_FIELDS = (
    "T_MeV",
    "muq_MeV",
    "muB_MeV",
    "xi",
    "alpha_T",
    "eta_over_s",
    "sigma_over_T",
    "zeta_over_s",
    "tau_u",
    "tau_d",
    "tau_s",
    "tau_ubar",
    "tau_dbar",
    "tau_sbar",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line.strip() and not line.startswith("#")]
    if not lines:
        return [], []
    reader = csv.DictReader(lines)
    if reader.fieldnames is None:
        raise ValueError(f"missing CSV header: {path}")
    return list(reader.fieldnames), list(reader)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_csv(path: Path, fieldnames: Iterable[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def git_head(repo_root: Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return ""


def as_float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if number == number and abs(number) != float("inf") else None


def find_run_manifest(source_root: Path, run_id: str) -> Path:
    candidates = sorted((source_root / "collector" / run_id).rglob("run_manifest.json"))
    if not candidates:
        candidates = sorted((source_root / f"run_{run_id}").rglob("run_manifest.json"))
    if len(candidates) != 1:
        raise FileNotFoundError(f"expected one run_manifest.json for {run_id}, found {len(candidates)}")
    return candidates[0]


def resolve_artifact(manifest_path: Path, artifact_path: str) -> Path | None:
    normalized = artifact_path.replace("\\", "/")
    candidates: list[Path] = []
    if "data/outputs/results/" in normalized:
        suffix = normalized.split("data/outputs/results/", 1)[1]
        candidates.append(manifest_path.parent.parent.parent.parent / suffix)
    if "data/outputs/" in normalized:
        suffix = normalized.split("data/outputs/", 1)[1]
        candidates.append(manifest_path.parent.parent.parent.parent / suffix)
    candidates.append(manifest_path.parent / Path(normalized).name)
    candidates.extend(sorted(manifest_path.parent.rglob(Path(normalized).name)))
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def artifact_hash_rows(manifest_path: Path) -> list[dict[str, Any]]:
    manifest = read_json(manifest_path)
    rows: list[dict[str, Any]] = []
    for artifact in manifest.get("artifacts", []):
        relative = str(artifact.get("path", ""))
        actual = resolve_artifact(manifest_path, relative)
        expected_hash = str(artifact.get("sha256", ""))
        actual_hash = sha256(actual) if actual else ""
        basename = Path(relative.replace("\\", "/")).name
        if actual is None:
            verdict = "missing_artifact"
            mismatch_kind = "missing"
        elif expected_hash == actual_hash:
            verdict = "hash_match"
            mismatch_kind = "none"
        elif basename in ALLOWED_HISTORICAL_MISMATCHES:
            verdict = "historical_known_mismatch"
            mismatch_kind = basename
        else:
            verdict = "unexpected_hash_mismatch"
            mismatch_kind = "unexpected"
        rows.append(
            {
                "manifest": str(manifest_path),
                "path": relative,
                "basename": basename,
                "expected_sha256": expected_hash,
                "actual_sha256": actual_hash,
                "verdict": verdict,
                "mismatch_kind": mismatch_kind,
            }
        )
    return rows


def scan_summary(manifest_path: Path, expected: dict[str, str]) -> dict[str, Any]:
    rows = artifact_hash_rows(manifest_path)
    scan_path = next((Path(row["manifest"]).parent / Path(row["path"].replace("\\", "/")).name for row in rows if row["basename"] == "phase_guided_transport_scan.csv"), None)
    if scan_path is None or not scan_path.is_file():
        scan_path = resolve_artifact(manifest_path, next((row["path"] for row in rows if row["basename"] == "phase_guided_transport_scan.csv"), ""))
    fields, scan_rows = read_csv(scan_path) if scan_path else ([], [])
    duplicates: set[tuple[str, ...]] = set()
    seen: set[tuple[str, ...]] = set()
    key_fields = ("muB_MeV", "xi", "alpha_T") if expected.get("mode", "").startswith("mode_a") else ("T_MeV", "muB_MeV", "xi")
    for row in scan_rows:
        key = tuple(row.get(field, "") for field in key_fields)
        if key in seen:
            duplicates.add(key)
        seen.add(key)
    nonfinite = []
    numeric_fields = tuple(
        field
        for field in NUMERIC_FIELDS
        if not (field == "alpha_T" and expected.get("mode", "").startswith("mode_b"))
    )
    for index, row in enumerate(scan_rows):
        for field in numeric_fields:
            if field in row and as_float(row[field]) is None:
                nonfinite.append((index, field))
    nonconverged = [row for row in scan_rows if row.get("converged", "").lower() != "true"]
    failed_path = next((resolve_artifact(manifest_path, row["path"]) for row in rows if row["basename"] == "failed_points.csv"), None)
    _, failed_rows = read_csv(failed_path) if failed_path else ([], [])
    return {
        "run_id": expected.get("run_id", ""),
        "label": expected.get("label", ""),
        "mode": expected.get("mode", ""),
        "manifest": str(manifest_path),
        "scan_rows": len(scan_rows),
        "expected_scan_rows": int(expected.get("scan_rows", 0) or 0),
        "scan_duplicate_keys": len(duplicates),
        "nonfinite_numeric_values": len(nonfinite),
        "allowed_alpha_T_nan_count": sum(
            str(row.get("alpha_T", "")).lower() == "nan"
            for row in scan_rows
            if expected.get("mode", "").startswith("mode_b")
        ),
        "nonconverged_rows": len(nonconverged),
        "failed_rows": len(failed_rows),
        "hash_mismatch_count": sum(row["verdict"] != "hash_match" for row in rows),
        "hash_unexpected_mismatch_count": sum(row["verdict"] == "unexpected_hash_mismatch" for row in rows),
        "hash_missing_count": sum(row["verdict"] == "missing_artifact" for row in rows),
        "hash_known_mismatch_files": ",".join(sorted({row["basename"] for row in rows if row["verdict"] == "historical_known_mismatch"})),
        "git_commit": expected.get("git_commit", ""),
        "workflow_sha_in_audit": expected.get("workflow_sha_in_audit", "").lower() == "true",
        "calc_sha_in_audit": expected.get("calc_sha_in_audit", "").lower() == "true",
        "hard_gate_ok": expected.get("hard_gate_ok", "").lower() == "true",
        "timing_ok": expected.get("timing_ok", "").lower() == "true",
    }


def load_summary(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    fields, rows = read_csv(path)
    return fields, rows


def normalized_key(row: dict[str, str], fields: tuple[str, ...]) -> tuple[str, ...]:
    values = []
    for field in fields:
        number = as_float(row.get(field, ""))
        values.append(f"{number:.12g}" if number is not None else row.get(field, ""))
    return tuple(values)


def candidate_legacy_review(aggregate_root: Path, aggregate_summary: dict[str, Any], output_root: Path) -> dict[str, Any]:
    numeric_rows = []
    comparison_path = aggregate_root / "candidate_legacy_comparison.csv"
    if comparison_path.is_file():
        _, comparison = read_csv(comparison_path)
        grouped: dict[tuple[str, str], list[dict[str, Any]]] = {}
        for row in comparison:
            abs_delta = as_float(row.get("abs_delta", ""))
            rel_delta = as_float(row.get("rel_delta", ""))
            if abs_delta is None:
                continue
            grouped.setdefault((row.get("mode", ""), row.get("column", "")), []).append({
                "mode": row.get("mode", ""),
                "column": row.get("column", ""),
                "key": row.get("key", ""),
                "abs_delta": abs_delta,
                "rel_delta": rel_delta if rel_delta is not None else 0.0,
            })
        for (mode, column), rows in sorted(grouped.items()):
            by_abs = max(rows, key=lambda item: item["abs_delta"])
            by_rel = max(rows, key=lambda item: item["rel_delta"])
            numeric_rows.append({
                "mode": mode,
                "column": column,
                "comparison_count": len(rows),
                "max_abs_delta": by_abs["abs_delta"],
                "max_abs_key": by_abs["key"],
                "max_rel_delta": by_rel["rel_delta"],
                "max_rel_key": by_rel["key"],
            })
    write_csv(output_root / "candidate_legacy_numeric_extrema.csv", numeric_rows[0].keys() if numeric_rows else ("mode", "column"), numeric_rows)

    phase_rows: list[dict[str, Any]] = []
    classification_fields = ("phase_reference_kind", "phase_structure", "quality_flag", "quality_reason")
    route_fields = ("phase_prev", "phase_curr", "phase_boundary_xi_used", "seed_source", "equilibrium_backend")
    phase_fields = classification_fields + route_fields
    mode_configs = {
        "mode_a_fixed_muB_phase_scaled": (("muB_MeV", "xi", "alpha_T"), aggregate_root / "mode_a_scan.csv"),
        "mode_b_fixed_T_sparse_muB": (("T_MeV", "muB_MeV", "xi"), aggregate_root / "mode_b_scan.csv"),
    }
    phase_summary: dict[str, Any] = {}
    for mode, (key_fields, candidate_path) in mode_configs.items():
        legacy_path = Path(str(aggregate_summary.get("candidate_legacy_summary", {}).get(mode, {}).get("legacy_path", "")))
        _, candidate_rows = read_csv(candidate_path) if candidate_path.is_file() else ([], [])
        _, legacy_rows = read_csv(legacy_path) if legacy_path.is_file() else ([], [])
        candidate_index = {normalized_key(row, key_fields): row for row in candidate_rows}
        legacy_index = {normalized_key(row, key_fields): row for row in legacy_rows}
        common = sorted(set(candidate_index) & set(legacy_index))
        classification_keys: set[tuple[str, ...]] = set()
        route_keys: set[tuple[str, ...]] = set()
        for key in common:
            candidate = candidate_index[key]
            legacy = legacy_index[key]
            for field in phase_fields:
                left = candidate.get(field, "")
                right = legacy.get(field, "")
                if left != right:
                    if field in classification_fields:
                        classification_keys.add(key)
                        classification = "phase_classification_mismatch_author_check"
                    else:
                        route_keys.add(key)
                        classification = "reference_route_difference_author_check"
                    phase_rows.append({
                        "mode": mode,
                        "key": "|".join(key),
                        "field": field,
                        "legacy": right,
                        "candidate": left,
                        "classification": classification,
                    })
        phase_summary[mode] = {
            "candidate_rows": len(candidate_rows),
            "legacy_rows": len(legacy_rows),
            "common_keys": len(common),
            "candidate_only_keys": len(set(candidate_index) - set(legacy_index)),
            "legacy_only_keys": len(set(legacy_index) - set(candidate_index)),
            "phase_classification_mismatch_keys": len(classification_keys),
            "route_difference_keys": len(route_keys),
            "legacy_exists": legacy_path.is_file(),
        }
    write_csv(output_root / "candidate_legacy_phase_mismatches.csv", phase_rows[0].keys() if phase_rows else ("mode", "key", "field"), phase_rows)
    review = {
        "schema_version": "issue130_rs_candidate_legacy_review_v1",
        "verdict": "author_check_required",
        "numeric_extrema_file": "candidate_legacy_numeric_extrema.csv",
        "phase_mismatches_file": "candidate_legacy_phase_mismatches.csv",
        "phase_summary": phase_summary,
        "non_goals": ["no automatic acceptance", "no solver rerun", "no production write"],
    }
    write_json(output_root / "candidate_legacy_review_summary.json", review)
    return review


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--aggregate-version", default="aggregate_replay_20260826_v4")
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[3])
    args = parser.parse_args()

    source_root = args.source_root.resolve()
    output_root = args.output_root.resolve()
    aggregate_root = source_root / args.aggregate_version
    summary_path = aggregate_root / "aggregate_summary.json"
    aggregate_manifest_path = aggregate_root / "manifest.json"
    per_shard_path = aggregate_root / "per_shard_summary.csv"
    for path in (summary_path, aggregate_manifest_path, per_shard_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    aggregate_summary = read_json(summary_path)
    aggregate_manifest = read_json(aggregate_manifest_path)
    manifest_rows = []
    for entry in aggregate_manifest.get("files", []):
        path = aggregate_root / str(entry["path"])
        actual_hash = sha256(path) if path.is_file() else ""
        manifest_rows.append(
            {
                "path": str(entry["path"]),
                "expected_bytes": int(entry.get("bytes", 0)),
                "actual_bytes": path.stat().st_size if path.is_file() else 0,
                "expected_sha256": str(entry.get("sha256", "")),
                "actual_sha256": actual_hash,
                "verdict": "match" if path.is_file() and actual_hash == entry.get("sha256") and path.stat().st_size == entry.get("bytes") else "mismatch",
            }
        )
    _, shard_rows = load_summary(per_shard_path)
    shard_audits = []
    all_artifact_rows = []
    for row in shard_rows:
        run_id = row.get("run_id", "")
        manifest_path = find_run_manifest(source_root, run_id)
        all_artifact_rows.extend(artifact_hash_rows(manifest_path))
        shard_audits.append(scan_summary(manifest_path, row))

    unexpected = [row for row in all_artifact_rows if row["verdict"] in {"unexpected_hash_mismatch", "missing_artifact"}]
    numerical_failures = [
        row
        for row in shard_audits
        if row["scan_rows"] != row["expected_scan_rows"]
        or row["scan_duplicate_keys"]
        or row["nonfinite_numeric_values"]
        or row["nonconverged_rows"]
        or row["failed_rows"]
        or not row["hard_gate_ok"]
        or not row["timing_ok"]
        or not row["calc_sha_in_audit"]
        or not row["workflow_sha_in_audit"]
    ]
    aggregate_manifest_ok = all(row["verdict"] == "match" for row in manifest_rows)
    aggregate_contract_ok = (
        aggregate_summary.get("selected_unique_run_count") == EXPECTED_RUNS
        and aggregate_summary.get("per_shard_hard_gate_all_ok") is True
        and aggregate_summary.get("actions_timing_all_ok") is True
        and aggregate_summary.get("all_hard_gates_including_timing") is True
        and aggregate_summary.get("direct_coexistence_gate", {}).get("has_minus_0003") is True
        and aggregate_summary.get("direct_coexistence_gate", {}).get("has_plus_0003") is True
        and aggregate_summary.get("direct_coexistence_gate", {}).get("has_zero") is False
    )
    verdict = "post_repair_audit_pass_diagnostic_only" if aggregate_manifest_ok and aggregate_contract_ok and not unexpected and not numerical_failures else "post_repair_audit_inconclusive"
    stop_reason = (
        "historical_sidecar_hash_defects_preserved_immutable; producer_contract_fixed_in_pr_269"
        if verdict == "post_repair_audit_pass_diagnostic_only"
        else "unexpected_or_numerical_audit_failure"
    )
    candidate_legacy = aggregate_summary.get("candidate_legacy_summary", {})
    quality_counts = aggregate_summary.get("quality_reason_counts", {})
    claims = [
        {
            "claim_id": "all_unique_shards_complete",
            "claim": "The selected 30 unique RS shards have complete finite converged scans and pass the recorded hard gates.",
            "status": "supported" if not numerical_failures and aggregate_contract_ok else "inconclusive",
            "evidence": "shard_audit.csv; aggregate_summary.json",
            "boundary": "solver-backed execution evidence remains diagnostic-only",
        },
        {
            "claim_id": "historical_manifest_defects_are_classified",
            "claim": "All historical shard manifest mismatches are limited to effective_config.json and failed_points.csv and are retained as immutable provenance evidence.",
            "status": "supported" if not unexpected else "inconclusive",
            "evidence": "artifact_hash_audit.csv; PR #269 merge SHA",
            "boundary": "the old manifests are not rewritten",
        },
        {
            "claim_id": "candidate_legacy_difference_requires_review",
            "claim": "Candidate/legacy phase and transport differences remain a separate author-review gate rather than an automatic promotion verdict.",
            "status": "author_check",
            "evidence": "aggregate_summary.json; candidate_legacy_comparison_summary.json",
            "boundary": "reference drift is not automatically classified as a solver regression",
        },
        {
            "claim_id": "no_solver_rerun",
            "claim": "This audit reads existing artifact and aggregate files only and does not call an equilibrium solver.",
            "status": "supported",
            "evidence": "audit_summary.json; command manifest",
            "boundary": "no production/reference write",
        },
    ]
    output_root.mkdir(parents=True, exist_ok=True)
    review = candidate_legacy_review(aggregate_root, aggregate_summary, output_root)
    write_csv(output_root / "aggregate_manifest_check.csv", manifest_rows[0].keys() if manifest_rows else ("path",), manifest_rows)
    write_csv(output_root / "artifact_hash_audit.csv", all_artifact_rows[0].keys() if all_artifact_rows else ("manifest",), all_artifact_rows)
    write_csv(output_root / "shard_audit.csv", shard_audits[0].keys() if shard_audits else ("run_id",), shard_audits)
    write_csv(output_root / "claim_ledger.csv", claims[0].keys(), claims)
    payload = {
        "schema_version": SCHEMA,
        "generated_utc": __import__("datetime").datetime.now(__import__("datetime").timezone.utc).isoformat(),
        "repo_head": git_head(args.repo_root.resolve()),
        "source_root": str(source_root),
        "aggregate_root": str(aggregate_root),
        "calculation_sha": CALCULATION_SHA,
        "workflow_head_sha": WORKFLOW_HEAD_SHA,
        "case_name": CASE_NAME,
        "repair_pr": 269,
        "repair_merge_sha": REPAIR_MERGE_SHA,
        "selected_unique_run_count": len(shard_rows),
        "duplicate_cancelled_run_ids": aggregate_summary.get("duplicate_cancelled_run_ids", []),
        "aggregate_manifest_ok": aggregate_manifest_ok,
        "aggregate_contract_ok": aggregate_contract_ok,
        "unexpected_artifact_failures": unexpected,
        "numerical_failures": numerical_failures,
        "known_historical_mismatch_counts": dict(Counter(row["basename"] for row in all_artifact_rows if row["verdict"] == "historical_known_mismatch")),
        "quality_reason_counts": quality_counts,
        "candidate_legacy_summary": candidate_legacy,
        "candidate_legacy_review": review,
        "verdict": verdict,
        "stop_reason": stop_reason,
        "solver_called": False,
        "production_write": False,
        "non_goals": ["no solver rerun", "no old artifact rewrite", "no RS promotion", "no legacy deletion"],
    }
    write_json(output_root / "audit_summary.json", payload)
    write_json(output_root / "command_manifest.json", {
        "command": "python scripts/analysis/relaxtime/audit_issue130_rs_sharded_provenance.py",
        "source_root": str(source_root),
        "aggregate_version": args.aggregate_version,
        "repo_head": payload["repo_head"],
        "solver_called": False,
    })
    output_files = [
        "aggregate_manifest_check.csv",
        "artifact_hash_audit.csv",
        "shard_audit.csv",
        "claim_ledger.csv",
        "candidate_legacy_numeric_extrema.csv",
        "candidate_legacy_phase_mismatches.csv",
        "candidate_legacy_review_summary.json",
        "audit_summary.json",
        "command_manifest.json",
    ]
    write_json(output_root / "manifest.json", {
        "schema_version": f"{SCHEMA}_manifest",
        "source_root": str(source_root),
        "aggregate_root": str(aggregate_root),
        "files": [
            {"path": name, "bytes": (output_root / name).stat().st_size, "sha256": sha256(output_root / name)}
            for name in output_files
        ],
    })
    print(json.dumps({"verdict": verdict, "selected_unique_run_count": len(shard_rows), "unexpected": len(unexpected), "numerical_failures": len(numerical_failures)}, ensure_ascii=False))
    return 0 if verdict == "post_repair_audit_pass_diagnostic_only" else 2


if __name__ == "__main__":
    raise SystemExit(main())
