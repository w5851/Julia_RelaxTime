#!/usr/bin/env python3
"""Collect and compare the paired Issue #130 RS numerical pilot.

The collector is intentionally diagnostic-only. It validates that the runtime
and explicit legacy rollback completed the same planned points, records hard
failures/non-finite/duplicate rows, and reports transport drift without
turning an arbitrary drift threshold into a production gate.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any, Iterable


REFERENCE_MODES = ("candidate_runtime", "legacy")
# The phase-reference choice can move the anchor temperature.  Pair rows by
# the requested point identity and report T/phase changes separately.
KEY_FIELDS = ("muB_MeV", "xi", "mode", "alpha_T")
COMPARE_FIELDS = (
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
REQUIRED_FIELDS = KEY_FIELDS + (
    "T_MeV",
    "converged",
    "quality_flag",
    "quality_reason",
    "phase_reference_kind",
    "phase_structure",
) + COMPARE_FIELDS


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--repo-sha", required=True)
    parser.add_argument("--calculation-sha", required=True)
    parser.add_argument("--run-id", default=os.environ.get("GITHUB_RUN_ID", ""))
    # 18 requested points plus one deterministic coexistence-side replacement
    # created by the existing direct_coexistence anchor route.
    parser.add_argument("--expected-points", type=int, default=19)
    parser.add_argument("--fail-on-hard-failure", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def read_data_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        return [], []
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line.strip() and not line.startswith("#")]
    if not lines:
        return [], []
    reader = csv.DictReader(lines)
    rows = list(reader)
    return list(reader.fieldnames or []), rows


def read_effective_config(path: Path) -> dict[str, Any]:
    payload = read_json(path)
    options = payload.get("options")
    return options if isinstance(options, dict) else payload


def bool_value(value: str | None) -> bool | None:
    if value is None:
        return None
    normalized = value.strip().lower()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no"}:
        return False
    return None


def as_float(value: str | None) -> float | None:
    if value is None or value.strip() == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def finite_number(value: str | None) -> bool:
    number = as_float(value)
    return number is not None and math.isfinite(number)


def key(row: dict[str, str]) -> tuple[str, ...]:
    return tuple(row.get(field, "") for field in KEY_FIELDS)


def duplicate_keys(rows: Iterable[dict[str, str]]) -> list[tuple[str, ...]]:
    seen: set[tuple[str, ...]] = set()
    duplicates: list[tuple[str, ...]] = []
    for row in rows:
        item = key(row)
        if item in seen:
            duplicates.append(item)
        seen.add(item)
    return duplicates


def mode_summary(root: Path, reference_mode: str, expected_points: int) -> dict[str, Any]:
    result = root / "results" / reference_mode / "phase_guided_transport_scan.csv"
    plan = root / "results" / reference_mode / "sampling_plan.csv"
    failed = root / "results" / reference_mode / "failed_points.csv"
    effective = root / "results" / reference_mode / "effective_config.json"
    manifest = root / "results" / reference_mode / "run_manifest.json"
    pilot_metadata = root / "results" / reference_mode / "pilot_metadata.json"
    status_path = root / "status" / f"{reference_mode}.txt"
    header, rows = read_data_rows(result)
    plan_header, plan_rows = read_data_rows(plan)
    missing_fields = [field for field in REQUIRED_FIELDS if field not in header]
    duplicate = duplicate_keys(rows)
    plan_duplicate = duplicate_keys(plan_rows)
    nonfinite: list[tuple[str, ...]] = []
    negative: list[tuple[str, ...]] = []
    nonconverged: list[tuple[str, ...]] = []
    invalid_quality: list[tuple[str, ...]] = []
    for row in rows:
        if any(not finite_number(row.get(field)) for field in COMPARE_FIELDS):
            nonfinite.append(key(row))
        if any((number := as_float(row.get(field))) is not None and number < 0 for field in COMPARE_FIELDS):
            negative.append(key(row))
        if bool_value(row.get("converged")) is not True:
            nonconverged.append(key(row))
        # scan_quality's flag is an issue marker: false means the point passed
        # the transport quality checks.
        if bool_value(row.get("quality_flag")) is not False:
            invalid_quality.append(key(row))
    _, failed_rows = read_data_rows(failed)
    status = status_path.read_text(encoding="utf-8").strip() if status_path.is_file() else "missing"
    config = read_effective_config(effective)
    run_manifest = read_json(manifest)
    metadata = read_json(pilot_metadata)
    return {
        "reference_mode": reference_mode,
        "result_path": str(result).replace("\\", "/"),
        "result_sha256": sha256(result) if result.is_file() else "",
        "effective_config_sha256": sha256(effective) if effective.is_file() else "",
        "run_manifest_sha256": sha256(manifest) if manifest.is_file() else "",
        "pilot_metadata_sha256": sha256(pilot_metadata) if pilot_metadata.is_file() else "",
        "status": status,
        "exit_success": status == "0",
        "row_count": len(rows),
        "expected_points": expected_points,
        "point_count_matches": len(rows) == expected_points,
        "plan_row_count": len(plan_rows),
        "plan_point_count_matches": len(plan_rows) == expected_points,
        "failed_point_count": len(failed_rows),
        "missing_fields": missing_fields,
        "duplicate_key_count": len(duplicate),
        "duplicate_keys": [list(item) for item in duplicate[:10]],
        "plan_duplicate_key_count": len(plan_duplicate),
        "plan_duplicate_keys": [list(item) for item in plan_duplicate[:10]],
        "nonfinite_count": len(nonfinite),
        "nonfinite_keys": [list(item) for item in nonfinite[:10]],
        "negative_transport_count": len(negative),
        "negative_transport_keys": [list(item) for item in negative[:10]],
        "nonconverged_count": len(nonconverged),
        "nonconverged_keys": [list(item) for item in nonconverged[:10]],
        "invalid_quality_count": len(invalid_quality),
        "invalid_quality_keys": [list(item) for item in invalid_quality[:10]],
        "solver_called": metadata.get("solver_called") is True,
        "pilot_metadata": metadata,
        "phase_reference_mode": config.get("phase_reference_mode", ""),
        "phase_reference_source": config.get("phase_reference_source", ""),
        "phase_reference_summary": config.get("phase_reference_summary", {}),
        "run_id": run_manifest.get("run_id", ""),
        "rows": rows,
        "plan_rows": plan_rows,
    }


def write_csv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def compare(runtime: dict[str, Any], legacy: dict[str, Any]) -> list[dict[str, Any]]:
    left = {key(row): row for row in runtime["rows"]}
    right = {key(row): row for row in legacy["rows"]}
    rows: list[dict[str, Any]] = []
    for item in sorted(set(left) & set(right)):
        left_row = left[item]
        right_row = right[item]
        base = dict(zip(KEY_FIELDS, item))
        for field in ("T_MeV", "muq_MeV", "phase_reference_kind", "phase_structure", "quality_reason"):
            base[f"{field}_candidate"] = left_row.get(field, "")
            base[f"{field}_legacy"] = right_row.get(field, "")
        for field in COMPARE_FIELDS:
            a = as_float(left_row.get(field))
            b = as_float(right_row.get(field))
            absolute = None if a is None or b is None else abs(a - b)
            denominator = max(abs(b), 1e-30) if b is not None else None
            relative = None if absolute is None or denominator is None else absolute / denominator
            base.update({f"{field}_candidate": a, f"{field}_legacy": b, f"{field}_abs_diff": absolute, f"{field}_rel_diff": relative})
        rows.append(base)
    return rows


def main() -> int:
    args = parse_args()
    output = args.output.resolve()
    runtime = mode_summary(args.pilot_root, REFERENCE_MODES[0], args.expected_points)
    legacy = mode_summary(args.pilot_root, REFERENCE_MODES[1], args.expected_points)
    comparison = compare(runtime, legacy)
    hard_failures = []
    for summary in (runtime, legacy):
        label = summary["reference_mode"]
        if not summary["exit_success"]:
            hard_failures.append(f"{label}:runner_exit_{summary['status']}")
        if not summary["point_count_matches"]:
            hard_failures.append(f"{label}:point_count_{summary['row_count']}/{summary['expected_points']}")
        if not summary["plan_point_count_matches"]:
            hard_failures.append(f"{label}:plan_point_count_{summary['plan_row_count']}/{summary['expected_points']}")
        if summary["missing_fields"]:
            hard_failures.append(f"{label}:missing_fields")
        if not summary["solver_called"]:
            hard_failures.append(f"{label}:solver_called_not_confirmed")
        if summary["failed_point_count"]:
            hard_failures.append(f"{label}:failed_points_{summary['failed_point_count']}")
        if summary["duplicate_key_count"]:
            hard_failures.append(f"{label}:duplicate_keys_{summary['duplicate_key_count']}")
        if summary["plan_duplicate_key_count"]:
            hard_failures.append(f"{label}:plan_duplicate_keys_{summary['plan_duplicate_key_count']}")
        if summary["nonfinite_count"]:
            hard_failures.append(f"{label}:nonfinite_{summary['nonfinite_count']}")
        if summary["negative_transport_count"]:
            hard_failures.append(f"{label}:negative_transport_{summary['negative_transport_count']}")
        if summary["nonconverged_count"]:
            hard_failures.append(f"{label}:nonconverged_{summary['nonconverged_count']}")
        if summary["invalid_quality_count"]:
            hard_failures.append(f"{label}:invalid_quality_{summary['invalid_quality_count']}")
        metadata = summary["pilot_metadata"]
        if metadata.get("reference_mode") != label:
            hard_failures.append(f"{label}:pilot_metadata_reference_mode")
        if metadata.get("calculation_sha") != args.calculation_sha:
            hard_failures.append(f"{label}:pilot_metadata_calculation_sha")
        expected_phase_mode = "runtime" if label == "candidate_runtime" else "legacy"
        expected_source = "candidate" if label == "candidate_runtime" else "legacy"
        if summary["phase_reference_mode"] != expected_phase_mode:
            hard_failures.append(f"{label}:effective_config_phase_reference_mode")
        if summary["phase_reference_source"] != expected_source:
            hard_failures.append(f"{label}:effective_config_phase_reference_source")
    plan_keys = [{key(row) for row in summary["plan_rows"]} for summary in (runtime, legacy)]
    if plan_keys[0] != plan_keys[1]:
        hard_failures.append("candidate_legacy_plan_key_mismatch")
    result_keys = [{key(row) for row in summary["rows"]} for summary in (runtime, legacy)]
    if result_keys[0] != result_keys[1]:
        hard_failures.append("candidate_legacy_result_key_mismatch")
    if not comparison:
        hard_failures.append("no_common_transport_rows")
    verdict = "pilot_pair_complete_diagnostic_only" if not hard_failures else "pilot_solver_or_curve_failure"
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / "transport_comparison.csv", list(comparison[0].keys()) if comparison else list(KEY_FIELDS), comparison)
    write_csv(
        output / "reference_summary.csv",
        ["reference_mode", "status", "row_count", "expected_points", "plan_row_count", "failed_point_count", "duplicate_key_count", "plan_duplicate_key_count", "nonfinite_count", "negative_transport_count", "nonconverged_count", "invalid_quality_count", "phase_reference_mode", "phase_reference_source", "solver_called"],
        [{field: summary[field] for field in ("reference_mode", "status", "row_count", "expected_points", "plan_row_count", "failed_point_count", "duplicate_key_count", "plan_duplicate_key_count", "nonfinite_count", "negative_transport_count", "nonconverged_count", "invalid_quality_count", "phase_reference_mode", "phase_reference_source", "solver_called")} for summary in (runtime, legacy)],
    )
    manifest = {
        "schema_version": "pnjl_issue130_rs_numerical_pilot_v1",
        "verdict": verdict,
        "diagnostic_only": True,
        "repo_sha": args.repo_sha,
        "calculation_sha": args.calculation_sha,
        "run_id": args.run_id,
        "solver_called": True,
        "nominal_requested_points": 18,
        "effective_planned_points": args.expected_points,
        "plan_expansion_reason": "direct_coexistence replaces xi=0 with certified +/- coexistence-side points for one first-order anchor" if args.expected_points != 18 else "none",
        "reference_modes": {summary["reference_mode"]: {key: value for key, value in summary.items() if key not in {"rows", "plan_rows"}} for summary in (runtime, legacy)},
        "common_row_count": len(comparison),
        "hard_failures": hard_failures,
        "non_goals": ["no production/reference write", "no legacy deletion", "no automatic numerical drift acceptance"],
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (output / "verdict.json").write_text(json.dumps({"verdict": verdict, "hard_failures": hard_failures}, indent=2) + "\n", encoding="utf-8")
    if args.fail_on_hard_failure and hard_failures:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
