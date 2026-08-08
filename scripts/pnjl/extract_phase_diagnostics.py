#!/usr/bin/env python3
"""Extract a small, provenance-rich telemetry payload from phase_summary.json.

The dense-reference shard intentionally keeps the complete phase summary in the
runner workspace.  This adapter copies only stable counters and classification
diagnostics into the shard artifact so that later merges can audit solver work
without publishing the full processed tree.  Missing historical summaries are
represented explicitly as ``telemetry_unavailable``; they are never converted
to zero-valued counters.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "pnjl_phase_shard_diagnostics_v1"
NORMALIZATION_VERSION = "phase_summary_projection_v1"

COUNTER_FIELDS = (
    "scan_total",
    "scan_success",
    "scan_failure",
    "point_requests",
    "cache_hits",
    "unique_solves",
    "targeted_additions",
    "failed_points",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--process-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--stage", required=True)
    parser.add_argument("--shard-id", required=True)
    parser.add_argument("--xi", required=True, type=float)
    parser.add_argument("--calculation-sha", required=True)
    parser.add_argument("--postprocess-sha", required=True)
    parser.add_argument("--source-workflow-run-id", default="")
    return parser.parse_args()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _finite_number(value: Any) -> int | float | None:
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    return None


def _counter(value: Any) -> int | float | None:
    number = _finite_number(value)
    if number is None:
        return None
    if number < 0:
        return None
    return number


def _dict(value: Any) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def _project_summary(summary: dict[str, Any]) -> dict[str, Any]:
    stats = _dict(summary.get("stats"))
    cache = _dict(stats.get("rho_support_cache"))
    records = stats.get("sweep_records")
    records = records if isinstance(records, list) else []

    counters = {
        "scan_total": _counter(stats.get("scan_total")),
        "scan_success": _counter(stats.get("scan_success")),
        "scan_failure": _counter(stats.get("scan_failure")),
        "point_requests": _counter(cache.get("point_requests")),
        "cache_hits": _counter(cache.get("cache_hits")),
        "unique_solves": _counter(cache.get("unique_solves")),
        "targeted_additions": _counter(cache.get("targeted_additions")),
        "failed_points": _counter(cache.get("failed_points")),
    }
    missing_fields = [name for name, value in counters.items() if value is None]

    status_counts: dict[str, int] = {}
    stage_counts: dict[str, int] = {}
    certificate_counts: dict[str, int] = {}
    geometry_missing_count = 0
    geometry_normalization_versions: set[str] = set()
    for record in records:
        if not isinstance(record, dict):
            continue
        for field, target in (
            ("status", status_counts),
            ("stage_used", stage_counts),
            ("hybrid_certificate_type", certificate_counts),
        ):
            value = record.get(field)
            if value is not None:
                key = str(value)
                target[key] = target.get(key, 0) + 1
        missing = record.get("geometry_missing_fields")
        if isinstance(missing, list) and missing:
            geometry_missing_count += 1
        version = record.get("geometry_normalization_version")
        if version not in (None, ""):
            geometry_normalization_versions.add(str(version))

    return {
        "counters": counters,
        "missing_counter_fields": missing_fields,
        "status_counts": dict(sorted(status_counts.items())),
        "stage_counts": dict(sorted(stage_counts.items())),
        "certificate_counts": dict(sorted(certificate_counts.items())),
        "geometry_missing_record_count": geometry_missing_count,
        "geometry_normalization_versions": sorted(geometry_normalization_versions),
        "sweep_record_count": len(records),
        "rho_unconverged_count": _counter(stats.get("rho_unconverged_count")),
        "temperature_unconverged_count": _counter(stats.get("temperature_unconverged_count")),
    }


def build_payload(args: argparse.Namespace) -> dict[str, Any]:
    summaries = sorted(args.process_root.rglob("phase_summary.json")) if args.process_root.is_dir() else []
    base: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "normalization_version": NORMALIZATION_VERSION,
        "tag": args.tag,
        "stage": args.stage,
        "shard_id": args.shard_id,
        "xi": args.xi,
        "calculation_sha": args.calculation_sha,
        "postprocess_sha": args.postprocess_sha,
        "source_workflow_run_id": args.source_workflow_run_id,
        "phase_summary_paths": [path.as_posix() for path in summaries],
        "phase_summary_sha256": {},
        "config_hash": None,
        "summary_schema_version": None,
        "availability": "telemetry_unavailable",
        "unavailable_reason": "phase_summary_missing",
        "missing_fields": list(COUNTER_FIELDS),
        "diagnostics": None,
    }
    if len(summaries) != 1:
        if len(summaries) > 1:
            base["unavailable_reason"] = "multiple_phase_summaries"
        return base

    summary_path = summaries[0]
    base["phase_summary_sha256"] = {summary_path.as_posix(): _sha256(summary_path)}
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        base["unavailable_reason"] = "phase_summary_invalid_json"
        return base
    if not isinstance(summary, dict):
        base["unavailable_reason"] = "phase_summary_not_object"
        return base

    projected = _project_summary(summary)
    base.update(
        {
            "availability": "available",
            "unavailable_reason": None,
            "config_hash": _dict(summary.get("config_snapshot")).get("config_hash")
            or summary.get("config_hash"),
            "summary_schema_version": summary.get("schema_version"),
            "missing_fields": projected["missing_counter_fields"],
            "diagnostics": projected,
        }
    )
    return base


def main() -> None:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    payload = build_payload(args)
    args.output.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"availability": payload["availability"], "output": args.output.as_posix()}))


if __name__ == "__main__":
    main()
