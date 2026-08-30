#!/usr/bin/env python3
"""Collect and contract-check an approved PNJL deep-oracle anchor set."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any


EXPECTED_POINTS_BY_SET = {
    "legacy_five": {
        (-0.5, 20.0),
        (0.0, 5.0),
        (0.0, 20.0),
        (0.5, 5.0),
        (0.5, 20.0),
    },
    "required_three": {
        (-0.5, 5.0),
        (-0.5, 20.0),
        (0.0, 5.0),
    },
}
SCHEMA = "cep_deep_oracle_v1"


def _float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def _finite(value: Any) -> bool:
    return math.isfinite(_float(value))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
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


def _job_dirs(input_dir: Path) -> list[Path]:
    return sorted(path.parent for path in input_dir.rglob("job_summary.json"))


def collect(
    input_dir: Path,
    output_dir: Path,
    calculation_sha: str,
    workflow_head_sha: str,
    tag: str,
    anchor_set: str = "legacy_five",
) -> dict[str, Any]:
    if anchor_set not in EXPECTED_POINTS_BY_SET:
        raise ValueError(f"unsupported anchor set: {anchor_set}")
    expected_points = EXPECTED_POINTS_BY_SET[anchor_set]
    output_dir.mkdir(parents=True, exist_ok=True)
    errors: list[str] = []
    jobs: list[dict[str, Any]] = []
    curve_rows: list[dict[str, Any]] = []
    slice_rows: list[dict[str, Any]] = []
    cost_rows: list[dict[str, Any]] = []
    accuracy_rows: list[dict[str, Any]] = []
    seen_curve_keys: set[tuple[float, float, float]] = set()
    seen_points: set[tuple[float, float]] = set()

    directories = _job_dirs(input_dir)
    for job_dir in directories:
        summary_path = job_dir / "job_summary.json"
        summary = _read_json(summary_path)
        xi = _float(summary.get("xi"))
        temperature = _float(summary.get("temperature_MeV"))
        point = (xi, temperature)
        jobs.append({"path": str(job_dir), "summary": summary})
        if summary.get("schema_version") != SCHEMA:
            errors.append(f"wrong schema in {summary_path}")
        if summary.get("calculation_sha") != calculation_sha:
            errors.append(f"calculation SHA mismatch in {summary_path}")
        if summary.get("workflow_head_sha") != workflow_head_sha:
            errors.append(f"workflow head SHA mismatch in {summary_path}")
        if summary.get("rho_coarse_step") != 0.003125 or summary.get("rho_fine_step") != 0.0015625:
            errors.append(f"deep rho steps mismatch in {summary_path}")
        if point in seen_points:
            errors.append(f"duplicate deep-oracle point xi={xi} T={temperature}")
        seen_points.add(point)

        for name, target in (
            ("curve_points.csv", curve_rows),
            ("slice_metrics.csv", slice_rows),
            ("method_costs.csv", cost_rows),
            ("cep_accuracy.csv", accuracy_rows),
        ):
            path = job_dir / name
            if not path.is_file():
                errors.append(f"missing {name} in {job_dir}")
                continue
            declared = summary.get("curve_file_sha256", {}).get(name)
            if name == "curve_points.csv" and declared and _sha256(path) != declared:
                errors.append(f"curve hash mismatch in {job_dir}")
            with path.open(newline="", encoding="utf-8") as handle:
                target.extend(csv.DictReader(handle))

        for row in curve_rows[-999999:]:
            if row.get("xi") != str(xi) or row.get("T_MeV") != str(temperature):
                continue
            rho = _float(row.get("rho"))
            mu = _float(row.get("muq_MeV"))
            residual = _float(row.get("residual_norm"))
            key = (xi, temperature, rho)
            if key in seen_curve_keys:
                errors.append(f"duplicate curve point key xi={xi} T={temperature} rho={rho}")
            seen_curve_keys.add(key)
            if not (_finite(rho) and _finite(mu) and _finite(residual)):
                errors.append(f"non-finite curve point xi={xi} T={temperature} rho={rho}")
            if str(row.get("converged", "")).lower() not in {"true", "1"}:
                errors.append(f"non-converged curve point xi={xi} T={temperature} rho={rho}")

    if set(seen_points) != expected_points:
        missing = sorted(expected_points - seen_points)
        extra = sorted(seen_points - expected_points)
        if missing:
            errors.append(f"missing deep-oracle points: {missing}")
        if extra:
            errors.append(f"unexpected deep-oracle points: {extra}")
    if len(jobs) != len(expected_points):
        errors.append(f"expected {len(expected_points)} jobs, found {len(jobs)}")

    _write_csv(output_dir / "curve_points.csv", curve_rows)
    _write_csv(output_dir / "slice_metrics.csv", slice_rows)
    _write_csv(output_dir / "method_costs.csv", cost_rows)
    _write_csv(output_dir / "cep_accuracy.csv", accuracy_rows)
    claims = [
        {
            "claim_id": "coverage",
            "claim": f"{len(expected_points)} approved deep-oracle anchors were independently recomputed",
            "status": "pass" if not errors else "workflow_failure",
            "boundary": "diagnostic only",
        },
        {
            "claim_id": "resolution",
            "claim": "0.003125 -> 0.0015625 independent rho refinement completed",
            "status": "pass" if not errors else "workflow_failure",
            "boundary": "does not itself authorize production",
        },
        {
            "claim_id": "physics",
            "claim": "deep-oracle statuses are evidence for author review, not automatic labels",
            "status": "not_claimed",
            "boundary": "author physical review required",
        },
    ]
    _write_csv(output_dir / "claim_ledger.csv", claims)

    files = {
        name: _sha256(output_dir / name)
        for name in (
            "curve_points.csv",
            "slice_metrics.csv",
            "method_costs.csv",
            "cep_accuracy.csv",
            "claim_ledger.csv",
        )
    }
    result = {
        "schema_version": SCHEMA,
        "tag": tag,
        "anchor_set": anchor_set,
        "calculation_sha": calculation_sha,
        "workflow_head_sha": workflow_head_sha,
        "expected_points": [list(point) for point in sorted(expected_points)],
        "observed_points": [list(point) for point in sorted(seen_points)],
        "job_count": len(jobs),
        "curve_point_count": len(curve_rows),
        "errors": errors,
        "status": "complete" if not errors else "workflow_failure",
        "physical_verdict": "author_review_required" if not errors else "not_available",
        "files": files,
        "provenance": {
            "calculation_sha": calculation_sha,
            "workflow_head_sha": workflow_head_sha,
            "reference_write": False,
            "solver_scope": f"{anchor_set}_only",
        },
    }
    (output_dir / "deep_oracle_summary.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    files["deep_oracle_summary.json"] = _sha256(output_dir / "deep_oracle_summary.json")
    manifest = {
        "schema_version": SCHEMA,
        "tag": tag,
        "anchor_set": anchor_set,
        "calculation_sha": calculation_sha,
        "workflow_head_sha": workflow_head_sha,
        "source_jobs": [
            {
                "path": job["path"],
                "xi": job["summary"].get("xi"),
                "temperature_MeV": job["summary"].get("temperature_MeV"),
                "summary_sha256": _sha256(Path(job["path"]) / "job_summary.json"),
            }
            for job in jobs
        ],
        "files": files,
        "status": result["status"],
        "physical_verdict": result["physical_verdict"],
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--calculation-sha", required=True)
    parser.add_argument("--workflow-head-sha", required=True)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--anchor-set", choices=tuple(EXPECTED_POINTS_BY_SET), default="legacy_five")
    args = parser.parse_args(argv)
    result = collect(
        args.input_dir,
        args.output_dir,
        args.calculation_sha,
        args.workflow_head_sha,
        args.tag,
        args.anchor_set,
    )
    print(json.dumps(result, sort_keys=True))
    return 0 if result["status"] == "complete" else 1


if __name__ == "__main__":
    sys.exit(main())
