#!/usr/bin/env python3
"""Validate the one-anchor Maxwell endpoint refinement artifact."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


EXPECTED_SCHEMA = "pnjl_maxwell_endpoint_refinement_v1"
EXPECTED_XI = -0.5
EXPECTED_T = 5.0
MAX_TARGETED = 12


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--calculation-sha", required=True)
    parser.add_argument("--workflow-head-sha", required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((args.input_dir / "manifest.json").read_text(encoding="utf-8"))
    errors: list[str] = []
    if manifest.get("schema_version") != EXPECTED_SCHEMA:
        errors.append("manifest schema mismatch")
    if manifest.get("calculation_sha") != args.calculation_sha:
        errors.append("manifest calculation SHA mismatch")
    if manifest.get("workflow_head_sha") != args.workflow_head_sha:
        errors.append("manifest workflow SHA mismatch")
    curve = rows(args.input_dir / "curve_points.csv")
    metrics = rows(args.input_dir / "slice_metrics.csv")
    frontier = rows(args.input_dir / "policy_frontier.csv")
    if not curve or not metrics or not frontier:
        errors.append("missing required non-empty CSV")
    keys = {(row["xi"], row["T_MeV"], row["rho"]) for row in curve}
    if len(keys) != len(curve):
        errors.append("duplicate curve key")
    for row in curve:
        if abs(float(row["xi"]) - EXPECTED_XI) > 1e-9 or abs(float(row["T_MeV"]) - EXPECTED_T) > 1e-9:
            errors.append("unexpected anchor")
        if row["calculation_sha"] != args.calculation_sha or row["workflow_head_sha"] != args.workflow_head_sha:
            errors.append("curve row provenance mismatch")
        if row["converged"].lower() != "true" or row["finite"].lower() != "true":
            errors.append("non-finite/non-converged curve row")
    if max(int(row["targeted_additions"]) for row in metrics) > MAX_TARGETED:
        errors.append("targeted cap exceeded")
    selected = json.loads((args.input_dir / "selected_policy.json").read_text(encoding="utf-8"))
    verdict = str(manifest.get("verdict", "workflow_failure"))
    if selected.get("verdict") != verdict:
        errors.append("selected policy verdict mismatch")
    aggregate = {
        "schema_version": EXPECTED_SCHEMA,
        "verdict": verdict if not errors else "workflow_failure",
        "workflow_contract_errors": errors,
        "calculation_sha": args.calculation_sha,
        "workflow_head_sha": args.workflow_head_sha,
        "curve_rows": len(curve),
        "metric_rows": len(metrics),
        "frontier_rows": len(frontier),
        "targeted_additions": max((int(row["targeted_additions"]) for row in metrics), default=0),
        "reference_write": False,
        "files": {path.name: sha256(path) for path in args.input_dir.iterdir() if path.is_file()},
    }
    (args.output_dir / "endpoint_refinement_summary.json").write_text(
        json.dumps(aggregate, indent=2) + "\n", encoding="utf-8"
    )
    (args.output_dir / "manifest.json").write_text(
        json.dumps(aggregate, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(aggregate))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
