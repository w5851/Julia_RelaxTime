#!/usr/bin/env python3
"""Freeze and validate one Issue #130 Maxwell expansion run.

This is a solver-free evidence utility.  It copies the numerical target and
aggregate artifacts into an external, versioned evidence directory and writes
one manifest containing every copied file hash.  It never writes data/reference
and never calls Julia or the PNJL solver.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def copy_tree(source: Path, destination: Path) -> None:
    if not source.is_dir():
        raise FileNotFoundError(source)
    if destination.exists():
        raise FileExistsError(f"refusing to overwrite frozen evidence: {destination}")
    shutil.copytree(source, destination)


def validate_target(target_dir: Path, aggregate_manifest: dict[str, Any]) -> dict[str, Any]:
    manifest = read_json(target_dir / "manifest.json")
    provenance = read_json(target_dir / "provenance.json")
    summary = read_json(target_dir / "target_summary.json")
    target_id = str(manifest.get("target_id", ""))
    if not target_id:
        raise ValueError(f"empty target_id: {target_dir}")
    for name, value in (
        ("manifest", manifest),
        ("provenance", provenance),
        ("summary", summary),
    ):
        if str(value.get("target_id", target_id)) != target_id:
            raise ValueError(f"target identity mismatch in {name}: {target_id}")
    expected_calc = str(aggregate_manifest.get("calculation_sha", ""))
    expected_workflow = str(aggregate_manifest.get("postprocess_sha", ""))
    if str(provenance.get("calculation_sha", "")) != expected_calc:
        raise ValueError(f"calculation SHA mismatch: {target_id}")
    if str(provenance.get("workflow_head_sha", "")) != expected_workflow:
        raise ValueError(f"workflow SHA mismatch: {target_id}")
    for value in (manifest, provenance, summary):
        if value.get("reference_write") is not False:
            raise ValueError(f"reference_write is not false: {target_id}")
        if value.get("oracle_labels_consumed") is not False:
            raise ValueError(f"oracle labels consumed: {target_id}")
    files = manifest.get("files")
    if not isinstance(files, dict):
        raise ValueError(f"missing file manifest: {target_id}")
    for relative, expected_hash in files.items():
        path = target_dir / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        actual = sha256(path)
        if actual != str(expected_hash):
            raise ValueError(f"file hash mismatch: {target_id}/{relative}")
    return {
        "target_id": target_id,
        "xi": summary.get("xi"),
        "T_MeV": summary.get("T_MeV"),
        "verdict": summary.get("verdict"),
        "final_status": summary.get("final_status"),
        "final_geometry_converged": summary.get("final_geometry_converged"),
        "finite_and_converged": summary.get("finite_and_converged"),
        "final_candidate_count": summary.get("final_candidate_count"),
        "final_candidate_mu_MeV": summary.get("final_candidate_mu_MeV"),
        "final_area_residual": summary.get("final_area_residual"),
        "targeted_additions": summary.get("targeted_additions"),
        "curve_points": summary.get("curve_points"),
    }


def file_records(root: Path) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for path in sorted(p for p in root.rglob("*") if p.is_file()):
        records.append(
            {
                "path": path.relative_to(root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    return records


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--calculation-sha", required=True)
    parser.add_argument("--postprocess-sha", required=True)
    parser.add_argument("--tag", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.run_id.isdigit():
        raise ValueError("run id must be numeric")
    if len(args.calculation_sha) != 40 or len(args.postprocess_sha) != 40:
        raise ValueError("SHA arguments must be 40 characters")
    source = args.source_root.resolve()
    aggregate_candidates = list(source.glob("*-aggregate"))
    targets_root = source / "targets"
    if len(aggregate_candidates) != 1:
        raise ValueError(f"expected one aggregate directory under {source}")
    aggregate = aggregate_candidates[0]
    aggregate_manifest = read_json(aggregate / "manifest.json")
    if aggregate_manifest.get("expected_target_count") != 276:
        raise ValueError("unexpected aggregate target count")
    if aggregate_manifest.get("materialized_target_count") != 276:
        raise ValueError("aggregate is not complete")
    if aggregate_manifest.get("missing_target_ids") or aggregate_manifest.get("unexpected_target_ids"):
        raise ValueError("aggregate target identity mismatch")
    if aggregate_manifest.get("calculation_sha") != args.calculation_sha:
        raise ValueError("aggregate calculation SHA mismatch")
    if aggregate_manifest.get("postprocess_sha") != args.postprocess_sha:
        raise ValueError("aggregate postprocess SHA mismatch")
    target_dirs = sorted(p for p in targets_root.iterdir() if p.is_dir())
    if len(target_dirs) != 276:
        raise ValueError(f"expected 276 target directories, found {len(target_dirs)}")
    summaries = [validate_target(path, aggregate_manifest) for path in target_dirs]
    ids = [row["target_id"] for row in summaries]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate target id")
    if any(row["final_candidate_count"] != 1 for row in summaries):
        raise ValueError("not every target has a unique final candidate")
    if any(row["final_geometry_converged"] is not True for row in summaries):
        raise ValueError("not every target has a final geometry certificate")
    if any(row["finite_and_converged"] is not True for row in summaries):
        raise ValueError("not every target curve is finite/converged")

    output = args.output_root.resolve()
    output.mkdir(parents=True, exist_ok=False)
    copy_tree(aggregate, output / "aggregate")
    copy_tree(targets_root, output / "targets")
    copied_files = file_records(output)
    manifest = {
        "schema_version": "pnjl_issue130_maxwell_cep_local_expansion_freeze_v1",
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "run_id": args.run_id,
        "tag": args.tag,
        "calculation_sha": args.calculation_sha,
        "postprocess_sha": args.postprocess_sha,
        "source_workflow_sha": aggregate_manifest.get("source_workflow_sha"),
        "aggregate_schema_version": aggregate_manifest.get("schema_version"),
        "target_count": len(summaries),
        "target_ids": sorted(ids),
        "reference_write": False,
        "oracle_labels_consumed": False,
        "solver_called": True,
        "aggregate_verdict": read_json(aggregate / "verdict.json").get("verdict"),
        "integrity": {
            "aggregate_manifest_errors": aggregate_manifest.get("errors", []),
            "files": copied_files,
            "file_count": len(copied_files),
            "total_bytes": sum(int(record["bytes"]) for record in copied_files),
        },
        "target_summary_digest": hashlib.sha256(
            json.dumps(summaries, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest(),
        "boundary": "diagnostic-only; strict_reference input candidate, not phase-reference promotion",
    }
    (output / "freeze_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    (output / "README.md").write_text(
        "# Issue #130 Maxwell expansion frozen evidence\n\n"
        f"Numerical source run: `{args.run_id}`.\n\n"
        f"Calculation SHA: `{args.calculation_sha}`. Postprocess SHA: `{args.postprocess_sha}`.\n\n"
        "This directory is an immutable diagnostic evidence copy of 276 Maxwell fixed-(xi,T) targets. "
        "All target manifests, per-file hashes, finite/converged flags, unique candidates and geometry "
        "certificates were validated before copying. `reference_write=false` and `oracle_labels_consumed=false`; "
        "the copy is an input candidate for strict_reference_v1, not a promotion.\n\n"
        "The aggregate manifest, target artifacts and complete rho-mu curves are retained byte-for-byte. "
        "See `freeze_manifest.json` for the complete hash inventory.\n",
        encoding="utf-8",
    )
    print(json.dumps({"output_root": str(output), "target_count": len(summaries), "file_count": len(copied_files)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
