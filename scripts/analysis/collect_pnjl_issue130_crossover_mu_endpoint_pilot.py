#!/usr/bin/env python3
"""Validate and aggregate Issue #130 crossover-mu endpoint pilot artifacts."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "pnjl_issue130_crossover_mu_endpoint_pilot_v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def target_index(path: Path) -> tuple[dict[str, dict[str, str]], str]:
    rows = read_csv(path)
    selected = {
        row["target_id"]: row
        for row in rows
        if row.get("pilot_selection") == "pilot_candidate"
    }
    if len(selected) != 8:
        raise ValueError(f"expected exactly 8 pilot targets, found {len(selected)}")
    if len(selected) != sum(row.get("pilot_selection") == "pilot_candidate" for row in rows):
        raise ValueError("duplicate pilot target_id")
    return selected, sha256(path)


def bool_value(value: Any) -> bool:
    return value is True or str(value).strip().lower() in {"true", "1", "yes"}


def aggregate(args: argparse.Namespace) -> int:
    target_path = Path(args.target_list).resolve()
    expected, target_hash = target_index(target_path)
    root = Path(args.input_dir).resolve()
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)

    summary_paths = sorted(root.rglob("target_summary.json"))
    errors: list[str] = []
    summaries: dict[str, dict[str, Any]] = {}
    response_rows: list[dict[str, Any]] = []
    artifact_hashes: dict[str, str] = {}

    for summary_path in summary_paths:
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            errors.append(f"invalid summary {summary_path}: {exc}")
            continue
        target_id = str(summary.get("target_id", ""))
        if target_id not in expected:
            errors.append(f"unexpected target_id={target_id}")
            continue
        if target_id in summaries:
            errors.append(f"duplicate artifact target_id={target_id}")
            continue
        if summary.get("schema_version") != SCHEMA_VERSION:
            errors.append(f"{target_id}: schema mismatch")
        if str(summary.get("calculation_sha", "")).lower() != args.calculation_sha.lower():
            errors.append(f"{target_id}: calculation SHA mismatch")
        if summary.get("target_list_sha256") != target_hash:
            errors.append(f"{target_id}: target list hash mismatch")
        target = summary.get("target") or {}
        if target.get("physical_side") != "crossover_mu_lt_CEP_proxy":
            errors.append(f"{target_id}: invalid physical side")
        solver = summary.get("solver") or {}
        if not bool_value(solver.get("solver_called")):
            errors.append(f"{target_id}: solver_called is false")
        provenance_path = summary_path.with_name("provenance.json")
        if not provenance_path.is_file():
            errors.append(f"{target_id}: provenance.json missing")
        else:
            provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
            if bool_value(provenance.get("reference_write")):
                errors.append(f"{target_id}: reference_write must be false")
            if str(provenance.get("calculation_sha", "")).lower() != args.calculation_sha.lower():
                errors.append(f"{target_id}: provenance calculation SHA mismatch")
            artifact_hashes[f"{target_id}/target_summary.json"] = sha256(summary_path)
            artifact_hashes[f"{target_id}/provenance.json"] = sha256(provenance_path)

        response_path = summary_path.with_name("temperature_response.csv")
        if not response_path.is_file():
            errors.append(f"{target_id}: temperature_response.csv missing")
        else:
            rows = read_csv(response_path)
            if not rows:
                errors.append(f"{target_id}: empty temperature response")
            for row in rows:
                response_rows.append({"target_id": target_id, **row})
                if not bool_value(row.get("finite")):
                    errors.append(f"{target_id}: non-finite supplemental response row")
            artifact_hashes[f"{target_id}/temperature_response.csv"] = sha256(response_path)
        summaries[target_id] = summary

    missing = sorted(set(expected) - set(summaries))
    errors.extend(f"missing target_id={target_id}" for target_id in missing)

    summary_rows: list[dict[str, Any]] = []
    for target_id in sorted(expected):
        summary = summaries.get(target_id)
        if summary is None:
            continue
        target = summary.get("target") or {}
        result = summary.get("result") or {}
        summary_rows.append(
            {
                "target_id": target_id,
                "xi": target.get("xi"),
                "target_mu_MeV": target.get("target_mu_MeV"),
                "mu_CEP_proxy_MeV": target.get("mu_CEP_proxy_MeV"),
                "T_crossover_MeV": result.get("T_crossover_MeV"),
                "rho_crossover": result.get("rho_crossover"),
                "derivative_peak": result.get("derivative_peak"),
                "found": result.get("found"),
                "status": result.get("status"),
                "finite_and_converged": summary.get("finite_and_converged"),
                "runner_seconds": summary.get("runner_seconds"),
                "estimated_detector_calls": (summary.get("solver") or {}).get("estimated_detector_calls"),
            }
        )

    if errors:
        verdict = "pilot_artifact_invalid" if any("missing" in item or "mismatch" in item or "schema" in item for item in errors) else "pilot_solver_or_curve_failure"
    elif len(summary_rows) != len(expected):
        verdict = "pilot_artifact_invalid"
    elif not all(bool_value(row["finite_and_converged"]) for row in summary_rows):
        verdict = "pilot_solver_or_curve_failure"
    elif not all(bool_value(row["found"]) and row["status"] == "crossover_candidate" for row in summary_rows):
        verdict = "pilot_endpoint_inconclusive"
    else:
        verdict = "pilot_candidate"

    write_csv(
        output / "pilot_summary.csv",
        summary_rows,
        [
            "target_id", "xi", "target_mu_MeV", "mu_CEP_proxy_MeV",
            "T_crossover_MeV", "rho_crossover", "derivative_peak", "found",
            "status", "finite_and_converged", "runner_seconds", "estimated_detector_calls",
        ],
    )
    write_csv(
        output / "temperature_response.csv",
        response_rows,
        ["target_id", "T_MeV", "derivative", "rho", "finite", "solver_ok", "error", "scan_derivative"],
    )

    claims = [
        {
            "claim_id": "pilot_targets_are_frozen",
            "claim": "The numerical pilot consumed only the eight representative target IDs from the immutable preflight list.",
            "status": "supported" if not missing else "failed",
            "evidence": "pilot_summary.csv; manifest.json",
            "boundary": "does not justify expansion to all xi",
        },
        {
            "claim_id": "pilot_does_not_write_reference",
            "claim": "The pilot artifacts record reference_write=false and do not replace C2 or phase-reference evidence.",
            "status": "supported" if not errors else "author_check",
            "evidence": "provenance.json; manifest.json",
            "boundary": "diagnostic-only",
        },
        {
            "claim_id": "pilot_verdict",
            "claim": f"Aggregate pilot verdict is {verdict}.",
            "status": "supported",
            "evidence": "verdict.json; pilot_summary.csv",
            "boundary": "not a phase-reference promotion decision",
        },
    ]
    write_csv(output / "claim_ledger.csv", claims, ["claim_id", "claim", "status", "evidence", "boundary"])

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "run_mode": args.run_mode,
        "source_run_id": args.source_run_id or None,
        "calculation_sha": args.calculation_sha.lower(),
        "postprocess_sha": args.postprocess_sha,
        "target_list": str(target_path),
        "target_list_sha256": target_hash,
        "expected_target_count": len(expected),
        "materialized_target_count": len(summary_rows),
        "missing_target_ids": missing,
        "solver_called": bool(summary_rows),
        "reference_write": False,
        "oracle_labels_consumed": False,
        "artifact_hashes": artifact_hashes,
        "errors": errors,
        "generated_at": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "verdict": verdict,
    }
    write_json(output / "manifest.json", manifest)
    write_json(output / "verdict.json", {"verdict": verdict, "errors": errors, "stop_expansion": verdict != "pilot_candidate"})
    (output / "AUDIT.md").write_text(
        "# Issue #130 crossover mu endpoint pilot\n\n"
        f"- verdict: `{verdict}`\n"
        f"- calculation SHA: `{args.calculation_sha}`\n"
        f"- expected targets: `{len(expected)}`\n"
        f"- materialized targets: `{len(summary_rows)}`\n"
        "- oracle labels consumed: `false`\n"
        "- reference write: `false`\n\n"
        "This is diagnostic-only evidence. A `pilot_candidate` permits a separately authorized all-xi expansion; it does not promote phase-reference.\n",
        encoding="utf-8",
    )
    print(json.dumps({"verdict": verdict, "errors": errors, "target_count": len(summary_rows)}, sort_keys=True))
    return 0 if verdict == "pilot_candidate" else 2


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--target-list", required=True)
    parser.add_argument("--calculation-sha", required=True)
    parser.add_argument("--postprocess-sha", required=True)
    parser.add_argument("--run-mode", choices=("numerical", "aggregate_replay"), default="numerical")
    parser.add_argument("--source-run-id", default="")
    return aggregate(parser.parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
