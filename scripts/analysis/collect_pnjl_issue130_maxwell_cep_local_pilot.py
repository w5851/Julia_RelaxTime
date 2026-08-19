#!/usr/bin/env python3
"""Validate and aggregate the Issue #130 Maxwell CEP-local pilot."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


SCHEMA_VERSION = "pnjl_issue130_maxwell_cep_local_pilot_v1"
RUNNER_SCHEMA = SCHEMA_VERSION
MAX_TARGETED = 12
MATERIALIZATION_CONTRACT_VERSION = "aggregate_identity_fallback_v1"
IDENTITY_FIELDS = ("calculation_sha", "workflow_head_sha")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def bool_value(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def finite(value: object) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def text_value(value: object) -> str:
    return "" if value is None else str(value).strip()


def target_rows(path: Path) -> dict[str, dict[str, str]]:
    rows = read_rows(path)
    selected = {
        row["target_id"]: row
        for row in rows
        if row.get("pilot_selection") == "pilot_candidate"
    }
    if len(selected) != 11:
        raise ValueError(f"expected exactly 11 pilot targets, found {len(selected)}")
    if any(row.get("target_kind") != "maxwell_fixed_xi_T" for row in selected.values()):
        raise ValueError("pilot list contains a non-fixed-(xi,T) target")
    return selected


def validate_identity(
    target_id: str,
    summary: dict,
    provenance: dict,
    manifest: dict,
    calculation_sha: str,
    workflow_head_sha: str,
) -> tuple[list[str], list[str]]:
    """Validate artifact identity and record compatible summary omissions.

    Early v1 numerical artifacts put immutable identities in provenance and
    manifest but omitted them from target_summary. Those two files remain the
    canonical sources; a missing summary field is accepted only when both agree
    with the requested inputs. A present but wrong summary value still fails.
    """

    errors: list[str] = []
    fallback_fields: list[str] = []
    expected_values = {
        "calculation_sha": calculation_sha,
        "workflow_head_sha": workflow_head_sha,
    }
    for field in IDENTITY_FIELDS:
        expected = expected_values[field].lower()
        source_values = {
            "provenance": text_value(provenance.get(field, "")),
            "manifest": text_value(manifest.get(field, "")),
        }
        for label, actual in source_values.items():
            if actual.lower() != expected:
                errors.append(f"{target_id}: {label} {field} mismatch")
        if source_values["provenance"].lower() != source_values["manifest"].lower():
            errors.append(f"{target_id}: provenance/manifest {field} mismatch")

        summary_value = text_value(summary.get(field, ""))
        sources_verified = all(
            source_values[label].lower() == expected
            for label in ("provenance", "manifest")
        )
        if not summary_value:
            if sources_verified:
                fallback_fields.append(field)
            else:
                errors.append(
                    f"{target_id}: summary {field} missing and verified identity unavailable"
                )
        elif summary_value.lower() != expected:
            errors.append(f"{target_id}: summary {field} mismatch")
    return errors, fallback_fields


def artifact_dirs(input_dir: Path) -> dict[str, Path]:
    found: dict[str, Path] = {}
    for summary_path in input_dir.rglob("target_summary.json"):
        summary = read_json(summary_path)
        target_id = str(summary.get("target_id", ""))
        if target_id:
            if target_id in found:
                raise ValueError(f"duplicate target artifact: {target_id}")
            found[target_id] = summary_path.parent
    return found


def validate_target(
    target_id: str,
    expected: dict[str, str],
    directory: Path,
    calculation_sha: str,
    workflow_head_sha: str,
) -> tuple[dict[str, object], dict[str, str], list[str], dict[str, str], list[str]]:
    errors: list[str] = []
    required = ["target_summary.json", "provenance.json", "manifest.json", "curve_points.csv", "slice_metrics.csv", "policy_frontier.csv", "method_costs.csv"]
    for name in required:
        if not (directory / name).is_file():
            errors.append(f"{target_id}: missing {name}")
    if errors:
        return {}, {}, errors, {}, []

    summary = read_json(directory / "target_summary.json")
    provenance = read_json(directory / "provenance.json")
    manifest = read_json(directory / "manifest.json")
    if summary.get("schema_version") != RUNNER_SCHEMA:
        errors.append(f"{target_id}: summary schema mismatch")
    if manifest.get("schema_version") != RUNNER_SCHEMA:
        errors.append(f"{target_id}: manifest schema mismatch")
    if summary.get("target_id") != target_id or provenance.get("target_id") != target_id:
        errors.append(f"{target_id}: target id mismatch")
    identity_errors, identity_fallback_fields = validate_identity(
        target_id, summary, provenance, manifest, calculation_sha, workflow_head_sha
    )
    errors.extend(identity_errors)
    if not bool_value(provenance.get("solver_called")):
        errors.append(f"{target_id}: solver_called must be true for numerical pilot")
    if bool_value(provenance.get("reference_write")) or bool_value(manifest.get("reference_write")):
        errors.append(f"{target_id}: reference_write must be false")
    if bool_value(provenance.get("oracle_labels_consumed")) or bool_value(manifest.get("oracle_labels_consumed")):
        errors.append(f"{target_id}: oracle labels must not be consumed")

    try:
        if abs(float(expected["xi"]) - float(summary["xi"])) > 1e-12:
            errors.append(f"{target_id}: xi mismatch")
        if abs(float(expected["T_MeV"]) - float(summary["T_MeV"])) > 1e-12:
            errors.append(f"{target_id}: temperature mismatch")
    except (KeyError, ValueError):
        errors.append(f"{target_id}: target coordinates are invalid")

    curve = read_rows(directory / "curve_points.csv")
    metrics = read_rows(directory / "slice_metrics.csv")
    frontier = read_rows(directory / "policy_frontier.csv")
    keys = {(row.get("target_id"), row.get("xi"), row.get("T_MeV"), row.get("rho")) for row in curve}
    if len(keys) != len(curve):
        errors.append(f"{target_id}: duplicate rho curve key")
    if not curve:
        errors.append(f"{target_id}: empty rho curve")
    for row in curve:
        if row.get("target_id") != target_id:
            errors.append(f"{target_id}: curve target mismatch")
        if row.get("converged", "").lower() != "true" or row.get("finite", "").lower() != "true":
            errors.append(f"{target_id}: non-finite/non-converged curve row")
        for field in ("xi", "T_MeV", "rho", "muq_MeV", "residual_norm"):
            if not finite(row.get(field)):
                errors.append(f"{target_id}: non-finite curve field {field}")
    if not metrics:
        errors.append(f"{target_id}: empty slice metrics")
    if metrics:
        try:
            max_targeted = max(int(row.get("targeted_additions", "-1")) for row in metrics)
            if max_targeted > MAX_TARGETED:
                errors.append(f"{target_id}: targeted cap exceeded")
        except ValueError:
            errors.append(f"{target_id}: invalid targeted_additions")
    final_status = str(summary.get("final_status", ""))
    try:
        candidate_count = int(summary.get("final_candidate_count", -1))
    except (TypeError, ValueError):
        candidate_count = -1
    geometry = bool_value(summary.get("final_geometry_converged"))
    finite_converged = bool_value(summary.get("finite_and_converged"))
    summary_row = {
        "target_id": target_id,
        "xi": expected["xi"],
        "T_MeV": expected["T_MeV"],
        "preflight_reason": expected.get("reason", ""),
        "verdict": summary.get("verdict", ""),
        "final_status": final_status,
        "final_reason": summary.get("final_reason", ""),
        "candidate_count": candidate_count,
        "final_candidate_mu_MeV": summary.get("final_candidate_mu_MeV", ""),
        "final_area_residual": summary.get("final_area_residual", ""),
        "final_geometry_converged": geometry,
        "finite_and_converged": finite_converged,
        "curve_points": len(curve),
        "metric_rows": len(metrics),
        "targeted_additions": summary.get("targeted_additions", ""),
        "identity_fallback_fields": ";".join(identity_fallback_fields),
    }
    curve_index = {
        "target_id": target_id,
        "curve_points": str(len(curve)),
        "curve_sha256": sha256(directory / "curve_points.csv"),
        "metrics_sha256": sha256(directory / "slice_metrics.csv"),
        "method_costs_sha256": sha256(directory / "method_costs.csv"),
    }
    hashes = {
        f"{target_id}/target_summary.json": sha256(directory / "target_summary.json"),
        f"{target_id}/provenance.json": sha256(directory / "provenance.json"),
        f"{target_id}/manifest.json": sha256(directory / "manifest.json"),
        f"{target_id}/curve_points.csv": curve_index["curve_sha256"],
        f"{target_id}/slice_metrics.csv": curve_index["metrics_sha256"],
        f"{target_id}/method_costs.csv": curve_index["method_costs_sha256"],
    }
    return summary_row, curve_index, errors, hashes, identity_fallback_fields


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--target-list", type=Path, required=True)
    parser.add_argument("--calculation-sha", required=True)
    parser.add_argument("--postprocess-sha", required=True)
    parser.add_argument("--run-mode", choices=("numerical", "aggregate_replay"), default="numerical")
    parser.add_argument("--source-run-id", default="")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    expected = target_rows(args.target_list)
    directories = artifact_dirs(args.input_dir)
    errors: list[str] = []
    summaries: list[dict[str, object]] = []
    curve_indices: list[dict[str, str]] = []
    hashes: dict[str, str] = {}
    identity_fallbacks: dict[str, list[str]] = {}
    for target_id, expected_row in expected.items():
        directory = directories.get(target_id)
        if directory is None:
            errors.append(f"{target_id}: numerical artifact missing")
            continue
        summary, curve_index, target_errors, target_hashes, fallback_fields = validate_target(
            target_id, expected_row, directory, args.calculation_sha, args.postprocess_sha
        )
        errors.extend(target_errors)
        if summary:
            summaries.append(summary)
        if curve_index:
            curve_indices.append(curve_index)
        hashes.update(target_hashes)
        if fallback_fields:
            identity_fallbacks[target_id] = fallback_fields
    unexpected = sorted(set(directories) - set(expected))
    errors.extend(f"unexpected target artifact: {target_id}" for target_id in unexpected)
    summaries.sort(key=lambda row: (float(row["xi"]), float(row["T_MeV"])))
    curve_indices.sort(key=lambda row: row["target_id"])
    write_csv(args.output_dir / "pilot_summary.csv", summaries, [
        "target_id", "xi", "T_MeV", "preflight_reason", "verdict", "final_status",
        "final_reason", "candidate_count", "final_candidate_mu_MeV", "final_area_residual",
        "final_geometry_converged", "finite_and_converged", "curve_points", "metric_rows",
        "targeted_additions", "identity_fallback_fields",
    ])
    write_csv(args.output_dir / "curve_index.csv", curve_indices, [
        "target_id", "curve_points", "curve_sha256", "metrics_sha256", "method_costs_sha256"
    ])
    all_valid = len(summaries) == len(expected) and not errors
    all_feasible = all(
        row["final_status"] == "first_order"
        and int(row["candidate_count"]) == 1
        and bool_value(row["final_geometry_converged"])
        and bool_value(row["finite_and_converged"])
        for row in summaries
    )
    if all_valid and all_feasible:
        verdict = "pilot_candidate"
    elif any("non-finite" in error or "solver" in error for error in errors) or any(
        row.get("finite_and_converged") is False for row in summaries
    ):
        verdict = "solver_or_curve_failure"
    else:
        verdict = "pilot_inconclusive"
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "verdict": verdict,
        "run_mode": args.run_mode,
        "source_run_id": args.source_run_id or None,
        "calculation_sha": args.calculation_sha,
        "postprocess_sha": args.postprocess_sha,
        "materialization_contract_version": MATERIALIZATION_CONTRACT_VERSION,
        "expected_target_count": len(expected),
        "materialized_target_count": len(summaries),
        "missing_target_ids": sorted(set(expected) - set(directories)),
        "unexpected_target_ids": unexpected,
        "errors": errors,
        "identity_fallbacks": identity_fallbacks,
        "identity_fallback_rule": "missing_summary_identity_uses_verified_provenance_and_manifest_v1",
        "solver_called": True,
        "reference_write": False,
        "oracle_labels_consumed": False,
        "target_list_sha256": sha256(args.target_list),
        "artifact_hashes": hashes,
    }
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    (args.output_dir / "materialization_diagnostics.json").write_text(
        json.dumps(
            {
                "materialization_contract_version": MATERIALIZATION_CONTRACT_VERSION,
                "identity_fallback_rule": "missing_summary_identity_uses_verified_provenance_and_manifest_v1",
                "target_identity_fallbacks": identity_fallbacks,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    (args.output_dir / "verdict.json").write_text(
        json.dumps({"verdict": verdict, "errors": errors, "stop_expansion": verdict != "pilot_candidate"}, indent=2) + "\n",
        encoding="utf-8",
    )
    write_csv(args.output_dir / "claim_ledger.csv", [
        {"claim_id": "target_matrix", "claim": "Exactly the 11 authorized Maxwell fixed-(xi,T) targets were materialized.", "status": "supported" if all_valid else "inconclusive", "evidence": "manifest.json; pilot_summary.csv", "boundary": "preflight pilot list only"},
        {"claim_id": "candidate_geometry", "claim": "Every target has a unique candidate and cross-level geometry convergence.", "status": "supported" if all_feasible else "inconclusive", "evidence": "pilot_summary.csv; curve_index.csv", "boundary": "diagnostic candidate, not production certificate"},
        {"claim_id": "identity_materialization", "claim": "Missing summary identities use only mutually verified provenance and manifest identities.", "status": "supported" if all_valid else "inconclusive", "evidence": "materialization_diagnostics.json; manifest.json", "boundary": "compatibility rule for v1 target artifacts"},
        {"claim_id": "reference_boundary", "claim": "The pilot does not write phase-reference and does not consume oracle labels.", "status": "supported", "evidence": "manifest.json; target provenance.json", "boundary": "phase-reference promotion remains separate"},
    ], ["claim_id", "claim", "status", "evidence", "boundary"])
    (args.output_dir / "README.md").write_text(
        f"# Issue #130 Maxwell CEP-local pilot\n\nverdict: `{verdict}`\n\n"
        "This artifact is diagnostic-only. It reruns the complete rho curve and the public strict candidate contract for the 11 authorized targets. It does not write phase-reference, read oracle labels, or alter C0/C1/C2 evidence.\n\n"
        f"Aggregate materialization contract: `{MATERIALIZATION_CONTRACT_VERSION}`. Missing identity fields in target summaries are accepted only when provenance and manifest agree with the requested immutable inputs; the fallback is recorded in manifest.json and materialization_diagnostics.json.\n",
        encoding="utf-8",
    )
    return 0 if verdict == "pilot_candidate" else 2


if __name__ == "__main__":
    raise SystemExit(main())
