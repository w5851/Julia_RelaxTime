#!/usr/bin/env python3
"""Audit p104-to-p128 convergence after the RS transport energy fix.

The input is a set of downloaded GitHub Actions artifacts plus a small JSON
manifest describing each run.  This script does not recompute physics and
does not promote artifacts.  It checks the convergence pair, verifies that
relaxation times remain unchanged relative to the superseded production case,
and records the expected transport-only semantic drift.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[3]
RESULT_ROOT = (
    PROJECT_ROOT
    / "data"
    / "outputs"
    / "results"
    / "relaxtime"
    / "transport"
    / "phase_guided"
)
MODE_A = "mode_a_fixed_muB_phase_scaled"
MODE_B = "mode_b_fixed_T_sparse_muB"

TAU_FIELDS = [
    "tau_u",
    "tau_d",
    "tau_s",
    "tau_ubar",
    "tau_dbar",
    "tau_sbar",
    "tauinv_u",
    "tauinv_d",
    "tauinv_s",
    "tauinv_ubar",
    "tauinv_dbar",
    "tauinv_sbar",
]
TRANSPORT_FIELDS = [
    "eta",
    "sigma",
    "zeta",
    "eta_over_s",
    "zeta_over_s",
    "sigma_over_T",
]
CONVERGENCE_FIELDS = TAU_FIELDS + TRANSPORT_FIELDS
DIAGNOSTIC_FIELDS = ["density", "rate", "contribution", "total", "tau_inv_species"]
NONNEGATIVE_SCAN_FIELDS = TAU_FIELDS + ["eta", "sigma", "zeta"]
NONNEGATIVE_DIAGNOSTIC_FIELDS = ["rate", "contribution", "total", "tau_inv_species"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download-root", required=True, type=Path)
    parser.add_argument("--run-manifest", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument(
        "--old-case",
        default="first_canonical_v1_p128_xi001_validated_anchored_prod_v1",
    )
    parser.add_argument("--convergence-rel-threshold", type=float, default=1e-2)
    parser.add_argument("--invariance-rel-threshold", type=float, default=1e-3)
    parser.add_argument("--abs-floor", type=float, default=1e-10)
    return parser.parse_args()


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    comments: list[str] = []
    data_lines: list[str] = []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for line in handle:
            if line.startswith("#"):
                comments.append(line.rstrip("\n"))
            elif line.strip():
                data_lines.append(line)
    if not data_lines:
        raise ValueError(f"CSV has no header: {path}")
    return comments, list(csv.DictReader(data_lines))


def write_csv(path: Path, header: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in header})


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_float(value: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError(f"non-finite numeric value: {value!r}")
    return parsed


def normalized(value: str | float) -> str:
    return f"{float(value):.10f}"


def point_key(mode: str, row: dict[str, str]) -> tuple[str, ...]:
    if mode == MODE_A:
        return (
            normalized(row["muB_MeV"]),
            normalized(row["alpha_T"]),
            normalized(row["xi"]),
        )
    if mode == MODE_B:
        return (
            normalized(row["T_MeV"]),
            normalized(row["muB_MeV"]),
            normalized(row["xi"]),
        )
    raise ValueError(f"unsupported mode: {mode}")


def rel_delta(new: float, old: float, abs_floor: float) -> float:
    return abs(new - old) / max(abs(new), abs(old), abs_floor)


def signed_rel_delta(new: float, old: float, abs_floor: float) -> float:
    return (new - old) / max(abs(old), abs_floor)


def one_file(root: Path, name: str) -> Path:
    matches = list(root.rglob(name))
    if len(matches) != 1:
        raise ValueError(f"expected one {name} under {root}, found {len(matches)}")
    return matches[0]


def load_run_artifacts(
    download_root: Path,
    manifest_rows: list[dict[str, Any]],
    source_commit: str,
) -> tuple[dict[tuple[str, str], list[dict[str, str]]], list[dict[str, Any]], list[str]]:
    scans: dict[tuple[str, str], list[dict[str, str]]] = {}
    audit_rows: list[dict[str, Any]] = []
    failures: list[str] = []

    for run in manifest_rows:
        run_id = str(run["run_id"])
        tier = str(run["tier"])
        mode = str(run["mode"])
        run_dir = download_root / run_id
        scan_path = one_file(run_dir, "phase_guided_transport_scan.csv")
        diag_path = one_file(run_dir, "channel_diagnostics.csv")
        failed_path = one_file(run_dir, "failed_points.csv")
        diagnostic_audit = one_file(run_dir, "DIAGNOSTIC_AUDIT.md")
        _, scan_rows = read_csv(scan_path)
        _, diag_rows = read_csv(diag_path)
        _, failed_rows = read_csv(failed_path)

        if str(run.get("head_sha", "")) != source_commit:
            failures.append(f"run {run_id}: head_sha mismatch")
        if str(run.get("conclusion", "")) != "success":
            failures.append(f"run {run_id}: conclusion is not success")
        if failed_rows:
            failures.append(f"run {run_id}: failed_points has {len(failed_rows)} rows")
        if not scan_rows:
            failures.append(f"run {run_id}: scan is empty")
            continue
        if {row.get("mode", "") for row in scan_rows} != {mode}:
            failures.append(f"run {run_id}: scan mode does not match manifest")
        if any(row.get("converged", "").lower() != "true" for row in scan_rows):
            failures.append(f"run {run_id}: unconverged scan row")

        nonfinite_scan = 0
        negative_scan = 0
        for row in scan_rows:
            for field in CONVERGENCE_FIELDS:
                try:
                    value = parse_float(row[field])
                except (KeyError, ValueError):
                    nonfinite_scan += 1
                    continue
                if field in NONNEGATIVE_SCAN_FIELDS and value < 0.0:
                    negative_scan += 1

        nonfinite_diag = 0
        negative_diag = 0
        for row in diag_rows:
            for field in DIAGNOSTIC_FIELDS:
                try:
                    value = parse_float(row[field])
                except (KeyError, ValueError):
                    nonfinite_diag += 1
                    continue
                if field in NONNEGATIVE_DIAGNOSTIC_FIELDS and value < 0.0:
                    negative_diag += 1

        if nonfinite_scan or negative_scan or nonfinite_diag or negative_diag:
            failures.append(
                f"run {run_id}: validity scan(nonfinite={nonfinite_scan}, negative={negative_scan}) "
                f"diagnostics(nonfinite={nonfinite_diag}, negative={negative_diag})"
            )

        key = (tier, mode)
        if key in scans:
            failures.append(f"duplicate manifest tier/mode pair: {tier}/{mode}")
        scans[key] = scan_rows
        audit_rows.append(
            {
                "run_id": run_id,
                "tier": tier,
                "mode": mode,
                "label": run.get("label", ""),
                "url": run.get("url", ""),
                "head_sha": run.get("head_sha", ""),
                "artifact_name": run.get("artifact_name", ""),
                "scan_rows": len(scan_rows),
                "diagnostic_rows": len(diag_rows),
                "failed_rows": len(failed_rows),
                "scan_nonfinite": nonfinite_scan,
                "scan_negative": negative_scan,
                "diagnostic_nonfinite": nonfinite_diag,
                "diagnostic_negative": negative_diag,
                "scan_sha256": sha256_file(scan_path),
                "diagnostics_sha256": sha256_file(diag_path),
                "failed_points_sha256": sha256_file(failed_path),
                "diagnostic_audit_sha256": sha256_file(diagnostic_audit),
            }
        )
    return scans, audit_rows, failures


def compare_convergence(
    scans: dict[tuple[str, str], list[dict[str, str]]],
    abs_floor: float,
) -> tuple[list[dict[str, Any]], list[str]]:
    rows: list[dict[str, Any]] = []
    failures: list[str] = []
    for mode in (MODE_A, MODE_B):
        low = {point_key(mode, row): row for row in scans.get(("p104", mode), [])}
        high = {point_key(mode, row): row for row in scans.get(("p128", mode), [])}
        if set(low) != set(high):
            failures.append(
                f"{mode}: p104/p128 point keys differ; "
                f"missing_low={len(set(high) - set(low))}, missing_high={len(set(low) - set(high))}"
            )
        for key in sorted(set(low) & set(high)):
            for field in CONVERGENCE_FIELDS:
                old = parse_float(low[key][field])
                new = parse_float(high[key][field])
                rows.append(
                    {
                        "mode": mode,
                        "key": "|".join(key),
                        "xi": high[key]["xi"],
                        "field": field,
                        "p104": old,
                        "p128": new,
                        "abs_delta": abs(new - old),
                        "rel_delta": rel_delta(new, old, abs_floor),
                        "signed_rel_delta": signed_rel_delta(new, old, abs_floor),
                    }
                )
    return rows, failures


def compare_old_new(
    scans: dict[tuple[str, str], list[dict[str, str]]],
    old_case: str,
    abs_floor: float,
) -> tuple[list[dict[str, Any]], list[str]]:
    rows: list[dict[str, Any]] = []
    failures: list[str] = []
    for mode in (MODE_A, MODE_B):
        old_path = RESULT_ROOT / mode / old_case / "phase_guided_transport_scan.csv"
        _, old_rows = read_csv(old_path)
        old_by_key = {point_key(mode, row): row for row in old_rows}
        for new_row in scans.get(("p128", mode), []):
            key = point_key(mode, new_row)
            old_row = old_by_key.get(key)
            if old_row is None:
                failures.append(f"{mode}: old case is missing {'|'.join(key)}")
                continue
            for field in CONVERGENCE_FIELDS:
                old = parse_float(old_row[field])
                new = parse_float(new_row[field])
                if field in TAU_FIELDS:
                    role = "tau_invariance"
                elif abs(float(new_row["xi"])) <= 1e-12:
                    role = "xi0_transport_invariance"
                else:
                    role = "transport_semantic_drift"
                rows.append(
                    {
                        "mode": mode,
                        "key": "|".join(key),
                        "xi": new_row["xi"],
                        "field": field,
                        "role": role,
                        "old": old,
                        "new": new,
                        "abs_delta": abs(new - old),
                        "rel_delta": rel_delta(new, old, abs_floor),
                        "signed_rel_delta": signed_rel_delta(new, old, abs_floor),
                    }
                )
    return rows, failures


def field_summary(rows: list[dict[str, Any]], field_key: str = "field") -> dict[str, Any]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[str(row[field_key])].append(row)
    summary: dict[str, Any] = {}
    for field, values in grouped.items():
        worst = max(values, key=lambda item: float(item["rel_delta"]))
        signed = [float(item.get("signed_rel_delta", 0.0)) for item in values]
        summary[field] = {
            "count": len(values),
            "max_rel": float(worst["rel_delta"]),
            "worst": worst,
            "min_signed_rel": min(signed),
            "max_signed_rel": max(signed),
        }
    return summary


def main() -> int:
    args = parse_args()
    if args.out_dir.exists():
        raise FileExistsError(f"output directory already exists: {args.out_dir}")
    manifest_rows = json.loads(args.run_manifest.read_text(encoding="utf-8-sig"))
    expected_pairs = {
        ("p104", MODE_A),
        ("p104", MODE_B),
        ("p128", MODE_A),
        ("p128", MODE_B),
    }
    manifest_pairs = {(str(row["tier"]), str(row["mode"])) for row in manifest_rows}
    failures: list[str] = []
    if manifest_pairs != expected_pairs:
        failures.append(f"run manifest pairs differ from expected: {sorted(manifest_pairs)}")

    scans, run_audit_rows, load_failures = load_run_artifacts(
        args.download_root, manifest_rows, args.source_commit
    )
    failures.extend(load_failures)
    convergence_rows, convergence_failures = compare_convergence(scans, args.abs_floor)
    failures.extend(convergence_failures)
    old_new_rows, old_new_failures = compare_old_new(
        scans, args.old_case, args.abs_floor
    )
    failures.extend(old_new_failures)

    convergence_worst = max(
        convergence_rows, key=lambda row: float(row["rel_delta"]), default=None
    )
    convergence_max = float(convergence_worst["rel_delta"]) if convergence_worst else math.inf
    if convergence_max > args.convergence_rel_threshold:
        failures.append(
            f"p104/p128 max relative drift {convergence_max} exceeds "
            f"{args.convergence_rel_threshold}"
        )

    tau_rows = [row for row in old_new_rows if row["field"] in TAU_FIELDS]
    tau_worst = max(tau_rows, key=lambda row: float(row["rel_delta"]), default=None)
    tau_max = float(tau_worst["rel_delta"]) if tau_worst else math.inf
    if tau_max > args.invariance_rel_threshold:
        failures.append(
            f"tau invariance max relative drift {tau_max} exceeds "
            f"{args.invariance_rel_threshold}"
        )

    xi0_transport_rows = [
        row
        for row in old_new_rows
        if row["field"] in TRANSPORT_FIELDS and abs(float(row["xi"])) <= 1e-12
    ]
    xi0_worst = max(
        xi0_transport_rows, key=lambda row: float(row["rel_delta"]), default=None
    )
    xi0_max = float(xi0_worst["rel_delta"]) if xi0_worst else math.inf
    if xi0_max > args.invariance_rel_threshold:
        failures.append(
            f"xi=0 transport invariance max relative drift {xi0_max} exceeds "
            f"{args.invariance_rel_threshold}"
        )

    args.out_dir.mkdir(parents=True, exist_ok=False)
    write_csv(
        args.out_dir / "download_manifest.csv",
        list(run_audit_rows[0].keys()) if run_audit_rows else [],
        run_audit_rows,
    )
    write_csv(
        args.out_dir / "p104_vs_p128_convergence.csv",
        [
            "mode",
            "key",
            "xi",
            "field",
            "p104",
            "p128",
            "abs_delta",
            "rel_delta",
            "signed_rel_delta",
        ],
        convergence_rows,
    )
    write_csv(
        args.out_dir / "old_vs_new_semantics_comparison.csv",
        [
            "mode",
            "key",
            "xi",
            "field",
            "role",
            "old",
            "new",
            "abs_delta",
            "rel_delta",
            "signed_rel_delta",
        ],
        old_new_rows,
    )

    status_counts = Counter(
        row.get("converged", "missing")
        for values in scans.values()
        for row in values
    )
    summary = {
        "schema": "phase_guided_energy_semantics_convergence_gate_v1",
        "verdict": "production-grade" if not failures else "blocked",
        "source_commit": args.source_commit,
        "old_case": args.old_case,
        "convergence_threshold": args.convergence_rel_threshold,
        "invariance_threshold": args.invariance_rel_threshold,
        "abs_floor": args.abs_floor,
        "run_count": len(run_audit_rows),
        "status_counts": dict(status_counts),
        "validity_counts": {
            "failed_rows": sum(int(row["failed_rows"]) for row in run_audit_rows),
            "scan_nonfinite": sum(int(row["scan_nonfinite"]) for row in run_audit_rows),
            "scan_negative": sum(int(row["scan_negative"]) for row in run_audit_rows),
            "diagnostic_nonfinite": sum(
                int(row["diagnostic_nonfinite"]) for row in run_audit_rows
            ),
            "diagnostic_negative": sum(
                int(row["diagnostic_negative"]) for row in run_audit_rows
            ),
        },
        "p104_vs_p128": {
            "comparison_cells": len(convergence_rows),
            "max_rel": convergence_max,
            "worst": convergence_worst,
            "field_summary": field_summary(convergence_rows),
        },
        "old_vs_new": {
            "comparison_cells": len(old_new_rows),
            "tau_invariance_max_rel": tau_max,
            "tau_invariance_worst": tau_worst,
            "xi0_transport_invariance_max_rel": xi0_max,
            "xi0_transport_invariance_worst": xi0_worst,
            "field_summary": field_summary(old_new_rows),
        },
        "failures": failures,
    }
    write_json(args.out_dir / "convergence_summary.json", summary)
    with (args.out_dir / "README.md").open(
        "w", encoding="utf-8", newline="\n"
    ) as handle:
        handle.write(
            "\n".join(
            [
                "# Issue #130 RS transport energy convergence gate",
                "",
                f"Verdict: `{summary['verdict']}`",
                "",
                f"- source commit: `{args.source_commit}`",
                f"- GitHub Actions runs: `{len(run_audit_rows)}`",
                f"- p104 -> p128 worst relative drift: `{convergence_max:.8g}`",
                f"- tau old/new invariance worst relative drift: `{tau_max:.8g}`",
                f"- xi=0 transport old/new invariance worst relative drift: `{xi0_max:.8g}`",
                f"- failed checks: `{len(failures)}`",
                "",
                "This gate covers the old convergence audit's most sensitive mode-A and mode-B panels",
                "at xi = -0.5, -0.1, 0.0, 0.35, 0.5.  It is a prerequisite for the full",
                "high-resolution production rerun; the Action artifacts themselves remain diagnostic-only.",
                "",
            ]
            )
        )
    write_json(
        args.out_dir / "manifest.json",
        {
            "format": "phase_guided_energy_semantics_convergence_manifest_v1",
            "verdict": summary["verdict"],
            "source_commit": args.source_commit,
            "files": [
                "README.md",
                "download_manifest.csv",
                "p104_vs_p128_convergence.csv",
                "old_vs_new_semantics_comparison.csv",
                "convergence_summary.json",
            ],
            "hashes": {
                name: sha256_file(args.out_dir / name)
                for name in (
                    "README.md",
                    "download_manifest.csv",
                    "p104_vs_p128_convergence.csv",
                    "old_vs_new_semantics_comparison.csv",
                    "convergence_summary.json",
                )
            },
        },
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
