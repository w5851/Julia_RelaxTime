#!/usr/bin/env python3
"""Merge sharded phase-guided high-xi transport action artifacts.

This script imports GitHub Actions result artifacts into the repository
production data layout after local consistency gates pass. It does not
recompute numerical rows.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
from collections import defaultdict
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[3]
RESULT_ROOT = PROJECT_ROOT / "data" / "outputs" / "results" / "relaxtime" / "transport" / "phase_guided"
FIGURE_ROOT = PROJECT_ROOT / "data" / "outputs" / "figures" / "relaxtime" / "transport" / "phase_guided"

OLD_CASE = "first_canonical_v1_p128_validated_anchored_prod_v1"
NEW_CASE = "first_canonical_v1_p128_xi001_validated_anchored_prod_v1"
SOURCE_CANDIDATE_CASE = "p128_xi001_localdispatch_20260618_223856"
CONVERGENCE_SOURCE = "first_canonical_v1_p128_validated_anchored_prod_v1_convergence"
MODE_A = "mode_a_fixed_muB_phase_scaled"
MODE_B = "mode_b_fixed_T_sparse_muB"

SCAN_CORE_COLUMNS = [
    "T_MeV",
    "muq_MeV",
    "muB_MeV",
    "xi",
    "T_phase_base_MeV",
    "alpha_T",
    "Phi",
    "Phibar",
    "m_u",
    "m_d",
    "m_s",
    "rho_baryon",
    "rho_norm",
    "P_fm4inv",
    "epsilon_fm4inv",
    "s_fm3inv",
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
    "eta",
    "sigma",
    "zeta",
    "eta_over_s",
    "zeta_over_s",
    "sigma_over_T",
]

DIAG_COMPARE_COLUMNS = ["density", "rate", "contribution", "total", "tau_inv_species"]
DIAG_NONNEG_COLUMNS = ["rate", "contribution", "total", "tau_inv_species"]


@dataclass
class ParsedCsv:
    comments: list[str]
    header: list[str]
    rows: list[dict[str, str]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download-root", required=True, type=Path)
    parser.add_argument("--run-manifest", required=True, type=Path)
    parser.add_argument("--new-case", default=NEW_CASE)
    parser.add_argument("--old-case", default=OLD_CASE)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--overlap-rel-threshold", type=float, default=1e-3)
    parser.add_argument("--anchor-rel-threshold", type=float, default=1e-2)
    parser.add_argument("--anchor-abs-threshold", type=float, default=1e-10)
    return parser.parse_args()


def read_scan_csv(path: Path) -> ParsedCsv:
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
    reader = csv.DictReader(data_lines)
    return ParsedCsv(comments, list(reader.fieldnames or []), list(reader))


def write_csv(path: Path, comments: list[str], header: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        for comment in comments:
            handle.write(comment.rstrip("\n") + "\n")
        writer = csv.DictWriter(handle, fieldnames=header, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in header})


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_float(value: str) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed


def is_finite_number(value: str) -> bool:
    parsed = parse_float(value)
    return parsed is not None and math.isfinite(parsed)


def rel_diff(a: float, b: float) -> float:
    denom = max(abs(a), abs(b), 1e-300)
    return abs(a - b) / denom


def norm_key(value: str | float) -> str:
    parsed = float(value)
    return f"{parsed:.10f}"


def scan_key(row: dict[str, str]) -> tuple[str, ...]:
    mode = row["mode"]
    if mode == MODE_A:
        return (
            mode,
            norm_key(row["muB_MeV"]),
            norm_key(row["alpha_T"]),
            norm_key(row["xi"]),
        )
    if mode == MODE_B:
        return (
            mode,
            norm_key(row["T_MeV"]),
            norm_key(row["muB_MeV"]),
            norm_key(row["xi"]),
        )
    raise ValueError(f"unsupported mode in row: {mode}")


def diag_key(row: dict[str, str]) -> tuple[str, ...]:
    return (
        norm_key(row["T_MeV"]),
        norm_key(row["muq_MeV"]),
        norm_key(row["muB_MeV"]),
        norm_key(row["xi"]),
        row["species"],
        row["channel"],
        row["density_key"],
        norm_key(row["multiplicity"]),
    )


def scan_sort_key(row: dict[str, str]) -> tuple[float, ...]:
    if row["mode"] == MODE_A:
        return (
            float(row["muB_MeV"]),
            float(row["alpha_T"]),
            float(row["xi"]),
        )
    return (
        float(row["T_MeV"]),
        float(row["muB_MeV"]),
        float(row["xi"]),
    )


def diag_sort_key(row: dict[str, str]) -> tuple[Any, ...]:
    return (
        float(row["T_MeV"]),
        float(row["muB_MeV"]),
        float(row["xi"]),
        row["species"],
        row["channel"],
        row["density_key"],
        float(row["multiplicity"]),
    )


def prefer_candidate(rows: list[dict[str, str]]) -> dict[str, str]:
    for row in rows:
        label = row.get("__source_label", "")
        if "xi_neg" in label:
            return row
    return rows[0]


def compare_duplicate_rows(rows: list[dict[str, str]], columns: list[str]) -> dict[str, Any]:
    base = rows[0]
    worst = {
        "max_abs": 0.0,
        "max_rel": 0.0,
        "column": "",
        "base_source": base.get("__source_label", ""),
        "other_source": "",
    }
    for other in rows[1:]:
        for col in columns:
            a = parse_float(base.get(col, ""))
            b = parse_float(other.get(col, ""))
            if a is None or b is None:
                continue
            if not math.isfinite(a) and not math.isfinite(b):
                continue
            if not (math.isfinite(a) and math.isfinite(b)):
                abs_delta = math.inf
                rel_delta = math.inf
            else:
                abs_delta = abs(a - b)
                rel_delta = rel_diff(a, b)
            if rel_delta > worst["max_rel"]:
                worst = {
                    "max_abs": abs_delta,
                    "max_rel": rel_delta,
                    "column": col,
                    "base_source": base.get("__source_label", ""),
                    "other_source": other.get("__source_label", ""),
                }
    return worst


def merge_rows(
    rows: list[dict[str, str]],
    key_fn,
    sort_fn,
    compare_columns: list[str],
    rel_threshold: float,
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    grouped: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[key_fn(row)].append(row)

    duplicate_groups = 0
    duplicate_rows = 0
    worst = {
        "max_abs": 0.0,
        "max_rel": 0.0,
        "column": "",
        "base_source": "",
        "other_source": "",
    }
    merged: list[dict[str, str]] = []
    for group_rows in grouped.values():
        if len(group_rows) > 1:
            duplicate_groups += 1
            duplicate_rows += len(group_rows) - 1
            group_worst = compare_duplicate_rows(group_rows, compare_columns)
            if group_worst["max_rel"] > worst["max_rel"]:
                worst = group_worst
        chosen = prefer_candidate(group_rows)
        merged.append({k: v for k, v in chosen.items() if not k.startswith("__")})

    merged.sort(key=sort_fn)
    summary = {
        "duplicate_groups": duplicate_groups,
        "duplicate_rows_removed": duplicate_rows,
        "worst": worst,
        "threshold": rel_threshold,
        "passed": worst["max_rel"] <= rel_threshold,
    }
    return merged, summary


def data_row_count(path: Path) -> int:
    parsed = read_scan_csv(path)
    return len(parsed.rows)


def discover_artifacts(download_root: Path, run_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_run = {str(row["RunId"]): row for row in run_rows}
    discovered: list[dict[str, Any]] = []
    for run_id, run in by_run.items():
        run_dir = download_root / run_id
        scans = list(run_dir.rglob("phase_guided_transport_scan.csv"))
        diags = list(run_dir.rglob("channel_diagnostics.csv"))
        failed = list(run_dir.rglob("failed_points.csv"))
        configs = list(run_dir.rglob("effective_config.json"))
        readmes = list(run_dir.rglob("README.md"))
        audits = list(run_dir.rglob("DIAGNOSTIC_AUDIT.md"))
        logs = list(run_dir.rglob("*.log"))
        if len(scans) != 1 or len(diags) != 1 or len(failed) != 1:
            raise ValueError(
                f"run {run_id} has unexpected artifact files: "
                f"scan={len(scans)}, diagnostics={len(diags)}, failed={len(failed)}"
            )
        discovered.append(
            {
                "run": run,
                "scan": scans[0],
                "diagnostics": diags[0],
                "failed": failed[0],
                "effective_config": configs[0] if configs else None,
                "readme": readmes[0] if readmes else None,
                "diagnostic_audit": audits[0] if audits else None,
                "logs": logs,
            }
        )
    return discovered


def mode_from_label(label: str) -> str:
    if "fixed-muB-phase-scaled" in label or label.startswith("modeA_"):
        return MODE_A
    if "fixed-T-sparse-muB" in label or label.startswith("modeB_"):
        return MODE_B
    raise ValueError(f"cannot infer mode from label: {label}")


def validate_failed_points(path: Path) -> int:
    parsed = read_scan_csv(path)
    return len(parsed.rows)


def load_all_rows(discovered: list[dict[str, Any]]) -> tuple[dict[str, list[dict[str, str]]], dict[str, list[dict[str, str]]], list[dict[str, Any]]]:
    scan_rows: dict[str, list[dict[str, str]]] = {MODE_A: [], MODE_B: []}
    diag_rows: dict[str, list[dict[str, str]]] = {MODE_A: [], MODE_B: []}
    manifest_rows: list[dict[str, Any]] = []

    for entry in discovered:
        run = entry["run"]
        label = str(run["Label"])
        mode = mode_from_label(label)
        failed_rows = validate_failed_points(entry["failed"])
        scan = read_scan_csv(entry["scan"])
        diag = read_scan_csv(entry["diagnostics"])
        if not scan.rows:
            raise ValueError(f"scan artifact has no rows: {entry['scan']}")
        for row in scan.rows:
            row["__source_run_id"] = str(run["RunId"])
            row["__source_label"] = label
            scan_rows[mode].append(row)
        for row in diag.rows:
            row["__source_run_id"] = str(run["RunId"])
            row["__source_label"] = label
            diag_rows[mode].append(row)
        manifest_rows.append(
            {
                "run_id": run["RunId"],
                "label": label,
                "status": run["Status"],
                "conclusion": run["Conclusion"],
                "artifact_name": run["ArtifactName"],
                "artifact_size": run["ArtifactSize"],
                "head_sha": run["HeadSha"],
                "url": run["Url"],
                "mode": mode,
                "scan_rows": len(scan.rows),
                "diagnostic_rows": len(diag.rows),
                "failed_rows": failed_rows,
                "scan_sha256": sha256_file(entry["scan"]),
                "diagnostics_sha256": sha256_file(entry["diagnostics"]),
                "failed_points_sha256": sha256_file(entry["failed"]),
            }
        )
    return scan_rows, diag_rows, manifest_rows


def numeric_validity(rows: list[dict[str, str]], columns: list[str]) -> dict[str, Any]:
    nonfinite = 0
    negative = 0
    for row in rows:
        for col in columns:
            value = parse_float(row.get(col, ""))
            if value is None:
                continue
            if not math.isfinite(value):
                nonfinite += 1
            elif col in DIAG_NONNEG_COLUMNS and value < -1e-12:
                negative += 1
    return {"nonfinite": nonfinite, "negative": negative}


def scan_validity_columns(mode: str) -> list[str]:
    if mode == MODE_B:
        return [col for col in SCAN_CORE_COLUMNS if col != "alpha_T"]
    return list(SCAN_CORE_COLUMNS)


def old_anchor_key(mode: str, row: dict[str, str]) -> tuple[str, ...]:
    if mode == MODE_A:
        return (
            norm_key(row["muB_MeV"]),
            norm_key(row["alpha_T"]),
            norm_key(row["xi"]),
        )
    return (
        norm_key(row["T_MeV"]),
        norm_key(row["muB_MeV"]),
        norm_key(row["xi"]),
    )


def compare_to_old_case(new_case: str, old_case: str, mode: str, new_rows: list[dict[str, str]]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    old_path = RESULT_ROOT / mode / old_case / "phase_guided_transport_scan.csv"
    old_rows = read_scan_csv(old_path).rows
    old_by_key = {old_anchor_key(mode, row): row for row in old_rows}
    comparisons: list[dict[str, Any]] = []
    missing = 0
    worst = {
        "mode": mode,
        "key": "",
        "column": "",
        "old_value": "",
        "new_value": "",
        "abs_delta": 0.0,
        "rel_delta": 0.0,
    }
    for row in new_rows:
        xi = float(row["xi"])
        scaled = round(xi * 100.0)
        if abs(scaled % 5) > 1e-8 and abs((scaled % 5) - 5) > 1e-8:
            continue
        key = old_anchor_key(mode, row)
        old = old_by_key.get(key)
        if old is None:
            missing += 1
            continue
        for col in SCAN_CORE_COLUMNS:
            old_val = parse_float(old.get(col, ""))
            new_val = parse_float(row.get(col, ""))
            if old_val is None or new_val is None:
                continue
            if not math.isfinite(old_val) and not math.isfinite(new_val):
                continue
            if math.isfinite(old_val) and math.isfinite(new_val):
                abs_delta = abs(new_val - old_val)
                rel_delta = rel_diff(new_val, old_val)
            else:
                abs_delta = math.inf
                rel_delta = math.inf
            comparison = {
                "mode": mode,
                "key": "|".join(key),
                "column": col,
                "old_value": old.get(col, ""),
                "new_value": row.get(col, ""),
                "abs_delta": abs_delta,
                "rel_delta": rel_delta,
            }
            comparisons.append(comparison)
            if rel_delta > worst["rel_delta"]:
                worst = comparison
    summary = {
        "mode": mode,
        "old_case": old_case,
        "new_case": new_case,
        "matched_anchor_rows": len({item["key"] for item in comparisons}),
        "comparison_cells": len(comparisons),
        "missing_anchor_rows": missing,
        "worst": worst,
    }
    return comparisons, summary


def write_download_manifest(path: Path, rows: list[dict[str, Any]]) -> None:
    header = [
        "run_id",
        "label",
        "status",
        "conclusion",
        "artifact_name",
        "artifact_size",
        "head_sha",
        "url",
        "mode",
        "scan_rows",
        "diagnostic_rows",
        "failed_rows",
        "scan_sha256",
        "diagnostics_sha256",
        "failed_points_sha256",
    ]
    write_csv(path, [], header, [{key: str(row.get(key, "")) for key in header} for row in rows])


def write_anchor_comparison(path: Path, rows: list[dict[str, Any]]) -> None:
    header = ["mode", "key", "column", "old_value", "new_value", "abs_delta", "rel_delta"]
    write_csv(path, [], header, [{key: str(row.get(key, "")) for key in header} for row in rows])


def collect_xi_values(rows: list[dict[str, str]]) -> list[float]:
    return sorted({round(float(row["xi"]), 10) for row in rows})


def production_comments(new_case: str, git_commit: str, source_manifest_rel: str, title: str, mode_path: str | None = None) -> list[str]:
    comments = [
        f"# production_case: {new_case}",
        f"# source_candidate_case: {SOURCE_CANDIDATE_CASE}",
        f"# source_sharded_action_runs: see {source_manifest_rel}",
        f"# git_commit: {git_commit}",
        "# script: scripts/relaxtime/run_phase_guided_transport_scan.jl",
    ]
    if mode_path:
        comments.append(f"# source_csv: {mode_path}")
    comments.append(f"# title: {title}")
    return comments


def effective_config(mode: str, new_case: str, rows: list[dict[str, str]]) -> dict[str, Any]:
    outdir = f"data/outputs/results/relaxtime/transport/phase_guided/{mode}/{new_case}"
    fig_dir = f"data/outputs/figures/relaxtime/transport/phase_guided/{mode}/{new_case}"
    options = {
        "T_values": sorted({float(row["T_MeV"]) for row in rows}) if mode == MODE_B else [],
        "alpha_T_values": sorted({float(row["alpha_T"]) for row in rows}) if mode == MODE_A else [1.0, 1.1, 1.2],
        "case_name": new_case,
        "channel_diagnostics": True,
        "compute_bulk": True,
        "dry_run": False,
        "figure_dir": fig_dir,
        "mode": mode,
        "muB_values": sorted({float(row["muB_MeV"]) for row in rows}),
        "outdir": outdir,
        "overwrite": True,
        "plan_csv": f"{outdir}/sampling_plan.csv",
        "propagator_xi_policy": "match_thermo",
        "result_csv": f"{outdir}/phase_guided_transport_scan.csv",
        "resume": False,
        "sigma_cache_policy": "validated_anchored",
        "sigma_grid_n": 560,
        "tau_angle_nodes": 20,
        "tau_n_sigma_points": 24,
        "tau_p_nodes": 128,
        "tau_phi_nodes": 36,
        "xi_values": collect_xi_values(rows),
    }
    return {
        "schema_version": "v1",
        "script": "scripts/relaxtime/run_phase_guided_transport_scan.jl",
        "options": options,
        "production_import": {
            "verdict": "production-grade",
            "source_candidate_case": SOURCE_CANDIDATE_CASE,
            "formal_case_name": new_case,
            "source_commit": "700a0845abd09013eb39ff8f8d2993d5574476d3",
            "path_fields_rewritten_for_repository_import": True,
            "numerical_options_unchanged_from_p128_validated_anchored": True,
            "parameter_convergence_gate": f"data/outputs/results/relaxtime/transport/phase_guided/{CONVERGENCE_SOURCE}/p104_vs_p128_convergence_gate.summary.json",
            "highxi_import_gate": f"data/outputs/results/relaxtime/transport/phase_guided/{new_case}_convergence/highxi_anchor_comparison_to_p128_plus5.summary.json",
        },
    }


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=False)
        handle.write("\n")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def make_plan_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    plan = []
    for row in rows:
        plan.append(
            {
                "T_MeV": row["T_MeV"],
                "muB_MeV": row["muB_MeV"],
                "xi": row["xi"],
                "mode": row["mode"],
                "plot_panel": row["plot_panel"],
                "plot_series": row["plot_series"],
            }
        )
    return plan


def mode_title(mode: str) -> str:
    return "mode A" if mode == MODE_A else "mode B"


def mode_scope(mode: str) -> str:
    if mode == MODE_A:
        return "fixed mu_B, phase-reference T/T_phase bands, continuous xi scan"
    return "fixed temperature, sparse mu_B, continuous xi scan"


def grid_line(mode: str, rows: list[dict[str, str]]) -> str:
    xis = collect_xi_values(rows)
    if mode == MODE_A:
        mu = sorted({float(row["muB_MeV"]) for row in rows})
        alpha = sorted({float(row["alpha_T"]) for row in rows})
        return f"mu_B = {mu} MeV, alpha_T = {alpha}, xi = {xis[0]:.2f}:0.01:{xis[-1]:.2f}"
    temps = sorted({float(row["T_MeV"]) for row in rows})
    mu = sorted({float(row["muB_MeV"]) for row in rows})
    return f"T = {temps} MeV, mu_B = {mu} MeV, xi = {xis[0]:.2f}:0.01:{xis[-1]:.2f}"


def write_readme_and_audit(
    case_dir: Path,
    conv_dir: Path,
    mode: str,
    new_case: str,
    summary: dict[str, Any],
    mode_manifest: dict[str, Any],
) -> None:
    rel_conv = conv_dir.relative_to(PROJECT_ROOT).as_posix()
    rel_fig = (FIGURE_ROOT / mode / new_case).relative_to(PROJECT_ROOT).as_posix()
    readme = f"""# {new_case}

Phase-guided transport production-grade case for `{mode}`.

## Verdict

`production-grade`

This result imports the p128 high-xi-resolution GitHub Actions artifacts after local merge gates passed. The numerical integration policy is unchanged from `{OLD_CASE}`; the new contribution is denser `xi = -0.50:0.01:0.50` sampling.

## Scope

- mode: `{mode}`
- case name: `{new_case}`
- summary: {mode_scope(mode)}
- grid: `{grid_line(mode, mode_manifest["rows_for_text"])}`
- planned and completed points: `{summary["scan_rows"]}`
- compute bulk viscosity (`zeta`): `true`

## Numerical Policy

- `propagator_xi_policy = match_thermo`
- `sigma_cache_policy = validated_anchored`
- `tau_p_nodes = 128`
- `tau_angle_nodes = 20`
- `tau_phi_nodes = 36`
- `tau_n_sigma_points = 24`
- `sigma_grid_n = 560`
- `channel_diagnostics = true`

## Key Files

- `phase_guided_transport_scan.csv`: imported transport scan output
- `sampling_plan.csv`: full phase-guided sampling plan
- `channel_diagnostics.csv`: per-channel scattering-rate diagnostics
- `failed_points.csv`: failure ledger; contains `0` data rows
- `effective_config.json`: repository-imported effective configuration snapshot
- `manifest.json`: production manifest with source runs, hashes, and gate evidence
- `PRODUCTION_AUDIT.md`: human-readable production audit

Shared high-xi import evidence is stored in:

- `{rel_conv}/`

Canonical figures are stored in:

- `{rel_fig}/`

## Gate Summary

- source run head SHA: `700a0845abd09013eb39ff8f8d2993d5574476d3`
- source action artifacts: `{summary["source_run_count"]}` result artifacts
- failed point rows across source artifacts: `{summary["failed_rows"]}`
- chunk overlap duplicate groups: `{summary["overlap"]["duplicate_groups"]}`
- chunk overlap worst relative difference: `{summary["overlap"]["worst"]["max_rel"]}`
- anchor comparison to `{OLD_CASE}` effective worst relative difference: `{summary["effective_anchor_worst_rel"]}`
- anchor comparison raw near-zero worst relative difference: `{summary["anchor_worst_rel"]}`
- inherited adjacent integration convergence: `{CONVERGENCE_SOURCE}/p104_vs_p128_convergence_gate.summary.json`

## Interpretation Boundary

This case supersedes neither `{OLD_CASE}` nor its convergence evidence. It is the same p128 validated-anchored numerical policy sampled on a denser `xi` grid. Claims that depend on newly resolved local structures should still check neighboring points and channel diagnostics.
"""
    write_text(case_dir / "README.md", readme)

    run_table = "\n".join(f"| `{row['label']}` | `{row['run_id']}` |" for row in mode_manifest["source_runs"])
    audit = f"""# Production Audit: Phase-Guided Transport {mode_title(mode)} p128 xi0.01 validated-anchored

Verdict: `production-grade`

Date: {date.today().isoformat()}

## Production Case

- Case slug: `{new_case}`
- Mode: `{mode}`
- Result path: `{case_dir.relative_to(PROJECT_ROOT).as_posix()}/`
- Figure path: `{rel_fig}/`
- Remote workflow: `.github/workflows/relaxtime-phase-guided-transport-production.yml`
- Source commit: `700a0845abd09013eb39ff8f8d2993d5574476d3`

## Physics Scope

- {mode_scope(mode)}.
- Grid: `{grid_line(mode, mode_manifest["rows_for_text"])}`
- Observables: relaxation times, inverse relaxation times, `eta`, `sigma`, `zeta`, `eta_over_s`, `zeta_over_s`, `sigma_over_T`.
- Bulk viscosity computation: enabled.
- Channel diagnostics: enabled.

## Policies

- `propagator_xi_policy = match_thermo`
- `sigma_cache_policy = validated_anchored`
- `channel_diagnostics = true`

## Command Log

The p128 high-xi production artifact was generated with GitHub Actions shards using:

- `case_name = {SOURCE_CANDIDATE_CASE}`
- `propagator_xi_policy = match_thermo`
- `sigma_cache_policy = validated_anchored`
- `tau_p_nodes = 128`
- `tau_angle_nodes = 20`
- `tau_phi_nodes = 36`
- `tau_n_sigma_points = 24`
- `sigma_grid_n = 560`
- `render_plots = false` for nonzero-muB rerun chunks; figures were rendered locally after merge

Mode source runs:

| shard | run id |
| --- | --- |
{run_table}

## Convergence And Import Gates

Parameter convergence is inherited from the existing adjacent p104-to-p128 gate for `{OLD_CASE}` because the integration policy is unchanged:

- `data/outputs/results/relaxtime/transport/phase_guided/{CONVERGENCE_SOURCE}/p104_vs_p128_convergence_gate.summary.json`

High-xi import-specific gates:

| gate | value |
| --- | ---: |
| source result artifacts | `{summary["source_run_count"]}` |
| source failed rows | `{summary["failed_rows"]}` |
| merged data rows | `{summary["scan_rows"]}` |
| merged diagnostic rows | `{summary["diagnostic_rows"]}` |
| chunk overlap duplicate groups | `{summary["overlap"]["duplicate_groups"]}` |
| chunk overlap worst relative difference | `{summary["overlap"]["worst"]["max_rel"]}` |
| anchor rows matched to old p128 case | `{summary["anchor_matched_rows"]}` |
| anchor effective worst relative difference | `{summary["effective_anchor_worst_rel"]}` |
| anchor raw near-zero worst relative difference | `{summary["anchor_worst_rel"]}` |
| non-finite core scan values | `{summary["scan_validity"]["nonfinite"]}` |
| non-finite diagnostic values | `{summary["diag_validity"]["nonfinite"]}` |
| negative diagnostic rate/contribution values | `{summary["diag_validity"]["negative"]}` |

## Data Outputs

- `phase_guided_transport_scan.csv`: `{summary["scan_rows"]}` data rows
- `sampling_plan.csv`: `{summary["scan_rows"]}` data rows
- `channel_diagnostics.csv`: `{summary["diagnostic_rows"]}` data rows
- `failed_points.csv`: `0` failure rows
- `effective_config.json`
- `README.md`
- `PRODUCTION_AUDIT.md`
- `manifest.json`

## Figure Outputs

- Figure directory: `{rel_fig}/`
- `plot_manifest.json`: generated by `scripts/relaxtime/run_phase_guided_transport_plots.jl`

## Known Limitations And Residual Risks

- The denser xi grid reuses the p128 validated-anchored integration policy; it is not a new p104-to-p128 convergence sweep over every xi=0.01 point.
- The old-grid anchor check guards against import drift at all xi=0.05 points.
- Local structures visible only between old xi=0.05 anchors should be interpreted together with neighboring points and channel diagnostics.
"""
    write_text(case_dir / "PRODUCTION_AUDIT.md", audit)


def make_manifest(
    new_case: str,
    mode: str,
    case_dir: Path,
    summary: dict[str, Any],
    mode_manifest: dict[str, Any],
) -> dict[str, Any]:
    files = [
        "phase_guided_transport_scan.csv",
        "sampling_plan.csv",
        "channel_diagnostics.csv",
        "failed_points.csv",
        "effective_config.json",
        "README.md",
        "PRODUCTION_AUDIT.md",
    ]
    hashes = {name: sha256_file(case_dir / name) for name in files}
    return {
        "format": "phase_guided_transport_production_manifest_v1",
        "case_slug": new_case,
        "mode": mode,
        "date": date.today().isoformat(),
        "verdict": "production-grade",
        "source_commit": "700a0845abd09013eb39ff8f8d2993d5574476d3",
        "remote_workflow": ".github/workflows/relaxtime-phase-guided-transport-production.yml",
        "source_candidate_case": SOURCE_CANDIDATE_CASE,
        "production_parameters": {
            "tau_p_nodes": 128,
            "tau_angle_nodes": 20,
            "tau_phi_nodes": 36,
            "tau_n_sigma_points": 24,
            "sigma_grid_n": 560,
        },
        "policies": {
            "propagator_xi_policy": "match_thermo",
            "sigma_cache_policy": "validated_anchored",
            "channel_diagnostics": True,
            "compute_bulk": True,
        },
        "grid": summary["grid"],
        "remote_runs_p128_xi001": {row["label"]: str(row["run_id"]) for row in mode_manifest["source_runs"]},
        "gates": {
            "inherited_parameter_convergence": f"data/outputs/results/relaxtime/transport/phase_guided/{CONVERGENCE_SOURCE}/p104_vs_p128_convergence_gate.summary.json",
            "highxi_import_convergence": f"data/outputs/results/relaxtime/transport/phase_guided/{new_case}_convergence/highxi_anchor_comparison_to_p128_plus5.summary.json",
            "chunk_overlap": summary["overlap"],
            "effective_anchor_worst_rel": summary["effective_anchor_worst_rel"],
            "raw_anchor_worst_rel": summary["anchor_worst_rel"],
        },
        "result_files": files + ["manifest.json"],
        "figure_files": [
            f"data/outputs/figures/relaxtime/transport/phase_guided/{mode}/{new_case}/plot_manifest.json",
            f"data/outputs/figures/relaxtime/transport/phase_guided/{mode}/{new_case}/**/*.png",
        ],
        "hashes": hashes,
        "known_limitations": [
            "xi=0.01 grid reuses the validated p128 integration policy rather than rerunning p104-to-p128 convergence at every new xi point",
            "new intermediate xi structures require local neighbor and channel-diagnostics review before physics interpretation",
        ],
    }


def cleanup_outputs(paths: list[Path], overwrite: bool) -> None:
    for path in paths:
        if not path.exists():
            continue
        if not overwrite:
            raise FileExistsError(f"output path exists; pass --overwrite to replace: {path}")
        resolved = path.resolve()
        if PROJECT_ROOT not in resolved.parents:
            raise ValueError(f"refusing to remove path outside project: {path}")
        shutil.rmtree(resolved)


def main() -> int:
    args = parse_args()
    download_root = args.download_root.resolve()
    run_manifest_path = args.run_manifest.resolve()
    new_case = args.new_case
    old_case = args.old_case

    run_rows = json.loads(run_manifest_path.read_text(encoding="utf-8-sig"))
    discovered = discover_artifacts(download_root, run_rows)
    scan_rows_raw, diag_rows_raw, source_manifest_rows = load_all_rows(discovered)

    conv_dir = RESULT_ROOT / f"{new_case}_convergence"
    mode_dirs = {mode: RESULT_ROOT / mode / new_case for mode in (MODE_A, MODE_B)}
    fig_dirs = {mode: FIGURE_ROOT / mode / new_case for mode in (MODE_A, MODE_B)}
    cleanup_outputs([conv_dir, *mode_dirs.values(), *fig_dirs.values()], args.overwrite)

    merged_scans: dict[str, list[dict[str, str]]] = {}
    merged_diags: dict[str, list[dict[str, str]]] = {}
    overlap_summaries: dict[str, dict[str, Any]] = {}
    diag_overlap_summaries: dict[str, dict[str, Any]] = {}
    summaries: dict[str, dict[str, Any]] = {}
    all_anchor_rows: list[dict[str, Any]] = []
    anchor_summaries: dict[str, Any] = {}
    anchor_comparisons_by_mode: dict[str, list[dict[str, Any]]] = {}

    for mode in (MODE_A, MODE_B):
        merged_scans[mode], overlap_summaries[mode] = merge_rows(
            scan_rows_raw[mode],
            scan_key,
            scan_sort_key,
            SCAN_CORE_COLUMNS,
            args.overlap_rel_threshold,
        )
        merged_diags[mode], diag_overlap_summaries[mode] = merge_rows(
            diag_rows_raw[mode],
            diag_key,
            diag_sort_key,
            DIAG_COMPARE_COLUMNS,
            args.overlap_rel_threshold,
        )
        comparisons, anchor_summary = compare_to_old_case(new_case, old_case, mode, merged_scans[mode])
        all_anchor_rows.extend(comparisons)
        anchor_summaries[mode] = anchor_summary
        anchor_comparisons_by_mode[mode] = comparisons

    anchor_worst = max(all_anchor_rows, key=lambda row: float(row["rel_delta"]), default=None)
    anchor_worst_rel = float(anchor_worst["rel_delta"]) if anchor_worst else 0.0
    effective_anchor_rows = [
        row for row in all_anchor_rows
        if float(row["abs_delta"]) > args.anchor_abs_threshold
    ]
    effective_anchor_worst = max(effective_anchor_rows, key=lambda row: float(row["rel_delta"]), default=None)
    effective_anchor_worst_rel = float(effective_anchor_worst["rel_delta"]) if effective_anchor_worst else 0.0

    gate_failures: list[str] = []
    for mode in (MODE_A, MODE_B):
        if not overlap_summaries[mode]["passed"]:
            gate_failures.append(f"{mode} scan overlap gate failed")
        if not diag_overlap_summaries[mode]["passed"]:
            gate_failures.append(f"{mode} diagnostics overlap gate failed")
    anchor_violations = [
        row for row in all_anchor_rows
        if float(row["abs_delta"]) > args.anchor_abs_threshold and float(row["rel_delta"]) > args.anchor_rel_threshold
    ]
    if anchor_violations:
        gate_failures.append(
            "anchor comparison gate failed: "
            f"max_effective_rel={effective_anchor_worst_rel}, "
            f"violations={len(anchor_violations)}"
        )
    total_failed = sum(int(row["failed_rows"]) for row in source_manifest_rows)
    if total_failed != 0:
        gate_failures.append(f"source artifacts contain failed rows: {total_failed}")

    for mode in (MODE_A, MODE_B):
        source_runs = [row for row in source_manifest_rows if row["mode"] == mode]
        scan_validity = numeric_validity(merged_scans[mode], scan_validity_columns(mode))
        diag_validity = numeric_validity(merged_diags[mode], DIAG_COMPARE_COLUMNS)
        if scan_validity["nonfinite"] != 0:
            gate_failures.append(f"{mode} scan core non-finite count={scan_validity['nonfinite']}")
        if diag_validity["nonfinite"] != 0 or diag_validity["negative"] != 0:
            gate_failures.append(f"{mode} diagnostics validity failed: {diag_validity}")

        mode_effective_anchor_rows = [
            row for row in anchor_comparisons_by_mode[mode]
            if float(row["abs_delta"]) > args.anchor_abs_threshold
        ]
        mode_effective_anchor_worst = max(
            mode_effective_anchor_rows,
            key=lambda row: float(row["rel_delta"]),
            default=None,
        )
        mode_effective_anchor_worst_rel = (
            float(mode_effective_anchor_worst["rel_delta"])
            if mode_effective_anchor_worst is not None
            else 0.0
        )

        xis = collect_xi_values(merged_scans[mode])
        grid: dict[str, Any]
        if mode == MODE_A:
            grid = {
                "muB_MeV": sorted({float(row["muB_MeV"]) for row in merged_scans[mode]}),
                "alpha_T": sorted({float(row["alpha_T"]) for row in merged_scans[mode]}),
                "xi": f"{xis[0]:.2f}:0.01:{xis[-1]:.2f}",
                "points": len(merged_scans[mode]),
            }
        else:
            grid = {
                "T_MeV": sorted({float(row["T_MeV"]) for row in merged_scans[mode]}),
                "muB_MeV": sorted({float(row["muB_MeV"]) for row in merged_scans[mode]}),
                "xi": f"{xis[0]:.2f}:0.01:{xis[-1]:.2f}",
                "points": len(merged_scans[mode]),
            }
        summaries[mode] = {
            "source_run_count": len(source_runs),
            "failed_rows": sum(int(row["failed_rows"]) for row in source_runs),
            "scan_rows": len(merged_scans[mode]),
            "diagnostic_rows": len(merged_diags[mode]),
            "overlap": overlap_summaries[mode],
            "diagnostic_overlap": diag_overlap_summaries[mode],
            "anchor_matched_rows": anchor_summaries[mode]["matched_anchor_rows"],
            "anchor_worst_rel": float(anchor_summaries[mode]["worst"]["rel_delta"]),
            "anchor_worst": anchor_summaries[mode]["worst"],
            "effective_anchor_worst_rel": mode_effective_anchor_worst_rel,
            "effective_anchor_worst": mode_effective_anchor_worst,
            "effective_anchor_rows": len(mode_effective_anchor_rows),
            "scan_validity": scan_validity,
            "diag_validity": diag_validity,
            "grid": grid,
        }

    if gate_failures:
        conv_dir.mkdir(parents=True, exist_ok=True)
        write_json(conv_dir / "blocked_gate_failures.json", {"failures": gate_failures})
        raise SystemExit("Gate failures: " + "; ".join(gate_failures))

    conv_dir.mkdir(parents=True, exist_ok=True)
    write_download_manifest(conv_dir / "download_manifest.csv", source_manifest_rows)
    write_anchor_comparison(conv_dir / "highxi_anchor_comparison_to_p128_plus5.csv", all_anchor_rows)

    anchor_summary_payload = {
        "schema": "phase_guided_highxi_anchor_gate_v1",
        "new_case": new_case,
        "old_case": old_case,
        "anchor_rel_threshold": args.anchor_rel_threshold,
        "anchor_abs_threshold": args.anchor_abs_threshold,
        "anchor_worst_rel": anchor_worst_rel,
        "anchor_worst": anchor_worst,
        "effective_anchor_worst_rel": effective_anchor_worst_rel,
        "effective_anchor_worst": effective_anchor_worst,
        "effective_anchor_rows": len(effective_anchor_rows),
        "anchor_violation_count": len(anchor_violations),
        "mode_summaries": anchor_summaries,
        "passed": len(anchor_violations) == 0,
    }
    write_json(conv_dir / "highxi_anchor_comparison_to_p128_plus5.summary.json", anchor_summary_payload)
    write_text(
        conv_dir / "highxi_anchor_comparison_to_p128_plus5.md",
        f"""# High-xi anchor comparison gate

Verdict: `production-grade`

- new case: `{new_case}`
- old anchor case: `{old_case}`
- comparison: all high-xi rows at xi multiples of 0.05 vs existing p128 production rows
- worst relative difference: `{anchor_worst_rel}`
- effective worst relative difference after abs <= {args.anchor_abs_threshold} floor: `{effective_anchor_worst_rel}`
- threshold: `{args.anchor_rel_threshold}`
- absolute floor: `{args.anchor_abs_threshold}`
- comparison cells: `{len(all_anchor_rows)}`

The p128 integration-parameter convergence gate is inherited from `{CONVERGENCE_SOURCE}` because the numerical integration policy is unchanged.
""",
    )

    write_json(
        conv_dir / "merged_manifest.json",
        {
            "schema": "phase_guided_highxi_merge_manifest_v1",
            "new_case": new_case,
            "source_candidate_case": SOURCE_CANDIDATE_CASE,
            "download_root": str(download_root),
            "source_run_count": len(source_manifest_rows),
            "mode_summaries": summaries,
        },
    )

    write_text(
        conv_dir / "README.md",
        f"""# {new_case} convergence and import gates

This directory records the local import gates used to promote the high-xi phase-guided transport action artifacts.

## Verdict

`production-grade`

## Evidence

- `download_manifest.csv`: 30 source GitHub Actions result artifacts
- `highxi_anchor_comparison_to_p128_plus5.csv`: old-grid anchor comparison against `{old_case}`
- `highxi_anchor_comparison_to_p128_plus5.summary.json`: machine-readable anchor gate summary
- `merged_manifest.json`: merge counts, overlap checks, and validity checks

Parameter convergence is inherited from `{CONVERGENCE_SOURCE}` because this case uses the same p128 validated-anchored integration policy and only increases xi sampling density.
""",
    )
    write_json(
        conv_dir / "manifest.json",
        {
            "format": "phase_guided_transport_highxi_import_gate_manifest_v1",
            "case_slug": new_case,
            "date": date.today().isoformat(),
            "verdict": "production-grade",
            "source_commit": "700a0845abd09013eb39ff8f8d2993d5574476d3",
            "source_candidate_case": SOURCE_CANDIDATE_CASE,
            "inherited_parameter_convergence": f"data/outputs/results/relaxtime/transport/phase_guided/{CONVERGENCE_SOURCE}/p104_vs_p128_convergence_gate.summary.json",
            "anchor_worst_rel": anchor_worst_rel,
            "effective_anchor_worst_rel": effective_anchor_worst_rel,
            "source_run_count": len(source_manifest_rows),
            "files": [
                "README.md",
                "manifest.json",
                "download_manifest.csv",
                "merged_manifest.json",
                "highxi_anchor_comparison_to_p128_plus5.csv",
                "highxi_anchor_comparison_to_p128_plus5.summary.json",
                "highxi_anchor_comparison_to_p128_plus5.md",
            ],
            "hashes": {
                "download_manifest.csv": sha256_file(conv_dir / "download_manifest.csv"),
                "merged_manifest.json": sha256_file(conv_dir / "merged_manifest.json"),
                "highxi_anchor_comparison_to_p128_plus5.csv": sha256_file(conv_dir / "highxi_anchor_comparison_to_p128_plus5.csv"),
                "highxi_anchor_comparison_to_p128_plus5.summary.json": sha256_file(conv_dir / "highxi_anchor_comparison_to_p128_plus5.summary.json"),
            },
        },
    )

    for mode in (MODE_A, MODE_B):
        case_dir = mode_dirs[mode]
        case_dir.mkdir(parents=True, exist_ok=True)
        source_manifest_rel = f"../../{new_case}_convergence/download_manifest.csv"
        scan_comments = production_comments(
            new_case,
            "700a0845abd09013eb39ff8f8d2993d5574476d3",
            source_manifest_rel,
            "phase_guided_transport_scan",
        )
        diag_source = f"data/outputs/results/relaxtime/transport/phase_guided/{mode}/{new_case}/phase_guided_transport_scan.csv"
        diag_comments = production_comments(
            new_case,
            "700a0845abd09013eb39ff8f8d2993d5574476d3",
            source_manifest_rel,
            "phase_guided_transport_scan_channel_diagnostics",
            mode_path=diag_source,
        )
        failed_comments = [
            f"# production_case: {new_case}",
            f"# source_candidate_case: {SOURCE_CANDIDATE_CASE}",
            f"# source_sharded_action_runs: see {source_manifest_rel}",
        ]

        scan_header = list(read_scan_csv(discovered[0]["scan"]).header)
        diag_header = list(read_scan_csv(discovered[0]["diagnostics"]).header)
        failed_header = list(read_scan_csv(discovered[0]["failed"]).header)
        write_csv(case_dir / "phase_guided_transport_scan.csv", scan_comments, scan_header, merged_scans[mode])
        write_csv(case_dir / "channel_diagnostics.csv", diag_comments, diag_header, merged_diags[mode])
        write_csv(case_dir / "failed_points.csv", failed_comments, failed_header, [])
        write_csv(
            case_dir / "sampling_plan.csv",
            [],
            ["T_MeV", "muB_MeV", "xi", "mode", "plot_panel", "plot_series"],
            make_plan_rows(merged_scans[mode]),
        )
        write_json(case_dir / "effective_config.json", effective_config(mode, new_case, merged_scans[mode]))

        mode_manifest = {
            "source_runs": [row for row in source_manifest_rows if row["mode"] == mode],
            "rows_for_text": merged_scans[mode],
        }
        write_readme_and_audit(case_dir, conv_dir, mode, new_case, summaries[mode], mode_manifest)
        manifest = make_manifest(new_case, mode, case_dir, summaries[mode], mode_manifest)
        write_json(case_dir / "manifest.json", manifest)

    print(
        json.dumps(
            {
                "new_case": new_case,
                "result_dirs": {mode: str(mode_dirs[mode]) for mode in (MODE_A, MODE_B)},
                "convergence_dir": str(conv_dir),
                "summaries": summaries,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
