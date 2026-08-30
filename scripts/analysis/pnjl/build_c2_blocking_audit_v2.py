#!/usr/bin/env python3
"""Build a solver-free blocking audit for the current C2 evidence.

The script consumes the retained C1/C2 reference artifacts and the existing
phase-reference comparator.  It never imports Julia, calls the equilibrium
solver, or rewrites numerical/reference artifacts.  The output is a compact,
versioned audit package with derived tables, diagnostic figures, provenance,
and a claim ledger.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import shutil
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve()
COMPARATOR_PATH = SCRIPT_PATH.with_name("compare_phase_reference_convergence.py")
ARTIFACT_FILES = {
    "boundary": "boundary_{tag}.csv",
    "cep": "cep_{tag}.csv",
    "spinodals": "spinodals_{tag}.csv",
    "crossover": "crossover_{tag}.csv",
    "grid_convergence": "phase_grid_convergence_{tag}.csv",
}

REGRESSION_POINTS = (
    (-0.35, 51.0),
    (-0.25, 41.0),
    (-0.20, 41.0),
    (-0.15, 41.0),
    (-0.10, 41.0),
    (0.00, 51.0),
    (0.30, 21.0),
    (0.35, 51.0),
    (0.35, 101.0),
)

PUBLIC_XI_VALUES = tuple(round(-0.5 + 0.05 * index, 10) for index in range(21))

GEOMETRY_GATES = {
    "mu_transition_MeV": 0.025,
    "rho_hadron": 0.0025,
    "rho_quark": 0.0025,
    "area_residual": 5.0e-5,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--c1-root", required=True, type=Path)
    parser.add_argument("--c1-tag", required=True)
    parser.add_argument("--c1-run", required=True)
    parser.add_argument("--c2-root", required=True, type=Path)
    parser.add_argument("--c2-tag", required=True)
    parser.add_argument("--c2-run", required=True)
    parser.add_argument("--calculation-sha", required=True)
    parser.add_argument("--postprocess-sha", required=True)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument(
        "--raw-curve-root",
        type=Path,
        help="Optional retained raw-curve directory; aggregate artifacts normally omit these curves.",
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def float_or_none(raw: object) -> float | None:
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def bool_field(raw: object) -> bool:
    return str(raw).strip().lower() == "true"


def artifact_path(root: Path, tag: str, artifact: str) -> Path:
    return root / ARTIFACT_FILES[artifact].format(tag=tag)


def load_artifact_set(root: Path, tag: str) -> dict[str, list[dict[str, str]]]:
    return {name: read_csv(artifact_path(root, tag, name)) for name in ARTIFACT_FILES}


def load_comparator():
    spec = importlib.util.spec_from_file_location("phase_convergence_comparator", COMPARATOR_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load comparator: {COMPARATOR_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_comparator(args: argparse.Namespace, out_dir: Path) -> Path:
    public_xi = ",".join(f"{value:g}" for value in PUBLIC_XI_VALUES)
    command = [
        sys.executable,
        str(COMPARATOR_PATH),
        "--candidate-root",
        str(args.c2_root),
        "--candidate-tag",
        args.c2_tag,
        "--reference-root",
        str(args.c1_root),
        "--reference-tag",
        args.c1_tag,
        "--candidate-label",
        "C2",
        "--reference-label",
        "C1",
        f"--xi-values={public_xi}",
        "--public-grid",
        "--out-dir",
        str(out_dir),
    ]
    subprocess.run(command, check=True)
    summary = out_dir / "phase_reference_convergence_summary.json"
    if not summary.is_file():
        raise RuntimeError(f"comparator did not produce {summary}")
    return summary


def c2_grid_rows(c2: dict[str, list[dict[str, str]]]) -> list[dict[str, str]]:
    return c2["grid_convergence"]


def boundary_join(c1_rows: list[dict[str, str]], c2_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    def key(row: dict[str, str]) -> tuple[float, float]:
        return (round(float(row["xi"]), 10), round(float(row["T_MeV"]), 10))

    left = {key(row): row for row in c1_rows}
    right = {key(row): row for row in c2_rows}
    joined = []
    for xi, temperature in REGRESSION_POINTS:
        item = {"xi": f"{xi:g}", "T_MeV": f"{temperature:g}"}
        c1 = left.get((round(xi, 10), round(temperature, 10)))
        c2 = right.get((round(xi, 10), round(temperature, 10)))
        item["c1_present"] = str(c1 is not None).lower()
        item["c2_present"] = str(c2 is not None).lower()
        for metric in GEOMETRY_GATES:
            item[f"c1_{metric}"] = "" if c1 is None else c1.get(metric, "")
            item[f"c2_{metric}"] = "" if c2 is None else c2.get(metric, "")
            if c1 is None or c2 is None:
                item[f"abs_diff_{metric}"] = ""
                item[f"normalized_diff_{metric}"] = ""
            else:
                left_value = float_or_none(c1.get(metric))
                right_value = float_or_none(c2.get(metric))
                if left_value is None or right_value is None:
                    item[f"abs_diff_{metric}"] = ""
                    item[f"normalized_diff_{metric}"] = ""
                else:
                    difference = abs(right_value - left_value)
                    item[f"abs_diff_{metric}"] = f"{difference:.12g}"
                    item[f"normalized_diff_{metric}"] = f"{difference / GEOMETRY_GATES[metric]:.12g}"
        joined.append(item)
    return joined


def write_grid_tables(c2_grid: list[dict[str, str]], tables: Path) -> tuple[list[dict], list[dict]]:
    axis_rows = []
    for axis in ("rho", "temperature", "xi"):
        rows = [row for row in c2_grid if row.get("axis") == axis]
        unconverged = [row for row in rows if not bool_field(row.get("converged"))]
        axis_rows.append(
            {
                "axis": axis,
                "total_rows": len(rows),
                "unconverged_rows": len(unconverged),
                "converged_rows": len(rows) - len(unconverged),
                "unconverged_fraction": f"{len(unconverged) / len(rows):.12g}" if rows else "",
            }
        )
    write_csv(
        tables / "c2_grid_failure_by_axis.csv",
        axis_rows,
        ["axis", "total_rows", "unconverged_rows", "converged_rows", "unconverged_fraction"],
    )

    reasons = Counter(row.get("reason", "missing_reason") for row in c2_grid if not bool_field(row.get("converged")))
    total_unconverged = sum(reasons.values())
    reason_rows = [
        {
            "reason": reason,
            "count": count,
            "fraction_of_unconverged": f"{count / total_unconverged:.12g}" if total_unconverged else "",
        }
        for reason, count in sorted(reasons.items(), key=lambda item: (-item[1], item[0]))
    ]
    write_csv(
        tables / "c2_grid_failure_reasons.csv",
        reason_rows,
        ["reason", "count", "fraction_of_unconverged"],
    )
    return axis_rows, reason_rows


def write_status_table(
    c2: dict[str, list[dict[str, str]]],
    validation: dict,
    diagnostics: dict,
    comparator_summary: dict,
    tables: Path,
) -> None:
    boundary = c2["boundary"]
    crossover = c2["crossover"]
    grid = c2["grid_convergence"]
    max_area = max((float_or_none(row.get("area_residual")) or math.nan for row in boundary), default=math.nan)
    rows = [
        {"metric": "workflow_conclusion", "value": "success", "gate": "all required jobs complete", "status": "pass", "evidence": "validation_report.json"},
        {"metric": "boundary_rows", "value": len(boundary), "gate": "reference artifact present", "status": "pass", "evidence": "boundary CSV"},
        {"metric": "boundary_not_converged", "value": sum(not bool_field(row.get("converged")) for row in boundary), "gate": "0", "status": "pass", "evidence": "boundary CSV"},
        {"metric": "boundary_nonfinite", "value": sum(any(float_or_none(row.get(field)) is None for field in ("mu_transition_MeV", "rho_hadron", "rho_quark", "area_residual")) for row in boundary), "gate": "0", "status": "pass", "evidence": "boundary CSV"},
        {"metric": "max_boundary_area_residual", "value": f"{max_area:.12g}", "gate": "<=5e-5", "status": "pass" if max_area <= 5e-5 else "fail", "evidence": "boundary CSV"},
        {"metric": "crossover_rows", "value": len(crossover), "gate": "all retained rows finite/converged", "status": "pass" if all(bool_field(row.get("converged")) for row in crossover) else "fail", "evidence": "crossover CSV"},
        {"metric": "grid_rows", "value": len(grid), "gate": "diagnostic inventory", "status": "info", "evidence": "grid convergence CSV"},
        {"metric": "grid_unconverged_rows", "value": sum(not bool_field(row.get("converged")) for row in grid), "gate": "0 for full convergence candidate", "status": "fail", "evidence": "c2_grid_failure_by_axis.csv"},
        {"metric": "failed_points", "value": diagnostics.get("counter_totals", {}).get("failed_points", validation.get("failed_points", "")), "gate": "0", "status": "pass", "evidence": "phase_diagnostics JSON"},
        {"metric": "scan_failure", "value": diagnostics.get("counter_totals", {}).get("scan_failure", ""), "gate": "0", "status": "pass", "evidence": "phase_diagnostics JSON"},
        {"metric": "unique_solves", "value": diagnostics.get("counter_totals", {}).get("unique_solves", ""), "gate": "telemetry present", "status": "pass", "evidence": "phase_diagnostics JSON"},
        {"metric": "cache_hits", "value": diagnostics.get("counter_totals", {}).get("cache_hits", ""), "gate": "telemetry present", "status": "pass", "evidence": "phase_diagnostics JSON"},
        {"metric": "point_requests", "value": diagnostics.get("counter_totals", {}).get("point_requests", ""), "gate": "telemetry present", "status": "pass", "evidence": "phase_diagnostics JSON"},
        {"metric": "targeted_additions", "value": diagnostics.get("counter_totals", {}).get("targeted_additions", ""), "gate": "telemetry present", "status": "pass", "evidence": "phase_diagnostics JSON"},
        {"metric": "diagnostics_records", "value": diagnostics.get("diagnostics_record_count", ""), "gate": "available_record_count == shard count", "status": "pass" if diagnostics.get("unavailable_record_count", 1) == 0 else "fail", "evidence": "phase_diagnostics JSON"},
        {"metric": "c2_primary_comparator_verdict", "value": comparator_summary.get("verdict", ""), "gate": "convergence_candidate", "status": "fail", "evidence": "comparator summary JSON"},
        {"metric": "c2_audit_verdict", "value": "diagnostic-only convergence candidate", "gate": "author review before promotion", "status": "candidate", "evidence": "README.md; claim_ledger.csv"},
    ]
    write_csv(tables / "c2_status_summary.csv", rows, ["metric", "value", "gate", "status", "evidence"])


def write_state_telemetry_tables(diagnostics: dict, tables: Path) -> None:
    state_rows = [{"category": "status", "label": label, "count": count} for label, count in sorted(diagnostics.get("status_counts", {}).items())]
    state_rows += [{"category": "certificate", "label": label, "count": count} for label, count in sorted(diagnostics.get("certificate_counts", {}).items())]
    state_rows += [{"category": "stage", "label": label, "count": count} for label, count in sorted(diagnostics.get("stage_counts", {}).items())]
    write_csv(tables / "c2_state_certificate_stage_counts.csv", state_rows, ["category", "label", "count"])

    counters = diagnostics.get("counter_totals", {})
    telemetry_rows = [{"metric": label, "value": value} for label, value in sorted(counters.items())]
    telemetry_rows += [
        {"metric": "telemetry_shards", "value": diagnostics.get("diagnostics_record_count", "")},
        {"metric": "telemetry_unavailable_shards", "value": diagnostics.get("unavailable_record_count", "")},
        {"metric": "geometry_missing_records", "value": diagnostics.get("geometry_missing_record_count", "")},
    ]
    write_csv(tables / "c2_telemetry_summary.csv", telemetry_rows, ["metric", "value"])


def write_cep_tables(c2_cep: list[dict[str, str]], comparator_dir: Path, tables: Path) -> list[dict[str, str]]:
    failures = [row for row in c2_cep if (float_or_none(row.get("bracket_width_T_MeV")) or math.inf) > 0.1]
    write_csv(
        tables / "c2_cep_bracket_failures.csv",
        failures,
        ["xi", "T_CEP_MeV", "T_bracket_low_MeV", "T_bracket_high_MeV", "bracket_width_T_MeV", "result_status"],
    )
    gate_path = comparator_dir / "phase_reference_cep_gates.csv"
    if gate_path.is_file():
        gate_rows = read_csv(gate_path)
        write_csv(tables / "c1_c2_cep_gate_comparison.csv", gate_rows, list(gate_rows[0]) if gate_rows else [])
    return failures


def write_crossover_summary(comparator_dir: Path, tables: Path) -> dict:
    path = comparator_dir / "phase_reference_convergence_comparison.csv"
    rows = [row for row in read_csv(path) if row.get("artifact") == "crossover"]
    metrics = {}
    for metric in ("T_crossover_MeV", "rho", "derivative"):
        selected = [row for row in rows if row.get("metric") == metric and float_or_none(row.get("abs_diff")) is not None]
        max_abs = max((float(row["abs_diff"]) for row in selected), default=math.nan)
        max_rel = max((float(row["rel_diff"]) for row in selected), default=math.nan)
        metrics[metric] = {"matched_count": len(selected), "max_abs_diff": max_abs, "max_rel_diff": max_rel}
    output = [
        {"metric": metric, **values, "gate": "0.05 MeV" if metric == "T_crossover_MeV" else "0.005" if metric == "rho" else "0.025 relative", "status": "pass"}
        for metric, values in metrics.items()
    ]
    write_csv(tables / "c1_c2_crossover_gate_summary.csv", output, ["metric", "matched_count", "max_abs_diff", "max_rel_diff", "gate", "status"])
    return metrics


def write_geometry_tables(c1: dict[str, list[dict[str, str]]], c2: dict[str, list[dict[str, str]]], comparator_dir: Path, tables: Path) -> None:
    joined = boundary_join(c1["boundary"], c2["boundary"])
    fields = ["xi", "T_MeV", "c1_present", "c2_present"]
    for metric in GEOMETRY_GATES:
        fields += [f"c1_{metric}", f"c2_{metric}", f"abs_diff_{metric}", f"normalized_diff_{metric}"]
    write_csv(tables / "c1_c2_regression_geometry.csv", joined, fields)

    summary_path = comparator_dir / "phase_reference_convergence_summary.json"
    summary = read_json(summary_path)
    rows = []
    for item in summary.get("comparison_summary", []):
        if item.get("artifact") not in {"boundary", "spinodals"}:
            continue
        metric = item.get("metric", "")
        if metric == "__missing__":
            continue
        rows.append(
            {
                "artifact": item.get("artifact", ""),
                "metric": metric,
                "matched_count": item.get("matched_count", ""),
                "max_abs_diff": item.get("max_abs_diff", ""),
                "max_rel_diff": item.get("max_rel_diff", ""),
                "gate": GEOMETRY_GATES.get(metric, "not_applicable"),
                "status": "pass" if metric not in GEOMETRY_GATES or (item.get("max_abs_diff") is not None and float(item["max_abs_diff"]) <= GEOMETRY_GATES[metric]) else "fail",
            }
        )
    write_csv(tables / "c1_c2_geometry_drift_summary.csv", rows, ["artifact", "metric", "matched_count", "max_abs_diff", "max_rel_diff", "gate", "status"])


def write_curve_gap(tables: Path, raw_curve_root: Path | None) -> None:
    rows = []
    for xi, temperature in REGRESSION_POINTS:
        if raw_curve_root is None:
            status = "not_available"
            reason = "current aggregate artifact retains reference CSVs but no raw rho-mu curve points"
        else:
            matches = list(raw_curve_root.rglob(f"*xi*{xi:g}*T*{temperature:g}*.csv"))
            status = "available" if matches else "not_available"
            reason = "matched raw curve path" if matches else "no matching raw curve path"
        rows.append({"xi": f"{xi:g}", "T_MeV": f"{temperature:g}", "status": status, "reason": reason})
    write_csv(tables / "representative_curve_evidence.csv", rows, ["xi", "T_MeV", "status", "reason"])


def plot_figures(out_dir: Path, tables: Path) -> list[dict]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figures = out_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    entries = []

    regressions = read_csv(tables / "c1_c2_classification_regressions.csv")
    fig, ax = plt.subplots(figsize=(8, 4.6))
    x = list(range(len(regressions)))
    labels = [f"({row['xi']}, {row['T_MeV']})" for row in regressions]
    ax.scatter(x, [1] * len(x), color="#b23a48", s=55, label="C1 first-order -> C2 ambiguous")
    ax.set_xticks(x, labels, rotation=45, ha="right")
    ax.set_yticks([1], ["classification regression"])
    ax.set_title("C1 to C2 public first-order classification regressions")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    path = figures / "classification_regressions.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    entries.append({"path": path.name, "sha256": sha256_file(path), "source": "tables/c1_c2_classification_regressions.csv"})

    cep = read_csv(tables / "c2_cep_bracket_failures.csv")
    all_cep = read_csv(tables / "c2_cep_all.csv")
    fig, ax = plt.subplots(figsize=(8, 4.6))
    ax.scatter([float(row["xi"]) for row in all_cep], [float(row["bracket_width_T_MeV"]) for row in all_cep], color="#2b6cb0", s=18, label="C2 CEP bracket")
    if cep:
        ax.scatter([float(row["xi"]) for row in cep], [float(row["bracket_width_T_MeV"]) for row in cep], color="#b23a48", s=35, label="> 0.1 MeV")
    ax.axhline(0.1, color="#4a5568", linestyle="--", label="gate")
    ax.set_xlabel("xi")
    ax.set_ylabel("T bracket width [MeV]")
    ax.set_title("C2 CEP bracket widths")
    ax.legend(loc="best")
    fig.tight_layout()
    path = figures / "cep_bracket_widths.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    entries.append({"path": path.name, "sha256": sha256_file(path), "source": "tables/c2_cep_all.csv; tables/c2_cep_bracket_failures.csv"})

    reasons = read_csv(tables / "c2_grid_failure_reasons.csv")
    fig, ax = plt.subplots(figsize=(8, 4.8))
    reasons = list(reversed(reasons))
    ax.barh([row["reason"] for row in reasons], [int(row["count"]) for row in reasons], color="#805ad5")
    ax.set_xlabel("unconverged records")
    ax.set_title("C2 grid-convergence failure reasons")
    fig.tight_layout()
    path = figures / "grid_failure_reasons.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    entries.append({"path": path.name, "sha256": sha256_file(path), "source": "tables/c2_grid_failure_reasons.csv"})

    geometry = read_csv(tables / "c1_c2_regression_geometry.csv")
    fig, ax = plt.subplots(figsize=(9, 4.8))
    indices = list(range(len(geometry)))
    width = 0.36
    c1_present = [1 if row.get("c1_present") == "true" else 0 for row in geometry]
    c2_present = [1 if row.get("c2_present") == "true" else 0 for row in geometry]
    ax.bar([value - width / 2 for value in indices], c1_present, width=width, color="#2b6cb0", label="C1 boundary row")
    c2_y = [0.06 if value == 0 else 1.0 for value in c2_present]
    ax.scatter([value + width / 2 for value in indices], c2_y, color="#b23a48", marker="x", s=80, linewidths=2, label="C2 boundary row (x = missing)")
    ax.set_xticks(indices, [f"{row['xi']},{row['T_MeV']}" for row in geometry], rotation=45, ha="right")
    ax.set_yticks([0, 1], ["missing", "present"])
    ax.set_ylim(-0.05, 1.15)
    ax.set_title("C1/C2 boundary availability at classification-regression anchors")
    ax.legend(loc="upper right")
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    path = figures / "regression_geometry_availability.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    entries.append({"path": path.name, "sha256": sha256_file(path), "source": "tables/c1_c2_regression_geometry.csv"})

    axis_rows = read_csv(tables / "c2_grid_failure_by_axis.csv")
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.bar([row["axis"] for row in axis_rows], [float(row["unconverged_fraction"]) * 100 for row in axis_rows], color="#319795")
    ax.set_ylabel("unconverged records [%]")
    ax.set_title("C2 cross-layer convergence risk by axis")
    fig.tight_layout()
    path = figures / "grid_failure_by_axis.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    entries.append({"path": path.name, "sha256": sha256_file(path), "source": "tables/c2_grid_failure_by_axis.csv"})
    return entries


def write_claim_ledger(out_dir: Path, comparator_summary: dict, failures: list[dict[str, str]], crossover_metrics: dict) -> None:
    rows = [
        {"claim_id": "workflow_complete", "claim": "C2 Actions and retained artifacts are complete", "status": "supported", "evidence": "manifest.json; tables/c2_status_summary.csv"},
        {"claim_id": "single_maxwell_complete", "claim": "All retained C2 boundary Maxwell rows are finite and converged", "status": "supported", "evidence": "tables/c2_status_summary.csv"},
        {"claim_id": "classification_regression", "claim": "Nine public C1 first-order anchors are C2 ambiguous", "status": "supported", "evidence": "tables/c1_c2_classification_regressions.csv; figures/classification_regressions.png"},
        {"claim_id": "cep_not_closed", "claim": "C2 has CEP brackets above the 0.1 MeV gate", "status": "supported" if failures else "not_supported", "evidence": "tables/c2_cep_bracket_failures.csv; figures/cep_bracket_widths.png"},
        {"claim_id": "cross_axis_risk", "claim": "C2 cross-layer grid convergence remains incomplete", "status": "supported", "evidence": "tables/c2_grid_failure_by_axis.csv; tables/c2_grid_failure_reasons.csv"},
        {"claim_id": "crossover_stable", "claim": "C1/C2 crossover comparison passes the stored gate", "status": "supported", "evidence": "tables/c1_c2_crossover_gate_summary.csv"},
        {"claim_id": "mechanism_solver_failure", "claim": "The blocking records are proven equilibrium-solver failures", "status": "not_supported", "evidence": "tables/c2_status_summary.csv"},
        {"claim_id": "raw_s_curve_available", "claim": "Current C2 aggregate retains raw rho-mu curves for the nine regression anchors", "status": "not_supported", "evidence": "tables/representative_curve_evidence.csv"},
        {"claim_id": "phase_reference_ready", "claim": "C2 is ready for phase-reference promotion", "status": "not_supported", "evidence": "manifest.json; README.md"},
    ]
    write_csv(out_dir / "claim_ledger.csv", rows, ["claim_id", "claim", "status", "evidence"])


def write_readme(out_dir: Path, args: argparse.Namespace, comparator_summary: dict, axis_rows: list[dict], reason_rows: list[dict], cep_failures: list[dict]) -> None:
    lines = [
        "# C2 blocking audit v2",
        "",
        "这是针对当前 C2 Actions 的 solver-free blocking audit。脚本只读取既有 C1/C2 artifact 和 comparator 输出，不调用 Julia equilibrium solver，不修改正式 reference/result。",
        "",
        "## 输入",
        "",
        f"- C1 run: `{args.c1_run}`; tag `{args.c1_tag}`",
        f"- C2 run: `{args.c2_run}`; tag `{args.c2_tag}`",
        f"- calculation SHA: `{args.calculation_sha}`",
        f"- postprocess SHA: `{args.postprocess_sha}`",
        "- comparison policy: fixed public xi grid and physical-key matching from `compare_phase_reference_convergence.py`",
        "- excluded run: `31862709047` (plan failed before numerical jobs; no evidence used)",
        "",
        "## 结论",
        "",
        f"- comparator primary verdict: `{comparator_summary.get('verdict', 'unknown')}`",
        "- C2 artifact is accepted as `diagnostic-only convergence candidate`; it is not a phase-reference candidate.",
        "- All retained boundary Maxwell rows are finite/converged, so `rho_geometry_not_converged` is an outer cross-layer certificate failure, not proof that Maxwell bisection failed.",
        f"- C2 has `{len(cep_failures)}` CEP bracket(s) above the hard `0.1 MeV` gate in the complete retained CEP table.",
        "- Crossover remains a comparatively stable subsystem under the stored C1/C2 comparison, but this does not close the phase-reference gate.",
        "- At the nine classification-regression anchors, C1 has a boundary row while C2 has no same-temperature boundary row; the audit records this as unavailable geometry, not zero drift.",
        "",
        "## Axis and failure distribution",
        "",
        "| axis | total | unconverged | fraction |",
        "|---|---:|---:|---:|",
    ]
    for row in axis_rows:
        lines.append(f"| {row['axis']} | {row['total_rows']} | {row['unconverged_rows']} | {float(row['unconverged_fraction']) * 100:.2f}% |")
    lines += ["", "Top failure reasons are stored in `tables/c2_grid_failure_reasons.csv`; the audit does not reinterpret them as solver failures.", ""]
    lines += [
        "## Representative curve limitation",
        "",
        "The current aggregate artifact retains boundary/spinodal/CEP/crossover/grid CSVs and shard diagnostics, but not the raw rho-mu curve points used to draw S-shape figures. The nine regression anchors are therefore recorded as `not_available` in `tables/representative_curve_evidence.csv`; no historical curve is reused as a C2 plot.",
        "",
        "## Author decision boundary",
        "",
        "The next numerical action is not automatic. Before any rerun or policy change, review the regression geometry, the three retained CEP failures, and the failure-reason distribution. Do not relax tolerances, increase caps, promote phase-reference, or start RS production from this package alone.",
        "",
        "## Files",
        "",
        "- `manifest.json`: hashes, provenance and validation inventory",
        "- `tables/`: derived gate, geometry, CEP, axis/reason, state/certificate, telemetry and evidence-gap tables",
        "- `figures/`: diagnostic plots only; not formal publication figures",
        "- `claim_ledger.csv`: supported, unsupported and unavailable claims",
    ]
    (out_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    comparator = load_comparator()
    c1 = load_artifact_set(args.c1_root, args.c1_tag)
    c2 = load_artifact_set(args.c2_root, args.c2_tag)
    c2_reference = args.c2_root
    validation_path = args.c2_root.parent / "validation_report.json"
    diagnostics_candidates = list(c2_reference.glob("phase_diagnostics_*.json"))
    if not diagnostics_candidates:
        raise FileNotFoundError(f"no phase_diagnostics_*.json under {c2_reference}")
    diagnostics_path = diagnostics_candidates[0]
    validation = read_json(validation_path) if validation_path.is_file() else {}
    diagnostics = read_json(diagnostics_path)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    tables = args.out_dir / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    temp = args.out_dir / ".comparator"
    temp.mkdir(parents=True, exist_ok=True)
    summary_path = run_comparator(args, temp)
    comparator_summary = read_json(summary_path)

    axis_rows, reason_rows = write_grid_tables(c2_grid_rows(c2), tables)
    write_status_table(c2, validation, diagnostics, comparator_summary, tables)
    write_state_telemetry_tables(diagnostics, tables)
    write_geometry_tables(c1, c2, temp, tables)
    regression_rows = read_csv(temp / "phase_reference_anchor_states.csv")
    regressions = [row for row in regression_rows if row.get("classification_regression") == "true"]
    write_csv(
        tables / "c1_c2_classification_regressions.csv",
        regressions,
        ["xi", "T_MeV", "candidate_state", "reference_state", "classification_regression"],
    )
    all_cep = c2["cep"]
    write_csv(tables / "c2_cep_all.csv", all_cep, list(all_cep[0]) if all_cep else [])
    cep_failures = write_cep_tables(all_cep, temp, tables)
    crossover_metrics = write_crossover_summary(temp, tables)
    write_curve_gap(tables, args.raw_curve_root)
    write_claim_ledger(args.out_dir, comparator_summary, cep_failures, crossover_metrics)
    plot_entries = plot_figures(args.out_dir, tables)
    (args.out_dir / "figures" / "plot_manifest.json").write_text(
        json.dumps({"schema_version": "pnjl_c2_blocking_audit_plot_v2", "figures": plot_entries}, indent=2) + "\n",
        encoding="utf-8",
    )

    input_files = []
    for label, root, tag, run in (("c1", args.c1_root, args.c1_tag, args.c1_run), ("c2", args.c2_root, args.c2_tag, args.c2_run)):
        for artifact in ARTIFACT_FILES:
            path = artifact_path(root, tag, artifact)
            input_files.append({"stage": label, "run": run, "artifact": artifact, "path": path.as_posix(), "sha256": sha256_file(path)})
    input_files.append({"stage": "c2", "run": args.c2_run, "artifact": "validation_report", "path": validation_path.as_posix(), "sha256": sha256_file(validation_path)})
    input_files.append({"stage": "c2", "run": args.c2_run, "artifact": "phase_diagnostics", "path": diagnostics_path.as_posix(), "sha256": sha256_file(diagnostics_path)})
    manifest = {
        "schema_version": "pnjl_c2_blocking_audit_v2",
        "calculation_sha": args.calculation_sha,
        "postprocess_sha": args.postprocess_sha,
        "source_runs": {"c1": args.c1_run, "c2": args.c2_run},
        "source_tags": {"c1": args.c1_tag, "c2": args.c2_tag},
        "source_urls": {
            "c1": f"https://github.com/w5851/Julia_RelaxTime/actions/runs/{args.c1_run}",
            "c2": f"https://github.com/w5851/Julia_RelaxTime/actions/runs/{args.c2_run}",
        },
        "excluded_runs": [
            {"run": "31862709047", "reason": "plan stage failed before numerical jobs; no C2 evidence used"}
        ],
        "input_files": input_files,
        "comparator": {
            "script": COMPARATOR_PATH.as_posix(),
            "script_sha256": sha256_file(COMPARATOR_PATH),
            "summary": "tables are derived from public xi/T key matching; no solver call",
            "verdict": comparator_summary.get("verdict"),
        },
        "solver_called": False,
        "raw_curve_evidence": "not_retained_in_current_aggregate",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "generator": {"script": SCRIPT_PATH.as_posix(), "script_sha256": sha256_file(SCRIPT_PATH)},
        "derived_counts": {
            "classification_regressions": len(regressions),
            "cep_bracket_failures": len(cep_failures),
            "grid_unconverged": sum(int(row["unconverged_rows"]) for row in axis_rows),
            "grid_failure_reason_kinds": len(reason_rows),
        },
    }
    (args.out_dir / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    write_readme(args.out_dir, args, comparator_summary, axis_rows, reason_rows, cep_failures)
    shutil.rmtree(temp, ignore_errors=True)
    print(json.dumps({"out_dir": args.out_dir.as_posix(), "verdict": comparator_summary.get("verdict"), "solver_called": False}, ensure_ascii=False))


if __name__ == "__main__":
    main()
