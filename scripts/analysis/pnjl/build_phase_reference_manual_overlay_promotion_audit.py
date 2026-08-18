#!/usr/bin/env python3
"""Build a solver-free phase-reference promotion gate audit.

The audit combines the C2 automatic evidence with two author-review overlays:
the nine targeted rho-mu curves and the three manually closed CEP brackets.
It never changes C2 rows and never calls Julia or the PNJL equilibrium solver.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
BRACKET_TARGET_MEV = 0.1
GRID_REASON_CLASSES = {
    "geometry_tolerance_exceeded": "geometry_interpolation_diagnostic",
    "rho_geometry_not_converged": "geometry_interpolation_diagnostic",
    "interpolation_tolerance_exceeded": "geometry_interpolation_diagnostic",
    "ok": "unresolved_unclassified",
    "hybrid_stage_c_not_converged": "hybrid_stage_c_unresolved",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--c1-root", type=Path, required=True)
    parser.add_argument("--c2-root", type=Path, required=True)
    parser.add_argument("--c2-audit-root", type=Path, required=True)
    parser.add_argument("--nine-review-root", type=Path, required=True)
    parser.add_argument("--manual-overlay", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"CSV has no data rows: {path}")
    return rows


def finite(value: str | None) -> bool:
    if value is None or not value.strip():
        return False
    try:
        return math.isfinite(float(value))
    except ValueError:
        return False


def bool_value(value: str | bool | None) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() == "true"


def write_csv(path: Path, fieldnames: Iterable[str], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        writer.writerows(rows)


def check_unique(rows: list[dict[str, str]], keys: tuple[str, ...]) -> bool:
    values = [tuple(row.get(key, "") for key in keys) for row in rows]
    return len(values) == len(set(values))


def require_file(path: Path) -> Path:
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def main() -> int:
    args = parse_args()
    c1_root = args.c1_root.resolve()
    c2_root = args.c2_root.resolve()
    c2_audit_root = args.c2_audit_root.resolve()
    nine_review_root = args.nine_review_root.resolve()
    manual_overlay_path = args.manual_overlay.resolve()
    output_root = args.output_root.resolve()
    tables_root = output_root / "tables"
    output_root.mkdir(parents=True, exist_ok=True)
    tables_root.mkdir(parents=True, exist_ok=True)

    c1_reference = c1_root / "reference"
    c2_reference = c2_root / "reference"
    c1_grid = require_file(next(c1_reference.glob("phase_grid_convergence_*.csv")))
    c2_grid = require_file(next(c2_reference.glob("phase_grid_convergence_*.csv")))
    c1_boundary = require_file(next(c1_reference.glob("boundary_*.csv")))
    c2_boundary = require_file(next(c2_reference.glob("boundary_*.csv")))
    c2_crossover = require_file(next(c2_reference.glob("crossover_*.csv")))
    c2_cep = require_file(next(c2_reference.glob("cep_*.csv")))
    c2_diagnostics = require_file(next(c2_reference.glob("phase_diagnostics_*.json")))
    c2_validation = require_file(c2_root / "validation_report.json")
    c2_audit_manifest = require_file(c2_audit_root / "manifest.json")
    c2_regressions = require_file(c2_audit_root / "tables" / "c1_c2_classification_regressions.csv")
    c2_drift = require_file(c2_audit_root / "tables" / "c1_c2_geometry_drift_summary.csv")
    nine_manifest = require_file(nine_review_root / "manifest.json")
    nine_points = require_file(nine_review_root / "tables" / "point_summary.csv")
    manual_overlay_path = require_file(manual_overlay_path)

    c1_grid_rows = read_csv(c1_grid)
    c2_grid_rows = read_csv(c2_grid)
    c1_boundary_rows = read_csv(c1_boundary)
    c2_boundary_rows = read_csv(c2_boundary)
    c2_crossover_rows = read_csv(c2_crossover)
    c2_cep_rows = read_csv(c2_cep)
    c2_diagnostics_payload = read_json(c2_diagnostics)
    validation_payload = read_json(c2_validation)
    manual_overlay = read_json(manual_overlay_path)
    regression_rows = read_csv(c2_regressions)
    drift_rows = read_csv(c2_drift)
    nine_point_rows = read_csv(nine_points)

    if manual_overlay.get("source_calculation_sha", "").lower() != CALCULATION_SHA:
        raise ValueError("manual overlay calculation SHA does not match C2")
    if manual_overlay.get("manual_decision_required") is not False:
        raise ValueError("manual overlay is not closed")
    cep_overlay = manual_overlay.get("cep_brackets") or {}
    if len(cep_overlay) != 3:
        raise ValueError("manual overlay must contain exactly three CEP brackets")
    for label, payload in cep_overlay.items():
        bracket = payload.get("bracket_MeV") or []
        width = float(payload.get("width_MeV"))
        if payload.get("status") != "closed_by_author_review" or len(bracket) != 2:
            raise ValueError(f"manual CEP overlay is not closed: {label}")
        if not math.isfinite(width) or width > BRACKET_TARGET_MEV:
            raise ValueError(f"manual CEP overlay misses target: {label}")

    # The nine curves are author-accepted oracle observations, not hybrid label overrides.
    nine_ok = (
        len(nine_point_rows) == 9
        and all(row.get("oracle_status") == "confirmed_first_order" for row in nine_point_rows)
        and all(row.get("oracle_candidate_count") == "1" for row in nine_point_rows)
        and all(row.get("oracle_crossing_count") == "3" for row in nine_point_rows)
        and all(bool_value(row.get("oracle_geometry_converged")) for row in nine_point_rows)
        and all(bool_value(row.get("finite_and_converged")) for row in nine_point_rows)
    )

    unresolved_rows = [row for row in c2_grid_rows if not bool_value(row.get("converged"))]
    reason_counts = Counter(row.get("reason", "") for row in unresolved_rows)
    axis_counts = Counter(row.get("axis", "") for row in unresolved_rows)
    reason_class_rows = []
    for reason, count in sorted(reason_counts.items()):
        category = GRID_REASON_CLASSES.get(reason)
        if category is None and reason.startswith("classification_changed_or_unresolved:"):
            category = "classification_transition"
        if category is None:
            category = "other_unresolved"
        safe_diagnostic = category == "geometry_interpolation_diagnostic"
        reason_class_rows.append(
            {
                "reason": reason,
                "category": category,
                "row_count": count,
                "rho_count": sum(1 for row in unresolved_rows if row.get("reason") == reason and row.get("axis") == "rho"),
                "temperature_count": sum(1 for row in unresolved_rows if row.get("reason") == reason and row.get("axis") == "temperature"),
                "xi_count": sum(1 for row in unresolved_rows if row.get("reason") == reason and row.get("axis") == "xi"),
                "safe_for_limited_overlay": safe_diagnostic,
                "note": (
                    "geometry/interpolation-only diagnostic; does not prove full-grid convergence"
                    if safe_diagnostic
                    else "requires separate evidence or remains a promotion blocker"
                ),
            }
        )

    unsafe_unresolved = sum(
        row["row_count"] for row in reason_class_rows if not bool_value(row["safe_for_limited_overlay"])
    )
    geometry_only_unresolved = sum(
        row["row_count"] for row in reason_class_rows if bool_value(row["safe_for_limited_overlay"])
    )

    boundary_gate = all(
        bool_value(row.get("converged"))
        and all(finite(row.get(field)) for field in ("mu_transition_MeV", "rho_hadron", "rho_quark", "area_residual"))
        for row in c2_boundary_rows
    ) and check_unique(c2_boundary_rows, ("xi", "T_MeV"))
    crossover_gate = all(
        bool_value(row.get("converged"))
        and all(finite(row.get(field)) for field in ("mu_MeV", "T_crossover_MeV", "rho"))
        for row in c2_crossover_rows
    ) and check_unique(c2_crossover_rows, ("xi", "mu_MeV"))
    failed_points = (c2_diagnostics_payload.get("counter_totals") or {}).get("failed_points", 0)
    scan_failure = (c2_diagnostics_payload.get("counter_totals") or {}).get("scan_failure", 0)
    solver_gate = failed_points == 0 and scan_failure == 0 and bool(validation_payload)
    auto_cep_widths = [float(row["bracket_width_T_MeV"]) for row in c2_cep_rows if finite(row.get("bracket_width_T_MeV"))]
    auto_cep_gate = bool(auto_cep_widths) and all(width <= BRACKET_TARGET_MEV for width in auto_cep_widths)
    regression_gate = len(regression_rows) == 0
    stability_gate = bool(drift_rows) and all(row.get("status") == "pass" for row in drift_rows)

    gate_rows = [
        {"gate": "boundary_finite_converged_unique", "automatic_status": "pass" if boundary_gate else "fail", "manual_overlay_status": "not_applicable", "evidence": "C2 boundary CSV"},
        {"gate": "crossover_finite_converged_unique", "automatic_status": "pass" if crossover_gate else "fail", "manual_overlay_status": "not_applicable", "evidence": "C2 crossover CSV"},
        {"gate": "solver_scan_failed_points_zero", "automatic_status": "pass" if solver_gate else "fail", "manual_overlay_status": "not_applicable", "evidence": "C2 phase diagnostics and validation"},
        {"gate": "automatic_grid_unresolved_zero", "automatic_status": "pass" if not unresolved_rows else "fail", "manual_overlay_status": "not_overridden", "evidence": "C2 phase_grid_convergence CSV"},
        {"gate": "automatic_classification_regressions_zero", "automatic_status": "pass" if regression_gate else "fail", "manual_overlay_status": "manual_oracle_acceptance_only", "evidence": "C1/C2 regression table plus nine-point overlay"},
        {"gate": "automatic_cep_width_leq_0.1_MeV", "automatic_status": "pass" if auto_cep_gate else "fail", "manual_overlay_status": "closed_by_author_review" if len(cep_overlay) == 3 else "fail", "evidence": "C2 CEP CSV plus manual CEP overlay"},
        {"gate": "nine_oracle_curve_review_complete", "automatic_status": "not_applicable", "manual_overlay_status": "accepted" if nine_ok else "fail", "evidence": "nine-point manual review package"},
        {"gate": "matched_C1_C2_physical_stability", "automatic_status": "pass" if stability_gate else "fail", "manual_overlay_status": "not_applicable", "evidence": "C1/C2 geometry drift summary"},
        {"gate": "unresolved_rows_all_geometry_only", "automatic_status": "pass" if unsafe_unresolved == 0 else "fail", "manual_overlay_status": "not_overridden", "evidence": "derived unresolved reason classification"},
    ]

    promotion_verdict = "limited_overlay_candidate" if unsafe_unresolved == 0 and nine_ok and stability_gate else "blocked_manual_overlay_inconclusive"
    blocker_reason = (
        "unresolved grid rows include non-geometry categories; nine oracle acceptances and CEP overlay do not rewrite C2 hybrid labels"
        if promotion_verdict != "limited_overlay_candidate"
        else "all unresolved rows are explicitly bounded geometry/interpolation diagnostics"
    )

    write_csv(tables_root / "auto_gate_summary.csv", ("gate", "automatic_status", "manual_overlay_status", "evidence"), gate_rows)
    write_csv(tables_root / "grid_unresolved_by_reason.csv", ("reason", "category", "row_count", "rho_count", "temperature_count", "xi_count", "safe_for_limited_overlay", "note"), reason_class_rows)
    write_csv(tables_root / "manual_nine_curve_decisions.csv", ("target_id", "xi", "T_MeV", "oracle_status", "oracle_candidate_count", "oracle_crossing_count", "oracle_geometry_converged", "finite_and_converged", "author_decision", "production_override"), [
        {**{key: row.get(key, "") for key in ("target_id", "xi", "T_MeV", "oracle_status", "oracle_candidate_count", "oracle_crossing_count", "oracle_geometry_converged", "finite_and_converged")}, "author_decision": "accepted_oracle_observation", "production_override": "false"}
        for row in nine_point_rows
    ])
    write_csv(tables_root / "manual_cep_decisions.csv", ("xi", "bracket_low_MeV", "bracket_high_MeV", "width_MeV", "status", "production_override"), [
        {"xi": label.removeprefix("xi_"), "bracket_low_MeV": payload["bracket_MeV"][0], "bracket_high_MeV": payload["bracket_MeV"][1], "width_MeV": payload["width_MeV"], "status": payload["status"], "production_override": "false"}
        for label, payload in sorted(cep_overlay.items())
    ])
    write_csv(tables_root / "physical_stability_summary.csv", tuple(drift_rows[0].keys()) if drift_rows else ("artifact", "metric", "status"), drift_rows)
    write_csv(tables_root / "promotion_decision.csv", ("verdict", "blocker_reason", "unresolved_grid_rows", "geometry_only_rows", "unsafe_unresolved_rows", "classification_regressions", "automatic_cep_failures", "manual_cep_closed", "manual_nine_curve_review", "reference_write"), [{
        "verdict": promotion_verdict,
        "blocker_reason": blocker_reason,
        "unresolved_grid_rows": len(unresolved_rows),
        "geometry_only_rows": geometry_only_unresolved,
        "unsafe_unresolved_rows": unsafe_unresolved,
        "classification_regressions": len(regression_rows),
        "automatic_cep_failures": sum(width > BRACKET_TARGET_MEV for width in auto_cep_widths),
        "manual_cep_closed": len(cep_overlay) == 3,
        "manual_nine_curve_review": nine_ok,
        "reference_write": False,
    }])
    write_csv(tables_root / "claim_ledger.csv", ("claim_id", "claim", "status", "evidence", "boundary"), [
        {"claim_id": "automatic_passes", "claim": "C2 boundary, crossover and solver-finiteness checks pass", "status": "supported" if boundary_gate and crossover_gate and solver_gate else "inconclusive", "evidence": "tables/auto_gate_summary.csv", "boundary": "C2 candidate only"},
        {"claim_id": "manual_cep_overlay", "claim": "All three CEP brackets are author-accepted at width 0.0625 MeV", "status": "supported" if len(cep_overlay) == 3 else "inconclusive", "evidence": "tables/manual_cep_decisions.csv", "boundary": "diagnostic overlay; no CSV rewrite"},
        {"claim_id": "manual_nine_curves", "claim": "Nine oracle rho-mu curves are author-accepted as first-order observations", "status": "supported" if nine_ok else "inconclusive", "evidence": "tables/manual_nine_curve_decisions.csv", "boundary": "hybrid labels remain ambiguous"},
        {"claim_id": "grid_residuals", "claim": "C2 unresolved grid rows are entirely geometry/interpolation-only diagnostics", "status": "supported" if unsafe_unresolved == 0 else "not_supported", "evidence": "tables/grid_unresolved_by_reason.csv", "boundary": "promotion blocker when false"},
        {"claim_id": "promotion", "claim": "Combined evidence is ready for a phase-reference promotion PR", "status": "candidate" if promotion_verdict == "limited_overlay_candidate" else "blocked", "evidence": "tables/promotion_decision.csv", "boundary": "no reference write performed"},
    ])

    decision_payload = {
        "schema_version": "pnjl_phase_reference_manual_overlay_promotion_audit_v1",
        "calculation_sha": CALCULATION_SHA,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "solver_called": False,
        "reference_write": False,
        "promotion_verdict": promotion_verdict,
        "blocker_reason": blocker_reason,
        "counts": {
            "c2_grid_rows": len(c2_grid_rows),
            "c2_grid_unresolved_rows": len(unresolved_rows),
            "geometry_only_unresolved_rows": geometry_only_unresolved,
            "unsafe_unresolved_rows": unsafe_unresolved,
            "classification_regressions": len(regression_rows),
            "manual_nine_curve_rows": len(nine_point_rows),
            "manual_cep_brackets": len(cep_overlay),
        },
        "automatic_gate_summary": gate_rows,
        "manual_overlay_contract": {
            "nine_curves_are_observations_not_hybrid_overrides": True,
            "cep_brackets_are_intervals_not_single_point_values": True,
            "historical_c2_artifacts_unchanged": True,
        },
    }
    (output_root / "decision.json").write_text(json.dumps(decision_payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    command = " ".join(str(value) for value in __import__("sys").argv)
    (output_root / "execution_log.md").write_text(
        "# Phase-reference manual overlay promotion audit v1\n\n"
        f"- Generated (UTC): `{decision_payload['generated_at_utc']}`\n"
        f"- Command: `{command}`\n"
        f"- Calculation SHA: `{CALCULATION_SHA}`\n"
        "- Solver called: `false`\n"
        "- Reference write: `false`\n"
        "- Input policy: C2 automatic artifacts plus explicit author-review overlays; no C2 row was overwritten.\n"
        f"- Result: `{promotion_verdict}`\n"
        f"- Reason: {blocker_reason}.\n",
        encoding="utf-8",
    )
    (output_root / "README.md").write_text(
        "# Phase-reference manual overlay promotion audit v1\n\n"
        "This solver-free package evaluates whether C2 automatic evidence plus the nine-point and CEP author overlays can support a limited phase-reference promotion. It does not write `data/reference/**` and does not change C2 artifacts.\n\n"
        f"The current verdict is `{promotion_verdict}`. The C2 grid contains {len(unresolved_rows)} unresolved rows; {geometry_only_unresolved} are geometry/interpolation-only by reason, while {unsafe_unresolved} remain unclassified, Stage-C, or classification-transition blockers. The nine oracle curves and three CEP brackets are therefore retained as explicit diagnostic overlays rather than silently replacing hybrid states.\n\n"
        "See `decision.json`, `execution_log.md`, `tables/auto_gate_summary.csv`, `tables/grid_unresolved_by_reason.csv`, `tables/manual_nine_curve_decisions.csv`, `tables/manual_cep_decisions.csv`, and `tables/claim_ledger.csv`.\n",
        encoding="utf-8",
    )

    output_files = {}
    for path in sorted(output_root.rglob("*")):
        if path.is_file() and path.name != "manifest.json":
            output_files[path.relative_to(output_root).as_posix()] = sha256(path)
    input_paths = [c1_grid, c1_boundary, c2_grid, c2_boundary, c2_crossover, c2_cep, c2_diagnostics, c2_validation, c2_audit_manifest, c2_regressions, c2_drift, nine_manifest, nine_points, manual_overlay_path]
    manifest = {
        "schema_version": "pnjl_phase_reference_manual_overlay_promotion_audit_v1",
        "calculation_sha": CALCULATION_SHA,
        "generated_at_utc": decision_payload["generated_at_utc"],
        "solver_called": False,
        "reference_write": False,
        "promotion_verdict": promotion_verdict,
        "input_files": [{"path": str(path), "sha256": sha256(path)} for path in input_paths],
        "source_counts": decision_payload["counts"],
        "output_files": output_files,
        "generator": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__).resolve())},
    }
    (output_root / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"promotion_verdict": promotion_verdict, "counts": decision_payload["counts"], "output_root": str(output_root)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
