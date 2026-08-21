#!/usr/bin/env python3
"""Run a solver-free promotion gate for the Issue #130 phase layers.

The gate validates the immutable strict/derived/render package after author
review.  It does not copy files into ``data/reference`` and it never calls the
PNJL solver.  A passing result means that a separate, versioned promotion PR
may be prepared; it does not authorize runtime consumption by itself.
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
from typing import Any, Iterable


CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
SOURCE_RUN_ID = "32354095831"
REPLAY_RUN_ID = "32451053476"
EXPECTED_XI_COUNT = 161
EXPECTED_XI_STEP = 0.00625
EXPECTED_MAXWELL_ROWS = 7162
EXPECTED_EXPANSION_ROWS = 276
EXPECTED_RENDER_COVERAGE_ROWS = 483
EXPECTED_REPLAY_TARGETS = 276
EPS = 1.0e-9


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--author-review-status", choices=("pending", "accepted"), default="pending")
    parser.add_argument("--review-pr", default="")
    parser.add_argument("--review-merge-sha", default="")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
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


def finite(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def as_bool(value: Any) -> bool:
    return str(value).strip().lower() == "true" if not isinstance(value, bool) else value


def check_unique(rows: Iterable[dict[str, Any]], keys: tuple[str, ...]) -> bool:
    values = [tuple(row.get(key, "") for key in keys) for row in rows]
    return len(values) == len(set(values))


def require(path: Path) -> Path:
    if not path.is_file():
        raise FileNotFoundError(path)
    return path


def expected_xi_grid() -> list[float]:
    return [round(-0.5 + index * EXPECTED_XI_STEP, 8) for index in range(EXPECTED_XI_COUNT)]


def xi_values(rows: Iterable[dict[str, str]]) -> list[float]:
    return sorted({round(float(row["xi"]), 8) for row in rows})


def gate(name: str, passed: bool, evidence: str, details: str = "") -> dict[str, str]:
    return {
        "gate": name,
        "status": "pass" if passed else "fail",
        "evidence": evidence,
        "details": details,
    }


def _all_finite(rows: Iterable[dict[str, str]], fields: tuple[str, ...]) -> bool:
    return all(finite(row.get(field)) for row in rows for field in fields)


def validate_strict(package_root: Path) -> tuple[list[dict[str, str]], dict[str, Any]]:
    root = package_root / "strict_reference_v1" / "tables"
    maxwell = read_csv(require(root / "maxwell_surface_strict_reference_v1.csv"))
    crossover = read_csv(require(root / "crossover_surface_strict_reference_v1.csv"))
    spinodal = read_csv(require(root / "spinodal_surface_strict_reference_v1.csv"))
    cep = read_csv(require(root / "cep_boundary_strict_reference_v1.csv"))

    maxwell_ok = (
        len(maxwell) == EXPECTED_MAXWELL_ROWS
        and check_unique(maxwell, ("xi", "T_MeV"))
        and _all_finite(maxwell, ("xi", "T_MeV", "mu_MeV", "rho_hadron", "rho_quark", "area_residual"))
        and all(row.get("surface") == "maxwell" for row in maxwell)
        and all(row.get("layer") == "strict_reference_v1" for row in maxwell)
        and all(row.get("reference_write") == "False" for row in maxwell)
        and all(row.get("oracle_labels_consumed") == "False" for row in maxwell)
    )
    expansion_rows = [row for row in maxwell if row.get("status") == "expansion_native_geometry_converged"]
    expansion_ok = (
        len(expansion_rows) == EXPECTED_EXPANSION_ROWS
        and all(row.get("source_run_id") == SOURCE_RUN_ID for row in expansion_rows)
        and all(row.get("source_target_id") for row in expansion_rows)
        and all(as_bool(row.get("geometry_converged")) for row in expansion_rows)
        and all(as_bool(row.get("finite_and_converged")) for row in expansion_rows)
    )
    unresolved = [row for row in maxwell if as_bool(row.get("grid_unresolved"))]
    unresolved_statuses = Counter(row.get("status", "") for row in unresolved)

    crossover_ok = (
        len(crossover) == 1343
        and check_unique(crossover, ("xi", "mu_MeV"))
        and _all_finite(crossover, ("xi", "mu_MeV", "T_MeV", "rho", "mu_CEP_proxy_MeV"))
        and all(row.get("layer") == "strict_reference_v1" for row in crossover)
        and all(row.get("reference_write") == "False" for row in crossover)
    )
    spinodal_ok = (
        len(spinodal) == 6886
        and check_unique(spinodal, ("xi", "T_MeV"))
        and _all_finite(spinodal, ("xi", "T_MeV", "mu_spinodal_hadron_MeV", "mu_spinodal_quark_MeV"))
        and all(row.get("layer") == "strict_reference_v1" for row in spinodal)
        and all(row.get("reference_write") == "False" for row in spinodal)
    )
    cep_ok = (
        len(cep) == 93
        and check_unique(cep, ("xi",))
        and _all_finite(cep, ("xi", "mu_CEP_proxy_MeV", "T_low_MeV", "T_high_MeV", "T_midpoint_MeV"))
        and all(float(row["T_low_MeV"]) <= float(row["T_midpoint_MeV"]) <= float(row["T_high_MeV"]) for row in cep)
        and all(row.get("layer") == "strict_reference_v1" for row in cep)
        and all(row.get("reference_write") == "False" for row in cep)
    )
    gates = [
        gate("strict_maxwell_finite_unique", maxwell_ok, "strict_reference_v1/tables/maxwell_surface_strict_reference_v1.csv"),
        gate("strict_maxwell_expansion_provenance", expansion_ok, "strict_reference_v1/tables/maxwell_surface_strict_reference_v1.csv", f"expansion_rows={len(expansion_rows)}"),
        gate("strict_unresolved_preserved", bool(unresolved) and maxwell_ok, "strict_reference_v1/tables/maxwell_surface_strict_reference_v1.csv", f"unresolved_rows={len(unresolved)}"),
        gate("strict_crossover_finite_unique", crossover_ok, "strict_reference_v1/tables/crossover_surface_strict_reference_v1.csv"),
        gate("strict_spinodal_finite_unique", spinodal_ok, "strict_reference_v1/tables/spinodal_surface_strict_reference_v1.csv"),
        gate("strict_cep_brackets_preserved", cep_ok, "strict_reference_v1/tables/cep_boundary_strict_reference_v1.csv"),
    ]
    summary = {
        "maxwell_rows": len(maxwell),
        "maxwell_expansion_rows": len(expansion_rows),
        "maxwell_unresolved_rows": len(unresolved),
        "maxwell_unresolved_statuses": dict(sorted(unresolved_statuses.items())),
        "crossover_rows": len(crossover),
        "spinodal_rows": len(spinodal),
        "cep_rows": len(cep),
    }
    return gates, summary


def validate_derived(package_root: Path) -> tuple[list[dict[str, str]], dict[str, Any]]:
    root = package_root / "derived_reference_v1" / "tables"
    tables = {
        "maxwell": read_csv(require(root / "maxwell_surface_derived_reference_v1.csv")),
        "crossover": read_csv(require(root / "crossover_surface_derived_reference_v1.csv")),
        "spinodal": read_csv(require(root / "spinodal_surface_derived_reference_v1.csv")),
        "cep": read_csv(require(root / "cep_boundary_derived_reference_v1.csv")),
        "coverage": read_csv(require(root / "surface_coverage_mask.csv")),
    }
    target = expected_xi_grid()
    surface_requirements = {
        "maxwell": (("xi", "T_MeV"), ("xi", "T_MeV", "mu_MeV", "rho_hadron", "rho_quark", "area_residual"), 12537),
        "crossover": (("xi", "mu_MeV"), ("xi", "mu_MeV", "T_MeV", "rho", "mu_CEP_proxy_MeV"), 3135),
        "spinodal": (("xi", "T_MeV"), ("xi", "T_MeV", "mu_spinodal_hadron_MeV", "mu_spinodal_quark_MeV"), 11989),
    }
    surface_ok = {}
    for name, (keys, fields, expected_rows) in surface_requirements.items():
        rows = tables[name]
        surface_ok[name] = (
            len(rows) == expected_rows
            and xi_values(rows) == target
            and check_unique(rows, keys)
            and _all_finite(rows, fields)
            and all(row.get("layer") in {"strict_reference_v1", "interpolated_noncertified"} for row in rows)
            and all(row.get("reference_write") == "False" for row in rows)
        )
    cep = tables["cep"]
    cep_ok = (
        len(cep) == EXPECTED_XI_COUNT
        and xi_values(cep) == target
        and check_unique(cep, ("xi",))
        and _all_finite(cep, ("xi", "mu_CEP_proxy_MeV", "T_low_MeV", "T_high_MeV", "T_midpoint_MeV"))
        and all(float(row["T_low_MeV"]) <= float(row["T_midpoint_MeV"]) <= float(row["T_high_MeV"]) for row in cep)
        and all(row.get("layer") in {"strict_reference_v1", "interpolated_noncertified"} for row in cep)
        and all(row.get("reference_write") == "False" for row in cep)
    )
    coverage = tables["coverage"]
    coverage_ok = (
        len(coverage) == EXPECTED_RENDER_COVERAGE_ROWS
        and check_unique(coverage, ("surface", "xi"))
        and set(row.get("surface") for row in coverage) == {"crossover", "maxwell", "spinodal"}
        and all(xi_values([row]) == [xi] for row in coverage for xi in [round(float(row["xi"]), 8)])
        and all(row.get("coverage_status") in {"native_support", "interpolated_common_support"} for row in coverage)
    )
    no_extrapolation = all(
        row.get("coverage_status") in {"native_support", "interpolated_common_support"}
        and row.get("axis_low", "") != ""
        and row.get("axis_high", "") != ""
        for row in coverage
    )
    gates = [
        gate("derived_uniform_xi_grid", all(surface_ok.values()) and cep_ok, "derived_reference_v1/tables/*.csv", f"xi_count={EXPECTED_XI_COUNT}"),
        gate("derived_common_support_only", coverage_ok and no_extrapolation, "derived_reference_v1/tables/surface_coverage_mask.csv"),
        gate("derived_noncertified_rows_explicit", all(row.get("layer") != "interpolated_noncertified" or row.get("status") != "native_uniform_xi" for rows in tables.values() for row in rows), "derived_reference_v1/tables/*.csv"),
    ]
    return gates, {"surface_rows": {name: len(rows) for name, rows in tables.items()}, "coverage_rows": len(coverage)}


def validate_render(package_root: Path) -> tuple[list[dict[str, str]], dict[str, Any]]:
    root = package_root / "phase_surface_render_v1"
    manifest_path = require(root / "manifest.json")
    manifest = read_json(manifest_path)
    plot_manifest = read_json(require(root / "figures" / "plot_manifest.json"))
    figure = require(root / "figures" / "phase_surface_render_mu_xi_T.png")
    input_records = plot_manifest.get("inputs") or {}
    semantics = plot_manifest.get("render_semantics") or {}
    project_root = package_root.parents[4]
    input_checks = []
    for key in ("derived_manifest", "crossover_table", "maxwell_table", "cep_table", "coverage_table"):
        record = input_records.get(key) or {}
        path = Path(record.get("path", ""))
        if not path.is_absolute():
            path = project_root / path
        if path.is_file():
            input_checks.append(sha256(path) == record.get("sha256"))
        else:
            input_checks.append(False)
    render_manifest_hash = manifest.get("figures", {}).get("plot_manifest.json")
    plot_manifest_hash_ok = render_manifest_hash == sha256(root / "figures" / "plot_manifest.json")
    figure_ok = (
        plot_manifest.get("figure", {}).get("sha256") == sha256(figure)
        and bool(input_checks)
        and all(input_checks)
        and plot_manifest_hash_ok
        and semantics.get("triangulation") is False
        and semantics.get("masked_cells_connected") is False
        and semantics.get("derived_rows_are_certified") is False
    )
    gates = [
        gate("render_no_triangulation", figure_ok, "phase_surface_render_v1/manifest.json; figures/plot_manifest.json"),
        gate("render_inputs_hash_bound", all(input_checks) and plot_manifest_hash_ok, "phase_surface_render_v1/figures/plot_manifest.json"),
    ]
    return gates, {"figure_sha256": sha256(figure), "figure_bytes": figure.stat().st_size, "render_input_hashes": input_checks}


def validate_common(package_root: Path) -> tuple[list[dict[str, str]], dict[str, Any]]:
    root_manifest = read_json(require(package_root / "manifest.json"))
    decision = read_json(require(package_root / "decision.json"))
    replay = read_json(require(package_root / "replay_provenance.json"))
    inventory_ok = True
    for relative, expected in (root_manifest.get("package_files") or {}).items():
        path = package_root / relative
        inventory_ok = inventory_ok and path.is_file() and sha256(path) == expected
    for child, expected in (root_manifest.get("children") or {}).items():
        path = package_root / child / "manifest.json"
        inventory_ok = inventory_ok and path.is_file() and sha256(path) == expected
    render_manifest = package_root / "phase_surface_render_v1" / "manifest.json"
    plot_manifest = package_root / "phase_surface_render_v1" / "figures" / "plot_manifest.json"
    inventory_ok = inventory_ok and render_manifest.is_file() and plot_manifest.is_file()
    common_ok = (
        root_manifest.get("calculation_sha") == CALCULATION_SHA
        and root_manifest.get("maxwell_source_run_id") == SOURCE_RUN_ID
        and root_manifest.get("reference_write") is False
        and root_manifest.get("solver_called") is False
        and root_manifest.get("oracle_labels_consumed") is False
        and decision.get("reference_write") is False
        and replay.get("run_mode") == "aggregate_replay"
        and replay.get("source_run_id") == SOURCE_RUN_ID
        and replay.get("replay_run_id") == REPLAY_RUN_ID
        and replay.get("calculation_sha") == CALCULATION_SHA
        and replay.get("solver_called") is False
        and replay.get("materialized_target_count") == EXPECTED_REPLAY_TARGETS
        and replay.get("expected_target_count") == EXPECTED_REPLAY_TARGETS
    )
    gates = [
        gate("package_and_replay_provenance", common_ok, "manifest.json; decision.json; replay_provenance.json"),
        gate("package_manifest_hash_inventory", inventory_ok, "manifest.json; child manifests"),
    ]
    return gates, {"calculation_sha": root_manifest.get("calculation_sha"), "source_run_id": replay.get("source_run_id"), "replay_run_id": replay.get("replay_run_id")}


def write_csv(path: Path, fieldnames: tuple[str, ...], rows: Iterable[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def audit_package(package_root: Path, output_root: Path, *, author_review_status: str = "pending", review_pr: str = "", review_merge_sha: str = "") -> dict[str, Any]:
    package_root = package_root.resolve()
    output_root = output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    (output_root / "tables").mkdir(parents=True, exist_ok=True)
    gates: list[dict[str, str]] = []
    summary: dict[str, Any] = {}
    for validator in (validate_common, validate_strict, validate_derived, validate_render):
        current_gates, current_summary = validator(package_root)
        gates.extend(current_gates)
        summary.update(current_summary)
    review_ok = author_review_status == "accepted" and bool(review_pr) and len(review_merge_sha) == 40
    gates.append(gate("author_review_recorded", review_ok, "promotion gate invocation", f"status={author_review_status}; pr={review_pr}"))
    verdict = "promotion_candidate" if all(row["status"] == "pass" for row in gates) else "promotion_blocked"
    review_record = {
        "schema_version": "pnjl_issue130_phase_reference_author_review_v1",
        "status": author_review_status,
        "review_pr": review_pr,
        "review_merge_sha": review_merge_sha,
        "package_manifest_sha256": sha256(package_root / "manifest.json"),
        "scope": ["strict_reference_v1", "derived_reference_v1", "phase_surface_render_v1", "coverage", "provenance"],
        "reference_write_authorized": False,
        "runtime_consumption_authorized": False,
    }
    write_json(output_root / "author_review_record_v1.json", review_record)
    write_csv(output_root / "tables" / "gate_summary.csv", ("gate", "status", "evidence", "details"), gates)
    unresolved_rows = [{"status": key, "row_count": value} for key, value in sorted(summary.get("maxwell_unresolved_statuses", {}).items())]
    write_csv(output_root / "tables" / "unresolved_summary.csv", ("status", "row_count"), unresolved_rows)
    target = {
        "planned_reference_root": "data/reference/pnjl/issue130_phase_reference_v1",
        "planned_layers": {"strict": "strict", "derived": "derived", "render_data": "render"},
        "planned_figure_root": "data/outputs/figures/pnjl/phase_reference/issue130_phase_reference_v1",
        "write_in_this_gate": False,
        "runtime_consumption_in_this_gate": False,
        "legacy_files_untouched": ["boundary.csv", "cep.csv", "spinodals.csv", "crossover_dense.csv", "phase_reference_dense_manifest.json"],
    }
    write_json(output_root / "reference_target.json", target)
    claim_rows = [
        {"claim_id": "package_integrity", "claim": "The reviewed three-layer package is hash-bound and replay-complete", "status": "supported" if all(row["status"] == "pass" for row in gates if row["gate"] != "author_review_recorded") else "blocked", "evidence": "tables/gate_summary.csv; manifest.json", "boundary": "solver-free gate only"},
        {"claim_id": "unresolved_semantics", "claim": "Strict unresolved geometry rows remain explicit and are not certified by derived interpolation", "status": "supported" if summary.get("maxwell_unresolved_rows", 0) > 0 else "blocked", "evidence": "tables/unresolved_summary.csv; strict_reference_v1/tables/maxwell_surface_strict_reference_v1.csv", "boundary": "downstream policy must choose strict versus derived explicitly"},
        {"claim_id": "versioned_promotion", "claim": "A separate versioned promotion PR may be prepared", "status": "candidate" if verdict == "promotion_candidate" else "blocked", "evidence": "decision.json; reference_target.json", "boundary": "no data/reference write or runtime switch in this gate"},
    ]
    write_csv(output_root / "claim_ledger.csv", ("claim_id", "claim", "status", "evidence", "boundary"), claim_rows)
    decision = {
        "schema_version": "pnjl_issue130_phase_reference_promotion_gate_v1",
        "verdict": verdict,
        "promotion_status": "eligible_for_versioned_promotion" if verdict == "promotion_candidate" else "blocked",
        "reference_write": False,
        "runtime_consumption": False,
        "package_root": str(package_root),
        "calculation_sha": CALCULATION_SHA,
        "source_run_id": SOURCE_RUN_ID,
        "replay_run_id": REPLAY_RUN_ID,
        "summary": summary,
        "next_action": "prepare a separate versioned reference import PR" if verdict == "promotion_candidate" else "repair the failed gate without changing numerical evidence",
    }
    write_json(output_root / "decision.json", decision)
    audit_lines = [
        "# Issue #130 phase-reference promotion gate v1",
        "",
        f"- Verdict: `{verdict}`",
        f"- Calculation SHA: `{CALCULATION_SHA}`",
        f"- Numerical source run: `{SOURCE_RUN_ID}`; aggregate replay: `{REPLAY_RUN_ID}`",
        f"- Strict Maxwell unresolved rows retained: `{summary.get('maxwell_unresolved_rows', 0)}`",
        "- This gate is solver-free. It does not write `data/reference` and does not switch runtime consumers.",
        "- A passing gate authorizes preparation of a separate versioned import PR only.",
        "",
        "## Boundaries",
        "",
        "- `strict_reference_v1` remains the evidence/source layer; unresolved geometry status is preserved.",
        "- `derived_reference_v1` is an explicit downstream derived layer; interpolated rows are not certificates.",
        "- `phase_surface_render_v1` is a reproducible render/data projection and cannot change physical status.",
    ]
    (output_root / "AUDIT.md").write_text("\n".join(audit_lines) + "\n", encoding="utf-8")
    manifest_files = {}
    for path in sorted(output_root.rglob("*")):
        if path.is_file() and path.name != "manifest.json":
            manifest_files[path.relative_to(output_root).as_posix()] = sha256(path)
    manifest = {
        "schema_version": "pnjl_issue130_phase_reference_promotion_gate_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "verdict": verdict,
        "calculation_sha": CALCULATION_SHA,
        "source_run_id": SOURCE_RUN_ID,
        "replay_run_id": REPLAY_RUN_ID,
        "package_manifest_sha256": sha256(package_root / "manifest.json"),
        "reference_write": False,
        "runtime_consumption": False,
        "files": manifest_files,
    }
    write_json(output_root / "manifest.json", manifest)
    return decision


def main() -> int:
    args = parse_args()
    decision = audit_package(
        args.package_root,
        args.output_root,
        author_review_status=args.author_review_status,
        review_pr=args.review_pr,
        review_merge_sha=args.review_merge_sha,
    )
    print(json.dumps({"verdict": decision["verdict"], "output_root": str(args.output_root)}, ensure_ascii=False))
    return 0 if decision["verdict"] == "promotion_candidate" else 1


if __name__ == "__main__":
    raise SystemExit(main())
