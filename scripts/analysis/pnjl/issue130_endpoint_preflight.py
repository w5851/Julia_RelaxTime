#!/usr/bin/env python3
"""Create solver-free Issue #130 endpoint target lists and cost envelopes."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCHEMA = "pnjl_issue130_endpoint_refinement_preflight_v1"
CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
SOURCE_RUN = "31862752226"
RHO_STEP = 0.00625
RHO_MAX = 4.0
RHO_KEYS = int(round(RHO_MAX / RHO_STEP)) + 1
REPRESENTATIVE_XI = (-0.5, 0.0, 0.2875, 0.5)
CEP_FRACTIONS = (0.5, 0.75)
MAXWELL_WINDOW = 8.0
MAXWELL_PILOT_PER_XI = 3


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"empty CSV: {path}")
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def number(row: dict[str, str], field: str, context: str) -> float:
    try:
        value = float(row[field])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"invalid {field} in {context}: {row.get(field)!r}") from exc
    if not math.isfinite(value):
        raise ValueError(f"non-finite {field} in {context}")
    return value


def optional_number(row: dict[str, str], field: str) -> float | str:
    raw = row.get(field, "")
    if raw is None or not str(raw).strip():
        return ""
    value = float(raw)
    if not math.isfinite(value):
        raise ValueError(f"non-finite optional {field}")
    return value


def truth(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def xi_key(value: float) -> float:
    return round(value, 10)


def key(xi: float, temperature: float) -> tuple[float, float]:
    return xi_key(xi), round(temperature, 10)


def slug(value: float) -> str:
    return f"{value:g}".replace("-", "m").replace(".", "p")


def v5_inputs(root: Path) -> tuple[dict[str, Any], dict[str, Path]]:
    manifest_path = root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("calculation_sha") != CALCULATION_SHA:
        raise ValueError("unexpected v5 calculation SHA")
    if manifest.get("solver_called") is not False or manifest.get("reference_write") is not False:
        raise ValueError("v5 input is not solver-free/read-only")
    table_root = root / "tables"
    paths = {
        "gap": table_root / "crossover_cep_sampling_gap.csv",
        "maxwell": table_root / "maxwell_surface_point_status.csv",
        "grid": table_root / "grid_unresolved_diagnostics.csv",
        "separation": table_root / "crossover_maxwell_endpoint_separation.csv",
    }
    if any(not path.is_file() for path in paths.values()):
        raise ValueError("v5 input table missing")
    return {
        "manifest": manifest,
        "gap": read_csv(paths["gap"]),
        "maxwell": read_csv(paths["maxwell"]),
        "grid": read_csv(paths["grid"]),
        "separation": read_csv(paths["separation"]),
    }, paths


def build_crossover(rows: list[dict[str, str]]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    targets: list[dict[str, Any]] = []
    seen: set[float] = set()
    reps = {xi_key(value) for value in REPRESENTATIVE_XI}
    invalid = 0
    for line, row in enumerate(rows, start=2):
        context = f"gap:{line}"
        xi = number(row, "xi", context)
        xkey = xi_key(xi)
        if xkey in seen:
            raise ValueError(f"duplicate crossover xi {xi}")
        seen.add(xkey)
        cep_mu = number(row, "muq_CEP_proxy_MeV", context)
        last_mu = number(row, "last_retained_muq_MeV", context)
        low = number(row, "T_CEP_bracket_low_MeV", context)
        high = number(row, "T_CEP_bracket_high_MeV", context)
        last_T = number(row, "last_retained_T_MeV", context)
        gap = number(row, "native_mu_gap_MeV", context)
        if not last_mu < cep_mu:
            invalid += 1
            continue
        for index, fraction in enumerate(CEP_FRACTIONS, start=1):
            target_mu = last_mu + fraction * (cep_mu - last_mu)
            if not last_mu < target_mu < cep_mu:
                raise ValueError("crossover target is not strictly below CEP proxy")
            targets.append({
                "target_id": f"crossover_xi_{slug(xi)}_fraction_{index}",
                "target_kind": "crossover_mu_endpoint",
                "xi": xi,
                "target_mu_MeV": target_mu,
                "mu_lower_bound_MeV": last_mu,
                "mu_upper_bound_MeV": cep_mu,
                "fraction_from_last_retained": fraction,
                "last_retained_T_MeV": last_T,
                "T_CEP_bracket_low_MeV": min(low, high),
                "T_CEP_bracket_high_MeV": max(low, high),
                "native_mu_gap_MeV": gap,
                "physical_side": "crossover_mu_lt_CEP_proxy",
                "representative_xi": xkey in reps,
                "pilot_selection": "pilot_candidate" if xkey in reps else "full_only",
                "selection_reason": "deterministic_internal_fraction",
            })
    if not targets:
        raise ValueError("no valid crossover target interval")
    return targets, {
        "source_rows": len(rows),
        "valid_xi_rows": len(seen) - invalid,
        "invalid_interval_rows": invalid,
        "target_count": len(targets),
        "representative_pilot_target_count": sum(row["pilot_selection"] == "pilot_candidate" for row in targets),
        "fractions": list(CEP_FRACTIONS),
    }


def reason_rank(status: str, diagnostics: list[dict[str, str]]) -> tuple[int, str, str]:
    text = "|".join([status, *[row.get("reason", "") for row in diagnostics]]).lower()
    rho = "rho" in text or "geometry" in text
    temperature = "temperature" in text
    xi = "xi" in text
    if rho and (temperature or xi):
        return 1, "mixed_geometry_interpolation", "rho plus interpolation status"
    if rho:
        return 1, "rho_geometry", "rho geometry status"
    if temperature:
        return 2, "temperature_interpolation", "temperature interpolation status"
    if xi:
        return 3, "xi_interpolation", "xi interpolation status"
    return 4, "input_incomplete", "no actionable unresolved component"


def build_maxwell(
    gap_rows: list[dict[str, str]],
    maxwell_rows: list[dict[str, str]],
    grid_rows: list[dict[str, str]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    reps = {xi_key(value) for value in REPRESENTATIVE_XI}
    cep: dict[float, dict[str, float]] = {}
    for line, row in enumerate(gap_rows, start=2):
        context = f"gap:{line}"
        xi = number(row, "xi", context)
        cep[xi_key(xi)] = {
            "low": min(number(row, "T_CEP_bracket_low_MeV", context), number(row, "T_CEP_bracket_high_MeV", context)),
            "high": max(number(row, "T_CEP_bracket_low_MeV", context), number(row, "T_CEP_bracket_high_MeV", context)),
        }
    grid_by_key: defaultdict[tuple[float, float], list[dict[str, str]]] = defaultdict(list)
    unkeyable = 0
    for line, row in enumerate(grid_rows, start=2):
        context = f"grid:{line}"
        xi = number(row, "xi", context)
        raw_T = row.get("T_MeV", "")
        if raw_T is None or not str(raw_T).strip():
            unkeyable += 1
            continue
        grid_by_key[key(xi, number(row, "T_MeV", context))].append(row)
    maxwell_by_key: dict[tuple[float, float], dict[str, str]] = {}
    for line, row in enumerate(maxwell_rows, start=2):
        context = f"maxwell:{line}"
        item_key = key(number(row, "xi", context), number(row, "T_MeV", context))
        if item_key in maxwell_by_key:
            raise ValueError(f"duplicate Maxwell key {item_key}")
        maxwell_by_key[item_key] = row

    targets: list[dict[str, Any]] = []
    outside = 0
    non_unresolved = 0
    for item_key, row in maxwell_by_key.items():
        if not truth(row.get("grid_unresolved", "false")):
            non_unresolved += 1
            continue
        xi, temperature = item_key
        bracket = cep.get(xi)
        if bracket is None:
            continue
        distance = bracket["low"] - temperature
        if temperature > bracket["low"] or distance > MAXWELL_WINDOW:
            outside += 1
            continue
        diagnostics = grid_by_key.get(item_key, [])
        priority, reason, detail = reason_rank(row.get("grid_status", ""), diagnostics)
        targets.append({
            "target_id": f"maxwell_xi_{slug(xi)}_T_{slug(temperature)}",
            "target_kind": "maxwell_fixed_xi_T",
            "xi": xi,
            "T_MeV": temperature,
            "target_mu_MeV": optional_number(row, "mu_MeV"),
            "cep_T_low_MeV": bracket["low"],
            "cep_T_high_MeV": bracket["high"],
            "distance_below_CEP_MeV": distance,
            "grid_status": row.get("grid_status", ""),
            "grid_unresolved": True,
            "boundary_row_present": True,
            "priority": priority,
            "reason": reason,
            "reason_detail": detail,
            "diagnostic_count": len(diagnostics),
            "representative_xi": xi in reps,
            "pilot_selection": "full_candidate",
        })

    # A grid diagnostic without a boundary row is retained as an audit target,
    # but it cannot enter numerical pilot because its source payload is incomplete.
    for item_key, diagnostics in grid_by_key.items():
        if item_key in maxwell_by_key or item_key[0] not in cep:
            continue
        bracket = cep[item_key[0]]
        distance = bracket["low"] - item_key[1]
        unresolved = [row for row in diagnostics if row.get("reason", "") != "ok" or not truth(row.get("boundary_row_present", "false"))]
        if item_key[1] > bracket["low"] or distance > MAXWELL_WINDOW or not unresolved:
            continue
        priority, _reason, detail = reason_rank("|".join(row.get("reason", "") for row in unresolved), unresolved)
        targets.append({
            "target_id": f"maxwell_xi_{slug(item_key[0])}_T_{slug(item_key[1])}",
            "target_kind": "maxwell_fixed_xi_T",
            "xi": item_key[0],
            "T_MeV": item_key[1],
            "target_mu_MeV": "",
            "cep_T_low_MeV": bracket["low"],
            "cep_T_high_MeV": bracket["high"],
            "distance_below_CEP_MeV": distance,
            "grid_status": "|".join(sorted({row.get("reason", "") for row in unresolved})),
            "grid_unresolved": True,
            "boundary_row_present": False,
            "priority": priority,
            "reason": "input_incomplete",
            "reason_detail": f"missing boundary row; {detail}",
            "diagnostic_count": len(diagnostics),
            "representative_xi": item_key[0] in reps,
            "pilot_selection": "input_incomplete",
        })

    grouped: defaultdict[float, list[dict[str, Any]]] = defaultdict(list)
    for row in targets:
        grouped[xi_key(float(row["xi"]))].append(row)
    for rows in grouped.values():
        rows.sort(key=lambda row: (int(row["priority"]), float(row["distance_below_CEP_MeV"]), float(row["T_MeV"]), row["target_id"]))
        eligible = 0
        for rank, row in enumerate(rows, start=1):
            row["rank_within_xi"] = rank
            if row["representative_xi"] and row["boundary_row_present"] and eligible < MAXWELL_PILOT_PER_XI:
                row["pilot_selection"] = "pilot_candidate"
                eligible += 1
            elif row["boundary_row_present"]:
                row["pilot_selection"] = "full_only"
    targets.sort(key=lambda row: (float(row["xi"]), int(row["rank_within_xi"])))
    reasons: defaultdict[str, int] = defaultdict(int)
    for row in targets:
        reasons[str(row["reason"])] += 1
    return targets, {
        "source_maxwell_rows": len(maxwell_rows),
        "source_grid_rows": len(grid_rows),
        "unkeyable_grid_rows": unkeyable,
        "candidate_count": len(targets),
        "eligible_existing_boundary_targets": sum(row["boundary_row_present"] for row in targets),
        "input_incomplete_targets": sum(not row["boundary_row_present"] for row in targets),
        "representative_pilot_target_count": sum(row["pilot_selection"] == "pilot_candidate" for row in targets),
        "outside_window_unresolved_rows": outside,
        "non_unresolved_maxwell_rows": non_unresolved,
        "cep_window_below_bracket_low_MeV": MAXWELL_WINDOW,
        "reason_counts": dict(sorted(reasons.items())),
    }


def cost(route: str, pilot: int, eligible: int) -> list[dict[str, Any]]:
    return [
        {
            "route": route,
            "scope": scope,
            "target_count": count,
            "rho_base_step": RHO_STEP,
            "rho_max": RHO_MAX,
            "base_rho_keys_per_target": RHO_KEYS,
            "expected_base_rho_keys": count * RHO_KEYS,
            "targeted_refinement_points": "not_budgeted_in_preflight",
            "expected_unique_solves": count * RHO_KEYS,
            "wall_time_status": "not_measured_solver_free_preflight",
            "cost_unit": "base_rho_curve_point",
            "stop_condition": "no_numerical_dispatch_from_preflight",
        }
        for scope, count in (("representative_xi_pilot", pilot), ("all_eligible_candidates", eligible))
    ]


def files_hashes(root: Path) -> dict[str, str]:
    return {
        path.relative_to(root).as_posix(): digest(path)
        for path in sorted(root.rglob("*"))
        if path.is_file() and path.relative_to(root).as_posix() != "manifest.json"
    }


def route_package(root: Path, name: str, targets: list[dict[str, Any]], summary: dict[str, Any], costs: list[dict[str, Any]], source_hashes: dict[str, str]) -> dict[str, Any]:
    route_root = root / name
    table_root = route_root / "tables"
    write_csv(table_root / "target_list.csv", targets)
    write_csv(table_root / "cost_envelope.csv", costs)
    write_json(table_root / "preflight_summary.json", summary)
    claims = [
        {"claim_id": f"{name}_solver_free", "claim": "target selection is solver-free", "status": "supported", "evidence": "manifest.json; tables/target_list.csv", "boundary": "not a physical certificate"},
        {"claim_id": f"{name}_pilot_pending", "claim": "numerical pilot is not authorized", "status": "supported", "evidence": "README.md; tables/cost_envelope.csv", "boundary": "requires separate author authorization"},
    ]
    write_csv(table_root / "claim_ledger.csv", claims)
    (route_root / "README.md").write_text(
        f"# Issue #130 {name} preflight\n\n"
        "Solver-free target selection from immutable C2 v5 evidence. No Julia, "
        "equilibrium, Maxwell, or oracle call is made. This is a candidate list, "
        "not a physical certificate.\n\n"
        f"- candidate count: `{summary.get('target_count', summary.get('candidate_count'))}`\n"
        f"- representative pilot count: `{summary['representative_pilot_target_count']}`\n"
        f"- base rho grid: `{RHO_STEP}` to `{RHO_MAX}` ({RHO_KEYS} points)\n"
        "- wall time: not measured in solver-free preflight\n",
        encoding="utf-8",
    )
    manifest = {
        "schema_version": SCHEMA,
        "route": name,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "calculation_sha": CALCULATION_SHA,
        "source_run_id": SOURCE_RUN,
        "solver_called": False,
        "reference_write": False,
        "run_mode": "solver_free_preflight",
        "verdict": f"{name}_target_list_ready_cost_unmeasured",
        "policy": {"representative_xi": list(REPRESENTATIVE_XI), "rho_base_step": RHO_STEP, "rho_max": RHO_MAX, "base_rho_keys_per_target": RHO_KEYS},
        "source_hashes": source_hashes,
        "summary": summary,
        "files": {},
    }
    manifest["files"] = files_hashes(route_root)
    write_json(route_root / "manifest.json", manifest)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v5-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    data, paths = v5_inputs(args.v5_root.resolve())
    source_hashes = {"v5_manifest.json": digest(args.v5_root / "manifest.json")}
    source_hashes.update({f"tables/{path.name}": digest(path) for path in paths.values()})
    crossover, crossover_summary = build_crossover(data["gap"])
    maxwell, maxwell_summary = build_maxwell(data["gap"], data["maxwell"], data["grid"])
    crossover_pilot = sum(row["pilot_selection"] == "pilot_candidate" for row in crossover)
    maxwell_pilot = sum(row["pilot_selection"] == "pilot_candidate" for row in maxwell)
    root = args.output_root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    crossover_manifest = route_package(root, "crossover_mu", crossover, crossover_summary, cost("crossover_mu_endpoint", crossover_pilot, len(crossover)), source_hashes)
    maxwell_manifest = route_package(root, "maxwell_local", maxwell, maxwell_summary, cost("maxwell_cep_local", maxwell_pilot, sum(row["boundary_row_present"] for row in maxwell)), source_hashes)
    (root / "README.md").write_text(
        "# Issue #130 endpoint refinement preflight v1\n\n"
        "Two solver-free routes are recorded: crossover internal mu targets below "
        "the CEP proxy, and fixed (xi,T) Maxwell unresolved targets near CEP. "
        "Neither route is a certificate; numerical pilot authorization remains pending.\n\n"
        f"- calculation SHA: `{CALCULATION_SHA}`\n- source run: `{SOURCE_RUN}`\n"
        f"- crossover candidates/pilot: `{len(crossover)}` / `{crossover_pilot}`\n"
        f"- Maxwell candidates/eligible/pilot: `{len(maxwell)}` / `{sum(row['boundary_row_present'] for row in maxwell)}` / `{maxwell_pilot}`\n"
        f"- base rho grid: `{RHO_STEP}` to `{RHO_MAX}` ({RHO_KEYS} points)\n",
        encoding="utf-8",
    )
    manifest = {
        "schema_version": SCHEMA,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "calculation_sha": CALCULATION_SHA,
        "source_run_id": SOURCE_RUN,
        "workflow_head_sha": None,
        "repo_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=Path(__file__).resolve().parents[3], text=True).strip(),
        "solver_called": False,
        "reference_write": False,
        "run_mode": "solver_free_preflight",
        "verdict": "two_route_target_lists_ready_cost_unmeasured",
        "pilot_authorization": "pending",
        "routes": {"crossover_mu": crossover_manifest["summary"], "maxwell_local": maxwell_manifest["summary"]},
        "source_hashes": source_hashes,
        "files": {},
    }
    manifest["files"] = files_hashes(root)
    write_json(root / "manifest.json", manifest)
    print(json.dumps({"crossover": crossover_summary, "maxwell": maxwell_summary}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
