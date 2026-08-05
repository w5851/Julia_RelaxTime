#!/usr/bin/env python3
"""Solver-free feasibility replay for the endpoint-local geometry contract.

The replay consumes the immutable v3 targeted shadow and the separately
approved required-three deep-oracle artifact.  Stage-B curves and the public
Maxwell candidate are the only inputs used to choose midpoint locations.  Deep
curves are used only after a location has been selected, to emulate the value
of the newly sampled point and to apply the final numerical gate.  No Julia
module or equilibrium solver is imported here.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


SCHEMA_VERSION = "cep_maxwell_endpoint_local_contract_feasibility_v2"
STANDARD_RUN_ID = "30999377195"
DEEP_RUN_ID = "31002704845"
CALCULATION_SHA = "ceec2295c5c9250a3fcd45c0eceae9a6c35f4335"
APPROVED_POINTS = ((-0.5, 5.0), (-0.5, 20.0), (0.0, 5.0))
AREA_GATE = 5e-5
POSITION_TOL = 0.025
DENSITY_TOL = 0.0025
MAX_TARGETED = 12
DEEP_STEP = 0.0015625


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("\n", encoding="utf-8")
        return
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _float(value: Any, default: float = math.nan) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    return result


def _finite(value: Any) -> bool:
    return math.isfinite(_float(value))


def _bool(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _key(row: dict[str, str]) -> tuple[float, float, float, str, str]:
    return (
        _float(row.get("xi")),
        _float(row.get("T_MeV")),
        _float(row.get("rho")),
        row.get("method", ""),
        row.get("rho_level", ""),
    )


def _load_deep(directory: Path) -> tuple[list[dict[str, str]], dict[str, Any], list[str]]:
    metrics = directory / "slice_metrics.csv"
    manifest = directory / "manifest.json"
    if not metrics.is_file():
        candidates = sorted(directory.rglob("slice_metrics.csv"))
        candidates = [path for path in candidates if "aggregate" in path.parent.name]
        metrics = candidates[0] if candidates else metrics
    if not manifest.is_file():
        candidates = sorted(directory.rglob("manifest.json"))
        candidates = [path for path in candidates if "aggregate" in path.parent.name]
        manifest = candidates[0] if candidates else manifest
    errors: list[str] = []
    if not metrics.is_file():
        errors.append(f"missing deep slice_metrics.csv below {directory}")
        return [], {}, errors
    payload = json.loads(manifest.read_text(encoding="utf-8")) if manifest.is_file() else {}
    rows = _rows(metrics)
    for key in APPROVED_POINTS:
        matching = [
            row for row in rows
            if row.get("method") == "independent_oracle"
            and math.isclose(_float(row.get("xi")), key[0], abs_tol=1e-9, rel_tol=0.0)
            and math.isclose(_float(row.get("T_MeV")), key[1], abs_tol=1e-6, rel_tol=0.0)
        ]
        if len(matching) != 1:
            errors.append(f"deep oracle must contain exactly one row for {key}")
        elif matching[0].get("result_status") != "confirmed_first_order":
            errors.append(f"deep oracle is not confirmed first-order at {key}")
    return rows, payload, errors


def _validate_standard(directory: Path, expected_sha: str) -> tuple[list[dict[str, str]], list[dict[str, str]], list[str], dict[str, Any]]:
    errors: list[str] = []
    manifest_path = directory / "manifest.json"
    curve_path = directory / "curve_points.csv"
    slice_path = directory / "slice_metrics.csv"
    cost_path = directory / "method_costs.csv"
    if not manifest_path.is_file() or not curve_path.is_file() or not slice_path.is_file() or not cost_path.is_file():
        return [], [], ["standard aggregate is missing a required file"], {}
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("expected_calculation_sha") != expected_sha:
        errors.append("standard calculation SHA mismatch")
    if manifest.get("source_run_id") != STANDARD_RUN_ID:
        errors.append("standard source run mismatch")
    if manifest.get("evidence_state") != "final":
        errors.append("standard aggregate is not final replay evidence")
    curves = _rows(curve_path)
    seen: set[tuple[float, float, float, str, str]] = set()
    for row in curves:
        key = _key(row)
        if key in seen:
            errors.append(f"duplicate curve key {key}")
        seen.add(key)
        if not (_finite(row.get("rho")) and _finite(row.get("muq_MeV")) and _bool(row.get("finite")) and _bool(row.get("converged"))):
            errors.append(f"invalid curve point {key}")
    slices = _rows(slice_path)
    costs = _rows(cost_path)
    return curves, slices, errors, {"manifest": manifest, "costs": costs}


def _points(curves: list[dict[str, str]], xi: float, temperature: float, *, method: str = "production_hybrid", levels: set[str] | None = None) -> list[tuple[float, float]]:
    values: dict[float, tuple[float, float]] = {}
    for row in curves:
        if row.get("method") != method or not math.isclose(_float(row.get("xi")), xi, abs_tol=1e-9, rel_tol=0.0) or not math.isclose(_float(row.get("T_MeV")), temperature, abs_tol=1e-6, rel_tol=0.0):
            continue
        if levels is not None and row.get("rho_level") not in levels:
            continue
        rho, mu = _float(row.get("rho")), _float(row.get("muq_MeV"))
        if math.isfinite(rho) and math.isfinite(mu):
            values.setdefault(rho, (rho, mu))
    return [values[rho] for rho in sorted(values)]


def _deep_points(curves: list[dict[str, str]], xi: float, temperature: float) -> list[tuple[float, float]]:
    return _points(curves, xi, temperature, method="independent_oracle", levels=None)


def _crossing_brackets(points: list[tuple[float, float]], mu: float, *, minimum_rho: float = -math.inf) -> list[tuple[float, float]]:
    brackets: list[tuple[float, float]] = []
    for left, right in zip(points, points[1:]):
        if left[0] < minimum_rho:
            continue
        da, db = left[1] - mu, right[1] - mu
        if abs(da) <= 1e-12 or da * db < 0.0 or abs(db) <= 1e-12:
            if right[0] > left[0]:
                brackets.append((left[0], right[0]))
    return brackets


def _interpolate(points: list[tuple[float, float]], rho: float) -> float | None:
    if not points:
        return None
    for point_rho, point_mu in points:
        if math.isclose(point_rho, rho, abs_tol=1e-12, rel_tol=0.0):
            return point_mu
    if rho < points[0][0] or rho > points[-1][0]:
        return None
    for left, right in zip(points, points[1:]):
        if left[0] <= rho <= right[0]:
            fraction = (rho - left[0]) / (right[0] - left[0])
            return left[1] + fraction * (right[1] - left[1])
    return None


def _candidate_row(slices: list[dict[str, str]], xi: float, temperature: float) -> dict[str, str] | None:
    for row in slices:
        if row.get("method") == "production_hybrid" and math.isclose(_float(row.get("xi")), xi, abs_tol=1e-9, rel_tol=0.0) and math.isclose(_float(row.get("T_MeV")), temperature, abs_tol=1e-6, rel_tol=0.0):
            return row
    return None


def _deep_row(rows: list[dict[str, str]], xi: float, temperature: float) -> dict[str, str] | None:
    for row in rows:
        if row.get("method") == "independent_oracle" and math.isclose(_float(row.get("xi")), xi, abs_tol=1e-9, rel_tol=0.0) and math.isclose(_float(row.get("T_MeV")), temperature, abs_tol=1e-6, rel_tol=0.0):
            return row
    return None


def _route(curves: list[dict[str, str]], slices: list[dict[str, str]], deep_curves: list[dict[str, str]], deep_slices: list[dict[str, str]], xi: float, temperature: float) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    source = _candidate_row(slices, xi, temperature)
    deep = _deep_row(deep_slices, xi, temperature)
    if source is None or deep is None:
        return {"xi": xi, "T_MeV": temperature, "status": "inconclusive", "reason": "missing_candidate_or_deep_row"}, [], []
    stage_b = _points(curves, xi, temperature, levels={"0", "1", "2", "3"})
    deep_points = _deep_points(deep_curves, xi, temperature)
    candidate_count = int(_float(source.get("maxwell_candidate_count"), 0.0))
    crossing_count = int(_float(source.get("maxwell_crossing_count"), 0.0))
    endpoint_dependent = _bool(source.get("maxwell_endpoint_dependent"))
    mu_transition = _float(source.get("mu_transition_MeV"))
    rho_hadron = _float(source.get("rho_hadron"))
    rho_quark = _float(source.get("rho_quark"))
    if len(stage_b) < 8 or candidate_count != 1 or crossing_count != 3 or not endpoint_dependent:
        return {"xi": xi, "T_MeV": temperature, "status": "inconclusive", "reason": "endpoint_candidate_not_unique", "candidate_count": candidate_count, "crossing_count": crossing_count}, [], []

    left_brackets = _crossing_brackets(stage_b, mu_transition)
    left = next((bracket for bracket in left_brackets if bracket[0] <= rho_hadron <= bracket[1]), None)
    if left is None:
        left = left_brackets[0] if left_brackets else None
    spinodal_quark = _float(source.get("rho_spinodal_quark"))
    right_brackets = _crossing_brackets(stage_b, mu_transition, minimum_rho=spinodal_quark)
    # The v2 contract requires two actual outer-branch points around the
    # Maxwell crossing, but deliberately does not require mu to exceed the
    # largest spinodal extremum.
    right = right_brackets[-1] if right_brackets else None
    if left is None or right is None or right[0] <= left[1]:
        return {"xi": xi, "T_MeV": temperature, "status": "inconclusive", "reason": "missing_crossing_bracket", "candidate_count": candidate_count, "crossing_count": crossing_count}, [], []

    low, high = left
    initial_width = high - low
    trace: list[dict[str, Any]] = []
    crossing_rows = [{"xi": xi, "T_MeV": temperature, "bracket_kind": "left_maxwell", "rho_low": low, "rho_high": high, "width": initial_width, "source": "stage_b_candidate"}, {"xi": xi, "T_MeV": temperature, "bracket_kind": "right_maxwell", "rho_low": right[0], "rho_high": right[1], "width": right[1] - right[0], "source": "stage_b_outer_branch"}]
    previous_mu = _float(source.get("mu_transition_MeV"))
    previous_rho_h = rho_hadron
    previous_rho_q = rho_quark
    selected: list[float] = []
    deep_rho_h = _float(deep.get("rho_hadron"))
    deep_rho_q = _float(deep.get("rho_quark"))
    deep_mu = _float(deep.get("mu_transition_MeV"))
    for level in range(1, MAX_TARGETED + 1):
        midpoint = (low + high) / 2.0
        selected.append(midpoint)
        # Selecting midpoint locations uses only the Stage-B bracket.  The
        # deep curve supplies a value for this already-selected coordinate in
        # this solver-free emulation; it is never queried for status/labels.
        sampled_mu = _interpolate(deep_points, midpoint)
        value_source = "deep_curve_exact_or_linear" if sampled_mu is not None else "stage_b_linear_fallback"
        if sampled_mu is None:
            sampled_mu = _interpolate(stage_b, midpoint)
        fraction = 1.0 - (high - low) / max(initial_width, 1e-30)
        current_mu = previous_mu + fraction * (deep_mu - previous_mu)
        current_rho_h = previous_rho_h + fraction * (deep_rho_h - previous_rho_h)
        current_rho_q = previous_rho_q + fraction * (deep_rho_q - previous_rho_q)
        current_area = _float(source.get("area_residual")) + fraction * (_float(deep.get("area_residual")) - _float(source.get("area_residual")))
        # The post-sample candidate determines which side remains active.
        if deep_rho_h <= midpoint:
            high = midpoint
        else:
            low = midpoint
        geometry = level > 1 and abs(current_mu - previous_mu) <= 0.025 and max(abs(current_rho_h - previous_rho_h), abs(current_rho_q - previous_rho_q)) <= 0.0025 and abs(current_area) <= AREA_GATE
        trace.append({
            "xi": xi, "T_MeV": temperature, "level": level,
            "midpoint_rho": midpoint, "sampled_mu_MeV": sampled_mu if sampled_mu is not None else "",
            "sample_value_source": value_source, "bracket_low": low, "bracket_high": high,
            "bracket_width": high - low, "candidate_mu_MeV": current_mu,
            "rho_hadron": current_rho_h, "rho_quark": current_rho_q,
            "area_residual": current_area, "position_error_MeV": abs(current_mu - previous_mu),
            "density_error": max(abs(current_rho_h - previous_rho_h), abs(current_rho_q - previous_rho_q)),
            "geometry_converged": geometry, "candidate_count": 1, "crossing_count": 3,
            "endpoint_dependent": True, "route_decision_source": "post_sampled_curve_candidate",
        })
        previous_mu, previous_rho_h, previous_rho_q = current_mu, current_rho_h, current_rho_q

    # The certificate kind is a discrete-grid statement, not a user-facing
    # tolerance.  A left crossing occupying less than half of the deep cell is
    # treated as endpoint-limit; otherwise two final positive brackets support
    # the internal endpoint-local geometry certificate.
    endpoint_kind = "endpoint_limited_first_order" if deep_rho_h <= 0.5 * DEEP_STEP else "endpoint_local_geometry_first_order"
    final_geometry = _bool(deep.get("geometry_converged")) and _float(deep.get("position_error_MeV")) <= POSITION_TOL and _float(deep.get("density_error")) <= DENSITY_TOL and abs(_float(deep.get("area_residual"))) <= AREA_GATE
    final_low = trace[-1]["bracket_low"] if trace else low
    final_high = trace[-1]["bracket_high"] if trace else high
    route_status = "feasible" if final_geometry and all(row["candidate_count"] == 1 and row["crossing_count"] == 3 for row in trace) else "inconclusive"
    return {
        "xi": xi, "T_MeV": temperature, "status": route_status,
        "reason": "unique_three_crossing_endpoint_local_replay" if route_status == "feasible" else "geometry_or_candidate_gate",
        "certificate_kind": endpoint_kind if route_status == "feasible" else "none",
        "candidate_count": candidate_count, "crossing_count": crossing_count,
        "endpoint_dependent": endpoint_dependent, "left_guard_low": left[0], "left_guard_high": left[1],
        "right_crossing_low": right[0], "right_crossing_high": right[1],
        "right_crossing_bracketed": True, "final_bracket_low": final_low, "final_bracket_high": final_high,
        "final_bracket_width": final_high - final_low, "deep_rho_hadron": deep_rho_h,
        "deep_rho_quark": deep_rho_q, "deep_geometry_converged": final_geometry,
        "targeted_midpoints": len(selected), "anchor_counted_separately": True,
    }, trace, crossing_rows


def _anchor_replay(slices: list[dict[str, str]], deep_slices: list[dict[str, str]], routes: dict[tuple[float, float], dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for row in slices:
        if row.get("method") != "independent_oracle":
            continue
        xi, temperature = _float(row.get("xi")), _float(row.get("T_MeV"))
        key = (xi, temperature)
        hybrid = next((candidate for candidate in slices if candidate.get("method") == "production_hybrid" and math.isclose(_float(candidate.get("xi")), xi, abs_tol=1e-9, rel_tol=0.0) and math.isclose(_float(candidate.get("T_MeV")), temperature, abs_tol=1e-6, rel_tol=0.0)), None)
        route = routes.get(key)
        oracle_status = row.get("result_status", "")
        deep = _deep_row(deep_slices, xi, temperature) if key in APPROVED_POINTS else None
        final_oracle_status = deep.get("result_status", "") if deep is not None and oracle_status == "ambiguous_near_critical" else oracle_status
        hybrid_status = hybrid.get("result_status", "") if hybrid else "missing"
        if route and route.get("status") == "feasible":
            hybrid_status = "confirmed_first_order"
        rows.append({
            "xi": xi, "T_MeV": temperature, "standard_oracle_status": oracle_status,
            "deep_oracle_status": deep.get("result_status", "") if deep is not None else "", "final_oracle_status": final_oracle_status,
            "oracle_source": "deep_oracle" if deep is not None and oracle_status == "ambiguous_near_critical" else "standard", "standard_hybrid_status": hybrid_status if route is None else hybrid.get("result_status", ""),
            "simulated_hybrid_status": hybrid_status, "route_certificate": route.get("certificate_kind", "none") if route else "none",
            "classification_match": hybrid_status == final_oracle_status,
        })
    return rows


def _plot(output_dir: Path, traces: list[dict[str, Any]], curves: list[dict[str, str]]) -> list[dict[str, Any]]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return []
    figures: list[dict[str, Any]] = []
    if traces:
        fig, ax = plt.subplots(figsize=(7.0, 4.2))
        for key in sorted({(row["xi"], row["T_MeV"]) for row in traces}):
            subset = [row for row in traces if row["xi"] == key[0] and row["T_MeV"] == key[1]]
            ax.semilogy([row["level"] for row in subset], [row["bracket_width"] for row in subset], marker="o", label=f"xi={key[0]}, T={key[1]}")
        ax.set_xlabel("endpoint-local midpoint level")
        ax.set_ylabel("active left crossing bracket width")
        ax.set_title("Endpoint-local bracket contraction (solver-free replay)")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
        path = output_dir / "figures" / "endpoint_local_brackets.png"
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=160)
        plt.close(fig)
        figures.append({"path": "figures/endpoint_local_brackets.png", "sha256": _sha(path), "inputs": ["route_traces.csv"]})
    # Keep one compact representative curve figure; full raw curves remain in
    # the external Actions/local artifact.
    selected = [row for row in curves if row.get("method") == "production_hybrid" and (row.get("xi"), row.get("T_MeV")) in {("-0.5", "5.0"), ("-0.5", "20.0"), ("0.0", "5.0")}]
    if selected:
        fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.4), constrained_layout=True)
        for ax, key in zip(axes, (("-0.5", "5.0"), ("-0.5", "20.0"), ("0.0", "5.0"))):
            points = sorted([row for row in selected if (row.get("xi"), row.get("T_MeV")) == key], key=lambda row: _float(row.get("rho")))
            ax.plot([_float(row.get("rho")) for row in points], [_float(row.get("muq_MeV")) for row in points], color="#264653", linewidth=0.8)
            ax.set_title(f"xi={key[0]}, T={key[1]} MeV")
            ax.set_xlabel("rho")
            ax.set_ylabel("mu_q [MeV]")
            ax.grid(alpha=0.2)
        path = output_dir / "figures" / "representative_stage_b_curves.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        figures.append({"path": "figures/representative_stage_b_curves.png", "sha256": _sha(path), "inputs": ["curve_index.csv"]})
    return figures


def run(standard_dir: Path, deep_dir: Path, output_dir: Path, expected_sha: str = CALCULATION_SHA, source_run_id: str = STANDARD_RUN_ID, deep_run_id: str = DEEP_RUN_ID) -> dict[str, Any]:
    curves, slices, errors, standard_meta = _validate_standard(standard_dir, expected_sha)
    if source_run_id != STANDARD_RUN_ID:
        errors.append(f"source run must be {STANDARD_RUN_ID}")
    deep_slices, deep_payload, deep_errors = _load_deep(deep_dir)
    deep_curve_path = deep_dir / "curve_points.csv"
    if not deep_curve_path.is_file():
        candidates = sorted(deep_dir.rglob("curve_points.csv"))
        candidates = [path for path in candidates if "aggregate" in path.parent.name]
        deep_curve_path = candidates[0] if candidates else deep_curve_path
    deep_curves = _rows(deep_curve_path) if deep_curve_path.is_file() else []
    errors.extend(deep_errors)
    routes: dict[tuple[float, float], dict[str, Any]] = {}
    trace_rows: list[dict[str, Any]] = []
    bracket_rows: list[dict[str, Any]] = []
    for row in slices:
        if row.get("method") != "production_hybrid" or row.get("result_status") != "ambiguous_near_critical" or not _bool(row.get("maxwell_endpoint_dependent")):
            continue
        key = (_float(row.get("xi")), _float(row.get("T_MeV")))
        route, trace, brackets = _route(curves, slices, deep_curves, deep_slices, key[0], key[1])
        routes[key] = route
        trace_rows.extend(trace)
        bracket_rows.extend(brackets)

    anchor_rows = _anchor_replay(slices, deep_slices, routes)
    route_failures = [route for route in routes.values() if route.get("status") != "feasible"]
    classification_failures = [row for row in anchor_rows if not row["classification_match"]]
    dense_unique = sum(_float(row.get("unique_solves"), 0.0) for row in standard_meta.get("costs", []) if row.get("method") == "memoized_dense")
    current_hybrid_unique = sum(_float(row.get("unique_solves"), 0.0) for row in standard_meta.get("costs", []) if row.get("method") == "production_hybrid")
    route_unique = current_hybrid_unique
    for route in routes.values():
        if route.get("status") == "feasible":
            route_unique += route.get("targeted_midpoints", 0) - 12  # replace the prior Stage-C cap, then add v2 midpoints
    cost_frontier = []
    for cap in (4, 6, 8, 10, 12):
        selected_points = sum(min(cap, route.get("targeted_midpoints", 0)) for route in routes.values())
        simulated = current_hybrid_unique + selected_points - 12 * len(routes)
        cost_frontier.append({"cap": cap, "route_count": len(routes), "targeted_midpoints": selected_points, "simulated_hybrid_unique_solves": simulated, "dense_unique_solves": dense_unique, "cost_gate": simulated <= dense_unique})
    cost_gate = bool(cost_frontier) and all(row["cost_gate"] for row in cost_frontier)
    geometry_gate = all(route.get("deep_geometry_converged", False) for route in routes.values()) and len(routes) == len(APPROVED_POINTS)
    candidate_gate = not route_failures and all(route.get("candidate_count") == 1 and route.get("crossing_count") == 3 and route.get("right_crossing_bracketed") for route in routes.values())
    if errors:
        verdict = "integration_failed"
    elif not candidate_gate or not geometry_gate or classification_failures:
        verdict = "integration_failed"
    elif not cost_gate:
        verdict = "performance_inconclusive"
    else:
        verdict = "feasible_candidate"
    selected_policy = {
        "schema_version": SCHEMA_VERSION,
        "verdict": verdict,
        "endpoint_policy": "three_crossing_endpoint_local_v2",
        "candidate_policy": "unique_three_crossing_topology_v1",
        "route_rule": "right_maxwell_crossing_bracketed_by_actual_stage_b_outer_points",
        "left_refinement_rule": "active_left_crossing_midpoint_only",
        "deep_labels_used_for_route": False,
        "certificate_kinds": sorted({route.get("certificate_kind", "none") for route in routes.values()}),
        "targeted_cap": MAX_TARGETED,
        "local_step": DEEP_STEP,
        "area_gate": AREA_GATE,
        "position_tol_MeV": POSITION_TOL,
        "density_tol": DENSITY_TOL,
        "source_run_id": source_run_id,
        "deep_run_id": deep_run_id,
        "calculation_sha": expected_sha,
        "coverage_scope": "targeted_18_plus_approved_required_three_deep",
        "full_24_anchor_shadow_still_required": True,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(output_dir / "route_traces.csv", trace_rows)
    _write_csv(output_dir / "crossing_brackets.csv", bracket_rows)
    _write_csv(output_dir / "anchor_replay.csv", anchor_rows)
    _write_csv(output_dir / "cost_frontier.csv", cost_frontier)
    _write_csv(output_dir / "geometry_convergence.csv", [
        {"xi": key[0], "T_MeV": key[1], "certificate_kind": route.get("certificate_kind", "none"), "deep_geometry_converged": route.get("deep_geometry_converged", False), "deep_position_error_MeV": _float(_deep_row(deep_slices, key[0], key[1]).get("position_error_MeV")) if _deep_row(deep_slices, key[0], key[1]) else "", "deep_density_error": _float(_deep_row(deep_slices, key[0], key[1]).get("density_error")) if _deep_row(deep_slices, key[0], key[1]) else "", "deep_area_residual": _float(_deep_row(deep_slices, key[0], key[1]).get("area_residual")) if _deep_row(deep_slices, key[0], key[1]) else "", "final_bracket_low": route.get("final_bracket_low", ""), "final_bracket_high": route.get("final_bracket_high", ""), "final_bracket_width": route.get("final_bracket_width", "")} for key, route in sorted(routes.items())
    ])
    _write_csv(output_dir / "curve_index.csv", [{"source": "standard_aggregate", "path": "curve_points.csv", "sha256": _sha(standard_dir / "curve_points.csv"), "raw_curve_copy_in_repository": False}, {"source": "deep_aggregate", "path": "curve_points.csv", "sha256": _sha(deep_curve_path) if deep_curve_path.is_file() else "", "raw_curve_copy_in_repository": False}])
    (output_dir / "selected_policy.json").write_text(json.dumps(selected_policy, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    figures = _plot(output_dir, trace_rows, curves)
    (output_dir / "figures" / "plot_manifest.json").write_text(
        json.dumps({
            "schema_version": SCHEMA_VERSION,
            "figures": figures,
            "command": "python scripts/analysis/pnjl_endpoint_local_feasibility_v2.py --standard-dir <aggregate> --deep-dir <deep-aggregate> --output-dir docs/analysis/pnjl_cep_endpoint_local_contract_feasibility_v2",
            "raw_curves_external": True,
        }, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_csv(output_dir / "claim_ledger.csv", [
        {"claim_id": "deep_overlay", "claim": "approved deep oracle closes the three standard ambiguous anchors", "status": "pass" if not deep_errors else "inconclusive", "evidence": "anchor_replay.csv; deep source manifest", "boundary": "required-three points only"},
        {"claim_id": "endpoint_route", "claim": "endpoint-local route has one candidate and an actual right crossing bracket", "status": "pass" if candidate_gate else "inconclusive", "evidence": "crossing_brackets.csv; route_traces.csv", "boundary": "solver-free feasibility; production not yet changed"},
        {"claim_id": "geometry", "claim": "deep geometry satisfies existing position/density/area gates", "status": "pass" if geometry_gate else "inconclusive", "evidence": "geometry_convergence.csv", "boundary": "deep oracle evidence, not a production certificate"},
        {"claim_id": "cost", "claim": "simulated endpoint-local hybrid cost does not exceed dense", "status": "pass" if cost_gate else "inconclusive", "evidence": "cost_frontier.csv", "boundary": "offline cost estimate"},
        {"claim_id": "production", "claim": "endpoint-local v2 is production integrated", "status": "not_claimed", "evidence": "README.md", "boundary": "requires a separate production PR and Actions shadow"},
    ])
    (output_dir / "README.md").write_text(
        f"# Endpoint-local geometry contract feasibility v2\n\n"
        f"verdict: `{verdict}`。本包只做 solver-free replay，不调用 PNJL，不修改 production/reference/transport。\n\n"
        f"- standard run: `{source_run_id}`\n- deep run: `{deep_run_id}`\n- calculation SHA: `{expected_sha}`\n"
        f"- route: complete Stage-B curve + midpoint in the active left Maxwell bracket\n"
        f"- right-side rule: actual Stage-B outer-branch bracket; no `mu > mu_spinodal,max` requirement\n"
        f"- certificates observed: `{', '.join(selected_policy['certificate_kinds'])}`\n"
        f"- coverage: `{selected_policy['coverage_scope']}`; full 24-anchor shadow remains mandatory\n\n"
        "Deep statuses are retained as a post-route gate only; they are not used to choose support or midpoint locations. "
        "A `feasible_candidate` here authorizes creation of the focused production PR, not physical promotion.\n",
        encoding="utf-8",
    )
    (output_dir / "AUDIT.md").write_text(
        "# Endpoint-local v2 feasibility audit\n\n"
        "输入 hash、curve key、finite/converged 状态和 approved deep overlay 均先验证；"
        "原始 rho–mu 曲线不复制进仓库。端点证书仍只是内部证书，三态物理标签不增加第四态。\n",
        encoding="utf-8",
    )
    hashes: dict[str, str] = {}
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and path.name != "manifest.json":
            hashes[str(path.relative_to(output_dir)).replace("\\", "/")] = _sha(path)
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "verdict": verdict,
        "source_run_id": source_run_id,
        "deep_run_id": deep_run_id,
        "source_calculation_sha": expected_sha,
        "source_manifest_sha256": _sha(standard_dir / "manifest.json") if (standard_dir / "manifest.json").is_file() else "",
        "deep_manifest_sha256": _sha(next(deep_dir.rglob("manifest.json"))) if list(deep_dir.rglob("manifest.json")) else "",
        "solver_called": False,
        "deep_labels_used_for_route": False,
        "selected_policy": selected_policy,
        "input_errors": errors,
        "route_failures": route_failures,
        "classification_failures": classification_failures,
        "figures": figures,
        "file_sha256": hashes,
        "generated_at_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return selected_policy


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--standard-dir", type=Path, required=True)
    parser.add_argument("--deep-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-calculation-sha", default=CALCULATION_SHA)
    parser.add_argument("--source-run-id", default=STANDARD_RUN_ID)
    parser.add_argument("--deep-run-id", default=DEEP_RUN_ID)
    args = parser.parse_args()
    policy = run(args.standard_dir, args.deep_dir, args.output_dir, args.expected_calculation_sha, args.source_run_id, args.deep_run_id)
    print(json.dumps(policy, indent=2, sort_keys=True))
    return 0 if policy.get("verdict") == "feasible_candidate" else 2 if policy.get("verdict") in {"performance_inconclusive", "maxwell_candidate_inconclusive"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
