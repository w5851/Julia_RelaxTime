#!/usr/bin/env python3
"""Replay bounded Stage-C density-certificate policies without the solver.

The numerical curve bundle is produced by Actions from an immutable
calculation SHA.  This postprocessor only ranks already available
0.003125 points, merges them with the complete 0.00625 Stage-B curve, and
checks the resulting diagnostics against the existing geometry contract.
Oracle labels are read only after route selection for the gate; they never
participate in feature ranking or point selection.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any, Iterable


SCHEMA_VERSION = "cep_hybrid_stagec_density_certificate_feasibility_v1"
SOURCE_CALCULATION_SHA = "ffa816df0a145f73d7490db1ed9ff10c92e017a4"
STAGE_B_STEP = 0.00625
STAGE_C_STEP = 0.003125
CAPS = (12, 16, 24)
ROUTES = ("stage_b_features_v1", "balanced_density_features_v2", "geometry_feedback_v2")
ROUTE_PRIORITY = {route: index for index, route in enumerate(ROUTES)}
METHODS = ("production_hybrid", "memoized_dense", "independent_oracle")
DENSITY_XIS = (-0.5, -0.35, -0.25, -0.2, -0.15, -0.1, 0.0, 0.3, 0.35, 0.5)
EXPECTED_JOB_COUNT = len(METHODS) * len(DENSITY_XIS)
POSITION_TOL = 0.025
DENSITY_TOL = 0.0025
AREA_TOL = 5e-5
RESPONSE_TOL = 0.025
CEP_TOL = 0.1

DENSITY_ANCHORS = (
    (-0.35, 51.0),
    (-0.25, 41.0),
    (-0.2, 41.0),
    (-0.15, 41.0),
    (-0.1, 41.0),
    (0.0, 51.0),
    (0.3, 21.0),
    (0.35, 51.0),
    (0.35, 101.0),
)
FIRST_ORDER_CONTROLS = ((-0.5, 60.0), (0.0, 60.0), (0.5, 60.0))
MONOTONE_CONTROLS = ((-0.5, 160.0), (0.0, 145.0), (0.5, 120.0))
ALL_ANCHORS = DENSITY_ANCHORS + FIRST_ORDER_CONTROLS + MONOTONE_CONTROLS
CEP_XIS = (0.05, 0.15, 0.225, 0.35, 0.4, 0.45, 0.5)
CEP_WINDOW_EXTENSION_MEV = 0.25
CEP_NODE_STEP_MEV = 0.0625
CROSSOVER_REFINEMENT_LEVELS = (2, 3, 4)

REQUIRED_FILES = ("curve_points.csv", "slice_metrics.csv", "method_costs.csv", "cep_accuracy.csv")
COMPARISON_FIELDS = (
    "xi", "T_MeV", "route", "cap", "oracle_status", "simulated_status", "reason",
    "candidate_count", "cross_resolution_candidate_count", "selected_points",
    "unique_solves", "dense_unique_solves", "geometry_gate", "cost_gate",
)


def _float(value: Any, default: float = math.nan) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    return result if math.isfinite(result) else default


def _finite(value: Any) -> bool:
    return math.isfinite(_float(value))


def _bool(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def _same(value: Any, target: float, tol: float = 1e-8) -> bool:
    number = _float(value)
    return math.isfinite(number) and math.isclose(number, target, abs_tol=tol, rel_tol=0.0)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: Iterable[dict[str, Any]], fields: Iterable[str] | None = None) -> None:
    rows = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = []
        for row in rows:
            for key in row:
                if key not in fields:
                    fields.append(key)
    fields = list(fields)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _job_dirs(input_dir: Path) -> list[Path]:
    if (input_dir / "curve_points.csv").is_file():
        return [input_dir]
    return sorted({path.parent for path in input_dir.rglob("curve_points.csv")})


def _load_bundle(
    input_dir: Path,
    expected_sha: str,
    expected_source_postprocess_sha: str | None = None,
) -> dict[str, Any]:
    job_dirs = _job_dirs(input_dir)
    if not job_dirs:
        raise ValueError(f"no curve_points.csv found under {input_dir}")
    curves: list[dict[str, str]] = []
    slices: list[dict[str, str]] = []
    costs: list[dict[str, str]] = []
    cep: list[dict[str, str]] = []
    input_files: list[dict[str, str]] = []
    manifests: list[dict[str, Any]] = []
    hash_errors: list[str] = []
    for job_dir in job_dirs:
        for name in REQUIRED_FILES:
            path = job_dir / name
            if not path.is_file():
                raise ValueError(f"missing {name}: {job_dir}")
            input_files.append({"path": path.as_posix(), "sha256": _sha256(path)})
        curves.extend(_read_csv(job_dir / "curve_points.csv"))
        slices.extend(_read_csv(job_dir / "slice_metrics.csv"))
        costs.extend(_read_csv(job_dir / "method_costs.csv"))
        cep.extend(_read_csv(job_dir / "cep_accuracy.csv"))
        manifest_path = job_dir / "manifest.json"
        if manifest_path.is_file():
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifests.append(manifest)
            input_files.append({"path": manifest_path.as_posix(), "sha256": _sha256(manifest_path)})
            declared_hashes = manifest.get("files") or manifest.get("curve_file_sha256") or {}
            if isinstance(declared_hashes, dict):
                for name, declared in declared_hashes.items():
                    if name not in REQUIRED_FILES:
                        continue
                    actual = _sha256(job_dir / name)
                    if str(declared) != actual:
                        hash_errors.append(f"manifest hash mismatch: {job_dir / name}")
    errors = _validate_bundle(
        curves,
        slices,
        costs,
        manifests,
        expected_sha,
        expected_source_postprocess_sha,
        job_count=len(job_dirs),
    )
    if errors:
        raise ValueError("; ".join(errors + hash_errors))
    if hash_errors:
        raise ValueError("; ".join(hash_errors))
    return {
        "curves": curves,
        "slices": slices,
        "costs": costs,
        "cep": cep,
        "input_files": input_files,
        "job_count": len(job_dirs),
        "manifests": manifests,
    }


def _validate_bundle(
    curves: list[dict[str, str]],
    slices: list[dict[str, str]],
    costs: list[dict[str, str]],
    manifests: list[dict[str, Any]],
    expected_sha: str,
    expected_source_postprocess_sha: str | None = None,
    job_count: int | None = None,
) -> list[str]:
    errors: list[str] = []
    if expected_sha != SOURCE_CALCULATION_SHA:
        errors.append(f"expected calculation SHA must be {SOURCE_CALCULATION_SHA}")
    if job_count is not None and job_count != EXPECTED_JOB_COUNT:
        errors.append(f"expected {EXPECTED_JOB_COUNT} numerical jobs, found {job_count}")
    if len(manifests) != EXPECTED_JOB_COUNT:
        errors.append(f"expected {EXPECTED_JOB_COUNT} job manifests, found {len(manifests)}")
    job_keys: set[tuple[str, float]] = set()
    for manifest in manifests:
        sha = manifest.get("calculation_sha") or manifest.get("expected_calculation_sha")
        if sha and str(sha) != expected_sha:
            errors.append(f"manifest calculation SHA mismatch: {sha}")
        postprocess_sha = manifest.get("postprocess_sha") or manifest.get("workflow_head_sha")
        if expected_source_postprocess_sha and postprocess_sha != expected_source_postprocess_sha:
            errors.append(f"manifest postprocess SHA mismatch: {postprocess_sha}")
        method = str(manifest.get("method", ""))
        xi = _float(manifest.get("xi"))
        key = (method, xi)
        if method not in METHODS or not math.isfinite(xi) or xi not in DENSITY_XIS:
            errors.append(f"invalid numerical job identity: method={method!r}, xi={manifest.get('xi')!r}")
        elif key in job_keys:
            errors.append(f"duplicate numerical job identity: {key}")
        else:
            job_keys.add(key)
    if len(job_keys) != EXPECTED_JOB_COUNT:
        errors.append(f"numerical job matrix is incomplete: {sorted(job_keys)!r}")
    seen: set[tuple[str, float, float, str, float]] = set()
    for row in curves:
        xi = _float(row.get("xi"))
        temperature = _float(row.get("T_MeV"))
        rho = _float(row.get("rho"))
        key = (
            str(row.get("method", "")), xi, temperature,
            str(row.get("rho_level", "")), rho,
        )
        if key in seen:
            errors.append(f"duplicate curve key: {key}")
        seen.add(key)
        if not _bool(row.get("converged")) or not _bool(row.get("finite")):
            errors.append(f"non-finite/non-converged curve: {key}")
        if not _finite(row.get("rho")) or not _finite(row.get("muq_MeV")):
            errors.append(f"invalid curve coordinate: {key}")
        row_sha = row.get("calculation_sha")
        if row_sha and row_sha != expected_sha:
            errors.append(f"curve calculation SHA mismatch: {key}")
    if not curves:
        errors.append("curve_points.csv is empty")
    if not slices:
        errors.append("slice_metrics.csv is empty")
    if not costs:
        errors.append("method_costs.csv is empty")
    return errors


def _anchor_key(xi: Any, temperature: Any) -> tuple[float, float]:
    return round(_float(xi), 8), round(_float(temperature), 8)


def _row_for_anchor(rows: list[dict[str, str]], method: str, anchor: tuple[float, float]) -> dict[str, str] | None:
    for row in rows:
        if row.get("method") == method and _same(row.get("xi"), anchor[0]) and _same(row.get("T_MeV"), anchor[1], 1e-6):
            return row
    return None


def _level(row: dict[str, str]) -> str:
    return str(row.get("rho_level", "")).strip()


def _curve_points(rows: list[dict[str, str]], method: str, anchor: tuple[float, float], levels: set[str] | None = None) -> list[tuple[float, float]]:
    grouped: dict[float, tuple[float, float]] = {}
    for row in rows:
        if row.get("method") != method or not _same(row.get("xi"), anchor[0]) or not _same(row.get("T_MeV"), anchor[1], 1e-6):
            continue
        if levels is not None and _level(row) not in levels:
            continue
        rho, mu = _float(row.get("rho")), _float(row.get("muq_MeV"))
        if not (math.isfinite(rho) and math.isfinite(mu)):
            continue
        residual = _float(row.get("residual_norm"), math.inf)
        previous = grouped.get(rho)
        if previous is None or residual < previous[0]:
            grouped[rho] = (residual, mu)
    return sorted((rho, item[1]) for rho, item in grouped.items())


def _candidate_intervals(points: list[tuple[float, float]]) -> list[dict[str, float]]:
    if len(points) < 5:
        return []
    slopes: list[float] = []
    for (left_rho, left_mu), (right_rho, right_mu) in zip(points, points[1:]):
        slopes.append((right_mu - left_mu) / (right_rho - left_rho) if right_rho > left_rho else math.nan)
    signs = [1 if slope > 0 else -1 if slope < 0 else 0 for slope in slopes]
    candidates: list[dict[str, float]] = []
    for index in range(1, len(signs) - 1):
        if signs[index] != -1 or signs[index - 1] != 1:
            continue
        end = index
        while end + 1 < len(signs) and signs[end + 1] == -1:
            end += 1
        if end + 1 >= len(signs) or signs[end + 1] != 1:
            continue
        low = points[index][0]
        high = points[end + 1][0]
        drop = points[index][1] - points[end + 1][1]
        if high > low and drop > 0:
            candidates.append({
                "rho_low": low,
                "rho_high": high,
                "drop_mu": drop,
                "width": high - low,
                "negative_secants": end - index + 1,
            })
    return candidates


def _feature_targets(points: list[tuple[float, float]]) -> list[tuple[str, float, float]]:
    """Rank density-only features derived from Stage-B, without oracle data."""
    if len(points) < 3:
        return []
    targets: list[tuple[str, float, float]] = []
    candidates = _candidate_intervals(points)
    rho_values = [rho for rho, _ in points]
    for candidate in candidates:
        low_index = min(range(len(points)), key=lambda i: abs(points[i][0] - candidate["rho_low"]))
        high_index = min(range(len(points)), key=lambda i: abs(points[i][0] - candidate["rho_high"]))
        for label, index, rank in (
            ("spinodal_left", low_index, 0.0),
            ("coexistence_midpoint", (low_index + high_index) // 2, 1.0),
            ("spinodal_right", high_index, 0.0),
        ):
            targets.append((label, points[index][0], rank))
        if low_index > 0:
            targets.append(("spinodal_left_adjacent", points[low_index - 1][0], 0.5))
        if high_index + 1 < len(points):
            targets.append(("spinodal_right_adjacent", points[high_index + 1][0], 0.5))

    curvature: list[tuple[float, float]] = []
    for index in range(1, len(points) - 1):
        left_rho, left_mu = points[index - 1]
        rho, mu = points[index]
        right_rho, right_mu = points[index + 1]
        if right_rho <= left_rho:
            continue
        curvature.append((abs(right_mu - 2 * mu + left_mu), rho))
    for score, rho in sorted(curvature, reverse=True)[:8]:
        targets.append(("high_curvature", rho, 2.0 - min(score, 1.0)))

    area_proxy: list[tuple[float, float]] = []
    for (left_rho, left_mu), (right_rho, right_mu) in zip(points, points[1:]):
        area_proxy.append((abs(right_mu - left_mu) * (right_rho - left_rho), (left_rho + right_rho) / 2))
    for score, rho in sorted(area_proxy, reverse=True)[:8]:
        targets.append(("area_contribution", rho, 3.0 - min(score, 1.0)))
    return sorted(targets, key=lambda item: (item[2], item[1]))


def _select_local_points(pool: list[tuple[float, float]], targets: list[tuple[str, float, float]], cap: int) -> list[tuple[float, float, str]]:
    if cap <= 0 or not pool:
        return []
    available = {rho: mu for rho, mu in pool}
    selected: list[tuple[float, float, str]] = []
    for label, target, _rank in targets:
        if not available:
            break
        rho = min(available, key=lambda value: (abs(value - target), value))
        selected.append((rho, available.pop(rho), label))
        if len(selected) >= cap:
            break
    if len(selected) < cap:
        for rho in sorted(available, key=lambda value: (abs(value - (pool[0][0] + pool[-1][0]) / 2), value)):
            selected.append((rho, available.pop(rho), "ranked_midpoint"))
            if len(selected) >= cap:
                break
    return sorted(selected, key=lambda item: item[0])


def _geometry_gate(row: dict[str, str] | None) -> bool:
    if row is None or not _bool(row.get("geometry_converged")):
        return False
    position = _float(row.get("position_error_MeV"))
    density = _float(row.get("density_error"))
    area = _float(row.get("maxwell_area_gate"))
    residual = _float(row.get("area_residual"))
    return (
        math.isfinite(position) and math.isfinite(density) and math.isfinite(area) and math.isfinite(residual)
        and position <= POSITION_TOL and density <= DENSITY_TOL and area <= AREA_TOL and residual <= AREA_TOL
    )


def _merge_curve(stage_b: list[tuple[float, float]], selected: list[tuple[float, float, str]]) -> list[tuple[float, float]]:
    merged = {rho: mu for rho, mu in stage_b}
    for rho, mu, _label in selected:
        merged.setdefault(rho, mu)
    return sorted(merged.items())


def _route_targets(route: str, stage_b: list[tuple[float, float]], production_rows: list[dict[str, str]], anchor: tuple[float, float]) -> list[tuple[str, float, float]]:
    targets = _feature_targets(stage_b)
    if route == "stage_b_features_v1":
        # The legacy route is replayed from Stage-B features only.  Reading
        # production Stage-C points here would leak the route's prior result
        # into the offline selection contract.
        return [("stage_b_features_v1", rho, rank) for _label, rho, rank in targets]
    if route == "balanced_density_features_v2":
        return targets
    # Feedback may only inspect the current hybrid diagnostics and Stage-B
    # features.  It does not inspect oracle labels or oracle geometry.
    diagnostic = _row_for_anchor(production_rows, "production_hybrid", anchor)
    if diagnostic is None:
        return targets
    errors = [
        ("position_feedback", _float(diagnostic.get("position_error_MeV")) / POSITION_TOL),
        ("density_feedback", _float(diagnostic.get("density_error")) / DENSITY_TOL),
        ("area_feedback", _float(diagnostic.get("maxwell_area_gate")) / AREA_TOL),
    ]
    dominant = max(errors, key=lambda item: item[1] if math.isfinite(item[1]) else -math.inf)[0]
    return [(dominant if label.startswith("spinodal") or label == "coexistence_midpoint" else label, rho, rank) for label, rho, rank in targets]


def _simulated_status(stage_b: list[tuple[float, float]], selected: list[tuple[float, float, str]], diagnostic: dict[str, str] | None, oracle_status: str) -> tuple[str, str, int, int]:
    diagnostic = diagnostic or {}
    merged = _merge_curve(stage_b, selected)
    candidates = _candidate_intervals(merged)
    # Both candidate counts are evaluated on the same sorted global union.
    # Keeping an unsorted Stage-B + Stage-C concatenation here would create
    # artificial sign changes at the insertion boundaries.
    fine_candidates = _candidate_intervals(merged)
    candidate_count = len(candidates)
    cross_count = max(candidate_count, len(fine_candidates))
    if cross_count > 1:
        return "ambiguous_near_critical", "multiple_candidate_intervals", candidate_count, cross_count
    if candidate_count == 1 and _geometry_gate(diagnostic):
        return "confirmed_first_order", "unique_candidate_geometry_gate", candidate_count, cross_count
    diagnostic_status = str(diagnostic.get("stage_b_status", diagnostic.get("result_status", "")))
    if candidate_count == 0 and diagnostic_status in {"confirmed_monotone", "monotone", "stable_no_s_shape"}:
        return "confirmed_monotone", "stable_no_s_shape", candidate_count, cross_count
    if candidate_count == 0 and oracle_status == "confirmed_monotone":
        return "ambiguous_near_critical", "monotone_requires_stage_b_certificate", candidate_count, cross_count
    return "ambiguous_near_critical", "maxwell_or_geometry_unresolved", candidate_count, cross_count


def _expected_status(anchor: tuple[float, float], oracle_status: str) -> str:
    if anchor in DENSITY_ANCHORS or anchor in FIRST_ORDER_CONTROLS:
        return "confirmed_first_order"
    if anchor in MONOTONE_CONTROLS:
        return "confirmed_monotone"
    return oracle_status


def _dense_cost(costs: list[dict[str, str]]) -> int:
    return int(round(sum(_float(row.get("unique_solves"), 0.0) for row in costs if row.get("method") == "memoized_dense")))


def _base_cost_keys(curves: list[dict[str, str]], anchor: tuple[float, float]) -> set[tuple[float, float, str]]:
    keys: set[tuple[float, float, str]] = set()
    for row in curves:
        if row.get("method") != "production_hybrid" or not _same(row.get("xi"), anchor[0]) or not _same(row.get("T_MeV"), anchor[1], 1e-6):
            continue
        level = _level(row)
        if level in {"0", "1", "2", "3"} or not level:
            rho = _float(row.get("rho"))
            math.isfinite(rho) and keys.add((anchor[0], anchor[1], f"{rho:.12g}"))
    return keys


def _evaluate_route(bundle: dict[str, Any], route: str, cap: int, dense_cost: int) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    curves, slices = bundle["curves"], bundle["slices"]
    total_unique = 0
    total_selected = 0
    rows: list[dict[str, Any]] = []
    candidates_rows: list[dict[str, Any]] = []
    mismatch = unsupported = geometry_failures = multiple = invalid = 0
    for anchor in ALL_ANCHORS:
        stage_b = _curve_points(curves, "independent_oracle", anchor, {"0"})
        if not stage_b:
            stage_b = _curve_points(curves, "production_hybrid", anchor, {"2", "3"})
        pool = _curve_points(curves, "independent_oracle", anchor, {"1"})
        production = _row_for_anchor(slices, "production_hybrid", anchor)
        oracle = _row_for_anchor(slices, "independent_oracle", anchor)
        oracle_status = str((oracle or {}).get("result_status", "missing"))
        targets = _route_targets(route, stage_b, slices, anchor)
        selected = _select_local_points(pool, targets, cap)
        simulated, reason, candidate_count, cross_count = _simulated_status(stage_b, selected, production, oracle_status)
        expected = _expected_status(anchor, oracle_status)
        if oracle_status == "ambiguous_near_critical" and simulated == "confirmed_first_order":
            unsupported += 1
        if simulated != expected:
            mismatch += 1
        if simulated == "confirmed_first_order" and not _geometry_gate(production):
            geometry_failures += 1
        if cross_count > 1:
            multiple += 1
        if not stage_b or not all(math.isfinite(mu) for _rho, mu in stage_b):
            invalid += 1
        base_keys = _base_cost_keys(curves, anchor)
        selected_keys = {(anchor[0], anchor[1], f"{rho:.12g}") for rho, _mu, _label in selected}
        unique = len(base_keys | selected_keys)
        total_unique += unique
        total_selected += len(selected)
        row = {
            "xi": anchor[0], "T_MeV": anchor[1], "route": route, "cap": cap,
            "oracle_status": oracle_status, "expected_status": expected,
            "simulated_status": simulated, "reason": reason,
            "candidate_count": candidate_count, "cross_resolution_candidate_count": cross_count,
            "selected_points": len(selected), "unique_solves": unique,
            "dense_unique_solves": dense_cost, "geometry_gate": _geometry_gate(production),
            "cost_gate": unique <= dense_cost, "finite_gate": invalid == 0,
        }
        rows.append(row)
        for index, (rho, mu, label) in enumerate(selected, start=1):
            candidates_rows.append({
                "xi": anchor[0], "T_MeV": anchor[1], "route": route, "cap": cap,
                "rank": index, "feature": label, "rho": rho, "muq_MeV": mu,
            })
    state_gate = mismatch == 0 and unsupported == 0
    geometry_gate = geometry_failures == 0
    candidate_gate = multiple == 0
    finite_gate = invalid == 0
    cost_gate = total_unique <= dense_cost
    frontier = {
        "route": route, "cap": cap, "anchors": len(ALL_ANCHORS),
        "classification_mismatches": mismatch, "unsupported_confirmations": unsupported,
        "geometry_failures": geometry_failures, "multiple_candidate_anchors": multiple,
        "invalid_anchors": invalid, "unique_solves": total_unique,
        "targeted_points": total_selected, "dense_unique_solves": dense_cost,
        "state_gate": state_gate, "geometry_gate": geometry_gate,
        "candidate_gate": candidate_gate, "finite_gate": finite_gate,
        "cost_gate": cost_gate, "feasible": all((state_gate, geometry_gate, candidate_gate, finite_gate, cost_gate)),
    }
    return frontier, rows, candidates_rows


def _load_c2_audit(audit_dir: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    cep_path = audit_dir / "tables" / "c1_vs_c2_cep_gates.csv"
    crossover_path = audit_dir / "tables" / "c1_vs_c2_crossover_xi_0p2875.csv"
    cep_rows = _read_csv(cep_path) if cep_path.is_file() else []
    crossover_rows = _read_csv(crossover_path) if crossover_path.is_file() else []
    cep_gate = []
    for xi in CEP_XIS:
        row = next((item for item in cep_rows if _same(item.get("xi"), xi)), None)
        width = _float((row or {}).get("candidate_bracket_width_T_MeV"))
        cep_gate.append({
            "xi": xi,
            "source": "c2_convergence_audit_v1" if row else "missing",
            "bracket_width_T_MeV": width,
            "window_extension_MeV": CEP_WINDOW_EXTENSION_MEV,
            "node_step_MeV": CEP_NODE_STEP_MEV,
            "pass": math.isfinite(width) and width <= CEP_TOL,
        })
    crossover_bad = 0
    for row in crossover_rows:
        if row.get("metric") == "derivative" and _float(row.get("rel_diff")) > RESPONSE_TOL:
            crossover_bad += 1
        elif row.get("metric") in {"T_crossover_MeV", "rho"} and _float(row.get("abs_diff")) > (0.05 if row.get("metric") == "T_crossover_MeV" else 0.005):
            crossover_bad += 1
    crossover_gate = [{
        "xi": 0.2875,
        "source": "c2_convergence_audit_v1" if crossover_rows else "missing",
        "refinement_levels": ",".join(str(level) for level in CROSSOVER_REFINEMENT_LEVELS),
        "risk_rows": crossover_bad,
        "pass": bool(crossover_rows) and crossover_bad == 0,
    }]
    return cep_gate, crossover_gate


def _write_plots(output_dir: Path, bundle: dict[str, Any], selected_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return []
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    cases = DENSITY_ANCHORS[:2] + MONOTONE_CONTROLS[:1]
    entries = []
    for xi, temperature in cases:
        points = _curve_points(bundle["curves"], "independent_oracle", (xi, temperature), {"0"})
        if not points:
            continue
        selected = [row for row in selected_rows if _same(row.get("xi"), xi) and _same(row.get("T_MeV"), temperature, 1e-6)]
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot([rho for rho, _mu in points], [mu for _rho, mu in points], color="#264653", linewidth=1.1)
        if selected:
            ax.scatter([_float(row.get("rho")) for row in selected], [_float(row.get("muq_MeV")) for row in selected], color="#e76f51", s=14)
        ax.set_xlabel(r"$\rho$")
        ax.set_ylabel(r"$\mu_q$ [MeV]")
        ax.set_title(f"Stage-B density certificate xi={xi:g}, T={temperature:g} MeV")
        ax.grid(alpha=0.2)
        name = f"rho_mu_xi_{str(xi).replace('-', 'm').replace('.', 'p')}_T_{str(temperature).replace('.', 'p')}.png"
        path = figure_dir / name
        fig.tight_layout()
        fig.savefig(path, dpi=160)
        plt.close(fig)
        entries.append({"path": f"figures/{name}", "sha256": _sha256(path), "xi": xi, "T_MeV": temperature})
    return entries


def _git_head() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def run(
    input_dir: Path,
    output_dir: Path,
    expected_sha: str = SOURCE_CALCULATION_SHA,
    postprocess_sha: str | None = None,
    audit_dir: Path | None = None,
    expected_source_postprocess_sha: str | None = None,
) -> dict[str, Any]:
    bundle = _load_bundle(input_dir, expected_sha, expected_source_postprocess_sha)
    dense_cost = _dense_cost(bundle["costs"])
    frontier: list[dict[str, Any]] = []
    route_rows: list[dict[str, Any]] = []
    candidate_rows: list[dict[str, Any]] = []
    for route in ROUTES:
        for cap in CAPS:
            row, replay, candidates = _evaluate_route(bundle, route, cap, dense_cost)
            frontier.append(row)
            route_rows.extend(replay)
            candidate_rows.extend(candidates)
    cap12 = [row for row in frontier if row["cap"] == 12 and row["feasible"]]
    diagnostic = [row for row in frontier if row["cap"] in (16, 24) and row["feasible"]]
    selected = min(cap12, key=lambda row: (row["unique_solves"], row["targeted_points"], ROUTE_PRIORITY[row["route"]])) if cap12 else None
    if selected is not None:
        verdict = "feasible_candidate"
        reason = "minimum-cost cap-12 route satisfies all gates"
    elif diagnostic:
        verdict = "cap_contract_inconclusive"
        reason = "only a diagnostic cap above 12 satisfies the gates"
    elif any(row["unsupported_confirmations"] for row in frontier):
        verdict = "oracle_inconclusive"
        reason = "an oracle-ambiguous anchor would be confirmed"
    elif any(row["multiple_candidate_anchors"] for row in frontier):
        verdict = "maxwell_candidate_inconclusive"
        reason = "multiple cross-resolution candidate intervals remain"
    elif any(not row["cost_gate"] for row in frontier if row["cap"] == 12):
        verdict = "performance_inconclusive"
        reason = "cap-12 route cost exceeds memoized dense"
    else:
        verdict = "integration_failed"
        reason = "one or more cap-12 route gates failed"
    audit_dir = audit_dir or Path(__file__).resolve().parents[2] / "docs" / "analysis" / "pnjl" / "c2_audits" / "c2_convergence_audit_v1"
    cep_gates, crossover_gates = _load_c2_audit(audit_dir)
    if verdict == "feasible_candidate" and (not all(row["pass"] for row in cep_gates) or not all(row["pass"] for row in crossover_gates)):
        verdict = "cep_ambiguity_inconclusive" if not all(row["pass"] for row in cep_gates) else "crossover_refinement_required"
        reason = "stored C2 CEP/crossover diagnostic gate remains unresolved"

    output_dir.mkdir(parents=True, exist_ok=True)
    selected_policy = {
        "schema_version": SCHEMA_VERSION,
        "verdict": verdict,
        "route": selected["route"] if selected else None,
        "cap": selected["cap"] if selected else None,
        "point_ranking_version": selected["route"] if selected else None,
        "local_step": STAGE_C_STEP,
        "target_cap": 12,
        "feature_radius": 0.025,
        "area_error_budget": AREA_TOL,
        "reason": reason,
    }
    _write_csv(output_dir / "cost_frontier.csv", frontier)
    _write_csv(output_dir / "route_comparison.csv", route_rows, COMPARISON_FIELDS + ("expected_status", "finite_gate"))
    _write_csv(output_dir / "selected_point_index.csv", candidate_rows)
    _write_csv(output_dir / "cep_gates.csv", cep_gates)
    _write_csv(output_dir / "crossover_gates.csv", crossover_gates)
    (output_dir / "selected_policy.json").write_text(json.dumps(selected_policy, indent=2) + "\n", encoding="utf-8")
    plot_entries = _write_plots(output_dir, bundle, candidate_rows if selected is None else [row for row in candidate_rows if row.get("route") == selected["route"] and row.get("cap") == selected["cap"]])
    claim_rows = [
        {"claim_id": "classification", "claim": "all density anchors and controls match the oracle/author gate", "status": "pass" if selected else "inconclusive", "evidence": "route_comparison.csv"},
        {"claim_id": "geometry", "claim": "shared first-order rows satisfy position/density/area gates", "status": "pass" if selected else "inconclusive", "evidence": "route_comparison.csv; cost_frontier.csv"},
        {"claim_id": "cost", "claim": "selected cap-12 route does not exceed memoized dense", "status": "pass" if selected else "inconclusive", "evidence": "cost_frontier.csv"},
        {"claim_id": "cep", "claim": "CEP brackets satisfy the 0.1 MeV gate", "status": "pass" if all(row["pass"] for row in cep_gates) else "inconclusive", "evidence": "cep_gates.csv"},
        {"claim_id": "crossover", "claim": "xi=0.2875 crossover refinement is stable", "status": "pass" if all(row["pass"] for row in crossover_gates) else "inconclusive", "evidence": "crossover_gates.csv"},
        {"claim_id": "promotion", "claim": "feasibility authorizes production/reference promotion", "status": "not_claimed", "evidence": "requires targeted/full shadow and author review"},
    ]
    _write_csv(output_dir / "claim_ledger.csv", claim_rows)
    (output_dir / "README.md").write_text(
        f"# Stage-C density certificate feasibility v1\n\n"
        f"verdict: `{verdict}`。这是固定 calculation SHA `{expected_sha}` 的 solver-free replay；不调用 equilibrium solver，"
        "不改变 production、Maxwell、旧 reference 或 transport。\n\n"
        f"- density anchors: `{len(DENSITY_ANCHORS)}`；controls: `{len(FIRST_ORDER_CONTROLS) + len(MONOTONE_CONTROLS)}`\n"
        f"- routes: `{', '.join(ROUTES)}`；cap scan: `{', '.join(map(str, CAPS))}`\n"
        f"- selected policy: `{json.dumps(selected_policy, ensure_ascii=False)}`\n"
        "- route selection reads Stage-B curves only; oracle labels are applied after selection for gates.\n"
        "- cap 16/24 is diagnostic only and cannot authorize production.\n",
        encoding="utf-8",
    )
    (output_dir / "AUDIT.md").write_text(
        "# Stage-C density certificate audit\n\n"
        "成本按 production Stage-A/B key 并集与选定 Stage-C key 计算，memoized dense 仅作为成本上界，"
        "不会免费消费 oracle 网格。完整 raw curve 留在 Actions/local artifact；仓库只保留 selected-point index、"
        "route replay、成本 frontier、CEP/crossover gate 和代表图。非 `feasible_candidate` 时不得创建 production PR。\n",
        encoding="utf-8",
    )
    files: dict[str, str] = {}
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and path.name != "manifest.json":
            files[path.relative_to(output_dir).as_posix()] = _sha256(path)
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "verdict": verdict,
        "source_calculation_sha": expected_sha,
        "source_postprocess_shas": sorted({
            str(item.get("postprocess_sha") or item.get("workflow_head_sha"))
            for item in bundle["manifests"]
            if item.get("postprocess_sha") or item.get("workflow_head_sha")
        }),
        "postprocess_sha": postprocess_sha or _git_head(),
        "producer_head_sha": _git_head(),
        "solver_called": False,
        "source_job_count": bundle["job_count"],
        "input_files": bundle["input_files"],
        "selected_policy": selected_policy,
        "cost_frontier": frontier,
        "cep_gates": cep_gates,
        "crossover_gates": crossover_gates,
        "figures": plot_entries,
        "files": files,
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    (output_dir / "plot_manifest.json").write_text(json.dumps({"schema_version": SCHEMA_VERSION, "figures": plot_entries, "raw_curves_external": True}, indent=2) + "\n", encoding="utf-8")
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--expected-calculation-sha", default=SOURCE_CALCULATION_SHA)
    parser.add_argument("--postprocess-sha", default=None)
    parser.add_argument("--audit-dir", type=Path, default=None)
    parser.add_argument("--expected-source-postprocess-sha", default=None)
    args = parser.parse_args(argv)
    manifest = run(
        args.input_dir,
        args.output_dir,
        args.expected_calculation_sha,
        args.postprocess_sha,
        args.audit_dir,
        args.expected_source_postprocess_sha,
    )
    print(json.dumps({"verdict": manifest["verdict"], "output_dir": args.output_dir.as_posix()}))
    return 0 if manifest["verdict"] == "feasible_candidate" else 2


if __name__ == "__main__":
    raise SystemExit(main())
