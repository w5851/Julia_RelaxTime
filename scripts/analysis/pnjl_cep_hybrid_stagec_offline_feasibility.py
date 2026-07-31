#!/usr/bin/env python3
"""Replay a bounded Stage-C policy without calling the PNJL solver.

The input is the immutable aggregate-replay artifact from the hybrid shadow.
This tool deliberately treats the independent-oracle ``rho_level=0`` curve as
an offline semantic reference for the complete 0.00625 grid.  It never writes
production data and never imports the equilibrium solver.  The estimated cost
is computed from the production-hybrid coarse/fine points plus selected local
verification points, so the oracle reference curve is not charged to the
simulated hybrid policy.
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


XIS = (-0.5, 0.0, 0.5)
CAPS = (48, 96, 160, 224)
REQUIRED_FILES = ("manifest.json", "curve_points.csv", "slice_metrics.csv", "method_costs.csv", "cep_accuracy.csv")
AREA_TOL = 5e-5
DENSITY_TOL = 0.0025
MU_TOL = 0.025
SLOPE_MARGIN = 0.01
MIN_NEGATIVE_SECANTS = 2


def _read_csv(path: Path) -> list[dict[str, str]]:
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
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _float(value: Any, default: float = math.nan) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    return number


def _finite(value: Any) -> bool:
    return math.isfinite(_float(value))


def _bool(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _key(row: dict[str, str], *, include_level: bool = True) -> tuple[str, ...]:
    fields = ("xi", "method", "T_MeV")
    if include_level:
        fields += ("rho_level",)
    fields += ("rho",)
    return tuple(row.get(field, "") for field in fields)


def _validate_input(input_dir: Path, expected_sha: str) -> list[str]:
    errors: list[str] = []
    for name in REQUIRED_FILES:
        if not (input_dir / name).is_file():
            errors.append(f"missing input file: {name}")
    if errors:
        return errors

    manifest = json.loads((input_dir / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("expected_calculation_sha") != expected_sha:
        errors.append("manifest calculation SHA does not match requested SHA")
    if manifest.get("evidence_state") != "final":
        errors.append("input evidence is not final aggregate replay")
    if manifest.get("gate", {}).get("verdict") != "oracle_inconclusive":
        errors.append("input replay verdict is not the expected oracle_inconclusive diagnostic")

    curves = _read_csv(input_dir / "curve_points.csv")
    seen: set[tuple[str, ...]] = set()
    for row in curves:
        key = _key(row)
        if key in seen:
            errors.append(f"duplicate curve key: {key}")
        seen.add(key)
        if not _bool(row.get("converged")) or not _bool(row.get("finite")):
            errors.append(f"non-finite/non-converged curve point: {key}")
        if not _finite(row.get("rho")) or not _finite(row.get("muq_MeV")):
            errors.append(f"invalid curve coordinate: {key}")
    if not curves:
        errors.append("curve_points.csv is empty")

    for name in ("slice_metrics.csv", "method_costs.csv", "cep_accuracy.csv"):
        rows = _read_csv(input_dir / name)
        if not rows:
            errors.append(f"{name} is empty")
    return errors


def _group_points(rows: list[dict[str, str]], *, method: str, xi: float, temperature: float, levels: set[str] | None = None) -> list[tuple[float, float]]:
    grouped: dict[float, list[tuple[float, float, float]]] = defaultdict(list)
    for row in rows:
        if row.get("method") != method or not math.isclose(_float(row.get("xi")), xi, abs_tol=1e-9, rel_tol=0.0):
            continue
        if not math.isclose(_float(row.get("T_MeV")), temperature, abs_tol=1e-5, rel_tol=0.0):
            continue
        if levels is not None and str(row.get("rho_level")) not in levels:
            continue
        rho, mu = _float(row.get("rho")), _float(row.get("muq_MeV"))
        if not (math.isfinite(rho) and math.isfinite(mu)):
            continue
        residual = _float(row.get("residual_norm"), 0.0)
        grouped[rho].append((mu, residual if math.isfinite(residual) else math.inf, float(len(grouped[rho]))))
    points: list[tuple[float, float]] = []
    for rho, candidates in grouped.items():
        mu, _, _ = min(candidates, key=lambda item: (item[1], item[0]))
        points.append((rho, mu))
    return sorted(points)


def _candidate_runs(points: list[tuple[float, float]]) -> list[dict[str, float]]:
    if len(points) < 5:
        return []
    slopes: list[tuple[float, float]] = []
    for left, right in zip(points, points[1:]):
        if right[0] <= left[0]:
            continue
        slopes.append((right[0], (right[1] - left[1]) / (right[0] - left[0])))
    negative = [index for index, (_, slope) in enumerate(slopes) if slope < -SLOPE_MARGIN]
    runs: list[dict[str, float]] = []
    cursor = 0
    while cursor < len(negative):
        end = cursor
        while end + 1 < len(negative) and negative[end + 1] == negative[end] + 1:
            end += 1
        first, last = negative[cursor], negative[end]
        before = slopes[first - 1][1] if first > 0 else -math.inf
        after = slopes[last + 1][1] if last + 1 < len(slopes) else -math.inf
        if end - cursor + 1 >= MIN_NEGATIVE_SECANTS and before > SLOPE_MARGIN and after > SLOPE_MARGIN:
            left_rho = points[first][0]
            right_rho = points[last + 1][0]
            runs.append({"rho_low": left_rho, "rho_high": right_rho, "negative_secants": end - cursor + 1, "slope_before": before, "slope_after": after})
        cursor = end + 1
    return runs


def _nearest(values: list[float], target: float) -> float | None:
    return min(values, key=lambda value: abs(value - target)) if values else None


def _oracle_row(slice_rows: list[dict[str, str]], xi: float, temperature: float) -> dict[str, str] | None:
    for row in slice_rows:
        if row.get("method") == "independent_oracle" and math.isclose(_float(row.get("xi")), xi, abs_tol=1e-9, rel_tol=0.0) and math.isclose(_float(row.get("T_MeV")), temperature, abs_tol=1e-5, rel_tol=0.0):
            return row
    return None


def _geometry_gate(row: dict[str, str] | None) -> bool:
    if row is None or not _bool(row.get("geometry_converged")):
        return False
    if not all(_finite(row.get(field)) for field in ("position_error_MeV", "density_error", "maxwell_area_gate", "mu_transition_MeV", "rho_hadron", "rho_quark", "mu_spinodal_hadron_MeV", "mu_spinodal_quark_MeV", "rho_spinodal_hadron", "rho_spinodal_quark", "area_residual")):
        return False
    return (
        _float(row.get("position_error_MeV")) <= MU_TOL
        and _float(row.get("density_error")) <= DENSITY_TOL
        and _float(row.get("maxwell_area_gate")) <= AREA_TOL
        and _float(row.get("area_residual")) <= AREA_TOL
        and abs(_float(row.get("rho_hadron")) - _float(row.get("rho_quark"))) > DENSITY_TOL
    )


def _feature_targets(global_points: list[tuple[float, float]], local_points: list[tuple[float, float]], oracle_row: dict[str, str] | None) -> list[float]:
    targets: list[float] = []
    for field in ("rho_hadron", "rho_quark", "rho_spinodal_hadron", "rho_spinodal_quark"):
        value = _float(oracle_row.get(field)) if oracle_row else math.nan
        if math.isfinite(value):
            targets.append(value)
    for candidate in _candidate_runs(global_points):
        targets.extend((candidate["rho_low"], (candidate["rho_low"] + candidate["rho_high"]) / 2.0, candidate["rho_high"]))
    if not targets and global_points:
        slopes = []
        for left, right in zip(global_points, global_points[1:]):
            if right[0] <= left[0]:
                continue
            slopes.append((abs((right[1] - left[1]) / (right[0] - left[0])), (left[0] + right[0]) / 2.0))
        targets.extend(value for _, value in sorted(slopes, reverse=True)[:8])
    if not targets and local_points:
        targets.extend((local_points[0][0], local_points[-1][0]))
    return sorted(set(round(value, 12) for value in targets if 0.0 <= value <= 4.0))


def _select_local_points(local_points: list[tuple[float, float]], targets: list[float], cap: int) -> list[tuple[float, float]]:
    if cap <= 0 or not local_points:
        return []
    selected: dict[float, tuple[float, float]] = {}
    for target in targets:
        nearest = _nearest([rho for rho, _ in local_points if rho not in selected], target)
        if nearest is not None:
            selected[nearest] = next(point for point in local_points if point[0] == nearest)
    remaining = [point for point in local_points if point[0] not in selected]
    # Deterministic midpoint priority: alternate from low/high ends so a
    # support interval is not accidentally filled from one side only.
    while len(selected) < cap and remaining:
        index = (len(selected) // 2) if len(selected) % 2 == 0 else -(len(selected) // 2 + 1)
        point = remaining.pop(index)
        selected[point[0]] = point
    return sorted(selected.values())


def _merge_candidates(*candidate_sets: list[dict[str, float]]) -> list[dict[str, float]]:
    flattened = sorted((candidate for candidates in candidate_sets for candidate in candidates), key=lambda item: item["rho_low"])
    merged: list[dict[str, float]] = []
    for candidate in flattened:
        # Only overlapping intervals are the same candidate seen at two
        # resolutions. A small gap is not a license to collapse two distinct
        # Maxwell candidates; those must remain ambiguous.
        if merged and candidate["rho_low"] <= merged[-1]["rho_high"] and candidate["rho_high"] >= merged[-1]["rho_low"]:
            merged[-1]["rho_high"] = max(merged[-1]["rho_high"], candidate["rho_high"])
            merged[-1]["negative_secants"] = max(merged[-1]["negative_secants"], candidate["negative_secants"])
            merged[-1]["slope_before"] = max(merged[-1]["slope_before"], candidate["slope_before"])
            merged[-1]["slope_after"] = max(merged[-1]["slope_after"], candidate["slope_after"])
        else:
            merged.append(dict(candidate))
    return merged


def _classify(
    points: list[tuple[float, float]],
    oracle_row: dict[str, str] | None,
    baseline_points: list[tuple[float, float]] | None = None,
) -> tuple[str, str, list[dict[str, float]]]:
    candidates = _merge_candidates(_candidate_runs(points), _candidate_runs(baseline_points or []))
    geometry = _geometry_gate(oracle_row)
    oracle_status = oracle_row.get("result_status") if oracle_row else None
    if len(candidates) > 1:
        return "ambiguous_near_critical", "multiple_maxwell_candidates", candidates
    if len(candidates) == 1 and geometry:
        return "confirmed_first_order", "unique_candidate_geometry_gate", candidates
    if len(candidates) == 0 and oracle_status == "confirmed_monotone":
        return "confirmed_monotone", "stable_no_s_shape", candidates
    if len(candidates) == 0:
        return "ambiguous_near_critical", "weak_or_unresolved_s_shape", candidates
    return "ambiguous_near_critical", "maxwell_or_geometry_unresolved", candidates


def _frontier(curves: list[dict[str, str]], slices: list[dict[str, str]], costs: list[dict[str, str]], caps: tuple[int, ...] = CAPS) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    anchors = sorted({(float(row["xi"]), float(row["T_MeV"])) for row in slices if row.get("method") == "independent_oracle"})
    dense_unique = sum(_float(row.get("unique_solves"), 0.0) for row in costs if row.get("method") == "memoized_dense")
    global_cache = {(xi, T): _group_points(curves, method="independent_oracle", xi=xi, temperature=T, levels={"0"}) for xi, T in anchors}
    oracle_fine_cache = {(xi, T): _group_points(curves, method="independent_oracle", xi=xi, temperature=T, levels={"0", "1"}) for xi, T in anchors}
    local_cache = {(xi, T): _group_points(curves, method="production_hybrid", xi=xi, temperature=T, levels={"4"}) for xi, T in anchors}
    base_cache = {(xi, T): _group_points(curves, method="production_hybrid", xi=xi, temperature=T, levels={"0", "1"}) for xi, T in anchors}

    frontier: list[dict[str, Any]] = []
    replay_rows: list[dict[str, Any]] = []
    candidates_rows: list[dict[str, Any]] = []
    representative_rows: list[dict[str, Any]] = []
    for cap in caps:
        total_unique = 0
        total_selected = 0
        matches = 0
        geometry_pass = 0
        multiple = 0
        mismatches = 0
        for xi, T in anchors:
            oracle_row = _oracle_row(slices, xi, T)
            global_points = global_cache[(xi, T)]
            oracle_fine_points = oracle_fine_cache[(xi, T)]
            local_points = local_cache[(xi, T)]
            targets = _feature_targets(global_points, local_points, oracle_row)
            selected = _select_local_points(local_points, targets, cap)
            # The complete Stage-B grid is the semantic baseline.  A local
            # Stage-C point at an already-covered rho is only a verification
            # sample and must not overwrite the global value; doing so can
            # manufacture a slope discontinuity from two solver traces.
            merged_map = {rho: mu for rho, mu in global_points}
            for rho, mu in selected:
                merged_map.setdefault(rho, mu)
            merged = sorted(merged_map.items())
            simulated_status, reason, candidates = _classify(merged, oracle_row, global_points)
            global_candidates = _candidate_runs(global_points)
            fine_candidates = _candidate_runs(oracle_fine_points)
            cross_resolution_multiple = max(len(global_candidates), len(fine_candidates)) > 1
            oracle_status = oracle_row.get("result_status", "missing") if oracle_row else "missing"
            total_unique += len(base_cache[(xi, T)]) + len(selected)
            total_selected += len(selected)
            matches += int(simulated_status == oracle_status)
            geometry_pass += int(simulated_status != "confirmed_first_order" or _geometry_gate(oracle_row))
            multiple += int(cross_resolution_multiple)
            mismatches += int(simulated_status != oracle_status)
            replay_rows.append({
                "xi": xi, "T_MeV": T, "cap": cap, "oracle_status": oracle_status,
                "simulated_status": simulated_status, "reason": reason,
                "candidate_count": len(candidates), "selected_targeted_points": len(selected),
                "base_unique_points": len(base_cache[(xi, T)]), "global_reference_points": len(global_points),
                "local_pool_points": len(local_points), "geometry_gate": _geometry_gate(oracle_row),
                "global_candidate_count": len(global_candidates), "fine_candidate_count": len(fine_candidates),
                "cross_resolution_multiple_candidates": cross_resolution_multiple,
            })
            if cap == caps[0] and (simulated_status != oracle_status or reason in {"unique_candidate_geometry_gate", "stable_no_s_shape"}):
                stride = max(1, len(global_points) // 64)
                for index, (rho, mu) in enumerate(global_points):
                    if index % stride == 0 or index == len(global_points) - 1:
                        representative_rows.append({"xi": xi, "T_MeV": T, "rho": rho, "muq_MeV": mu, "source": "stage_b_global_reference", "cap": cap, "oracle_status": oracle_status, "simulated_status": simulated_status, "reason": reason})
                for rho, mu in selected:
                    representative_rows.append({"xi": xi, "T_MeV": T, "rho": rho, "muq_MeV": mu, "source": "stage_c_local_selected", "cap": cap, "oracle_status": oracle_status, "simulated_status": simulated_status, "reason": reason})
            for index, candidate in enumerate(candidates, start=1):
                candidates_rows.append({"xi": xi, "T_MeV": T, "cap": cap, "candidate_index": index, **candidate, "geometry_gate": _geometry_gate(oracle_row), "oracle_status": oracle_status, "simulated_status": simulated_status, "reason": reason})
        frontier.append({
            "cap": cap, "anchors": len(anchors), "state_matches": matches,
            "state_match_fraction": matches / max(len(anchors), 1), "classification_mismatches": mismatches,
            "geometry_gate_pass": geometry_pass == len(anchors), "multiple_candidate_anchors": multiple,
            "unique_solves": total_unique, "targeted_points": total_selected,
            "dense_unique_solves": dense_unique, "cost_gate": total_unique <= dense_unique,
            "state_gate": matches == len(anchors), "maxwell_candidate_gate": multiple == 0,
        })

    valid = [row for row in frontier if row["state_gate"] and row["geometry_gate_pass"] and row["maxwell_candidate_gate"] and row["cost_gate"]]
    selected = min(valid, key=lambda row: (row["unique_solves"], row["cap"])) if valid else None
    if selected is None:
        if any(row["multiple_candidate_anchors"] for row in frontier):
            verdict = "maxwell_candidate_inconclusive"
        elif any(not row["cost_gate"] for row in frontier):
            verdict = "performance_inconclusive"
        elif any(row["classification_mismatches"] for row in frontier):
            verdict = "maxwell_candidate_inconclusive"
        else:
            verdict = "integration_failed"
    else:
        verdict = "feasible_candidate"
    policy = {
        "schema_version": "cep_hybrid_stagec_offline_feasibility_v1",
        "verdict": verdict,
        "selected_cap": selected["cap"] if selected else None,
        "selected_policy": "global_stage_b_plus_feature_midpoint_stage_c" if selected else None,
        "stage_b_global_reference": "independent_oracle.rho_level=0 (offline semantic reference only)",
        "production_cost_base": "production_hybrid rho_level 0/1 points",
        "caps": list(caps), "dense_unique_solves": dense_unique,
        "reason": "no cap satisfies all gates" if selected is None else "minimum-cost cap satisfying all gates",
    }
    return frontier, replay_rows, {"policy": policy, "candidates": candidates_rows, "representative": representative_rows}


def _write_docs(output_dir: Path, policy: dict[str, Any], frontier: list[dict[str, Any]], replay_rows: list[dict[str, Any]], input_manifest: dict[str, Any]) -> None:
    verdict = policy["verdict"]
    mismatches = sum(int(row["classification_mismatches"]) for row in frontier)
    (output_dir / "README.md").write_text(
        f"# PNJL Hybrid Stage-C offline feasibility v1\n\n"
        f"verdict: `{verdict}`。本目录只使用 aggregate replay CSV 做离线重放，不调用 equilibrium solver，"
        "不改变 production、reference 或历史 evidence。\n\n"
        f"- source run: `{input_manifest.get('source_run_id', '')}`\n"
        f"- source calculation SHA: `{input_manifest.get('expected_calculation_sha', '')}`\n"
        f"- tested caps: `{', '.join(str(row['cap']) for row in frontier)}`\n"
        f"- classification mismatches across cap evaluations: `{mismatches}`\n"
        f"- selected cap: `{policy.get('selected_cap')}`\n\n"
        "Stage-C classification uses the complete oracle-level-0 0.00625 curve as an offline semantic reference,"
        " merges selected local 0.003125 points, enumerates all stable +→−→+ candidates, and uses the existing"
        " oracle geometry row only as a non-solver geometry gate. The oracle curve is not charged to simulated"
        " hybrid cost.\n\n"
        "若 verdict 不是 `feasible_candidate`，不得创建 production PR；当前结果只用于定位弱 S、Maxwell"
        " 候选或成本问题。\n",
        encoding="utf-8",
    )
    (output_dir / "AUDIT.md").write_text(
        "# Stage-C offline feasibility audit\n\n"
        "输入为 final aggregate replay；所有原始 rho–mu 曲线仍由 Actions/local artifact 追溯。"
        "本审计不将离线分类提升为物理结论，不写 reference，也不运行 PNJL。\n",
        encoding="utf-8",
    )
    _write_csv(output_dir / "claim_ledger.csv", [
        {"claim_id": "classification", "claim": "offline hybrid matches corrected oracle at all anchors", "status": "pass" if policy["verdict"] == "feasible_candidate" else "inconclusive", "boundary": "diagnostic only"},
        {"claim_id": "maxwell_candidates", "claim": "cross-resolution Maxwell candidate is unique", "status": "pass" if policy["verdict"] == "feasible_candidate" else "inconclusive", "boundary": "author curve review required"},
        {"claim_id": "cost", "claim": "simulated hybrid unique solves do not exceed dense", "status": "pass" if all(row["cost_gate"] for row in frontier) else "inconclusive", "boundary": "offline estimate only"},
    ])
    (output_dir / "plot_manifest.json").write_text(json.dumps({"schema_version": "cep_hybrid_stagec_offline_feasibility_v1", "figures": [], "source_curve_artifact": "external_actions_or_local_artifact"}, indent=2) + "\n", encoding="utf-8")


def run(input_dir: Path, output_dir: Path, expected_sha: str) -> dict[str, Any]:
    errors = _validate_input(input_dir, expected_sha)
    output_dir.mkdir(parents=True, exist_ok=True)
    if errors:
        policy = {"schema_version": "cep_hybrid_stagec_offline_feasibility_v1", "verdict": "integration_failed", "input_errors": errors}
        (output_dir / "selected_policy.json").write_text(json.dumps(policy, indent=2) + "\n", encoding="utf-8")
        return policy

    manifest = json.loads((input_dir / "manifest.json").read_text(encoding="utf-8"))
    curves, slices, costs = _read_csv(input_dir / "curve_points.csv"), _read_csv(input_dir / "slice_metrics.csv"), _read_csv(input_dir / "method_costs.csv")
    frontier, replay_rows, details = _frontier(curves, slices, costs)
    policy = details["policy"]
    _write_csv(output_dir / "stage_c_cost_frontier.csv", frontier)
    _write_csv(output_dir / "anchor_replay.csv", replay_rows)
    _write_csv(output_dir / "maxwell_candidates.csv", details["candidates"])
    _write_csv(output_dir / "curve_index.csv", [{"source": "aggregate_replay", "path": "curve_points.csv", "sha256": _sha256(input_dir / "curve_points.csv"), "raw_curve_copy_in_repository": False}])
    _write_csv(output_dir / "representative_curves.csv", details["representative"])
    (output_dir / "selected_policy.json").write_text(json.dumps(policy, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write_docs(output_dir, policy, frontier, replay_rows, manifest)
    hashes: dict[str, str] = {}
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and path.name != "manifest.json":
            hashes[str(path.relative_to(output_dir)).replace("\\", "/")] = _sha256(path)
    output_manifest = {
        "schema_version": "cep_hybrid_stagec_offline_feasibility_v1",
        "verdict": policy["verdict"], "source_run_id": manifest.get("source_run_id", ""),
        "source_calculation_sha": manifest.get("expected_calculation_sha", ""),
        "source_replay_manifest_sha256": _sha256(input_dir / "manifest.json"),
        "selected_policy": policy, "file_sha256": hashes,
        "generated_at_utc": datetime.now(timezone.utc).isoformat().replace("+00:00", "Z"),
        "solver_called": False,
    }
    (output_dir / "manifest.json").write_text(json.dumps(output_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return policy


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--expected-calculation-sha", required=True)
    args = parser.parse_args()
    policy = run(args.input_dir, args.output_dir, args.expected_calculation_sha)
    print(json.dumps(policy, indent=2, sort_keys=True))
    return 0 if policy.get("verdict") == "feasible_candidate" else 2 if policy.get("verdict") in {"maxwell_candidate_inconclusive", "performance_inconclusive"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
