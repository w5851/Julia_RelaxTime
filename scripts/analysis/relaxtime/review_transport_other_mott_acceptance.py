#!/usr/bin/env python3
"""Audit completed other-Mott grids and explicitly track missing reference work."""
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import argparse
import importlib.util
import json
import math
import tomllib


def check_grid(observations, plan, require_all):
    by_id = {row["point_id"]: row for row in plan}
    actual = [row["point_id"] for row in observations]
    if len(actual) != len(set(actual)) or not set(actual) <= set(by_id):
        raise ValueError("Duplicate or unplanned observations")
    if require_all and set(actual) != set(by_id):
        raise ValueError("Incomplete grid")
    for row in observations:
        original = by_id[row["point_id"]]
        for field in ("curve_id", "window", "meson"):
            if row[field] != original[field]:
                raise ValueError("Observation identity mismatch: " + field)
        for field in ("xi", "xi_mott", "T_MeV", "muB_MeV"):
            if not math.isfinite(float(row[field])) or abs(float(row[field])-float(original[field])) > 1e-11:
                raise ValueError("Observation coordinate mismatch: " + field)
        if not math.isfinite(float(row["background_residual"])) or float(row["background_residual"]) > 1e-10:
            raise ValueError("Background residual exceeds diagnostic contract")
    return sorted(set(by_id)-set(actual))


def xi_comparisons(edges, references, observables, slope):
    results, missing = [], []
    for edge in edges:
        ids = [edge[name] for name in ("left_id", "midpoint_id", "right_id")]
        if not all(point in references for point in ids):
            missing.append(dict(edge, missing_ids=[point for point in ids if point not in references]))
            continue
        a, m, b = (references[point] for point in ids)
        xa, xm, xb = (float(row["xi"]) for row in (a, m, b))
        if not xa < xm < xb or any(row["curve_id"] != edge["curve_id"] for row in (a, m, b)):
            raise ValueError("Invalid xi refinement geometry")
        for observable in observables:
            coarse = slope(float(a[observable]), float(b[observable]), xb-xa)
            left = slope(float(a[observable]), float(m[observable]), xm-xa)
            right = slope(float(m[observable]), float(b[observable]), xb-xm)
            results.append(dict(curve_id=edge["curve_id"], observable=observable,
                left_id=ids[0], midpoint_id=ids[1], right_id=ids[2], xi_left=xa, xi_midpoint=xm, xi_right=xb,
                coarse_log_slope=coarse, refined_left_log_slope=left, refined_right_log_slope=right,
                max_split_slope_change=max(abs(left-coarse), abs(right-coarse))))
    return results, missing


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--source", type=Path, default=Path("D:/CodexCompute/trfull-20261026-v1/case/source"))
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()
    case = args.case.resolve()
    spec = importlib.util.spec_from_file_location("frozen_slopes", case/"review_transport_other_mott_slopes_v3.py")
    review = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(review)
    read = lambda path: json.loads(path.read_text(encoding="utf-8-sig"))
    source = read(case/"expected_source_hashes.json")
    assert len(source) == 340 and all(review.sha(args.source/name) == digest for name, digest in source.items())
    fine_runs = list(read(case/"fine_slopes_v1/summary.json")["provenance"])
    reference_runs = []
    for path in sorted(case.glob("tighter_*.status.json")):
        status = read(path)
        if status["status"] == "completed_diagnostic":
            reference_runs.append(path.name.removesuffix(".status.json"))
    for name in fine_runs+reference_runs:
        status = read(case/(name+".status.json"))
        assert status["status"] == "completed_diagnostic" and status["source_all_match"]
        assert status["source_files_checked_after"] == 340
    obs, rates, weights, channels, provenance = review.collect(case, fine_runs)
    assert len(obs) == 195
    fine_plan = review.read_csv(case/"locate_v1/transport_plan.csv")
    check_grid(obs, fine_plan, True)
    targets = {row["curve_id"]: row for row in read(case/"inputs/targets.json")["targets"]}
    assert all(targets[row["curve_id"]]["segment_left"] <= float(row["xi"]) <=
               targets[row["curve_id"]]["segment_right"] for row in obs)
    products, closure = review.analyze(obs, rates, weights, channels)
    observations = {row["point_id"]: row for row in obs}
    controls = read(case/"convergence_plan_v1/summary.json")
    expected_reference_ids = {point for job in controls["jobs"] for point in job["point_ids"]}
    whole_plan = review.read_csv(case/"convergence_plan_v1/transport_plan.csv")
    reference_plan = [row for row in whole_plan if row["point_id"] in expected_reference_ids]
    assert len(reference_plan) == len(expected_reference_ids) == controls["total_evaluations"]
    point_checks, reference_obs, slope_checks, xi_checks = [], [], [], []
    references = {}
    if reference_runs:
        reference_data = review.collect(case, reference_runs)
        reference_obs = reference_data[0]
        check_grid(reference_obs, reference_plan, args.require_complete)
        ref_products, ref_closure = review.analyze(*reference_data[:4])
        products.update({"reference_"+key: value for key, value in ref_products.items()})
        closure = max(closure, ref_closure)
        provenance.update(reference_data[4])
        references = {row["point_id"]: row for row in reference_obs}
        for row in reference_obs:
            if row["point_id"] not in observations:
                continue
            original = observations[row["point_id"]]
            for observable in review.OBSERVABLES:
                a, b = float(original[observable]), float(row[observable])
                point_checks.append(dict(curve_id=row["curve_id"], point_id=row["point_id"], observable=observable,
                    fine_value=a, tighter_value=b, relative_difference=(b-a)/b))
        slope_checks = review.compare_native_slopes(products["slopes"], references)
    missing_ids = sorted(expected_reference_ids-set(references))
    xi_checks, missing_xi = xi_comparisons(controls["refinements"], references, review.OBSERVABLES, review.log_slope)
    if args.require_complete and (missing_ids or missing_xi):
        raise ValueError("Reference acceptance requires every planned point and split interval")
    # Every profile must retain the tighter accuracy and all physical quadratures.
    profiles = {}
    standard = tomllib.loads((case/"inputs/tighter.toml").read_text(encoding="utf-8"))
    for variant in {row["variant"] for row in reference_obs}:
        profile_path = case/"inputs"/(variant+".toml")
        actual = tomllib.loads(profile_path.read_text(encoding="utf-8"))
        for key in ("algorithm", "background", "a_integral", "tau_integral", "terminal_integral"):
            assert actual[key] == standard[key]
        for key in ("rtol", "atol", "weighted_window", "seed_order", "max_rounds", "max_events"):
            assert actual["sampling"][key] == standard["sampling"][key]
        profiles[variant] = dict(sha256=review.sha(profile_path), sampling=actual["sampling"])
    interval_groups = defaultdict(list)
    for row in products["slopes"]:
        if abs(row["delta_xi_midpoint"]) <= .004:
            interval_groups[row["curve_id"], row["observable"]].append(row)
    metrics = []
    for (curve, observable), rows in sorted(interval_groups.items()):
        lower, upper = min(rows, key=lambda r: r["log_slope"]), max(rows, key=lambda r: r["log_slope"])
        precision = [row for row in slope_checks if row["curve_id"] == curve and row["observable"] == observable]
        splits = [row for row in xi_checks if row["curve_id"] == curve and row["observable"] == observable]
        span = upper["log_slope"]-lower["log_slope"]
        delta = max((row["absolute_difference"] for row in precision), default=None)
        metrics.append(dict(curve_id=curve, meson=upper["meson"], window=upper["window"], observable=observable,
            xi_mott=upper["xi_mott"], min_log_slope=lower["log_slope"], max_log_slope=upper["log_slope"],
            slope_range=span, min_left=lower["point_left"], min_right=lower["point_right"],
            max_left=upper["point_left"], max_right=upper["point_right"],
            min_delta_xi=lower["delta_xi_midpoint"], max_delta_xi=upper["delta_xi_midpoint"],
            matched_precision_intervals=len(precision), maximum_precision_slope_difference=delta,
            precision_difference_over_range=delta/span if delta is not None and span else None,
            split_intervals=len(splits), max_observed_split_change=max((row["max_split_slope_change"] for row in splits), default=None)))
    group_index = {(row["curve_id"], row["point_left"], row["point_right"], row["species"], row["process_group"]): row
                   for row in products["process_group_derivatives"]}
    modulation = []
    for row in metrics:
        if not row["observable"].startswith("tau_"):
            continue
        species = row["observable"].removeprefix("tau_")
        for group in ("K_s_processes", "other_processes"):
            a = group_index[row["curve_id"], row["min_left"], row["min_right"], species, group]
            b = group_index[row["curve_id"], row["max_left"], row["max_right"], species, group]
            modulation.append(dict(curve_id=row["curve_id"], species=species, process_group=group,
                total_slope_range=row["slope_range"], min_left=row["min_left"], min_right=row["min_right"],
                max_left=row["max_left"], max_right=row["max_right"],
                delta_slope_from_rate=b["log_slope_from_rate"]-a["log_slope_from_rate"],
                delta_slope_from_density=b["log_slope_from_density"]-a["log_slope_from_density"],
                delta_slope_contribution=b["log_slope_contribution"]-a["log_slope_contribution"]))
    neighbours = []
    for curve, target in targets.items():
        center = float(next(row for row in obs if row["curve_id"] == curve)["xi_mott"])
        for crossing in read(case/"survey_review_v2/summary.json")["root_crossings"]:
            if crossing["window"] != target["window"] or crossing["meson"] == target["meson"]:
                continue
            a, b = crossing["xi_left"], crossing["xi_right"]
            distance = max(a-center, center-b, 0.)
            if distance <= .025:
                neighbours.append(dict(curve_id=curve, target_meson=target["meson"], xi_mott=center,
                    neighbour_meson=crossing["meson"], cut=crossing["cut"], status=crossing["status"],
                    xi_left=a, xi_right=b, distance_to_bracket=distance,
                    overlaps_plus_minus_0008=(a <= center+.008 and b >= center-.008)))
    args.output.mkdir(parents=True, exist_ok=False)
    for key, rows in dict(products, point_convergence=point_checks, slope_convergence=slope_checks,
        xi_resolution=xi_checks, window_metrics=metrics, slope_modulation_by_process=modulation, nearby_mott_brackets=neighbours).items():
        review.write_csv(args.output/(key+".csv"), rows)
    summary = dict(status="completed_diagnostic_review", whole_task_complete=False,
        acceptance_scope="main grid and completed reference runs only; no meson-exclusivity or formal-promotion verdict",
        main_grid_complete=True, expected_main_points=195, points=len(obs), curves=len(targets), channels=len(rates),
        source_files_checked=340, source_all_match=True, exact_decomposition_max_absolute_closure_error=closure,
        expected_reference_evaluations=len(expected_reference_ids), reference_evaluations=len(reference_obs),
        reference_plan_complete=not missing_ids, missing_reference_ids=missing_ids,
        comparison_points=len({row["point_id"] for row in point_checks}),
        new_xi_evaluations=len(set(references)-set(observations)),
        matched_slope_comparisons=len(slope_checks), xi_resolution_checks=len(xi_checks), missing_xi_intervals=missing_xi,
        max_point_relative_difference=max((abs(row["relative_difference"]) for row in point_checks), default=None),
        max_slope_precision_difference=max((row["absolute_difference"] for row in slope_checks), default=None),
        fine_runs=fine_runs, reference_runs=reference_runs, provenance=provenance, reference_profiles=profiles,
        contrasts=products["contrasts"], timestamp_utc=datetime.now(timezone.utc).isoformat(),
        script_sha256=review.sha(Path(__file__)), slope_reviewer_sha256=review.sha(case/"review_transport_other_mott_slopes_v3.py"))
    (args.output/"summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps({key: summary[key] for key in ("main_grid_complete", "points", "reference_evaluations",
        "expected_reference_evaluations", "reference_plan_complete", "matched_slope_comparisons", "xi_resolution_checks",
        "max_point_relative_difference", "max_slope_precision_difference")}, allow_nan=False))


if __name__ == "__main__":
    main()
