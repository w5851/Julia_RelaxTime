#!/usr/bin/env python3
"""Compare symmetric Mott-window slope modulation on quote-remote.

Reuse the frozen slope reviewer's exact decomposition on three actual nodes.
Inserted reference nodes must not change the endpoints used for attribution.
This is a derived diagnostic; it performs no new physics solve or interpolation.
"""
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import argparse
import importlib.util
import json
import math
import shutil


HALF_WIDTHS = (.004, .002, .001, .0005, .0002, .0001)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--review", default="acceptance_review_complete_v1")
    args = parser.parse_args()
    case = args.case.resolve()
    review_path = case / args.review
    summary_path = review_path / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8-sig"))
    if not summary["main_grid_complete"] or not summary["reference_plan_complete"]:
        raise ValueError("Complete main and reference grids are required")
    reviewer_path = case / "review_transport_other_mott_slopes_v3.py"
    spec = importlib.util.spec_from_file_location("frozen_slopes", reviewer_path)
    review = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(review)
    if review.sha(reviewer_path) != summary["slope_reviewer_sha256"]:
        raise ValueError("Frozen reviewer identity mismatch")

    centers, process_rows, final_rows, unavailable = [], [], [], []
    input_hashes = {str(summary_path.relative_to(case)): review.sha(summary_path),
                    str(reviewer_path.relative_to(case)): review.sha(reviewer_path)}
    max_closure = max_frozen_difference = 0.0
    frozen_checks = 0
    for profile, run_key, frozen_name in (
        ("fine", "fine_runs", "centered_slopes.csv"),
        ("tighter", "reference_runs", "reference_centered_slopes.csv"),
    ):
        data = review.collect(case, summary[run_key])
        if any(digest != summary["provenance"][name] for name, digest in data[4].items()):
            raise ValueError("Run provenance differs from complete review")
        input_hashes.update({name + "/summary.json": digest for name, digest in data[4].items()})
        frozen_path = review_path / frozen_name
        input_hashes[str(frozen_path.relative_to(case))] = review.sha(frozen_path)
        frozen = {(r["curve_id"], float(r["half_width"]), r["observable"]): r
                  for r in review.read_csv(frozen_path)}
        by_curve = defaultdict(list)
        for row in data[0]:
            by_curve[row["curve_id"]].append(row)
        for curve, rows in sorted(by_curve.items()):
            xi_mott = float(rows[0]["xi_mott"])
            for h in HALF_WIDTHS:
                selected = []
                for x in (xi_mott-h, xi_mott, xi_mott+h):
                    matches = [r for r in rows if abs(float(r["xi"])-x) < 1e-12]
                    if len(matches) > 1:
                        raise ValueError("Ambiguous center-window coordinate")
                    selected.extend(matches)
                if len(selected) != 3:
                    unavailable.append(dict(profile=profile, curve_id=curve, half_width=h))
                    continue
                ids = [r["point_id"] for r in selected]
                subsets = [[r for r in table if r["point_id"] in ids] for table in data[1:4]]
                products, closure = review.analyze(selected, *subsets)
                max_closure = max(max_closure, closure)
                info = dict(profile=profile, curve_id=curve, meson=rows[0]["meson"],
                    window=rows[0]["window"], T_MeV=float(rows[0]["T_MeV"]),
                    muB_MeV=float(rows[0]["muB_MeV"]), xi_mott=xi_mott, half_width=h,
                    point_left=ids[0], point_center=ids[1], point_right=ids[2])
                center_by_observable = {}
                for observable in review.OBSERVABLES:
                    pair = [r for r in products["slopes"] if r["observable"] == observable]
                    left, right = pair
                    if (left["point_left"], left["point_right"], right["point_left"], right["point_right"]) != (ids[0], ids[1], ids[1], ids[2]):
                        raise ValueError("Attribution endpoints are not the selected three nodes")
                    delta = right["log_slope"] - left["log_slope"]
                    row = dict(info, observable=observable, left_log_slope=left["log_slope"],
                        right_log_slope=right["log_slope"], right_minus_left=delta,
                        estimated_s_interpolation_modulation_error=sum(
                            r["estimated_s_interpolation_slope_error"] for r in pair))
                    expected = frozen[curve, h, observable]
                    for key in ("left_log_slope", "right_log_slope", "right_minus_left"):
                        difference = abs(row[key]-float(expected[key]))
                        max_frozen_difference = max(max_frozen_difference, difference)
                        if not math.isclose(row[key], float(expected[key]), abs_tol=1e-8, rel_tol=1e-10):
                            raise ValueError("Symmetric slopes differ from frozen centered table")
                    frozen_checks += 1
                    centers.append(row)
                    center_by_observable[observable] = delta
                for product, target, keys, fields in (
                    ("process_group_derivatives", process_rows, ("species", "process_group"),
                     ("log_slope_from_rate", "log_slope_from_density", "log_slope_contribution")),
                    ("final_species_derivatives", final_rows, ("observable", "species"),
                     ("log_slope_from_tau", "log_slope_from_coefficient", "log_slope_contribution")),
                ):
                    pairs = defaultdict(list)
                    for row in products[product]:
                        pairs[tuple(row[k] for k in keys)].append(row)
                    sums = defaultdict(float)
                    for identity, (left, right) in pairs.items():
                        values = {"delta_"+field: right[field]-left[field] for field in fields}
                        target.append(dict(info, **dict(zip(keys, identity)), **values))
                        observable = "tau_"+identity[0] if product.startswith("process") else identity[0]
                        sums[observable] += values["delta_log_slope_contribution"]
                    for observable, total in sums.items():
                        difference = abs(total-center_by_observable[observable])
                        max_closure = max(max_closure, difference)
                        if difference >= 1e-7:
                            raise ValueError("Modulation decomposition does not close")

    lookup = {(r["profile"], r["curve_id"], r["half_width"], r["observable"]): r for r in centers}
    comparisons = []
    for row in centers:
        if row["profile"] != "fine":
            continue
        ref = lookup.get(("tighter", row["curve_id"], row["half_width"], row["observable"]))
        if ref is not None:
            comparisons.append(dict(curve_id=row["curve_id"], meson=row["meson"],
                half_width=row["half_width"], observable=row["observable"],
                fine_modulation=row["right_minus_left"], tighter_modulation=ref["right_minus_left"],
                absolute_difference=abs(row["right_minus_left"]-ref["right_minus_left"])))
    args.output.mkdir(parents=True, exist_ok=False)
    files = dict(centered_modulation=centers, process_modulation=process_rows,
                 final_species_modulation=final_rows, centered_precision=comparisons)
    for name, rows in files.items():
        review.write_csv(args.output/(name+".csv"), rows)
    (args.output/"source").mkdir()
    for path in (Path(__file__), reviewer_path):
        shutil.copyfile(path, args.output/"source"/path.name)
    report = dict(schema="other_mott_center_modulation_v1", status="completed_derived_diagnostic",
        timestamp_utc=datetime.now(timezone.utc).isoformat(), compute_host="quote-remote",
        input_review=str(summary_path), input_sha256=input_hashes,
        generator_sha256=review.sha(Path(__file__)),
        definition="J_y(h)=S_y[xi_M,xi_M+h]-S_y[xi_M-h,xi_M]; S_y=delta(log(y))/delta(xi)",
        units="dimensionless; xi is dimensionless", half_widths=list(HALF_WIDTHS),
        row_counts={name:len(rows) for name, rows in files.items()},
        frozen_centered_rows_verified=frozen_checks,
        max_difference_from_frozen_centered_table=max_frozen_difference,
        max_absolute_decomposition_closure_error=max_closure,
        unavailable_symmetric_triplets=unavailable,
        input_values_changed=False, interpolation=False, new_physics_solves=False,
        scope="Finite-window modulation, not the xi-continuous derivative limit. s estimates cover interpolation only. Full-process grouping is not an amplitude ablation.",
        output_sha256={name+".csv":review.sha(args.output/(name+".csv")) for name in files})
    (args.output/"summary.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps({k:report[k] for k in ("status", "row_counts", "frozen_centered_rows_verified",
        "max_difference_from_frozen_centered_table", "max_absolute_decomposition_closure_error")}))


if __name__ == "__main__":
    main()
