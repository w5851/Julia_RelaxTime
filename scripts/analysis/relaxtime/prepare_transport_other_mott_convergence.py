#!/usr/bin/env python3
"""Choose same-coordinate precision controls and resolved-peak xi midpoints."""
from collections import defaultdict
from pathlib import Path
import argparse
import csv
import hashlib
import json


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    original_plan = args.case / "locate_v1/transport_plan.csv"
    original = read_csv(original_plan)
    by_id = {r["point_id"]: r for r in original}
    grouped = defaultdict(list)
    for row in read_csv(args.case / "fine_slopes_v1/slopes.csv"):
        if abs(float(row["delta_xi_midpoint"])) <= .004:
            grouped[row["curve_id"], row["observable"]].append(row)
    chosen, edges = defaultdict(set), defaultdict(dict)
    for row in original:
        if row["kind"] == "mott_centered" and abs(float(row["delta_xi"])) in (0, .0002, .001):
            chosen[row["curve_id"]].add(row["point_id"])
    for (curve, _), rows in grouped.items():
        # Resolve both the maximum and minimum slope, including a negative response.
        for edge in (max(rows, key=lambda r: float(r["log_slope"])),
                     min(rows, key=lambda r: float(r["log_slope"]))):
            key = edge["point_left"], edge["point_right"]
            chosen[curve].update(key)
            edges[curve][key] = edge
    rows = [dict(row) for row in original]
    jobs, refinements = [], []
    for curve, ids in sorted(chosen.items()):
        extra = []
        for i, ((left, right), edge) in enumerate(sorted(edges[curve].items()), 1):
            a, b = by_id[left], by_id[right]
            x = (float(a["xi"])+float(b["xi"]))/2
            if any(r["curve_id"] == curve and abs(float(r["xi"])-x)<1e-12 for r in original):
                continue
            row = dict(a)
            pid = f"{curve}_x{i:02d}"
            row.update(point_id=pid, xi=x, delta_xi=x-float(a["xi_mott"]),
                       kind="xi_midpoint_refinement", original_point_id="")
            rows.append(row)
            extra.append(pid)
            refinements.append(dict(curve_id=curve, midpoint_id=pid, left_id=left, right_id=right,
                                    xi=x, original_width=float(b["xi"])-float(a["xi"])))
        jobs.append(dict(name="tighter_"+curve+"_v1", curve_id=curve, variant="tighter",
                         point_ids=sorted(ids)+extra, same_coordinate_points=len(ids), new_xi_points=len(extra)))
    args.output.mkdir(parents=True, exist_ok=False)
    with (args.output / "transport_plan.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(original[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary = dict(status="prepared_precision_and_xi_controls", jobs=jobs, refinements=refinements,
        original_plan_sha256=hashlib.sha256(original_plan.read_bytes()).hexdigest(),
        total_evaluations=sum(len(j["point_ids"]) for j in jobs),
        new_xi_points=sum(j["new_xi_points"] for j in jobs),
        selection="Center and +/-0.0002,+/-0.001, plus endpoints and midpoints of each observable's sampled max/min slope within +/-0.004.")
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
