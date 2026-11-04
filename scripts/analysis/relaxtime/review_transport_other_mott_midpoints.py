#!/usr/bin/env python3
"""Audit continuous threshold crossings after midpoint background solves."""
from collections import defaultdict
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    path = args.case / "midpoint_screen_v1/thresholds.csv"
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    original = json.loads((args.case / "survey_review_v2/summary.json").read_text(encoding="utf-8-sig"))
    known = [r for r in original["diagnostic_threshold_crossings"] if r["continuous_background"]
             and r["cut"] in ("pair", "uu")]
    groups = defaultdict(list)
    for row in rows:
        assert all(math.isfinite(float(row[k])) for k in ("inverse_re", "inverse_im", "inverse_abs")), row
        groups[row["continuation_segment"], row["meson"], row["cut"]].append(row)
    events, new_events, ranges = [], [], []
    max_spacing = 0
    for (segment, meson, cut), curve in sorted(groups.items()):
        curve.sort(key=lambda r: float(r["xi"]))
        closest = min(curve, key=lambda r: float(r["inverse_abs"]))
        ranges.append(dict(segment=segment, meson=meson, cut=cut, points=len(curve),
            min_inverse_abs=float(closest["inverse_abs"]), closest_xi=float(closest["xi"]),
            min_inverse_re=min(float(r["inverse_re"]) for r in curve),
            max_inverse_re=max(float(r["inverse_re"]) for r in curve)))
        for a, b in zip(curve, curve[1:]):
            xa, xb = float(a["xi"]), float(b["xi"])
            max_spacing = max(max_spacing, xb-xa)
            if float(a["inverse_re"])*float(b["inverse_re"]) >= 0:
                continue
            closed = cut in ("pair", "uu") and max(abs(float(r["inverse_im"])) for r in (a, b)) < 1e-7
            matched = any(k["window"] == a["window"] and k["meson"] == meson and k["cut"] == cut
                          and k["xi_left"]-.01 <= (xa+xb)/2 <= k["xi_right"]+.01 for k in known)
            item = dict(segment=segment, window=a["window"], meson=meson, cut=cut,
                xi_left=xa, xi_right=xb, closed_threshold_candidate=closed,
                matches_original_candidate=matched if closed else None,
                inverse_left=float(a["inverse_re"]), inverse_right=float(b["inverse_re"]),
                imaginary_left=float(a["inverse_im"]), imaginary_right=float(b["inverse_im"]))
            events.append(item)
            if closed and not matched:
                new_events.append(item)
    args.output.mkdir(parents=True, exist_ok=False)
    for name, values in (("crossings.csv", events), ("ranges.csv", ranges)):
        with (args.output/name).open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(values[0]))
            writer.writeheader()
            writer.writerows(values)
    summary = dict(status="completed_threshold_screen_review", threshold_rows=len(rows),
        max_continuous_xi_spacing=max_spacing,
        closed_threshold_crossings=sum(r["closed_threshold_candidate"] for r in events),
        new_closed_candidates=new_events, crossings=events,
        source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        limitation="Finite sample spacing; no claim excluding arbitrarily narrow or tangential unsampled zeros.")
    (args.output/"summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
