#!/usr/bin/env python3
"""Select confirmed closed-threshold events and two pi controls on quote-remote."""
import csv
import hashlib
import json
from pathlib import Path
import argparse


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    args = parser.parse_args()
    source_path = args.case / "inputs/source_rows.csv"
    report_path = args.case / "survey_review_v2/summary.json"
    report = json.loads(report_path.read_text(encoding="utf-8-sig"))
    with source_path.open(encoding="utf-8-sig", newline="") as stream:
        source = list(csv.DictReader(stream))
    by_id = {r["point_id"]: r for r in source}
    targets = []
    for event in report["root_crossings"]:
        if event["status"] != "coarse_confirmed_mott":
            continue
        control = event["meson"] == "pi" and event["window"] in (
            "mode_a:muB900.0_alpha1.0", "mode_a:muB900.0_alpha1.1")
        if event["meson"] in ("pi", "eta") and not control:
            continue
        segment = by_id[event["point_left"]]["continuation_segment"]
        xs = [float(r["xi"]) for r in source if r["continuation_segment"] == segment]
        targets.append(dict(curve_id=f"other{len(targets)+1:02d}", role="positive_control" if control else "primary",
            segment_left=min(xs), segment_right=max(xs),
            search_left=max(min(xs), event["xi_left"]-0.02),
            search_right=min(max(xs), event["xi_right"]+0.02), **event))
    assert len(targets) == 13
    assert sum(t["meson"] == "K" for t in targets) == 10
    assert sum(t["meson"] == "sigma" for t in targets) == 1
    config = dict(schema="transport_other_mott_targets_v1", targets=targets,
        xi_offsets=[-.008, -.004, -.002, -.001, -.0005, -.0002, 0,
                     .0002, .0005, .001, .002, .004, .008],
        source_rows_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
        survey_review_sha256=hashlib.sha256(report_path.read_bytes()).hexdigest(),
        endpoint_bisection_tolerance=1e-9,
        notes="Other high-cut crossings and unresolved branches are retained in the survey; they are not closed-threshold Mott events.")
    path = args.case / "inputs/targets.json"
    if path.exists():
        raise FileExistsError(path)
    path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(config))


if __name__ == "__main__":
    main()
