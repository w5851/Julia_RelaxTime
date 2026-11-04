#!/usr/bin/env python3
"""Freeze one completed location for transport while other locations finish."""
from pathlib import Path
import argparse
import csv
import hashlib
import json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("curve")
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    config = json.loads((args.case / "inputs/targets.json").read_text(encoding="utf-8-sig"))
    target = next(t for t in config["targets"] if t["curve_id"] == args.curve)
    raw = (args.case / "locate_v1/locations.json").read_bytes()
    location = next(t for t in json.loads(raw) if t["curve_id"] == args.curve)
    center = location["levels"]["base"]["xi"]
    assert all(abs(v["xi"]-center) < 1e-6 for v in location["levels"].values())
    with (args.case / "inputs/source_rows.csv").open(encoding="utf-8-sig", newline="") as stream:
        source = {r["point_id"]: r for r in csv.DictReader(stream)}
    points = []
    for i, offset in enumerate(config["xi_offsets"], 1):
        if target["segment_left"] <= center + offset <= target["segment_right"]:
            points.append(dict(point_id=f"{args.curve}_p{i:02d}", curve_id=args.curve,
                meson=target["meson"], cut=target["cut"], window=target["window"],
                T_MeV=target["T_MeV"], muB_MeV=target["muB_MeV"], xi=center+offset,
                xi_mott=center, delta_xi=offset, kind="mott_centered", original_point_id=""))
    for i, name in enumerate(("point_left", "point_right"), 1):
        point = source[target[name]]
        xi = float(point["xi"])
        points.append(dict(point_id=f"{args.curve}_c{i}", curve_id=args.curve,
            meson=target["meson"], cut=target["cut"], window=target["window"],
            T_MeV=target["T_MeV"], muB_MeV=target["muB_MeV"], xi=xi,
            xi_mott=center, delta_xi=xi-center, kind="original_control", original_point_id=point["point_id"]))
    args.output.mkdir(parents=True, exist_ok=False)
    with (args.output / "transport_plan.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(points[0]))
        writer.writeheader()
        writer.writerows(points)
    summary = dict(status="completed_location_reused_for_pilot", location=location,
                   live_location_snapshot_sha256=hashlib.sha256(raw).hexdigest(),
                   plan_sha256=hashlib.sha256((args.output / "transport_plan.csv").read_bytes()).hexdigest())
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
