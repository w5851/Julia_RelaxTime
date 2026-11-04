#!/usr/bin/env python3
"""Select peak and exact-threshold xi midpoints; execute on quote-remote."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def rows(path):
    with path.open(encoding="utf-8-sig",newline="") as stream:
        return list(csv.DictReader(stream))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case",type=Path)
    parser.add_argument("review",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    review=json.loads((args.review/"review.json").read_text(encoding="utf-8-sig"))
    assert review["status"]=="completed_diagnostic_review"
    for file,digest in review["output_sha256"].items():
        assert sha(args.review/file)==digest,file
    points={r["point_id"]:r for r in rows(args.case/"locate_v3/transport_plan.csv")}
    # The first precision batch contains the two closest points on each side
    # of the crossing. Reuse any already scheduled peak endpoint.
    scheduled={curve+suffix for curve in ("a_mu0","a_mu450","b_mu0") for suffix in ("_p01","_p02")}
    extra=[]
    selections=[]
    for curve in review["curves"]:
        peak=curve["observables"]["tau_ubar"]
        left,right=points[peak["point_left"]],points[peak["point_right"]]
        lo,hi=float(left["xi"]),float(right["xi"])
        mid=(lo+hi)/2
        assert lo<mid<hi and math.isfinite(mid)
        for row in (left,right):
            if row["point_id"] not in scheduled:
                extra.append(row.copy())
        center=float(left["xi_mott"])
        new=dict(left,point_id=curve["curve_id"]+"_peak_mid",xi=mid,delta_xi=mid-center,
                 kind="peak_interval_midpoint",original_point_id="")
        extra.append(new)
        selections.append(dict(curve_id=curve["curve_id"],refinement_kind="peak_interval",
            point_left=left["point_id"],point_right=right["point_id"],
            midpoint_id=new["point_id"],xi_left=lo,xi_right=hi,xi_midpoint=mid,
            fine_secant_log_slope=peak["fine_max_log_slope"],observable="tau_ubar"))
        threshold_left=points[curve["curve_id"]+"_p01"]
        threshold_right=points[curve["curve_id"]+"_p02"]
        threshold_lo,threshold_hi=float(threshold_left["xi"]),float(threshold_right["xi"])
        assert threshold_lo<center<threshold_hi
        assert abs(center-(threshold_lo+threshold_hi)/2)<1e-12
        threshold=dict(threshold_left,point_id=curve["curve_id"]+"_mott_mid",xi=center,delta_xi=0.0,
            kind="threshold_interval_midpoint",original_point_id="")
        extra.append(threshold)
        selections.append(dict(curve_id=curve["curve_id"],refinement_kind="threshold_interval",
            point_left=threshold_left["point_id"],point_right=threshold_right["point_id"],
            midpoint_id=threshold["point_id"],xi_left=threshold_lo,xi_right=threshold_hi,xi_midpoint=center,
            fine_secant_log_slope=None,observable="all_observables"))
    args.output.mkdir(parents=True,exist_ok=False)
    with (args.output/"transport_plan.csv").open("w",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(extra[0]))
        writer.writeheader();writer.writerows(extra)
    summary=dict(status="completed_diagnostic_selection",selection=selections,transport_plan_points=len(extra),
        purpose="Halve the steepest observed tau_ubar interval and the central crossing interval; all endpoints and midpoints compared at rtol=1e-6",
        source_review_sha256=sha(args.review/"review.json"),
        parent_location_summary_sha256=sha(args.case/"locate_v3/summary.json"),
        analysis_script_sha256=sha(Path(__file__)),output_sha256={"transport_plan.csv":sha(args.output/"transport_plan.csv")})
    (args.output/"summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
