#!/usr/bin/env python3
"""Compact remote diagnostics without changing or replacing any retained results."""
from pathlib import Path
from collections import defaultdict
import csv
import json
import math
import sys


def read_csv(path):
    with path.open(encoding="utf-8-sig",newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    case=Path(sys.argv[1]);mode=sys.argv[2]
    if mode=="events":
        rows=[]
        plan=read_csv(case/"inputs/transport_plan.csv")
        folders=[case/"inputs/previous",case/"left_s6_v1",case/"right_s6_v1"]
        for folder in folders:
            obs=read_csv(folder/"observables.csv")
            rates=read_csv(folder/"rates.csv")
            for point in obs:
                ident=next((p["point_id"] for p in plan if abs(float(p["xi"])-float(point["xi"]))<1e-12),None)
                if ident is None:continue
                for row in rates:
                    if row["point_id"]!=point["point_id"] or row["channel"] not in ("uu_to_uu","ud_to_ud","us_to_us","uubar_to_uubar"):continue
                    rows.append(dict(point_id=ident,xi=float(point["xi"]),mu=float(point["m_u"]),
                        channel=row["channel"],events=int(row["events"]),unique=int(row["unique_events"]),
                        rate=float(row["rate"]),error=float(row["estimated_error"])))
        print(json.dumps(sorted(rows,key=lambda r:(r["channel"],r["xi"]))))
    elif mode=="peak":
        paths=(case/"spectrum_v2/peak_areas.csv",case/"thermal_audit_v1/thermal_peak_areas.csv")
        areas=[r for r in read_csv(paths[0]) if r["channel"]=="uubar_to_uubar"]
        indexed={(r["point_id"],r["order"]):r for r in areas}
        checks=[]
        for r in read_csv(paths[1]):
            if r["channel"]!="uubar_to_uubar":continue
            a,b=indexed[r["point_id"],"16"],indexed[r["point_id"],"32"]
            checks.append(dict(point_id=r["point_id"],x=float(r["x_peak_inv_fm2"]),
                area_order_relative=float(a["mixed_square_area_dimensionless"])/float(b["mixed_square_area_dimensionless"])-1,
                thermal_discrete_to_continuum=float(r["discrete_exact_mixed_square_rate_fm2"])/float(r["continuous_mixed_square_rate_fm2"]),
                peak_interpolation_l1_fm2=float(r["discrete_interpolation_l1_fm2"])))
        print(json.dumps(checks))
    elif mode=="background":
        spectra={r["point_id"]:r for r in read_csv(case/"spectrum_v2/background.csv")}
        plan=read_csv(case/"inputs/transport_plan.csv")
        checks=[]
        for folder in (case/"inputs/previous",case/"left_s6_v1",case/"right_s6_v1"):
            for row in read_csv(folder/"observables.csv"):
                point=next((p for p in plan if abs(float(p["xi"])-float(row["xi"]))<1e-12),None)
                if point is None:continue
                spec=spectra[point["point_id"]]
                difference={k:float(row[k])-float(spec[k]) for k in ("m_u","m_d","m_s","Phi","Phibar")}
                checks.append(dict(point_id=point["point_id"],differences=difference,
                    max_absolute=max(abs(v) for v in difference.values()),
                    relative_pair_threshold_difference=(float(row["m_u"])**2/float(spec["m_u"])**2-1)))
        print(json.dumps(checks))
    elif mode=="live_rates":
        native={}
        plan=read_csv(case/"inputs/transport_plan.csv")
        for folder in (case/"inputs/previous",case/"left_s6_v1",case/"right_s6_v1"):
            points=read_csv(folder/"observables.csv")
            ids={r["point_id"]:next((p["point_id"] for p in plan if abs(float(p["xi"])-float(r["xi"]))<1e-12),None) for r in points}
            for r in read_csv(folder/"rates.csv"):
                if ids.get(r["point_id"]):native[ids[r["point_id"]],r["channel"]]=r
        reference={}
        for name in ("reference_center_s7_v1","reference_left_s7_v1","reference_right_s7_v1"):
            path=case/name/"rates.csv"
            if path.exists():
                for r in read_csv(path):reference[r["point_id"],r["channel"]]=r
        channels={k[1] for k in reference}
        result=[]
        for channel in sorted(channels):
            ids=("sig_m020","sig_c000","sig_p020")
            if not all((i,channel) in reference for i in ids):continue
            item=dict(channel=channel)
            for tag,source in (("s6",native),("s7",reference)):
                values=[float(source[i,channel]["rate"]) for i in ids]
                sl=math.log1p((values[1]-values[0])/values[0])/2e-5
                sr=math.log1p((values[2]-values[1])/values[1])/2e-5
                item[tag+"_left"]=sl;item[tag+"_right"]=sr;item[tag+"_jump"]=sr-sl
            item["maximum_point_relative_change"]=max(abs(float(reference[i,channel]["rate"])/float(native[i,channel]["rate"])-1) for i in ids)
            result.append(item)
        print(json.dumps(dict(available_reference_channels=len(reference),matched_triplets=result)))
    else:
        raise ValueError(mode)


if __name__=="__main__":main()
