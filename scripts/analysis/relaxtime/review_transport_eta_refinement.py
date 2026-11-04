#!/usr/bin/env python3
"""Compare raw coarse/fine xi transport curves. Execute on quote-remote."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

OBSERVABLES = ("tau_u", "tau_d", "tau_s", "tau_ubar", "tau_dbar", "tau_sbar",
               "eta_over_s", "zeta_over_s", "sigma_over_T")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def number(row, key):
    value = float(row[key])
    if not math.isfinite(value):
        raise ValueError(f"Nonfinite {key} in {row.get('point_id', '')}")
    return value


def intervals(rows, label, curve_id):
    output = []
    for left, right in zip(rows, rows[1:]):
        lo, hi = number(left, "xi"), number(right, "xi")
        if hi <= lo:
            raise ValueError("Duplicate or unordered xi values")
        item = dict(curve_id=curve_id, grid=label, point_left=left["point_id"],
                    point_right=right["point_id"], xi_left=lo, xi_right=hi,
                    xi_midpoint=(lo+hi)/2, delta_xi=hi-lo,
                    xi_mott=number(left, "xi_mott"),
                    delta_xi_midpoint=(lo+hi)/2-number(left,"xi_mott"))
        for key in OBSERVABLES:
            a, b = number(left, key), number(right, key)
            if min(a, b) <= 0:
                raise ValueError("Log slopes require positive observables")
            item[key+"_log_slope"] = math.log(b/a)/(hi-lo)
            item[key+"_change_percent"] = 100*(b/a-1)
        output.append(item)
    return output


def collect_runs(case, names, preview):
    rows, channels, hashes = [], [], {}
    for name in names:
        folder = case/name
        summary_path = folder/"summary.json"
        partial = False
        if not summary_path.exists() and (folder/"partial_summary.json").exists():
            summary_path=folder/"partial_summary.json"
            partial=True
        if summary_path.exists():
            summary = read_json(summary_path)
            assert summary["status"] == ("interrupted_with_complete_points" if partial else "completed_diagnostic")
            for filename, sha in summary["output_sha256"].items():
                assert digest(folder/filename) == sha, (name, filename)
            hashes[str(summary_path.relative_to(case))] = digest(summary_path)
        elif not preview:
            raise ValueError(f"Incomplete run: {name}")
        obs_path = folder/"observables.csv"
        if not obs_path.exists():
            continue
        local_rows = read_csv(obs_path)
        if partial:
            accepted=set(summary["accepted_point_ids"])
            local_rows=[r for r in local_rows if r["point_id"] in accepted]
        ids = {r["point_id"] for r in local_rows}
        local_channels = [r for r in read_csv(folder/"rates.csv") if r["point_id"] in ids]
        for row in local_rows:
            point_channels = [r for r in local_channels if r["point_id"] == row["point_id"]]
            assert len(point_channels) == len({r["channel"] for r in point_channels}) == 17
            assert all(number(r, "rate") > 0 for r in point_channels)
            assert number(row, "background_residual") <= 1e-10
            for key in OBSERVABLES:
                assert number(row, key) > 0
        rows.extend(local_rows)
        channels.extend(local_channels)
        hashes[str(obs_path.relative_to(case))] = digest(obs_path)
    assert len(rows) == len({r["point_id"] for r in rows}), "duplicate points across runs"
    return rows, channels, hashes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("--runs", required=True, help="Comma-separated fine runs, including optional extra points")
    parser.add_argument("--locate", default="locate_v3")
    parser.add_argument("--inputs", default="inputs_v2")
    parser.add_argument("--comparison-runs", default="")
    parser.add_argument("--peak-plan", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--preview", action="store_true")
    args = parser.parse_args()
    case = args.case.resolve()
    rows, channels, hashes = collect_runs(case, args.runs.split(","), args.preview)
    locations = read_json(case/args.locate/"locations.json")
    originals = {r["point_id"]: r for r in read_csv(case/args.inputs/"source_rows.csv")}
    comparisons, chords, fine_slopes, coarse_slopes, summaries = [], [], [], [], []
    for loc in locations:
        curve_id = loc["curve_id"]
        curve = sorted((r for r in rows if r["curve_id"] == curve_id), key=lambda r: number(r, "xi"))
        controls = [r for r in curve if r["original_point_id"]]
        fs, cs = intervals(curve, "fine", curve_id), intervals(controls, "coarse_same_precision", curve_id)
        fine_slopes.extend(fs)
        coarse_slopes.extend(cs)
        for row in controls:
            old = originals[row["original_point_id"]]
            assert abs(number(row, "xi")-number(old, "xi")) < 1e-12
            for key in OBSERVABLES:
                comparisons.append(dict(curve_id=curve_id, point_id=row["point_id"],
                    original_point_id=old["point_id"], xi=number(row, "xi"), observable=key,
                    original=number(old, key), recomputed=number(row, key),
                    relative_change_percent=100*(number(row, key)/number(old, key)-1)))
        for row in curve:
            if row["original_point_id"]:
                continue
            x = number(row, "xi")
            for left, right in zip(controls, controls[1:]):
                lo, hi = number(left, "xi"), number(right, "xi")
                if lo <= x <= hi:
                    fraction = (x-lo)/(hi-lo)
                    for key in OBSERVABLES:
                        # Explicit diagnostic chord; never substituted for data.
                        chord = (1-fraction)*number(left, key)+fraction*number(right, key)
                        chords.append(dict(curve_id=curve_id, point_id=row["point_id"], xi=x,
                            observable=key, value=number(row, key), coarse_chord=chord,
                            deviation_percent=100*(number(row, key)/chord-1)))
                    break
        peaks = {}
        for key in OBSERVABLES:
            field = key+"_log_slope"
            if fs:
                peak = max(fs, key=lambda r: r[field])
                coarse_peak = max((r[field] for r in cs), default=None)
                deviation = max((abs(r["deviation_percent"]) for r in chords
                                 if r["curve_id"] == curve_id and r["observable"] == key), default=None)
                peaks[key] = dict(fine_max_log_slope=peak[field], xi_left=peak["xi_left"],
                    xi_right=peak["xi_right"], point_left=peak["point_left"], point_right=peak["point_right"],
                    recomputed_coarse_max_log_slope=coarse_peak,
                    peak_over_coarse=(peak[field]/coarse_peak if coarse_peak else None),
                    max_abs_coarse_chord_deviation_percent=deviation)
        summaries.append(dict(curve_id=curve_id, points=len(curve), coarse_controls=len(controls),
            xi_mott=loc["levels"]["base"]["xi"],
            A64_xi_shift=loc["levels"]["A64"]["xi"]-loc["levels"]["base"]["xi"],
            background_xi_shift=loc["levels"]["background_ref"]["xi"]-loc["levels"]["A64"]["xi"],
            pair_quadrature_xi_shift=loc["levels"]["pair_ref"]["difference"], observables=peaks))
    tighter_rows, tighter_comparisons, tighter_slopes, slope_checks, xi_checks = [], [], [], [], []
    if args.comparison_runs:
        tighter_rows, _, additional_hashes = collect_runs(case, args.comparison_runs.split(","), args.preview)
        hashes.update(additional_hashes)
        fine_by_id = {r["point_id"]: r for r in rows}
        for row in tighter_rows:
            if row["point_id"] not in fine_by_id:
                assert row["point_kind"] in ("peak_interval_midpoint","threshold_interval_midpoint"), "unknown new reference point"
                continue
            fine = fine_by_id[row["point_id"]]
            assert abs(number(row,"xi")-number(fine,"xi")) < 1e-12
            for key in OBSERVABLES:
                tighter_comparisons.append(dict(curve_id=row["curve_id"],point_id=row["point_id"],
                    xi=number(row,"xi"),observable=key,fine=number(fine,key),tighter=number(row,key),
                    relative_change_percent=100*(number(row,key)/number(fine,key)-1)))
        for loc in locations:
            curve_id=loc["curve_id"]
            curve=sorted((r for r in tighter_rows if r["curve_id"]==curve_id),key=lambda r:number(r,"xi"))
            ts=intervals(curve,"tighter",curve_id)
            tighter_slopes.extend(ts)
            matched=[r for r in curve if r["point_id"] in fine_by_id]
            fs=intervals([fine_by_id[r["point_id"]] for r in matched],"fine_matched",curve_id)
            mt=intervals(matched,"tighter_matched",curve_id)
            for fine,tighter in zip(fs,mt):
                for key in OBSERVABLES:
                    a,b=fine[key+"_log_slope"],tighter[key+"_log_slope"]
                    slope_checks.append(dict(curve_id=curve_id,point_left=fine["point_left"],
                        point_right=fine["point_right"],xi_left=fine["xi_left"],xi_right=fine["xi_right"],
                        observable=key,fine_log_slope=a,tighter_log_slope=b,absolute_difference=b-a,
                        relative_difference_percent=100*(b/a-1) if a else None))
        if args.peak_plan:
            selection=read_json(args.peak_plan/"summary.json")
            hashes[str(args.peak_plan/"summary.json")]=digest(args.peak_plan/"summary.json")
            by_id={r["point_id"]:r for r in tighter_rows}
            for item in selection["selection"]:
                left,mid,right=(by_id[item[k]] for k in ("point_left","midpoint_id","point_right"))
                lo,middle,hi=(number(r,"xi") for r in (left,mid,right))
                assert abs(middle-(lo+hi)/2)<1e-12
                for key in OBSERVABLES:
                    a,b,c=(number(r,key) for r in (left,mid,right))
                    whole=math.log(c/a)/(hi-lo)
                    half_left=math.log(b/a)/(middle-lo)
                    half_right=math.log(c/b)/(hi-middle)
                    assert abs((half_left*(middle-lo)+half_right*(hi-middle))/(hi-lo)-whole)<1e-9
                    xi_checks.append(dict(curve_id=item["curve_id"],
                        refinement_kind=item.get("refinement_kind","peak_interval"),
                        observable=key,xi_left=lo,xi_midpoint=middle,
                        xi_right=hi,parent_log_slope=whole,left_half_log_slope=half_left,right_half_log_slope=half_right,
                        peak_half_over_parent=max(half_left,half_right)/whole if whole else None,
                        midpoint_chord_deviation_percent=100*(b/((a+c)/2)-1)))
    summary=dict(status="partial_preview" if args.preview else "completed_diagnostic_review",
        points=len(rows),channels=len(channels),curves=summaries,input_hashes=hashes,
        location_summary_sha256=digest(case/args.locate/"summary.json"),
        script_sha256=digest(Path(__file__)),comparison_points=len(tighter_rows),
        definitions=dict(log_slope="ln(y_right/y_left)/(xi_right-xi_left)",
            chord="explicit straight chord between recomputed original nodes, used only for comparison",
            peak="largest observed secant slope on supplied grid, not a proven continuous supremum"),
        limits="No smoothing or formal-data replacement; eta association is not an amplitude ablation")
    if not args.preview:
        assert args.output is not None
        args.output.mkdir(parents=True,exist_ok=False)
        expected={r["point_id"] for r in read_csv(case/args.locate/"transport_plan.csv")}
        assert expected <= {r["point_id"] for r in rows}, "missing initial plan points"
        for name,table in (("observables",rows),("channels",channels),("fine_slopes",fine_slopes),
                ("coarse_slopes",coarse_slopes),("original_controls",comparisons),("coarse_chord",chords),
                ("tighter_observables",tighter_rows),("tighter_comparison",tighter_comparisons),
                ("tighter_slopes",tighter_slopes),("slope_convergence",slope_checks),("xi_refinement",xi_checks)):
            write_csv(args.output/(name+".csv"),table)
        summary["output_sha256"]={p.name:digest(p) for p in args.output.iterdir() if p.is_file()}
        (args.output/"review.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(summary,indent=2,allow_nan=False))


if __name__ == "__main__":
    main()
