"""Summarize controlled sampling changes without modifying historical curves."""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import statistics

from analyze_threshold_cache_upgrade import audit_run


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(line for line in stream if not line.startswith("#")))


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path,data):
    path.write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")


def write_rows(path,data):
    if not data:
        return
    with path.open("w",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(dict.fromkeys(key for row in data for key in row)))
        writer.writeheader();writer.writerows(data)


def curvature(xs,ys):
    residuals=[]
    for i in range(1,len(xs)-1):
        t=(xs[i]-xs[i-1])/(xs[i+1]-xs[i-1])
        residuals.append(abs(math.log(ys[i])-(1-t)*math.log(ys[i-1])-t*math.log(ys[i+1])))
    return max(residuals,default=0.0)


def normalized_record(row,run):
    result=dict(row,run=run)
    for side,legacy,level in (("coarse","dense192_rate",192),("fine","dense384_rate",384)):
        field="reference_"+side+"_rate"
        if field not in result:
            result[field]=row[legacy]
            result["reference_"+side+"_level"]=str(level)
    return result


def rate_value(row,column):
    if column in row:return float(row[column])
    return float(row[{"reference_coarse_rate":"dense192_rate","reference_fine_rate":"dense384_rate"}[column]])


def reconstruct_species(point, contributions, records, species):
    selected=[r for r in contributions if r["species"]==species]
    required=[r for r in selected if r["recompute"]=="true"]
    available=all((r["point_id"],r["channel"]) in records for r in required)
    original=float(point["tau_"+species])
    result=dict(point_id=point["point_id"],window_ids=point["window_ids"],
        muB_MeV=float(point["muB_MeV"]),alpha_T=float(point["alpha_T"]),
        xi=float(point["xi"]),species=species,source_tau=original,
        selected_channels=len(required),complete_selected_channels=available)
    for flag in ("original_replayed","reference_converged","adaptive_validated"):
        result[flag]=available and all(records[r["point_id"],r["channel"]][flag]=="true" for r in required)
    methods={"replayed":"replayed_rate","adaptive":"adaptive_rate", "reference_coarse":"reference_coarse_rate","reference":"reference_fine_rate"}
    for method,column in methods.items():
        if not available:
            # An unfinished channel set must not look like a completed update.
            result[method+"_tau"]=None
            continue
        inverse=1/original+sum(float(r["multiplicity"])*float(r["density"])*
            (rate_value(records[r["point_id"],r["channel"]],column)-float(r["rate"])) for r in required)
        result[method+"_tau"]=1/inverse if inverse>0 else None
    return result


def main(root,out,run_names,figure_out=None,refinements=()):
    root,out=root.resolve(),out.resolve()
    figure_out=out if figure_out is None else figure_out.resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError("Use a new empty output directory")
    if figure_out!=out and figure_out.exists() and any(figure_out.iterdir()):
        raise ValueError("Use a new empty figure directory")
    selection=read_json(root/"anomaly_selection_v1/selection.json")
    points=rows(root/"anomaly_selection_v1/review_points.csv")
    plan=rows(root/"anomaly_selection_v1/channel_plan.csv")
    audits=[];records={};failures=[];initial_records={};superseded=[]
    selected_input_hashes={r["sha256"] for r in selection["outputs"]}
    for name,is_refinement in [(name,False) for name in run_names]+[(name,True) for name in refinements]:
        run=root/name
        audit=audit_run(run,allow_interrupted=True)
        if not selected_input_hashes.issubset(audit["input_hashes"]):
            raise ValueError(f"Run inputs differ from the selected historical backgrounds: {name}")
        audits.append(audit)
        if (run/"rates.csv").exists():
            for r in rows(run/"rates.csv"):
                key=(r["point_id"],r["channel"])
                if is_refinement:
                    if key not in records:
                        raise ValueError(f"Refinement has no original completed case: {key}")
                    old=records[key]
                    if old["run"]!=r["superseded_run"]:
                        raise ValueError(f"Refinement donor does not match the active result: {key}")
                    if any(float(old[k])!=float(r[k]) for k in ("source_rate","muB_MeV","alpha_T","xi")):
                        raise ValueError(f"Refinement changes the retained state or source rate: {key}")
                    superseded.append(dict(old,replaced_by=name))
                elif key in records:
                    raise ValueError(f"Duplicate case requires an explicit refinement policy: {key}")
                records[key]=normalized_record(r,name)
                if not is_refinement:initial_records[key]=records[key]
        if (run/"failures.csv").exists():
            failures.extend(dict(r,run=name) for r in rows(run/"failures.csv"))
    expected={(r["point_id"],r["channel"]) for r in plan if r["recompute"]=="true"}
    missing=sorted(expected-set(records))
    contribution_index=collections.defaultdict(list)
    for r in plan:contribution_index[r["point_id"]].append(r)
    tau_rows=[]
    for p in points:
        pd=contribution_index[p["point_id"]]
        for species in ("u","d","s","ubar","dbar","sbar"):
            tau_rows.append(reconstruct_species(p,pd,records,species))
    windows=[]
    for w in selection["windows"]:
        ids={p["point_id"] for p in points if w["window_id"] in p["window_ids"].split(";")}
        selected=[r for r in records.values() if r["point_id"] in ids]
        by_species={q:sorted([r for r in tau_rows if r["point_id"] in ids and r["species"]==q],key=lambda r:r["xi"])
                    for q in ("u","s","ubar","sbar")}
        diagnostic=max(("u","s"),key=lambda q:curvature([r["xi"] for r in by_species[q]],[r["source_tau"] for r in by_species[q]]))
        rs=by_species[diagnostic];xs=[r["xi"] for r in rs]
        complete=all(r["complete_selected_channels"] and r["reference_tau"] is not None for r in rs)
        old_curvature=curvature(xs,[r["source_tau"] for r in rs])
        new_curvature=curvature(xs,[r["reference_tau"] for r in rs]) if complete else None
        expected_jobs=sum(k[0] in ids for k in expected)
        all_jobs_complete=len(selected)==expected_jobs
        windows.append(dict(window_id=w["window_id"],muB=w["panel"],alpha_T=w["series"],
            left_xi=min(xs),right_xi=max(xs),diagnostic_species=diagnostic,
            jobs=len(selected),expected_jobs=expected_jobs,
            replay_failures=sum(r["original_replayed"]!="true" for r in selected),
            reference_unresolved=sum(r["reference_converged"]!="true" for r in selected),
            adaptive_unresolved=sum(r["adaptive_validated"]!="true" for r in selected),
            selected_channels_complete=all_jobs_complete,
            selected_channels_validated=all_jobs_complete and all(r["original_replayed"]==r["reference_converged"]==r["adaptive_validated"]=="true" for r in selected),
            raw_log_curvature=old_curvature,recomputed_log_curvature=new_curvature,
            max_tau_relative_change=max(abs(r["reference_tau"]/r["source_tau"]-1) for r in rs) if complete else None,
            interpretation="rapid-change control" if w["control"] else "assess convergence before physical attribution"))
    performance=[]
    for name in ("thermal_performance_v1","native_performance_v1","grouped_performance_v2"):
        if not (root/name/"status.json").exists():continue
        if read_json(root/name/"status.json").get("complete") is not True:continue
        audits.append(audit_run(root/name))
        groups=collections.defaultdict(list)
        for r in rows(root/name/"timings.csv"):
            groups[r.get("muB_MeV","0.0"),r["xi"],r.get("channel",r.get("process")),r["method"]].append(r)
        for k,rs in groups.items():
            performance.append(dict(run=name,muB_MeV=k[0],xi=k[1],channel=k[2],method=k[3],n=len(rs),
                wall_median=statistics.median(float(r["wall_seconds"]) for r in rs),
                cpu_median=statistics.median(float(r["cpu_seconds"]) for r in rs),
                max_relative_error=max(float(r["relative_error"]) for r in rs)))
    out.mkdir(parents=True,exist_ok=True)
    write_rows(out/"rates.csv",list(records.values()))
    write_rows(out/"tau_reconstruction.csv",tau_rows)
    write_rows(out/"window_summary.csv",windows)
    write_rows(out/"performance.csv",performance)
    write_rows(out/"failures.csv",failures)
    write_rows(out/"superseded_rates.csv",superseded)
    write_json(out/"audit.json",audits)
    summary=dict(status="analysis_only",expected_jobs=len(expected),successful_jobs=len(records),missing_jobs=missing,
        backgrounds=len(points),refined_cases=len({(r["point_id"],r["channel"]) for r in superseded}),
        initial_reference_unresolved=sum(r["reference_converged"]!="true" for r in initial_records.values()),
        initial_adaptive_unresolved=sum(r["adaptive_validated"]!="true" for r in initial_records.values()),
        original_replay_failures=sum(r["original_replayed"]!="true" for r in records.values()),
        reference_unresolved=sum(r["reference_converged"]!="true" for r in records.values()),
        adaptive_unresolved=sum(r["adaptive_validated"]!="true" for r in records.values()),windows=windows,
        max_adaptive_vs_reference_tau=max((abs(r["adaptive_tau"]/r["reference_tau"]-1) for r in tau_rows
            if r["adaptive_tau"] is not None and r["reference_tau"] is not None),default=None),
        source_unchanged=all(sha(Path(r["path"]))==r["sha256"] for r in selection["sources"]),
        scope="mode A selected windows; partial channel replacement at retained backgrounds, no smoothing or background solve")
    write_json(out/"summary.json",summary)
    figure_out.mkdir(parents=True,exist_ok=True)
    make_plot(figure_out,tau_rows,windows)
    antiquark_windows=[]
    for window in windows:
        selected=[r for r in tau_rows if window["window_id"] in r["window_ids"].split(";")]
        species=max(("ubar","sbar"),key=lambda q:curvature(
            [r["xi"] for r in sorted(selected,key=lambda r:r["xi"]) if r["species"]==q],
            [r["source_tau"] for r in sorted(selected,key=lambda r:r["xi"]) if r["species"]==q]))
        antiquark_windows.append(dict(window,diagnostic_species=species))
    make_plot(figure_out,tau_rows,antiquark_windows,basename="anomaly_antiquark_windows")
    (out/"source").mkdir()
    for p in (Path(__file__),Path(__file__).with_name("analyze_threshold_cache_upgrade.py")):
        shutil.copyfile(p,out/"source"/p.name)
    write_json(out/"manifest.json",dict(status="analysis_only",runs=run_names,refinements=list(refinements),
        sources=[dict(path=str(p),sha256=sha(p)) for p in (Path(__file__),root/"anomaly_selection_v1/selection.json")],
        figures=[dict(path=str(figure_out/f"{name}.{suffix}"),sha256=sha(figure_out/f"{name}.{suffix}"))
            for name in ("anomaly_windows","anomaly_antiquark_windows") for suffix in ("png","pdf")],
        outputs=[dict(path=p.name,sha256=sha(p)) for p in sorted(out.iterdir()) if p.is_file()]))
    print(json.dumps({k:v for k,v in summary.items() if k!="windows"},ensure_ascii=False))
    for w in windows:print(w)


def make_plot(out,data,windows,basename="anomaly_windows"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(math.ceil(len(windows)/3),3,figsize=(12,10),constrained_layout=True)
    for index,(ax,w) in enumerate(zip(axes.flat,windows)):
        rs=sorted([r for r in data if w["window_id"] in r["window_ids"].split(";") and r["species"]==w["diagnostic_species"]],key=lambda r:r["xi"])
        for col,label,style in (("source_tau","Raw production","o--"),("adaptive_tau","Adaptive cache",".-"),("reference_tau","Fixed fine reference","x:")):
            ax.plot([r["xi"] for r in rs],[r[col] for r in rs],style,label=label,lw=1.2,ms=4)
        unresolved=[r for r in rs if r["reference_tau"] is not None and not r["reference_converged"]]
        if unresolved:
            ax.scatter([r["xi"] for r in unresolved],[r["reference_tau"] for r in unresolved],
                s=35,facecolors="none",edgecolors="red",marker="s",label="Reference unresolved",zorder=5)
        mu=float(w["muB"].removeprefix("muB"));alpha=float(w["alpha_T"].removeprefix("alpha"))
        symbol={"u":"u","s":"s","ubar":r"\bar u","sbar":r"\bar s"}[w["diagnostic_species"]]
        ax.set_title(f"({chr(97+index)}) "+fr"$\mu_B={mu:g}\,\mathrm{{MeV}},\ \alpha_T={alpha:g},\ q={symbol}$",fontsize=9)
        ax.set(xlabel=r"$\xi$",ylabel=r"$\tau$ [fm]")
        ax.tick_params(direction="in",top=True,right=True)
        ax.legend(fontsize=7,frameon=False)
    for ax in list(axes.flat)[len(windows):]:ax.set_visible(False)
    for suffix in ("png","pdf"):fig.savefig(out/f"{basename}.{suffix}",dpi=180)
    plt.close(fig)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_root",type=Path);parser.add_argument("new_output_directory",type=Path)
    parser.add_argument("runs",nargs="+")
    parser.add_argument("--figure-output-directory",type=Path)
    parser.add_argument("--refinements",nargs="+",default=[])
    args=parser.parse_args();main(args.run_root,args.new_output_directory,args.runs,args.figure_output_directory,args.refinements)
