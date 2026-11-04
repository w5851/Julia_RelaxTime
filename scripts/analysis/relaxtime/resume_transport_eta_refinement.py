#!/usr/bin/env python3
"""Freeze complete points from interrupted diagnostic workers; quote-remote only."""
import csv
import hashlib
import json
import math
from pathlib import Path
import sys


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    with path.open(encoding="utf-8-sig",newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    case=Path(sys.argv[1]).resolve()
    destination=case/"resume_v2.json"
    if destination.exists():
        raise FileExistsError(destination)
    plan=rows(case/"locate_v3/transport_plan.csv")
    originals={r["point_id"]:r for r in plan}
    jobs=[]
    partial_runs=[]
    fields=("tau_u","tau_d","tau_s","tau_ubar","tau_dbar","tau_sbar","eta_over_s","zeta_over_s","sigma_over_T")
    for variant in ("fine","tighter"):
        for curve in ("a_mu0","a_mu450","b_mu0"):
            name=f"{variant}_{curve}_v1"
            folder=case/name
            accepted=[]
            rejected=[]
            if (folder/"observables.csv").exists():
                current=rows(folder/"observables.csv")
                rates=rows(folder/"rates.csv")
                for row in current:
                    point=row.get("point_id","")
                    try:
                        channels=[r for r in rates if r.get("point_id")==point]
                        assert len(channels)==len({r["channel"] for r in channels})==17
                        assert all(math.isfinite(float(r["rate"])) and float(r["rate"])>0 for r in channels)
                        assert all(math.isfinite(float(row[k])) and float(row[k])>0 for k in fields)
                        assert float(row["background_residual"])<=1e-10
                        assert point in originals and abs(float(row["xi"])-float(originals[point]["xi"]))<1e-12
                    except (ValueError,TypeError,KeyError,AssertionError) as err:
                        rejected.append(dict(point_id=point,reason=repr(err)))
                    else:
                        accepted.append(point)
                assert len(accepted)==len(set(accepted)),"duplicate accepted points"
                summary=dict(status="interrupted_with_complete_points",accepted_point_ids=accepted,
                    rejected_rows=rejected,variant=variant,curve_id=curve,
                    interruption_record_sha256=sha(case/"interruption_v1.json"),
                    analysis_script_sha256=sha(case/"run_transport_eta_refinement.jl"),
                    input_config_sha256=sha(case/"inputs_v2/config.json"),
                    transport_plan_sha256=sha(case/"locate_v3/transport_plan.csv"),
                    output_sha256={p.name:sha(p) for p in folder.iterdir() if p.is_file()})
                (folder/"partial_summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
                partial_runs.append(name)
            wanted=[r["point_id"] for r in plan if r["curve_id"]==curve and
                    (variant=="fine" or r["point_id"] in (curve+"_p01",curve+"_p02"))]
            missing=[point for point in wanted if point not in accepted]
            if missing:
                jobs.append(dict(name=f"{variant}_{curve}_v2",variant=variant,curve_id=curve,point_ids=missing))
    result=dict(status="ready",jobs=jobs,partial_runs=partial_runs,
        change="transport integrals unchanged; reuse independent locator instead of per-point multistart pole probes",
        plan_sha256=sha(case/"locate_v3/transport_plan.csv"),script_sha256=sha(Path(__file__)))
    destination.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
