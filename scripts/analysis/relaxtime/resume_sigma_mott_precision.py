#!/usr/bin/env python3
"""Complete the sigma fixed-background precision comparison after its pilot."""
from pathlib import Path
from types import SimpleNamespace
from concurrent.futures import ThreadPoolExecutor
import argparse
import importlib.util
import json
import math
import os
import subprocess
import sys
import time


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case",type=Path)
    args=parser.parse_args()
    case=args.case.resolve()
    spec=importlib.util.spec_from_file_location("sigma_runner",case/"run_sigma_mott_refinement_v2.py")
    runner=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    path=case/"reference_controls_v2.status.json"
    with path.open("x",encoding="utf-8") as stream:
        json.dump(dict(status="waiting_for_pilot",pid=os.getpid(),started_utc=runner.now()),stream)
    state=dict(status="waiting_for_pilot",pid=os.getpid(),started_utc=runner.now(),completed=[])
    try:
        while True:
            pilot=runner.read_json(case/"reference_center_s7_v1.status.json")
            if pilot["status"]=="failed":
                raise RuntimeError("Center precision reference failed")
            rate_file=case/"reference_center_s7_v1"/"rates.csv"
            partial=runner.read_csv(rate_file) if rate_file.exists() else []
            if pilot["status"]=="completed_diagnostic" or len(partial)>=2:
                break
            state.update(updated_utc=runner.now(),free_GiB=runner.free_gib())
            runner.write_json(path,state)
            time.sleep(20)
        old=runner.read_csv(case/"center_s7_v2"/"rates.csv")
        new={r["channel"]:r for r in runner.read_csv(case/"reference_center_s7_v1"/"rates.csv")}
        checks=[]
        for row in old:
            if row["channel"] not in new:
                continue
            r=new[row["channel"]]
            delta=float(r["rate"])-float(row["rate"])
            allowance=float(r["estimated_error"])+float(row["estimated_error"])+1e-12
            passed=abs(delta)<=3*allowance
            checks.append(dict(channel=row["channel"],native_rate=float(row["rate"]),fresh_s_rate=float(r["rate"]),
                relative_difference=delta/float(r["rate"]),combined_estimated_error=allowance,passed=passed))
        assert len(checks)>=2 and all(r["passed"] for r in checks),"Cache-lifecycle reference disagrees with completed native channels"
        runner.write_json(case/"cache_lifecycle_comparison_pilot_v1.json",dict(status="passed",checks=checks,
            interpretation="Difference checked against sum of interpolation estimators; not a proof of arbitrary cache equivalence.",
            original_summary="center_s7_v2_resource_stop.json",reference_summary="reference_center_s7_v1/summary.json"))
        while runner.free_gib()<6.5:
            state.update(status="waiting_for_memory",updated_utc=runner.now(),free_GiB=runner.free_gib())
            runner.write_json(path,state);time.sleep(20)
        state.update(status="running_side_references",cache_comparison="passed",updated_utc=runner.now())
        runner.write_json(path,state)
        jobs=[SimpleNamespace(name=name,stage="spectrum",profile="s7",points=point,
                    script="run_sigma_mott_precision_reference.jl") for name,point in
              (("reference_left_s7_v1","sig_m020"),("reference_right_s7_v1","sig_p020"))]
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures=[pool.submit(runner.job,case,j) for j in jobs]
            for j,future in zip(jobs,futures):
                future.result()
                state["completed"].append(j.name)
                runner.write_json(path,state)
        while True:
            pilot=runner.read_json(case/"reference_center_s7_v1.status.json")
            if pilot["status"]=="failed":
                raise RuntimeError("Center precision reference failed")
            if pilot["status"]=="completed_diagnostic":
                break
            time.sleep(20)
        new={r["channel"]:r for r in runner.read_csv(case/"reference_center_s7_v1"/"rates.csv")}
        final_checks=[]
        for row in old:
            r=new[row["channel"]]
            delta=float(r["rate"])-float(row["rate"])
            allowance=float(r["estimated_error"])+float(row["estimated_error"])+1e-12
            final_checks.append(dict(channel=row["channel"],native_rate=float(row["rate"]),fresh_s_rate=float(r["rate"]),
                relative_difference=delta/float(r["rate"]),combined_estimated_error=allowance,passed=abs(delta)<=3*allowance))
        assert all(r["passed"] for r in final_checks),"Final lifecycle checks failed"
        runner.write_json(case/"cache_lifecycle_comparison_v1.json",dict(status="passed",checks=final_checks,
            interpretation="Original retained five-channel pilot versus fresh-s reference; errors cover interpolation only."))
        # The analysis itself also runs remotely. Preserve exact stdout/stderr.
        env=os.environ.copy();env["SIGMA_REVIEW_SCRIPT"]=str(case/"review_sigma_mott_refinement_v3.py")
        with (case/"analysis_unit_checks_final.stdout.log").open("xb") as out,(case/"analysis_unit_checks_final.stderr.log").open("xb") as err:
            subprocess.run([sys.executable,"-m","unittest","-v","test_sigma_mott_analysis_v2"],cwd=case,
                env=env,stdout=out,stderr=err,check=True)
        with (case/"review_final_v1.stdout.log").open("xb") as out,(case/"review_final_v1.stderr.log").open("xb") as err:
            subprocess.run([sys.executable,str(case/"review_sigma_mott_refinement_v3.py"),str(case),
                str(case/"review_final_v1")],cwd=case,stdout=out,stderr=err,check=True)
        runner.verify_source(case)
        state.update(status="completed_diagnostic",finished_utc=runner.now(),review="review_final_v1/summary.json")
        runner.write_json(path,state)
    except BaseException as exc:
        state.update(status="failed",error=str(exc),finished_utc=runner.now())
        runner.write_json(path,state)
        raise


if __name__=="__main__":
    main()
