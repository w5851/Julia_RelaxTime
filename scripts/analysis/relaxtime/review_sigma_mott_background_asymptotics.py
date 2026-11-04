#!/usr/bin/env python3
"""quote-remote only: reconcile smaller stencils, background grids and coefficient limits."""
from pathlib import Path
import csv
import importlib.util
import json
import math
import os
import subprocess
import sys
import zipfile


def main():
    case=Path(sys.argv[1])
    spec=importlib.util.spec_from_file_location("sigma_asym",case/"review_sigma_mott_asymptotics_v2.py")
    review=importlib.util.module_from_spec(spec);spec.loader.exec_module(review)
    out=case/"asym_background_review_v1";out.mkdir(exist_ok=False)
    original=case/"asym_final_review_v1"
    source=case/"asym_background_checks_v1"
    def checked_summary(folder):
        meta=json.loads((folder/"summary.json").read_text(encoding="utf-8-sig"))
        for name,expected in meta["files"].items():
            if review.digest(folder/name)!=expected:raise ValueError("Frozen file changed: "+name)
        return meta
    base_meta=checked_summary(original);meta=checked_summary(source)
    state=json.loads((case/"asym_background_checks_v1.status.json").read_text(encoding="utf-8-sig"))
    assert state["status"]=="completed_diagnostic" and state["exit_code"]==0
    assert state["source_files_checked_before"]==state["source_files_checked_after"]==340
    base=review.read_csv(original/"side_curvatures.csv")
    by_key={(r["side"],float(r["h"])):r for r in base}
    refs=review.read_csv(source/"background_reference.csv")
    variants=("same_background_small_step","fine_background")
    all_curves=[]
    for variant in variants:
        values=review.derive([r for r in refs if r["variant"]==variant],meta["steps"],fractions=(.125,.0625))
        all_curves.extend(dict(variant=variant,**row) for row in values)
    smaller={(r["side"],r["h"]):r for r in all_curves if r["variant"]==variants[0]}
    comparison=[]
    for row in all_curves:
        prev=by_key[row["side"],row["h"]]
        same=smaller[row["side"],row["h"]]
        comparison.append(dict(variant=row["variant"],side=row["side"],h=row["h"],
            curvature=row["curvature"],scaled_curvature=row["scaled_curvature"],
            change_from_original=row["curvature"]-float(prev["curvature"]),
            relative_change_from_original=row["curvature"]/float(prev["curvature"])-1,
            relative_change_from_same_step=row["curvature"]/same["curvature"]-1,
            quadrature_defect=row["quadrature_defect"],step_defect=row["step_defect"],
            effective_exponent=row["effective_exponent"]))
    # For a logarithmic right-side divergence, curvature gain per log refinement
    # tends to a finite coefficient. A power-law exponent is not the best test.
    log_rows=[]
    for variant,rows in [("original",base),*[(v,[r for r in all_curves if r["variant"]==v]) for v in variants]]:
        right=sorted([r for r in rows if r["side"]=="right"],key=lambda r:float(r["h"]),reverse=True)
        for previous,current in zip(right[:-1],right[1:]):
            gain=(float(current["curvature"])-float(previous["curvature"]))/math.log(float(previous["h"])/float(current["h"]))
            log_rows.append(dict(variant=variant,h=float(current["h"]),previous_h=float(previous["h"]),
                log_coefficient_estimate=gain,analytical_coefficient=meta["right_log_coefficient"],
                relative_difference=gain/meta["right_log_coefficient"]-1))
    def write_csv(name,rows):
        with (out/name).open("w",encoding="utf-8",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    write_csv("reference_curvatures.csv",all_curves)
    write_csv("background_step_comparison.csv",comparison)
    write_csv("right_log_coefficient.csv",log_rows)
    env=os.environ.copy();env["SIGMA_ASYM_REVIEW_SCRIPT"]=str(case/"review_sigma_mott_asymptotics_v2.py")
    tests=subprocess.run([sys.executable,str(case/"test_sigma_mott_asymptotics_v4.py")],env=env,capture_output=True,text=True)
    (out/"tests.stdout.log").write_text(tests.stdout,encoding="utf-8")
    (out/"tests.stderr.log").write_text(tests.stderr,encoding="utf-8")
    if tests.returncode:raise RuntimeError("Derivative review tests failed")
    summary=dict(status="completed_background_and_step_review",tests=dict(passed=True,count=7),
        execution=dict(exit_code=state["exit_code"],source_files_checked_before=340,source_files_checked_after=340),
        original_center=meta["original_center"],fine_center=meta["fine_center"],fine_center_shift=meta["fine_center_shift"],
        independent_left_coefficient=meta["left_sqrt_h_coefficient"],independent_right_log_coefficient=meta["right_log_coefficient"],
        scope=base_meta["scope"],comparison= comparison,right_log_comparison=log_rows,
        error_semantics="Empirical comparisons, not rigorous enclosure. Each grid is compared about its own Mott center. Smaller steps also include strict kernel and thermal quadrature.",
        files={p.name:review.digest(p) for p in out.iterdir() if p.is_file()})
    (out/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    paths=[p for folder in (source,out) for p in folder.rglob("*") if p.is_file()]
    paths.extend(case/name for name in ("check_sigma_mott_background_asymptotics_v1.jl",
        "review_sigma_mott_background_asymptotics.py","review_sigma_mott_asymptotics_v2.py","test_sigma_mott_asymptotics_v4.py",
        "asym_background_checks_v1.status.json","asym_background_checks_v1.stdout.log","asym_background_checks_v1.stderr.log"))
    manifest=case/"asym_background_return_manifest_v1.json"
    manifest.write_text(json.dumps(dict(files={p.relative_to(case).as_posix():dict(bytes=p.stat().st_size,sha256=review.digest(p))
        for p in paths}),indent=2)+"\n",encoding="utf-8")
    package=case/"packages/asymptotics_background_return_v1.zip"
    with zipfile.ZipFile(package,"x",zipfile.ZIP_DEFLATED) as archive:
        for p in paths:archive.write(p,p.relative_to(case).as_posix())
        archive.write(manifest,manifest.name)
    print(json.dumps(dict(status=summary["status"],package=str(package),bytes=package.stat().st_size,
        sha256=review.digest(package),files=len(paths),comparison=comparison,right_log_comparison=log_rows),allow_nan=False))


if __name__=="__main__":main()
