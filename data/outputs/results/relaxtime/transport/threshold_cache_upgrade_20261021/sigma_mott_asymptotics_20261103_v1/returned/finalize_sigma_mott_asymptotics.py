#!/usr/bin/env python3
"""quote-remote numerical review, execution checks and immutable evidence return."""
from pathlib import Path
import csv
import json
import math
import os
import subprocess
import sys
import zipfile
import review_sigma_mott_asymptotics as review


def write_csv(path,rows):
    with path.open("w",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def main():
    case=Path(sys.argv[1]);out=case/"asym_final_review_v1";out.mkdir(exist_ok=False)
    baseline=review.read_csv(case/"asym_full_v2/integrals.csv")
    refs=review.read_csv(case/"asym_checks_v1/reference_integrals.csv")
    source_meta=json.loads((case/"asym_full_v2/summary.json").read_text(encoding="utf-8-sig"))
    check_meta=json.loads((case/"asym_checks_v1/summary.json").read_text(encoding="utf-8-sig"))
    for folder in ("asym_full_v2","asym_checks_v1"):
        summary=json.loads((case/folder/"summary.json").read_text(encoding="utf-8-sig"))
        for name,expected in summary["files"].items():
            if review.digest(case/folder/name)!=expected:raise ValueError("Changed numerical input")
    main_rows=review.derive(baseline,source_meta["steps"])
    by_key={(r["side"],r["h"]):r for r in main_rows}
    comparisons=[]
    for variant in check_meta["variants"]:
        derived=review.derive([r for r in refs if r["variant"]==variant],check_meta["steps"])
        for row in derived:
            source=by_key[row["side"],row["h"]]
            comparisons.append(dict(variant=variant,side=row["side"],h=row["h"],
                reference_curvature=row["curvature"],baseline_curvature=source["curvature"],
                difference=row["curvature"]-source["curvature"],
                relative_difference=row["curvature"]/source["curvature"]-1,
                scaled_curvature=row["scaled_curvature"],
                step_defect=row["step_defect"],quadrature_defect=row["quadrature_defect"]))
    write_csv(out/"side_curvatures.csv",main_rows)
    write_csv(out/"precision_comparison.csv",comparisons)
    # Endpoint normal-form coefficients are independently estimated from numerical kernel
    # coefficients by a separate Julia probe when present; never used to generate rates.
    env=os.environ.copy();env["SIGMA_ASYM_REVIEW_SCRIPT"]=str(case/"review_sigma_mott_asymptotics.py")
    run=subprocess.run([sys.executable,str(case/"test_sigma_mott_asymptotics_v3.py")],env=env,capture_output=True,text=True)
    (out/"tests.stdout.log").write_text(run.stdout,encoding="utf-8")
    (out/"tests.stderr.log").write_text(run.stderr,encoding="utf-8")
    if run.returncode:raise RuntimeError("Derivative diagnostic unit tests failed")
    statuses={}
    for name in ("asym_full_v2","asym_checks_v1"):
        state=json.loads((case/(name+".status.json")).read_text(encoding="utf-8-sig"))
        if state["status"]!="completed_diagnostic" or state["source_files_checked_after"]!=340:
            raise ValueError("Incomplete or unverified remote job: "+name)
        statuses[name]=dict(exit_code=state["exit_code"],source_files_checked_before=state["source_files_checked_before"],
            source_files_checked_after=state["source_files_checked_after"])
    info=dict(status="completed_side_curvature_diagnostic",source=str(case/"asym_full_v2/summary.json"),
        source_sha256=review.digest(case/"asym_full_v2/summary.json"),scope=source_meta["scope"],xi_mott=source_meta["xi_mott"],
        tests=dict(passed=True,count=5),executions=statuses,
        derivative="Second xi derivative of continuous near-threshold mixed-scalar positive-square rate, not complete tau or terminal transport.",
        error_semantics="Numerical comparison envelope, not a certified bound. Loop/representation and continuous thermal-weight references at four distances on both sides.",
        guards="The guarded variant applies the existing determinant threshold exactly, including its discontinuous activation.",
        left_last3=[r for r in main_rows if r["side"]=="left"][-3:],
        right_last3=[r for r in main_rows if r["side"]=="right"][-3:],
        maximum_reference_relative_change={v:max(abs(r["relative_difference"]) for r in comparisons if r["variant"]==v)
            for v in check_meta["variants"]},
        files={p.name:review.digest(p) for p in out.iterdir() if p.is_file()})
    (out/"summary.json").write_text(json.dumps(info,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    folders=("asym_pilot_v1","asym_pilot_v2","asym_full_v1","asym_full_v2","asym_checks_v1","asym_final_review_v1")
    sources=("run_sigma_mott_asymptotics_v1.jl","run_sigma_mott_asymptotics_v2.jl","run_sigma_mott_asymptotics.jl",
        "run_sigma_mott_asymptotics_v3.jl","check_sigma_mott_asymptotics_v1.jl","review_sigma_mott_asymptotics.py",
        "finalize_sigma_mott_asymptotics.py","test_sigma_mott_asymptotics_v3.py","run_sigma_mott_thermal_audit.jl",
        "run_sigma_mott_spectrum_v2.jl","run_transport_other_mott_refinement_library_v2.jl",
        "run_transport_mott_probe.jl","transport_mott_thresholds.jl","transport_mott_kaon_endpoint.jl",
        "run_sigma_mott_refinement_v2.py","expected_source_hashes.json","contract.json")
    paths=[case/name for name in sources]
    for folder in folders:
        paths.extend(p for p in (case/folder).rglob("*") if p.is_file())
        for suffix in (".status.json",".stdout.log",".stderr.log"):
            if (case/(folder+suffix)).is_file():paths.append(case/(folder+suffix))
    paths.extend(p for p in (case/"inputs").rglob("*") if p.is_file())
    files={p.relative_to(case).as_posix():dict(bytes=p.stat().st_size,sha256=review.digest(p)) for p in paths}
    manifest=case/"asym_return_manifest_v1.json"
    manifest.write_text(json.dumps(dict(files=files),indent=2)+"\n",encoding="utf-8")
    package=case/"packages/asymptotics_return_v1.zip"
    with zipfile.ZipFile(package,"x",zipfile.ZIP_DEFLATED) as archive:
        for p in paths:archive.write(p,p.relative_to(case).as_posix())
        archive.write(manifest,manifest.name)
    print(json.dumps(dict(status=info["status"],package=str(package),bytes=package.stat().st_size,
        sha256=review.digest(package),files=len(files),reference_changes=info["maximum_reference_relative_change"])))


if __name__=="__main__":main()
