#!/usr/bin/env python3
"""Freeze a compact completed quote-remote sigma evidence package, without promotion."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import zipfile


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case",type=Path)
    parser.add_argument("name")
    parser.add_argument("--review",default="review_uniform_with_cutoff_v1")
    args=parser.parse_args()
    case=args.case.resolve()
    if not args.name.replace("_","").isalnum():
        raise ValueError("Invalid package name")
    package=case/"packages"/args.name
    package.mkdir(parents=True,exist_ok=False)
    complete=tuple(dict.fromkeys(("left_s6_v1","right_s6_v1","spectrum_v2","thermal_audit_v1","weight_checks_v1",
              "reference_center_s7_v1","reference_left_s7_v1","reference_right_s7_v1","cutoff_audit_v1",
              "review_final_v1","review_with_cutoff_v1","review_uniform_v1",args.review)))
    for name in complete:
        summary=read_json(case/name/"summary.json")
        assert summary["status"].startswith("completed"),(name,summary.get("status"))
        for file,digest in summary["files"].items():
            assert sha(case/name/file)==digest,(name,file)
        if not name.startswith("review_"):
            status=read_json(case/(name+".status.json"))
            assert status["status"]=="completed_diagnostic" and status["exit_code"]==0
            assert status["source_files_checked_before"]==status["source_files_checked_after"]==340
            assert status["source_all_match"] is True
    assert read_json(case/"cache_lifecycle_comparison_v1.json")["status"]=="passed"
    assert read_json(case/"final_acceptance_analysis_v1.json")["status"]=="passed"
    assert read_json(case/args.review/"cutoff_event_validation.json")["status"]=="passed"
    contract=read_json(case/"contract.json")
    for name,digest in contract["input_sha256"].items():
        assert sha(case/name)==digest,("frozen input",name)
    candidates=[]
    for path in case.iterdir():
        if path.is_file() and path.suffix in (".jl",".py",".ps1",".json",".log"):
            candidates.append(path)
    for name in (*complete,"inputs"):
        candidates.extend(p for p in (case/name).rglob("*") if p.is_file() and "__pycache__" not in p.parts)
    candidates=sorted(set(candidates),key=lambda p:p.relative_to(case).as_posix())
    entries=[dict(path=p.relative_to(case).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in candidates]
    manifest=dict(schema="sigma_mott_diagnostic_package_v1",generated_utc=datetime.now(timezone.utc).isoformat(),
        source_host="quote-remote",remote_case=str(case),completed_runs=list(complete),
        authoritative_review=args.review,
        historical_review_scope="review_final_v1 and review_with_cutoff_v1 preserve the initial mixed-tier hybrid controls. Use the uniform review for cross-window comparisons.",
        calculation_sha="2e5c411c716d5cbc524c3fc6cfb335fe05b39444",files=entries,
        preserved_failed_runs=["spectrum_v1","center_s7_v1","center_s7_v2"],
        failed_run_scope="Root logs/status retained. Interrupted center_s7_v2 completed rates retained separately.",
        production_kernel_changed=False,manuscript_eligible=False,author_accepted=False)
    # Interrupted native rates are important independent reference evidence.
    interrupted=case/"center_s7_v2/rates.csv"
    if interrupted.exists():
        candidates.append(interrupted)
        manifest["files"].append(dict(path=interrupted.relative_to(case).as_posix(),bytes=interrupted.stat().st_size,sha256=sha(interrupted)))
    (package/"package_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    archive=package/(args.name+".zip")
    with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,compresslevel=6) as bundle:
        for path in candidates:
            bundle.write(path,path.relative_to(case).as_posix())
        bundle.write(package/"package_manifest.json","package_manifest.json")
    receipt=dict(status="packaged_completed_diagnostic",archive=archive.name,bytes=archive.stat().st_size,
        archive_sha256=sha(archive),manifest_sha256=sha(package/"package_manifest.json"),files=len(candidates))
    (package/"archive_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(receipt),flush=True)


if __name__=="__main__":main()
