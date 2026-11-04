#!/usr/bin/env python3
"""Verify and package the completed eta diagnostic on quote-remote."""
import hashlib
import json
from pathlib import Path
import sys
import zipfile
from datetime import datetime, timezone


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def main():
    case=Path(sys.argv[1]).resolve()
    source=Path(sys.argv[2]).resolve()
    archive=case/"eta_refinement_results_v1.zip"
    if archive.exists() or (case/"package_manifest.json").exists():
        raise FileExistsError("retained package already exists")
    for stage in ("locate_v3.status.json","transport_queue_v2.status.json","peak_refinement_v1.status.json"):
        status=read(case/stage)
        assert status["status"]=="finished",stage
        assert all(r["exit_code"]==0 for r in status["workers"]),stage
    expected=read(case/"expected_source_hashes.json")
    assert len(expected)==340
    assert all(sha(source/name)==digest for name,digest in expected.items())
    config=read(case/"inputs_v2/config.json")
    assert sha(case/"inputs_v2/source_rows.csv")==config["source_rows_sha256"]
    for name,profile in config["profiles"].items():
        assert sha(case/"inputs_v2"/profile["path"])==profile["sha256"],name
    review=read(case/"review_final_v1/review.json")
    assert review["status"]=="completed_diagnostic_review" and review["points"]==51
    assert review["channels"]==51*17
    assert review["comparison_points"]>=9
    for file,digest in review["output_sha256"].items():
        assert sha(case/"review_final_v1"/file)==digest,file
    receipt=dict(status="completed_diagnostic",calculation_sha=config["calculation_sha"],
        source_files_checked=340,all_source_hashes_match=True,points=review["points"],
        channels=review["channels"],comparison_points=review["comparison_points"],
        review_sha256=sha(case/"review_final_v1/review.json"),
        validation_scope="completeness, positivity, selected xi/precision comparisons, immutable source identity; no amplitude ablation or production promotion",
        timestamp_utc=datetime.now(timezone.utc).isoformat())
    (case/"validation.json").write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    excluded={"package_manifest.json","archive_receipt.json","finish_eta_v1.status.json"}
    files=[p for p in sorted(case.rglob("*")) if p.is_file() and p.suffix.lower()!=".zip"
           and p.name not in excluded and "__pycache__" not in p.parts]
    manifest=dict(schema="transport_eta_refinement_package_v1",**receipt,
        files=[dict(path=p.relative_to(case).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in files])
    manifest_path=case/"package_manifest.json"
    manifest_path.write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    with zipfile.ZipFile(archive,"x",compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in [*files,manifest_path]:
            z.write(p,p.relative_to(case).as_posix())
    transfer=dict(archive=archive.name,bytes=archive.stat().st_size,sha256=sha(archive),
        package_manifest_sha256=sha(manifest_path),manifest_files=len(files))
    (case/"archive_receipt.json").write_text(json.dumps(transfer,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(transfer,indent=2))


if __name__=="__main__":
    main()
