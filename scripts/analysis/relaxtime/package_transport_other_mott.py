#!/usr/bin/env python3
"""Package completed diagnostic phases without copying live numerical outputs."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import zipfile


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("source", type=Path)
    parser.add_argument("name")
    parser.add_argument("--runs", required=True)
    parser.add_argument("--reviews", default="survey_review_v2")
    args = parser.parse_args()
    case, source = args.case.resolve(), args.source.resolve()
    runs = args.runs.split(",")
    for name in runs:
        status = json.loads((case / (name+".status.json")).read_text(encoding="utf-8-sig"))
        assert status["status"] == "completed_diagnostic", (name, status)
        assert status["source_all_match"] and status["source_files_checked_after"] == 340
    expected = json.loads((case / "expected_source_hashes.json").read_text(encoding="utf-8-sig"))
    assert len(expected) == 340 and all(sha(source / name) == h for name, h in expected.items())
    files = {p for p in case.iterdir() if p.is_file() and p.suffix.lower() in (".jl", ".py", ".ps1")}
    files.add(case / "expected_source_hashes.json")
    for name in ["inputs", *runs, *args.reviews.split(",")]:
        folder = case / name
        assert folder.resolve().is_relative_to(case) and folder.is_dir()
        files.update(p for p in folder.rglob("*") if p.is_file() and "__pycache__" not in p.parts)
    for name in runs:
        files.update(p for p in case.glob(name+".*") if p.is_file())
    for name in ("slope_unit_v1.log", "slope_unit_v2.log", "slope_unit_v3.stdout.log", "slope_unit_v3.stderr.log", "pipeline_syntax_v1.log",
                 "slope_unit_v4.stdout.log", "slope_unit_v4.stderr.log",
                 "pipeline_v1.config.json", "pipeline_v2.config.json", "pipeline_v1_superseded.json",
                 "fine_queue_plan_v1.json", "location_review_v1.json",
                 "review_survey_v1.stdout.log", "review_survey_v2.stdout.log"):
        if (case / name).is_file():
            files.add(case / name)
    output = case / "packages" / args.name
    assert output.resolve().is_relative_to(case)
    output.mkdir(parents=True, exist_ok=False)
    validation = dict(status="completed_diagnostic_phases", completed_runs=runs,
        calculation_sha="2e5c411c716d5cbc524c3fc6cfb335fe05b39444", source_files_checked=340,
        all_source_hashes_match=True, whole_task_complete=False, production_data_changed=False,
        timestamp_utc=datetime.now(timezone.utc).isoformat())
    validation_path = output / "phase_validation.json"
    validation_path.write_text(json.dumps(validation, indent=2) + "\n", encoding="utf-8")
    entries = [(p, p.relative_to(case).as_posix()) for p in sorted(files)] + [(validation_path, validation_path.name)]
    manifest = dict(schema="transport_other_mott_phase_package_v1", validation=validation,
        files=[dict(path=name, bytes=p.stat().st_size, sha256=sha(p)) for p, name in entries])
    manifest_path = output / "package_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    archive = output / (args.name + ".zip")
    with zipfile.ZipFile(archive, "x", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as bundle:
        for path, name in [*entries, (manifest_path, manifest_path.name)]:
            bundle.write(path, name)
    receipt = dict(archive=archive.name, bytes=archive.stat().st_size, sha256=sha(archive),
                   manifest_files=len(entries), package_manifest_sha256=sha(manifest_path))
    (output / "archive_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt))


if __name__ == "__main__":
    main()
