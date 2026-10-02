#!/usr/bin/env python3
"""Losslessly consolidate v10/v11 manifest sidecars; never render or solve."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import tomllib
import zipfile

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.plotting.plot_bundle import build_bundle, expand_bundle, load_chart_records
from scripts.plotting.plot_manifest import input_record, sha256_file, write_manifest
from scripts.plotting.plot_provenance import is_code_record, read_manifest_record

BASE = ROOT / "data/outputs/figures/relaxtime/transport/phase_guided"
ANALYSIS = ROOT / "docs/analysis/relaxtime/phase_guided_transport"
DESTINATION = ANALYSIS / "publication_plot_manifest_migration_v1"
ARCHIVE = DESTINATION / "original_manifest_graph.zip"
REPORT = DESTINATION / "manifest.json"
CASES = (
    ("publication_clean_v10_png_review", "plot_manifest.json"),
    ("publication_clean_v11_png_review", "plot_manifest.json"),
    ("publication_clean_v11_pdf_review", "pdf_review_index.json"),
)


def relative(path: Path) -> str:
    return path.relative_to(ROOT).as_posix()


def archive_originals() -> dict:
    files = sorted(path for case, _ in CASES for path in (BASE / case).rglob("*.json"))
    if len(files) != 225:
        raise ValueError("expected exactly 225 original manifests/indexes before migration")
    original = {relative(path): path.read_bytes() for path in files}
    if ARCHIVE.exists():
        with zipfile.ZipFile(ARCHIVE) as archive:
            if set(archive.namelist()) != set(original):
                raise ValueError("existing archive has a different manifest set")
            if any(archive.read(name) != payload for name, payload in original.items()):
                raise ValueError("existing archive differs from original manifest bytes")
    else:
        ARCHIVE.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(ARCHIVE, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
            for name, payload in original.items():
                info = zipfile.ZipInfo(name, date_time=(2026, 10, 2, 0, 0, 0))
                info.compress_type = zipfile.ZIP_DEFLATED
                archive.writestr(info, payload, compresslevel=9)
        with zipfile.ZipFile(ARCHIVE) as archive:
            assert all(archive.read(name) == payload for name, payload in original.items())
    return input_record(ARCHIVE, role="immutable_original_manifest_graph")


def protected_records() -> list[dict]:
    paths = {path for case, _ in CASES for path in (BASE / case).rglob("*")
             if path.is_file() and path.suffix != ".json"}
    paths.add(ANALYSIS / "publication_clean_current.json")
    paths.add(ANALYSIS / "publication_clean_v11_stage_acceptance_v1.json")
    for version in (10, 11):
        directory = ANALYSIS / f"phase_guided_transport_publication_clean_v{version}_png_review"
        paths.update(path for path in directory.rglob("*") if path.is_file())
        package = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
        for record in package["inputs"]:
            path = ROOT / record["path"]
            if (path.is_file() and not is_code_record(path, ROOT)
                    and not (path.is_relative_to(BASE) and path.suffix == ".json")):
                paths.add(path)
    pdf_analysis = ANALYSIS / "phase_guided_transport_publication_clean_v11_pdf_review"
    paths.update(path for path in pdf_analysis.rglob("*") if path.is_file())
    return [input_record(path, role="unchanged_accepted_evidence") for path in sorted(paths)]


def apply_migration() -> dict:
    if REPORT.exists():
        raise FileExistsError("migration is immutable; use --check after applying it")
    archived = archive_originals()
    registry = tomllib.loads((ROOT / "config/plotting/historical_snapshots.toml").read_text(encoding="utf-8"))
    if not any(entry["archive"] == archived["path"] and entry["sha256"] == archived["sha256"]
               for entry in registry.get("manifest_archives", [])):
        raise ValueError("register the verified archive before applying migration")
    protected = protected_records()
    planned = []
    for case, old_name in CASES:
        source = BASE / case / old_name
        index, pairs = load_chart_records(source, root=ROOT)
        if len(pairs) != 74:
            raise ValueError(f"{case} must contain 74 charts")
        bundle = build_bundle(index, [record for _, record in pairs])
        bundle["manifest_migration"] = {
            "source_index": input_record(source, role="original_frozen_index"),
            "original_manifest_archive": archived,
            "transform": "lossless shared-field factoring; no rendering or numerical changes",
        }
        if [record for _, record in expand_bundle(bundle)] != [record for _, record in pairs]:
            raise ValueError("consolidation lost figure evidence")
        planned.append((source, pairs, bundle))
    results = []
    for source, pairs, bundle in planned:
        target = source.parent / "plot_manifest.json"
        sources = [ROOT / chart["manifest"] for chart, _ in pairs]
        if any(not path.resolve().is_relative_to(source.parent.resolve()) for path in sources):
            raise ValueError("sidecar removal escaped the intended case directory")
        old_bytes = source.stat().st_size + sum(path.stat().st_size for path in sources)
        write_manifest(target, bundle, overwrite=target == source)
        for path in sources:
            path.unlink()
        if source != target:
            source.unlink()
        results.append({**input_record(target, role="consolidated_figure_manifest"),
                        "source_index": bundle["manifest_migration"]["source_index"],
                        "figure_count": len(pairs), "original_manifest_bytes": old_bytes})
    report = {
        "schema": "publication_plot_manifest_migration_v1", "task_classification": "independent",
        "migration_kind": "metadata_storage_only", "source_code_commit": "0725c0c4f87bbcbced4b7641c0a086e65c6fe42b",
        "solver_called": False, "figures_rendered": False, "canonical_data_modified": False,
        "manuscript_eligibility_changed": False, "current_publication_layer_changed": False,
        "original_manifest_count": 225, "figure_record_count": 222,
        "original_manifest_archive": archived, "bundles": results, "protected": protected,
    }
    write_manifest(REPORT, report)
    check_migration()
    return report


def check_migration() -> dict:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    for record in [report["original_manifest_archive"], *report["bundles"], *report["protected"]]:
        path = ROOT / record["path"]
        if sha256_file(path) != record["sha256"] or path.stat().st_size != record["bytes"]:
            raise ValueError(f"migration integrity failure: {record['path']}")
    archive = ROOT / report["original_manifest_archive"]["path"]
    with zipfile.ZipFile(archive) as snapshot:
        if len(snapshot.namelist()) != report["original_manifest_count"]:
            raise ValueError("retired manifest inventory changed")
        for result in report["bundles"]:
            path = ROOT / result["path"]
            expected_assets = {record["path"] for record in report["protected"]
                               if (ROOT / record["path"]).is_relative_to(path.parent)}
            actual_assets = {relative(asset) for asset in path.parent.rglob("*")
                             if asset.is_file() and asset.suffix != ".json"}
            if actual_assets != expected_assets:
                raise ValueError("figure asset inventory changed")
            if sorted(p.name for p in path.parent.rglob("*.json")) != ["plot_manifest.json"]:
                raise ValueError("figure package still contains manifest sidecars")
            bundle = json.loads(path.read_text(encoding="utf-8"))
            original = read_manifest_record(result["source_index"], root=ROOT)
            pairs = expand_bundle(bundle)
            if len(pairs) != result["figure_count"]:
                raise ValueError("bundle figure count changed")
            if bundle["package"] != {key: value for key, value in original.items() if key != "charts"}:
                raise ValueError("index metadata was lost")
            for old_chart, (chart, record) in zip(original["charts"], pairs):
                expected = json.loads(snapshot.read(old_chart["manifest"]))
                if record != expected:
                    raise ValueError("consolidated figure record differs from original")
                if chart != {key: value for key, value in old_chart.items()
                             if key not in {"manifest", "manifest_sha256"}}:
                    raise ValueError("figure summary was lost")
                if hashlib.sha256(snapshot.read(old_chart["manifest"])).hexdigest() != old_chart["manifest_sha256"]:
                    raise ValueError("original manifest hash changed")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--archive-only", action="store_true")
    action.add_argument("--apply", action="store_true")
    action.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.archive_only:
        print(json.dumps(archive_originals(), indent=2))
        return
    report = apply_migration() if args.apply else check_migration()
    print(json.dumps({"status": "verified", "figures": report["figure_record_count"],
                      "protected_files": len(report["protected"]), "bundle_count": len(report["bundles"]),
                      "original_bytes": sum(item["original_manifest_bytes"] for item in report["bundles"]),
                      "consolidated_bytes": sum(item["bytes"] for item in report["bundles"])}, indent=2))


if __name__ == "__main__":
    main()
