#!/usr/bin/env python3
"""Verify selected values and archive the augmented grid on quote-remote."""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import importlib.util
from pathlib import Path
import zipfile


def verify_and_package(case, source_root):
    spec = importlib.util.spec_from_file_location("assembler", case / "merge_transport_eta_refinement_v2.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    m.check_case_inputs(case)
    review = m.read_json(case / "grid/review.json")
    for relative, expected in review["output_sha256"].items():
        m.require(m.sha(case / "grid" / relative) == expected, "changed assembled output")
    selected = {r["point_id"]: r for r in m.read_csv(case / "selection/points.csv")}
    new_rates = m.read_csv(case / "selection/channels.csv")
    original_values_checked, original_rates_checked, selected_values_checked = 0, 0, 0
    counts = Counter()
    all_ids = []
    for mode in ("mode_a", "mode_b"):
        base = case / "inputs/base" / mode
        grid = case / "grid" / mode
        originals = {m.key(r): r for r in m.read_csv(base / "phase_guided_transport_scan.csv")}
        rows = m.read_csv(grid / "phase_guided_transport_grid.csv")
        provenance = [r for r in m.read_csv(case / "grid/point_provenance.csv") if r["mode"] == mode]
        m.require(len(rows) == len(provenance), "provenance count mismatch")
        for index, (row, origin) in enumerate(zip(rows, provenance), 2):
            all_ids.append(row["point_id"])
            counts[row["source_kind"]] += 1
            m.require(origin["point_id"] == row["point_id"] and int(origin["grid_csv_line"]) == index,
                      "provenance line does not identify the output point")
            if row["source_kind"] == "v14_retained":
                original = originals[m.key(row)]
                m.require(all(row[k] == v for k, v in original.items()), "changed retained original field")
                original_values_checked += 1
            else:
                point = selected[row["source_point_id"]]
                for field in (*m.OBSERVABLES, "T_fm", "muq_fm", "Phi", "Phibar", "m_u", "m_d", "m_s", "s_fm3inv",
                              *("n_" + q for q in m.SPECIES), *("tauinv_" + q for q in m.SPECIES)):
                    m.require(row[field] == point[field], "selected result was transformed")
                m.require(row["iterations"] == "", "unrecorded solver iteration count was invented")
                if row["source_kind"] == "xi_inserted":
                    m.require(row["phase_curr"] == "unknown", "new point was assigned an uncomputed phase label")
                selected_values_checked += 1
        rates = m.read_csv(grid / "phase_guided_transport_grid_sampling.csv")
        base_rates = {(m.rate_key(r), r["channel"]): r for r in m.read_csv(base / "phase_guided_transport_scan_sampling.csv")}
        point_map = {m.rate_key(r): r for r in rows}
        selection_rates = {(r["point_id"], r["channel"]): r for r in new_rates}
        m.require(len(rates) == 17 * len(rows), "wrong total channel count")
        for rate in rates:
            point = point_map[m.rate_key(rate)]
            if point["source_kind"] == "v14_retained":
                m.require(rate == base_rates[m.rate_key(rate), rate["channel"]], "changed retained channel data")
                original_rates_checked += 1
            else:
                chosen = selection_rates[point["source_point_id"], rate["channel"]]
                m.require(all(rate[k] == chosen[k] for k in rate if k in chosen), "changed selected channel data")
    m.require(dict(counts) == dict(v14_retained=1810, recomputed_original=9, xi_inserted=48), "wrong source counts")
    m.require(len(all_ids) == len(set(all_ids)) == 1867, "duplicate or missing point IDs")
    expected_sources = m.read_json(case / "inputs/refinement/expected_source_hashes.json")
    for relative, expected in expected_sources.items():
        m.require(m.sha(source_root / relative) == expected, "frozen numerical source changed during replay")
    seams = m.read_csv(case / "grid/boundary_secants.csv")
    summary = dict(schema="transport_eta_assembly_verification_v1", status="passed",
        host="quote-remote", generated_utc=datetime.now(timezone.utc).isoformat(),
        unchanged_original_rows_checked=original_values_checked, unchanged_original_channel_rows_checked=original_rates_checked,
        selected_rows_checked=selected_values_checked, point_counts=dict(counts), points=len(all_ids), channels=len(all_ids) * 17,
        grid_review_sha256=m.sha(case / "grid/review.json"), source_files_checked=len(expected_sources),
        all_source_hashes_match=True, script_sha256=m.sha(Path(__file__)),
        boundary_intervals_checked=len(seams), boundary_secants=seams,
        boundary_scope="Raw log secants with adjacent-interval context; no sign reversals. This is not a bound on derivative quadrature error.",
        gap_scope="Two original first-order coordinate pairs and phase labels preserved; no new phase classification")
    m.write_json(case / "assembly_verification.json", summary)
    print({k: v for k, v in summary.items() if k != "boundary_secants"}, flush=True)
    for row in seams:
        if row["observable"] == "sigma_over_T":
            print("BOUNDARY", row, flush=True)
    manifest_path, archive = case / "package_manifest.json", case / "eta_augmented_results_v1.zip"
    m.require(not manifest_path.exists() and not archive.exists(), "archive already exists")
    files = [dict(path=p.relative_to(case).as_posix(), bytes=p.stat().st_size, sha256=m.sha(p))
             for p in sorted(case.rglob("*")) if p.is_file() and "__pycache__" not in p.parts]
    m.write_json(manifest_path, dict(schema="transport_eta_augmented_package_v1", status="numerically_accepted_augmented_grid",
        calculation_sha=m.CALCULATION_SHA, points=1867, channels=31739, grid_review_sha256=m.sha(case / "grid/review.json"),
        source_files_checked=len(expected_sources), all_source_hashes_match=True, files=files))
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as out:
        for item in files:
            out.write(case / item["path"], item["path"])
        out.write(manifest_path, manifest_path.name)
    receipt = dict(archive=archive.name, bytes=archive.stat().st_size, sha256=m.sha(archive),
                   package_manifest_sha256=m.sha(manifest_path), manifest_files=len(files))
    m.write_json(case / "archive_receipt.json", receipt)
    print(receipt, flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("source_root", type=Path)
    args = parser.parse_args()
    verify_and_package(args.case.resolve(), args.source_root.resolve())
