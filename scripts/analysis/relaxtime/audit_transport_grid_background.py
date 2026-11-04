"""Compare reviewed transport backgrounds and catalogue adjacent-grid changes.

Run on quote-remote. The background target is the existing 1 percent mass/entropy
qualification contract. The 50 percent neighbour screen only selects diagnostics;
it is neither an integration-error estimate nor a tolerance adjustment.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import math
from pathlib import Path
import sys
import zipfile

from review_transport_full_grid import digest, point_key, read, read_csv, require, write_csv

OBSERVABLES = ("tau_u", "tau_d", "tau_s", "tau_ubar", "tau_dbar", "tau_sbar",
               "eta_over_s", "zeta_over_s", "sigma_over_T")


def load_background(archive, checksum, mode):
    require(digest(archive) == checksum, "background archive differs from retained evidence")
    with zipfile.ZipFile(archive) as bundle:
        names = bundle.namelist()
        csv_names = [name for name in names if name.endswith("/background.csv")]
        review_names = [name for name in names if name.endswith("reviewed_result.json")]
        require(len(csv_names) == len(review_names) == 1, "ambiguous background evidence")
        reviewed = json.loads(bundle.read(review_names[0]).decode("utf-8-sig"))
        require(reviewed["phase_plan_passed" if mode == "A" else "mode_b_background_qualified"],
                "background grid is not qualified")
        content = bundle.read(csv_names[0])
        rows = list(csv.DictReader(io.StringIO(content.decode("utf-8-sig"))))
    mapped = {point_key(row, mode): row for row in rows}
    require(len(mapped) == len(rows) == reviewed["points"] and all(row["passed"] == "true" for row in rows),
            "incomplete or unqualified background points")
    return mapped, dict(archive_sha256=checksum, background_csv_sha256=hashlib.sha256(content).hexdigest(),
                       prior_review=reviewed)


def compare_background(row, reference, mode):
    require(abs(float(row["T_MeV"]) - float(reference["T_MeV"])) <= 1e-5,
            "background temperature differs from the qualified grid")
    changes = {name: abs(float(row[name]) / float(reference[name]) - 1) for name in ("m_u", "m_s")}
    if "m_d" in reference:
        changes["m_d"] = abs(float(row["m_d"]) / float(reference["m_d"]) - 1)
    changes["entropy"] = abs(float(row["s_fm3inv"]) / float(reference["entropy"]) - 1)
    distinct = reference["distinct_branches"] == "true"
    same_branch = not distinct or row["phase_curr"] == reference["stable_phase"]
    passed = all(math.isfinite(value) and value <= 0.01 for value in changes.values()) and same_branch
    return dict(mode=mode, T_MeV=row["T_MeV"], muB_MeV=row["muB_MeV"], xi=row["xi"],
                max_mass_relative_change=max(value for key, value in changes.items() if key != "entropy"),
                entropy_relative_change=changes["entropy"], distinct_reference_branches=distinct,
                expected_phase=reference["stable_phase"], actual_phase=row["phase_curr"],
                same_distinct_branch=same_branch, passed=passed)


def neighbouring_changes(rows, backgrounds, expected):
    curves = defaultdict(list)
    expected_curves = defaultdict(list)
    for key in expected:
        expected_curves[key[:-1]].append(key[-1])
    for curve in expected_curves:
        expected_curves[curve].sort()
    for key, row in rows.items():
        curves[key[:-1]].append((key, row))
    changes = []
    for curve, points in curves.items():
        points.sort(key=lambda item: item[0][-1])
        order = {xi: index for index, xi in enumerate(expected_curves[curve])}
        for (left_key, left), (right_key, right) in zip(points, points[1:]):
            if order[right_key[-1]] - order[left_key[-1]] != 1:
                continue
            a, b = backgrounds[left_key], backgrounds[right_key]
            boundary = (a["distinct_branches"] == b["distinct_branches"] == "true" and
                        a["stable_phase"] != b["stable_phase"])
            for observable in OBSERVABLES:
                x, y = float(left[observable]), float(right[observable])
                change = abs(x - y) / max(abs(x), abs(y), 1e-300)
                changes.append(dict(mode=curve[0], muB_MeV=curve[1], curve_parameter=curve[2],
                    xi_left=left_key[-1], xi_right=right_key[-1], observable=observable, left=x, right=y,
                    relative_to_larger=change, reference_first_order_boundary=boundary,
                    mass_change=abs(float(right["m_u"]) / float(left["m_u"]) - 1),
                    entropy_change=abs(float(right["s_fm3inv"]) / float(left["s_fm3inv"]) - 1),
                    selected_for_review=change > 0.5))
    return changes


def audit(case, grid_review, background_inputs, output):
    require(not output.exists(), "use a new audit output directory")
    report = read(grid_review / "reviewed_grid.json")
    require(not report["errors"], "resolve shard artifact errors before background audit")
    specification = read(case / "dispatch_manifest.json")
    require(report["input_metadata_sha256"]["queue_plan.json"] == digest(case / "queue_plan.json"),
            "grid review belongs to another case")
    inputs = read(background_inputs / "background_inputs.json")
    backgrounds, origins = {}, {}
    for mode, item in inputs["modes"].items():
        mapped, identity = load_background(background_inputs / item["archive"], item["sha256"], mode)
        backgrounds.update(mapped)
        origins[mode] = identity
    entries = {entry["run_id"]: entry for entry in specification["runs"]}
    entries.update({entry["run_id"]: dict(entry, directory="pilot_runs/" + str(entry["run_id"]))
                    for entry in specification["reused_pilot_runs"]})
    actual, comparisons = {}, []
    for result in report["runs"]:
        entry = entries[result["run_id"]]
        paths = list((case / entry["directory"] / "artifact").rglob("phase_guided_transport_scan.csv"))
        require(len(paths) == 1, "ambiguous reviewed CSV")
        for row in read_csv(paths[0])[1]:
            key = point_key(row, entry["mode"])
            require(key in backgrounds and key not in actual, "background point missing or duplicated")
            actual[key] = row
            comparisons.append(compare_background(row, backgrounds[key], entry["mode"]))
    require(len(actual) == report["completed_points"], "reviewed point count changed")
    changes = neighbouring_changes(actual, backgrounds, backgrounds)
    screened = [row for row in changes if row["selected_for_review"]]
    output.mkdir(parents=True)
    write_csv(output / "background_comparisons.csv", list(comparisons[0]), comparisons)
    if changes:
        write_csv(output / "adjacent_changes.csv", list(changes[0]), changes)
    if screened:
        write_csv(output / "screened_changes.csv", list(screened[0]), screened)
    summary = dict(host="quote-remote", created_utc=datetime.now(timezone.utc).isoformat(),
        points=len(actual), full_grid=len(actual) == len(backgrounds),
        background_target=0.01, all_backgrounds_passed=all(row["passed"] for row in comparisons),
        maximum_mass_relative_change=max(row["max_mass_relative_change"] for row in comparisons),
        maximum_entropy_relative_change=max(row["entropy_relative_change"] for row in comparisons),
        branch_mismatches=[row for row in comparisons if not row["same_distinct_branch"]],
        screened_observable_pairs=len(screened),
        screened_nonboundary_pairs=[row for row in screened if not row["reference_first_order_boundary"]],
        known_boundary_pairs=[row for row in screened if row["reference_first_order_boundary"]],
        top_nonboundary_changes=sorted((row for row in changes if not row["reference_first_order_boundary"]),
                                     key=lambda row: row["relative_to_larger"], reverse=True)[:20],
        quality_flags_retained=len(report["quality_flags"]),
        grid_review_sha256=digest(grid_review / "reviewed_grid.json"),
        background_inputs_sha256=digest(background_inputs / "background_inputs.json"),
        background_sources=origins, reviewer_sha256=digest(Path(__file__)),
        files={p.name: digest(p) for p in output.iterdir() if p.is_file()},
        scope="Background replay and adjacent-grid screening; representative convergence evidence remains separate")
    (output / "background_audit.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in summary.items() if key not in
          ("background_sources", "screened_nonboundary_pairs", "known_boundary_pairs", "top_nonboundary_changes")}))
    return 0 if summary["all_backgrounds_passed"] else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("case", "grid_review", "background_inputs", "output"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    sys.exit(audit(args.case.resolve(), args.grid_review.resolve(), args.background_inputs.resolve(), args.output.resolve()))
