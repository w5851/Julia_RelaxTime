#!/usr/bin/env python3
"""Assemble an eta-refined transport grid from frozen results on quote-remote."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from decimal import Decimal
from pathlib import Path
import shutil
import sys

CALCULATION_SHA = "2e5c411c716d5cbc524c3fc6cfb335fe05b39444"
CURVES = ("a_mu0", "a_mu450", "b_mu0")
SPECIES = ("u", "d", "s", "ubar", "dbar", "sbar")
OBSERVABLES = tuple("tau_" + s for s in SPECIES) + ("eta", "zeta", "sigma", "eta_over_s", "zeta_over_s", "sigma_over_T")
FINE_RUNS = ("fine_a_mu0_v1", "fine_a_mu450_v1", "fine_a_mu0_v2", "fine_a_mu450_v2", "fine_b_mu0_v2")
TIGHT_RUNS = tuple("tighter_" + c + "_v2" for c in CURVES) + tuple("tighter_peak_" + c + "_v1" for c in CURVES)
META_FIELDS = ("mode", "phase_reference_kind", "scan_group", "group_label", "plot_panel", "plot_panel_label",
               "plot_series", "plot_series_label", "T_phase_base_MeV", "alpha_T")
PHASE_FIELDS = ("phase_prev", "phase_curr", "phase_structure", "phase_boundary_xi_used")
EXTRA_FIELDS = ("point_id", "source_kind", "sampling_profile", "source_point_id", "source_run", "source_csv_line")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def write_csv(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields or list(rows[0]), extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def key(row):
    return row["plot_panel"], row["plot_series"], Decimal(row["xi"])


def rate_key(row):
    # Sampling CSV temperatures occasionally differ from the main CSV by one ulp.
    return round(float(row["T_MeV"]), 8), round(float(row["muB_MeV"]), 8), Decimal(row["xi"])


def check_case_inputs(case):
    manifest = read_json(case / "input_manifest.json")
    require(manifest["calculation_sha"] == CALCULATION_SHA, "unexpected calculation source")
    for item in manifest["files"]:
        path = case / item["path"]
        require(path.stat().st_size == item["bytes"] and sha(path) == item["sha256"], f"changed input: {item['path']}")
    base = case / "inputs/base"
    acceptance = read_json(base / "acceptance.json")
    promotion = read_json(base / "promotion.json")
    require(acceptance["full_grid_numerical_acceptance_passed"] is True and acceptance["points"] == 1819,
            "base grid has not passed numerical acceptance")
    require(promotion["formal_production_promoted"] is True and promotion["calculation_sha"] == CALCULATION_SHA,
            "base grid has not been promoted")
    for relative, expected in acceptance["merged_files"].items():
        require(sha(base / relative) == expected, "base CSV differs from its acceptance")
    require(sha(base / "background_audit.json") == acceptance["background_audit_sha256"], "base background audit changed")
    require(sha(base / "reviewed_grid.json") == acceptance["grid_review_sha256"], "base grid review changed")
    return manifest


def select_results(fine, tighter):
    """Prefer the tighter result only for an identical point and coordinate."""
    selected = {}
    for rows, expected_variant in ((fine, "fine"), (tighter, "tighter")):
        seen = set()
        for row in rows:
            identity = row["point_id"]
            require(identity not in seen, "duplicate result within one precision level")
            seen.add(identity)
            require(row["variant"] == expected_variant, "precision/source mismatch")
            if identity in selected:
                previous = selected[identity]
                require(all(row[k] == previous[k] for k in ("curve_id", "xi", "T_MeV", "muB_MeV", "original_point_id")),
                        "tighter point does not match its fine coordinate")
            selected[identity] = row
    points = list(selected.values())
    require(len({(p["curve_id"], Decimal(p["xi"])) for p in points}) == len(points), "duplicate physical coordinate")
    return sorted(points, key=lambda p: (p["curve_id"], Decimal(p["xi"])))


def prepare(case, eta_case):
    check_case_inputs(case)
    output = case / "selection"
    require(not output.exists(), "selection already exists")
    sys.path.insert(0, str(eta_case))
    import review_transport_eta_refinement as review
    package = read_json(eta_case / "package_manifest.json")
    require(package["calculation_sha"] == CALCULATION_SHA and package["all_source_hashes_match"] is True,
            "unverified refinement source")
    require(sha(eta_case / "package_manifest.json") == read_json(case / "inputs/refinement/local_return_verification.json")["archive"]["package_manifest_sha256"],
            "refinement package differs from the verified return")
    for item in package["files"]:
        path = eta_case / item["path"]
        require(path.stat().st_size == item["bytes"] and sha(path) == item["sha256"], f"changed refinement input: {item['path']}")
    report = read_json(eta_case / "review_final_v1/review.json")
    require(report["points"] == 51 and report["comparison_points"] == 18, "incomplete refinement review")
    require(sha(eta_case / "review_final_v1/review.json") == package["review_sha256"], "refinement review changed")
    for name, digest in report["output_sha256"].items():
        require(sha(eta_case / "review_final_v1" / name) == digest, "changed refinement review")
    levels = {}
    all_rates = {}
    source_records = []
    for variant, names in (("fine", FINE_RUNS), ("tighter", TIGHT_RUNS)):
        values = []
        for name in names:
            rows, rates, _ = review.collect_runs(eta_case, [name], False)
            actual = read_csv(eta_case / name / "observables.csv")
            lines = {r["point_id"]: i for i, r in enumerate(actual, 2)}
            summary_name = "summary.json" if (eta_case / name / "summary.json").exists() else "partial_summary.json"
            for filename in ("observables.csv", "rates.csv", summary_name):
                source_records.append(dict(path=f"{name}/{filename}", sha256=sha(eta_case / name / filename)))
            for row in rows:
                values.append(dict(row, variant=variant, source_run=name,
                                   source_csv_line=str(lines[row["point_id"]]),
                                   source_observables_sha256=sha(eta_case / name / "observables.csv"),
                                   source_rates_sha256=sha(eta_case / name / "rates.csv")))
            all_rates[name] = rates
        levels[variant] = values
    points = select_results(levels["fine"], levels["tighter"])
    require(len(points) == 57, "expected 57 distinct selected points")
    require(sum(p["variant"] == "tighter" for p in points) == 18, "expected 18 tighter and 39 fine points")
    for curve in CURVES:
        part = [p for p in points if p["curve_id"] == curve]
        require(len(part) == 19 and sum(bool(p["original_point_id"]) for p in part) == 3,
                "expected 16 new points and 3 original-coordinate replays per curve")
    rates = []
    for row in points:
        selected = [dict(r, variant=row["variant"], source_run=row["source_run"]) for r in all_rates[row["source_run"]]
                    if r["point_id"] == row["point_id"]]
        require(len(selected) == len({r["channel"] for r in selected}) == 17, "incomplete channel selection")
        rates.extend(selected)
    output.mkdir()
    write_csv(output / "points.csv", points)
    write_csv(output / "channels.csv", rates)
    write_json(output / "selection.json", dict(status="selected_frozen_results", points=57, added=48, recomputed=9,
        fine_points=sum(p["variant"] == "fine" for p in points), tighter_points=sum(p["variant"] == "tighter" for p in points),
        calculation_sha=CALCULATION_SHA, eta_case=str(eta_case), source_records=source_records,
        input_review_sha256=sha(eta_case / "review_final_v1/review.json"),
        points_sha256=sha(output / "points.csv"), channels_sha256=sha(output / "channels.csv"), script_sha256=sha(Path(__file__))))
    print(json.dumps(dict(status="selected_frozen_results", points=57, added=48, recomputed=9)))


def assembled_row(point, original, template, thermo, mode):
    """Do not invent missing historical solver-iteration or phase-tracker fields."""
    result = {name: "" for name in template}
    result.update({name: template[name] for name in META_FIELDS})
    for name in ("T_MeV", "muq_MeV", "muB_MeV"):
        result[name] = template[name]
    result.update({name: point[name] for name in point if name in result})
    result["T_MeV"], result["muq_MeV"], result["muB_MeV"] = (template[name] for name in ("T_MeV", "muq_MeV", "muB_MeV"))
    result.update({name: thermo[name] for name in thermo if name in result})
    # These source values are already frozen; the thermodynamic replay checks them.
    for name in (*OBSERVABLES, "T_fm", "muq_fm", "Phi", "Phibar", "m_u", "m_d", "m_s", "s_fm3inv",
                 *("n_" + s for s in SPECIES), *("tauinv_" + s for s in SPECIES)):
        result[name] = point[name]
    result.update(converged="true", iterations="", residual_norm=point["background_residual"],
                  equilibrium_backend="models_multi", seed_source="reused_eta_refinement_stable_multiseed",
                  run_id=f"mott_eta_fine_20261029_v1/{point['source_run']}")
    if original:
        result.update({name: original[name] for name in PHASE_FIELDS})
    else:
        result.update(phase_prev="unknown", phase_curr="unknown", phase_structure="unknown", phase_boundary_xi_used="")
    result.update(point_id=point["original_point_id"] or "eta_" + point["point_id"],
                  source_kind="recomputed_original" if original else "xi_inserted",
                  sampling_profile="eta_xi_" + point["variant"] + "_v1", source_point_id=point["point_id"],
                  source_run=point["source_run"], source_csv_line=point["source_csv_line"])
    require(result["mode"] == mode, "mode mapping changed")
    return result


def slopes(points, observable):
    return [dict(left=a["point_id"], right=b["point_id"], xi_left=float(a["xi"]), xi_right=float(b["xi"]),
                 slope=math.log(float(b[observable]) / float(a[observable])) / (float(b["xi"]) - float(a["xi"])))
            for a, b in zip(points, points[1:])]


def merge(case, eta_case):
    check_case_inputs(case)
    output = case / "grid"
    require(not output.exists(), "grid already exists")
    selection = read_json(case / "selection/selection.json")
    require(sha(case / "selection/points.csv") == selection["points_sha256"], "changed selected points")
    require(sha(case / "selection/channels.csv") == selection["channels_sha256"], "changed selected channels")
    thermal_review = read_json(case / "thermo/summary.json")
    require(thermal_review["status"] == "passed" and thermal_review["points"] == 57, "thermodynamic replay failed")
    require(sha(case / "thermo/rows.csv") == thermal_review["rows_sha256"], "changed thermodynamic replay")
    points = read_csv(case / "selection/points.csv")
    thermal = {r["point_id"]: r for r in read_csv(case / "thermo/rows.csv")}
    require(len(thermal) == len(points) == 57 and set(thermal) == {r["point_id"] for r in points}, "thermo point identities differ")
    for point in points:
        require(all(Decimal(point[k]) == Decimal(thermal[point["point_id"]][k]) for k in ("T_MeV", "muB_MeV", "xi")),
                "thermo coordinate mismatch")
    enriched = {r["point_id"]: r for r in read_csv(eta_case / "inputs_v2/source_rows.csv")}
    by_mode = defaultdict(list)
    for row in enriched.values():
        by_mode[row["mode_key"]].append(row)
    channels = read_csv(case / "selection/channels.csv")
    counts, evidence, changes, seam_rows, provenance = {}, {}, [], [], []
    output.mkdir()
    for mode, old_count, added in (("mode_a", 910, 32), ("mode_b", 909, 16)):
        base = case / "inputs/base" / mode
        originals = read_csv(base / "phase_guided_transport_scan.csv")
        fields = list(originals[0])
        require(len(originals) == old_count, "wrong base point count")
        original_by_key = {key(r): r for r in originals}
        enriched_keys = {key(r): r for r in by_mode[mode]}
        require(len(original_by_key) == len(originals) and len(enriched_keys) == len(by_mode[mode]), "duplicate base coordinates")
        require(set(original_by_key) == set(enriched_keys), "eta reference and base coordinates differ")
        for identity, row in original_by_key.items():
            require(all(row[n] == enriched_keys[identity][n] for n in fields), "eta reference differs from v14 base row")
        mapped = {identity: dict(row, point_id=enriched_keys[identity]["point_id"], source_kind="v14_retained",
                                 sampling_profile="thermal_weighted_v1", source_point_id=enriched_keys[identity]["point_id"],
                                 source_run=row["run_id"], source_csv_line=enriched_keys[identity]["source_csv_line"])
                  for identity, row in original_by_key.items()}
        updates = [p for p in points if ("mode_a" if p["curve_id"].startswith("a_") else "mode_b") == mode]
        affected = set()
        for point in updates:
            candidates = [r for r in originals if math.isclose(float(r["T_MeV"]), float(point["T_MeV"]), rel_tol=0, abs_tol=1e-10)
                          and float(r["muB_MeV"]) == float(point["muB_MeV"])]
            require(candidates, "selected point has no parent curve")
            template = candidates[0]
            identity = (template["plot_panel"], template["plot_series"], Decimal(point["xi"]))
            previous = original_by_key.get(identity)
            require(bool(previous) == bool(point["original_point_id"]), "insert/replay identity mismatch")
            if previous:
                require(enriched_keys[identity]["point_id"] == point["original_point_id"], "wrong original-coordinate replay")
                for name in OBSERVABLES:
                    changes.append(dict(mode=mode, point_id=point["original_point_id"], xi=point["xi"], observable=name,
                        original=previous[name], selected=point[name], change_percent=100 * (float(point[name]) / float(previous[name]) - 1)))
            result = assembled_row(point, previous, template, thermal[point["point_id"]], template["mode"])
            mapped[identity] = result
            affected.add(identity[:2])
        rows = [mapped[k] for k in sorted(mapped)]
        require(len(rows) == old_count + added, "wrong merged point count")
        for identity, original in original_by_key.items():
            if mapped[identity]["source_kind"] == "v14_retained":
                require(all(mapped[identity][f] == original[f] for f in fields), "unrelated original value changed")
        for row in rows:
            require(all(math.isfinite(float(row[k])) and float(row[k]) > 0 for k in OBSERVABLES), "invalid plotted value")
            allowed_empty = {"iterations", "phase_boundary_xi_used"}
            require(all(row[k] != "" or k in allowed_empty or originals[0][k] == "" for k in fields),
                    "missing complete-grid field")
        write_csv(output / mode / "phase_guided_transport_grid.csv", rows, fields + list(EXTRA_FIELDS))
        # Retain source sampling rows verbatim; update only the selected coordinates.
        old_rates = read_csv(base / "phase_guided_transport_scan_sampling.csv")
        rate_fields = list(old_rates[0])
        update_keys = {rate_key(p) for p in updates}
        merged_rates = [r for r in old_rates if rate_key(r) not in update_keys]
        for point in updates:
            for r in channels:
                if r["point_id"] != point["point_id"]:
                    continue
                v = {k: r.get(k, "") for k in rate_fields}
                v.update(converged="true", run_id=f"mott_eta_fine_20261029_v1/{point['source_run']}")
                require(float(v["rate"]) > 0 and v["reason"] == "thermal_defect_and_two_shifted_audits", "unaccepted channel")
                require(rate_key(v) == rate_key(point) and r["variant"] == point["variant"] and r["source_run"] == point["source_run"],
                        "selected channel belongs to another point or precision")
                merged_rates.append(v)
        rates_by_point = defaultdict(list)
        for r in merged_rates:
            rates_by_point[rate_key(r)].append(r)
        require(len(rates_by_point) == len(rows), "sampling/main coordinate mismatch")
        for row in rows:
            rs = rates_by_point[rate_key(row)]
            require(len(rs) == len({r["channel"] for r in rs}) == 17, "missing or duplicate sampling channel")
        write_csv(output / mode / "phase_guided_transport_grid_sampling.csv", merged_rates, rate_fields)
        for panel, series in sorted(affected):
            group = [r for r in rows if r["plot_panel"] == panel and r["plot_series"] == series]
            for observable in ("tau_ubar", "eta_over_s", "zeta_over_s", "sigma_over_T"):
                by_id = {r["point_id"]: r for r in group}
                local_slopes = slopes(group, observable)
                for index, interval in enumerate(local_slopes):
                    a, b = by_id[interval["left"]], by_id[interval["right"]]
                    if (a["source_kind"] == "v14_retained") != (b["source_kind"] == "v14_retained"):
                        seam_rows.append(dict(mode=mode, panel=panel, series=series, observable=observable,
                            **interval, left_source=a["source_kind"], right_source=b["source_kind"],
                            preceding_slope=local_slopes[index - 1]["slope"] if index else "",
                            following_slope=local_slopes[index + 1]["slope"] if index + 1 < len(local_slopes) else ""))
        counts[mode] = len(rows)
        for index, row in enumerate(rows, 2):
            provenance.append(dict(mode=mode, grid_csv_line=index, point_id=row["point_id"], T_MeV=row["T_MeV"],
                muB_MeV=row["muB_MeV"], xi=row["xi"], source_kind=row["source_kind"],
                sampling_profile=row["sampling_profile"], source_run=row["source_run"],
                source_point_id=row["source_point_id"], source_csv_line=row["source_csv_line"]))
        evidence[mode] = dict(retained=old_count - len(updates) + added, recomputed=len(updates) - added, added=added,
                             quality_flags=sum(r["quality_flag"] == "true" for r in rows))
    require(sum(counts.values()) == 1867, "expected 1867 merged points")
    require(sum(c["quality_flags"] for c in evidence.values()) == 470, "original tau ratio flags changed")
    require(all(r["slope"] > 0 for r in seam_rows), "unexpected sign change at refinement boundary")
    # The two audited first-order gaps belong to untouched curves.
    audit = read_json(case / "inputs/base/background_audit.json")
    require(len(audit["known_boundary_pairs"]) == 2, "unexpected first-order boundary count")
    for pair in audit["known_boundary_pairs"]:
        mode = "mode_a" if pair["mode"] == "A" else "mode_b"
        grid = read_csv(output / mode / "phase_guided_transport_grid.csv")
        parameter = float(pair["curve_parameter"])
        curve = [r for r in grid if float(r["muB_MeV"]) == float(pair["muB_MeV"])
                 and float(r["alpha_T"] if mode == "mode_a" else r["T_MeV"]) == parameter]
        matches = [(a, b) for a, b in zip(curve, curve[1:]) if
                   (float(a["xi"]), float(b["xi"])) == (float(pair["xi_left"]), float(pair["xi_right"]))]
        require(len(matches) == 1 and all(r["source_kind"] == "v14_retained" for r in matches[0])
                and {r["phase_curr"] for r in matches[0]} == {"quark", "hadron"}, "changed first-order boundary")
    shutil.copyfile(case / "inputs/base/background_audit.json", output / "base_background_audit.json")
    write_csv(output / "point_provenance.csv", provenance)
    write_csv(output / "recomputed_original_changes.csv", changes)
    write_csv(output / "boundary_secants.csv", seam_rows)
    files = {p.relative_to(output).as_posix(): sha(p) for p in output.rglob("*") if p.is_file()}
    review = dict(schema="transport_eta_augmented_grid_v1", status="numerically_accepted_augmented_grid",
        calculation_sha=CALCULATION_SHA, workflow_ref="full-rs-eta-refinement-20261029-v1", points=1867,
        original_coordinates=1819, retained_original_values=1810, recomputed_original_values=9, added_points=48,
        channels=1867 * 17, mode_counts=counts, modes=evidence, quality_flags_retained=470,
        known_first_order_pairs=audit["known_boundary_pairs"], boundary_secants_positive=True,
        max_recomputed_observable_change_percent=max(abs(r["change_percent"]) for r in changes),
        full_grid_numerical_acceptance_passed=True, artifact_contract_passed=True,
        formal_production_promoted=False, publication_use_approved=False,
        source_selection=selection, thermodynamic_replay=thermal_review, output_sha256=files,
        missing_metadata=dict(iterations="not retained by the diagnostic driver; empty for its 57 selected points",
            phase_tracker="old labels retained at 9 original coordinates; unknown/empty at 48 new coordinates; no new phase classification claimed"),
        retained_provenance="inputs/base/point_provenance.csv retains the GitHub run and artifact identities for original points",
        numerical_policy="raw result selection with recorded stricter profiles; no interpolation, smoothing or value blending",
        limits="inherits qualified production settings and local xi/s evidence; no new all-channel outer/t/B0 convergence or eta amplitude ablation",
        generated_utc=datetime.now(timezone.utc).isoformat(), assembler_sha256=sha(Path(__file__)))
    write_json(output / "review.json", review)
    print(json.dumps(dict(status=review["status"], points=1867, added=48, recomputed=9, mode_counts=counts,
                          max_recomputed_observable_change_percent=review["max_recomputed_observable_change_percent"])))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "merge"))
    parser.add_argument("case", type=Path)
    parser.add_argument("eta_case", type=Path)
    args = parser.parse_args()
    (prepare if args.stage == "prepare" else merge)(args.case.resolve(), args.eta_case.resolve())


if __name__ == "__main__":
    main()
