"""Review collected full-RS shards on quote-remote and merge a complete grid.

The original CSV values and artifacts remain unchanged. Merged tables are only
written when all requested points, including reused pilots, pass the artifact
contract. Partial reviews retain errors and never certify a complete grid.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(line for line in stream if line.strip() and not line.startswith("#"))
        rows = list(reader)
        return reader.fieldnames, rows


def write_csv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def job_points(job, mode):
    args = job["inputs"]
    mu = [float(v) for v in args["muB_list"].split(",")]
    second = [float(v) for v in args["alpha_t_list" if mode == "A" else "t_list"].split(",")]
    xis = [float(v) for v in args["xi_list"].split(",")]
    require(len(mu) == len(second) == 1, "shard must follow exactly one curve")
    points = [(mode, mu[0], second[0], value) for value in xis]
    require(len(points) == job["expected_points"] == len(set(points)), "duplicate or missing planned point")
    return points


def coverage(bound, entries, *, partial=False):
    expected, requested = {}, set()
    for mode, original in bound.items():
        for job in original["provisional_full_jobs"]:
            for point in job_points(job, mode):
                require(point not in expected, "duplicate point in bound full grid")
                expected[point] = float(job["expected_T_MeV"])
        require(sum(point[0] == mode for point in expected) == original["full_points"],
                "bound full-grid count differs")
    for entry in entries:
        for point in job_points(entry["job"], entry["mode"]):
            require(point not in requested, "overlap between shards or reused pilots")
            require(point in expected and abs(float(entry["job"]["expected_T_MeV"]) - expected[point]) <= 1e-10,
                    "shard changes the bound grid or its phase anchor")
            requested.add(point)
    require(requested <= expected.keys() if partial else requested == expected.keys(),
            "requested shards do not cover the bound full grid")
    return expected


def point_key(row, mode):
    mu = float(row["muB_MeV"])
    second = float(row["alpha_T"] if mode == "A" else row["T_MeV"])
    xi = float(row["xi"])
    require(all(math.isfinite(v) for v in (mu, second, xi)), "nonfinite point coordinate")
    return mode, mu, second, xi


def verify_quality_flag(row):
    # scan_quality.jl marks a finite particle/antiparticle ratio above 6.0.
    # This heuristic is retained for review; it is not a convergence failure.
    values = (float(row["tau_u"]), float(row["tau_ubar"]))
    require(all(math.isfinite(value) and value > 0 for value in values), "invalid quality ratio input")
    ratio = max(values) / min(values)
    metric = float(row["quality_metric"])
    high = ratio > 6.0
    require(row["quality_flag"].lower() == str(high).lower() and
            row["quality_reason"] == ("tau_u_ubar_ratio_high" if high else "ok") and
            math.isfinite(metric) and math.isclose(metric, ratio, rel_tol=1e-12, abs_tol=1e-12),
            "quality flag differs from the recorded finite-ratio contract")


def verify_rows(rows, diagnostics, entry, expected, seen):
    mode = entry["mode"]
    require(len(diagnostics) == len(rows) * 17, "extra or missing channel diagnostics")
    local = set()
    for row in rows:
        point = point_key(row, mode)
        require(point in expected and abs(float(row["T_MeV"]) - expected[point]) <= 1e-5,
                "actual point or temperature differs from bound grid")
        require(point not in seen and point not in local, "duplicate actual point across runs")
        verify_quality_flag(row)
        residual = float(row["residual_norm"])
        require(math.isfinite(residual) and 0 <= residual <= 1e-8, "background residual exceeds qualified tolerance")
        local.add(point)
    require(local == set(job_points(entry["job"], mode)), "actual curve coordinates differ from shard")
    return local


def review(directory, output, pilot_reviewer, *, partial=False):
    require(not output.exists(), "use a new review output directory")
    output.mkdir(parents=True)
    pilot = load_module("full_grid_pilot", pilot_reviewer)
    specification = read(directory / "dispatch_manifest.json")
    plan = read(directory / "queue_plan.json")
    require(all(specification[k] == value for k, value in plan.items() if k != "runs"), "frozen queue plan changed")
    for relative, checksum in plan["frozen_files"].items():
        require(digest(pilot.child(directory, relative)) == checksum, "frozen input changed: " + relative)
    source = read(directory / "source_manifest.json")
    source_root = directory / "source"
    require(source["calculation_sha"] == specification["calculation_sha"] == specification["workflow_sha"],
            "source archive does not match requested calculation")
    pilot.verified_zip(directory / "calculation-source.zip", source_root, source["sha256"])
    old_spec = read(directory / "inputs/pilot_dispatch_manifest.json")
    old_review = read(directory / "inputs/reviewed_pilot.json")
    require(old_review["pilot_artifacts_passed"] and old_review["dispatch_manifest_sha256"] ==
            digest(directory / "inputs/pilot_dispatch_manifest.json"), "pilot reuse identity differs")
    require(old_spec["runs"] == specification["reused_pilot_runs"], "reused pilot entries changed")
    bound = {mode: read(directory / f"inputs/mode_{mode.lower()}_runner_inputs.json") for mode in ("A", "B")}
    for mode in bound:
        name = f"mode_{mode.lower()}_runner_inputs.json"
        require(digest(directory / "inputs" / name) == old_spec["bound_file_hashes"][name], "bound grid changed")
        batches = read(directory / f"inputs/mode_{mode.lower()}_batches.json")
        require([entry["job"] for entry in specification["planned_runs"] if entry["mode"] == mode] == batches["jobs"],
                "planned jobs differ from reviewed batch requests")
    expected = coverage(bound, specification["planned_runs"] + old_spec["runs"])
    planned = {entry["job"]["label"]: entry for entry in specification["planned_runs"]}
    entries = list(specification["runs"])
    require(len({entry["run_id"] for entry in entries}) == len(entries), "duplicate full-grid run ID")
    for entry in entries:
        frozen = planned.get(entry["job"]["label"])
        require(frozen is not None and all(entry[key] == value for key, value in frozen.items()),
                "recorded run differs from its frozen request")
    entries += [dict(entry, directory="pilot_runs/" + str(entry["run_id"])) for entry in old_spec["runs"]]
    require(len({entry["run_id"] for entry in entries}) == len(entries), "pilot and new run IDs overlap")
    seen, reports, errors, pending, provenance, quality_flags = set(), [], [], [], [], []
    merged, samples, fields, sample_fields = defaultdict(list), defaultdict(list), {}, {}
    for entry in entries:
        run_root = pilot.child(directory, entry["directory"])
        if not (run_root / "archive_manifest.json").exists():
            pending.append(entry["run_id"])
            continue
        try:
            result = pilot.review_run(run_root, entry, specification, source_root)
            manifests = list((run_root / "artifact").rglob("run_manifest.json"))
            require(len(manifests) == 1, "ambiguous artifact result")
            result_root = manifests[0].parent
            csv_path = result_root / "phase_guided_transport_scan.csv"
            header, rows = read_csv(csv_path)
            sample_header, diagnostics = read_csv(result_root / "phase_guided_transport_scan_sampling.csv")
            mode = entry["mode"]
            require(mode not in fields or fields[mode] == header, "result CSV headers differ between shards")
            require(mode not in sample_fields or sample_fields[mode] == sample_header, "sampling CSV headers differ")
            require(not read_csv(result_root / "failed_points.csv")[1], "failed attempts need an explicit audit")
            keys = verify_rows(rows, diagnostics, entry, expected, seen)
            seen.update(keys)
            fields[mode], sample_fields[mode] = header, sample_header
            merged[mode].extend(rows)
            samples[mode].extend(diagnostics)
            for row in rows:
                if row["quality_flag"].lower() == "true":
                    quality_flags.append(dict(mode=mode, run_id=entry["run_id"],
                         **{key: row[key] for key in ("T_MeV", "muB_MeV", "xi", "quality_reason", "quality_metric")}))
                provenance.append(dict(mode=mode, muB_MeV=row["muB_MeV"], T_MeV=row["T_MeV"],
                      alpha_T=row["alpha_T"], xi=row["xi"], github_run_id=entry["run_id"],
                      github_run_attempt=entry["run_attempt"], github_run_url=entry["run_url"],
                      artifact_sha256=result["artifact_sha256"], source_csv=csv_path.relative_to(directory).as_posix(),
                      source_csv_sha256=digest(csv_path), result_run_id=row["run_id"]))
            reports.append(result)
        except (ValueError, KeyError, OSError) as error:
            errors.append(dict(run_id=entry["run_id"], label=entry["job"]["label"], error=str(error)))
        print(json.dumps(dict(reviewed=len(reports), points=len(seen), errors=len(errors), run_id=entry["run_id"])), flush=True)
    complete = seen == expected.keys() and len(specification["runs"]) == len(specification["planned_runs"])
    passed = complete and not errors and not pending
    merged_hashes = {}
    if passed:
        for mode in bound:
            rows = sorted(merged[mode], key=lambda row: point_key(row, mode))
            destination = output / f"mode_{mode.lower()}/phase_guided_transport_scan.csv"
            write_csv(destination, fields[mode], rows)
            diagnostics = sorted(samples[mode], key=lambda row: (float(row["muB_MeV"]), float(row["T_MeV"]),
                                  float(row["xi"]), row["channel"]))
            sample_path = destination.with_name("phase_guided_transport_scan_sampling.csv")
            write_csv(sample_path, sample_fields[mode], diagnostics)
            for path in (destination, sample_path):
                merged_hashes[path.relative_to(output).as_posix()] = digest(path)
        provenance.sort(key=lambda row: (row["mode"], float(row["muB_MeV"]), float(row["T_MeV"]), float(row["xi"])))
        destination = output / "point_provenance.csv"
        write_csv(destination, list(provenance[0]), provenance)
        merged_hashes[destination.name] = digest(destination)
    report = dict(created_utc=datetime.now(timezone.utc).isoformat(), numerical_review_host="quote-remote",
             calculation_sha=specification["calculation_sha"], workflow_ref=specification["workflow_ref"],
             grid_artifact_contract_passed=passed, whole_grid_transport_calculated=complete,
             formal_production_promoted=False, publication_use_approved=False,
             expected_points=len(expected), completed_points=len(seen), reviewed_runs=len(reports),
             qualified_channels=sum(item["qualified_channels"] for item in reports),
             modes={mode: dict(expected=bound[mode]["full_points"], completed=len(merged[mode])) for mode in bound},
             missing_points=[list(point) for point in sorted(expected.keys() - seen)],
             errors=errors, pending_run_ids=pending, runs=reports, merged_files=merged_hashes,
             quality_flags=quality_flags, quality_ratio_threshold=6.0,
             quality_flag_semantics="Finite tau_u/tau_ubar asymmetry heuristic; retain and audit with the qualified representative evidence, not a numerical-error estimator",
             reviewer_sha256=digest(Path(__file__)), pilot_reviewer_sha256=digest(pilot_reviewer),
             input_metadata_sha256={name: digest(directory / name) for name in
                   ("dispatch_manifest.json", "queue_plan.json", "source_manifest.json")},
             supplemental_background_and_local_jump_audit_required=passed)
    with (output / "reviewed_grid.json").open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    print(json.dumps({k: v for k, v in report.items() if k not in ("runs", "missing_points", "quality_flags")}), flush=True)
    return 0 if passed or (partial and not errors) else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--pilot-reviewer", type=Path, default=Path(__file__).with_name("review_transport_runner_pilot.py"))
    parser.add_argument("--partial", action="store_true")
    args = parser.parse_args()
    sys.exit(review(args.directory.resolve(), args.output.resolve(), args.pilot_reviewer.resolve(), partial=args.partial))
