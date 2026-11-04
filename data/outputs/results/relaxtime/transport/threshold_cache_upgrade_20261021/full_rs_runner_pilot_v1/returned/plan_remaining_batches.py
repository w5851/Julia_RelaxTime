"""Prepare complete, nonoverlapping batches after the actual runner pilot review."""
from __future__ import annotations

import argparse
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def group_key(job, mode):
    args = job["inputs"]
    mu = args["muB_list"].split(",")
    second = args["alpha_t_list" if mode == "A" else "t_list"].split(",")
    require(len(mu) == len(second) == 1, "a batch must follow one scan curve")
    return (float(mu[0]), float(second[0]))


def plan(root):
    specification = read(root / "dispatch_manifest.json")
    review = read(root / "reviewed_pilot.json")
    require(review["pilot_artifacts_passed"] and review["runner_cost_qualified"], "pilot is not qualified")
    require(review["dispatch_manifest_sha256"] == digest(root / "dispatch_manifest.json"), "dispatch changed after review")
    output = root / "remaining_batches"
    require(not output.exists(), "use a new batch-plan directory")
    reviewed_runs = {run["run_id"]: run for run in review["runs"]}
    require(len(reviewed_runs) == len(specification["runs"]), "incomplete pilot review")
    prepared = {}
    for mode in ("A", "B"):
        filename = specification["bound_inputs"][mode]
        require(digest(root / filename) == specification["bound_file_hashes"][filename], "bound grid changed")
        original = read(root / filename)
        curves, all_points = {}, set()
        for job in original["provisional_full_jobs"]:
            group = group_key(job, mode)
            curve = curves.setdefault(group, dict(template=deepcopy(job), xi={}))
            require(curve["template"]["expected_T_MeV"] == job["expected_T_MeV"], "curve temperature changed")
            for token in job["inputs"]["xi_list"].split(","):
                point = (*group, float(token))
                require(point not in all_points, "duplicate full-grid point")
                all_points.add(point)
                curve["xi"][float(token)] = token
        require(len(all_points) == original["full_points"], "full-grid count differs")
        completed, reuse = set(), []
        for entry in specification["runs"]:
            if entry["mode"] != mode:
                continue
            job = entry["job"]
            run = reviewed_runs[entry["run_id"]]
            group = group_key(job, mode)
            require(run["calculation_sha"] == specification["calculation_sha"], "pilot calculation differs")
            require(len(run["channel_times"]) == run["completed_points"] == job["expected_points"],
                    "actual reviewed pilot point set is incomplete")
            points = {(*group, row["xi"]) for row in run["channel_times"]}
            requested = {(*group, float(xi)) for xi in job["inputs"]["xi_list"].split(",")}
            require(points == requested and len(points) == job["expected_points"] and
                    not (points & completed) and points <= all_points, "pilot reuse point set differs")
            completed.update(points)
            reuse.append(dict(run_id=run["run_id"], run_url=run["run_url"], label=run["label"],
                              artifact_sha256=run["artifact_sha256"], points=run["completed_points"],
                              actual_coordinates=[{k: row[k] for k in ("T_MeV", "muB_MeV", "xi")}
                                                  for row in run["channel_times"]]))
        budget = review["modes"][mode]
        width = budget["max_points_per_shard"]
        require(isinstance(width, int) and width >= 1, "invalid accepted batch width")
        remaining, jobs = set(), []
        for group, curve in sorted(curves.items()):
            pending = [xi for xi in sorted(curve["xi"]) if (*group, xi) not in completed]
            for index, start in enumerate(range(0, len(pending), width)):
                chunk = pending[start:start+width]
                job = deepcopy(curve["template"])
                mu, second = group
                label = (f"full_modeA_muB{mu:g}_alpha{second:g}_{index:03d}" if mode == "A"
                         else f"full_modeB_T{second:g}_muB{mu:g}_{index:03d}")
                job["label"] = label
                job["inputs"]["shard_label"] = label
                job["inputs"]["xi_list"] = ",".join(curve["xi"][xi] for xi in chunk)
                job["expected_points"] = len(chunk)
                job["planned_job_seconds"] = (budget["reserved_overhead_seconds"] +
                                               len(chunk) * budget["per_point_budget_seconds"])
                require(job["planned_job_seconds"] <= budget["job_target_seconds"] + 1e-8,
                        "batch exceeds the reviewed planning budget")
                require(job["inputs"]["calculation_sha"] == specification["calculation_sha"] and
                        job["inputs"]["run_tier"] == "production" and job["inputs"]["render_plots"] is False,
                        "batch changed the calculation contract")
                points = {(*group, xi) for xi in chunk}
                require(not (points & remaining), "duplicate remaining batch point")
                remaining.update(points)
                jobs.append(job)
        require(not (remaining & completed) and remaining | completed == all_points,
                "remaining and reused points do not exactly cover the full grid")
        require(len({job["label"] for job in jobs}) == len(jobs), "duplicate batch labels")
        prepared[mode] = dict(mode=mode, case_name=original["case_name"], sampling_profile=original["sampling_profile"],
                              calculation_sha=specification["calculation_sha"], workflow_sha=specification["workflow_sha"],
                              workflow_ref=specification["workflow_ref"], workflow_path=specification["workflow_path"],
                              required_total_points=len(all_points), reused_pilot_points=len(completed),
                              remaining_points=len(remaining), reuse_pilot_runs=reuse, budget=budget, jobs=jobs,
                              full_dispatch_ready=True, runner_cost_qualified=True, whole_grid_transport_calculated=False,
                              full_jobs_dispatched=False, formal_production_promoted=False,
                              review_sha256=digest(root / "reviewed_pilot.json"),
                              coverage_sha256=hashlib.sha256(json.dumps(sorted(all_points)).encode()).hexdigest())
    output.mkdir()
    written = {}
    for mode, payload in prepared.items():
        path = output / f"mode_{mode.lower()}_batches.json"
        path.write_text(json.dumps(payload, indent=2)+"\n", encoding="utf-8")
        written[path.name] = digest(path)
        for job in payload["jobs"]:
            request = dict(ref=specification["workflow_ref"],
                           inputs={key: str(value).lower() if isinstance(value, bool) else value
                                   for key, value in job["inputs"].items()})
            path = output / (job["label"] + ".request.json")
            path.write_text(json.dumps(request, indent=2)+"\n", encoding="utf-8")
            written[path.name] = digest(path)
    manifest = dict(created_utc=datetime.now(timezone.utc).isoformat(), calculation_sha=specification["calculation_sha"],
                    controller_sha256=digest(Path(__file__)), review_sha256=digest(root / "reviewed_pilot.json"),
                    modes={mode: {key: payload[key] for key in ("required_total_points", "reused_pilot_points", "remaining_points")}
                           | dict(jobs=len(payload["jobs"]), max_points_per_shard=payload["budget"]["max_points_per_shard"])
                           for mode, payload in prepared.items()}, files=written, full_jobs_dispatched=False,
                    merge_contract="Preserve each run artifact and its actual coordinates; verify each curve/xi key exactly once across reused pilot and remaining runs before any promotion")
    (output / "batches_manifest.json").write_text(json.dumps(manifest, indent=2)+"\n", encoding="utf-8")
    return {key: value for key, value in manifest.items() if key != "files"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    print(json.dumps(plan(args.directory.resolve()), indent=2))
