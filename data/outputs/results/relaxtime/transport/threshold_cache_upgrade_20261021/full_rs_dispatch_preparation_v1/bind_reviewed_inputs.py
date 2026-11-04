"""Bind the received mode-A/B runner drafts to the reviewed calculation commit.

This batch-specific controller checks hashes and Git content only. Numerical
qualification is supplied by the retained, independently reviewed remote runs.
It creates new input files and does not dispatch workflows or promote results.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys
import zipfile

sys.dont_write_bytecode = True
EVIDENCE = Path("data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021")
WORKFLOW = ".github/workflows/relaxtime-phase-guided-transport-production.yml"
TIMING_ENTRYPOINT = "benchmark/relaxtime/transport_runner_timing.jl"
TIMING_CASE = "full_rs_runner_timing_preflight_setup_v1"
SPECS = {
    "A": {
        "case": "full_rs_mode_a_runner_inputs_v1",
        "draft": "runner_inputs.json",
        "review": "preparation_manifest.json",
        "production": "full_rs_production_replay_v3/returned/production_replay",
        "full_points": 910,
        "scan_mode": "fixed-muB-phase-scaled",
    },
    "B": {
        "case": "full_rs_mode_b_qualification_setup_v1",
        "draft": "runner_inputs/runner_inputs.json",
        "review": "reviewed_result.json",
        "production": "full_rs_mode_b_qualification_setup_v1/returned/production",
        "full_points": 909,
        "scan_mode": "fixed-T-sparse-muB",
    },
}


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def received_case(folder):
    declared = read(folder / "returned_archive_manifest.json")
    receipt = read(folder / "local_receipt.json")
    archive = folder / "returned_evidence.zip"
    require(receipt["archive"] == declared, "receipt and archive identity differ")
    require(digest(archive) == declared["sha256"] and archive.stat().st_size == declared["bytes"],
            "received archive changed")
    unpacked = (folder / "returned").resolve()
    with zipfile.ZipFile(archive) as bundle:
        names = bundle.namelist()
        require(len(names) == len(set(names)) == declared["files"], "invalid archive member set")
        for name in names:
            target = (unpacked / name).resolve()
            require(target.is_relative_to(unpacked), "archive member leaves the received directory")
            require(target.read_bytes() == bundle.read(name), f"received evidence changed: {name}")
    return declared


def job_points(jobs, draft, mode):
    coordinates = set()
    labels = set()
    for job in jobs:
        args = job["inputs"]
        require(job["label"] == args["shard_label"] and job["label"] not in labels,
                "duplicate or inconsistent shard label")
        labels.add(job["label"])
        require(args["case_name"] == draft["case_name"] and
                args["sampling_profile"] == draft["sampling_profile"] and
                args["mode"] == SPECS[mode]["scan_mode"] and
                args["run_tier"] == "production" and args["render_plots"] is False,
                "draft job differs from the reviewed calculation scope")
        mu = [float(v) for v in args["muB_list"].split(",")]
        second = [float(v) for v in args["alpha_t_list" if mode == "A" else "t_list"].split(",")]
        xi = [float(v) for v in args["xi_list"].split(",")]
        points = {(m, s, x) for m in mu for s in second for x in xi}
        require(len(points) == job["expected_points"] and not (coordinates & points),
                "duplicate or incomplete draft point set")
        coordinates.update(points)
    return coordinates


def bind(root, commit, workflow_ref, outdir):
    require(re.fullmatch(r"[0-9a-f]{40}", commit) is not None, "use a full calculation commit")
    require(not outdir.exists(), "use a new output directory")
    require(subprocess.check_output(["git", "-C", str(root), "rev-parse", workflow_ref],
                                    text=True).strip() == commit,
            "local workflow ref does not identify the requested calculation commit")
    execution_files = {}
    for relative in (WORKFLOW, TIMING_ENTRYPOINT):
        committed = subprocess.check_output(["git", "-C", str(root), "show", f"{commit}:{relative}"])
        current = (root / relative).read_bytes()
        require(committed.replace(b"\r\n", b"\n") == current.replace(b"\r\n", b"\n"),
                f"execution file differs from the selected commit: {relative}")
        execution_files[relative] = dict(committed_sha256=hashlib.sha256(committed).hexdigest(),
                                        worktree_sha256=hashlib.sha256(current).hexdigest())
    helper_path = root / "scripts/analysis/relaxtime/prepare_transport_production_inputs.py"
    spec = importlib.util.spec_from_file_location("transport_inputs", helper_path)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    base = root / EVIDENCE
    timing_folder = base / TIMING_CASE
    timing_archive = received_case(timing_folder)
    timing_review_path = timing_folder / "returned/reviewed_result.json"
    timing_review = read(timing_review_path)
    timing_package = read(timing_folder / "returned/package_manifest.json")
    require(timing_review["timing_contract_passed"] and
            timing_review["actual_dry_run_planned_points"] == 2 and
            timing_review["unit_assertions"] >= 25,
            "runner timing preflight is incomplete")
    require(timing_review["overlay_archive_sha256"] == timing_package["overlay_archive_sha256"],
            "timing review and source overlay differ")
    for relative, checksum in timing_package["overlay_files"].items():
        require(digest(root / relative) == checksum,
                f"timing preflight no longer covers the current file: {relative}")
    prepared, provenance, common_sources = {}, {}, None
    for mode, entry in SPECS.items():
        folder = base / entry["case"]
        archive = received_case(folder)
        draft_path = folder / "returned" / entry["draft"]
        review = read(folder / "returned" / entry["review"])
        draft = read(draft_path)
        require(digest(draft_path) == review["runner_inputs_sha256"], "unreviewed runner input")
        require(draft["calculation_sha"] is None and not draft["pilot_input_ready"] and
                not draft["full_dispatch_ready"], "expected an unbound provisional draft")
        if mode == "B":
            require(review["mode_b_background_qualified"] and
                    review["mode_b_transport_representatives_qualified"] and
                    review["layer_targets_passed"], "mode-B qualification is incomplete")
        production = base / entry["production"]
        manifest = read(production / "manifest.json")
        sources = helper.verified_sources(production, manifest, root)
        helper.verified_calculation_commit(commit, root, sources)
        require(common_sources is None or sources == common_sources,
                "mode-A and mode-B use different numerical code")
        common_sources = sources
        require(all(timing_package["source_files"].get(path) == checksum
                    for path, checksum in sources.items()),
                "runner timing was checked against different numerical code")
        require(hashlib.sha256(json.dumps(sources, sort_keys=True).encode()).hexdigest() ==
                draft["numerical_source_sha256"], "draft numerical source identity differs")
        profile = root / "config/models/pnjl/transport_sampling" / (draft["sampling_profile"] + ".toml")
        require(digest(profile) == draft["sampling_profile_sha256"], "draft profile changed")
        pilot = job_points(draft["pilot_jobs"], draft, mode)
        full = job_points(draft["provisional_full_jobs"], draft, mode)
        require(pilot and pilot <= full and len(full) == draft["full_points"] == entry["full_points"],
                "pilot and full plan are inconsistent")
        payload = deepcopy(draft)
        payload.update(mode=mode, calculation_sha=commit, workflow_ref=workflow_ref,
                       workflow_sha=commit, pilot_input_ready=True, full_dispatch_ready=False,
                       runner_cost_qualified=False, remote_ref_verified=False,
                       remaining_before_dispatch=["publish and verify the exact workflow ref on GitHub"],
                       unbound_runner_inputs_sha256=digest(draft_path),
                       received_archive_sha256=archive["sha256"])
        for key in ("pilot_jobs", "provisional_full_jobs"):
            for job in payload[key]:
                job["inputs"]["calculation_sha"] = commit
        prepared[mode] = payload
        provenance[mode] = dict(received_archive=archive, draft_sha256=digest(draft_path),
                                source_manifest_sha256=digest(production / "manifest.json"),
                                numerical_sources=len(sources), pilot_points=len(pilot),
                                pilot_jobs=len(payload["pilot_jobs"]), full_points=len(full))
    outdir.mkdir(parents=True)
    written = {}
    for mode, payload in prepared.items():
        path = outdir / f"mode_{mode.lower()}_runner_inputs.json"
        path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        written[path.name] = digest(path)
        for job in payload["pilot_jobs"]:
            path = outdir / f"{job['label']}.inputs.json"
            path.write_text(json.dumps(job["inputs"], indent=2) + "\n", encoding="utf-8")
            written[path.name] = digest(path)
    manifest = dict(created_utc=datetime.now(timezone.utc).isoformat(), calculation_sha=commit,
                    workflow=WORKFLOW, workflow_ref=workflow_ref, workflow_sha=commit,
                    workflow_file_sha256=digest(root / WORKFLOW),
                    execution_files=execution_files,
                    timing_preflight=dict(received_archive=timing_archive,
                                          review_sha256=digest(timing_review_path),
                                          unit_assertions=timing_review["unit_assertions"],
                                          actual_dry_run_planned_points=2),
                    controller_sha256=digest(Path(__file__)), evidence=provenance, files=written,
                    actions_dispatched=False, remote_ref_verified=False,
                    full_dispatch_ready=False, runner_cost_qualified=False,
                    check_scope="Received bytes, source/profile hashes, Git content and job metadata; no numerical execution")
    (outdir / "binding_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_root", type=Path)
    parser.add_argument("calculation_sha")
    parser.add_argument("workflow_ref")
    parser.add_argument("outdir", type=Path)
    args = parser.parse_args()
    print(json.dumps(bind(args.source_root.resolve(), args.calculation_sha, args.workflow_ref,
                          args.outdir.resolve()), indent=2))
