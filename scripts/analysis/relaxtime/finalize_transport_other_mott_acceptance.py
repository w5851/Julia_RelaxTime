#!/usr/bin/env python3
"""Run the frozen complete-reference audit after the quote-remote recovery queue."""
from datetime import datetime, timezone
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time
import traceback

STAGE = "final_acceptance_v1"
REVIEW = "acceptance_review_complete_v1"
PACKAGE = "complete_acceptance_snapshot_v1"


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ready_state(state):
    if state.get("status") == "failed" or state.get("failed_runs"):
        raise ValueError("Recovery failed; preserve its existing evidence")
    if state.get("status") == "running_reference_recovery":
        return False
    if state.get("status") != "ready_for_agent_review":
        raise ValueError("Unexpected recovery state")
    if state.get("pending_jobs") != 0 or state.get("active_jobs") or state.get("external_jobs"):
        raise ValueError("Recovery reports ready while workers remain")
    references = state.get("reference_runs", [])
    if not references or len(references) != len(set(references)):
        raise ValueError("Invalid completed reference run identities")
    if set(references) != set(state.get("finished_runs", [])):
        raise ValueError("Final and finished reference runs disagree")
    return True


def validate_complete_summary(summary, expected_references):
    required = dict(status="completed_diagnostic_review", main_grid_complete=True,
        points=195, curves=13, channels=3315, source_files_checked=340, source_all_match=True,
        reference_evaluations=expected_references, expected_reference_evaluations=expected_references,
        reference_plan_complete=True, whole_task_complete=False)
    for key, expected in required.items():
        actual = summary.get(key)
        if actual != expected or isinstance(expected, bool) and actual is not expected:
            raise ValueError("Incomplete or ineligible acceptance summary: " + key)
    if summary.get("missing_reference_ids") != [] or summary.get("missing_xi_intervals") != []:
        raise ValueError("Missing reference points or xi split intervals")


def verify_frozen(case, config):
    for name, digest in config["files_sha256"].items():
        path = (case / name).resolve()
        if not path.is_relative_to(case) or sha(path) != digest:
            raise ValueError("Frozen postprocessing input changed: " + name)
    source = Path(config["source"]).resolve()
    hashes = read_json(case / "expected_source_hashes.json")
    if len(hashes) != 340 or any(sha(source / name) != digest for name, digest in hashes.items()):
        raise ValueError("Frozen numerical source changed")
    recovery = read_json(case / "recovery_v1.config.json")
    ids = recovery["expected_reference_ids"]
    if len(ids) != len(set(ids)) or len(ids) != config["expected_reference_evaluations"]:
        raise ValueError("Frozen reference plan coverage changed")


def run_stage(case, argv, log_path):
    with log_path.open("xb") as stream:
        subprocess.run(argv, cwd=case, stdin=subprocess.DEVNULL, stdout=stream,
            stderr=subprocess.STDOUT, creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0), check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("--wait", action="store_true", help="Wait without starting any numerical worker")
    args = parser.parse_args()
    case = args.case.resolve()
    config_path = case / "final_acceptance_setup_v1/config.json"
    config = read_json(config_path)
    config_digest = sha(config_path)
    state_path = case / (STAGE + ".status.json")
    lock_path = case / (STAGE + ".lock")
    with lock_path.open("x", encoding="utf-8") as lock:
        lock.write(str(os.getpid()))
        lock.flush()
        state = dict(status="waiting_for_reference_recovery", pid=os.getpid(),
            started_utc=datetime.now(timezone.utc).isoformat(), whole_task_complete=False,
            config_sha256=config_digest)

        def checkpoint(**changes):
            state.update(changes, updated_utc=datetime.now(timezone.utc).isoformat())
            state_path.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")

        try:
            if state_path.exists():
                raise FileExistsError("Retained postprocessing status already exists")
            checkpoint(status="verifying_frozen_inputs")
            verify_frozen(case, config)
            checkpoint(status="waiting_for_reference_recovery")
            while True:
                recovery = read_json(case / "recovery_v1.status.json")
                if ready_state(recovery):
                    break
                if not args.wait:
                    raise ValueError("Recovery is still running; use --wait or run after completion")
                time.sleep(30)
            verify_frozen(case, config)
            if sha(config_path) != config_digest:
                raise ValueError("Frozen postprocessing configuration changed")
            output = case / REVIEW
            package = case / "packages" / PACKAGE
            evidence = case / "final_acceptance_evidence_v1"
            if any(path.exists() for path in (output, package, evidence)):
                raise FileExistsError("Refusing to overwrite retained complete-reference evidence")
            evidence.mkdir()
            checkpoint(status="auditing_complete_reference_grid")
            audit = [sys.executable, str(case / config["acceptance_script"]), str(case), str(output),
                "--source", config["source"], "--require-complete"]
            run_stage(case, audit, evidence / "acceptance.log")
            summary = read_json(output / "summary.json")
            validate_complete_summary(summary, config["expected_reference_evaluations"])
            if set(summary["reference_runs"]) != set(recovery["reference_runs"]):
                raise ValueError("Accepted runs differ from the completed recovery queue")
            verify_frozen(case, config)
            completion = dict(schema="other_mott_complete_reference_execution_v1", audit_command=audit,
                summary_sha256=sha(output / "summary.json"), reference_grid_complete=True,
                reference_evaluations=summary["reference_evaluations"], whole_task_complete=False,
                scientific_verdict="pending_agent_review", production_data_changed=False,
                timestamp_utc=datetime.now(timezone.utc).isoformat())
            (evidence / "execution.json").write_text(json.dumps(completion, indent=2) + "\n", encoding="utf-8")
            runs = list(dict.fromkeys(["survey_v1", "endpoint_checks_v1", "locate_v1", "midpoint_screen_v1",
                *summary["fine_runs"], *summary["reference_runs"]]))
            reviews = ["survey_review_v2", "midpoint_review_v1", "fine_slopes_v1", "convergence_plan_v1",
                "final_slopes_recovery_v1", "pilot_plan_other01_v1", "recovery_evidence_v1", REVIEW,
                "final_acceptance_setup_v1", "final_acceptance_evidence_v1"]
            if (case / "acceptance_evidence_v1").is_dir():
                reviews.append("acceptance_evidence_v1")
            checkpoint(status="packaging_complete_reference_grid")
            run_stage(case, [sys.executable, str(case / config["packager_script"]), str(case), config["source"],
                PACKAGE, "--runs", ",".join(runs), "--reviews", ",".join(reviews)], case / (STAGE + ".package.log"))
            receipt = read_json(package / "archive_receipt.json")
            archive = package / receipt["archive"]
            if archive.stat().st_size != receipt["bytes"] or sha(archive) != receipt["sha256"]:
                raise ValueError("Complete-reference archive identity mismatch")
            checkpoint(status="ready_for_agent_review", reference_grid_complete=True,
                reference_evaluations=summary["reference_evaluations"], review=REVIEW + "/summary.json",
                package=str((package / "archive_receipt.json").relative_to(case)).replace("\\", "/"),
                finished_utc=datetime.now(timezone.utc).isoformat())
            print(json.dumps(state), flush=True)
        except BaseException as error:
            # The current run owns this status only after its first checkpoint.
            if state.get("updated_utc"):
                checkpoint(status="failed", error=str(error), traceback=traceback.format_exc())
            raise
        finally:
            lock.close()
            lock_path.unlink()


if __name__ == "__main__":
    main()
