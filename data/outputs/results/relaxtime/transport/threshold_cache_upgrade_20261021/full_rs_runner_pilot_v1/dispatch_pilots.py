"""Dispatch only this batch's seven prepared pilots, retaining each request attempt."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parent
manifest_path = ROOT / "dispatch_manifest.json"
spec = json.loads(manifest_path.read_text(encoding="utf-8"))
assert spec["remote_ref_verified"] and not spec["runs"]
assert len(spec["planned_runs"]) == 7
assert sum(entry["job"]["expected_points"] for entry in spec["planned_runs"]) == 11
assert spec["calculation_sha"] == spec["workflow_sha"] == "2e5c411c716d5cbc524c3fc6cfb335fe05b39444"
route = f"repos/{spec['repository']}/actions/workflows/{spec['workflow_id']}"


def matching_runs(title):
    raw = subprocess.check_output(["gh", "api", route + "/runs?event=workflow_dispatch&per_page=100"])
    return [run for run in json.loads(raw)["workflow_runs"] if
            run["head_sha"] == spec["workflow_sha"] and run["display_title"] == title]


for entry in spec["planned_runs"]:
    job = entry["job"]
    label = job["label"]
    request = ROOT / entry["request_file"]
    assert hashlib.sha256(request.read_bytes()).hexdigest() == entry["request_sha256"]
    payload = json.loads(request.read_text(encoding="utf-8"))
    assert payload["ref"] == spec["workflow_ref"] and payload["inputs"]["run_tier"] == "production"
    assert payload["inputs"]["calculation_sha"] == spec["calculation_sha"]
    title = f"transport {job['inputs']['case_name']} / {label} / production"
    attempt_path = ROOT / (label + ".dispatch_attempt.json")
    assert not attempt_path.exists(), "An attempt already exists; recover its run instead of dispatching again"
    assert not matching_runs(title), "A matching workflow already exists"
    attempt = dict(label=label, request_file=entry["request_file"], request_sha256=entry["request_sha256"],
                   started_utc=datetime.now(timezone.utc).isoformat(), request_may_have_been_accepted=True)
    attempt_path.write_text(json.dumps(attempt, indent=2)+"\n", encoding="utf-8")
    print("DISPATCH " + label, flush=True)
    result = subprocess.run(["gh", "api", "--method", "POST", route + "/dispatches", "--input", str(request)],
                            capture_output=True, encoding="utf-8", errors="replace")
    attempt.update(returncode=result.returncode, stdout=result.stdout, stderr=result.stderr,
                   response_utc=datetime.now(timezone.utc).isoformat())
    attempt_path.write_text(json.dumps(attempt, indent=2)+"\n", encoding="utf-8")
    if result.returncode:
        raise RuntimeError("Dispatch outcome requires inspection: " + label + " " + result.stderr)
    found = []
    for _ in range(24):
        found = matching_runs(title)
        if found:
            break
        time.sleep(2)
    assert len(found) == 1, "Recover the accepted request before any further dispatch"
    run = found[0]
    assert run["head_branch"] == spec["workflow_ref"]
    attempt.update(run_id=run["id"], run_url=run["html_url"], run_attempt=run["run_attempt"])
    attempt_path.write_text(json.dumps(attempt, indent=2)+"\n", encoding="utf-8")
    spec["runs"].append(dict(mode=entry["mode"], job=job, run_id=run["id"], run_attempt=run["run_attempt"],
                              run_url=run["html_url"], directory=f"runs/{run['id']}",
                              request_sha256=entry["request_sha256"]))
    manifest_path.write_text(json.dumps(spec, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(dict(label=label, run_id=run["id"], status=run["status"], url=run["html_url"])), flush=True)
print("Seven pilots dispatched; full-grid dispatch remains unstarted", flush=True)
