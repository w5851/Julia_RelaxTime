#!/usr/bin/env python3
"""Wait for a free remote numerical slot, then audit moving-cutoff membership."""
from pathlib import Path
from types import SimpleNamespace
import importlib.util
import json
import subprocess
import sys
import time

case=Path(sys.argv[1]).resolve()
spec=importlib.util.spec_from_file_location("sigma_runner",case/"run_sigma_mott_refinement_v2.py")
runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)
path=case/"cutoff_queue_v1.status.json"
with path.open("x",encoding="utf-8") as stream:
    json.dump(dict(status="waiting_for_slot",started_utc=runner.now()),stream)
while True:
    processes=subprocess.run(["tasklist","/FI","IMAGENAME eq julia.exe","/FO","CSV","/NH"],
        capture_output=True,text=True,check=True).stdout
    count=sum(line.lower().startswith('"julia.exe"') for line in processes.splitlines())
    if count<3 and runner.free_gib()>3.5:break
    time.sleep(20)
try:
    runner.job(case,SimpleNamespace(name="cutoff_audit_v1",stage="spectrum",profile="s6",points="all",
                                   script="audit_sigma_mott_cutoff.jl"))
    runner.write_json(path,dict(status="completed_diagnostic",finished_utc=runner.now()))
except BaseException as exc:
    runner.write_json(path,dict(status="failed",error=str(exc),finished_utc=runner.now()))
    raise
