"""Coordinate the finite Andronic study with two already-started Julia workers.

Each numerical request calls the unchanged density runner. Design, assessment,
resolution review and plotting remain separate frozen stages with their hashes.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
from scripts.plotting.plot_manifest import sha256_file


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix",default="andronic_envelope_20261018")
    parser.add_argument("--experiment",type=Path,
        default=ROOT/"data/inputs/meson_freezeout_experiment_20261018/experimental_ratios.csv")
    args=parser.parse_args()
    base=ROOT/"data/outputs/results/relaxtime/meson_density"
    prefix=base/args.prefix
    case=lambda suffix:Path(str(prefix)+"_"+suffix)
    execution=case("execution")
    if not execution.is_dir():raise ValueError("prepare the initial design and start both worker processes first")
    log=(execution/"pipeline.log").open("a",encoding="utf-8",buffering=1)
    sequence=[1,1]
    for part in (0,1):
        while (execution/f"worker_{part}_{sequence[part]:03d}.json").exists():sequence[part]+=1

    def say(value):
        text=json.dumps(value,ensure_ascii=False)
        print(text,flush=True);log.write(text+"\n")

    def run(script,*argv):
        command=[sys.executable,str(ROOT/"scripts/analysis/relaxtime"/script),*map(str,argv)]
        say({"command":command})
        process=subprocess.run(command,cwd=ROOT,text=True,encoding="utf-8",capture_output=True)
        log.write(process.stdout+process.stderr)
        if process.stdout:print(process.stdout,flush=True)
        if process.returncode:
            if process.stderr:print(process.stderr,file=sys.stderr,flush=True)
            raise RuntimeError(f"analysis stage failed: {script}")
        return process.stdout

    def wait(paths,label):
        started=time.monotonic();last=-30.
        while not all(p.is_file() for p in paths):
            elapsed=time.monotonic()-started
            if elapsed-last>=30:
                say({"waiting_for":label,"elapsed_seconds":round(elapsed),"pending":[str(p) for p in paths if not p.exists()]})
                last=elapsed
            time.sleep(5)
        say({"completed":label})

    def publish(part,**job):
        path=execution/f"worker_{part}_{sequence[part]:03d}.json"
        sequence[part]+=1
        value=dict(schema="freezeout_envelope_worker_request_v1",**job)
        temporary=path.with_suffix(".pending")
        temporary.write_text(json.dumps(value,indent=2)+"\n",encoding="utf-8")
        temporary.replace(path)
        return path.with_name(path.stem+"_done.json")

    def dispatch_design(design_path,stage):
        design=json.loads(design_path.read_text(encoding="utf-8"));jobs=[];runs=[]
        for part,names in enumerate(design["parts"]):
            if not names:continue
            output=case(f"{stage}_{part}");runs.append(output)
            jobs.append(publish(part,action="scan",output=str(output),profiles=names,
                energies_GeV=design["design"]["screening_energies_GeV"],resolution="screening",
                profile_dir=str(design_path.parent/"profiles")))
        wait(jobs,stage)
        return runs

    def assess(designs,runs,stage):
        output=case(stage)
        run("andronic_parameter_envelope.py","assess","--designs",*designs,"--runs",*runs,
            "--experiment",args.experiment,"--output",output)
        return output/"assessment.json"

    initial=case("design")/"design.json"
    runs=[case(f"screen_{part}") for part in (0,1)]
    wait([p/"manifest.json" for p in runs],"initial 46-curve screening")
    designs=[initial]
    assessment=assess(designs,runs,"initial_assessment")
    local=case("local_design")
    run("andronic_parameter_envelope.py","refine","--assessment",assessment,"--designs",*designs,"--output",local)
    runs+=dispatch_design(local/"design.json","local_screen");designs.append(local/"design.json")
    assessment=assess(designs,runs,"local_assessment")
    mapping=case("mapping_design")
    output=run("andronic_parameter_envelope.py","horn-map","--assessment",assessment,"--designs",*designs,"--output",mapping)
    mapping_decision=dict(checked_utc=datetime.now(timezone.utc).isoformat(),output=output,
        design_created=(mapping/"design.json").exists(),policy="Exact T(muB)-preserving joint A,E rescaling; reject outside the selected box")
    (execution/"mapping_decision.json").write_text(json.dumps(mapping_decision,indent=2)+"\n",encoding="utf-8")
    if mapping_decision["design_created"]:
        runs+=dispatch_design(mapping/"design.json","mapping_screen");designs.append(mapping/"design.json")
    assessment=assess(designs,runs,"all_screen_assessment")
    selection=case("selected_design")
    run("andronic_parameter_envelope.py","select","--assessment",assessment,"--designs",*designs,"--output",selection)
    selected=json.loads((selection/"design.json").read_text(encoding="utf-8"))
    names=[c["name"] for c in selected["candidates"]]
    plan=json.loads((execution/"numerical_plan.json").read_text(encoding="utf-8"))
    jobs=[]
    jobs.append(publish(0,action="scan",output=str(case("standard")),profiles=names,
        energies_GeV=plan["standard_energies_GeV"],resolution="standard",profile_dir=str(selection/"profiles")))
    jobs.append(publish(0,action="scan",output=str(case("fine_0")),profiles=names[:1],
        energies_GeV=plan["fine_energies_GeV"],resolution="fine",profile_dir=str(selection/"profiles")))
    jobs.append(publish(1,action="scan",output=str(case("fine_1")),profiles=names[1:],
        energies_GeV=plan["fine_energies_GeV"],resolution="fine",profile_dir=str(selection/"profiles")))
    wait(jobs,"selected standard/fine refinement")
    for part in (0,1):publish(part,action="stop")
    old_summary=base/"freezeout_envelope_20261018_summary/summary.json"
    old=json.loads(old_summary.read_text(encoding="utf-8"))
    run("summarize_freezeout_envelope.py","--screen-assessment",assessment,
        "--design",selection/"design.json","--standard-runs",case("standard"),
        "--fine-runs",case("fine_0"),case("fine_1"),
        "--baseline-runs",*[r["directory"] for r in old["baseline_inputs"]],
        "--experiment",args.experiment,"--other-family-summary",old_summary,"--output",case("summary"))
    run("plot_freezeout_thermodynamic_envelope.py","prepare","--summary-dir",case("summary"),
        "--output",case("Tmu_prepared"))
    snapshot=execution/Path(__file__).name
    if not snapshot.exists():shutil.copy2(__file__,snapshot)
    say({"status":"numerical_and_analytic_tables_complete","summary":str(case("summary")/"summary.json"),
        "pipeline_sha256":sha256_file(snapshot),"plotting_pending":True})


if __name__=="__main__":main()
