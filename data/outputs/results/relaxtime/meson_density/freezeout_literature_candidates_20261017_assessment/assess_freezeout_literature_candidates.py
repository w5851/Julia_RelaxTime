"""Assess frozen literature-default trajectories on a declared Au+Au subset.

The eight-point comparison excludes 7.7 GeV because the FRG paper excluded it
from its parameter fit. The separately retained 7.7 GeV predictions do not
establish a horn position. No interpolation, fit, or figure production occurs.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.analysis.relaxtime.plot_freezeout_comparison import (
    compare_resolutions, load_calculation, load_experiment, residual_summary,
)
from scripts.plotting.plot_manifest import sha256_file

ENERGIES = (11.5, 14.5, 19.6, 27.0, 39.0, 62.4, 130.0, 200.0)
FRG_PROFILES = {"frg2024_mudep", "frg2024_muindep"}


def assess(models, experiments):
    observations = [r for r in experiments if r["primary_auau"] and r["energy"] in ENERGIES]
    if sorted(r["energy"] for r in observations) != list(ENERGIES):
        raise ValueError("the declared eight Au+Au experimental energies are required")
    common, excluded = [], []
    for row in observations:
        energy = row["energy"]
        if any(energy not in points for points in models.values()):
            raise ValueError(f"missing directly calculated experimental energy: {energy}")
        failed = [name for name, points in models.items() if not points[energy]["passed"]]
        if failed:
            excluded.append({"id": row["id"], "energy_GeV": energy, "failed_profiles": failed})
        else:
            common.append(row)
    if not common:
        raise ValueError("no common valid observations")
    summaries, residuals = {}, []
    for name, points in models.items():
        pairs = [(points[r["energy"]]["ratio"], r["ratio"]) for r in common]
        summary = residual_summary(pairs)
        summary["ratio_at_7p7"] = points[7.7]["ratio"] if points[7.7]["passed"] else None
        summary["ratio_at_200"] = points[200.0]["ratio"] if points[200.0]["passed"] else None
        summary["7p7_outside_source_fit"] = name in FRG_PROFILES
        summaries[name] = summary
        for row in common:
            point = points[row["energy"]]
            residuals.append({"profile": name, "experiment_id": row["id"],
                "energy_GeV": row["energy"], "T_MeV": point["T_MeV"],
                "muB_MeV": point["muB_MeV"], "model_ratio": point["ratio"],
                "experimental_ratio": row["ratio"], "residual": point["ratio"] - row["ratio"]})
    baseline = summaries["andronic2018"]["rmse"]
    for summary in summaries.values():
        summary["rmse_reduction_vs_andronic_fraction"] = 1 - summary["rmse"] / baseline
    return {"energies_GeV": [r["energy"] for r in common], "excluded_failures": excluded,
        "profiles": summaries, "ranking": sorted(summaries, key=lambda n: summaries[n]["rmse"])}, residuals


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("standard-run", "fine-run", "baseline-standard-run", "baseline-fine-run", "output"):
        parser.add_argument("--" + key, type=Path, required=True)
    parser.add_argument("--experiment", type=Path, default=ROOT / "data/inputs/meson_freezeout_experiment_20261017/experimental_ratios.csv")
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("output exists; use a new assessment directory")
    experiments = load_experiment(args.experiment)
    batches = {}
    for key in ("standard_run", "fine_run", "baseline_standard_run", "baseline_fine_run"):
        batches[key] = load_calculation(getattr(args, key))
    standard, _, state = batches["standard_run"]
    fine, _, refined_state = batches["fine_run"]
    for key in ("source_hashes", "config", "background_contract", "coordinate_contract", "points"):
        if state[key] != refined_state[key]:
            raise ValueError(f"candidate batches changed {key}")
    for key in ("qmax", "omega_lower_inv_fm", "endpoint_policy"):
        if state["settings"][key] != refined_state["settings"][key]:
            raise ValueError(f"integration support changed: {key}")
    verified_shared = {}
    for label in ("standard", "fine"):
        current = batches[label + "_run"][2]
        old = batches["baseline_" + label + "_run"][2]
        for key in ("settings", "config", "background_contract", "coordinate_contract"):
            if current[key] != old[key]:
                raise ValueError(f"baseline differs in {key}: {label}")
        shared = set(current["source_hashes"]) & set(old["source_hashes"])
        runner = "scripts/analysis/relaxtime/run_charged_gbu_freezeout_comparison.jl"
        # The earlier frozen runner predates the corrected DEFAULT_ENERGIES.
        # Both experiments explicitly supplied energies; preserve both snapshots.
        changed = [p for p in shared if current["source_hashes"][p] != old["source_hashes"][p]]
        if any(p != runner for p in changed):
            raise ValueError(f"shared numerical input changed: {changed}")
        verified_shared[label] = {"identical_source_paths": sorted(shared - set(changed)),
            "changed_runner_path": changed,
            "runner_difference": "earlier frozen runner predates CLI default-energy correction; both used explicit energy lists"}
    results = {}
    fine_residuals = []
    for label in ("standard", "fine"):
        models = {**batches["baseline_" + label + "_run"][0], **batches[label + "_run"][0]}
        results[label], residuals = assess(models, experiments)
        if label == "fine":
            fine_residuals = residuals
    result = {"schema": "literature_freezeout_candidate_assessment_v1",
        "metric": "equal-weight RMSE on common Au+Au observations at 11.5-200 GeV",
        "scope": "restricted-energy diagnostic; not a full-energy ranking or horn-shape validation",
        "parameter_adjustment": False, "density_route": "q0_lambda_reference",
        "experimental_csv_sha256": sha256_file(args.experiment),
        "comparisons": results, "resolution_comparison": compare_resolutions(fine, standard),
        "baseline_source_checks": verified_shared,
        "source_limitations": {"frg2024": "source uses muI=muS=0, without strangeness neutrality; 7.7 GeV excluded from fit",
            "functional_qcd2026": "preprint; T0 and kappa2 constrained to source chiral-line values; uses mean row only",
            "comparison_background": "all imported T,muB trajectories evaluated with project quark-only BQS nS=0,nQ/nB=0.4"},
        "inputs": {key: {"directory": str(getattr(args, key).resolve()),
            "manifest_sha256": sha256_file(getattr(args, key) / "manifest.json")} for key in batches}}
    args.output.mkdir(parents=True)
    (args.output / "assessment.json").write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    with (args.output / "residuals.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fine_residuals[0]))
        writer.writeheader()
        writer.writerows(fine_residuals)
    shutil.copy2(__file__, args.output / Path(__file__).name)
    print(json.dumps(result["comparisons"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
