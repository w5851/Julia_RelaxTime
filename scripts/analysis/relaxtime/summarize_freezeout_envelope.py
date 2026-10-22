"""Freeze a matched-resolution review of selected exploratory curves."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from scripts.analysis.relaxtime.freezeout_parameter_envelope import (
    DEFAULT_EXPERIMENT, POLYNOMIAL_FAMILY, parameter_keys, assess_models, complete_curve_band, load_designs,
    trajectory, write_csv, write_json,
)
from scripts.analysis.relaxtime.plot_freezeout_comparison import (
    compare_resolutions, load_calculation, load_experiment,
)
from scripts.plotting.plot_manifest import sha256_file


def checked_runs(paths):
    merged, states, inputs = {}, [], []
    for path in paths:
        models, manifest, state = load_calculation(path)
        for relative, expected in state["source_hashes"].items():
            target = path / "source_snapshot" / relative
            if relative.startswith("../"):
                target = path / "source_snapshot/external_profiles" / Path(relative).name
            if sha256_file(target) != expected:
                raise ValueError(f"source snapshot hash mismatch: {target}")
        if set(merged) & set(models):
            raise ValueError("duplicate curve across run parts")
        merged.update(models)
        states.append(state)
        inputs.append(dict(directory=str(path.resolve()), scan_identity=manifest["scan_identity"],
            manifest_sha256=sha256_file(path / "manifest.json"),
            point_count=manifest["point_count"], successful_points=manifest["successful_points"],
            warning_points=manifest["warning_points"],
            snapshot_hashes_checked=len(state["source_hashes"])))
    return merged, states, inputs


def compatible_physics(reference, other, *, same_resolution):
    for field in ("config", "background_contract", "coordinate_contract", "density_route", "git_head"):
        if reference[field] != other[field]:
            raise ValueError(f"physics context changed: {field}")
    settings = ("qmax", "omega_lower_inv_fm", "endpoint_policy")
    if same_resolution:
        settings = reference["settings"]
    for key in settings:
        if reference["settings"][key] != other["settings"][key]:
            raise ValueError(f"numerical support changed: {key}")
    common = set(reference["source_hashes"]) & set(other["source_hashes"])
    changed = [k for k in common if reference["source_hashes"][k] != other["source_hashes"][k]]
    runner = "scripts/analysis/relaxtime/run_charged_gbu_freezeout_comparison.jl"
    if any(k != runner for k in changed):
        raise ValueError(f"shared physical source changed: {changed}")
    return dict(identical_source_count=len(common)-len(changed), changed=changed,
        allowed_runner_difference="explicit parameter-policy metadata and safe profile snapshot paths; older defaults also had a CLI default-energy change; all runs supply energies explicitly")


def curve_rows(models, minimum, maximum):
    return [dict(profile=name, energy_GeV=e, T_MeV=p["T_MeV"], muB_MeV=p["muB_MeV"],
        ratio=p["ratio"] if p["passed"] else None, passed=p["passed"])
        for name, points in models.items() for e, p in sorted(points.items()) if minimum <= e <= maximum]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--screen-assessment", type=Path, required=True)
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument("--standard-runs", type=Path, nargs="+", required=True)
    parser.add_argument("--fine-runs", type=Path, nargs="+", required=True)
    parser.add_argument("--baseline-runs", type=Path, nargs="+", required=True)
    parser.add_argument("--experiment", type=Path, default=DEFAULT_EXPERIMENT)
    parser.add_argument("--other-family-summary", type=Path,
        help="Optional frozen review to compare using the same fine physics and fixed observations")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("summary output exists")
    designs, candidates = load_designs([args.design])
    design = designs[0]
    family = design.get("family", POLYNOMIAL_FAMILY)
    keys = parameter_keys(family)
    experiments = load_experiment(args.experiment)
    standard, std_states, std_inputs = checked_runs(args.standard_runs)
    fine, fine_states, fine_inputs = checked_runs(args.fine_runs)
    baselines, old_states, old_inputs = checked_runs(args.baseline_runs)
    if set(standard) != set(fine):
        raise ValueError("standard/fine curve names differ")
    checks = [compatible_physics(fine_states[0], s, same_resolution=True) for s in fine_states+old_states]
    checks += [compatible_physics(fine_states[0], s, same_resolution=False) for s in std_states]
    checks += [compatible_physics(std_states[0], s, same_resolution=True) for s in std_states]
    for name, points in fine.items():
        for energy, point in points.items():
            t, mu = trajectory(candidates[name]["values"], energy, family)
            if abs(t*1000-point["T_MeV"]) > 1e-9 or abs(mu*1000-point["muB_MeV"]) > 1e-9:
                raise ValueError("fine coordinates differ from selected parameters")
    combined_meta = {**candidates, **{n:dict(role="published_default") for n in baselines}}
    cfg = design["design"]
    metrics = assess_models({**baselines, **fine}, experiments, combined_meta,
        cfg["primary_energies_GeV"], cfg["near_best_rmse_factor"])
    standard_metrics = assess_models(standard, experiments, candidates,
        cfg["primary_energies_GeV"], cfg["near_best_rmse_factor"])
    fine_ranking = [n for n in metrics["ranking"] if n in fine]
    baseline_ranking = [n for n in metrics["ranking"] if n in baselines]
    best = fine_ranking[0]
    score = metrics["profiles"][best]["rmse"]
    resolution = compare_resolutions(fine, standard)
    max_relative = max((r["max_relative_difference"]["relative_to_reference"] for r in resolution["profiles"].values()
        if r["max_relative_difference"] is not None), default=None)
    screen = json.loads(args.screen_assessment.read_text(encoding="utf-8"))
    for name, expected in json.loads((args.screen_assessment.parent / "checksums.json").read_text(encoding="utf-8")).items():
        if sha256_file(args.screen_assessment.parent / name) != expected:
            raise ValueError("screening assessment output changed")
    for record in screen["run_inputs"]:
        if sha256_file(Path(record["directory"]) / "manifest.json") != record["manifest_sha256"]:
            raise ValueError("screening assessment input changed")
    for record in screen["design_inputs"]:
        if sha256_file(Path(record["path"])) != record["sha256"]:
            raise ValueError("screening design changed")
    _, sampled_candidates = load_designs([Path(r["path"]) for r in screen["design_inputs"]])
    for name in fine:
        if sampled_candidates[name]["values"] != candidates[name]["values"]:
            raise ValueError("selected and screening parameters differ")
    near_parameters = [dict(profile=name, **dict(zip(keys, sampled_candidates[name]["values"])),
        screening_rmse=screen["profiles"][name]["rmse"])
        for name in screen["near_best_band_complete_profiles"]]
    fine_grid = sorted(next(iter(fine.values())))
    if any(sorted(p) != fine_grid for p in fine.values()):
        raise ValueError("selected fine curves do not share the dense grid")
    fine_band, complete = complete_curve_band(fine, list(fine), fine_grid)
    analysis_sources = [Path(__file__).resolve(),
        ROOT / "scripts/analysis/relaxtime/freezeout_parameter_envelope.py",
        ROOT / "scripts/analysis/relaxtime/plot_freezeout_comparison.py"]
    result = dict(schema="freezeout_envelope_review_v1", family=family, parameter_order=list(keys), best_selected=best,
        calculation_sha=fine_states[0]["git_head"],
        selected_parameters=dict(zip(keys, candidates[best]["values"])),
        selected_fine_ranking=fine_ranking, selected_standard_ranking=standard_metrics["ranking"],
        selected_ranking_stable=fine_ranking == standard_metrics["ranking"],
        selected_count=len(fine), fine_complete_band_members=complete,
        improvement_vs_andronic_fraction=1-score/metrics["profiles"]["andronic2018"]["rmse"],
        best_previous_default=baseline_ranking[0],
        improvement_vs_best_previous_default_fraction=1-score/metrics["profiles"][baseline_ranking[0]]["rmse"],
        maximum_standard_fine_relative_difference=max_relative,
        metrics=metrics, standard_metrics=standard_metrics, physics_checks=checks,
        standard_inputs=std_inputs, fine_inputs=fine_inputs, baseline_inputs=old_inputs,
        screen_assessment=dict(path=str(args.screen_assessment.resolve()), sha256=sha256_file(args.screen_assessment)),
        design_input=dict(path=str(args.design.resolve()), sha256=sha256_file(args.design)),
        sampled_design_inputs=screen["design_inputs"],
        near_best_parameter_sets=near_parameters,
        experiment=dict(path=str(args.experiment.resolve()), sha256=sha256_file(args.experiment)),
        analysis_sources=[dict(path=p.relative_to(ROOT).as_posix(), sha256=sha256_file(p)) for p in analysis_sources],
        limitations=["best among sampled/selected curves; no global optimum guarantee",
            "screening band has not been uniformly recomputed at fine resolution",
            "source-specific covariance relaxed by author-approved exploration; no confidence level",
            "low-energy extensions do not reproduce each source paper's validity range",
            "fixed finite integration window; no endpoint/tail limit certification",
            "no final-state decay or rapidity-acceptance model"])
    if args.other_family_summary is not None:
        other_path = args.other_family_summary.resolve()
        for name, expected in json.loads((other_path.parent / "checksums.json").read_text(encoding="utf-8")).items():
            if sha256_file(other_path.parent / name) != expected:
                raise ValueError("other-family frozen review changed")
        other = json.loads(other_path.read_text(encoding="utf-8"))
        if (other["experiment"]["sha256"] != result["experiment"]["sha256"] or
                other["metrics"]["primary_energies_GeV"] != cfg["primary_energies_GeV"]):
            raise ValueError("other-family observations differ")
        other_models, other_states, other_inputs = checked_runs([Path(r["directory"]) for r in other["fine_inputs"]])
        other_checks = [compatible_physics(fine_states[0],s,same_resolution=True) for s in other_states]
        other_metrics = assess_models(other_models,experiments,
            {n:dict(role="other_family_selected") for n in other_models},
            cfg["primary_energies_GeV"],cfg["near_best_rmse_factor"])
        other_best = other["best_selected"]
        other_rmse = other_metrics["profiles"][other_best]["rmse"]
        if other_rmse != other["metrics"]["profiles"][other_best]["rmse"]:
            raise ValueError("other-family fixed-point score cannot be reproduced")
        result["other_family_comparison"] = dict(summary_path=str(other_path),sha256=sha256_file(other_path),
            family=other.get("family",POLYNOMIAL_FAMILY),best_selected=other_best,
            selected_parameters=other["selected_parameters"],rmse=other_rmse,
            this_family_rmse=score,relative_rmse_change=score/other_rmse-1,
            fine_inputs=other_inputs,physics_checks=other_checks,
            interpretation="Matched fine resolution and fixed nine observations; finite sampled family minima, not global optima")
    args.output.mkdir(parents=True)
    write_json(args.output / "summary.json", result)
    write_json(args.output / "resolution_comparison.json", resolution)
    write_csv(args.output / "fine_curves.csv", curve_rows(fine, cfg["minimum_energy_GeV"], cfg["maximum_energy_GeV"]))
    write_csv(args.output / "baseline_curves.csv", curve_rows(baselines, cfg["minimum_energy_GeV"], cfg["maximum_energy_GeV"]))
    write_csv(args.output / "fine_band.csv", fine_band)
    shutil.copy2(args.screen_assessment.parent / "sampled_band.csv", args.output / "screen_band.csv")
    write_csv(args.output / "selected_parameters.csv", [dict(profile=n, **dict(zip(keys,candidates[n]["values"]))) for n in fine])
    write_csv(args.output / "near_best_parameters.csv", near_parameters)
    for path in analysis_sources:
        target = args.output / "source_snapshot" / path.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
    write_json(args.output / "checksums.json", {p.relative_to(args.output).as_posix():sha256_file(p)
        for p in sorted(args.output.rglob("*")) if p.is_file()})
    print(json.dumps({k:result[k] for k in ("best_selected","selected_parameters","improvement_vs_andronic_fraction",
        "best_previous_default","improvement_vs_best_previous_default_fraction","maximum_standard_fine_relative_difference",
        "selected_ranking_stable")}, indent=2))


if __name__ == "__main__":
    main()
