"""Prepare and assess a reproducible cross-literature freeze-out exploration.

The pooled coordinate box is an exploratory prior, not a confidence region.
Every candidate is one fixed five-parameter curve over the whole energy window.
The separate Julia runner computes all ratio points directly using q0 extrapolation.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import random
import shutil
import sys
import tomllib

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts.analysis.relaxtime.plot_freezeout_comparison import (
    load_calculation, load_experiment, residual_summary,
)
from scripts.plotting.plot_manifest import sha256_file

KEYS = ("a_GeV", "b_GeV_inv1", "c_GeV_inv3", "d_GeV", "e_GeV_inv1")
DEFAULT_CATALOG = ROOT / "config/physics/freezeout/parameter_domains/literature_envelope_v1.toml"
DEFAULT_EXPERIMENT = ROOT / "data/inputs/meson_freezeout_experiment_20261017/experimental_ratios.csv"


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def write_csv(path, rows):
    if not rows:
        raise ValueError("cannot write an empty table")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def source_parameters(source):
    values = list(map(float, source["values"]))
    if len(values) != 5 or not all(map(math.isfinite, values)):
        raise ValueError("source requires five finite coefficients")
    if source["kind"] == "normalized_quartic":
        t0, k2, k4, d, e = values
        if t0 <= 0:
            raise ValueError("T0 must be positive")
        return [t0, k2 / t0, k4 / t0**3, d, e]
    if source["kind"] != "polynomial":
        raise ValueError("unsupported source temperature family")
    return values


def pooled_bounds(catalog):
    lows, highs = [[] for _ in KEYS], [[] for _ in KEYS]
    for source in catalog["sources"]:
        values = source_parameters(source)
        errors = source.get("errors", [0.] * 5)
        if len(errors) != 5 or any(not math.isfinite(e) or e < 0 for e in errors):
            raise ValueError("source errors must be five finite nonnegative values")
        if source["kind"] != "polynomial" and any(errors):
            raise ValueError("normalized errors require a joint conversion, not coefficient substitution")
        for j, (v, error) in enumerate(zip(values, errors)):
            lows[j].append(v - error)
            highs[j].append(v + error)
    for partial in catalog.get("partial_constraints", []):
        lo, hi = partial["a_interval_GeV"]
        if not 0 < lo <= hi or not all(map(math.isfinite, (lo, hi))):
            raise ValueError("invalid partial a interval")
        lows[0].append(lo)
        highs[0].append(hi)
        if partial.get("b_zero_supported", False):
            lows[1].append(0.)
    bounds = [[min(lo), max(hi)] for lo, hi in zip(lows, highs)]
    if any(not math.isfinite(x) for pair in bounds for x in pair):
        raise ValueError("nonfinite pooled bound")
    return bounds


def interior_bounds(bounds, fraction):
    if not 0 <= fraction < .5:
        raise ValueError("trim fraction must be in [0, 0.5)")
    return [[lo + fraction*(hi-lo), hi - fraction*(hi-lo)] for lo, hi in bounds]


def trajectory(parameters, energy):
    if len(parameters) != 5 or not all(map(math.isfinite, parameters)):
        raise ValueError("five finite parameters required")
    if not math.isfinite(energy) or energy <= 0:
        raise ValueError("energy must be positive and finite")
    a, b, c, d, e = parameters
    if 1 + e*energy <= 0:
        raise ValueError("nonpositive chemical-potential denominator")
    mu = d / (1 + e*energy)
    return a - b*mu**2 - c*mu**4, mu


def trajectory_valid(parameters, domain):
    if len(parameters) != 5 or not all(map(math.isfinite, parameters)):
        return False
    a, b, c, d, e = parameters
    if not (a > 0 and b >= 0 and c >= 0 and d > 0 and e > 0):
        return False
    # For these signs mu decreases and T increases with energy. Endpoints
    # therefore certify positive T throughout this continuous interval.
    return all(trajectory(parameters, x)[0] > 0 for x in domain)


def in_box(parameters, bounds):
    return all(lo <= x <= hi for x, (lo, hi) in zip(parameters, bounds))


def verify_sources(catalog):
    records = {}
    for source in catalog["sources"] + catalog.get("partial_constraints", []) + catalog.get("context_constraints", []):
        path = Path(catalog["source_root"]) / source["pdf"]
        if sha256_file(path) != source["sha256"]:
            raise ValueError(f"source PDF hash mismatch: {path}")
        records[str(path)] = source["sha256"]
    return records


def initial_candidates(catalog):
    cfg = catalog["design"]
    bounds = pooled_bounds(catalog)
    inner = interior_bounds(bounds, cfg["interior_trim_fraction"])
    domain = [cfg["minimum_energy_GeV"], cfg["maximum_energy_GeV"]]
    candidates = []
    for source in catalog["sources"]:
        candidates.append(dict(name=source["id"], role="source_anchor", values=source_parameters(source),
            provenance=[source["id"]]))
    sources = {r["name"]: r["values"] for r in candidates}
    for a in (.147, .150):
        for parent in ("ihrg2019", "frg2024_muindep", "functional_qcd2026_mean"):
            p = [a, *sources[parent][1:]]
            candidates.append(dict(name=f"hybrid_a{round(a*1000)}_{parent}", role="explicit_hybrid",
                values=p, provenance=[parent, "bazavov2016_STAR20", "cross-source exploration"]))
    rng = random.Random(cfg["seed"])
    size = cfg["proposal_count"]
    columns = []
    for lo, hi in inner:
        strata = list(range(size))
        rng.shuffle(strata)
        columns.append([lo + (hi-lo)*(i+rng.random())/size for i in strata])
    rejected = []
    accepted = 0
    for index, values in enumerate(zip(*columns)):
        if not trajectory_valid(values, domain):
            rejected.append(dict(proposal=index, values=list(values), reason="nonpositive trajectory temperature"))
            continue
        accepted += 1
        candidates.append(dict(name=f"envelope_{accepted:03d}", role="interior_sample", values=list(values),
            provenance=["pooled coordinate intervals; stratified proposals with positivity rejection"]))
        if accepted == cfg["interior_sample_count"]:
            break
    if accepted != cfg["interior_sample_count"]:
        raise ValueError("insufficient valid candidates; increase proposal_count")
    for item in candidates:
        if not in_box(item["values"], bounds) or not trajectory_valid(item["values"], domain):
            raise ValueError(f"invalid source anchor or hybrid: {item['name']}")
    return candidates, bounds, inner, rejected


def save_design(output, catalog_path, catalog, candidates, bounds, inner, rejected, parent=None):
    if output.exists():
        raise ValueError("output exists; choose a new design directory")
    source_hashes = verify_sources(catalog)
    cfg = catalog["design"]
    output.mkdir(parents=True)
    profile_dir = output / "profiles"
    profile_dir.mkdir()
    for item in candidates:
        name = item["name"]
        tag = "literature_anchor:" + name if item["role"] == "source_anchor" else "cross_literature_exploration:" + name
        lines = [f'profile_name = "{name}"', 'family = "cleymans_polynomial_rational_v1"',
            f'source_tag = "{tag}"', '', '[coefficients]']
        lines += [f"{key} = {value!r}" for key, value in zip(KEYS, item["values"])]
        path = profile_dir / (name + ".toml")
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        item["profile_sha256"] = sha256_file(path)
        item["in_interior_box"] = in_box(item["values"], inner)
    result = dict(schema="freezeout_parameter_design_v1", created_utc=datetime.now(timezone.utc).isoformat(),
        catalog_sha256=sha256_file(catalog_path), source_pdf_hashes=source_hashes,
        parameter_order=list(KEYS), pooled_bounds=bounds, interior_bounds=inner,
        design=cfg, candidates=candidates, rejected_proposals=rejected,
        parent=parent, parameter_policy="literature_envelope_exploration",
        band_semantics="finite sampled parameter envelope; no statistical confidence level",
        physical_prior="a,d,e>0; b,c>=0; T>0 throughout declared domain; source-specific negative curvature not imported",
        generator_sha256=sha256_file(Path(__file__)),
        parts=[[c["name"] for c in candidates[i::2]] for i in range(2)])
    write_json(output / "design.json", result)
    write_csv(output / "parameters.csv", [dict(profile=c["name"], role=c["role"],
        **dict(zip(KEYS, c["values"])), in_interior_box=c["in_interior_box"]) for c in candidates])
    shutil.copy2(catalog_path, output / catalog_path.name)
    shutil.copy2(__file__, output / Path(__file__).name)
    print(json.dumps(dict(output=str(output), candidates=len(candidates), rejected=len(rejected),
        pooled_bounds=bounds, interior_bounds=inner), indent=2))


def load_designs(paths):
    designs, candidates = [], {}
    for path in paths:
        design = json.loads(path.read_text(encoding="utf-8"))
        designs.append(design)
        for c in design["candidates"]:
            if c["name"] in candidates:
                raise ValueError("duplicate profile across designs")
            if sha256_file(path.parent / "profiles" / (c["name"] + ".toml")) != c["profile_sha256"]:
                raise ValueError("generated profile changed")
            candidates[c["name"]] = c
    if any(d["catalog_sha256"] != designs[0]["catalog_sha256"] for d in designs):
        raise ValueError("parameter catalog changed across designs")
    return designs, candidates


def assess_models(models, experiments, candidates, energies, factor):
    observations = [r for r in experiments if r["primary_auau"] and r["energy"] in energies]
    if sorted(r["energy"] for r in observations) != sorted(energies):
        raise ValueError("exact fixed experimental energy set required")
    summaries = {}
    for name, points in models.items():
        failed = [e for e in energies if e not in points or not points[e]["passed"]]
        pairs = [(points[r["energy"]]["ratio"], r["ratio"]) for r in observations] if not failed else []
        summary = residual_summary(pairs)
        eight = [(points[r["energy"]]["ratio"], r["ratio"]) for r in observations if r["energy"] != 7.7] if not failed else []
        high = [(points[r["energy"]]["ratio"], r["ratio"]) for r in observations if r["energy"] >= 39] if not failed else []
        summary.update(role=candidates[name]["role"], failed_required_energies=failed,
            complete_primary=not failed, eight_point_rmse=residual_summary(eight)["rmse"],
            high_energy=residual_summary(high),
            ratio_7p7=points[7.7]["ratio"] if 7.7 in points and points[7.7]["passed"] else None,
            ratio_200=points[200.]["ratio"] if 200. in points and points[200.]["passed"] else None,
            failed_all_energies=[e for e, p in points.items() if not p["passed"]])
        horn_grid = sorted(e for e in points if 4.87 <= e <= 14.5)
        horn_complete = len(horn_grid) >= 3 and all(points[e]["passed"] for e in horn_grid)
        peak = max(horn_grid, key=lambda e: points[e]["ratio"]) if horn_complete else None
        summary["sampled_horn"] = dict(energy_GeV=peak,
            ratio=points[peak]["ratio"] if peak is not None else None,
            peak_at_grid_boundary=peak in (horn_grid[0], horn_grid[-1]) if peak is not None else None,
            complete=horn_complete, support_count=len(horn_grid))
        anchors = [4.87, 7.7, 11.5]
        actual = {r["energy"]: r["ratio"] for r in experiments if r["primary_auau"]}
        if all(e in points and points[e]["passed"] and e in actual for e in anchors):
            predicted = [points[e]["ratio"] for e in anchors]
            dr = [predicted[j+1]-predicted[j] for j in (0, 1)]
            observed = [actual[anchors[j+1]]-actual[anchors[j]] for j in (0, 1)]
            summary["horn_three_point_difference_rmse"] = math.sqrt(sum((x-y)**2 for x,y in zip(dr,observed))/2)
            summary["horn_rises_then_falls"] = dr[0] > 0 and dr[1] < 0
        else:
            summary["horn_three_point_difference_rmse"] = None
            summary["horn_rises_then_falls"] = None
        summaries[name] = summary
    ranking = sorted((n for n, s in summaries.items() if s["complete_primary"]), key=lambda n: (summaries[n]["rmse"], n))
    if not ranking:
        raise ValueError("no candidate has all fixed scoring observations")
    best = ranking[0]
    near = [n for n in ranking if summaries[n]["rmse"] <= factor*summaries[best]["rmse"]]
    return dict(metric="equal-weight RMSE at nine Au+Au observations, 7.7-200 GeV",
        primary_energies_GeV=energies, observation_count=len(energies), profiles=summaries,
        ranking=ranking, best_sampled=best, near_best_profiles=near, near_best_factor=factor,
        failures_policy="incomplete candidate is ineligible; never shrink the experimental set",
        optimization_scope="finite sampled exploration, not proof of a global optimum")


def complete_curve_band(models, names, energies):
    complete = [n for n in names if all(e in models[n] and models[n][e]["passed"] for e in energies)]
    rows = []
    for energy in energies:
        vals = [(models[n][energy]["ratio"], n) for n in complete]
        rows.append(dict(energy_GeV=energy, minimum=min(vals)[0] if vals else None,
            maximum=max(vals)[0] if vals else None, minimum_profile=min(vals)[1] if vals else None,
            maximum_profile=max(vals)[1] if vals else None, complete_curve_count=len(complete)))
    return rows, complete


def select_refined_curves(assessment, count=3):
    """Retain the RMSE leader and a near-best horn candidate before filling by RMSE."""
    if not isinstance(count, int) or count < 1:
        raise ValueError("refined curve count must be a positive integer")
    ranking = assessment["ranking"]
    selected = [ranking[0]]
    shape = [n for n in assessment["near_best_profiles"]
        if assessment["profiles"][n].get("horn_rises_then_falls")]
    if shape and count > 1:
        name = min(shape, key=lambda n: (assessment["profiles"][n]["horn_three_point_difference_rmse"], n))
        if name not in selected:
            selected.append(name)
    for name in ranking:
        if len(selected) >= count:
            break
        if name not in selected:
            selected.append(name)
    return selected


def remap_peak(parameters, original_energy, target_energy, new_e):
    """Map a sampled thermodynamic peak to another energy without changing T(muB)."""
    if not math.isfinite(new_e) or new_e <= 0 or not math.isfinite(target_energy) or target_energy <= 0:
        raise ValueError("new e and target energy must be positive and finite")
    _, mu_peak = trajectory(parameters, original_energy)
    return [*parameters[:3], mu_peak*(1+new_e*target_energy), new_e]


def assess_runs(args):
    if args.output.exists():
        raise ValueError("assessment output exists")
    designs, candidates = load_designs(args.designs)
    models, run_inputs, states = {}, [], []
    for run in args.runs:
        data, manifest, state = load_calculation(run)
        if state.get("parameter_policy") != "literature_envelope_exploration":
            raise ValueError("run does not declare parameter exploration")
        for name, points in data.items():
            if name in models or name not in candidates:
                raise ValueError("duplicate run profile or profile absent from design")
            c = candidates[name]
            for energy, point in points.items():
                t, mu = trajectory(c["values"], energy)
                if not math.isclose(1000*t, point["T_MeV"], abs_tol=1e-9, rel_tol=0) or not math.isclose(1000*mu, point["muB_MeV"], abs_tol=1e-9, rel_tol=0):
                    raise ValueError("calculated coordinates do not match the fixed parameter vector")
            models[name] = points
        states.append(state)
        run_inputs.append(dict(directory=str(run.resolve()), manifest_sha256=sha256_file(run / "manifest.json")))
    for state in states[1:]:
        for key in ("settings", "config", "background_contract", "coordinate_contract"):
            if state[key] != states[0][key]:
                raise ValueError(f"cannot pool runs with different {key}")
        shared = set(state["source_hashes"]) & set(states[0]["source_hashes"])
        if any(state["source_hashes"][p] != states[0]["source_hashes"][p] for p in shared):
            raise ValueError("shared numerical source changed")
    cfg = designs[0]["design"]
    result = assess_models(models, load_experiment(args.experiment), candidates,
        cfg["primary_energies_GeV"], cfg["near_best_rmse_factor"])
    band_energies = sorted(next(iter(models.values())))
    if any(sorted(points) != band_energies for points in models.values()):
        raise ValueError("all band candidates require the same directly calculated energy grid")
    band, complete = complete_curve_band(models, list(models), band_energies)
    near_band, near_complete = complete_curve_band(models, result["near_best_profiles"], band_energies)
    result.update(schema="freezeout_parameter_assessment_v1", run_inputs=run_inputs,
        design_inputs=[dict(path=str(p.resolve()), sha256=sha256_file(p)) for p in args.designs],
        experimental_csv=str(args.experiment.resolve()), experimental_csv_sha256=sha256_file(args.experiment),
        band_energies_GeV=band_energies, band_complete_profiles=complete,
        near_best_band_complete_profiles=near_complete, resolution=states[0]["settings"],
        band_semantics="min/max of the declared finite set of complete trajectories; not confidence bounds")
    args.output.mkdir(parents=True)
    write_json(args.output / "assessment.json", result)
    write_csv(args.output / "sampled_band.csv", [{**r, "near_minimum": n["minimum"], "near_maximum": n["maximum"],
        "near_complete_curve_count": n["complete_curve_count"]} for r, n in zip(band, near_band)])
    rows = []
    for name, points in models.items():
        for energy, p in sorted(points.items()):
            rows.append(dict(profile=name, energy_GeV=energy, T_MeV=p["T_MeV"], muB_MeV=p["muB_MeV"],
                ratio=p["ratio"] if p["passed"] else None, passed=p["passed"]))
    write_csv(args.output / "curves.csv", rows)
    shutil.copy2(__file__, args.output / Path(__file__).name)
    write_json(args.output / "checksums.json", {p.name: sha256_file(p) for p in sorted(args.output.iterdir()) if p.is_file()})
    print(json.dumps({"best_sampled":result["best_sampled"], "top": [(n,result["profiles"][n]["rmse"]) for n in result["ranking"][:8]],
        "near_best_count":len(result["near_best_profiles"]), "complete_band_count":len(complete)}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("prepare", "refine", "horn-map"):
        sub = commands.add_parser(name)
        sub.add_argument("--catalog", type=Path, default=DEFAULT_CATALOG)
        sub.add_argument("--output", type=Path, required=True)
        if name != "prepare":
            sub.add_argument("--assessment", type=Path, required=True)
            sub.add_argument("--designs", type=Path, nargs="+", required=True)
        if name == "horn-map":
            sub.add_argument("--target-energy", type=float, default=7.7)
            sub.add_argument("--map-e", type=float, default=.273)
    sub = commands.add_parser("assess")
    sub.add_argument("--designs", type=Path, nargs="+", required=True)
    sub.add_argument("--runs", type=Path, nargs="+", required=True)
    sub.add_argument("--experiment", type=Path, default=DEFAULT_EXPERIMENT)
    sub.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "assess":
        assess_runs(args)
        return
    catalog = tomllib.loads(args.catalog.read_text(encoding="utf-8"))
    candidates, bounds, inner, rejected = initial_candidates(catalog)
    parent = None
    if args.command != "prepare":
        designs, existing = load_designs(args.designs)
        assessment = json.loads(args.assessment.read_text(encoding="utf-8"))
        best = assessment["best_sampled"]
        base = existing[best]["values"]
        domain = [catalog["design"]["minimum_energy_GeV"], catalog["design"]["maximum_energy_GeV"]]
        candidates, rejected = [], []
        if args.command == "refine":
            for j, (lo, hi) in enumerate(bounds):
                for sign in (-1, 1):
                    values = base.copy()
                    values[j] = min(hi, max(lo, values[j] + sign*catalog["design"]["local_step_fraction"]*(hi-lo)))
                    if values == base:
                        continue
                    name = f"local_{KEYS[j][0]}_{'minus' if sign < 0 else 'plus'}"
                    item = dict(name=name, role="local_exploration", values=values, provenance=[best, "one predefined axial refinement round"])
                    if trajectory_valid(values, domain):
                        candidates.append(item)
                    else:
                        rejected.append(item)
        else:
            seen = set()
            for source_name in assessment["ranking"]:
                source_values = existing[source_name]["values"]
                abc = tuple(source_values[:3])
                peak = assessment["profiles"][source_name]["sampled_horn"]
                if abc in seen or not peak["complete"] or peak["peak_at_grid_boundary"]:
                    continue
                seen.add(abc)
                values = remap_peak(source_values, peak["energy_GeV"], args.target_energy, args.map_e)
                item = dict(name=f"hornmap_{len(candidates)+1:03d}",role="energy_mapping_exploration",values=values,
                    provenance=[source_name, "d=mu_at_sampled_peak*(1+e*target_energy); e is Cleymans central value"],
                    original_peak_energy_GeV=peak["energy_GeV"], target_energy_GeV=args.target_energy,
                    map_e_GeV_inv1=args.map_e)
                if in_box(values,bounds) and trajectory_valid(values,domain):
                    candidates.append(item)
                else:
                    rejected.append(item)
                if len(candidates) == 3:
                    break
            if not candidates:
                raise ValueError("no admissible energy-mapping candidates inside the declared envelope")
        parent = dict(assessment=str(args.assessment.resolve()), sha256=sha256_file(args.assessment), best=best)
    save_design(args.output, args.catalog, catalog, candidates, bounds, inner, rejected, parent)


if __name__ == "__main__":
    main()
