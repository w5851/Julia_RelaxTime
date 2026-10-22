"""Design the Andronic-family counterpart of the frozen polynomial exploration.

Scoring, missing-point policy and band membership reuse the shared analysis.
This script prepares inputs only; the existing Julia runner evaluates densities.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import random
import shutil
import sys
import tomllib

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from scripts.analysis.relaxtime.freezeout_parameter_envelope import (
    ANDRONIC_FAMILY, ANDRONIC_KEYS, DEFAULT_EXPERIMENT, assess_runs,
    in_box, interior_bounds, load_designs, select_refined_curves,
    trajectory, trajectory_valid, write_csv, write_json,
)
from scripts.plotting.plot_manifest import sha256_file

DEFAULT_CONFIG = ROOT / "config/physics/freezeout/parameter_domains/andronic_exploration_v1.toml"


def read_config(path):
    config = tomllib.loads(path.read_text(encoding="utf-8"))
    domain_path = ROOT / config["domain"]
    if sha256_file(domain_path) != config["domain_sha256"]:
        raise ValueError("selected Andronic domain changed")
    domain = tomllib.loads(domain_path.read_text(encoding="utf-8"))
    if domain["family"] != ANDRONIC_FAMILY or domain["parameter_order"] != list(ANDRONIC_KEYS):
        raise ValueError("Andronic family and parameter order required")
    if config["design"]["primary_energies_GeV"] != domain["comparison"]["primary_energies_GeV"]:
        raise ValueError("fixed comparison energies changed")
    return config, domain


def initial_candidates(config, domain):
    cfg = config["design"]
    bounds = [domain["bounds"][k] for k in ANDRONIC_KEYS]
    inner = interior_bounds(bounds, cfg["interior_trim_fraction"])
    candidates, seen = [], set()

    def add(name, role, values, provenance):
        if tuple(values) in seen:
            return
        if not in_box(values, bounds) or not trajectory_valid(values,
                [cfg["minimum_energy_GeV"], cfg["maximum_energy_GeV"]], ANDRONIC_FAMILY):
            raise ValueError(f"invalid Andronic candidate: {name}")
        seen.add(tuple(values))
        candidates.append(dict(name=name, role=role, values=list(values), provenance=provenance))

    baseline = next(s["values"] for s in domain["sources"] if s["id"] == "andronic2018")
    for source in domain["sources"]:
        if source["values"] == baseline:
            continue  # Reuse existing 2018 fine/standard results as the reference.
        add(source["id"], "source_anchor", source["values"], [source["id"]])
    rng = random.Random(cfg["seed"])
    size = cfg["interior_sample_count"]
    columns = []
    for lo, hi in inner:
        strata = list(range(size))
        rng.shuffle(strata)
        columns.append([lo+(hi-lo)*(i+rng.random())/size for i in strata])
    for i, values in enumerate(zip(*columns), 1):
        add(f"andronic_sample_{i:03d}", "interior_sample", values,
            ["stratified sampling in the middle 80% of each selected coordinate interval"])
    if cfg["boundary_axis_samples"]:
        center = [(lo+hi)/2 for lo, hi in bounds]
        for j, (key, pair) in enumerate(zip(ANDRONIC_KEYS, bounds)):
            for side, value in zip(("low", "high"), pair):
                values = center.copy()
                values[j] = value
                add(f"andronic_axis_{key}_{side}", "axis_boundary_sample", values,
                    [f"{key} at {side} bound; other four coordinates at box center"])
    return candidates, bounds, inner


def rescale_energy(parameters, original_energy, target_energy):
    """Preserve T(muB) exactly while mapping one energy to another.

    k=x_old/x_new, E'=k E, A'=A-ln(k)/w; Tlim,w,D stay fixed.
    Changing D alone, as in the polynomial family, would change this geometry.
    """
    if not all(math.isfinite(x) and x > 0 for x in (original_energy, target_energy)):
        raise ValueError("original and target energies must be positive and finite")
    tlim, offset, width, d, e = parameters
    if width <= 0:
        raise ValueError("log width must be positive")
    factor = original_energy/target_energy
    return [tlim, offset-math.log(factor)/width, width, d, factor*e]


def source_evidence(domain):
    inputs, pdfs = [], {}
    for field in ("source_manifest", "previous_source_manifest"):
        path = ROOT / domain[field]
        inputs.append(dict(path=str(path.resolve()), sha256=sha256_file(path)))
    path = ROOT / domain["source_manifest"]
    manifest = json.loads(path.read_text(encoding="utf-8"))
    for source in manifest["references"]:
        pdf = path.parent / source["local_pdf"]
        if pdf.stat().st_size != source["bytes"] or sha256_file(pdf) != source["sha256"]:
            raise ValueError(f"source PDF changed: {pdf}")
        pdfs[str(pdf.resolve())] = source["sha256"]
    return inputs, pdfs


def save_design(output, config_path, config, domain, candidates, parent=None, rejected=None):
    if output.exists() or not candidates:
        raise ValueError("choose a new output directory and at least one candidate")
    bounds = [domain["bounds"][k] for k in ANDRONIC_KEYS]
    inner = interior_bounds(bounds, config["design"]["interior_trim_fraction"])
    evidence, pdfs = source_evidence(domain)
    output.mkdir(parents=True)
    profiles = output / "profiles"
    profiles.mkdir()
    candidates = deepcopy(candidates)
    for item in candidates:
        name = item["name"]
        text = [f'profile_name = "{name}"', f'family = "{ANDRONIC_FAMILY}"',
            f'source_tag = "andronic_literature_envelope_v2:{name}"', "", "[coefficients]"]
        text += [f"{key} = {value!r}" for key, value in zip(ANDRONIC_KEYS, item["values"])]
        path = profiles / (name+".toml")
        path.write_text("\n".join(text)+"\n", encoding="utf-8")
        item["profile_sha256"] = sha256_file(path)
        item["in_interior_box"] = in_box(item["values"], inner)
    result = dict(schema="freezeout_parameter_design_v1", family=ANDRONIC_FAMILY,
        created_utc=datetime.now(timezone.utc).isoformat(), catalog_sha256=sha256_file(config_path),
        domain_input=dict(path=str((ROOT/config["domain"]).resolve()), sha256=config["domain_sha256"]),
        source_evidence=evidence, source_pdf_hashes=pdfs,
        parameter_order=list(ANDRONIC_KEYS), pooled_bounds=bounds, interior_bounds=inner,
        design=config["design"], candidates=candidates, rejected_proposals=rejected or [],
        reused_reference="andronic2018; Elfner 2023 duplicates it; JHEP 2022 provides no D/E",
        parent=parent, parameter_policy="literature_envelope_exploration",
        band_semantics="finite sampled complete-curve envelope; no statistical confidence level",
        physical_prior="Tlim,w,D,E>0; analytic T rises and muB decreases on the declared positive-energy interval",
        generator_sha256=sha256_file(Path(__file__)),
        parts=[[c["name"] for c in candidates[i::2]] for i in range(2)])
    write_json(output / "design.json", result)
    write_csv(output / "parameters.csv", [dict(profile=c["name"], role=c["role"],
        **dict(zip(ANDRONIC_KEYS,c["values"])), in_interior_box=c["in_interior_box"]) for c in candidates])
    for path in [config_path, ROOT/config["domain"], Path(__file__),
            ROOT/"scripts/analysis/relaxtime/freezeout_parameter_envelope.py"]:
        shutil.copy2(path, output/path.name)
    print(json.dumps(dict(output=str(output), candidates=len(candidates), rejected=len(rejected or []))))
    return result


def checked_parent(args, config_path):
    designs, existing = load_designs(args.designs)
    if any(d.get("family") != ANDRONIC_FAMILY or d["catalog_sha256"] != sha256_file(config_path) for d in designs):
        raise ValueError("parent design uses another family or catalog")
    assessment = json.loads(args.assessment.read_text(encoding="utf-8"))
    for name, expected in json.loads((args.assessment.parent/"checksums.json").read_text(encoding="utf-8")).items():
        if sha256_file(args.assessment.parent/name) != expected:
            raise ValueError("parent assessment changed")
    supplied = {str(p.resolve()):sha256_file(p) for p in args.designs}
    if supplied != {r["path"]:r["sha256"] for r in assessment["design_inputs"]}:
        raise ValueError("assessment and parent design set differ")
    return assessment, existing


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("prepare", "refine", "horn-map", "select"):
        sub = commands.add_parser(name)
        sub.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
        sub.add_argument("--output", type=Path, required=True)
        if name != "prepare":
            sub.add_argument("--assessment", type=Path, required=True)
            sub.add_argument("--designs", type=Path, nargs="+", required=True)
    sub = commands.add_parser("assess")
    sub.add_argument("--designs", type=Path, nargs="+", required=True)
    sub.add_argument("--runs", type=Path, nargs="+", required=True)
    sub.add_argument("--experiment", type=Path, default=DEFAULT_EXPERIMENT)
    sub.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "assess":
        assess_runs(args)
        return
    config, domain = read_config(args.config)
    initial, bounds, _ = initial_candidates(config, domain)
    if args.command == "prepare":
        save_design(args.output,args.config,config,domain,initial)
        return
    assessment, existing = checked_parent(args,args.config)
    best = assessment["best_sampled"]
    parent = dict(assessment=str(args.assessment.resolve()), sha256=sha256_file(args.assessment), best=best)
    candidates, rejected = [], []
    seen = {tuple(c["values"]) for c in existing.values()}
    if args.command == "select":
        names = select_refined_curves(assessment,config["design"]["refined_curve_count"])
        candidates = [existing[n] for n in names]
    elif args.command == "refine":
        base = existing[best]["values"]
        for j, (lo,hi) in enumerate(bounds):
            for sign in (-1,1):
                values = base.copy()
                values[j] = min(hi,max(lo,values[j]+sign*config["design"]["local_step_fraction"]*(hi-lo)))
                if tuple(values) in seen:
                    continue
                seen.add(tuple(values))
                candidates.append(dict(name=f"andronic_local_{ANDRONIC_KEYS[j]}_{'minus' if sign<0 else 'plus'}",
                    role="local_exploration", values=values, provenance=[best,"one axial refinement round"]))
    else:
        for name in assessment["ranking"]:
            peak = assessment["profiles"][name]["sampled_horn"]
            if not peak["complete"] or peak["peak_at_grid_boundary"]:
                continue
            values = rescale_energy(existing[name]["values"],peak["energy_GeV"],7.7)
            if tuple(values) in seen:
                continue
            seen.add(tuple(values))
            item = dict(name=f"andronic_hornmap_{len(candidates)+1:03d}",role="energy_mapping_exploration",
                values=values,provenance=[name,"E'=k*E; A'=A-ln(k)/w; k=x_peak/7.7; exact T(muB) preservation"],
                original_peak_energy_GeV=peak["energy_GeV"],target_energy_GeV=7.7)
            if in_box(values,bounds):
                candidates.append(item)
            else:
                rejected.append(item)
            if len(candidates)==config["design"]["maximum_mapping_candidates"]:
                break
        if not candidates:
            print(json.dumps(dict(mapping_candidates=0,reason="No exact energy rescaling remains in the authorized five-dimensional box",rejected=rejected)))
            return
    save_design(args.output,args.config,config,domain,candidates,parent,rejected)


if __name__ == "__main__":
    main()
