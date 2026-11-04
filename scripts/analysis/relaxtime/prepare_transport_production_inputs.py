"""Prepare reviewable runner pilot and provisional mode-A shards after local qualification."""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess

OBSERVABLES = {"tau_u", "tau_d", "tau_s", "tau_ubar", "tau_dbar", "tau_sbar",
               "eta_over_s", "zeta_over_s", "sigma_over_T"}
SOURCE_ROOT = Path(__file__).resolve().parents[3]
NUMERICAL_PREFIXES = ("src/", "config/models/pnjl/", "config/physics/")


def is_numerical_source(relative):
    path = Path(relative)
    return ((relative.startswith(NUMERICAL_PREFIXES) and path.suffix in (".jl", ".toml"))
            or (path.parent.as_posix() == "scripts/relaxtime" and path.suffix == ".jl")
            or relative in ("Project.toml", "Manifest.toml"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def verified_sources(directory, manifest, source_root):
    """Require phase and replay evidence to describe the code being prepared."""
    sources = {}
    for entry in manifest["source_files"]:
        relative = entry["path"].replace("\\", "/")
        # Qualification snapshots include the shared scan scripts directly in
        # scripts/relaxtime; its separate workflow subdirectories are unused.
        if not is_numerical_source(relative):
            continue
        if relative in sources:
            raise ValueError("duplicate numerical source in qualification")
        frozen = (directory / "source" / relative).resolve()
        current = (source_root / relative).resolve()
        if (not frozen.is_relative_to((directory / "source").resolve()) or not frozen.is_file()
                or digest(frozen) != entry["sha256"]):
            raise ValueError(f"qualification source snapshot changed: {relative}")
        if not current.is_relative_to(source_root) or not current.is_file() or digest(current) != entry["sha256"]:
            raise ValueError(f"qualification numerical source differs from current code: {relative}")
        sources[relative] = entry["sha256"]
    expected = {p.relative_to(source_root).as_posix()
                for prefix in NUMERICAL_PREFIXES for p in (source_root / prefix).rglob("*")
                if p.is_file() and p.suffix in (".jl", ".toml")}
    expected.update(p.relative_to(source_root).as_posix()
                    for p in (source_root / "scripts/relaxtime").glob("*.jl") if p.is_file())
    expected.update(name for name in ("Project.toml", "Manifest.toml") if (source_root / name).is_file())
    if not sources or set(sources) != expected:
        raise ValueError("incomplete numerical source identity")
    return sources


def verified_calculation_commit(calculation_sha, source_root, sources):
    """Bind dispatch to the qualified code, allowing Git's normal EOL conversion.

    An arbitrary well-formed SHA is not a calculation identity. Check both the
    numerical file set and its committed content; unrelated docs are excluded.
    """
    def git(*args):
        try:
            return subprocess.check_output(["git", "-C", str(source_root), *args],
                encoding="utf-8", stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as error:
            raise ValueError("calculation commit is unavailable in this checkout") from error

    if git("cat-file", "-t", calculation_sha).strip() != "commit":
        raise ValueError("calculation SHA must identify a commit")
    committed = {p for p in git("ls-tree", "-r", "-z", "--name-only", calculation_sha).split("\0")
                 if p and is_numerical_source(p)}
    if committed != set(sources):
        raise ValueError("calculation commit does not contain the qualified numerical file set")
    changed = {p for p in git("diff", "--no-ext-diff", "--no-renames", "--name-only", "-z",
        calculation_sha, "--", "src", "config/models/pnjl", "config/physics",
        "scripts/relaxtime", "Project.toml", "Manifest.toml").split("\0")
        if p and is_numerical_source(p)}
    if changed:
        raise ValueError("calculation commit differs from qualified numerical code: " +
                         ", ".join(sorted(changed)[:5]))


def verified_reference(directory, manifest):
    path = directory / "inputs" / "qualified_comparisons.csv"
    checksum = digest(path)
    if checksum != manifest["parameters"]["reference_comparisons_sha256"]:
        raise ValueError("qualified comparison snapshot changed")
    targets = {}
    for row in read_csv(path):
        key = (row["point_id"], row["observable"])
        values = tuple(float(row[k]) for k in ("candidate", "reference"))
        if key in targets or key[1] not in OBSERVABLES or not all(math.isfinite(v) and v > 0 for v in values):
            raise ValueError("invalid qualified comparison set")
        targets[key] = values
    points = {point for point, _ in targets}
    if len(points) != 12 or set(targets) != {(point, name) for point in points for name in OBSERVABLES}:
        raise ValueError("qualified reference must contain twelve complete points")
    return checksum, targets


def point_key(row):
    values = tuple(float(row[k]) for k in ("T_MeV", "muB_MeV", "xi"))
    if not all(math.isfinite(v) for v in values):
        raise ValueError("nonfinite point coordinate")
    return values


def verified_phase_plan(directory):
    status = json.loads((directory / "status.json").read_text(encoding="utf-8"))
    if not status.get("phase_plan_passed") or not status.get("complete") or status.get("failures") != 0:
        raise ValueError("new-background phase plan has not passed")
    plan = read_csv(directory / "sampling_plan.csv")
    expected = {point_key(r) for r in plan}
    if not plan or len(expected) != len(plan) or status.get("points") != len(plan) or status.get("completed_points") != len(plan):
        raise ValueError("incomplete or duplicate phase plan")
    anchors = read_csv(directory / "anchors.csv")
    if {float(r["muB_MeV"]) for r in anchors} != {float(r["muB_MeV"]) for r in plan}:
        raise ValueError("phase anchor set does not match the plan")
    for row in anchors:
        delta = float(row["refinement_delta_MeV"])
        if row["certified"].lower() != "true" or not math.isfinite(delta) or abs(delta) > 0.1:
            raise ValueError("uncertified phase anchor")
    background = read_csv(directory / "background.csv")
    if len(background) != len(expected) or {point_key(r) for r in background} != expected:
        raise ValueError("background qualification does not cover the plan")
    for row in background:
        mass_delta, entropy_delta = float(row["max_mass_relative_change"]), float(row["entropy_relative_change"])
        same_phase = (row["distinct_branches"].lower() != "true" or row["refined_distinct_branches"].lower() != "true"
                      or row["stable_phase"] == row["refined_stable_phase"])
        if (row["passed"].lower() != "true" or not same_phase
                or not math.isfinite(mass_delta) or not 0 <= mass_delta <= 0.01
                or not math.isfinite(entropy_delta) or abs(entropy_delta) > 0.01):
            raise ValueError("background numerical qualification failed")
    return plan


def verify_comparison(row):
    if row["observable"] not in OBSERVABLES or row["passed"].lower() != "true":
        raise ValueError("invalid replay comparison")
    production, candidate, reference = (float(row[k]) for k in ("production", "qualified_candidate", "qualified_reference"))
    if not all(math.isfinite(v) and v > 0 for v in (production, candidate, reference)):
        raise ValueError("invalid replay observable")
    replay, refined = production / candidate - 1, production / reference - 1
    if abs(replay) > 1e-4 or abs(refined) > 0.01:
        raise ValueError("replay exceeds numerical acceptance limits")
    for key, actual in (("replay_relative_change", replay), ("reference_relative_change", refined)):
        recorded = float(row[key])
        if not math.isfinite(recorded) or not math.isclose(actual, recorded, rel_tol=1e-10, abs_tol=1e-14):
            raise ValueError("recorded replay error differs from actual values")


def prepare(phase_dir, replay_dirs, outdir, case_name, calculation_sha=None, shard_points=20, *, source_root=None):
    if case_name in (".", "..") or not re.fullmatch(r"[A-Za-z0-9_.=-]+", case_name):
        raise ValueError("invalid case name")
    if calculation_sha is not None and not re.fullmatch(r"[0-9a-f]{40}", calculation_sha):
        raise ValueError("calculation SHA must be full, immutable and lowercase")
    if shard_points < 1:
        raise ValueError("shard_points must be positive")
    plan = verified_phase_plan(phase_dir)
    phase_manifest = json.loads((phase_dir / "manifest.json").read_text(encoding="utf-8"))
    profile = phase_manifest["parameters"]["profile"]
    root = (source_root or SOURCE_ROOT).resolve()
    phase_sources = verified_sources(phase_dir, phase_manifest, root)
    compared = set()
    replay_evidence = []
    reference_checksum, reference_targets = None, None
    for directory in replay_dirs:
        status = json.loads((directory / "status.json").read_text(encoding="utf-8"))
        if not status.get("production_replay_passed") or not status.get("complete") or status.get("failures") != 0:
            raise ValueError("production replay has not passed")
        replay_manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
        replay_profile = replay_manifest["parameters"]["profile"]
        if (replay_profile["id"], replay_profile["source_sha256"]) != (profile["id"], profile["source_sha256"]):
            raise ValueError("replay and phase qualification use different sampling profiles")
        if verified_sources(directory, replay_manifest, root) != phase_sources:
            raise ValueError("replay and phase qualification use different numerical sources")
        checksum, targets = verified_reference(directory, replay_manifest)
        if reference_checksum is not None and checksum != reference_checksum:
            raise ValueError("replay shards use different qualified references")
        reference_checksum, reference_targets = checksum, targets
        comparisons = read_csv(directory / "comparisons.csv")
        points = {row["point_id"] for row in comparisons}
        if (not comparisons or status.get("comparisons") != len(comparisons)
                or status.get("completed_points") != len(points) or status.get("points") != len(points)):
            raise ValueError("invalid replay comparison")
        for row in comparisons:
            verify_comparison(row)
            key = (row["point_id"], row["observable"])
            if key not in targets or (float(row["qualified_candidate"]), float(row["qualified_reference"])) != targets[key]:
                raise ValueError("replay point or target differs from the qualified reference")
            if key in compared:
                raise ValueError("duplicate replay comparison")
            compared.add(key)
        replay_evidence.append({"directory": str(directory.resolve()), "comparisons_sha256": digest(directory / "comparisons.csv")})
    point_ids = {p for p, _ in compared}
    if len(point_ids) != 12 or any({name for p, name in compared if p == point} != OBSERVABLES for point in point_ids):
        raise ValueError("all twelve replay points and 108 comparisons are required")
    if compared != set(reference_targets):
        raise ValueError("replay does not cover the qualified reference points")
    if calculation_sha is not None:
        verified_calculation_commit(calculation_sha, root, phase_sources)
    groups = defaultdict(list)
    for row in plan:
        groups[(float(row["muB_MeV"]), float(row["alpha_T"]))].append(row)
    for rows in groups.values():
        rows.sort(key=lambda r: float(r["xi"]))

    def jobs_for(selected, stage):
        grouped = defaultdict(list)
        for row in selected:
            grouped[(float(row["muB_MeV"]), float(row["alpha_T"]))].append(row)
        jobs = []
        for (mu, alpha), rows in sorted(grouped.items()):
            rows.sort(key=lambda r: float(r["xi"]))
            for start in range(0, len(rows), shard_points):
                chunk = rows[start:start + shard_points]
                label = f"{stage}_muB{mu:g}_alpha{alpha:g}_{start // shard_points:02d}"
                inputs = {"case_name": case_name, "shard_label": label,
                          "sampling_profile": profile["id"], "mode": "fixed-muB-phase-scaled",
                          "muB_list": f"{mu:g}", "alpha_t_list": f"{alpha:g}",
                          "xi_list": ",".join(r["xi"] for r in chunk), "t_list": "200",
                          "run_tier": "production", "render_plots": False}
                if calculation_sha is not None:
                    inputs["calculation_sha"] = calculation_sha
                jobs.append({"label": label, "inputs": inputs, "expected_points": len(chunk),
                             "expected_T_MeV": float(chunk[0]["T_MeV"])})
        return jobs

    pilot_keys = [(0, 1, -.22), (0, 1, .37), (450, 1.1, .35), (900, 1, -.04), (0, 1.2, -.5)]
    selected = []
    for mu, alpha, xi in pilot_keys:
        matches = [r for r in groups[(mu, alpha)] if abs(float(r["xi"]) - xi) < 1e-10]
        if len(matches) != 1:
            raise ValueError(f"pilot point missing: {mu}, {alpha}, {xi}")
        selected.append(matches[0])
    side = [r for r in groups[(900, 1)] if r["coexistence_side"] == "quark_side"]
    if len(side) != 1:
        raise ValueError("certified first-order side point missing")
    selected.append(side[0])
    payload = {"schema_version": 1, "case_name": case_name, "sampling_profile": profile["id"],
               "sampling_profile_sha256": profile["source_sha256"], "calculation_sha": calculation_sha,
               "numerical_source_sha256": hashlib.sha256(json.dumps(phase_sources, sort_keys=True).encode()).hexdigest(),
               "reference_comparisons_sha256": reference_checksum,
               "phase_source": str(phase_dir.resolve()), "sampling_plan_sha256": digest(phase_dir / "sampling_plan.csv"),
               "replay_evidence": replay_evidence, "pilot_jobs": jobs_for(selected, "pilot"),
               "provisional_full_jobs": jobs_for(plan, "full"), "full_points": len(plan),
               "pilot_input_ready": calculation_sha is not None, "full_dispatch_ready": False,
               "remaining_before_dispatch": ["select the commit containing the reviewed implementation"] if calculation_sha is None else [],
               "full_shard_policy": "provisional; calibrate startup, complete point time and peak memory on the six-point runner pilot before full dispatch"}
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("use a new empty output directory")
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "runner_inputs.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return payload


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("phase_dir", type=Path)
    p.add_argument("outdir", type=Path)
    p.add_argument("--replay-dir", type=Path, action="append", required=True)
    p.add_argument("--case-name", required=True)
    p.add_argument("--calculation-sha")
    p.add_argument("--shard-points", type=int, default=20)
    args = p.parse_args()
    result = prepare(args.phase_dir, args.replay_dir, args.outdir, args.case_name, args.calculation_sha, args.shard_points)
    print(json.dumps({k: result[k] for k in ("full_points", "pilot_input_ready", "full_dispatch_ready")}, indent=2))


if __name__ == "__main__":
    main()
