"""Select mode-B transport checks from a qualified, immutable background grid."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import tomllib

ROOT = Path(__file__).resolve().parents[3]
TEMPERATURES = (120, 160, 200)
CHEMICAL_POTENTIALS = (0, 450, 900)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def key(row):
    return tuple(float(row[name]) for name in ("T_MeV", "muB_MeV", "xi"))


def expected_keys():
    return {(t, mu, i / 100) for t in TEMPERATURES for mu in CHEMICAL_POTENTIALS for i in range(-50, 51)}


def validate_background(plan, background, status):
    expected = expected_keys()
    if (not status.get("complete") or not status.get("phase_plan_passed")
            or not status.get("cold_seed_audits_passed") or status.get("failures") != 0
            or status.get("points") != 909 or status.get("completed_points") != 909):
        raise ValueError("mode-B background qualification is incomplete")
    if (len(plan) != 909 or len(background) != 909 or {key(r) for r in plan} != expected
            or {key(r) for r in background} != expected
            or any(r["mode"] != "mode_b_fixed_T_sparse_muB" for r in plan)):
        raise ValueError("mode-B background must cover the exact requested 909-point grid")
    for row in background:
        if row["passed"].lower() != "true" or row["coexistence_resolved"].lower() != "true":
            raise ValueError("failed or unresolved mode-B background point")
        masses = [float(row[name]) for name in ("m_u", "m_d", "m_s", "entropy")]
        mass_change, entropy_change = float(row["max_mass_relative_change"]), float(row["entropy_relative_change"])
        if (not all(math.isfinite(v) and v > 0 for v in masses)
                or not math.isfinite(mass_change) or not 0 <= mass_change <= .01
                or not math.isfinite(entropy_change) or abs(entropy_change) > .01):
            raise ValueError("mode-B background numerical check failed")
        if (row["distinct_branches"] == row["refined_distinct_branches"] == "true"
                and row["stable_phase"] != row["refined_stable_phase"]):
            raise ValueError("mode-B background phase choice changed with refinement")


def select_points(plan, boundaries):
    """Cover all nine curves, anisotropy extremes, isotropic controls and phase brackets."""
    available = {key(row) for row in plan}
    selected = {}

    def add(point, reason):
        if point not in available:
            raise ValueError(f"requested qualification point is missing: {point}")
        selected.setdefault(point, []).append(reason)

    for t in TEMPERATURES:
        for mu in CHEMICAL_POTENTIALS:
            add((t, mu, -.5), "negative anisotropy endpoint on each requested curve")
    for t, mu in ((120, 0), (120, 900), (200, 900)):
        add((t, mu, .5), "opposite anisotropy endpoint at a domain corner")
    for t, mu in ((120, 900), (160, 450), (200, 0)):
        add((t, mu, 0), "isotropic control across temperature and chemical potential")
    seen = set()
    for row in boundaries:
        t, mu = float(row["T_MeV"]), float(row["muB_MeV"])
        left, right = float(row["bracket_lower_xi"]), float(row["bracket_upper_xi"])
        identity = t, mu, left, right
        if identity in seen or not left < float(row["xi"]) < right or not left < float(row["refined_xi"]) < right:
            raise ValueError("invalid or duplicated mode-B phase boundary bracket")
        seen.add(identity)
        for xi in (left, right):
            add((t, mu, xi), "qualified grid endpoint bracketing first-order coexistence")
    return selected


def point_id(point):
    t, mu, xi = point
    return f"T{t:.1f}_muB{mu:.1f}_xi{xi:+.4f}"


def stable_trace(forward, backward, background):
    high, low = sorted((forward, backward), key=lambda row: float(row["m_u"]), reverse=True)
    stable = high if float(high["omega"]) <= float(low["omega"]) else low
    if tuple(float(stable[name]) for name in ("m_u", "m_d", "m_s")) != tuple(
            float(background[name]) for name in ("m_u", "m_d", "m_s")):
        raise ValueError("stable trace does not reproduce the qualified background masses")
    if not all(math.isfinite(float(stable[name])) for name in ("Phi", "Phibar", "omega")):
        raise ValueError("nonfinite retained background seed")
    return stable


def write_csv(path, data):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)


def prepare(directory, outdir, *, source_root=ROOT):
    source_root = source_root.resolve()
    spec = importlib.util.spec_from_file_location("production_inputs", source_root / "scripts/analysis/relaxtime/prepare_transport_production_inputs.py")
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    manifest, status = read_json(directory / "manifest.json"), read_json(directory / "status.json")
    source_hashes = helper.verified_sources(directory, manifest, source_root)
    plan, background = rows(directory / "sampling_plan.csv"), rows(directory / "background.csv")
    validate_background(plan, background, status)
    boundary_path = directory / "phase_boundaries.csv"
    boundaries = rows(boundary_path) if boundary_path.exists() else []
    if len(boundaries) != status["first_order_boundaries"]:
        raise ValueError("mode-B boundary count does not match the qualified background")
    chosen = select_points(plan, boundaries)
    background_index = {key(row): row for row in background}
    trace_files, trace_index = {}, {}
    for group in {row["scan_group"] for row in plan}:
        for direction in ("forward", "backward"):
            path = directory / "traces" / f"{group}_candidate_{direction}.csv"
            data = rows(path)
            if len(data) != 101 or len({key(row) for row in data}) != 101:
                raise ValueError("candidate continuation trace is incomplete")
            trace_files[str(path.resolve())] = digest(path)
            trace_index[group, direction] = {key(row): row for row in data}
    hbarc = tomllib.loads((source_root / "config/physics/default.toml").read_text(encoding="utf-8"))["physical"]["hbarc"]
    review = []
    for point, reasons in sorted(chosen.items()):
        bg = background_index[point]
        group = bg["scan_group"]
        seed = stable_trace(trace_index[group, "forward"][point], trace_index[group, "backward"][point], bg)
        t, mu, xi = point
        review.append(dict(point_id=point_id(point), T_MeV=t, muB_MeV=mu, T_fm=t/hbarc, muq_fm=mu/(3*hbarc), xi=xi,
                           **{name: float(seed[name]) for name in ("m_u", "m_d", "m_s", "Phi", "Phibar")},
                           scan_group=group, stable_phase=bg["stable_phase"], selection_reason="; ".join(reasons)))
    kernel = [(point_id((120, 0, .5)), "ssbar_to_ssbar"),
              (point_id((120, 900, -.5)), "uubar_to_uubar"),
              (point_id((200, 900, -.5)), "usbar_to_usbar")]
    for row in boundaries:
        kernel.extend((point_id((float(row["T_MeV"]), float(row["muB_MeV"]), float(row[side]))), "uubar_to_uubar")
                      for side in ("bracket_lower_xi", "bracket_upper_xi"))
    kernel = sorted(set(kernel))
    outer = [(point_id(point), "udbar_to_udbar") for point in ((120, 0, .5), (120, 900, -.5), (200, 900, -.5))]
    selection = dict(mode="B", points=[row["point_id"] for row in review],
                     pilot_point_id=point_id((120, 900, -.5)), kernel_cases=kernel, outer_cases=outer,
                     background_directory=str(directory.resolve()), background_manifest_sha256=digest(directory / "manifest.json"),
                     background_status_sha256=digest(directory / "status.json"), sampling_plan_sha256=digest(directory / "sampling_plan.csv"),
                     background_csv_sha256=digest(directory / "background.csv"), boundary_sha256=digest(boundary_path) if boundaries else None,
                     trace_files=trace_files, numerical_sources=source_hashes,
                     preparation_script_sha256=digest(Path(__file__)),
                     source_archive=manifest.get("source_archive"),
                     scope="Selected mode-B transport points; the 909-point background grid alone does not qualify transport",
                     transport_qualified=False, formal_production_promoted=False)
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("use a new empty input directory")
    outdir.mkdir(parents=True, exist_ok=True)
    write_csv(outdir / "review_points.csv", review)
    # The calculation uses the canonical native channel list; this file records
    # only selected diagnostic channels and does not redefine that list.
    channel_plan = [dict(point_id=point, channel=channel, stage=stage)
                    for stage, cases in (("kernel", kernel), ("outer", outer)) for point, channel in cases]
    write_csv(outdir / "channel_plan.csv", channel_plan)
    selection["review_points_sha256"] = digest(outdir / "review_points.csv")
    selection["channel_plan_sha256"] = digest(outdir / "channel_plan.csv")
    (outdir / "selection.json").write_text(json.dumps(selection, indent=2) + "\n", encoding="utf-8")
    return selection


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("background_dir", type=Path)
    parser.add_argument("outdir", type=Path)
    parser.add_argument("--source-root", type=Path, default=ROOT)
    args = parser.parse_args()
    result = prepare(args.background_dir, args.outdir, source_root=args.source_root)
    print(json.dumps({name: result[name] for name in ("points", "pilot_point_id", "kernel_cases", "outer_cases")}, indent=2))
