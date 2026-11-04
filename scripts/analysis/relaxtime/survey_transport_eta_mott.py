#!/usr/bin/env python3
"""Prepare and review an eta-pole survey on the accepted grid, on quote-remote."""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path

GRID_HASHES = {
    "mode_a": "4d812371ddd69678606fb0920813a274395865c1cd94c28f88cd26ef4a5c68c4",
    "mode_b": "5cddbeeb7944271071d08ce0a1ec6e5e32068bf09eaa3a9db473767df02ab2f8",
}
OBSERVABLES = ("tau_u", "tau_d", "tau_s", "tau_ubar", "tau_dbar", "tau_sbar",
               "eta_over_s", "zeta_over_s", "sigma_over_T")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def gap_after(row):
    """Only the two already reviewed first-order gaps split continuation."""
    x = float(row["xi"])
    return (row["mode_key"] == "mode_a" and float(row["muB_MeV"]) == 900
            and float(row["alpha_T"]) == 1 and x == -0.003) or (
                row["mode_key"] == "mode_b" and float(row["T_MeV"]) == 120
                and float(row["muB_MeV"]) == 900 and x == -0.13)


def prepare(args):
    args.output.mkdir(parents=True, exist_ok=False)
    rows, sources = [], []
    for mode, expected in GRID_HASHES.items():
        source = args.grid / mode / "phase_guided_transport_scan.csv"
        assert digest(source) == expected, f"changed accepted grid: {source}"
        sources.append(dict(path=str(source), sha256=expected, bytes=source.stat().st_size))
        groups = defaultdict(list)
        for line, row in enumerate(read_csv(source), 2):
            assert row["converged"] == "true", (source, line)
            assert all(math.isfinite(float(row[name])) and float(row[name]) > 0 for name in OBSERVABLES)
            groups[row["scan_group"]].append(dict(mode_key=mode, source_csv_line=line, **row))
        count = 0
        for group, curve in sorted(groups.items()):
            curve.sort(key=lambda row: float(row["xi"]))
            assert len({float(row["xi"]) for row in curve}) == len(curve)
            segment = 0
            for row in curve:
                count += 1
                window = f"{mode}:{group}"
                rows.append(dict(point_id=f"{mode[-1]}{count:04d}", window=window,
                                 continuation_segment=f"{window}:segment{segment}", **row))
                if gap_after(row):
                    segment += 1
    assert len(rows) == 1819
    assert len({r["window"] for r in rows}) == 18
    assert len({r["continuation_segment"] for r in rows}) == 20
    write_csv(args.output / "source_rows.csv", rows)
    manifest = dict(schema="transport_eta_mott_inputs_v1", host="quote-remote", sources=sources,
                    original_points=len(rows), curves=18, continuation_segments=20,
                    first_order_gaps=[dict(window=r["window"], xi_left=float(r["xi"])) for r in rows if gap_after(r)],
                    source_values_unchanged=True, background_recomputed=False,
                    source_rows_sha256=digest(args.output / "source_rows.csv"),
                    preparation_script_sha256=digest(Path(__file__)))
    write_json(args.output / "input_manifest.json", manifest)
    print(json.dumps(manifest))


def changes(left, right):
    result = dict(point_left=left["point_id"], point_right=right["point_id"],
                  xi_left=float(left["xi"]), xi_right=float(right["xi"]))
    for name in OBSERVABLES:
        result[name + "_change_percent"] = 100 * (float(right[name]) / float(left[name]) - 1)
    return result


def read_probe(probe, input_hash):
    summary = json.loads((probe / "summary.json").read_text(encoding="utf-8-sig"))
    assert summary["input_rows_sha256"] == input_hash
    for name, expected in summary["output_sha256"].items():
        assert digest(probe / name) == expected, (probe, name)
    poles = read_csv(probe / "poles.csv")
    assert len(poles) == summary["pole_rows"] == 2 * summary["root_points"]
    paired = defaultdict(dict)
    for row in poles:
        assert row["meson"] == "eta"
        assert row["direction"] not in paired[row["point_id"]]
        paired[row["point_id"]][row["direction"]] = row
    checks, trusted = [], {}
    for point, pair in paired.items():
        assert set(pair) == {"forward", "reverse"}
        f, r = pair["forward"], pair["reverse"]
        qualified = f["qualified"] == r["qualified"] == "true"
        dm = abs(float(f["mass_MeV"]) - float(r["mass_MeV"]))
        dg = abs(float(f["gamma_MeV"]) - float(r["gamma_MeV"]))
        dm = dm if math.isfinite(dm) else None
        dg = dg if math.isfinite(dg) else None
        multiple = max(int(f["unique_qualified_roots"]), int(r["unique_qualified_roots"])) > 1
        stable = qualified and dm is not None and dg is not None and dm < 1e-3 and dg < 1e-3 and not multiple
        checks.append(dict(point_id=point, window=f["window"], xi=float(f["xi"]),
                           both_qualified=qualified, mass_difference_MeV=dm,
                           gamma_difference_MeV=dg, multiple_roots=multiple, trusted=stable))
        if stable:
            trusted[point] = f
    return summary, poles, checks, trusted


def review(args):
    args.output.mkdir(parents=True, exist_ok=False)
    input_hash = digest(args.inputs / "source_rows.csv")
    manifest = json.loads((args.inputs / "input_manifest.json").read_text(encoding="utf-8-sig"))
    assert manifest["source_rows_sha256"] == input_hash
    rows = read_csv(args.inputs / "source_rows.csv")
    summary, poles, checks, trusted = read_probe(args.probe, input_hash)
    assert summary["root_points"] == len(rows)
    groups = defaultdict(list)
    for row in rows:
        groups[row["window"]].append(row)
    by_point = {r["point_id"]: r for r in rows}
    events, adjacency, local_steps, curve_ranges = [], [], [], []
    for window, curve in sorted(groups.items()):
        curve.sort(key=lambda row: float(row["xi"]))
        valid = [trusted[r["point_id"]] for r in curve if r["point_id"] in trusted]
        curve_ranges.append(dict(window=window, T_MeV=float(curve[0]["T_MeV"]),
            muB_MeV=float(curve[0]["muB_MeV"]), points=len(curve), trusted_points=len(valid),
            xi_min=float(curve[0]["xi"]), xi_max=float(curve[-1]["xi"]),
            light_gap_min_MeV=min((float(r["gap_light_MeV"]) for r in valid), default=None),
            light_gap_max_MeV=max((float(r["gap_light_MeV"]) for r in valid), default=None)))
        for index, (left, right) in enumerate(zip(curve, curve[1:])):
            is_gap = left["continuation_segment"] != right["continuation_segment"]
            adjacency.append(dict(window=window, first_order_gap=is_gap, **changes(left, right)))
            a, b = trusted.get(left["point_id"]), trusted.get(right["point_id"])
            if a is None or b is None:
                continue
            for cut in ("light", "strange"):
                gl, gr = float(a[f"gap_{cut}_MeV"]), float(b[f"gap_{cut}_MeV"])
                if gl * gr >= 0:
                    continue
                kind = "reverse" if gl > 0 else "forward"
                event_id = f"event{len(events)+1:02d}"
                event = dict(event_id=event_id, window=window, mode=left["mode_key"],
                    T_MeV=float(left["T_MeV"]), muB_MeV=float(left["muB_MeV"]),
                    alpha_T=float(left["alpha_T"]) if left["mode_key"] == "mode_a" else None,
                    cut=cut, direction=kind, first_order_gap=is_gap,
                    gap_left_MeV=gl, gap_right_MeV=gr,
                    gamma_left_MeV=float(a["gamma_MeV"]), gamma_right_MeV=float(b["gamma_MeV"]),
                    **changes(left, right))
                end = curve[min(index + 2, len(curve)-1)]
                if end["continuation_segment"] == left["continuation_segment"]:
                    event["two_step_xi_right"] = float(end["xi"])
                    for name in OBSERVABLES:
                        event["two_step_" + name + "_change_percent"] = 100 * (float(end[name]) / float(left[name]) - 1)
                else:
                    event["two_step_xi_right"] = None
                    for name in OBSERVABLES:
                        event["two_step_" + name + "_change_percent"] = None
                events.append(event)
                for j in range(max(0, index-4), min(len(curve)-1, index+6)):
                    l, r = curve[j], curve[j+1]
                    local_steps.append(dict(event_id=event_id, window=window, relative_step=j-index,
                        first_order_gap=l["continuation_segment"] != r["continuation_segment"], **changes(l, r)))
    # A signed real part at an open cut is a screen, not proof of a named Mott crossing.
    threshold_groups = defaultdict(list)
    threshold_nonfinite = []
    for row in read_csv(args.probe / "threshold_inverses.csv"):
        if not all(math.isfinite(float(row[k])) for k in ("inverse_re", "inverse_im")):
            threshold_nonfinite.append(row)
        if float(row["offset_inv_fm"]) == -1e-6:
            threshold_groups[row["window"], row["meson"], row["cut"]].append(row)
    screening = []
    for (window, meson, cut), values in sorted(threshold_groups.items()):
        values.sort(key=lambda r: float(r["xi"]))
        for a, b in zip(values, values[1:]):
            if float(a["inverse_re"]) * float(b["inverse_re"]) < 0:
                screening.append(dict(window=window, meson=meson, cut=cut,
                    point_left=a["point_id"], point_right=b["point_id"],
                    xi_left=float(a["xi"]), xi_right=float(b["xi"]),
                    real_left=float(a["inverse_re"]), real_right=float(b["inverse_re"]),
                    max_abs_imag=max(abs(float(a["inverse_im"])), abs(float(b["inverse_im"]))),
                    first_order_gap=by_point[a["point_id"]]["continuation_segment"] !=
                                    by_point[b["point_id"]]["continuation_segment"]))
    refinement = []
    if args.refinement:
        refined_summary, refined_poles, refined_checks, refined_trusted = read_probe(args.refinement, input_hash)
        for point, row in refined_trusted.items():
            if point not in trusted:
                continue
            base = trusted[point]
            refinement.append(dict(point_id=point, window=base["window"], xi=float(base["xi"]),
                mass_difference_MeV=float(row["mass_MeV"])-float(base["mass_MeV"]),
                gamma_difference_MeV=float(row["gamma_MeV"])-float(base["gamma_MeV"]),
                light_gap_sign_unchanged=float(row["gap_light_MeV"])*float(base["gap_light_MeV"]) > 0))
        assert len(refinement) == refined_summary["root_points"], "unresolved A-refinement roots"
    write_csv(args.output / "bidirectional_root_check.csv", checks)
    write_csv(args.output / "root_threshold_crossings.csv", events)
    write_csv(args.output / "adjacent_transport_changes.csv", adjacency)
    write_csv(args.output / "local_transport_changes.csv", local_steps)
    write_csv(args.output / "curve_gap_ranges.csv", curve_ranges)
    write_csv(args.output / "threshold_screening.csv", screening)
    write_csv(args.output / "threshold_nonfinite.csv", threshold_nonfinite)
    write_csv(args.output / "A_refinement.csv", refinement)
    report = dict(schema="transport_eta_mott_survey_review_v1", host="quote-remote",
        status="diagnostic_only", input_rows_sha256=input_hash,
        probe_summary_sha256=digest(args.probe / "summary.json"),
        original_points=len(rows), curves=len(groups), pole_rows=len(poles),
        trusted_points=len(trusted), unresolved_or_ambiguous=[c for c in checks if not c["trusted"]],
        threshold_nonfinite_count=len(threshold_nonfinite),
        max_bidirectional_mass_difference_MeV=max(c["mass_difference_MeV"] for c in checks if c["mass_difference_MeV"] is not None),
        max_bidirectional_gamma_difference_MeV=max(c["gamma_difference_MeV"] for c in checks if c["gamma_difference_MeV"] is not None),
        events=events, curve_ranges=curve_ranges,
        reverse_light_continuous_events=[e for e in events if e["cut"] == "light" and
                                         e["direction"] == "reverse" and not e["first_order_gap"]],
        A_refinement=refinement, transport_data="unchanged accepted full-RS rows; no new transport integrals",
        scope="eta algebraic light eigenbranch on original xi grid; does not establish exclusive eta causation or new integral convergence",
        output_sha256={p.name: digest(p) for p in args.output.iterdir() if p.is_file()})
    write_json(args.output / "review.json", report)
    print(json.dumps({key: report[key] for key in ("original_points", "curves", "trusted_points",
          "threshold_nonfinite_count", "events", "unresolved_or_ambiguous")}, allow_nan=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("prepare")
    p.add_argument("grid", type=Path)
    p.add_argument("output", type=Path)
    p.set_defaults(func=prepare)
    p = commands.add_parser("review")
    p.add_argument("inputs", type=Path)
    p.add_argument("probe", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--refinement", type=Path)
    p.set_defaults(func=review)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
