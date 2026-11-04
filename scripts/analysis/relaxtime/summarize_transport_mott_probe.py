"""Audit retained local Mott probes and existing channel rates on quote-remote."""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def analyse_probe(probe):
    summary = json.loads((probe / "summary.json").read_text(encoding="utf-8-sig"))
    for name, expected in summary["output_sha256"].items():
        assert digest(probe / name) == expected, (probe, name)
    threshold_rows = read_csv(probe / "threshold_inverses.csv")
    curves = defaultdict(list)
    nonfinite_endpoints = []
    for row in threshold_rows:
        finite = all(math.isfinite(float(row[k])) for k in ("inverse_re", "inverse_im", "inverse_abs"))
        if not finite:
            assert float(row["offset_inv_fm"]) == 0, "nonfinite result away from exact branch point"
            nonfinite_endpoints.append({key:row[key] for key in ("point_id", "meson", "cut", "offset_inv_fm")})
        # The bare width helper can be undefined at the exact branch point.
        # Preserve those rows; use a fixed one-sided evaluation for all points.
        if float(row["offset_inv_fm"]) == -1e-6:
            assert finite
            curves[row["window"], row["meson"], row["cut"]].append(row)
    crossings = []
    for (window, meson, cut), rows in sorted(curves.items()):
        rows.sort(key=lambda row: float(row["xi"]))
        for left, right in zip(rows, rows[1:]):
            fl, fr = float(left["inverse_re"]), float(right["inverse_re"])
            if fl * fr < 0:
                crossings.append(dict(window=window, meson=meson, cut=cut,
                    xi_left=float(left["xi"]), xi_right=float(right["xi"]),
                    point_left=left["point_id"], point_right=right["point_id"],
                    inverse_re_left=fl, inverse_re_right=fr,
                    evaluation_offset_inv_fm=-1e-6,
                    max_endpoint_inverse_im=max(abs(float(left["inverse_im"])), abs(float(right["inverse_im"]))),
                    classification="threshold_real_part_zero_bracket"))
    poles = read_csv(probe / "poles.csv")
    by_key = {(r["point_id"], r["meson"], r["direction"]): r for r in poles}
    comparisons = []
    for (point, meson, direction), forward in by_key.items():
        if direction != "forward":
            continue
        reverse = by_key[point, meson, "reverse"]
        comparisons.append(dict(point_id=point, window=forward["window"], xi=float(forward["xi"]),
            meson=meson, both_qualified=forward["qualified"] == reverse["qualified"] == "true",
            mass_difference_MeV=abs(float(forward["mass_MeV"])-float(reverse["mass_MeV"])),
            gamma_difference_MeV=abs(float(forward["gamma_MeV"])-float(reverse["gamma_MeV"])),
            forward_roots=int(forward["unique_qualified_roots"]),
            reverse_roots=int(reverse["unique_qualified_roots"])))
    root_curves = defaultdict(list)
    for r in poles:
        if r["direction"] == "forward":
            root_curves[r["window"],r["meson"]].append(r)
    root_crossings = []
    for (window,meson), curve in root_curves.items():
        curve.sort(key=lambda r:float(r["xi"]))
        cuts = (("light", "gap_light_MeV"),("strange","gap_strange_MeV")) if meson in (
            "eta","eta_prime","sigma","sigma_prime") else (("pair","gap_MeV"),)
        for cut,key in cuts:
            for a,b in zip(curve,curve[1:]):
                if a["qualified"] == b["qualified"] == "true" and float(a[key])*float(b[key]) < 0:
                    root_crossings.append(dict(window=window,meson=meson,cut=cut,
                        xi_left=float(a["xi"]),xi_right=float(b["xi"]),
                        gap_left_MeV=float(a[key]),gap_right_MeV=float(b[key]),
                        gamma_left_MeV=float(a["gamma_MeV"]),gamma_right_MeV=float(b["gamma_MeV"])))
    return summary, crossings, comparisons, nonfinite_endpoints, root_crossings


def channel_changes(probe):
    rows = read_csv(probe / "channel_contributions.csv")
    indexed = {(r["point_id"], r["species"], r["channel"]): r for r in rows}
    changes = []
    for left, right in (("p004", "p005"), ("p014", "p015")):
        for species in ("u", "ubar", "s", "sbar"):
            for key, a in indexed.items():
                if key[:2] != (left, species):
                    continue
                b = indexed[right, species, key[2]]
                changes.append(dict(window=a["window"], species=species, channel=a["channel"],
                    xi_left=float(a["xi"]), xi_right=float(b["xi"]),
                    rate_left_fm2=float(a["rate_fm2"]), rate_right_fm2=float(b["rate_fm2"]),
                    rate_relative_change=float(b["rate_fm2"])/float(a["rate_fm2"])-1,
                    contribution_left_inv_fm=float(a["inverse_tau_contribution_inv_fm"]),
                    contribution_right_inv_fm=float(b["inverse_tau_contribution_inv_fm"]),
                    contribution_change_inv_fm=float(b["inverse_tau_contribution_inv_fm"])-float(a["inverse_tau_contribution_inv_fm"]),
                    total_change_inv_fm=float(b["inverse_tau_inv_fm"])-float(a["inverse_tau_inv_fm"]),
                    fraction_of_total_change=(float(b["inverse_tau_contribution_inv_fm"])-float(a["inverse_tau_contribution_inv_fm"])) /
                        (float(b["inverse_tau_inv_fm"])-float(a["inverse_tau_inv_fm"]))))
    return sorted(changes, key=lambda r: (r["window"], r["species"], r["contribution_change_inv_fm"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("probe", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--comparison-probe", type=Path)
    args = parser.parse_args()
    assert not args.output.exists(), "refusing to overwrite existing review"
    args.output.mkdir(parents=True)
    summary, crossings, comparisons, nonfinite, root_crossings = analyse_probe(args.probe)
    write_csv(args.output / "threshold_crossing_brackets.csv", crossings)
    write_csv(args.output / "bidirectional_root_check.csv", comparisons)
    write_csv(args.output / "exact_endpoint_nonfinite.csv", nonfinite)
    write_csv(args.output / "root_threshold_crossings.csv", root_crossings)
    changes = channel_changes(args.probe)
    write_csv(args.output / "channel_changes.csv", changes)
    reviewed = dict(probe=str(args.probe), probe_summary=summary,
        status="analysis_only", core_or_production_changes=False,
        threshold_crossings=crossings,
        root_threshold_crossings=root_crossings,
        exact_endpoint_nonfinite_count=len(nonfinite),
        threshold_evaluation="one-sided M_threshold - 1e-6 fm^-1; exact nonfinite rows retained",
        ambiguous_or_disagreeing_roots=[r for r in comparisons if not r["both_qualified"] or
            r["mass_difference_MeV"] > 1e-3 or r["gamma_difference_MeV"] > 1e-3 or
            max(r["forward_roots"],r["reverse_roots"]) > 1],
        top_channel_changes=[r for r in changes if r["species"] in ("u", "ubar") and
            abs(r["contribution_change_inv_fm"]) >= 0.01*abs(r["total_change_inv_fm"])])
    forward = [r for r in read_csv(args.probe/"poles.csv") if r["direction"]=="forward"]
    ranges = []
    for window in sorted({r["window"] for r in forward}):
        for meson in ("eta","eta_prime","sigma","sigma_prime"):
            rows = [r for r in forward if r["window"]==window and r["meson"]==meson]
            ranges.append(dict(window=window,meson=meson,
                min_light_gap_MeV=min(float(r["gap_light_MeV"]) for r in rows),
                max_light_gap_MeV=max(float(r["gap_light_MeV"]) for r in rows),
                min_strange_gap_MeV=min(float(r["gap_strange_MeV"]) for r in rows),
                max_strange_gap_MeV=max(float(r["gap_strange_MeV"]) for r in rows)))
    write_csv(args.output/"mixed_threshold_gap_ranges.csv",ranges)
    reviewed["mixed_threshold_gap_ranges"] = ranges
    reviewed["max_bidirectional_mass_difference_MeV"] = max(r["mass_difference_MeV"] for r in comparisons)
    reviewed["max_bidirectional_gamma_difference_MeV"] = max(r["gamma_difference_MeV"] for r in comparisons)
    if args.comparison_probe:
        other_summary, other_crossings, other_comparisons, other_nonfinite, other_roots = analyse_probe(args.comparison_probe)
        primary = {(r["point_id"],r["meson"],r["direction"]):r for r in read_csv(args.probe/"poles.csv")}
        refinement = []
        for r in read_csv(args.comparison_probe/"poles.csv"):
            b = primary.get((r["point_id"],r["meson"],r["direction"]))
            if b is None or r["direction"] != "forward":
                continue
            refinement.append(dict(point_id=r["point_id"],meson=r["meson"],
                both_qualified=b["qualified"]==r["qualified"]=="true",
                mass_change_MeV=float(r["mass_MeV"])-float(b["mass_MeV"]),
                gamma_change_MeV=float(r["gamma_MeV"])-float(b["gamma_MeV"])))
        write_csv(args.output/"A_refinement_poles.csv",refinement)
        reviewed["A_refinement"]=dict(summary=other_summary,threshold_crossings=other_crossings,
            pole_changes=refinement,
            max_mass_change_MeV=max(abs(r["mass_change_MeV"]) for r in refinement),
            max_gamma_change_MeV=max(abs(r["gamma_change_MeV"]) for r in refinement),
            threshold_brackets_unchanged=
                [(r["window"],r["meson"],r["cut"],r["xi_left"],r["xi_right"]) for r in crossings]==
                [(r["window"],r["meson"],r["cut"],r["xi_left"],r["xi_right"]) for r in other_crossings])
    reviewed["files"] = {p.name:digest(p) for p in args.output.iterdir() if p.is_file()}
    (args.output/"review.json").write_text(json.dumps(reviewed,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(dict(threshold_crossings=crossings,
        qualified_poles=summary["qualified_pole_rows"],pole_rows=summary["pole_rows"],
        ambiguous_or_disagreeing=len(reviewed["ambiguous_or_disagreeing_roots"]))))


if __name__ == "__main__":
    main()
