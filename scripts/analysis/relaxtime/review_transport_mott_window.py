#!/usr/bin/env python3
"""Review a frozen Mott window and channel changes on quote-remote."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from summarize_transport_mott_probe import analyse_probe, digest, read_csv, write_csv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("probe", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--left", required=True)
    parser.add_argument("--right", required=True)
    parser.add_argument("--comparison-probe", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    summary, crossings, comparisons, nonfinite, root_crossings = analyse_probe(args.probe)
    channels = read_csv(args.probe / "channel_contributions.csv")
    indexed = {(r["point_id"], r["species"], r["channel"]): r for r in channels}
    changes = []
    for species in ("u", "ubar", "s", "sbar"):
        for (point, q, channel), a in indexed.items():
            if point != args.left or q != species:
                continue
            b = indexed[args.right, species, channel]
            delta = float(b["inverse_tau_contribution_inv_fm"]) - float(a["inverse_tau_contribution_inv_fm"])
            total_delta = float(b["inverse_tau_inv_fm"]) - float(a["inverse_tau_inv_fm"])
            changes.append(dict(species=species, channel=channel, xi_left=float(a["xi"]), xi_right=float(b["xi"]),
                                rate_left_fm2=float(a["rate_fm2"]), rate_right_fm2=float(b["rate_fm2"]),
                                rate_change_percent=100 * (float(b["rate_fm2"]) / float(a["rate_fm2"]) - 1),
                                inverse_tau_change_inv_fm=delta, fraction_of_inverse_tau_change=delta / total_delta,
                                fraction_left=float(a["fraction"]), fraction_right=float(b["fraction"])))
    issues = [r for r in comparisons if not r["both_qualified"] or r["mass_difference_MeV"] > 1e-3
              or r["gamma_difference_MeV"] > 1e-3 or max(r["forward_roots"], r["reverse_roots"]) > 1]
    result = dict(schema="transport_mott_window_review_v1", status="diagnostic_only", host="quote-remote",
                  probe_summary_sha256=digest(args.probe / "summary.json"), root_points=summary["root_points"],
                  qualified_pole_rows=summary["qualified_pole_rows"], pole_rows=summary["pole_rows"],
                  threshold_crossings=crossings, root_threshold_crossings=root_crossings,
                  ambiguous_or_disagreeing_roots=issues, exact_endpoint_nonfinite_count=len(nonfinite),
                  max_bidirectional_mass_difference_MeV=max(r["mass_difference_MeV"] for r in comparisons),
                  max_bidirectional_gamma_difference_MeV=max(r["gamma_difference_MeV"] for r in comparisons),
                  top_ubar_channel_changes=sorted([r for r in changes if r["species"] == "ubar"],
                                                 key=lambda r: abs(r["inverse_tau_change_inv_fm"]), reverse=True)[:6])
    args.output.mkdir(parents=True)
    write_csv(args.output / "threshold_crossings.csv", crossings)
    write_csv(args.output / "root_threshold_crossings.csv", root_crossings)
    write_csv(args.output / "bidirectional_check.csv", comparisons)
    write_csv(args.output / "exact_endpoint_nonfinite.csv", nonfinite)
    write_csv(args.output / "channel_changes.csv", changes)
    if args.comparison_probe:
        other_summary, other_crossings, _, _, _ = analyse_probe(args.comparison_probe)
        baseline = {(r["point_id"], r["meson"], r["direction"]): r for r in read_csv(args.probe / "poles.csv")}
        refined = []
        for row in read_csv(args.comparison_probe / "poles.csv"):
            key = (row["point_id"], row["meson"], row["direction"])
            if key not in baseline:
                continue
            original = baseline[key]
            refined.append(dict(point_id=row["point_id"], meson=row["meson"], direction=row["direction"],
                                both_qualified=row["qualified"] == original["qualified"] == "true",
                                mass_change_MeV=float(row["mass_MeV"]) - float(original["mass_MeV"]),
                                gamma_change_MeV=float(row["gamma_MeV"]) - float(original["gamma_MeV"])))
        write_csv(args.output / "A_refinement_poles.csv", refined)
        result["A_refinement"] = dict(effective_A=other_summary["effective_A"],
                                      all_qualified=all(r["both_qualified"] for r in refined),
                                      max_mass_change_MeV=max(abs(r["mass_change_MeV"]) for r in refined),
                                      max_gamma_change_MeV=max(abs(r["gamma_change_MeV"]) for r in refined),
                                      threshold_brackets_unchanged=[(r["meson"], r["cut"], r["xi_left"], r["xi_right"]) for r in crossings]
                                      == [(r["meson"], r["cut"], r["xi_left"], r["xi_right"]) for r in other_crossings])
    (args.output / "review.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
