#!/usr/bin/env python3
"""Audit the eta-window convergence evidence; numerical analysis runs remotely."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from review_transport_mott_sampling_convergence import digest, finite, read_csv, read_json, verify_job, write_csv

POINTS = ("x04p013", "x04p014", "x04p015")
CHANNELS = ("uubar_to_uubar", "uubar_to_ddbar")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    expected = {(p, c) for p in POINTS for c in CHANNELS}
    studies = {}
    for run in ("s_v1", "outer_v1"):
        path = args.case / run
        summary = read_json(path / "summary.json")
        if digest(path / "channels.csv") != summary["channels_sha256"]:
            raise ValueError(f"aggregate hash differs: {run}")
        rows = read_csv(path / "channels.csv")
        if len(rows) != len(expected) or {(r["point_id"], r["channel"]) for r in rows} != expected:
            raise ValueError(f"incomplete or duplicated point/channel coverage: {run}")
        for row in rows:
            verify_job(path / row["job_id"])
        studies[run] = rows
    s_rows, outer = studies["s_v1"], studies["outer_v1"]
    reference_passed = all(max(finite(r, "reference_previous_l1"), finite(r, "reference_shifted_l1")) <= 2e-5
                           and abs(finite(r, "reproduction_relative")) <= 1e-6 for r in s_rows)
    outer_passed = all(r["outer_converged"] == "true" and abs(finite(r, "joint2_relative_to_joint1")) <= 1e-3
                       and finite(r, "max_s_l1") <= 2e-5 for r in outer)
    path = args.case / "propagation_v1"
    ps = read_json(path / "summary.json")
    for name, value in ps["files"].items():
        if digest(path / name) != value:
            raise ValueError(f"propagation hash differs: {name}")
    for run, value in ps["input_hashes"].items():
        if digest(args.case / run / "channels.csv") != value:
            raise ValueError(f"propagation input differs: {run}")
    replay = read_csv(path / "formal_replay.csv")
    if len(replay) != len(POINTS) or any(r["passed"] != "true" for r in replay):
        raise ValueError("formal replay failed")
    tau = read_csv(path / "relaxation_times.csv")
    observables = read_csv(path / "observables.csv")
    steps = []
    maxima = {}
    for configuration in ("original", "reference", "outer_joint1", "outer_joint2"):
        curve = sorted([r for r in tau if r["configuration"] == configuration and r["species"] == "ubar"],
                       key=lambda r: finite(r, "xi"))
        if len(curve) != len(POINTS):
            raise ValueError("incomplete tau propagation")
        for left, right in zip(curve, curve[1:]):
            steps.append(dict(configuration=configuration, xi_left=finite(left, "xi"), xi_right=finite(right, "xi"),
                              tau_ubar_left_fm=finite(left, "tau_fm"), tau_ubar_right_fm=finite(right, "tau_fm"),
                              increase_percent=100 * (finite(right, "tau_fm") / finite(left, "tau_fm") - 1)))
        tr = [r for r in tau if r["configuration"] == configuration]
        ob = [r for r in observables if r["configuration"] == configuration]
        if len(tr) != len(POINTS) * 6 or len(ob) != len(POINTS):
            raise ValueError("incomplete final observable coverage")
        maxima[configuration] = dict(tau_percent=100 * max(abs(finite(r, "relative_to_formal")) for r in tr),
            sigma_percent=100 * max(abs(finite(r, "sigma_relative_to_formal")) for r in ob),
            eta_percent=100 * max(abs(finite(r, "eta_relative_to_formal")) for r in ob))
    s_passed = all(finite(r, "base_reference_l1") <= 1e-3 for r in s_rows)
    final_passed = all(value <= 1 for value in maxima["outer_joint2"].values())
    result = dict(schema="transport_mott_window_convergence_review_v1", status="completed_diagnostic_review",
                  computation_host="quote-remote", source_sha=ps["source_sha"], points=list(POINTS), channels=list(CHANNELS),
                  s_jobs=len(s_rows), outer_jobs=len(outer), formal_replay_passed=True,
                  reference_checks_passed=reference_passed, all_base_s_l1_within_production_budget=s_passed,
                  all_outer_converged=outer_passed, all_selected_final_changes_within_one_percent=final_passed,
                  max_s_signed_percent=100 * max(abs(finite(r, "base_reference_signed")) for r in s_rows),
                  max_s_l1_percent=100 * max(finite(r, "base_reference_l1") for r in s_rows),
                  max_outer_base_change_percent=100 * max(abs(finite(r, "joint2_relative_to_base")) for r in outer),
                  max_outer_last_step_percent=100 * max(abs(finite(r, "joint2_relative_to_joint1")) for r in outer),
                  final_change_maxima_percent=maxima, tau_ubar_steps=steps,
                  local_replacement_indicated=not(reference_passed and s_passed and outer_passed and final_passed),
                  original_grid_and_figures_changed=False,
                  limits=["Three frozen backgrounds, two selected mixed-pseudoscalar-sensitive complete channels.",
                          "This is not a pure-eta amplitude ablation or a new all-channel t/B0/background convergence claim.",
                          "No new xi support points or bulk-derivative replay; original densities and terminal integration are retained."],
                  input_hashes={run: digest(args.case / run / "channels.csv") for run in studies},
                  propagation_summary_sha256=digest(path / "summary.json"), review_script_sha256=digest(Path(__file__)))
    args.output.mkdir(parents=True)
    write_csv(args.output / "tau_ubar_step_comparison.csv", steps)
    (args.output / "review.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
