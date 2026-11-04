#!/usr/bin/env python3
"""Summarize observed coarse/fine slope ranges on quote-remote."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("review", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    review = json.loads((args.review/"review.json").read_text(encoding="utf-8-sig"))
    assert review["status"] == "completed_diagnostic_review"
    for name, digest in review["output_sha256"].items():
        assert sha(args.review/name) == digest, name
    fine, coarse = (rows(args.review/(name+".csv")) for name in ("fine_slopes", "coarse_slopes"))
    output = []
    for curve in review["curves"]:
        cs = [r for r in coarse if r["curve_id"] == curve["curve_id"]]
        assert len(cs) >= 2
        lo = min(float(r["xi_left"]) for r in cs)
        hi = max(float(r["xi_right"]) for r in cs)
        fs = [r for r in fine if r["curve_id"] == curve["curve_id"]
              and float(r["xi_left"]) >= lo-1e-12 and float(r["xi_right"]) <= hi+1e-12]
        assert fs
        for key, peak in curve["observables"].items():
            field = key+"_log_slope"
            assert all(math.isfinite(float(r[field])) for r in fs+cs)
            fmin, fmax = min(fs, key=lambda r: float(r[field])), max(fs, key=lambda r: float(r[field]))
            output.append(dict(curve_id=curve["curve_id"], observable=key,
                xi_window_left=lo, xi_window_right=hi,
                fine_min_log_slope=float(fmin[field]), fine_max_log_slope=float(fmax[field]),
                fine_min_xi_left=float(fmin["xi_left"]), fine_min_xi_right=float(fmin["xi_right"]),
                fine_max_xi_left=float(fmax["xi_left"]), fine_max_xi_right=float(fmax["xi_right"]),
                coarse_min_log_slope=min(float(r[field]) for r in cs),
                coarse_max_log_slope=max(float(r[field]) for r in cs),
                all_fine_secants_positive=all(float(r[field]) > 0 for r in fs),
                max_abs_coarse_chord_deviation_percent=peak["max_abs_coarse_chord_deviation_percent"]))
    args.output.mkdir(parents=True, exist_ok=False)
    table = args.output/"slope_ranges.csv"
    with table.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    summary = dict(status="completed_observed_shape_summary",
        review_sha256=sha(args.review/"review.json"), script_sha256=sha(Path(__file__)),
        output_sha256={table.name: sha(table)},
        scope="Only secants wholly inside the shared original-node window; no smoothing or derivative fit",
        limits="Observed extrema are not proven continuous derivative bounds or meson-amplitude contributions",
        primary_observables=[r for r in output if r["observable"] in ("tau_ubar", "sigma_over_T")])
    (args.output/"summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
