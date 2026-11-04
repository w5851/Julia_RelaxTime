#!/usr/bin/env python3
"""Review poles and threshold screens on quote-remote; no exclusivity assumption."""
from collections import Counter, defaultdict
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import re


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    if rows:
        with path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def finite(value):
    value = float(value)
    return value if math.isfinite(value) else None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    probe = args.case / "survey_v1"
    inputs = args.case / "inputs/source_rows.csv"
    raw_summary = (probe / "summary.json").read_text(encoding="utf-8-sig")
    # JSON3.pretty in the frozen environment did not escape Windows Dict keys.
    # Repair only the known candidate-file key prefix in a separate review copy.
    normalized_summary, repaired_keys = re.subn(r'"candidates\\([^"\n]+)"\s*:',
                                               r'"candidates/\1":', raw_summary)
    assert repaired_keys in (0, 320), repaired_keys
    report = json.loads(normalized_summary)
    (args.output / "survey_summary_readable.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")
    assert sha(inputs) == report["input_rows_sha256"]
    for name, digest in report["output_sha256"].items():
        assert sha(probe / name) == digest, name
    source = read_csv(inputs)
    by_id = {r["point_id"]: r for r in source}
    groups = defaultdict(list)
    for row in source:
        groups[row["window"]].append(row)
    pairs = defaultdict(dict)
    for row in read_csv(probe / "poles.csv"):
        pairs[row["meson"], row["point_id"]][row["direction"]] = row
    selected, checks = {}, []
    for (meson, point), pair in pairs.items():
        assert set(pair) == {"forward", "reverse"}
        a, b = pair["forward"], pair["reverse"]
        qualified = a["qualified"] == b["qualified"] == "true"
        dm = finite(abs(float(a["mass_MeV"]) - float(b["mass_MeV"])))
        dg = finite(abs(float(a["gamma_MeV"]) - float(b["gamma_MeV"])))
        roots = max(int(a["unique_qualified_roots"]), int(b["unique_qualified_roots"]))
        consistent = qualified and dm is not None and dg is not None and max(dm, dg) < 1e-3
        trusted = consistent and roots == 1
        checks.append(dict(meson=meson, point_id=point, window=a["window"], xi=float(a["xi"]),
                           both_qualified=qualified, bidirectional_consistent=consistent,
                           mass_difference_MeV=dm, gamma_difference_MeV=dg,
                           max_distinct_roots=roots, trusted_unique=trusted))
        selected[meson, point] = dict(a, trusted_unique=trusted, bidirectional_consistent=consistent)
    write_csv(args.output / "root_checks.csv", checks)
    events, ranges = [], []
    mixed = {"eta", "eta_prime", "sigma", "sigma_prime"}
    for window, curve in sorted(groups.items()):
        curve.sort(key=lambda r: float(r["xi"]))
        for meson in report["mesons"]:
            cuts = (("uu", "gap_light_MeV"), ("ss", "gap_strange_MeV")) if meson in mixed else (("pair", "gap_MeV"),)
            for cut, gap_key in cuts:
                roots = [selected[meson, r["point_id"]] for r in curve]
                valid = [r for r in roots if r["trusted_unique"]]
                gaps = [float(r[gap_key]) for r in valid]
                ranges.append(dict(window=window, meson=meson, cut=cut,
                    points=len(roots), trusted_points=len(valid),
                    min_gap_MeV=min(gaps, default=None), max_gap_MeV=max(gaps, default=None),
                    unresolved_points=len(roots)-len(valid)))
                for left, right in zip(curve, curve[1:]):
                    a, b = selected[meson, left["point_id"]], selected[meson, right["point_id"]]
                    ga, gb = float(a[gap_key]), float(b[gap_key])
                    if not (math.isfinite(ga) and math.isfinite(gb) and ga * gb < 0):
                        continue
                    continuous = left["continuation_segment"] == right["continuation_segment"]
                    trusted = a["trusted_unique"] and b["trusted_unique"]
                    bound = a if ga < 0 else b
                    width_bound = float(bound["gamma_MeV"])
                    if not continuous:
                        status = "first_order_background_jump"
                    elif not trusted:
                        status = "unresolved_spectral_branch"
                    elif cut == "ss":
                        status = "higher_cut_crossing"
                    elif abs(width_bound) < 1e-3:
                        status = "coarse_confirmed_mott"
                    else:
                        status = "mass_crossing_with_finite_width"
                    events.append(dict(window=window, meson=meson, cut=cut, status=status,
                        T_MeV=float(left["T_MeV"]), muB_MeV=float(left["muB_MeV"]),
                        point_left=left["point_id"], point_right=right["point_id"],
                        xi_left=float(left["xi"]), xi_right=float(right["xi"]),
                        gap_left_MeV=ga, gap_right_MeV=gb,
                        gamma_left_MeV=float(a["gamma_MeV"]), gamma_right_MeV=float(b["gamma_MeV"]),
                        continuous_background=continuous, both_trusted_unique=trusted,
                        direction="reverse_mott" if ga > 0 else "mott"))
    write_csv(args.output / "root_crossings.csv", events)
    write_csv(args.output / "root_ranges.csv", ranges)
    threshold_groups = defaultdict(list)
    for row in read_csv(probe / "thresholds.csv"):
        threshold_groups[row["window"], row["meson"], row["cut"], row["evaluator"]].append(row)
    inverse_events, inverse_ranges = [], []
    for (window, meson, cut, evaluator), rows in sorted(threshold_groups.items()):
        rows.sort(key=lambda r: float(r["xi"]))
        valid = [r for r in rows if finite(r["inverse_abs"]) is not None]
        closest = min(valid, key=lambda r: float(r["inverse_abs"]), default=None)
        inverse_ranges.append(dict(window=window, meson=meson, cut=cut, evaluator=evaluator,
            points=len(rows), finite_points=len(valid),
            min_inverse_re=min((float(r["inverse_re"]) for r in valid), default=None),
            max_inverse_re=max((float(r["inverse_re"]) for r in valid), default=None),
            min_inverse_abs=None if closest is None else float(closest["inverse_abs"]),
            xi_closest=None if closest is None else float(closest["xi"])))
        for a, b in zip(rows, rows[1:]):
            fa, fb = float(a["inverse_re"]), float(b["inverse_re"])
            if not (math.isfinite(fa) and math.isfinite(fb) and fa * fb < 0):
                continue
            continuous = a["continuation_segment"] == b["continuation_segment"]
            closed = max(abs(float(a["inverse_im"])), abs(float(b["inverse_im"]))) < 1e-8
            left = by_id[a["point_id"]]
            inverse_events.append(dict(window=window, meson=meson, cut=cut, evaluator=evaluator,
                T_MeV=float(left["T_MeV"]), muB_MeV=float(left["muB_MeV"]),
                point_left=a["point_id"], point_right=b["point_id"],
                xi_left=float(a["xi"]), xi_right=float(b["xi"]),
                inverse_left=fa, inverse_right=fb,
                imaginary_left=float(a["inverse_im"]), imaginary_right=float(b["inverse_im"]),
                continuous_background=continuous, negligible_imaginary=closed,
                status="candidate_only" if continuous else "first_order_background_jump"))
    write_csv(args.output / "threshold_crossings.csv", inverse_events)
    write_csv(args.output / "threshold_ranges.csv", inverse_ranges)
    summary = dict(schema="transport_other_mott_survey_review_v1", status="completed_diagnostic",
        input_rows_sha256=sha(inputs), survey_summary_sha256=sha(probe / "summary.json"),
        reviewer_sha256=sha(Path(__file__)), points=len(source), curves=len(groups),
        root_status_counts={m: dict(Counter("trusted" if r["trusted_unique"] else "unresolved"
            for r in checks if r["meson"] == m)) for m in report["mesons"]},
        root_crossing_counts=dict(Counter(r["status"] for r in events)), root_crossings=events,
        diagnostic_threshold_crossings=[r for r in inverse_events if r["evaluator"] == "diagnostic_endpoint"],
        scope="Finite-grid survey, not proof excluding arbitrarily narrow crossings; high mixed cut is separate from bound-state Mott.",
        transport_recomputed=False, production_kernel_changed=False,
        metadata_repair=dict(original_retained=True, candidate_path_keys_normalized=repaired_keys,
                             numerical_values_changed=False))
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, allow_nan=False))


if __name__ == "__main__":
    main()
