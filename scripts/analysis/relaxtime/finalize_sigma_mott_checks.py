#!/usr/bin/env python3
"""Combine uniform-precision continuum-peak and cutoff-membership diagnostics."""
from pathlib import Path
import argparse
import importlib.util
import json
import math
from copy import deepcopy

WIDTHS = ((2e-4, "200"), (1e-4, "100"), (5e-5, "050"), (2e-5, "020"))
CUTOFF_CHANNELS = {"uu_to_uu", "ud_to_ud"}


def validate_cutoff_events(effects, events):
    """Reconstruct every side correction from signed, directly evaluated events."""
    expected = {(h, channel) for h, _ in WIDTHS for channel in CUTOFF_CHANNELS}
    assert len(effects) == len(expected)
    assert {(float(r["h"]), r["channel"]) for r in effects} == expected
    checked, maximum_error = 0, 0.0
    for row in effects:
        h, channel = float(row["h"]), row["channel"]
        suffix = dict(WIDTHS)[h]
        for side, tag in (("left", "m"), ("right", "p")):
            selected = [r for r in events if float(r["h"]) == h and r["channel"] == channel
                        and r["point_id"] == "sig_"+tag+suffix]
            assert len(selected) == int(row[side+"_changed_events"])
            for event in selected:
                accepted = str(event["accepted_here"]).lower() == "true"
                center = str(event["accepted_center"]).lower() == "true"
                assert accepted != center
                weight, cross_section = float(event["thermal_weight"]), float(event["direct_sigma_fm2"])
                assert weight > 0 and cross_section >= 0
                value = weight*cross_section*(int(center)-int(accepted))
                assert math.isclose(value, float(event["fixed_membership_delta_rate_fm2"]),
                                    rel_tol=1e-12, abs_tol=1e-30)
            summed = math.fsum(float(r["fixed_membership_delta_rate_fm2"]) for r in selected)
            expected_sum = float(row[side+"_delta_rate_fm2"])
            error = abs(summed-expected_sum)
            maximum_error = max(maximum_error, error)
            assert math.isclose(summed, expected_sum, rel_tol=1e-11, abs_tol=1e-25)
            checked += len(selected)
    assert checked == len(events), "Unassigned cutoff event"
    return dict(status="passed",cases=len(effects),events=checked,
                max_sum_absolute_error_fm2=maximum_error)


def apply_fixed_membership(review, hybrid, points, ids, deltas):
    """Apply rate changes before inverting tau, then rebuild terminal quantities."""
    corrected = deepcopy(hybrid)
    for index, side in ((0, "left"), (2, "right")):
        ident = ids[index]
        data, delta = points[ident], deltas[side]
        assert set(delta) == CUTOFF_CHANNELS
        for channel, change in delta.items():
            raw = float(next(r["rate"] for r in data["rates"] if r["channel"] == channel))
            assert raw+change > 0
        for species in review.SPECIES:
            inverse_change = math.fsum(float(c["collision_density_weight_inv_fm3"])*delta[c["channel"]]
                for c in data["contributions"] if c["species"] == species and c["channel"] in delta)
            inverse = float(hybrid[ident]["tauinv_"+species])+inverse_change
            assert math.isfinite(inverse) and inverse > 0
            corrected[ident]["tauinv_"+species] = inverse
            corrected[ident]["tau_"+species] = 1/inverse
        for key in review.FINAL:
            corrected[ident][key] = math.fsum(float(w[key+"_coefficient"])*corrected[ident]["tau_"+w["species"]]
                                             for w in data["weights"])
    return corrected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--base-review", default="review_uniform_v1")
    parser.add_argument("--review-script", default="review_sigma_mott_refinement_v4.py")
    args = parser.parse_args()
    case, out = args.case.resolve(), args.output.resolve()
    spec = importlib.util.spec_from_file_location("sigma_review", case/args.review_script)
    review = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(review)
    base, cutoff = case/args.base_review, case/"cutoff_audit_v1"
    summary = review.read_json(base/"summary.json")
    cutoff_summary = review.read_json(cutoff/"summary.json")
    assert summary["status"] == "completed_diagnostic_review"
    assert cutoff_summary["status"] == "completed_diagnostic"
    assert "precision_policy" in summary
    for folder, report in ((base, summary), (cutoff, cutoff_summary),
                           (case/"thermal_audit_v1", review.read_json(case/"thermal_audit_v1/summary.json"))):
        for name, digest in report["files"].items():
            assert review.sha(folder/name) == digest
    plan = review.read_csv(case/"inputs/transport_plan.csv")
    s6 = review.load_run(case/"inputs/previous", plan)
    for name in ("left_s6_v1", "right_s6_v1"):
        extra = review.load_run(case/name, plan)
        assert not set(extra).intersection(s6)
        s6.update(extra)
    s7 = {}
    for name in ("reference_center_s7_v1", "reference_left_s7_v1", "reference_right_s7_v1"):
        extra = review.load_run(case/name, plan)
        assert not set(extra).intersection(s7)
        s7.update(extra)
    assert len(s6) == 9 and set(s7) == {"sig_m020", "sig_c000", "sig_p020"}
    # Cutoff events depend on the frozen background, not the s-table tolerance.
    for ident, data in s7.items():
        fields = ("xi", "m_u", "m_d", "m_s", "Phi", "Phibar")+tuple("n_"+s for s in review.SPECIES)
        for key in fields:
            assert float(data["observable"][key]) == float(s6[ident]["observable"][key])
        old_weights = {w["species"]:w for w in s6[ident]["weights"]}
        for new in data["weights"]:
            old = old_weights[new["species"]]
            for key in review.FINAL:
                assert float(old[key+"_coefficient"]) == float(new[key+"_coefficient"])
    bands = review.read_csv(case/"thermal_audit_v1/thermal_band_comparison.csv")
    effects = review.read_csv(cutoff/"cutoff_effect.csv")
    audit = validate_cutoff_events(effects, review.read_csv(cutoff/"cutoff_events.csv"))
    new_slopes, corrections = [], []
    for precision, points in (("s6", s6), ("s7", s7)):
        for h, suffix in WIDTHS:
            ids = ("sig_m"+suffix, "sig_c000", "sig_p"+suffix)
            if not set(ids) <= set(points):
                continue
            deltas = {side: {r["channel"]: float(r[side+"_delta_rate_fm2"])
                            for r in effects if float(r["h"]) == h} for side in ("left", "right")}
            for xmax in (0.001, 0.01):
                hybrid, _ = review.hybrid_points({i: points[i] for i in ids}, bands, xmax)
                corrected = apply_fixed_membership(review, hybrid, points, ids, deltas)
                label = f"combined_{precision}_x{xmax:g}"
                new_slopes.extend(r for r in review.centered_slopes(label, corrected) if r["h"] == h)
                for index, side in ((0, "left"), (2, "right")):
                    for channel, change in deltas[side].items():
                        corrections.append(dict(precision=precision,h=h,point_id=ids[index],
                            near_cutoff_inv_fm2=xmax,channel=channel,fixed_membership_delta_rate_fm2=change))
    out.mkdir(exist_ok=False)
    for path in base.iterdir():
        if path.is_file() and path.name!="summary.json":(out/path.name).write_bytes(path.read_bytes())
    original=review.read_csv(out/"centered_slopes.csv")
    review.write_csv(out/"centered_slopes.csv",original+new_slopes)
    review.write_csv(out/"cutoff_membership_rate_corrections.csv",corrections)
    for name in ("cutoff_effect.csv","cutoff_events.csv"):
        (out/name).write_bytes((cutoff/name).read_bytes())
    (out/"cutoff_event_validation.json").write_text(json.dumps(audit,indent=2)+"\n",encoding="utf-8")
    all_slopes = review.read_csv(out/"centered_slopes.csv")
    summary.update(finalization="Uniform s6 curves and an independent uniform inner s7 triplet; cutoff control added",
        cutoff_audit_sha256=review.sha(cutoff/"summary.json"),base_review=args.base_review,
        original_review_sha256=review.sha(base/"summary.json"),finalizer_sha256=review.sha(Path(__file__)),
        cutoff_event_validation=audit,
        cutoff_interpretation="The physical cutoff is retained. Holding discrete event membership at the center attributes quadrature staircase terms; this is not a new continuum solution.",
        comparison_scope="Near-threshold two-channel replacement and fixed-cutoff counterfactual are applied independently within each precision tier. All other rates, densities and terminal coefficients retain their frozen values.",
        selected_inner_slopes=[r for r in all_slopes if float(r["h"]) == 2e-5
            and r["observable"] in ("tau_ubar", "sigma_over_T")
            and r["variant"] in ("s6", "s7", "hybrid_s6_fine_x0.001", "hybrid_s7_fine_x0.001", "combined_s6_x0.001", "combined_s7_x0.001")],
        files={p.name:review.sha(p) for p in out.iterdir() if p.is_file()})
    (out/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(dict(status=summary["status"],cutoff_event_validation=audit,
        source_s6_points=len(s6),source_s7_points=len(s7),files=len(summary["files"]))),flush=True)


if __name__=="__main__":main()
