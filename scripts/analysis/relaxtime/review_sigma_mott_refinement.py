#!/usr/bin/env python3
"""Review sigma center/energy checks on quote-remote; creates a new evidence view."""
from pathlib import Path
from collections import defaultdict
import argparse
import csv
import hashlib
import json
import math

SPECIES = ("u", "d", "s", "ubar", "dbar", "sbar")
FINAL = ("eta_over_s", "zeta_over_s", "sigma_over_T")
OBSERVABLES = tuple("tau_" + s for s in SPECIES) + FINAL


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_csv(path, rows):
    if rows:
        with path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def slope(a, b, width):
    if min(a, b, width) <= 0:
        raise ValueError("Positive values and interval width required")
    return math.log1p((b-a)/a)/width


def centered_slopes(variant, points):
    rows = []
    if "sig_c000" not in points:
        return rows
    center = points["sig_c000"]
    for h, suffix in ((2e-4, "200"), (1e-4, "100"), (5e-5, "050"), (2e-5, "020")):
        left, right = points.get("sig_m"+suffix), points.get("sig_p"+suffix)
        if left is None or right is None:
            continue
        assert abs(float(center["xi"])-float(left["xi"])-h) < 1e-13
        assert abs(float(right["xi"])-float(center["xi"])-h) < 1e-13
        for key in OBSERVABLES:
            sl = slope(float(left[key]), float(center[key]), h)
            sr = slope(float(center[key]), float(right[key]), h)
            rows.append(dict(variant=variant, h=h, observable=key, left_slope=sl,
                             right_slope=sr, jump=sr-sl, jump_over_h=(sr-sl)/h))
    return rows


def interpolation_error_estimates(points):
    estimates={}
    for ident,data in points.items():
        obs=data["observable"]
        ri={r["channel"]:r for r in data["rates"]}
        relative={}
        absolute_tau={}
        for sp in SPECIES:
            inverse=float(obs["tauinv_"+sp])
            err=sum(float(c["collision_density_weight_inv_fm3"])*float(ri[c["channel"]]["estimated_error"])
                    for c in data["contributions"] if c["species"]==sp)
            assert 0<=err<inverse
            absolute_tau[sp]=err/(inverse*(inverse-err))
            relative["tau_"+sp]=absolute_tau[sp]/float(obs["tau_"+sp])
        for key in FINAL:
            relative[key]=sum(abs(float(w[key+"_coefficient"]))*absolute_tau[w["species"]]
                              for w in data["weights"])/float(obs[key])
        estimates[ident]=relative
    return estimates


def validate_point(row, rates, weights, contributions):
    if len(rates) != 17 or len(weights) != 6 or len(contributions) != 102:
        raise ValueError("Incomplete point decomposition: " + row["point_id"])
    assert len({r["channel"] for r in rates}) == 17
    assert float(row["background_residual"]) <= 1e-10
    assert all(float(row[k]) > 0 and math.isfinite(float(row[k])) for k in OBSERVABLES)
    for rate in rates:
        assert rate["reason"] == "thermal_defect_and_two_shifted_audits" and int(rate["guard_phase"]) == 2
        assert float(rate["rate"]) > 0 and 0 <= float(rate["estimated_error"])
    for species in SPECIES:
        value = sum(float(r["inverse_tau_contribution_inv_fm"]) for r in contributions if r["species"] == species)
        assert math.isclose(value, float(row["tauinv_"+species]), rel_tol=1e-10)
        assert math.isclose(value*float(row["tau_"+species]), 1.0, rel_tol=1e-10)
    for key in FINAL:
        value = sum(float(w[key+"_coefficient"])*float(row["tau_"+w["species"]]) for w in weights)
        assert math.isclose(value, float(row[key]), rel_tol=1e-10)


def load_run(folder, plan, preview=False):
    if not (folder/"observables.csv").exists():
        return {}
    summary_path = folder/"summary.json"
    if summary_path.exists():
        summary = read_json(summary_path)
        assert summary["status"] == "completed_diagnostic"
        for name, digest in summary["files"].items():
            assert sha(folder/name) == digest, (folder.name, name)
    elif not preview:
        return {}
    obs = read_csv(folder/"observables.csv")
    groups = {}
    for filename in ("rates", "species_weights", "channel_contributions"):
        by_point = defaultdict(list)
        for row in read_csv(folder/(filename+".csv")):
            by_point[row["point_id"]].append(row)
        groups[filename] = by_point
    output = {}
    for row in obs:
        match = [p for p in plan if abs(float(p["xi"])-float(row["xi"])) < 1e-12]
        if not match:
            continue
        assert len(match) == 1
        old_id = row["point_id"]
        data = dict(observable=dict(row, point_id=match[0]["point_id"]),
            rates=groups["rates"][old_id], weights=groups["species_weights"][old_id],
            contributions=groups["channel_contributions"][old_id], source=str(folder), original_point_id=old_id)
        validate_point(row, data["rates"], data["weights"], data["contributions"])
        output[match[0]["point_id"]] = data
    return output


def fit_log_power(xs, ys):
    xx, yy = [math.log(x) for x in xs], [math.log(y) for y in ys]
    mx, my = sum(xx)/len(xx), sum(yy)/len(yy)
    return sum((x-mx)*(y-my) for x, y in zip(xx, yy))/sum((x-mx)**2 for x in xx)


def logmean(a, b):
    return a if a == b else (b-a)/math.log1p((b-a)/a)


def slope_decomposition(variant, points):
    channel_rows, terminal_rows = [], []
    if "sig_c000" not in points:
        return channel_rows, terminal_rows
    for h, suffix in ((2e-4,"200"),(1e-4,"100"),(5e-5,"050"),(2e-5,"020")):
        if not {"sig_m"+suffix,"sig_p"+suffix} <= set(points):
            continue
        triplet = [points[key] for key in ("sig_m"+suffix,"sig_c000","sig_p"+suffix)]
        parts = []
        terminals = []
        for left, right in zip(triplet[:-1],triplet[1:]):
            a, b = left["observable"], right["observable"]
            ai = {(r["species"],r["channel"]):r for r in left["contributions"]}
            bi = {(r["species"],r["channel"]):r for r in right["contributions"]}
            aw = {r["species"]:r for r in left["weights"]}
            bw = {r["species"]:r for r in right["weights"]}
            by_channel, by_species = {}, {}
            for sp, channel in ai:
                ac, bc = ai[sp,channel], bi[sp,channel]
                n0,n1 = float(ac["collision_density_weight_inv_fm3"]),float(bc["collision_density_weight_inv_fm3"])
                w0,w1 = float(ac["rate_fm2"]),float(bc["rate_fm2"])
                denominator = h*logmean(float(a["tauinv_"+sp]),float(b["tauinv_"+sp]))
                rate_part = -(n0+n1)/2*(w1-w0)/denominator
                density_part = -(w0+w1)/2*(n1-n0)/denominator
                by_channel[sp,channel] = rate_part,density_part
            for sp in SPECIES:
                total = sum(sum(v) for (s,c),v in by_channel.items() if s == sp)
                assert abs(total-slope(float(a["tau_"+sp]),float(b["tau_"+sp]),h)) < 1e-7
            for key in FINAL:
                denominator = h*logmean(float(a[key]),float(b[key]))
                for sp in SPECIES:
                    c0,c1=float(aw[sp][key+"_coefficient"]),float(bw[sp][key+"_coefficient"])
                    t0,t1=float(a["tau_"+sp]),float(b["tau_"+sp])
                    by_species[key,sp]=((c0+c1)/2*(t1-t0)/denominator,
                                       (t0+t1)/2*(c1-c0)/denominator)
                total=sum(sum(v) for (k,sp),v in by_species.items() if k==key)
                assert abs(total-slope(float(a[key]),float(b[key]),h)) < 1e-7
            parts.append(by_channel);terminals.append(by_species)
        for key in parts[0]:
            l,r=parts[0][key],parts[1][key]
            channel_rows.append(dict(variant=variant,h=h,species=key[0],channel=key[1],
                jump_from_rate=r[0]-l[0],jump_from_density=r[1]-l[1],
                jump_contribution=sum(r)-sum(l)))
        for key in terminals[0]:
            l,r=terminals[0][key],terminals[1][key]
            terminal_rows.append(dict(variant=variant,h=h,observable=key[0],species=key[1],
                jump_from_tau=r[0]-l[0],jump_from_coefficient=r[1]-l[1],
                jump_contribution=sum(r)-sum(l)))
    return channel_rows,terminal_rows


def spectral_summary(case):
    folder = case/"spectrum_v2"
    if not (folder/"summary.json").exists():
        return {}, []
    summary = read_json(folder/"summary.json")
    for name, digest in summary["files"].items():
        assert sha(folder/name) == digest
    roots = [r for r in read_csv(folder/"real_axis_roots.csv") if r["backend"] == "native"]
    area = {(r["point_id"], r["channel"]): r for r in read_csv(folder/"peak_areas.csv") if int(r["order"]) == 32}
    spectra = defaultdict(list)
    for row in read_csv(folder/"cross_section_spectrum.csv"):
        spectra[row["point_id"], row["channel"]].append(row)
    result = []
    for root in roots:
        for channel in ("uubar_to_uubar", "uubar_to_ddbar"):
            key = root["point_id"], channel
            high = max(spectra[key], key=lambda r: float(r["s_mixed_scalar_square_fm2"]))
            result.append(dict(point_id=root["point_id"], channel=channel,
                delta_xi=float(root["delta_xi"]), peak_x_inv_fm2=float(root["x_inv_fm2"]),
                half_width_inv_fm2=float(root["half_width_inv_fm2"]),
                mixed_peak_height_fm2=float(high["s_mixed_scalar_square_fm2"]),
                full_peak_height_fm2=float(high["native_sigma_fm2"]),
                mixed_area_dimensionless=float(area[key]["mixed_square_area_dimensionless"])))
    selected = [r for r in result if r["channel"] == "uubar_to_uubar"]
    x = [r["peak_x_inv_fm2"] for r in selected]
    return dict(width_power_vs_peak_x=fit_log_power(x, [r["half_width_inv_fm2"] for r in selected]),
        height_power_vs_peak_x=fit_log_power(x, [r["mixed_peak_height_fm2"] for r in selected]),
        area_power_vs_peak_x=fit_log_power(x, [r["mixed_area_dimensionless"] for r in selected]),
        max_amplitude_reconstruction_relative=summary["max_reconstruction_relative"]), result


def hybrid_points(points, bands, cutoff, weight_variant="fine"):
    result, corrections = {}, []
    for ident, data in points.items():
        relevant = [r for r in bands if r["point_id"] == ident and float(r["x_right_inv_fm2"]) <= cutoff*(1+1e-12)]
        if not relevant:
            continue
        delta = defaultdict(float)
        field = "continuous_full_rate_fm2" if weight_variant == "fine" else "base_weight_full_rate_fm2"
        for row in relevant:
            delta[row["channel"]] += float(row[field])-float(row["discrete_full_rate_fm2"])
        assert set(delta) == {"uubar_to_uubar", "uubar_to_ddbar"}
        old = data["observable"]
        corrected = dict(old)
        for species in SPECIES:
            inverse_change = sum(float(c["collision_density_weight_inv_fm3"])*delta[c["channel"]]
                for c in data["contributions"] if c["species"] == species and c["channel"] in delta)
            corrected["tauinv_"+species] = float(old["tauinv_"+species])+inverse_change
            corrected["tau_"+species] = 1/corrected["tauinv_"+species]
            assert corrected["tau_"+species] > 0
        for key in FINAL:
            corrected[key] = sum(float(w[key+"_coefficient"])*corrected["tau_"+w["species"]] for w in data["weights"])
        result[ident] = corrected
        for channel, change in delta.items():
            original = float(next(r["rate"] for r in data["rates"] if r["channel"] == channel))
            assert original+change > 0
            corrections.append(dict(point_id=ident, cutoff_inv_fm2=cutoff, weight_variant=weight_variant,
                channel=channel, raw_rate_fm2=original, delta_rate_fm2=change,
                relative_change=change/original, diagnostic_rate_fm2=original+change))
    return result, corrections


def thermal_controls(precision_grids, bands, cutoffs=(1e-3, 2e-3, 5e-3, 1e-2)):
    """Keep all three secant endpoints at the same s interpolation tolerance."""
    slopes, corrections = [], []
    for precision, points in precision_grids.items():
        for cutoff in cutoffs:
            for weight_variant in ("fine", "base"):
                hybrid, changes = hybrid_points(points, bands, cutoff, weight_variant)
                corrections.extend(dict(precision=precision, **row) for row in changes)
                slopes.extend(centered_slopes(f"hybrid_{precision}_{weight_variant}_x{cutoff:g}", hybrid))
    return slopes, corrections


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--preview", action="store_true")
    args = parser.parse_args()
    case, out = args.case.resolve(), args.output.resolve()
    out.mkdir(exist_ok=False)
    plan = read_csv(case/"inputs"/"transport_plan.csv")
    s6 = load_run(case/"inputs"/"previous", plan)
    for name in ("left_s6_v1", "right_s6_v1"):
        extra = load_run(case/name, plan, args.preview)
        assert not set(extra).intersection(s6)
        s6.update(extra)
    s7 = {}
    for folder in sorted(case.glob("*s7*")):
        if folder.is_dir():
            extra = load_run(folder, plan, args.preview)
            assert not set(extra).intersection(s7)
            s7.update(extra)
    if not args.preview:
        assert len(s6) == 9, "Full nine-point center grid required"
        assert {"sig_m020", "sig_c000", "sig_p020"} <= set(s7), "Inner equal-endpoint precision controls required"
    slopes = []
    channel_decomposition,terminal_decomposition=[],[]
    error_rows=[]
    for variant, points in (("s6", s6), ("s7", s7)):
        slopes.extend(centered_slopes(variant, {k: v["observable"] for k, v in points.items()}))
        c,t=slope_decomposition(variant,points)
        channel_decomposition.extend(c);terminal_decomposition.extend(t)
        estimates=interpolation_error_estimates(points)
        if "sig_c000" in estimates:
            for h,suffix in ((2e-4,"200"),(1e-4,"100"),(5e-5,"050"),(2e-5,"020")):
                ids=("sig_m"+suffix,"sig_c000","sig_p"+suffix)
                if not set(ids)<=set(estimates):
                    continue
                for key in OBSERVABLES:
                    errors=[-math.log1p(-estimates[i][key]) for i in ids]
                    error_rows.append(dict(variant=variant,h=h,observable=key,
                        estimated_s_interpolation_jump_error=(errors[0]+2*errors[1]+errors[2])/h,
                        scope="s interpolation only; estimator, not total certified error"))
    precision = []
    low = {(r["h"], r["observable"]): r for r in slopes if r["variant"] == "s6"}
    for row in slopes:
        if row["variant"] == "s7" and (row["h"], row["observable"]) in low:
            previous = low[row["h"], row["observable"]]
            precision.append(dict(h=row["h"], observable=row["observable"], s6_jump=previous["jump"],
                s7_jump=row["jump"], jump_difference=row["jump"]-previous["jump"],
                max_side_slope_difference=max(abs(row["left_slope"]-previous["left_slope"]),
                                            abs(row["right_slope"]-previous["right_slope"]))))
    powers, peaks = spectral_summary(case)
    corrections = []
    audit = case/"thermal_audit_v1"
    thermal_summary = {}
    if (audit/"summary.json").exists():
        report = read_json(audit/"summary.json")
        for name, digest in report["files"].items():
            assert sha(audit/name) == digest
        bands = read_csv(audit/"thermal_band_comparison.csv")
        thermal_slopes, corrections = thermal_controls({"s6": s6, "s7": s7}, bands)
        slopes.extend(thermal_slopes)
        peak_rates = [r for r in read_csv(audit/"thermal_peak_areas.csv") if r["channel"] == "uubar_to_uubar"]
        thermal_summary = dict(continuum_peak_power_vs_x=fit_log_power(
            [float(r["x_peak_inv_fm2"]) for r in peak_rates],
            [float(r["continuous_mixed_square_rate_fm2"]) for r in peak_rates]),
            peak_direct_event_checks=sum(int(r["events"]) for r in read_csv(audit/"thermal_peak_areas.csv")),
            maximum_weight_quadrature_relative=max(abs(float(r["quadrature_relative"])) for r in read_csv(audit/"continuous_weight.csv")))
    elif not args.preview:
        raise AssertionError("Continuous thermal audit is required")
    write_csv(out/"centered_slopes.csv", slopes)
    write_csv(out/"precision_comparison.csv", precision)
    write_csv(out/"peak_scaling.csv", peaks)
    write_csv(out/"near_band_rate_corrections.csv", corrections)
    write_csv(out/"centered_channel_decomposition.csv",channel_decomposition)
    write_csv(out/"centered_terminal_decomposition.csv",terminal_decomposition)
    write_csv(out/"interpolation_error_estimates.csv",error_rows)
    points = [dict(variant=tag, point_id=ident, source=data["source"], original_point_id=data["original_point_id"],
        xi=float(data["observable"]["xi"]), **{k: float(data["observable"][k]) for k in OBSERVABLES})
        for tag, values in (("s6", s6), ("s7", s7)) for ident, data in sorted(values.items())]
    write_csv(out/"observables.csv", points)
    summary = dict(status="partial_diagnostic" if args.preview else "completed_diagnostic_review",
        source_s6_points=len(s6), source_s7_points=len(s7), spectral=powers, thermal=thermal_summary,
        precision_policy="Every raw or hybrid slope uses three endpoints from one precision tier; missing s7 endpoints are not filled with s6.",
        precision_tiers=dict(s6=dict(s_rtol=1e-6,points=len(s6)),s7=dict(s_rtol=1e-7,points=len(s7))),
        physical_kernel_changed=False, manuscript_eligible=False,
        interpretation="Finite-window slopes and a local integration-variable control; not a proof of differentiability. Hybrid observables replace only the indicated near-threshold bands of two light annihilation channels and are diagnostic.",
        script_sha256=sha(Path(__file__)), files={p.name: sha(p) for p in out.iterdir() if p.is_file()})
    (out/"summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
