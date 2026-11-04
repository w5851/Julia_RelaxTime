#!/usr/bin/env python3
"""Local log-secant and exact finite-difference attribution, on quote-remote."""
from collections import defaultdict
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import tomllib

SPECIES = ("u", "d", "s", "ubar", "dbar", "sbar")
FINAL = ("eta_over_s", "zeta_over_s", "sigma_over_T")
OBSERVABLES = tuple("tau_" + s for s in SPECIES) + FINAL
K_S_CHANNELS = {"usbar_to_usbar", "subar_to_subar"}


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


def logmean(a, b):
    if a <= 0 or b <= 0:
        raise ValueError("logarithmic mean requires positive endpoints")
    return a if a == b else (b-a) / math.log1p((b-a)/a)


def log_slope(a, b, dx):
    if dx <= 0 or a <= 0 or b <= 0:
        raise ValueError("positive values and increasing xi required")
    return math.log1p((b-a)/a) / dx


def product_difference(a0, a1, b0, b1):
    """Exact midpoint identity: delta(ab)=mean(a)delta(b)+mean(b)delta(a)."""
    return (a0+a1)/2 * (b1-b0), (b0+b1)/2 * (a1-a0)


def compare_native_slopes(native_slopes, references):
    """Compare identical endpoints even when the reference adds interior nodes."""
    checks = []
    for row in native_slopes:
        left, right, observable = row["point_left"], row["point_right"], row["observable"]
        if left not in references or right not in references:
            continue
        a, b = references[left], references[right]
        if any(abs(float(r["xi"])-float(row[k])) >= 1e-12
               for r, k in ((a,"xi_left"),(b,"xi_right"))):
            raise ValueError("Slope comparison requires identical endpoint coordinates")
        fine = row["log_slope"]
        tight = log_slope(float(a[observable]), float(b[observable]), row["interval_width"])
        checks.append(dict(curve_id=row["curve_id"],point_left=left,point_right=right,observable=observable,
            fine_log_slope=fine,tighter_log_slope=tight,absolute_difference=abs(tight-fine),
            relative_difference=None if abs(tight)<1e-8 else abs((tight-fine)/tight)))
    return checks


def collect(case, names, preview=False):
    obs, rates, weights, channels, provenance = [], [], [], [], {}
    for name in names:
        folder = case / name
        summary_path = folder / "summary.json"
        if not summary_path.exists() and preview:
            if not (folder / "observables.csv").exists():
                continue
        else:
            summary = json.loads(summary_path.read_text(encoding="utf-8-sig"))
            assert summary["status"] == "completed_diagnostic"
            for file, digest in summary["files"].items():
                assert sha(folder / file) == digest, (name, file)
            provenance[name] = sha(summary_path)
        rows = read_csv(folder / "observables.csv")
        ids = {r["point_id"] for r in rows}
        obs.extend(rows)
        local_rates = [r for r in read_csv(folder / "rates.csv") if r["point_id"] in ids]
        variants = {r["variant"] for r in rows}
        assert len(variants) == 1
        settings = tomllib.loads((case/"inputs"/(next(iter(variants))+".toml")).read_text(encoding="utf-8"))["sampling"]
        assert all(r["reason"] == "thermal_defect_and_two_shifted_audits" and int(r["guard_phase"]) == 2 for r in local_rates)
        assert all(0 <= float(r["estimated_error"]) <= (settings["atol"]+settings["rtol"]*abs(float(r["rate"]))) * (1+1e-10)
                   and float(r["rate"]) > 0 for r in local_rates)
        for target, file in ((rates, "rates.csv"), (weights, "species_weights.csv"),
                             (channels, "channel_contributions.csv")):
            target.extend(r for r in read_csv(folder / file) if r["point_id"] in ids)
    assert len(obs) == len({r["point_id"] for r in obs}), "duplicate point results"
    assert len(rates) == 17*len(obs) and len(weights) == 6*len(obs) and len(channels) == 102*len(obs)
    assert all(float(r[key]) > 0 and math.isfinite(float(r[key])) for r in obs for key in OBSERVABLES)
    return obs, rates, weights, channels, provenance


def analyze(obs, rates, weights, channels):
    ci = {(r["point_id"], r["species"], r["channel"]): r for r in channels}
    wi = {(r["point_id"], r["species"]): r for r in weights}
    ri = {(r["point_id"], r["channel"]): r for r in rates}
    processes = sorted({r["channel"] for r in rates})
    assert len(processes) == 17
    estimates = {}
    for row in obs:
        ident = row["point_id"]
        tau_errors = {}
        for sp in SPECIES:
            inverse = float(row["tauinv_"+sp])
            reconstructed = sum(float(ci[ident, sp, c]["inverse_tau_contribution_inv_fm"]) for c in processes)
            assert math.isclose(reconstructed, inverse, rel_tol=1e-10)
            inverse_error = sum(float(ci[ident, sp, c]["collision_density_weight_inv_fm3"])
                                * float(ri[ident, c]["estimated_error"]) for c in processes)
            assert 0 <= inverse_error < inverse
            tau_errors[sp] = inverse_error / (inverse*(inverse-inverse_error))
            estimates[ident, "tau_"+sp] = tau_errors[sp] / float(row["tau_"+sp])
        for key in FINAL:
            coefficients = {sp: float(wi[ident, sp][key+"_coefficient"]) for sp in SPECIES}
            rebuilt = sum(coefficients[sp]*float(row["tau_"+sp]) for sp in SPECIES)
            assert math.isclose(rebuilt, float(row[key]), rel_tol=1e-10)
            estimates[ident, key] = sum(abs(coefficients[sp])*tau_errors[sp] for sp in SPECIES) / float(row[key])
    groups = defaultdict(list)
    for row in obs:
        groups[row["curve_id"]].append(row)
    slopes, channel_derivatives, final_derivatives, group_derivatives, contrasts, centered = [], [], [], [], [], []
    closure_error = 0.0
    for curve, rows in groups.items():
        rows.sort(key=lambda r: float(r["xi"]))
        curve_slopes = []
        for a, b in zip(rows, rows[1:]):
            ida, idb = a["point_id"], b["point_id"]
            xa, xb = float(a["xi"]), float(b["xi"])
            dx = xb-xa
            info = dict(curve_id=curve, meson=a["meson"], window=a["window"],
                T_MeV=float(a["T_MeV"]), muB_MeV=float(a["muB_MeV"]), point_left=ida, point_right=idb,
                xi_left=xa, xi_right=xb, xi_midpoint=(xa+xb)/2, delta_xi_midpoint=(xa+xb)/2-float(a["xi_mott"]),
                xi_mott=float(a["xi_mott"]), interval_width=dx)
            for key in OBSERVABLES:
                s = log_slope(float(a[key]), float(b[key]), dx)
                s_error = (-math.log1p(-estimates[ida, key])-math.log1p(-estimates[idb, key])) / dx
                item = dict(info, observable=key, log_slope=s, estimated_s_interpolation_slope_error=s_error)
                slopes.append(item)
                curve_slopes.append(item)
            for sp in SPECIES:
                denominator = dx*logmean(float(a["tauinv_"+sp]), float(b["tauinv_"+sp]))
                values = []
                for ch in processes:
                    ca, cb = ci[ida, sp, ch], ci[idb, sp, ch]
                    rate_part, density_part = product_difference(
                        float(ca["collision_density_weight_inv_fm3"]), float(cb["collision_density_weight_inv_fm3"]),
                        float(ca["rate_fm2"]), float(cb["rate_fm2"]))
                    item = dict(info, species=sp, channel=ch,
                        log_slope_from_rate=-rate_part/denominator,
                        log_slope_from_density=-density_part/denominator,
                        log_slope_contribution=-(rate_part+density_part)/denominator)
                    values.append(item)
                    channel_derivatives.append(item)
                expected = log_slope(float(a["tau_"+sp]), float(b["tau_"+sp]), dx)
                difference = abs(sum(r["log_slope_contribution"] for r in values)-expected)
                closure_error = max(closure_error, difference)
                assert difference < 1e-7, (curve, sp, ida, idb, difference)
                for name, subset in (("K_s_processes", [r for r in values if r["channel"] in K_S_CHANNELS]),
                                     ("other_processes", [r for r in values if r["channel"] not in K_S_CHANNELS])):
                    group_derivatives.append(dict(info, species=sp, process_group=name,
                        **{key: sum(r[key] for r in subset) for key in (
                            "log_slope_from_rate", "log_slope_from_density", "log_slope_contribution")}))
            for key in FINAL:
                denominator = dx*logmean(float(a[key]), float(b[key]))
                values = []
                for sp in SPECIES:
                    relaxation, coefficient = product_difference(
                        float(wi[ida, sp][key+"_coefficient"]), float(wi[idb, sp][key+"_coefficient"]),
                        float(a["tau_"+sp]), float(b["tau_"+sp]))
                    values.append((relaxation+coefficient)/denominator)
                    final_derivatives.append(dict(info, observable=key, species=sp,
                        log_slope_from_tau=relaxation/denominator,
                        log_slope_from_coefficient=coefficient/denominator,
                        log_slope_contribution=(relaxation+coefficient)/denominator))
                difference = abs(sum(values)-log_slope(float(a[key]), float(b[key]), dx))
                closure_error = max(closure_error, difference)
                assert difference < 1e-7, (curve, key, difference)
        center = float(rows[0]["xi_mott"])
        central = next((r for r in rows if abs(float(r["xi"])-center)<1e-12), None)
        if central:
            for h in (.004, .002, .001, .0005, .0002, .0001):
                a = next((r for r in rows if abs(float(r["xi"])-(center-h))<1e-12), None)
                b = next((r for r in rows if abs(float(r["xi"])-(center+h))<1e-12), None)
                if a and b:
                    for key in OBSERVABLES:
                        left = log_slope(float(a[key]), float(central[key]), h)
                        right = log_slope(float(central[key]), float(b[key]), h)
                        centered.append(dict(curve_id=curve, meson=rows[0]["meson"], observable=key,
                            half_width=h, xi_mott=center, left_log_slope=left, right_log_slope=right,
                            right_minus_left=right-left))
        for key in OBSERVABLES:
            values = [r for r in curve_slopes if r["observable"] == key]
            if not values:
                continue
            central_values = [r for r in values if abs(r["delta_xi_midpoint"]) <= .002]
            outer_values = [r for r in values if abs(r["delta_xi_midpoint"]) >= .0025]
            contrasts.append(dict(curve_id=curve, meson=rows[0]["meson"], window=rows[0]["window"],
                observable=key, intervals=len(values), min_slope=min(r["log_slope"] for r in values),
                max_slope=max(r["log_slope"] for r in values),
                absolute_slope_range=max(r["log_slope"] for r in values)-min(r["log_slope"] for r in values),
                central_max_abs_slope=max((abs(r["log_slope"]) for r in central_values), default=None),
                outer_max_abs_slope=max((abs(r["log_slope"]) for r in outer_values), default=None),
                max_estimated_s_interpolation_slope_error=max(r["estimated_s_interpolation_slope_error"] for r in values)))
    return dict(slopes=slopes, channel_derivatives=channel_derivatives,
                final_species_derivatives=final_derivatives, process_group_derivatives=group_derivatives,
                contrasts=contrasts, centered_slopes=centered), closure_error


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--runs", required=True)
    parser.add_argument("--comparison-runs", default="")
    parser.add_argument("--preview", action="store_true")
    args = parser.parse_args()
    obs, rates, weights, channels, provenance = collect(args.case, args.runs.split(","), args.preview)
    if not obs:
        raise ValueError("No complete point yet")
    products, closure = analyze(obs, rates, weights, channels)
    point_checks, slope_checks, xi_checks = [], [], []
    reference_evaluations = new_xi_evaluations = 0
    if args.comparison_runs:
        references = collect(args.case, args.comparison_runs.split(","), args.preview)
        reference_evaluations = len(references[0])
        ref_products, ref_closure = analyze(*references[:4])
        closure = max(closure, ref_closure)
        provenance.update(references[4])
        by_id = {r["point_id"]: r for r in obs}
        ref_by_id = {r["point_id"]: r for r in references[0]}
        new_xi_evaluations = len(set(ref_by_id)-set(by_id))
        for row in references[0]:
            if row["point_id"] not in by_id:
                continue
            original = by_id[row["point_id"]]
            assert all(row[k] == original[k] for k in ("curve_id", "meson", "window"))
            assert all(abs(float(row[k])-float(original[k])) < 1e-12
                       for k in ("xi", "xi_mott", "T_MeV", "muB_MeV")), row["point_id"]
            for key in OBSERVABLES:
                a, b = float(by_id[row["point_id"]][key]), float(row[key])
                point_checks.append(dict(curve_id=row["curve_id"], point_id=row["point_id"], observable=key,
                    fine_value=a, tighter_value=b, relative_difference=(b-a)/b))
        slope_checks = compare_native_slopes(products["slopes"], ref_by_id)
        refinement_path = args.case/"convergence_plan_v1/summary.json"
        if refinement_path.exists():
            for edge in json.loads(refinement_path.read_text(encoding="utf-8-sig"))["refinements"]:
                a, m, b = (ref_by_id[edge[name]] for name in ("left_id", "midpoint_id", "right_id"))
                xa, xm, xb = (float(r["xi"]) for r in (a, m, b))
                for key in OBSERVABLES:
                    coarse = log_slope(float(a[key]), float(b[key]), xb-xa)
                    left = log_slope(float(a[key]), float(m[key]), xm-xa)
                    right = log_slope(float(m[key]), float(b[key]), xb-xm)
                    xi_checks.append(dict(curve_id=edge["curve_id"], observable=key,
                        left_id=edge["left_id"], midpoint_id=edge["midpoint_id"], right_id=edge["right_id"],
                        xi_left=xa, xi_midpoint=xm, xi_right=xb, coarse_log_slope=coarse,
                        refined_left_log_slope=left, refined_right_log_slope=right,
                        max_split_slope_change=max(abs(left-coarse), abs(right-coarse))))
        products.update({"reference_"+k: v for k, v in ref_products.items()})
    args.output.mkdir(parents=True, exist_ok=False)
    for key, rows in dict(products, point_convergence=point_checks, slope_convergence=slope_checks,
                         xi_resolution=xi_checks).items():
        write_csv(args.output / (key+".csv"), rows)
    summary = dict(status="partial_diagnostic_preview" if args.preview else "completed_diagnostic_review",
        points=len(obs), curves=len({r["curve_id"] for r in obs}), channels=len(rates),
        exact_decomposition_max_absolute_closure_error=closure, provenance=provenance,
        contrasts=products["contrasts"], comparison_points=len({r["point_id"] for r in point_checks}),
        reference_evaluations=reference_evaluations, new_xi_evaluations=new_xi_evaluations,
        matched_slope_comparisons=len(slope_checks),
        xi_resolution_checks=len(xi_checks),
        max_point_relative_difference=max((abs(r["relative_difference"]) for r in point_checks), default=None),
        max_slope_precision_difference=max((r["absolute_difference"] for r in slope_checks), default=None),
        statement="Finite log secants. Estimated bounds cover s interpolation only. Process groups are not meson amplitude ablations.",
        formula="S_tau=-sum_c Delta(n_eff*w_c)/(Delta_xi*logmean(R_left,R_right)); exact midpoint product decomposition.",
        script_sha256=sha(Path(__file__)))
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, allow_nan=False))


if __name__ == "__main__":
    main()
