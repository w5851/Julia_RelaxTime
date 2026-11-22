"""Audit and freeze the remotely computed Mott windows for the appendix."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path

CORE = "096199b96d34045bf8e2d9611c2d786d6f92fb8e"
PROFILE = "a26d470cf09730ba79d4a440bf9e3d7bc17c202f59323591edafe4a47d249313"
MESONS = ("pi", "K", "sigma")
TIERS = ("fine", "reference")
RATE_CONVERGENCE = "thermal_defect_and_two_shifted_audits"
# Units of 0.000125 in xi. These mirror the driver and require complete grids.
STEPS = {
    "fine": sorted(set(range(-64, 65, 8)) | set(range(-16, 17, 2)) | {-3, -1, 1, 3}),
    "reference": [-64, -56, -8, -2, -1, 0, 1, 2, 8, 56, 64],
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return [dict((k, float(v)) for k, v in row.items()) for row in csv.DictReader(stream)]


def write_rows(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, data):
    path.write_bytes((json.dumps(data, indent=2, allow_nan=False) + "\n").encode("utf-8"))


def check_meta(meta):
    p = meta["profile"]
    if not (meta["core_commit"] == CORE and p["source_sha256"] == PROFILE
            and meta["julia_version"] == "1.12.5"
            and p["a_polarization"] == p["a_coupling"] == "ward_matched_RS"
            and p["b0_direction"] == "sphere_average"
            and meta["runtime_sampling"]["relative_velocity"] == "medium_moller"
            and p["background"]["p_nodes"] == 96 and p["background"]["cos_nodes"] == 24
            and p["background"]["thermal_p_max_inv_fm"] == 25.0):
        raise ValueError("Numerical record is not from the pinned Ward/Møller prescription")


def secant(left, right):
    a, b = left["eta_over_s"], right["eta_over_s"]
    dx = right["delta_xi"] - left["delta_xi"]
    if not (math.isfinite(a) and math.isfinite(b) and min(a, b, dx) > 0):
        raise ValueError("Secants require positive eta/s and increasing abscissae")
    return math.log1p((b-a)/a)/dx


def line_fit(xs, ys):
    if len(xs) != len(ys) or len(xs) < 2 or not all(map(math.isfinite, [*xs, *ys])):
        raise ValueError("Incomplete linear model")
    xbar, ybar = math.fsum(xs)/len(xs), math.fsum(ys)/len(ys)
    slope = math.fsum((x-xbar)*(y-ybar) for x, y in zip(xs, ys))/math.fsum((x-xbar)**2 for x in xs)
    return ybar-slope*xbar, slope


def collect_transport(root, centers, sources=None):
    grids = {(m, t): {} for m in MESONS for t in TIERS}
    hashes = {}
    for path in sorted(root.rglob("*.json")):
        if not path.stem.isdecimal() or path.parent.name not in MESONS:
            continue
        row = read_json(path)
        meta = row["metadata"]
        check_meta(meta)
        if sources is not None:
            prefix = path.relative_to(root).parts[0]
            source = sources.get(prefix, {})
            if (source.get("driver_sha256") != meta["driver_sha256"]
                    or (source.get("host") and source["host"] != meta["run_host"])):
                raise ValueError(f"Transport node does not match the recorded source: {path}")
        meson, tier, index = row["meson"], row["tier"], row["index"]
        grid = grids[meson, tier]
        if index in grid:
            raise ValueError(f"Duplicate node: {meson} {tier} {index}")
        if not (1 <= index <= len(STEPS[tier]) and
                math.isclose(row["delta_xi"], STEPS[tier][index-1]*0.000125, abs_tol=1e-13) and
                math.isclose(row["xi"], row["xi_mott"]+row["delta_xi"], abs_tol=1e-13)):
            raise ValueError(f"Unexpected scan coordinates: {path}")
        loc = centers[meson]["location"]
        if not (abs(row["xi_mott"]-loc["xi"]) < 1e-11 and row["T_MeV"] == loc["T_MeV"]
                and row["muB_MeV"] == loc["muB_MeV"] and row["residual"] <= 1e-11):
            raise ValueError(f"Background mismatch: {path}")
        rates, diagnostics = row["rates"], row["rate_diagnostics"]
        cfg = meta["runtime_sampling"]
        recheck = (tier == "reference" and cfg["rtol"] == 1e-6 and cfg["max_evaluations"] == 300_000
                   and cfg["max_nodes"] == 80_000 and cfg["max_rounds"] == 90)
        if not ((recheck or cfg["rtol"] == {"fine": 2e-6, "reference": 2e-7}[tier])
                and len(rates) == len(diagnostics) == 17
                and set(rates) == {d["channel"] for d in diagnostics}
                and all(math.isfinite(v) and v > 0 for v in rates.values())
                and all(d["reason"] == RATE_CONVERGENCE and d["guard_phase"] == 2
                        and d["rate"] == rates[d["channel"]]
                        and 0 <= d["estimated_error"] <= cfg["atol"] + cfg["rtol"]*d["rate"]
                        for d in diagnostics)):
            raise ValueError(f"Incomplete rate set or wrong sampling tier: {path}")
        if abs(row["m_u_MeV"]-loc["center_masses_inv_fm"][0]*197.3269804) > 20.0:
            raise ValueError(f"Unexpected phase jump: {path}")
        grid[index] = {k:row[k] for k in ("meson","tier","index","T_MeV","muB_MeV","xi","xi_mott",
            "delta_xi","eta_over_s","m_u_MeV","m_s_MeV","residual")}
        grid[index]["driver_sha256"] = meta["driver_sha256"]
        grid[index]["sampling_rtol"] = cfg["rtol"]
        grid[index]["sampling_recheck"] = recheck
        grid[index]["run_host"] = meta["run_host"]
        if not recheck and grid[index]["driver_sha256"] != centers[meson]["metadata"]["driver_sha256"]:
            raise ValueError(f"Scan and center used different drivers: {path}")
        hashes[path.relative_to(root).as_posix()] = digest(path)
    missing = {f"{m}/{t}": sorted(set(range(1, len(STEPS[t])+1))-set(grids[m,t]))
               for m,t in grids if len(grids[m,t]) != len(STEPS[t])}
    if missing:
        raise ValueError(f"Missing independent transport nodes: {missing}")
    return {mt:list(dict(sorted(grid.items())).values()) for mt,grid in grids.items()}, hashes


def local_slopes(grids):
    slopes, checks = [], {}
    for meson in MESONS:
        fine, ref = [grids[meson,tier] for tier in TIERS]
        fine_at = {round(r["delta_xi"],12):r for r in fine}
        eta_defect = max(abs(r["eta_over_s"] / fine_at[round(r["delta_xi"], 12)]["eta_over_s"] - 1)
                         for r in ref)
        sampling = []
        for a,b in zip(ref[:-1],ref[1:]):
            d = secant(a,b)
            native = secant(fine_at[round(a["delta_xi"],12)],fine_at[round(b["delta_xi"],12)])
            sampling.append({"u_left":a["delta_xi"],"u_right":b["delta_xi"],"slope":native,
                             "reference_slope":d,"absolute_defect":abs(d-native)})
        selected = []
        for a,b in zip(fine[:-1],fine[1:]):
            aa, bb = a["delta_xi"],b["delta_xi"]
            selected.append({"meson":meson,"u_left":aa,"u_right":bb,"u_midpoint":(aa+bb)/2,
                             "slope":secant(a,b),"linear_background_residual":0.0})
        around = []
        for sign in [-1,1]:
            z = fine_at[0.0]
            outer = fine_at[sign*0.00025]
            inner = fine_at[sign*0.000125]
            ordered = (lambda a,b:(a,b)) if sign>0 else (lambda a,b:(b,a))
            d = secant(*ordered(z,outer))
            d2 = secant(*ordered(z,inner))
            around.append({"side":sign,"interval":0.00025,"slope":d,"half_step_slope":d2,
                           "difference":d2-d})
        entry = {"same_endpoint_sampling":sampling,"adjacent_to_mott":around,
                 "max_shared_node_eta_relative_defect":eta_defect,
                 "max_shared_endpoint_slope_defect":max(r["absolute_defect"] for r in sampling),
                 "slope_min":min(r["slope"] for r in selected),
                 "slope_max":max(r["slope"] for r in selected)}
        if meson == "sigma":
            sidebands = [r for r in selected if .003 <= abs(r["u_midpoint"]) <= .008]
            a,b = line_fit([r["u_midpoint"] for r in sidebands],[r["slope"] for r in sidebands])
            for r in selected:
                r["linear_background_residual"] = r["slope"]-a-b*r["u_midpoint"]
            # Show the central region covered entirely by the refined xi grid.
            half_width = .0005
            inset = [r["linear_background_residual"] for r in selected if abs(r["u_midpoint"])<=half_width]
            entry["inset"] = {"fit_range":[.003,.008],"intercept":a,"gradient":b,"display_half_width":half_width,
                              "residual_min":min(inset),"residual_max":max(inset)}
        checks[meson] = entry
        slopes += selected
    return slopes, checks


def scalar_results(root, center):
    tiers, provenance, stencils = {}, {}, {}
    for t in TIERS:
        folder = root/f"scalar_{t}"
        p = folder/"summary.json"
        meta = read_json(p)
        check_meta(meta["metadata"])
        if not (meta["status"] == "complete" and meta["determinant_floor"] is False
                and "strict_center_bracket" in meta
                and meta["lim"]["resonance_sign"] == -1 and abs(meta["center"]-center)<1e-9):
            raise ValueError("Unexpected scalar integration provenance")
        for filename, expected in meta["files"].items():
            if digest(folder/filename) != expected:
                raise ValueError(f"Changed scalar output: {filename}")
        tiers[t] = read_rows(folder/"curvatures.csv")
        stencils[t] = read_rows(folder/"integrals.csv")
        if not all(math.isfinite(v) for rows in (tiers[t],stencils[t]) for row in rows for v in row.values()):
            raise ValueError(f"Nonfinite scalar data in {t}")
        meta["max_background_residual"] = max(r["background_residual"] for r in stencils[t])
        if meta["max_background_residual"] >= 1e-12:
            raise ValueError(f"Unresolved scalar mean-field stencil in {t}")
        verification = folder/"source_verification.json"
        if verification.is_file():
            verified = read_json(verification)
            if verified["core_commit"] != CORE or verified["runtime_files_matched"] <= 0:
                raise ValueError(f"Unverified remote scalar core in {t}")
            meta["remote_source_verification"] = verified
        provenance[t] = {"sha256":digest(p), **meta}
    lim = provenance["fine"]["lim"]["expected"]
    if not (lim > 0 and math.isclose(lim, provenance["reference"]["lim"]["expected"], rel_tol=1e-5)
            and provenance["fine"]["scalar_driver_sha256"] == provenance["reference"]["scalar_driver_sha256"]):
        raise ValueError("Scalar tiers must share the same driver and independent limit")
    references = {r["h"]:r for r in tiers["reference"]}
    comparisons = []
    segments, consecutive = [], []
    for r in sorted(tiers["fine"], key=lambda r: r["h"], reverse=True):
        h = r["h"]
        row = {**r, "analytic_limit":lim, "reference_scaled_curvature":"", "reference_relative_defect":""}
        if h in references:
            s = references[h]
            row["reference_scaled_curvature"] = s["scaled_curvature"]
            row["reference_relative_defect"] = abs(r["scaled_curvature"]-s["scaled_curvature"])/lim
        error = (r["step_defect"]+r["quadrature_defect"])*math.sqrt(h)/lim
        row["relative_numerical_defect"] = error
        row["usable"] = bool(r["xr"] > 0 and error <= .025
                             and (row["reference_relative_defect"] == "" or row["reference_relative_defect"] <= .025))
        comparisons.append(row)
        if row["usable"] and abs(row["relative_to_limit"]) <= .05:
            consecutive.append(row)
        else:
            segments.append(consecutive)
            consecutive = []
    plateau = max([*segments, consecutive], key=len)
    if len(plateau)<3 or max(r["h"] for r in plateau) < 4*min(r["h"] for r in plateau):
        raise ValueError("The independent scalar curvature plateau is not resolved over a factor four in h")
    if len([r for r in plateau if r["reference_relative_defect"] != ""])<2:
        raise ValueError("The scalar plateau needs two common precision points")
    curvature_exponent = line_fit([math.log(r["h"]) for r in plateau],
                                  [math.log(r["curvature"]) for r in plateau])[1]
    area_exponent = line_fit([math.log(r["xr"]) for r in plateau],
                             [math.log(r["peak"]) for r in plateau])[1]
    if abs(curvature_exponent+.5)>.06 or abs(area_exponent-1.5)>.06:
        raise ValueError("Plateau and peak-area exponents are inconsistent with the threshold law")
    checks = {"analytic_limit_fm2":lim,"h_min":min(r["h"] for r in plateau),
        "h_max":max(r["h"] for r in plateau),"fit_points":len(plateau),
        "curvature_exponent":curvature_exponent,"peak_area_exponent":area_exponent,
        "max_relative_to_limit":max(abs(r["relative_to_limit"]) for r in plateau),
        "nearest_scaled_curvature_fm2":plateau[-1]["scaled_curvature"],
        "nearest_relative_to_limit":plateau[-1]["relative_to_limit"],
        "max_numerical_defect":max(r["relative_numerical_defect"] for r in plateau),
        "max_reference_defect":max(r["reference_relative_defect"] for r in plateau
                                   if r["reference_relative_defect"] != "")}
    return comparisons, checks, provenance, stencils


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan-root",required=True,type=Path)
    parser.add_argument("--scalar-root",required=True,type=Path)
    parser.add_argument("--out",required=True,type=Path)
    parser.add_argument("--action-url",required=True)
    parser.add_argument("--action-head-sha",required=True)
    parser.add_argument("--source-index",type=Path,
                        help="Source records keyed by top-level folder when combining runs")
    args = parser.parse_args()
    if not re.fullmatch(r"[0-9a-f]{40}", args.action_head_sha):
        parser.error("--action-head-sha must be a full commit hash")
    centers = {}
    for meson in MESONS:
        path, = args.scan_root.rglob(f"center_{meson}.json")
        centers[meson] = {"sha256":digest(path), **read_json(path)}
        check_meta(centers[meson]["metadata"])
    sources = read_json(args.source_index) if args.source_index else None
    grids, node_hashes = collect_transport(args.scan_root, centers, sources)
    slopes, slope_checks = local_slopes(grids)
    scalar, scalar_checks, scalar_provenance, stencils = scalar_results(
        args.scalar_root, centers["sigma"]["location"]["xi"])
    args.out.mkdir(parents=True, exist_ok=False)
    tables = {"transport_nodes.csv":[r for g in grids.values() for r in g],
              "local_slopes.csv":slopes,"scalar_curvature.csv":scalar,
              **{f"scalar_stencil_{tier}.csv": rows for tier, rows in stencils.items()}}
    for name,rows in tables.items():
        write_rows(args.out/name,rows)
    result = {"status":"frozen_analysis","core_commit":CORE,"profile_sha256":PROFILE,
              "action_url":args.action_url,"action_head_sha":args.action_head_sha,
              "transport_sources":sources,
              "centers":centers,"transport_checks":slope_checks,
              "scalar_checks":scalar_checks,"scalar_provenance":scalar_provenance,
              "transport_artifact_sha256":node_hashes,
              "postprocess_sha256":digest(Path(__file__)),
              "files":{name:digest(args.out/name) for name in tables}}
    write_json(args.out/"numerical_manifest.json",result)
    print(json.dumps(scalar_checks,indent=2))


if __name__ == "__main__":
    main()
