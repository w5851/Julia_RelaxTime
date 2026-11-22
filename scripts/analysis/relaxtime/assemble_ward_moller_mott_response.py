"""Audit and freeze the remotely computed Mott windows for the appendix."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

CORE = "096199b96d34045bf8e2d9611c2d786d6f92fb8e"
PROFILE = "a26d470cf09730ba79d4a440bf9e3d7bc17c202f59323591edafe4a47d249313"
MESONS = ("pi", "K", "sigma")
TIERS = ("fine", "reference")
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


def collect_transport(root, centers):
    grids = {(m, t): {} for m in MESONS for t in TIERS}
    hashes = {}
    for path in sorted(root.rglob("*.json")):
        if not path.stem.isdecimal() or path.parent.name not in MESONS:
            continue
        row = read_json(path)
        meta = row["metadata"]
        check_meta(meta)
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
        if not (meta["runtime_sampling"]["rtol"] == {"fine": 2e-6, "reference": 2e-7}[tier]
                and len(row["rates"]) == len(row["rate_diagnostics"]) == 17
                and all(math.isfinite(v) and v > 0 for v in row["rates"].values())):
            raise ValueError(f"Incomplete rate set or wrong sampling tier: {path}")
        if abs(row["m_u_MeV"]-loc["center_masses_inv_fm"][0]*197.3269804) > 20.0:
            raise ValueError(f"Unexpected phase jump: {path}")
        grid[index] = {k:row[k] for k in ("meson","tier","index","T_MeV","muB_MeV","xi","xi_mott",
            "delta_xi","eta_over_s","m_u_MeV","m_s_MeV","residual")}
        grid[index]["driver_sha256"] = meta["driver_sha256"]
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
                 "slope_min":min(r["slope"] for r in selected),
                 "slope_max":max(r["slope"] for r in selected)}
        if meson == "sigma":
            sidebands = [r for r in selected if .003 <= abs(r["u_midpoint"]) <= .008]
            a,b = line_fit([r["u_midpoint"] for r in sidebands],[r["slope"] for r in sidebands])
            for r in selected:
                r["linear_background_residual"] = r["slope"]-a-b*r["u_midpoint"]
            inset = [r["linear_background_residual"] for r in selected if abs(r["u_midpoint"])<=.001]
            entry["inset"] = {"fit_range":[.003,.008],"intercept":a,"gradient":b,"display_half_width":.001,
                              "residual_min":min(inset),"residual_max":max(inset)}
        checks[meson] = entry
        slopes += selected
    return slopes, checks


def scalar_results(root, center):
    tiers, provenance = {}, {}
    for t in TIERS:
        folder = root/f"scalar_{t}"
        p = folder/"summary.json"
        meta = read_json(p)
        check_meta(meta["metadata"])
        if not (meta["status"] == "complete" and meta["determinant_floor"] is False
                and meta["lim"]["resonance_sign"] == -1 and abs(meta["center"]-center)<1e-9):
            raise ValueError("Unexpected scalar integration provenance")
        for filename, expected in meta["files"].items():
            if digest(folder/filename) != expected:
                raise ValueError(f"Changed scalar output: {filename}")
        tiers[t] = read_rows(folder/"curvatures.csv")
        provenance[t] = {"sha256":digest(p), **meta}
    lim = provenance["fine"]["lim"]["expected"]
    references = {r["h"]:r for r in tiers["reference"]}
    comparisons = []
    for r in tiers["fine"]:
        h = r["h"]
        row = {**r, "analytic_limit":lim, "reference_relative_defect":""}
        if h in references:
            s = references[h]
            row["reference_relative_defect"] = abs(r["scaled_curvature"]-s["scaled_curvature"])/lim
        error = (r["step_defect"]+r["quadrature_defect"])*math.sqrt(h)/lim
        row["usable"] = bool(r["xr"] > 0 and error <= .025
                             and (row["reference_relative_defect"] == "" or row["reference_relative_defect"] <= .025))
        comparisons.append(row)
    plateau = [r for r in comparisons if r["usable"] and abs(r["relative_to_limit"]) <= .05]
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
        "max_reference_defect":max(r["reference_relative_defect"] for r in plateau
                                   if r["reference_relative_defect"] != "")}
    return comparisons, checks, provenance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan-root",required=True,type=Path)
    parser.add_argument("--scalar-root",required=True,type=Path)
    parser.add_argument("--out",required=True,type=Path)
    parser.add_argument("--action-url",required=True)
    args = parser.parse_args()
    centers = {}
    for meson in MESONS:
        path, = args.scan_root.rglob(f"center_{meson}.json")
        centers[meson] = {"sha256":digest(path), **read_json(path)}
        check_meta(centers[meson]["metadata"])
    grids, node_hashes = collect_transport(args.scan_root, centers)
    slopes, slope_checks = local_slopes(grids)
    scalar, scalar_checks, scalar_provenance = scalar_results(args.scalar_root,centers["sigma"]["location"]["xi"])
    args.out.mkdir(parents=True, exist_ok=False)
    tables = {"transport_nodes.csv":[r for g in grids.values() for r in g],
              "local_slopes.csv":slopes,"scalar_curvature.csv":scalar}
    for name,rows in tables.items():
        write_rows(args.out/name,rows)
    result = {"status":"frozen_analysis","core_commit":CORE,"profile_sha256":PROFILE,
              "action_url":args.action_url,"centers":centers,"transport_checks":slope_checks,
              "scalar_checks":scalar_checks,"scalar_provenance":scalar_provenance,
              "transport_artifact_sha256":node_hashes,
              "postprocess_sha256":digest(Path(__file__)),
              "files":{name:digest(args.out/name) for name in tables}}
    write_json(args.out/"numerical_manifest.json",result)
    print(json.dumps(scalar_checks,indent=2))


if __name__ == "__main__":
    main()
