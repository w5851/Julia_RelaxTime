"""Normalize retained published midrapidity observations; no network or fitting."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INPUT = ROOT / "data/inputs/meson_freezeout_experiment_20261017"


def errors_quadrature(errors, *, exclude=()):
    """Keep asymmetric signs out of variances; do not double-add subsets."""
    lower = upper = 0.0
    for error in errors:
        if error.get("hide") or error.get("label", "") in exclude:
            continue
        if "symerror" in error:
            lo = hi = abs(float(error["symerror"]))
        else:
            lo = abs(float(error["asymerror"]["minus"]))
            hi = abs(float(error["asymerror"]["plus"]))
        lower += lo * lo
        upper += hi * hi
    return math.sqrt(lower), math.sqrt(upper)


def prepare(input_root: Path):
    source_path = input_root / "sources.json"
    sources = json.loads(source_path.read_text(encoding="utf-8"))
    for source in sources["sources"]:
        actual = hashlib.sha256((input_root / source["path"]).read_bytes()).hexdigest()
        if actual != source["sha256"]:
            raise ValueError(f"source hash mismatch: {source['path']}")
    read = lambda name: json.loads((input_root / "raw" / name).read_text(encoding="utf-8"))
    rows = []

    def add(identifier, energy, ratio, errors, *, system, experiment, centrality,
            rapidity, source, charge="K+/pi+", method="published_ratio",
            error_method="published_stat_sys_quadrature", primary=False,
            horn=False, notes=""):
        if not (math.isfinite(energy) and energy > 0 and math.isfinite(ratio) and ratio >= 0):
            raise ValueError("invalid experimental energy/ratio")
        rows.append(dict(id=identifier, sqrt_s_NN_GeV=energy, ratio=ratio,
            error_minus=errors[0], error_plus=errors[1], system=system,
            experiment=experiment, centrality=centrality, rapidity=rapidity,
            charge=charge, method=method, error_method=error_method,
            primary_auau=primary, auxiliary_horn=horn, source=source, notes=notes))

    star = read("star2017_kpi_hepdata.json")
    groups = star["qualifiers"]["Sqrt{s_{NN}} [GeV]"]
    plus = next(q["group"] for q in groups if q["value"] == "ka+/pi+")
    for row in star["values"]:
        value = next(y for y in row["y"] if y["group"] == plus)
        energy = float(row["x"][0]["value"])
        add(f"STAR_BES_{energy:g}", energy, float(value["value"]), errors_quadrature(value["errors"]),
            system="Au+Au", experiment="STAR", centrality="0-5%", rapidity="|y|<0.1",
            source=star["doi"], primary=True)

    na49 = read("na492008_kpi_hepdata.json")
    if na49["headers"][1]["name"] != "MULT(Q=K+)/MULT(Q=PI+)":
        raise ValueError("NA49 first dependent variable must be positive-charge ratio")
    for row in na49["values"]:
        energy = float(row["x"][0]["value"])
        if energy >= 100:  # replace historical RHIC compilation with STAR table
            continue
        value = next(y for y in row["y"] if y["group"] == 0)
        ags = energy < 6
        add(f"{'AGS' if ags else 'NA49'}_{energy:g}", energy, float(value["value"]),
            errors_quadrature(value["errors"]), system="Au+Au" if ags else "Pb+Pb",
            experiment="AGS" if ags else "NA49", centrality="0-5%" if ags else "0-7.2%",
            rapidity="midrapidity" if ags else "K: |y|<0.1; pi: 0<y<0.2",
            source=na49["doi"], primary=ags, horn=not ags,
            error_method="published_combined" if ags else "published_stat_sys_quadrature",
            notes="compiled AGS point; published energy retained" if ags else
                "shared systematic errors across SPS energies; no covariance matrix")

    high = read("star2009_table_ix.json")
    if high["observable"] != "Kplus_over_piplus":
        raise ValueError("wrong STAR charge column")
    for row in high["rows"]:
        energy = row["sqrt_s_NN_GeV"]
        add(f"STAR_{energy:g}", energy, row["ratio"], (row["error"], row["error"]),
            system="Au+Au", experiment="STAR", centrality=row["centrality"], rapidity="|y|<0.1",
            source="10.1103/PhysRevC.79.034909 Table IX", primary=True,
            error_method="published_combined", method="published_table_transcription")

    pion = read("star145_piplus_hepdata.json")["values"][0]
    kaon = read("star145_kplus_hepdata.json")["values"][0]
    if pion["x"][0]["value"] != "338" or kaon["x"][0]["value"] != "338":
        raise ValueError("14.5 GeV yield rows must be the central Npart=338 rows")
    pi_y, k_y = pion["y"][0], kaon["y"][0]
    pi_n, k_n = float(pi_y["value"]), float(k_y["value"])
    ratio = k_n / pi_n
    pi_err, k_err = errors_quadrature(pi_y["errors"])[0], errors_quadrature(k_y["errors"])[0]
    err = ratio * math.hypot(pi_err / pi_n, k_err / k_n)
    add("STAR_14.5", 14.5, ratio, (err, err), system="Au+Au", experiment="STAR",
        centrality="0-5%", rapidity="|y|<0.1", source="10.1103/PhysRevC.101.024905 dN/dy",
        method="ratio_of_published_yields", error_method="zero_covariance_propagation",
        primary=True, notes="26.7/141; yield covariance unavailable; uncertainty is derived")

    for name, energy in (("alice2013_kpi_hepdata.json", 2760.), ("alice2020_kpi_hepdata.json", 5020.)):
        data = read(name)
        # Group 0 is charged multiplicity, not the ratio. Group 1 is charge sum.
        value = next(y for y in data["values"][0]["y"] if y["group"] == 1)
        add(f"ALICE_{energy:g}", energy, float(value["value"]),
            errors_quadrature(value["errors"], exclude=("syst. uncorr.",)),
            system="Pb+Pb", experiment="ALICE", centrality="0-5%", rapidity="|y|<0.5",
            source=data["doi"], charge="(K++K-)/(pi++pi-)",
            notes="auxiliary charge-sum comparison near muB=0; separately listed uncorrelated systematic retained in raw source, not added again")

    with (input_root / "hades_digitized/points.csv").open(encoding="utf-8", newline="") as handle:
        point = next(csv.DictReader(handle))
    if point["status"] != "measured":
        raise ValueError("HADES point was not resolved")
    ratio = float(point["y"])
    add("HADES_2.4", 2.4, ratio, (ratio-float(point["y_lower"]), float(point["y_upper"])-ratio),
        system="Au+Au", experiment="HADES", centrality="central (source compilation)", rapidity="y approx 0",
        source="arXiv:2205.06726v5 Figure 1 left, SIS", method="digitized",
        error_method="pixel_reading_interval_not_experimental_error",
        notes="display only; nominal energy 2.4 GeV; excluded from numerical ranking")
    rows.sort(key=lambda r: (r["sqrt_s_NN_GeV"], r["id"]))
    if len({row["id"] for row in rows}) != len(rows):
        raise ValueError("duplicate experimental ids")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rows = prepare(args.input_root)
    output = args.output or args.input_root / "experimental_ratios.csv"
    if output.exists():
        raise FileExistsError(f"refusing to overwrite experimental data: {output}")
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    print(json.dumps({"rows": len(rows), "primary_AuAu": sum(r["primary_auau"] for r in rows),
        "auxiliary_NA49": sum(r["auxiliary_horn"] for r in rows), "output": str(output)}))


if __name__ == "__main__":
    main()
