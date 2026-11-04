"""Check the common thermal-rate normalization independently of Julia quadrature."""
from __future__ import annotations

import argparse
import csv
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import shutil

import numpy as np

SPECIES = ("u", "d", "s", "ubar", "dbar", "sbar")


def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@lru_cache(None)
def grid(order, scale, semi=True):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    t, w = (nodes + 1) / 2, weights / 2
    if not semi:
        return t * scale, w * scale
    mask = t < 0.9999  # Same negligible far-tail exclusion as the runtime mapping.
    t, w = t[mask], w[mask]
    return scale * t / (1 - t), w * scale / (1 - t) ** 2


def density(row, species, order, angles, scale, *, semi=True, xi=None):
    p, pw = grid(order, scale, semi)
    c, cw = grid(angles, 1.0, False)  # Even RS integrand: cos(theta) in [0,1].
    xi = float(row["xi"]) if xi is None else xi
    energy = np.sqrt(float(row[f"m_{species[0]}"]) ** 2 + p[:, None] ** 2 * (1 + xi * c[None, :] ** 2))
    anti = species.endswith("bar")
    phi, phibar = float(row["Phi"]), float(row["Phibar"])
    if anti:
        phi, phibar = phibar, phi
    x = (energy + (1 if anti else -1) * float(row["muq_fm"])) / float(row["T_fm"])
    # Stable rational form on each sign of (E -/+ mu)/T.
    z = np.exp(-np.abs(x))
    f_positive = (phi*z + 2*phibar*z*z + z**3) / (1 + 3*phi*z + 3*phibar*z*z + z**3)
    f_negative = (phi*z*z + 2*phibar*z + 1) / (z**3 + 3*phi*z*z + 3*phibar*z + 1)
    occupation = np.where(x >= 0, f_positive, f_negative)
    return float(6 / (2 * np.pi**2) * np.sum(pw[:, None] * p[:, None] ** 2 * cw[None, :] * occupation))


def verify(base, outdir):
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError("Use a new empty output directory")
    inputs_path = base / "convergence_selection_v1/review_points.csv"
    inputs = {r["point_id"]: r for r in rows(inputs_path)}
    historical_path = base / "density_qualification_v1/densities.csv"
    historical_errors = []
    for row in rows(historical_path):
        if row["configuration"] == "semi32_scale10_c20":
            value = density(inputs[row["point_id"]], row["species"], 32, 20, 10.0)
        elif row["configuration"] == "finite10_p64_c8":
            value = density(inputs[row["point_id"]], row["species"], 64, 8, 10.0, semi=False)
        else:
            continue
        historical_errors.append(abs(value / float(row["density"]) - 1))
    if not historical_errors or max(historical_errors) > 1e-10:
        raise ValueError("Independent Python density does not reproduce the stored Julia checks")

    sources = [inputs_path, historical_path]
    states = {}
    for name, reference_only in (("joint_native_pilot_p48_v1", False), ("joint_native_pilot_v2", True),
                                ("joint_native_representatives_p48_s0_recovered_v1", False),
                                ("joint_native_representatives_p48_s1_recovered_v1", False)):
        path = base / name / "background.csv"
        sources.append(path)
        for row in rows(path):
            if reference_only and row["configuration"] != "reference":
                continue
            key = row["point_id"], row["configuration"]
            if key in states:
                raise ValueError(f"Duplicate background: {key}")
            states[key] = dict(inputs[row["point_id"]], **row)
    if len(states) != 24 or len({p for p, _ in states}) != 12:
        raise ValueError("Expected both backgrounds at all 12 joint points")
    results, bounds = [], []
    for (point, configuration), row in states.items():
        relative = []
        for species in SPECIES:
            actual = density(row, species, 32, 36 if configuration == "candidate" else 48, 10.0)
            reference = density(row, species, 256, 1, 5.0, xi=0.0) / math.sqrt(1 + float(row["xi"]))
            check = density(row, species, 128, 1, 5.0, xi=0.0) / math.sqrt(1 + float(row["xi"]))
            if not all(math.isfinite(v) and v > 0 for v in (actual, reference, check)):
                raise ValueError("Invalid independent density")
            reference_change = check / reference - 1
            if abs(reference_change) > 1e-9:
                raise ValueError("Independent radial density reference is not converged")
            relative.append(actual / reference - 1)
            results.append(dict(point_id=point, configuration=configuration, species=species,
                                actual_density=actual, reference_density=reference,
                                relative_change=relative[-1], reference_relative_change=reference_change))
        # With a fixed positive numerator, these denominators bound every
        # incoming pair's normalization error and hence any positive rate sum.
        pair_bound = max(abs(1 / ((1+a)*(1+b)) - 1) for a in relative for b in relative)
        bounds.append(dict(point_id=point, configuration=configuration, maximum_pair_rate_error=pair_bound))
    worst = max(abs(r["relative_change"]) for r in results)
    bound = max(r["maximum_pair_rate_error"] for r in bounds)
    summary = dict(analysis_only=True, backgrounds=len(states), comparisons=len(results),
                   reference="n(xi)=n(0)/sqrt(1+xi), radial semi-infinite p256 scale5, checked against p128",
                   historical_julia_comparisons=len(historical_errors),
                   historical_julia_max_relative=max(historical_errors),
                   maximum_density_relative=worst, maximum_pair_normalization_relative=bound,
                   normalization_budget=0.001, passed=bound <= 0.001,
                   sources=[dict(path=str(p.resolve()), sha256=digest(p)) for p in sources], numpy=np.__version__)
    outdir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(__file__, outdir / Path(__file__).name)
    for name, data in (("densities.csv", results), ("pair_bounds.csv", bounds)):
        with (outdir / name).open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("base", type=Path)
    parser.add_argument("outdir", type=Path)
    args = parser.parse_args()
    summary = verify(args.base, args.outdir)
    print(json.dumps({k: v for k, v in summary.items() if k != "sources"}, indent=2))
    raise SystemExit(0 if summary["passed"] else 1)
