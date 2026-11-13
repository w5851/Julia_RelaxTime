"""Condense the PR3 numerical runs, preserving failures and comparison inputs."""

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
BASE_COMMIT = "abe6dd2b83587c587b0e87deee82aa4526e8608b"
SOURCE_FILES = (
    "src/models/pnjl_physics/PNJLIntegrals.jl",
    "src/models/derivatives/PNJLChiBTaylorDiff.jl",
    "src/models/derivatives/TaylorDiffForwardDiffCompat.jl",
    "src/models/derivatives/ThermoDerivatives.jl",
    "src/models/derivatives/ConservedChargeSusceptibilities.jl",
    "scripts/analysis/pnjl/low_temperature_consumers.jl",
)


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def errors(value, path=""):
    if isinstance(value, dict):
        if value.get("status") == "error":
            yield {"path": path, "error": value.get("error")}
        for key, item in value.items():
            yield from errors(item, f"{path}/{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            yield from errors(item, f"{path}/{index}")


def maximum_scaled_difference(value, reference):
    if isinstance(reference, dict):
        return max((maximum_scaled_difference(value[k], v) for k, v in reference.items()
                    if k not in {"elapsed_seconds", "jacobian_condition", "residual_norm"}), default=0)
    if isinstance(reference, list):
        return max((maximum_scaled_difference(a, b) for a, b in zip(value, reference)), default=0)
    if isinstance(reference, (float, int)) and not isinstance(reference, bool):
        if not isinstance(value, (float, int)):
            return math.inf
        return abs(value - reference) / max(abs(reference), 1e-9)
    return 0 if value == reference else math.inf


def coordinates(row):
    return {k: row[k] for k in ("T_MeV", "mu_q_inv_fm", "xi")}


def consumer(value):
    bulk = value["bulk_result"]["value"]
    mass = value["mass_result"]["value"]
    thermo = value["thermo_result"]["value"]
    return {
        "state": value["state"], "masses": value["masses"],
        "pressure_mu_0_to_4": value["pressure_mu"]["value"],
        "density": value["density"], "dM_dT": mass["dM_dT"],
        "dM_dmu_q": mass["dM_dmu"],
        "chi2_B": value["chi_result"]["value"], "chi4_B": value["chi4_result"]["value"],
        "thermo": {k: thermo[k] for k in ("dP_dT", "dP_dmu", "dEpsilon_dT",
                   "dEpsilon_dmu", "dn_dT", "dn_dmu")},
        "bulk": {k: v for k, v in bulk.items() if k in
                 ("v_n_sq", "dμB_dT_sigma", "ds_dT", "ds_dμB", "dn_dT", "dn_dμB",
                  "c_p", "s", "n_B")},
        "residual_norm": value["residual_norm"],
        "jacobian_condition": value["jacobian_condition"],
    }


def summarize_fixed(row):
    methods = row["methods"]
    adaptive = [m["result"]["value"] for m in methods if m["method"] == "adaptive"]
    high = [m["result"]["value"] for m in methods if m["method"] == "bigfloat_window_oracle"]
    windows = [m["result"]["value"] for m in methods if m["method"] == "energy_window_oracle"]
    reference = high[0] if high else windows[-1]["chemical"] if windows else row["zero_limit"]
    compared = [(a, b) for a, b in zip(adaptive[-1]["chemical"], reference) if b is not None]
    return {**coordinates(row), "at_fermi_threshold": row["at_fermi_threshold"],
            "adaptive": adaptive[-1], "zero_limit": row["zero_limit"],
            "chemical_reference": reference,
            "reference_method": "192-bit q window" if high else "96-node q window" if windows else "analytic sea",
            "max_scaled_tolerance_spread": max(maximum_scaled_difference(v, adaptive[-1]) for v in adaptive),
            "max_scaled_chemical_reference_difference": max((abs(a-b)/max(abs(b), 1e-9) for a,b in compared), default=0),
            "window_64_to_96_max_scaled_difference": maximum_scaled_difference(windows[0], windows[-1]) if windows else None}


def summarize_consumers(row):
    values = [m["result"]["value"] for m in row["methods"]]
    summaries = [consumer(v) for v in values]
    final = summaries[-1]
    return {**coordinates(row), **final,
            "max_scaled_tolerance_spread": max(maximum_scaled_difference(v, final) for v in summaries),
            "max_residual_across_tolerances": max(v["residual_norm"] for v in values),
            "zero_fixed_state": values[-1]["zero_fixed_state"],
            "response_interpretation": "vanishing_density_ratios_unqualified" if abs(final["bulk"]["n_B"]) < 1e-8 else "finite_density_branch",
            "oracle": row["oracle"]["value"] if row["oracle"] else None}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    names = ("fixed_final.json", "consumers_full.json", "audit_final.json", "fixed_initial.json")
    reports = {name: read(args.input_dir / name) for name in names}
    failures = [dict(source=name, **error) for name, report in reports.items()
                if name != "fixed_initial.json" for error in errors(report)]
    if failures:
        raise RuntimeError(json.dumps(failures, ensure_ascii=False))
    fixed = [summarize_fixed(r) for r in reports["fixed_final.json"]["rows"]]
    consumers = [summarize_consumers(r) for r in reports["consumers_full.json"]["rows"]]
    audit = []
    for row in reports["audit_final.json"]["rows"]:
        value = row["result"]["value"]
        if row["check"] == "total_derivatives":
            value["original"] = consumer(value["original"])
            value["strict"] = consumer(value["strict"])
        audit.append({"check": row["check"], **value})
    counterexample = next(r for r in reports["fixed_initial.json"]["rows"]
                          if r["T_MeV"] == 0.01 and r["mu_q_inv_fm"] == 1.8 and r["xi"] == 0)
    counterexample["methods"] = [m for m in counterexample["methods"] if m["method"] == "adaptive"]
    summary = {
        "schema": "pnjl_low_temperature_pr3_summary_v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "execution": {"host": "quote-remote", "julia": "1.12.5", "threads": 1,
                      "openblas_threads": 1, "model_profile": "PNJL default",
                      "source_base_commit": BASE_COMMIT,
                      "reproduction_source_sha256_lf": {
                          p: hashlib.sha256((ROOT/p).read_bytes().replace(b"\r\n", b"\n")).hexdigest()
                          for p in SOURCE_FILES}},
        "controls": {"rtol": [1e-6, 1e-8, 1e-10], "atol_over_rtol": 1e-3,
                     "maxevals": 200000, "base_residual_tol": 1e-9,
                     "strict_base_residual_tol": 1e-13,
                     "vanishing_density_reporting_threshold_inv_fm3": 1e-8,
                     "fixed_masses_inv_fm": [1.5, 1.5, 2.5], "fixed_phi": 0.1, "fixed_phibar": 0.1,
                     "audit_oracle_phi": 0.12, "audit_oracle_phibar": 0.28,
                     "oracle_bits": 192, "oracle_nodes": [96, 128], "oracle_tail_in_T": [56, 72]},
        "units": {"T_MeV": "MeV", "mu_q_inv_fm": "fm^-1",
                  "state": "[phi_u,phi_d,phi_s,Phi,PhiBar]; condensates fm^-3, loops dimensionless",
                  "pressure_mu_0_to_4": "d^n P / d mu_q^n, n=0:4, fm^(n-4)",
                  "chemical": "same as pressure_mu_0_to_4; fixed mean fields",
                  "thermal": "[dP/dT, d2P/dT2], fixed mean fields, fm^-3 and fm^-2",
                  "temperature_derivatives": "with respect to T in fm^-1",
                  "nonfinite": "NaN and +/-Inf are strings",
                  "scaled_difference": "abs(candidate-reference)/max(abs(reference),1e-9), componentwise; not a per-component relative-error guarantee"},
        "counts": {"fixed_points": len(fixed), "self_consistent_points": len(consumers),
                   "audit_checks": len(audit), "capture_failures": len(failures)},
        "legacy_counterexample": {"source_commit": BASE_COMMIT, **counterexample},
        "fixed": fixed, "consumers": consumers, "audit": audit,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(summary["counts"]))
    print("fixed max reference difference", max(r["max_scaled_chemical_reference_difference"] for r in fixed))
    print("consumer max tolerance spread", max(r["max_scaled_tolerance_spread"] for r in consumers))


if __name__ == "__main__":
    main()
