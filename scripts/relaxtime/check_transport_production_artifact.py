"""Check actual profile, point completeness and channel budgets of a runner artifact."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

RATE_KEYS = {
    "uu_to_uu", "ud_to_ud", "us_to_us", "udbar_to_udbar", "dubar_to_dubar",
    "uubar_to_uubar", "uubar_to_ddbar", "usbar_to_usbar", "subar_to_subar",
    "uubar_to_ssbar", "ss_to_ss", "ssbar_to_uubar", "ssbar_to_ssbar",
    "ubardbar_to_ubardbar", "ubarubar_to_ubarubar", "ubarsbar_to_ubarsbar", "sbarsbar_to_sbarsbar",
}
OBSERVABLES = ("tau_u", "tau_d", "tau_s", "tau_ubar", "tau_dbar", "tau_sbar",
               "eta_over_s", "zeta_over_s", "sigma_over_T")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(line for line in stream if line.strip() and not line.startswith("#")))


def key(row: dict[str, str]) -> tuple[float, float, float]:
    return tuple(round(float(row[k]), 10) for k in ("T_MeV", "muB_MeV", "xi"))


def check_artifact(directory: Path, profile: Path, calculation_sha: str, tier: str,
                   source_root: Path | None = None) -> dict:
    manifest = json.loads((directory / "run_manifest.json").read_text(encoding="utf-8"))
    cfg = json.loads((directory / "effective_config.json").read_text(encoding="utf-8"))["options"]
    if manifest["git_commit"] != calculation_sha:
        raise ValueError("artifact calculation commit differs from the requested source")
    sampling = cfg["sampling_profile"]
    if sampling["algorithm"] != "thermal_weighted" or sampling["source_sha256"] != hashlib.sha256(profile.read_bytes()).hexdigest():
        raise ValueError("artifact does not use the requested sampling profile")
    plan = read_rows(directory / "sampling_plan.csv")
    expected = {key(row) for row in plan}
    if not expected or len(expected) != len(plan):
        raise ValueError("empty or duplicate sampling plan")
    if tier == "dry-run":
        if not cfg["dry_run"] or not manifest["summary"]["dry_run"]:
            raise ValueError("dry-run flag does not match the artifact")
        return {"dry_run": True, "planned_points": len(expected)}
    if tier != "production" or cfg["dry_run"]:
        raise ValueError("invalid production tier")
    contract = json.loads((directory / "phase_guided_transport_scan_sampling_contract.json").read_text(encoding="utf-8"))["identity"]
    if contract["profile_sha256"] != sampling["source_sha256"] or contract["profile_id"] != sampling["id"]:
        raise ValueError("result sampling contract differs from its effective configuration")
    if contract["propagator_xi_policy"] != cfg["propagator_xi_policy"] or contract["compute_bulk"] != cfg["compute_bulk"]:
        raise ValueError("result physical contract differs from its effective configuration")
    root = (source_root or Path(__file__).resolve().parents[2]).resolve()
    sources = contract["sources"]
    if not sources or len({s["path"] for s in sources}) != len(sources):
        raise ValueError("missing or duplicate calculation source identity")
    for source in sources:
        path = (root / source["path"]).resolve()
        if not path.is_relative_to(root) or not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != source["sha256"]:
            raise ValueError(f"artifact numerical source differs from the checked source: {source['path']}")
    rows = read_rows(directory / "phase_guided_transport_scan.csv")
    actual = {key(row) for row in rows}
    if actual != expected or len(rows) != len(actual):
        raise ValueError("result point set does not match the plan")
    summary = manifest["summary"]
    if summary["error_count"] or summary["success_count"] + summary["skipped_count"] != len(expected):
        raise ValueError("scan did not complete every planned point")
    unresolved = [row for row in read_rows(directory / "failed_points.csv") if key(row) not in actual]
    if unresolved:
        raise ValueError("unresolved failed points remain")
    diagnostics = read_rows(directory / "phase_guided_transport_scan_sampling.csv")
    for row in rows:
        if row.get("converged", "").strip().lower() != "true":
            raise ValueError(f"unconverged result at {key(row)}")
        if not row["run_id"].strip():
            raise ValueError("result has no run identity")
        for name in OBSERVABLES:
            value = float(row[name])
            if not math.isfinite(value) or value < 0 or (name != "zeta_over_s" and value == 0):
                raise ValueError(f"invalid final observable {name} at {key(row)}")
        selected = [d for d in diagnostics if key(d) == key(row) and d["run_id"] == row["run_id"]]
        if len(selected) != len(RATE_KEYS) or {d["channel"] for d in selected} != RATE_KEYS:
            raise ValueError(f"missing or duplicate channel diagnostics at {key(row)}")
        for d in selected:
            rate, error = float(d["rate"]), float(d["estimated_error"])
            budget = sampling["sampling"]["atol"] + sampling["sampling"]["rtol"] * abs(rate)
            if (d["converged"].lower() != "true" or not math.isfinite(rate) or rate < 0
                    or not math.isfinite(error) or error < 0 or error > budget * (1 + 1e-12)):
                raise ValueError(f"channel budget failed at {key(row)}: {d['channel']}")
    return {"dry_run": False, "completed_points": len(rows), "qualified_channels": len(rows) * len(RATE_KEYS)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--profile", type=Path, required=True)
    parser.add_argument("--calculation-sha", required=True)
    parser.add_argument("--tier", choices=("dry-run", "production"), required=True)
    args = parser.parse_args()
    print(json.dumps(check_artifact(args.directory, args.profile, args.calculation_sha, args.tier), indent=2))


if __name__ == "__main__":
    main()
