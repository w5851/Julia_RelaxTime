"""Promote a reviewed twelve-point transport candidate into a versioned baseline.

The candidate comes from
`scripts/dev/export_phase_guided_twelve_point_baseline.jl`; this step records it
with the provenance the regression reads back (`baseline_sha256`,
`profile_sha256`) plus what produced it.

The `qualified_candidate` / `qualified_reference` columns are filled from the same
runs, so they are **not** a fresh independent qualification: the historical
external qualification directory recorded by v3 was not re-run.  The metadata says
so explicitly instead of implying a requalification that did not happen.
"""

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
BASELINE_DIR = REPO / "tests" / "baselines" / "relaxtime"
DEFAULT_PROFILE = REPO / "config" / "models" / "pnjl" / "transport_sampling" / "thermal_weighted_v1.toml"
REFERENCE_PROFILE = REPO / "config" / "models" / "pnjl" / "transport_sampling" / "thermal_weighted_reference_v1.toml"

SOURCE_FILES = (
    "src/relaxtime/AverageScatteringRate.jl",
    "src/relaxtime/ThermalWeightedCache.jl",
    "src/relaxtime/ThermalWeightedSampling.jl",
    "src/relaxtime/RelaxationTime.jl",
    "src/models/workflow_apps/TransportWorkflow.jl",
    "src/models/workflow_apps/TransportSamplingProfiles.jl",
    "config/physics/default.toml",
    "Manifest.toml",
)

OBSERVABLES = {"tau_u", "tau_d", "tau_s", "tau_ubar", "tau_dbar", "tau_sbar",
               "eta_over_s", "zeta_over_s", "sigma_over_T"}
EXPECTED_HEADER = ["point_id", "T_MeV", "muB_MeV", "xi", "observable",
                   "candidate", "reference", "qualified_candidate"]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("version", help="target version, for example v4")
    parser.add_argument("--reason", required=True, help="why the numeric commitment changed")
    parser.add_argument("--supersedes", default="baseline_thermal_weighted_transport_v3.csv")
    args = parser.parse_args(argv)

    candidate: Path = args.candidate
    if not candidate.is_file():
        parser.error(f"candidate not found: {candidate}")
    lines = candidate.read_text(encoding="utf-8").splitlines()
    if lines[0].split(",") != EXPECTED_HEADER:
        parser.error(f"candidate header mismatch: {lines[0].split(',')}")

    rows = [line.split(",") for line in lines[1:] if line.strip()]
    if len(rows) != 108:
        parser.error(f"expected 108 rows (12 points x 9 observables), found {len(rows)}")
    if len({(row[0], row[4]) for row in rows}) != len(rows):
        parser.error("duplicate (point_id, observable) row")
    if {row[4] for row in rows} != OBSERVABLES:
        parser.error(f"unexpected observable set: {sorted({row[4] for row in rows})}")
    if len({row[0] for row in rows}) != 12:
        parser.error("expected 12 distinct points")
    for row in rows:
        for value in row[5:8]:
            if not value or value.lower() in {"nan", "inf", "-inf"}:
                parser.error(f"non-finite value in {row[0]}/{row[4]}")

    target_csv = BASELINE_DIR / f"baseline_thermal_weighted_transport_{args.version}.csv"
    target_json = target_csv.with_suffix(".json")
    for path in (target_csv, target_json):
        if path.exists():
            parser.error(f"refusing to overwrite existing baseline: {path}")

    previous = json.loads((BASELINE_DIR / "baseline_thermal_weighted_transport_v3.json").read_text(encoding="utf-8"))
    commit = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    julia_version = subprocess.run(["julia", "--version"], capture_output=True, text=True).stdout.strip().replace("julia version ", "")

    target_csv.write_text(candidate.read_text(encoding="utf-8"), encoding="utf-8", newline="\n")
    record = {
        "schema_version": previous.get("schema_version", 1),
        "purpose": previous.get("purpose"),
        "profile_sha256": sha256(DEFAULT_PROFILE),
        "reference_profile": str(REFERENCE_PROFILE.relative_to(REPO)).replace("\\", "/"),
        "reference_profile_sha256": sha256(REFERENCE_PROFILE),
        "qualification_status": "numeric_commitment_refresh_without_requalification",
        "qualification_note": (
            "v3 was qualified against an external joint-validation directory. That round was not "
            "repeated here: `qualified_candidate`/`qualified_reference` are filled from the same "
            "candidate/reference runs recorded in this refresh, so they must not be read as a new "
            "independent qualification."
        ),
        "numeric_commitment_change": args.reason,
        "supersedes": args.supersedes,
        "regenerated_by": "scripts/dev/export_phase_guided_twelve_point_baseline.jl",
        "host": "quote-remote",
        "julia_version": julia_version,
        "threads": 1,
        "source_commit": commit,
        "source_files": {rel: sha256(REPO / rel) for rel in SOURCE_FILES},
        "configurations": {
            "candidate": {"profile": "thermal_weighted_v1", "profile_sha256": sha256(DEFAULT_PROFILE)},
            "reference": {"profile": "thermal_weighted_reference_v1",
                          "profile_sha256": sha256(REFERENCE_PROFILE),
                          "scope": "tighter quadrature nodes only; default sampling tolerance"},
        },
        "points": 12,
        "comparisons": len(rows),
        "regression_rtol": previous.get("regression_rtol", 1e-6),
        "regression_atol": previous.get("regression_atol", 1e-8),
        "independent_relative_target": previous.get("independent_relative_target", 0.01),
        "baseline_semantics": previous.get("baseline_semantics", "production_entrypoint"),
        "baseline_sha256": sha256(target_csv),
    }
    target_json.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print(f"  promoted {target_csv.relative_to(REPO)}")
    print(f"  metadata {target_json.relative_to(REPO)}  sha256={record['baseline_sha256'][:16]}")
    print(f"  profile  {record['profile_sha256'][:16]} (must match the loaded default profile)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
