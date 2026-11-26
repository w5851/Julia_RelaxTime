"""Promote a reviewed kappa/lambda candidate into a versioned regression baseline.

The candidate itself comes from
`scripts/dev/export_kappa_lambda_thermal_weighted_baseline.jl`; this step only
records it together with the provenance the regression reads back:
`baseline_sha256` (so a hand edit is caught), the source-file fingerprints and
the configuration that produced both columns.  It refuses to overwrite an
existing version, so the previous baseline stays available for comparison.
"""

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
BASELINE_DIR = REPO / "tests" / "baselines" / "relaxtime"

SOURCE_FILES = (
    "src/relaxtime/AverageScatteringRate.jl",
    "src/relaxtime/ThermalWeightedCache.jl",
    "src/relaxtime/ThermalWeightedSampling.jl",
    "src/relaxtime/RelaxationTime.jl",
    "src/models/workflow_apps/TransportWorkflow.jl",
    "config/physics/default.toml",
    "Manifest.toml",
)

EXPECTED_HEADER = ["record_id", "T_MeV", "muB_MeV", "field", "expected", "refined"]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate", type=Path, help="reviewed candidate CSV")
    parser.add_argument("version", help="target version, for example v3")
    parser.add_argument("--reason", required=True, help="why the numeric commitment changed")
    args = parser.parse_args(argv)

    candidate: Path = args.candidate
    if not candidate.is_file():
        parser.error(f"candidate not found: {candidate}")
    header = candidate.read_text(encoding="utf-8").splitlines()[0].split(",")
    if header != EXPECTED_HEADER:
        parser.error(f"candidate header mismatch: {header}")

    target_csv = BASELINE_DIR / f"baseline_kappa_lambda_thermal_weighted_{args.version}.csv"
    target_json = target_csv.with_suffix(".json")
    for path in (target_csv, target_json):
        if path.exists():
            parser.error(f"refusing to overwrite existing baseline: {path}")

    rows = [line.split(",") for line in candidate.read_text(encoding="utf-8").splitlines()[1:] if line.strip()]
    if len(rows) != 14:
        parser.error(f"expected 14 records (4 points x 3 kappa + 2 lambda), found {len(rows)}")
    record_ids = [row[0] for row in rows]
    if len(set(record_ids)) != len(record_ids):
        parser.error("duplicate record_id in candidate")
    for row in rows:
        for value in row[4:6]:
            if not value or value.lower() in {"nan", "inf", "-inf"}:
                parser.error(f"non-finite value in record {row[0]}")

    previous = json.loads((BASELINE_DIR / "baseline_kappa_lambda_thermal_weighted_v2.json").read_text(encoding="utf-8"))
    commit = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    julia_version = subprocess.run(["julia", "--version"], capture_output=True, text=True).stdout.strip().replace("julia version ", "")

    target_csv.write_text(candidate.read_text(encoding="utf-8"), encoding="utf-8", newline="\n")
    record = {
        "schema_version": previous.get("schema_version", 1),
        "purpose": previous.get("purpose"),
        "date": subprocess.run(["git", "-C", str(REPO), "log", "-1", "--format=%ad", "--date=short"],
                               capture_output=True, text=True).stdout.strip(),
        "julia_version": julia_version or previous.get("julia_version"),
        "host": "quote-remote",
        "threads": 8,
        "source_commit": commit[:40],
        "numerical_source_commit": previous.get("numerical_source_commit"),
        "legacy_sampling_source_commit": previous.get("legacy_sampling_source_commit"),
        "legacy_target": previous.get("legacy_target"),
        "legacy_target_sha256": previous.get("legacy_target_sha256"),
        "baseline_sha256": sha256(target_csv),
        "source_files": {rel: sha256(REPO / rel) for rel in SOURCE_FILES},
        "configurations": previous.get("configurations"),
        "scope": previous.get("scope"),
        "regression_rtol": previous.get("regression_rtol", 1e-6),
        "regression_atol": previous.get("regression_atol", 1e-8),
        "refinement_relative_limit": previous.get("refinement_relative_limit", 1e-3),
        "supersedes": f"baseline_kappa_lambda_thermal_weighted_{previous.get('schema_version', 1) and 'v2'}.csv",
        "numeric_commitment_change": args.reason,
        "reproduction": (
            "Run scripts/dev/export_kappa_lambda_thermal_weighted_baseline.jl to a fresh path on the "
            "recorded julia_version and candidate/refined configurations, review the point-by-point "
            "difference against the superseded baseline, then promote with "
            "scripts/dev/promote_kappa_lambda_thermal_weighted_baseline.py."
        ),
    }
    target_json.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
    print(f"  promoted {target_csv.relative_to(REPO)}")
    print(f"  metadata {target_json.relative_to(REPO)}  sha256={record['baseline_sha256'][:16]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
