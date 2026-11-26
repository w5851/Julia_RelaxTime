"""Compare a kappa/lambda candidate against the baseline it would supersede.

Prints a point-by-point table (absolute and relative change for both recorded
columns) so a numeric commitment change is reviewed per record instead of being
summarised by a single correlation number.
"""

import argparse
import csv
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def load(path: Path) -> dict[tuple[str, str], dict[str, float]]:
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    table = {}
    for row in rows:
        key = (row["record_id"], row["field"])
        if key in table:
            raise SystemExit(f"duplicate record in {path.name}: {key}")
        table[key] = {
            "T": float(row["T_MeV"]), "muB": float(row["muB_MeV"]),
            "expected": float(row["expected"]), "refined": float(row["refined"]),
        }
    return table


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--baseline", type=Path,
                        default=REPO / "tests" / "baselines" / "relaxtime" / "baseline_kappa_lambda_thermal_weighted_v2.csv")
    args = parser.parse_args()

    candidate, baseline = load(args.candidate), load(args.baseline)
    missing = sorted(set(baseline) - set(candidate))
    added = sorted(set(candidate) - set(baseline))
    if missing:
        print(f"  missing records: {missing}")
    if added:
        print(f"  unexpected records: {added}")

    print(f"  {'record':22} {'T':>6} {'muB':>5} {'column':9} {'baseline':>14} {'candidate':>14} {'rel%':>9}")
    worst = 0.0
    for key in sorted(set(candidate) & set(baseline)):
        new, old = candidate[key], baseline[key]
        for column in ("expected", "refined"):
            relative = (new[column] / old[column] - 1) * 100
            worst = max(worst, abs(relative))
            print(f"  {key[0]:22} {new['T']:6.0f} {new['muB']:5.0f} {column:9} "
                  f"{old[column]:14.8f} {new[column]:14.8f} {relative:9.4f}")
    print(f"\n  worst |relative| = {worst:.4f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
