"""Compare the twelve-point transport candidate against the baseline it would supersede.

Reports the per-observable drift and the candidate/reference agreement band, and
flags any point where the pair leaves the recorded 1% independent-reference
tolerance.
"""

import argparse
import csv
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
LIMIT = 0.01


def load(path: Path) -> dict[tuple[str, str], dict[str, float]]:
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    table = {}
    for row in rows:
        key = (row["point_id"], row["observable"])
        if key in table:
            raise SystemExit(f"duplicate row: {key}")
        table[key] = {name: float(row[name]) for name in ("candidate", "reference", "qualified_candidate")}
    return table


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--baseline", type=Path,
                        default=REPO / "tests" / "baselines" / "relaxtime" / "baseline_thermal_weighted_transport_v3.csv")
    args = parser.parse_args()

    new, old = load(args.candidate), load(args.baseline)
    missing = sorted(set(old) - set(new))
    added = sorted(set(new) - set(old))
    if missing:
        print(f"  missing rows ({len(missing)}): {missing[:4]}")
    if added:
        print(f"  unexpected rows ({len(added)}): {added[:4]}")

    worst_pair = (0.0, "")
    worst_drift = (0.0, "")
    over_limit = []
    print(f"  {'point_id':34} {'observable':14} {'v3':>13} {'v4':>13} {'rel%':>9} {'pair%':>8}")
    for key in sorted(set(new) & set(old)):
        value_new, value_old = new[key], old[key]
        drift = (value_new["candidate"] / value_old["candidate"] - 1) * 100
        pair = abs(value_new["candidate"] / value_new["reference"] - 1) * 100
        if abs(drift) > abs(worst_drift[0]):
            worst_drift = (drift, f"{key[0]}/{key[1]}")
        if pair > worst_pair[0]:
            worst_pair = (pair, f"{key[0]}/{key[1]}")
        if pair / 100 > LIMIT:
            over_limit.append((key, pair))
        if abs(drift) >= 0.1:
            print(f"  {key[0]:34} {key[1]:14} {value_old['candidate']:13.7f} {value_new['candidate']:13.7f} "
                  f"{drift:9.4f} {pair:8.4f}")

    print(f"\n  worst |drift|      = {worst_drift[0]:.4f}%  at {worst_drift[1]}")
    print(f"  worst |pair|       = {worst_pair[0]:.4f}%  at {worst_pair[1]}  (limit {LIMIT * 100:.0f}%)")
    print(f"  rows outside band = {len(over_limit)}")
    for key, pair in over_limit[:5]:
        print(f"    {key[0]} {key[1]} pair={pair:.4f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
