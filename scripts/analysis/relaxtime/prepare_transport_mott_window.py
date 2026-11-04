#!/usr/bin/env python3
"""Select an original-grid Mott diagnostic window; execute on quote-remote."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def identity(row):
    return tuple(float(row[key]) for key in ("T_MeV", "muB_MeV", "xi"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--muB", type=float, required=True)
    parser.add_argument("--alpha", type=float, required=True)
    parser.add_argument("--xi-min", type=float, required=True)
    parser.add_argument("--xi-max", type=float, required=True)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise FileExistsError(args.output_dir)
    paths = {name: args.grid_dir / name for name in (
        "phase_guided_transport_scan.csv", "phase_guided_transport_scan_sampling.csv")}
    expected = {
        "phase_guided_transport_scan.csv": "4d812371ddd69678606fb0920813a274395865c1cd94c28f88cd26ef4a5c68c4",
        "phase_guided_transport_scan_sampling.csv": "5990c286aa667ca7c62e750f193984ac44056c2eb3d4c5a90018eaa7b0a9ec66",
    }
    sources = []
    for name, path in paths.items():
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != expected[name]:
            raise ValueError(f"frozen full-grid input differs: {path}")
        sources.append(dict(path=str(path), sha256=digest, bytes=path.stat().st_size))
    window = f"muB{args.muB:g}_alpha{args.alpha:g}_xi{args.xi_min:g}_{args.xi_max:g}"
    selected = []
    for line, row in enumerate(read_csv(paths["phase_guided_transport_scan.csv"]), 2):
        if (float(row["muB_MeV"]) == args.muB and float(row["alpha_T"]) == args.alpha
                and args.xi_min <= float(row["xi"]) <= args.xi_max):
            selected.append(dict(source_csv_line=line, **row))
    selected.sort(key=lambda row: float(row["xi"]))
    if len(selected) < 3 or len({identity(row) for row in selected}) != len(selected):
        raise ValueError("window must contain at least three distinct original points")
    rows = [dict(point_id=f"{args.prefix}{i:03d}", window=window, **row)
            for i, row in enumerate(selected, 1)]
    if any(row["converged"] != "true" for row in rows):
        raise ValueError("unconverged original background")
    keys = {identity(row): row for row in rows}
    rates = [dict(point_id=keys[identity(row)]["point_id"], window=window, **row)
             for row in read_csv(paths["phase_guided_transport_scan_sampling.csv"])
             if identity(row) in keys]
    for key, source_row in keys.items():
        channels = [row for row in rates if identity(row) == key]
        if len(channels) != 17 or len({row["channel"] for row in channels}) != 17:
            raise ValueError(f"incomplete channel coverage: {key}")
        if any(row["run_id"] != source_row["run_id"] or row["converged"] != "true" for row in channels):
            raise ValueError(f"channel provenance/acceptance mismatch: {key}")
    changes = []
    for left, right in zip(rows, rows[1:]):
        item = dict(left_id=left["point_id"], right_id=right["point_id"],
                    xi_left=float(left["xi"]), xi_right=float(right["xi"]),
                    phase_left=left["phase_curr"], phase_right=right["phase_curr"])
        for name in ("tau_u", "tau_s", "tau_ubar", "tau_sbar", "eta_over_s", "zeta_over_s", "sigma_over_T"):
            item[name + "_change_percent"] = 100 * (float(right[name]) / float(left[name]) - 1)
        changes.append(item)
    args.output_dir.mkdir(parents=True)
    write_csv(args.output_dir / "source_rows.csv", rows)
    write_csv(args.output_dir / "source_rates.csv", rates)
    write_csv(args.output_dir / "adjacent_changes.csv", changes)
    summary = dict(schema="transport_mott_window_inputs_v1", window=window,
                   computation_host="quote-remote", sources=sources, rows=len(rows), channels=len(rates),
                   source_values_unchanged=True, solver_called=False,
                   points=[{key: row[key] for key in ("point_id", "xi", "T_MeV", "phase_curr", "tau_ubar")}
                           for row in rows],
                   largest_tau_ubar_steps=sorted(changes, key=lambda row: abs(row["tau_ubar_change_percent"]), reverse=True)[:6])
    (args.output_dir / "input_manifest.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
