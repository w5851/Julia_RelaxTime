#!/usr/bin/env python3
"""Freeze v15 input rows for the other-Mott diagnostic (quote-remote)."""
from pathlib import Path
from collections import defaultdict
import argparse
import csv
import hashlib
import json
import shutil


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("grid", type=Path)
    parser.add_argument("profiles", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    rows, sources = [], []
    for mode in ("mode_a", "mode_b"):
        path = args.grid / mode / "phase_guided_transport_grid.csv"
        sources.append(dict(path=str(path), bytes=path.stat().st_size, sha256=sha(path)))
        with path.open(encoding="utf-8-sig", newline="") as stream:
            groups = defaultdict(list)
            for row in csv.DictReader(stream):
                groups[row["scan_group"]].append(row)
        for group, curve in sorted(groups.items()):
            segment = 0
            for row in sorted(curve, key=lambda r: float(r["xi"])):
                window = f"{mode}:{group}"
                rows.append(dict(window=window, continuation_segment=f"{window}:segment{segment}",
                                 mode_key=mode, **row))
                first_order_left = (
                    mode == "mode_a" and row["muB_MeV"] == "900.0"
                    and row["alpha_T"] == "1.0" and float(row["xi"]) == -0.003
                ) or (
                    mode == "mode_b" and float(row["T_MeV"]) == 120
                    and float(row["muB_MeV"]) == 900 and float(row["xi"]) == -0.13
                )
                segment += int(first_order_left)
    assert len(rows) == len({r["point_id"] for r in rows}) == 1867
    assert len({r["window"] for r in rows}) == 18
    assert len({r["continuation_segment"] for r in rows}) == 20
    with (args.output / "source_rows.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    profiles = {}
    for name in ("base", "fine", "tighter"):
        path = args.profiles / f"{name}.toml"
        shutil.copyfile(path, args.output / path.name)
        profiles[name] = dict(source=str(path), sha256=sha(path))
    manifest = dict(schema="transport_other_mott_inputs_v1", host="quote-remote",
                    calculation_sha="2e5c411c716d5cbc524c3fc6cfb335fe05b39444",
                    sources=sources, profiles=profiles, points=len(rows), curves=18, segments=20,
                    source_rows_sha256=sha(args.output / "source_rows.csv"),
                    preparation_script_sha256=sha(Path(__file__)),
                    production_data_changed=False, pdf_deferred=True)
    (args.output / "input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest))


if __name__ == "__main__":
    main()
