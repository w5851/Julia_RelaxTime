#!/usr/bin/env python3
"""Convert WebPlotDigitizer wide tau CSV into validation long-table format."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


SERIES_META = {
    "tau_s_muB800": {"muB_MeV": "800.0", "flavor": "s", "antiparticle": "false"},
    "tau_sbar_muB800": {"muB_MeV": "800.0", "flavor": "s", "antiparticle": "true"},
    "tau_u_muB0": {"muB_MeV": "0.0", "flavor": "u", "antiparticle": "false"},
    "tau_s_muB0": {"muB_MeV": "0.0", "flavor": "s", "antiparticle": "false"},
    "tau_u_muB800": {"muB_MeV": "800.0", "flavor": "u", "antiparticle": "false"},
    "tau_ubar_muB800": {"muB_MeV": "800.0", "flavor": "u", "antiparticle": "true"},
}


def _find_project_root() -> Path:
    script_dir = Path(__file__).resolve().parent
    candidates = [script_dir, script_dir.parent, script_dir.parent.parent]
    candidates.append(Path.cwd())

    for start in candidates:
        current = start
        for _ in range(5):
            if (current / "Project.toml").exists() or (current / ".git").exists():
                return current
            parent = current.parent
            if parent == current:
                break
            current = parent
    return Path.cwd()


PROJECT_ROOT = _find_project_root()
DEFAULT_INPUT = PROJECT_ROOT / "docs" / "reference" / "wpd_datasets (2).csv"
DEFAULT_OUTPUT = PROJECT_ROOT / "tests" / "validation" / "data" / "relaxtime_tau_literature_digitized_longtable_v1.csv"


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def _normalize(input_path: Path, output_path: Path) -> int:
    with input_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle))

    if len(rows) < 2:
        raise ValueError(f"invalid WPD CSV: {input_path}")

    header = rows[0]
    output_path.parent.mkdir(parents=True, exist_ok=True)

    count = 0
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "point_id",
            "series",
            "muB_MeV",
            "flavor",
            "antiparticle",
            "T_MeV",
            "tau_fm",
            "source",
        ])

        for column in range(0, len(header), 2):
            series = header[column].strip()
            if not series:
                continue
            if series not in SERIES_META:
                raise KeyError(f"unsupported series: {series}")

            meta = SERIES_META[series]
            series_index = 0
            for row in rows[2:]:
                if column + 1 >= len(row):
                    continue
                x_raw = row[column].strip()
                y_raw = row[column + 1].strip()
                if not x_raw or not y_raw:
                    continue
                series_index += 1
                count += 1
                writer.writerow([
                    f"{series}_{series_index}",
                    series,
                    meta["muB_MeV"],
                    meta["flavor"],
                    meta["antiparticle"],
                    x_raw,
                    y_raw,
                    "wpd_datasets (2).csv",
                ])

    return count


def main() -> None:
    args = _parse_args()
    count = _normalize(args.input, args.output)
    print(f"Wrote {count} rows to {args.output}")


if __name__ == "__main__":
    main()