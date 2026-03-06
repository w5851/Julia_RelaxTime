#!/usr/bin/env python3
"""Convert WebPlotDigitizer wide sigma CSVs into a validation long-table."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


SERIES_META = {
    "T200MeV": {"muB_MeV": "0.0", "process": "ud_to_ud", "T_MeV": "200.0"},
    "T160MeV": {"muB_MeV": "0.0", "process": "ud_to_ud", "T_MeV": "160.0"},
    "T250MeV": {"muB_MeV": "0.0", "process": "ud_to_ud", "T_MeV": "250.0"},
    "T160MeV_ustous": {"muB_MeV": "0.0", "process": "us_to_us", "T_MeV": "160.0"},
    "T250MeV_ustous": {"muB_MeV": "0.0", "process": "us_to_us", "T_MeV": "250.0"},
    "T210MeV_udbartoudbar": {"muB_MeV": "0.0", "process": "udbar_to_udbar", "T_MeV": "210.0"},
    "T250MeV_udbartoudbar": {"muB_MeV": "0.0", "process": "udbar_to_udbar", "T_MeV": "250.0"},
    "T250MeV_usbartousbar": {"muB_MeV": "0.0", "process": "usbar_to_usbar", "T_MeV": "250.0"},
    "T210MeV_usbartousbar": {"muB_MeV": "0.0", "process": "usbar_to_usbar", "T_MeV": "210.0"},
    "T100_udtoud": {"muB_MeV": "800.0", "process": "ud_to_ud", "T_MeV": "100.0"},
    "T180_udtoud": {"muB_MeV": "800.0", "process": "ud_to_ud", "T_MeV": "180.0"},
    "T180_ustous": {"muB_MeV": "800.0", "process": "us_to_us", "T_MeV": "180.0"},
    "T100_ustous": {"muB_MeV": "800.0", "process": "us_to_us", "T_MeV": "100.0"},
    "T180_udbartoudbar": {"muB_MeV": "800.0", "process": "udbar_to_udbar", "T_MeV": "180.0"},
    "T100_udbartoudbar": {"muB_MeV": "800.0", "process": "udbar_to_udbar", "T_MeV": "100.0"},
    "T150_udbartoudbar": {"muB_MeV": "800.0", "process": "udbar_to_udbar", "T_MeV": "150.0"},
    "T100_usbartousbar": {"muB_MeV": "800.0", "process": "usbar_to_usbar", "T_MeV": "100.0"},
    "T180_usbartousbar": {"muB_MeV": "800.0", "process": "usbar_to_usbar", "T_MeV": "180.0"},
    "T150_usbartousbar": {"muB_MeV": "800.0", "process": "usbar_to_usbar", "T_MeV": "150.0"},
}


def _find_project_root() -> Path:
    script_dir = Path(__file__).resolve().parent
    candidates = [script_dir, script_dir.parent, script_dir.parent.parent, Path.cwd()]

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
DEFAULT_INPUT_MUB0 = PROJECT_ROOT / "docs" / "reference" / "sigma(mb)-sqrt(s)(GeV)_muB0.csv"
DEFAULT_INPUT_MUB800 = PROJECT_ROOT / "docs" / "reference" / "sigma(mb)-sqrt(s)(GeV)_muB800.csv"
DEFAULT_OUTPUT = PROJECT_ROOT / "tests" / "validation" / "data" / "relaxtime_sigma_literature_digitized_longtable_v1.csv"

INPUT_SERIES = {
    DEFAULT_INPUT_MUB0.name: [
        "T200MeV",
        "T160MeV",
        "T250MeV",
        "T160MeV_ustous",
        "T250MeV_ustous",
        "T210MeV_udbartoudbar",
        "T250MeV_udbartoudbar",
        "T250MeV_usbartousbar",
        "T210MeV_usbartousbar",
    ],
    DEFAULT_INPUT_MUB800.name: [
        "T100_udtoud",
        "T180_udtoud",
        "T180_ustous",
        "T100_ustous",
        "T180_udbartoudbar",
        "T100_udbartoudbar",
        "T150_udbartoudbar",
        "T100_usbartousbar",
        "T180_usbartousbar",
        "T150_usbartousbar",
    ],
}


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-mub0", type=Path, default=DEFAULT_INPUT_MUB0)
    parser.add_argument("--input-mub800", type=Path, default=DEFAULT_INPUT_MUB800)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def _normalize_one(input_path: Path, keep_series: list[str], writer: csv.writer, count: int) -> int:
    with input_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle))

    if len(rows) < 2:
        raise ValueError(f"invalid WPD CSV: {input_path}")

    header = rows[0]
    allowed = set(keep_series)

    for column in range(0, len(header), 2):
        series = header[column].strip()
        if not series or series not in allowed:
            continue

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
                meta["process"],
                meta["T_MeV"],
                meta["muB_MeV"],
                x_raw,
                y_raw,
                input_path.name,
            ])

    return count


def main() -> None:
    args = _parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    count = 0
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "point_id",
            "series",
            "process",
            "T_MeV",
            "muB_MeV",
            "sqrt_s_GeV",
            "sigma_mb",
            "source",
        ])

        count = _normalize_one(args.input_mub0, INPUT_SERIES[args.input_mub0.name], writer, count)
        count = _normalize_one(args.input_mub800, INPUT_SERIES[args.input_mub800.name], writer, count)

    print(f"Wrote {count} rows to {args.output}")


if __name__ == "__main__":
    main()