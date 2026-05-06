#!/usr/bin/env python3
"""
Normalize a WebPlotDigitizer multi-dataset CSV export into project literature
target CSV/meta files for meson-density validation.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class DatasetSpec:
    raw_label: str
    slug: str
    quantity: str
    mu_pi_mev: str
    curve_label: str
    anomalous_mode: str
    notes: str


DATASETS: tuple[DatasetSpec, ...] = (
    DatasetSpec(
        raw_label="mu_100_K/pi+with_anomalous",
        slug="blaschke2019col_kplus_piplus_mu_pi_100_fig4_right_with_anomalous",
        quantity="K^+ / pi^+",
        mu_pi_mev="100",
        curve_label="right-panel black thin curve with anomalous-mode contribution",
        anomalous_mode="included",
        notes="This is one boundary of the K^+/pi^+ band-like result layer.",
    ),
    DatasetSpec(
        raw_label="mu_100_K/pi+no_anomalous",
        slug="blaschke2019col_kplus_piplus_mu_pi_100_fig4_right_no_anomalous",
        quantity="K^+ / pi^+",
        mu_pi_mev="100",
        curve_label="right-panel black thin curve without anomalous-mode contribution",
        anomalous_mode="excluded",
        notes="This is the lower K^+/pi^+ reference curve under the same mu_pi.",
    ),
    DatasetSpec(
        raw_label="mu_100_K/pi-",
        slug="blaschke2019col_kminus_piminus_mu_pi_100_fig4_right",
        quantity="K^- / pi^-",
        mu_pi_mev="100",
        curve_label="right-panel red thin dotted curve",
        anomalous_mode="not_applicable",
        notes="Recommended first single-curve charged/freeze-out validation target.",
    ),
    DatasetSpec(
        raw_label="mu_134p5_K/pi-",
        slug="blaschke2019col_kminus_piminus_mu_pi_134p5_fig4_right",
        quantity="K^- / pi^-",
        mu_pi_mev="134.5",
        curve_label="right-panel red thick dashed curve",
        anomalous_mode="not_applicable",
        notes="Second K^-/pi^- curve; close to the mu_pi=100 line in the source figure.",
    ),
)


def _parse_float(value: str) -> float | None:
    text = value.strip()
    if not text:
        return None
    return float(text)


def _clamp_small_negative_to_zero(value: float) -> tuple[float, bool]:
    if value < 0.0:
        return 0.0, True
    return value, False


def _longest_increasing_prefix(points: Iterable[tuple[float, float, int]]) -> tuple[list[tuple[float, float]], int]:
    cleaned: list[tuple[float, float]] = []
    dropped = 0
    last_x: float | None = None
    for x_value, y_value, _row_index in points:
        if last_x is not None and x_value <= last_x:
            dropped += 1
            continue
        cleaned.append((x_value, y_value))
        last_x = x_value
    return cleaned, dropped


def normalize_wpd_csv(source_csv: Path, output_dir: Path) -> list[dict[str, object]]:
    with source_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle))

    if len(rows) < 2:
        raise ValueError("WPD CSV must contain at least dataset and axis header rows.")

    dataset_row = rows[0]
    axis_row = rows[1]
    if len(dataset_row) % 2 != 0:
        raise ValueError("Expected paired X/Y columns in WPD CSV.")

    output_dir.mkdir(parents=True, exist_ok=True)
    reports: list[dict[str, object]] = []

    label_to_spec = {spec.raw_label: spec for spec in DATASETS}
    for column in range(0, len(dataset_row), 2):
        raw_label = dataset_row[column].strip()
        if not raw_label:
            continue
        x_axis = axis_row[column].strip()
        y_axis = axis_row[column + 1].strip()
        if x_axis != "X" or y_axis != "Y":
            raise ValueError(f"Unexpected axis labels for {raw_label}: {x_axis}, {y_axis}")
        spec = label_to_spec.get(raw_label)
        if spec is None:
            raise KeyError(f"Unmapped dataset label: {raw_label}")

        raw_points: list[tuple[float, float, int]] = []
        skipped_blank = 0
        clamped_negative = 0
        for row_index, row in enumerate(rows[2:], start=3):
            if column >= len(row):
                skipped_blank += 1
                continue
            raw_x = row[column].strip() if column < len(row) else ""
            raw_y = row[column + 1].strip() if column + 1 < len(row) else ""
            if not raw_x or not raw_y:
                skipped_blank += 1
                continue
            x_value = _parse_float(raw_x)
            y_value = _parse_float(raw_y)
            if x_value is None or y_value is None:
                skipped_blank += 1
                continue
            y_value, was_clamped = _clamp_small_negative_to_zero(y_value)
            if was_clamped:
                clamped_negative += 1
            raw_points.append((x_value, y_value, row_index))

        cleaned_points, dropped_nonmonotonic = _longest_increasing_prefix(raw_points)
        if not cleaned_points:
            raise ValueError(f"No valid points extracted for dataset {raw_label}")

        csv_path = output_dir / f"{spec.slug}.csv"
        with csv_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["x_value", "y_value"])
            writer.writerows(cleaned_points)

        meta_path = output_dir / f"{spec.slug}.meta.md"
        meta_path.write_text(
            "\n".join(
                [
                    f"# {spec.slug}",
                    "",
                    f"- source_paper: `Blaschke:2019col`",
                    f"- figure: `Figure 4`",
                    f"- panel: `right`",
                    f"- raw_dataset_label: `{spec.raw_label}`",
                    f"- physical_quantity: `{spec.quantity}`",
                    f"- x_axis: `sqrt(s_NN) [GeV]`",
                    f"- y_axis: `ratio`",
                    f"- mu_pi_MeV: `{spec.mu_pi_mev}`",
                    f"- curve_label: `{spec.curve_label}`",
                    f"- anomalous_mode: `{spec.anomalous_mode}`",
                    f"- path_context: `freeze-out / critical-line phenomenology; exact stitched path not yet reconstructed parametrically`",
                    f"- extraction_method: `manual digitization via WebPlotDigitizer, then normalized by project helper`",
                    f"- raw_points: `{len(raw_points)}`",
                    f"- kept_points: `{len(cleaned_points)}`",
                    f"- dropped_nonmonotonic_tail_points: `{dropped_nonmonotonic}`",
                    f"- clamped_negative_points_to_zero: `{clamped_negative}`",
                    f"- skipped_blank_or_incomplete_rows: `{skipped_blank}`",
                    f"- notes: {spec.notes}",
                ]
            )
            + "\n",
            encoding="utf-8",
        )

        reports.append(
            {
                "slug": spec.slug,
                "raw_label": spec.raw_label,
                "raw_points": len(raw_points),
                "kept_points": len(cleaned_points),
                "dropped_nonmonotonic": dropped_nonmonotonic,
                "clamped_negative_to_zero": clamped_negative,
                "skipped_blank": skipped_blank,
                "csv_path": str(csv_path),
                "meta_path": str(meta_path),
            }
        )

    manifest_path = output_dir / "manifest_wpd_import_report.csv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "slug",
                "raw_label",
                "raw_points",
                "kept_points",
                "dropped_nonmonotonic",
                "clamped_negative_to_zero",
                "skipped_blank",
                "csv_path",
                "meta_path",
            ],
        )
        writer.writeheader()
        writer.writerows(reports)

    return reports


def main() -> None:
    repo_root = Path(__file__).resolve().parents[3]
    source_csv = Path(r"D:\Desktop\wpd_datasets.csv")
    output_dir = repo_root / "data" / "outputs" / "results" / "relaxtime" / "literature" / "meson_density_targets"
    reports = normalize_wpd_csv(source_csv, output_dir)
    for report in reports:
        print(
            f"{report['slug']}: kept {report['kept_points']}/{report['raw_points']} "
            f"(dropped_nonmonotonic={report['dropped_nonmonotonic']}, "
            f"clamped_negative_to_zero={report['clamped_negative_to_zero']}, "
            f"skipped_blank={report['skipped_blank']})"
        )


if __name__ == "__main__":
    main()
