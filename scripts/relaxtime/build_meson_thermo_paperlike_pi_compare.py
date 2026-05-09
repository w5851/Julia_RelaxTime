#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _read_scan(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as f:
        header = None
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if header is None:
                header = [x.strip() for x in line.split(",")]
                continue
            cols = [x.strip() for x in line.split(",")]
            if len(cols) != len(header):
                raise ValueError(f"invalid row in {path}: {raw!r}")
            rows.append(dict(zip(header, cols)))
    if not rows:
        raise ValueError(f"no data rows in {path}")
    return rows


def _relative_posix(path: Path) -> str:
    try:
        return path.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def _series_from_scan(scan_csv: Path) -> dict[str, list[float]]:
    rows = _read_scan(scan_csv)
    xs: list[float] = []
    pi_qp_t4: list[float] = []
    pi_ld_t4: list[float] = []
    for row in rows:
        t_mev = float(row["T_MeV"])
        t_fm = float(row["T_fm"])
        t4 = t_fm**4
        xs.append(t_mev)
        pi_qp_t4.append(float(row["P_primary_qp"]) / t4)
        pi_ld_t4.append(float(row["P_primary_ld"]) / t4)
    return {
        "T_MeV": xs,
        "P_pi_qp_over_T4": pi_qp_t4,
        "P_pi_ld_over_T4": pi_ld_t4,
    }


def build_paperlike_compare(source_lambda: Path, source_2lambda: Path, out_dir: Path) -> dict[str, Path]:
    scan_lambda = source_lambda / "scan.csv"
    scan_2lambda = source_2lambda / "scan.csv"
    if not scan_lambda.is_file():
        raise FileNotFoundError(f"scan not found: {scan_lambda}")
    if not scan_2lambda.is_file():
        raise FileNotFoundError(f"scan not found: {scan_2lambda}")

    out_dir.mkdir(parents=True, exist_ok=True)

    data_lambda = _series_from_scan(scan_lambda)
    data_2lambda = _series_from_scan(scan_2lambda)

    fig, axes = plt.subplots(2, 1, figsize=(8.8, 7.0), constrained_layout=True, sharex=True)
    axes[0].plot(
        data_lambda["T_MeV"],
        data_lambda["P_pi_qp_over_T4"],
        linewidth=2.0,
        color="#1f77b4",
        label=r"$P_{\pi}^{QP}/T^4,\ \Lambda_{LD}=\Lambda$",
    )
    axes[0].plot(
        data_2lambda["T_MeV"],
        data_2lambda["P_pi_qp_over_T4"],
        linewidth=2.0,
        color="#1f77b4",
        linestyle="--",
        label=r"$P_{\pi}^{QP}/T^4,\ \Lambda_{LD}=2\Lambda$",
    )
    axes[0].set_ylabel(r"$P/T^4$")
    axes[0].set_title("Paper-like pion-only QP/LD comparison")
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(frameon=False)

    axes[1].plot(
        data_lambda["T_MeV"],
        data_lambda["P_pi_ld_over_T4"],
        linewidth=2.0,
        color="#d62728",
        label=r"$P_{\pi}^{LD}/T^4,\ \Lambda_{LD}=\Lambda$",
    )
    axes[1].plot(
        data_2lambda["T_MeV"],
        data_2lambda["P_pi_ld_over_T4"],
        linewidth=2.0,
        color="#d62728",
        linestyle="--",
        label=r"$P_{\pi}^{LD}/T^4,\ \Lambda_{LD}=2\Lambda$",
    )
    axes[1].set_xlabel("T [MeV]")
    axes[1].set_ylabel(r"$P/T^4$")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(frameon=False)

    compare_png = out_dir / "paperlike_pi_qp_ld_over_t4_compare.png"
    fig.savefig(compare_png, dpi=160)
    plt.close(fig)

    readme = out_dir / "README.md"
    readme.write_text(
        "\n".join(
            [
                "# Paper-like pion-only compare",
                "",
                "This directory compares only the pion-channel contributions extracted from canonical meson thermo scans.",
                "",
                "Inputs:",
                f"- lambda case: `{_relative_posix(scan_lambda)}`",
                f"- 2lambda case: `{_relative_posix(scan_2lambda)}`",
                "",
                "Outputs:",
                "- `paperlike_pi_qp_ld_over_t4_compare.png`",
                "- `README.md`",
                "",
                "Notes:",
                "- This is a paper-like comparison asset, not a validation artifact.",
                "- The source workflow still solves the project canonical meson-thermo point; this compare only extracts the pion-channel primary fields.",
                "- Curves are plotted as `P_pi^QP/T^4` and `P_pi^LD/T^4` for `Λ_LD = Λ` and `2Λ`.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    return {
        "compare_png": compare_png,
        "readme": readme,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-lambda", type=Path, required=True)
    parser.add_argument("--source-2lambda", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    outputs = build_paperlike_compare(
        args.source_lambda.resolve(),
        args.source_2lambda.resolve(),
        args.out_dir.resolve(),
    )
    for key, path in outputs.items():
        print(f"{key}={path}")


if __name__ == "__main__":
    main()
