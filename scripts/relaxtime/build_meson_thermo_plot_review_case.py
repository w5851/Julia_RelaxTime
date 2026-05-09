#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE_DIR = (
    PROJECT_ROOT
    / "data"
    / "outputs"
    / "results"
    / "relaxtime"
    / "meson_thermo"
    / "canonical_muB0_phase_shift_current_pi_sigma_pi"
)
DEFAULT_OUT_DIR = (
    PROJECT_ROOT
    / "data"
    / "outputs"
    / "results"
    / "relaxtime"
    / "meson_thermo"
    / "plot_review"
    / "canonical_muB0_phase_shift_current_pi_sigma_pi"
)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        header = None
        lines = []
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if header is None:
                header = [x.strip() for x in line.split(",")]
                continue
            cols = [x.strip() for x in line.split(",")]
            if len(cols) != len(header):
                raise ValueError(f"invalid CSV row in {path}: {raw!r}")
            lines.append(dict(zip(header, cols)))
        return lines


def _write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _to_float(row: dict[str, str], key: str) -> float:
    return float(row[key])


def _relative_posix(path: Path) -> str:
    try:
        return path.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def build_plot_review_case(source_dir: Path, out_dir: Path) -> dict[str, Path]:
    scan_csv = source_dir / "scan.csv"
    readme_src = source_dir / "README.md"
    if not scan_csv.is_file():
        raise FileNotFoundError(f"source scan.csv not found: {scan_csv}")
    if not readme_src.is_file():
        raise FileNotFoundError(f"source README.md not found: {readme_src}")

    rows = _read_csv(scan_csv)
    if not rows:
        raise ValueError(f"source scan.csv has no data rows: {scan_csv}")

    out_dir.mkdir(parents=True, exist_ok=True)

    workflow_csv = out_dir / "workflow_scan.csv"
    shutil.copyfile(scan_csv, workflow_csv)

    summary_rows: list[dict[str, object]] = []
    for row in rows:
        p_total = _to_float(row, "P_total")
        p_q = _to_float(row, "P_quark_meanfield")
        p_m = _to_float(row, "P_meson")
        p_qp = _to_float(row, "P_meson_qp")
        p_ld = _to_float(row, "P_meson_ld")
        trace = _to_float(row, "trace_anomaly")
        qp_share = float("nan") if p_m == 0.0 else p_qp / p_m
        ld_share = float("nan") if p_m == 0.0 else p_ld / p_m
        summary_rows.append(
            {
                "T_MeV": _to_float(row, "T_MeV"),
                "P_total": p_total,
                "P_quark_meanfield": p_q,
                "P_meson": p_m,
                "P_meson_qp": p_qp,
                "P_meson_ld": p_ld,
                "trace_anomaly": trace,
                "qp_share": qp_share,
                "ld_share": ld_share,
                "equilibrium_converged": row["equilibrium_converged"],
                "thermo_derivation_mode": row["thermo_derivation_mode"],
            }
        )

    summary_csv = out_dir / "plot_review_summary.csv"
    _write_csv(
        summary_csv,
        summary_rows,
        [
            "T_MeV",
            "P_total",
            "P_quark_meanfield",
            "P_meson",
            "P_meson_qp",
            "P_meson_ld",
            "trace_anomaly",
            "qp_share",
            "ld_share",
            "equilibrium_converged",
            "thermo_derivation_mode",
        ],
    )

    xs = [float(row["T_MeV"]) for row in summary_rows]
    p_total = [float(row["P_total"]) for row in summary_rows]
    p_q = [float(row["P_quark_meanfield"]) for row in summary_rows]
    p_m = [float(row["P_meson"]) for row in summary_rows]
    p_qp = [float(row["P_meson_qp"]) for row in summary_rows]
    p_ld = [float(row["P_meson_ld"]) for row in summary_rows]
    trace = [float(row["trace_anomaly"]) for row in summary_rows]
    qp_share = [float(row["qp_share"]) for row in summary_rows]
    ld_share = [float(row["ld_share"]) for row in summary_rows]

    fig, ax = plt.subplots(figsize=(8.8, 5.4), constrained_layout=True)
    ax.plot(xs, p_total, label="P_total", linewidth=2.1, color="#1f77b4")
    ax.plot(xs, p_q, label="P_quark_meanfield", linewidth=2.0, color="#7f7f7f")
    ax.plot(xs, p_m, label="P_meson", linewidth=2.0, color="#d62728")
    ax.set_xlabel("T [MeV]")
    ax.set_ylabel("pressure [fm^-4]")
    ax.set_title("Canonical mu_B=0 meson thermo: pressure decomposition")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    pressure_png = out_dir / "pressure_overlay.png"
    fig.savefig(pressure_png, dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.8, 5.4), constrained_layout=True)
    ax.plot(xs, trace, linewidth=2.0, color="#9467bd")
    ax.axhline(0.0, linewidth=1.0, color="#666666", alpha=0.6)
    ax.set_xlabel("T [MeV]")
    ax.set_ylabel("(epsilon - 3P) / T^4")
    ax.set_title("Canonical mu_B=0 meson thermo: trace anomaly")
    ax.grid(True, alpha=0.25)
    trace_png = out_dir / "trace_anomaly_overlay.png"
    fig.savefig(trace_png, dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(8.8, 6.9), constrained_layout=True, sharex=True)
    axes[0].plot(xs, p_qp, linewidth=2.0, color="#2ca02c", label="P_meson_qp")
    axes[0].plot(xs, p_ld, linewidth=2.0, color="#ff7f0e", label="P_meson_ld")
    axes[0].set_ylabel("pressure [fm^-4]")
    axes[0].set_title("Canonical mu_B=0 meson thermo: QP/LD split")
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(frameon=False)
    axes[1].plot(xs, qp_share, linewidth=2.0, color="#2ca02c", label="QP share")
    axes[1].plot(xs, ld_share, linewidth=2.0, color="#ff7f0e", label="LD share")
    axes[1].set_xlabel("T [MeV]")
    axes[1].set_ylabel("share of P_meson")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(frameon=False)
    split_png = out_dir / "qp_ld_split.png"
    fig.savefig(split_png, dpi=160)
    plt.close(fig)

    max_total = max(p_total)
    min_trace = min(trace)
    max_ld_share = max(ld_share)
    min_qp_share = min(qp_share)

    readme = out_dir / "README.md"
    readme.write_text(
        "\n".join(
            [
                "# Plot Review Case: canonical_muB0_phase_shift_current_pi_sigma_pi",
                "",
                "Current role:",
                "",
                "- plot-review asset for manual trend inspection",
                "- not an external validation gate",
                "- not a replacement for fixedpoint/path regression",
                "",
                "Inputs:",
                "",
                f"- workflow source: `{_relative_posix(scan_csv)}`",
                f"- source README: `{_relative_posix(readme_src)}`",
                "",
                "Outputs:",
                "",
                "- `workflow_scan.csv`",
                "- `plot_review_summary.csv`",
                "- `pressure_overlay.png`",
                "- `trace_anomaly_overlay.png`",
                "- `qp_ld_split.png`",
                "- `README.md`",
                "",
                "Manual review checklist:",
                "",
                "- compare whether `P_total` remains above `P_quark_meanfield` across the path",
                "- compare whether `P_meson` decreases with temperature in the current short path",
                "- compare whether `trace_anomaly` remains negative and approaches zero from below",
                "- compare whether the QP share grows while the LD share shrinks along the path",
                "- compare whether no nonphysical jumps or sign flips appear after code updates",
                "",
                "Summary metrics:",
                "",
                f"- points: `{len(summary_rows)}`",
                f"- max P_total: `{max_total:.6f}`",
                f"- min trace anomaly: `{min_trace:.6f}`",
                f"- max LD share: `{max_ld_share:.6f}`",
                f"- min QP share: `{min_qp_share:.6f}`",
                "",
                "Interpretation note:",
                "",
                "- This directory supports canonical trend review and regression fallback only.",
                "- It does not promote the current project-internal curve to `tests/validation/`.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    return {
        "workflow_csv": workflow_csv,
        "summary_csv": summary_csv,
        "pressure_png": pressure_png,
        "trace_png": trace_png,
        "split_png": split_png,
        "readme": readme,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()
    outputs = build_plot_review_case(args.source_dir.resolve(), args.out_dir.resolve())
    for key, path in outputs.items():
        print(f"{key}={path}")


if __name__ == "__main__":
    main()
