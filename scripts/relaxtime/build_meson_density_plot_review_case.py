#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SOURCE_DIR = (
    PROJECT_ROOT
    / "data"
    / "outputs"
    / "results"
    / "relaxtime"
    / "meson_density"
    / "freezeout_validation"
    / "blaschke2019col_kminus_piminus_mu_pi_100_phase_shift_gbu_default"
)
OUT_DIR = (
    PROJECT_ROOT
    / "data"
    / "outputs"
    / "results"
    / "relaxtime"
    / "meson_density"
    / "plot_review"
    / "freezeout_kminus_piminus_mu_pi_100"
)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def _write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _to_float(row: dict[str, str], key: str) -> float:
    return float(row[key])


def build_plot_review_case() -> dict[str, Path]:
    comparison_rows = _read_csv(SOURCE_DIR / "comparison_vs_target.csv")
    workflow_rows = _read_csv(SOURCE_DIR / "workflow_scan.csv")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    plot_rows: list[dict[str, object]] = []
    for row in comparison_rows:
        target_y = _to_float(row, "target_y")
        model_y = _to_float(row, "model_y")
        abs_diff = _to_float(row, "abs_diff")
        rel_diff_raw = row["rel_diff"].strip()
        rel_diff = "" if rel_diff_raw == "NaN" else float(rel_diff_raw)
        plot_rows.append(
            {
                "x_value": _to_float(row, "x_value"),
                "target_y": target_y,
                "model_y": model_y,
                "abs_diff": abs_diff,
                "rel_diff": rel_diff,
                "T_MeV": _to_float(row, "T_MeV"),
                "muB_MeV": _to_float(row, "muB_MeV"),
            }
        )

    plot_csv = OUT_DIR / "plot_review_comparison.csv"
    _write_csv(
        plot_csv,
        plot_rows,
        ["x_value", "target_y", "model_y", "abs_diff", "rel_diff", "T_MeV", "muB_MeV"],
    )

    workflow_csv = OUT_DIR / "workflow_scan.csv"
    workflow_csv.write_text((SOURCE_DIR / "workflow_scan.csv").read_text(encoding="utf-8"), encoding="utf-8")

    comparison_csv = OUT_DIR / "comparison_vs_target.csv"
    comparison_csv.write_text((SOURCE_DIR / "comparison_vs_target.csv").read_text(encoding="utf-8"), encoding="utf-8")

    xs = [row["x_value"] for row in plot_rows]
    target = [row["target_y"] for row in plot_rows]
    model = [row["model_y"] for row in plot_rows]
    abs_diff = [row["abs_diff"] for row in plot_rows]
    rel_diff = [row["rel_diff"] if row["rel_diff"] != "" else None for row in plot_rows]

    fig, ax = plt.subplots(figsize=(8.6, 5.2), constrained_layout=True)
    ax.plot(xs, target, label="literature target", linewidth=2.0, color="#1f77b4")
    ax.plot(xs, model, label="workflow model", linewidth=2.0, color="#d62728")
    ax.set_xlabel("sqrt(s_NN) [GeV]")
    ax.set_ylabel("K- / pi-")
    ax.set_title("Freeze-out plot review: kminus/piminus, mu_pi=100 MeV")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    overlay_png = OUT_DIR / "overlay_kminus_piminus_mu_pi_100.png"
    fig.savefig(overlay_png, dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(8.6, 6.8), constrained_layout=True, sharex=True)
    axes[0].plot(xs, abs_diff, linewidth=2.0, color="#9467bd")
    axes[0].set_ylabel("abs diff")
    axes[0].grid(True, alpha=0.25)
    rel_x = [x for x, y in zip(xs, rel_diff) if y is not None]
    rel_y = [y for y in rel_diff if y is not None]
    axes[1].plot(rel_x, rel_y, linewidth=2.0, color="#2ca02c")
    axes[1].set_xlabel("sqrt(s_NN) [GeV]")
    axes[1].set_ylabel("rel diff")
    axes[1].grid(True, alpha=0.25)
    residual_png = OUT_DIR / "residual_kminus_piminus_mu_pi_100.png"
    fig.savefig(residual_png, dpi=160)
    plt.close(fig)

    max_abs_diff = max(abs_diff)
    max_rel_diff = max(y for y in rel_y) if rel_y else float("nan")
    readme = OUT_DIR / "README.md"
    readme.write_text(
        "\n".join(
            [
                "# Plot Review Case: freezeout_kminus_piminus_mu_pi_100",
                "",
                "Current role:",
                "",
                "- plot-review asset for manual trend inspection",
                "- not an external validation gate",
                "- not a regression truth source by itself",
                "",
                "Inputs:",
                "",
                f"- workflow source: `{SOURCE_DIR.relative_to(PROJECT_ROOT) / 'workflow_scan.csv'}`",
                f"- comparison source: `{SOURCE_DIR.relative_to(PROJECT_ROOT) / 'comparison_vs_target.csv'}`",
                "",
                "Outputs:",
                "",
                "- `workflow_scan.csv`",
                "- `comparison_vs_target.csv`",
                "- `plot_review_comparison.csv`",
                "- `overlay_kminus_piminus_mu_pi_100.png`",
                "- `residual_kminus_piminus_mu_pi_100.png`",
                "",
                "Manual review checklist:",
                "",
                "- compare target/model monotonic rise",
                "- compare whether the workflow remains systematically below the literature target",
                "- compare whether no nonphysical reversals or spikes appear",
                "- compare whether the relative-difference band shape changes after code updates",
                "",
                "Summary metrics:",
                "",
                f"- points: `{len(plot_rows)}`",
                f"- max abs diff: `{max_abs_diff:.6f}`",
                f"- max rel diff (finite only): `{max_rel_diff:.6f}`",
                "",
                "Interpretation note:",
                "",
                "- This directory supports regression fallback and manual plot review only.",
                "- It does not promote the literature curve to `tests/validation/`.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    return {
        "plot_csv": plot_csv,
        "workflow_csv": workflow_csv,
        "comparison_csv": comparison_csv,
        "overlay_png": overlay_png,
        "residual_png": residual_png,
        "readme": readme,
    }


def main() -> None:
    outputs = build_plot_review_case()
    for key, path in outputs.items():
        print(f"{key}={path}")


if __name__ == "__main__":
    main()
