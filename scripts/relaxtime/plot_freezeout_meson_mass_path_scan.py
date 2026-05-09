#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CASE = "default_baseline_freezeout_xi0_loggrid_1to200_n30"
DEFAULT_CSV = (
    PROJECT_ROOT
    / "data"
    / "outputs"
    / "results"
    / "relaxtime"
    / "meson_mass"
    / "path_scan"
    / "freezeout"
    / f"{DEFAULT_CASE}.csv"
)

FREEZEOUT_DEFAULT_PROFILE = {
    "a_GeV": 0.166,
    "b_GeV_inv1": 0.139,
    "c_GeV_inv3": 0.053,
    "d_GeV": 1.308,
    "e_GeV_inv1": 0.273,
}


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


def _relative_posix(path: Path) -> str:
    return path.relative_to(PROJECT_ROOT).as_posix()


def _git(cmd: list[str]) -> str:
    try:
        return subprocess.check_output(cmd, cwd=PROJECT_ROOT, text=True, encoding="utf-8").strip()
    except Exception:
        return "unknown"


def _freezeout_temperature_mev(sqrt_s_gev: float) -> float:
    coeff = FREEZEOUT_DEFAULT_PROFILE
    mu_b = coeff["d_GeV"] / (1.0 + coeff["e_GeV_inv1"] * sqrt_s_gev)
    t_gev = coeff["a_GeV"] - coeff["b_GeV_inv1"] * mu_b**2 - coeff["c_GeV_inv3"] * mu_b**4
    return 1000.0 * t_gev


def _solve_zero_temperature_sqrts() -> float:
    lo = 1.2
    hi = 1.5
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if _freezeout_temperature_mev(mid) > 0.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def build_case_artifacts(csv_path: Path) -> dict[str, Path]:
    rows = _read_csv(csv_path)
    if not rows:
        raise ValueError(f"CSV is empty: {csv_path}")

    out_dir = csv_path.parent
    case_name = csv_path.stem
    figures_dir = out_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    valid_rows: list[dict[str, str]] = []
    invalid_rows: list[dict[str, str]] = []
    for row in rows:
        is_valid = row["equilibrium_converged"].strip().lower() == "true" and _to_float(row, "T_MeV") > 0.0
        if is_valid:
            valid_rows.append(row)
        else:
            invalid_rows.append(row)
    if not valid_rows:
        raise ValueError(f"No valid rows available for plotting: {csv_path}")

    valid_plot_rows: list[dict[str, object]] = []
    for row in valid_rows:
        valid_plot_rows.append(
            {
                "path_point_index": int(row["path_point_index"]),
                "sqrt_s_NN_GeV": _to_float(row, "sqrt_s_NN_GeV"),
                "T_MeV": _to_float(row, "T_MeV"),
                "muB_MeV": _to_float(row, "muB_MeV"),
                "M_pi": _to_float(row, "M_pi"),
                "M_K": _to_float(row, "M_K"),
                "M_eta": _to_float(row, "M_eta"),
                "M_eta_prime": _to_float(row, "M_eta_prime"),
                "M_sigma_pi": _to_float(row, "M_sigma_pi"),
                "M_sigma_K": _to_float(row, "M_sigma_K"),
                "M_sigma": _to_float(row, "M_sigma"),
                "M_sigma_prime": _to_float(row, "M_sigma_prime"),
            }
        )

    valid_csv = out_dir / f"{case_name}.valid_points.csv"
    _write_csv(
        valid_csv,
        valid_plot_rows,
        [
            "path_point_index",
            "sqrt_s_NN_GeV",
            "T_MeV",
            "muB_MeV",
            "M_pi",
            "M_K",
            "M_eta",
            "M_eta_prime",
            "M_sigma_pi",
            "M_sigma_K",
            "M_sigma",
            "M_sigma_prime",
        ],
    )

    xs = [row["sqrt_s_NN_GeV"] for row in valid_plot_rows]
    Ts = [row["T_MeV"] for row in valid_plot_rows]
    muBs = [row["muB_MeV"] for row in valid_plot_rows]

    fig, axes = plt.subplots(2, 1, figsize=(8.8, 7.0), constrained_layout=True, sharex=True)
    axes[0].plot(xs, Ts, linewidth=2.0, color="#d62728")
    axes[0].set_ylabel("T [MeV]")
    axes[0].set_title(f"Freeze-out path overview: {case_name}")
    axes[0].grid(True, alpha=0.25)
    axes[1].plot(xs, muBs, linewidth=2.0, color="#1f77b4")
    axes[1].set_xlabel("sqrt(s_NN) [GeV]")
    axes[1].set_ylabel("mu_B [MeV]")
    axes[1].grid(True, alpha=0.25)
    path_png = figures_dir / f"{case_name}_path_overview.png"
    fig.savefig(path_png, dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(8.8, 7.4), constrained_layout=True, sharex=True)
    axes[0].plot(xs, [row["M_pi"] for row in valid_plot_rows], label="pi", linewidth=2.0)
    axes[0].plot(xs, [row["M_K"] for row in valid_plot_rows], label="K", linewidth=2.0)
    axes[0].set_ylabel("mass [fm^-1]")
    axes[0].set_title(f"Meson masses on freeze-out path: {case_name}")
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(frameon=False, ncol=2)
    axes[1].plot(xs, [row["M_eta"] for row in valid_plot_rows], label="eta", linewidth=2.0)
    axes[1].plot(xs, [row["M_eta_prime"] for row in valid_plot_rows], label="eta_prime", linewidth=2.0)
    axes[1].set_xlabel("sqrt(s_NN) [GeV]")
    axes[1].set_ylabel("mass [fm^-1]")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(frameon=False, ncol=2)
    pseudoscalar_png = figures_dir / f"{case_name}_pseudoscalar_masses.png"
    fig.savefig(pseudoscalar_png, dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(8.8, 7.4), constrained_layout=True, sharex=True)
    axes[0].plot(xs, [row["M_sigma_pi"] for row in valid_plot_rows], label="sigma_pi", linewidth=2.0)
    axes[0].plot(xs, [row["M_sigma_K"] for row in valid_plot_rows], label="sigma_K", linewidth=2.0)
    axes[0].set_ylabel("mass [fm^-1]")
    axes[0].set_title(f"Scalar-channel masses on freeze-out path: {case_name}")
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(frameon=False, ncol=2)
    axes[1].plot(xs, [row["M_sigma"] for row in valid_plot_rows], label="sigma", linewidth=2.0)
    axes[1].plot(xs, [row["M_sigma_prime"] for row in valid_plot_rows], label="sigma_prime", linewidth=2.0)
    axes[1].set_xlabel("sqrt(s_NN) [GeV]")
    axes[1].set_ylabel("mass [fm^-1]")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(frameon=False, ncol=2)
    scalar_png = figures_dir / f"{case_name}_scalar_masses.png"
    fig.savefig(scalar_png, dpi=160)
    plt.close(fig)

    zero_t_threshold = _solve_zero_temperature_sqrts()
    recommended_physical_lower_bound = 1.45
    invalid_summary = [
        {
            "path_point_index": int(row["path_point_index"]),
            "sqrt_s_NN_GeV": _to_float(row, "sqrt_s_NN_GeV"),
            "T_MeV": _to_float(row, "T_MeV"),
            "muB_MeV": _to_float(row, "muB_MeV"),
            "message": row["message"],
        }
        for row in invalid_rows
    ]

    provenance = {
        "schema_version": "v1",
        "case_name": case_name,
        "generated_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "source_csv": _relative_posix(csv_path),
        "valid_points_csv": _relative_posix(valid_csv),
        "figures": [
            _relative_posix(path_png),
            _relative_posix(pseudoscalar_png),
            _relative_posix(scalar_png),
        ],
        "git_commit": _git(["git", "rev-parse", "HEAD"]),
        "git_branch": _git(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
        "points_total": len(rows),
        "points_plotted": len(valid_rows),
        "points_excluded": len(invalid_rows),
        "exclusion_rule": "exclude rows with equilibrium_converged != true or T_MeV <= 0.0",
        "known_invalid_points": invalid_summary,
        "freezeout_default_zero_temperature_threshold_sqrt_s_NN_GeV": zero_t_threshold,
        "recommended_physical_lower_bound_sqrt_s_NN_GeV": recommended_physical_lower_bound,
        "note": "Formal production CSV remains [1, 200] GeV as requested; future physical freezeout scans should use sqrt(s_NN) >= 1.45 GeV for the default baseline profile.",
    }
    provenance_json = out_dir / f"{case_name}.provenance.json"
    provenance_json.write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    readme = out_dir / f"{case_name}.README.md"
    readme.write_text(
        "\n".join(
            [
                f"# Freeze-out Meson Mass Case: {case_name}",
                "",
                "Current role:",
                "",
                "- formal freezeout meson-mass path output",
                "- xi = 0.0",
                "- sqrt(s_NN) sampled on a log10-uniform grid over [1, 200] GeV with 30 points",
                "",
                "Known path-validity note:",
                "",
                "- The formal production CSV is intentionally kept unchanged on `[1, 200] GeV`.",
                "- For the `default` freezeout baseline profile, the first two sampled points produce nonpositive `T_MeV` and are therefore retained as explicit failure rows in the CSV.",
                f"- The default-profile `T=0` threshold is approximately `sqrt(s_NN) = {zero_t_threshold:.6f} GeV`.",
                f"- For future physical freezeout scans, use `sqrt(s_NN) >= {recommended_physical_lower_bound:.2f} GeV` for this profile family.",
                "",
                "Plotting policy:",
                "",
                "- Figures in `figures/` exclude rows with `equilibrium_converged != true` or `T_MeV <= 0`.",
                f"- Excluded points in this case: `{len(invalid_rows)}`.",
                "",
                "Files:",
                "",
                f"- source CSV: `{_relative_posix(csv_path)}`",
                f"- valid-points CSV: `{_relative_posix(valid_csv)}`",
                f"- provenance: `{_relative_posix(provenance_json)}`",
                f"- path overview: `{_relative_posix(path_png)}`",
                f"- pseudoscalar masses: `{_relative_posix(pseudoscalar_png)}`",
                f"- scalar masses: `{_relative_posix(scalar_png)}`",
                "",
                "Quick summary:",
                "",
                f"- total points: `{len(rows)}`",
                f"- plotted valid points: `{len(valid_rows)}`",
                f"- excluded invalid points: `{len(invalid_rows)}`",
                f"- first plotted point: `sqrt(s_NN) = {valid_plot_rows[0]['sqrt_s_NN_GeV']:.6f} GeV`",
                f"- last plotted point: `sqrt(s_NN) = {valid_plot_rows[-1]['sqrt_s_NN_GeV']:.6f} GeV`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    return {
        "valid_csv": valid_csv,
        "path_png": path_png,
        "pseudoscalar_png": pseudoscalar_png,
        "scalar_png": scalar_png,
        "provenance_json": provenance_json,
        "readme": readme,
    }


def main() -> None:
    outputs = build_case_artifacts(DEFAULT_CSV)
    for key, path in outputs.items():
        print(f"{key}={path}")


if __name__ == "__main__":
    main()
