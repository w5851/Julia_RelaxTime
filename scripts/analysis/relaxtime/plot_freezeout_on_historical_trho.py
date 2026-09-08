"""Overlay saved BQS coordinates, not new ratios, on the retained historical grid."""
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as effects
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
import numpy as np

import render_combined_meson_density_fig3_like as historical

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "data/outputs/results/relaxtime/analysis/charged_rpa_phase_backend"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    source = Path(os.environ.get("FREEZEOUT_TRHO_OUTPUT", BASE / "freezeout_on_historical_trho_20260905"))
    out = source / "plots"
    if out.exists():
        raise RuntimeError(f"refusing to overwrite {out}")
    coords_path = source / "freezeout_trho_coordinates.csv"
    manifest = json.loads((source / "manifest.json").read_text())
    assert sha(coords_path) == manifest["output_hashes"][coords_path.name]
    with coords_path.open(newline="", encoding="utf-8") as f:
        points = list(csv.DictReader(f))
    assert len(points) == 10 and all(p["passed"] == "true" for p in points)
    historical_csv = ROOT / "data/outputs/results/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/combined_meson_density_scan.csv"
    plot_manifest = ROOT / "data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/plot_manifest.json"
    original_png = plot_manifest.parent / "combined_meson_density_scan.png"
    original_svg = plot_manifest.parent / "combined_meson_density_scan.svg"
    protected = {str(p): sha(p) for p in (historical_csv, plot_manifest, original_png, original_svg)}
    old = historical.read_rows(historical_csv)
    pm = json.loads(plot_manifest.read_text())
    spec = next(p for p in pm["figures"] if p["format"] == "png")
    norm = LogNorm(vmin=spec["color_vmin"], vmax=spec["color_vmax"])
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("0.7")
    x = np.array([float(p["rho_norm"]) for p in points])
    y = np.array([float(p["T_MeV"]) for p in points])
    energy = np.array([float(p["sqrt_s_NN_GeV"]) for p in points])
    xmax = max(1.03, 1.15 * x.max())
    ymin = min(115., y.min() - 8.)

    def panel(ax, regime, zoom=False, detailed=False):
        xs, ts, values = historical.matrix(old, regime, "kpi_ratio", "rho_target")
        xe, te = historical.centers_to_edges(xs), historical.centers_to_edges(ts)
        # Do not let outer half-cells suggest data below T=120 or rho/rho0=0.05.
        xe[0], xe[-1], te[0], te[-1] = xs[0], xs[-1], ts[0], ts[-1]
        mat = np.ma.masked_where(~np.isfinite(values) | (values <= 0), values)
        im = ax.pcolormesh(xe, te, mat, cmap=cmap, norm=norm, shading="flat", rasterized=True)
        ax.set_facecolor("#f8f8f8")
        ax.add_patch(Rectangle((.05, 120), .95, 100, fill=False, edgecolor="#666666", linestyle="--", linewidth=1))
        ax.plot(x, y, color="white", lw=4.5, zorder=4)
        ax.plot(x, y, "-o", color="#111111", mfc="white", mec="#111111", ms=5, lw=1.5,
                zorder=5, label="Current BQS freezeout (10 saved states)")
        outside = np.array([p["coverage"] != "inside_historical_grid" for p in points])
        ax.scatter(x[outside], y[outside], marker="D", facecolors="#ff8a3d", edgecolors="#222222", s=43, zorder=6)
        ax.set_xlim(0, max(.065, x.max()*1.65) if zoom else xmax)
        ax.set_ylim(ymin, 181 if zoom else 225)
        ax.set_xlabel(r"Net baryon density $\rho_B/\rho_0$  ($\rho_0=0.16$ fm$^{-3}$)")
        ax.set_ylabel("T [MeV]")
        ax.set_title(("Trajectory detail" if zoom else regime.replace("_", " ")), fontsize=12)
        if detailed:
            # A fixed label rail keeps the close high-energy points distinguishable.
            for idx, yt in zip(np.argsort(-energy), np.linspace(176, ymin+8, len(points))):
                ax.annotate(f"{energy[idx]:g} GeV", (x[idx], y[idx]), (ax.get_xlim()[1]*.70, yt),
                            fontsize=9, va="center", bbox=dict(fc="white", ec="none", alpha=.85, pad=1),
                            arrowprops=dict(arrowstyle="-", color="#555555", lw=.65), zorder=8)
        else:
            for e, delta in ((3, (9, 0)), (7.7, (10, -15)), (200, (10, 13))):
                i = int(np.where(energy == e)[0][0])
                txt = ax.annotate(f"{e:g} GeV", (x[i], y[i]), xytext=delta, textcoords="offset points", fontsize=9, zorder=8)
                txt.set_path_effects([effects.withStroke(linewidth=3, foreground="white")])
        return im

    out.mkdir()
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 7.6), layout="constrained")
    panel(axes[0], "phase_shift_gbu_reference")
    im = panel(axes[1], "phase_shift_gbu_reference", zoom=True, detailed=True)
    axes[0].legend(loc="upper right", fontsize=8)
    fig.colorbar(im, ax=axes, label="Historical K+/pi+ (log scale; unchanged limits)", shrink=.82)
    fig.suptitle("Historical GBU heatmap + current chemical-freezeout trajectory\n"
                 "Color = OLD ratio; line = background coordinates only, NOT the new ratio", fontsize=14)
    fig.supxlabel("Old: quark-only rho_u/rho_d=0.876 | Current: BQS rho_Q/rho_B=0.4 (u/d=0.875); rho_S=0\n"
                  "Orange diamonds: outside historical sampled rectangle. Pale area: no historical data. Lines only guide the eye.", fontsize=9)
    fig.savefig(out / "historical_gbu_with_freezeout.png", dpi=180)
    fig.savefig(out / "historical_gbu_with_freezeout.svg", metadata={"Date": None})
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 11), layout="constrained")
    for ax, regime in zip(axes.ravel(), ("stable", "strict_bw_stage1", "phase_shift_current", "phase_shift_gbu_reference")):
        im = panel(ax, regime)
    fig.colorbar(im, ax=axes.ravel().tolist(), label="Historical K+/pi+ (log scale)", shrink=.8)
    fig.suptitle("Historical T-rho panels with the SAME current BQS freezeout trajectory\n"
                 "Background-coordinate overlay only; no meson calculation or historical-grid interpolation", fontsize=13)
    fig.supxlabel("Old u/d=0.876; current u/d=0.875; both quark-only with rho_S=0. Orange: outside sampled rectangle.\n"
                  "Gray cells: historical missing/invalid ratios. Pale exterior: no historical data. Color limits retained.", fontsize=9)
    fig.savefig(out / "historical_all_regimes_with_freezeout.png", dpi=180)
    plt.close(fig)
    assert all(sha(p) == h for p, h in protected.items())
    report = {"status": "diagnostic_coordinate_overlay", "source_coordinate_sha256": sha(coords_path),
              "source_manifest_sha256": sha(source / "manifest.json"), "historical_input_hashes": protected,
              "script_sha256": sha(__file__), "renderer_sha256": sha(historical.__file__),
              "color_vmin": norm.vmin, "color_vmax": norm.vmax, "points": len(points),
              "inside_historical_grid": sum(p["coverage"] == "inside_historical_grid" for p in points),
              "outside_historical_grid": sum(p["coverage"] != "inside_historical_grid" for p in points),
              "historical_ratio_interpolation": False, "new_ratio_overlay": False, "production_authorized": False,
              "cell_boundary_policy": "outer half-cells clipped to the sampled centre rectangle; no extrapolation",
              "output_hashes": {p.name: sha(p) for p in out.iterdir() if p.is_file()}}
    (out / "manifest.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
