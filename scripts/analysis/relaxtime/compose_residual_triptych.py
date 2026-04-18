#!/usr/bin/env python3

import argparse
import pathlib

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _load_grid(df: pd.DataFrame):
    m_vals = np.sort(df["m_fm"].unique())
    g_vals = np.sort(df["g_fm"].unique())
    re_grid = df.pivot(index="g_fm", columns="m_fm", values="residual_re").values
    im_grid = df.pivot(index="g_fm", columns="m_fm", values="residual_im").values
    n_grid = df.pivot(index="g_fm", columns="m_fm", values="residual_norm").values
    return m_vals, g_vals, re_grid, im_grid, n_grid


def _symlog_norm(z: np.ndarray):
    vmax = np.nanmax(np.abs(z))
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = 1.0
    linthresh = max(vmax * 1e-3, 1e-12)
    return mcolors.SymLogNorm(linthresh=linthresh, vmin=-vmax, vmax=vmax, base=10)


def _draw_triptych(surface_csv: pathlib.Path, out_png: pathlib.Path, title_prefix: str):
    df = pd.read_csv(surface_csv)
    m_vals, g_vals, re_grid, im_grid, n_grid = _load_grid(df)
    m_mesh, g_mesh = np.meshgrid(m_vals, g_vals)

    fig, axes = plt.subplots(1, 3, figsize=(20, 6), constrained_layout=True)

    norm_re = _symlog_norm(re_grid)
    im_re = axes[0].imshow(
        re_grid,
        origin="lower",
        aspect="auto",
        extent=[m_vals.min(), m_vals.max(), g_vals.min(), g_vals.max()],
        cmap="RdBu_r",
        norm=norm_re,
    )
    axes[0].set_title("Re residual (symlog)")
    axes[0].set_xlabel("mass m (fm^-1)")
    axes[0].set_ylabel("gamma (fm^-1)")
    cbar_re = fig.colorbar(im_re, ax=axes[0], fraction=0.046, pad=0.04)
    cbar_re.set_label("residual_re")

    norm_im = _symlog_norm(im_grid)
    im_im = axes[1].imshow(
        im_grid,
        origin="lower",
        aspect="auto",
        extent=[m_vals.min(), m_vals.max(), g_vals.min(), g_vals.max()],
        cmap="RdBu_r",
        norm=norm_im,
    )
    axes[1].set_title("Im residual (symlog)")
    axes[1].set_xlabel("mass m (fm^-1)")
    axes[1].set_ylabel("gamma (fm^-1)")
    cbar_im = fig.colorbar(im_im, ax=axes[1], fraction=0.046, pad=0.04)
    cbar_im.set_label("residual_im")

    n_plot = np.log10(np.maximum(n_grid, 1e-16))
    im_n = axes[2].imshow(
        n_plot,
        origin="lower",
        aspect="auto",
        extent=[m_vals.min(), m_vals.max(), g_vals.min(), g_vals.max()],
        cmap="viridis",
    )
    fig.colorbar(im_n, ax=axes[2], fraction=0.046, pad=0.04, label="log10 ||F||")
    levels = np.linspace(np.nanmin(n_plot), np.nanmax(n_plot), 12)
    axes[2].contour(m_mesh, g_mesh, n_plot, levels=levels, colors="white", linewidths=0.35, alpha=0.45)
    axes[2].contour(m_mesh, g_mesh, re_grid, levels=[0.0], colors=["crimson"], linewidths=2.0)
    axes[2].contour(m_mesh, g_mesh, im_grid, levels=[0.0], colors=["royalblue"], linewidths=2.0)

    idx = np.unravel_index(np.argmin(n_grid), n_grid.shape)
    g0 = g_vals[idx[0]]
    m0 = m_vals[idx[1]]
    axes[2].plot([m0], [g0], marker="o", color="black", ms=4)
    axes[2].set_title("Heatmap + Re/Im zero contours")
    axes[2].set_xlabel("mass m (fm^-1)")
    axes[2].set_ylabel("gamma (fm^-1)")

    fig.suptitle(title_prefix)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Compose Re/Im/contour triptych from residual surface CSV")
    parser.add_argument("--surface", required=True, help="Path to residual_surface.csv")
    parser.add_argument("--out", required=True, help="Output PNG")
    parser.add_argument("--title", default="Residual triptych", help="Figure title")
    args = parser.parse_args()

    _draw_triptych(pathlib.Path(args.surface), pathlib.Path(args.out), args.title)
    print(f"Wrote triptych: {args.out}")


if __name__ == "__main__":
    main()
