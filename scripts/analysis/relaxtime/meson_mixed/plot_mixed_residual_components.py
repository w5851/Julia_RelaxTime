#!/usr/bin/env python3

import argparse
import pathlib

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _pivot(df: pd.DataFrame, value_col: str):
    m_vals = np.sort(df["m_fm"].unique())
    g_vals = np.sort(df["g_fm"].unique())
    z = df.pivot(index="g_fm", columns="m_fm", values=value_col).values
    return m_vals, g_vals, z


def _symlog_norm(z: np.ndarray):
    vmax = np.nanmax(np.abs(z))
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = 1.0
    linthresh = max(vmax * 1e-3, 1e-12)
    return mcolors.SymLogNorm(linthresh=linthresh, vmin=-vmax, vmax=vmax, base=10)


def _plot_component(df: pd.DataFrame, comp_col: str, out_png: pathlib.Path, title: str):
    m_vals, g_vals, z = _pivot(df, comp_col)
    fig, ax = plt.subplots(figsize=(8, 6))
    norm = _symlog_norm(z)
    im = ax.imshow(
        z,
        origin="lower",
        aspect="auto",
        extent=[m_vals.min(), m_vals.max(), g_vals.min(), g_vals.max()],
        cmap="RdBu_r",
        norm=norm,
    )

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(f"{comp_col} (symlog, sign-preserving)")

    ax.set_xlabel("mass m (fm^-1)")
    ax.set_ylabel("gamma (fm^-1)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Plot signed residual Re/Im heatmaps")
    parser.add_argument("--surface", required=True, help="Path to residual_surface.csv")
    parser.add_argument("--outdir", required=True, help="Output figure directory")
    parser.add_argument("--title-prefix", default="mixed residual", help="Figure title prefix")
    args = parser.parse_args()

    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.surface)

    _plot_component(
        df,
        "residual_re",
        outdir / "residual_re_heatmap_symlog.png",
        f"{args.title_prefix}: Re residual",
    )
    _plot_component(
        df,
        "residual_im",
        outdir / "residual_im_heatmap_symlog.png",
        f"{args.title_prefix}: Im residual",
    )

    print(f"Wrote component heatmaps under: {outdir}")


if __name__ == "__main__":
    main()
