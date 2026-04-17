#!/usr/bin/env python3

import argparse
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _plot_heatmap(df: pd.DataFrame, out_png: pathlib.Path, title: str):
    m_vals = np.sort(df["m_fm"].unique())
    g_vals = np.sort(df["g_fm"].unique())

    pivot = df.pivot(index="g_fm", columns="m_fm", values="log10_residual_norm")
    z = pivot.values

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.imshow(
        z,
        origin="lower",
        aspect="auto",
        extent=[m_vals.min(), m_vals.max(), g_vals.min(), g_vals.max()],
        cmap="viridis",
    )
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("log10 ||F||")

    ax.set_xlabel("mass m (fm^-1)")
    ax.set_ylabel("gamma (fm^-1)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    plt.close(fig)


def _plot_slice(df: pd.DataFrame, x_col: str, y_col: str, out_png: pathlib.Path, xlabel: str, title: str):
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(df[x_col].to_numpy(), df[y_col].to_numpy(), lw=1.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("log10 ||F||")
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_png, dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Plot mixed residual surface and slices")
    parser.add_argument("--surface", required=True, help="Path to residual_surface.csv")
    parser.add_argument("--slice-m", required=True, help="Path to residual_slice_mass.csv")
    parser.add_argument("--slice-g", required=True, help="Path to residual_slice_gamma.csv")
    parser.add_argument("--outdir", required=True, help="Output figure directory")
    parser.add_argument("--title", default="Mixed Residual Surface", help="Plot title prefix")
    args = parser.parse_args()

    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df_surface = pd.read_csv(args.surface)
    df_slice_m = pd.read_csv(args.slice_m)
    df_slice_g = pd.read_csv(args.slice_g)

    _plot_heatmap(
        df_surface,
        outdir / "residual_surface_heatmap.png",
        f"{args.title}: log10 ||F(m, gamma)||",
    )
    _plot_slice(
        df_slice_m,
        "m_fm",
        "log10_residual_norm",
        outdir / "residual_slice_mass.png",
        "mass m (fm^-1)",
        f"{args.title}: slice at fixed gamma",
    )
    _plot_slice(
        df_slice_g,
        "g_fm",
        "log10_residual_norm",
        outdir / "residual_slice_gamma.png",
        "gamma (fm^-1)",
        f"{args.title}: slice at fixed mass",
    )

    print(f"Wrote figures under: {outdir}")


if __name__ == "__main__":
    main()
