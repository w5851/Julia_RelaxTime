#!/usr/bin/env python3

import argparse
import pathlib

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


def _plot_overlay(df: pd.DataFrame, out_png: pathlib.Path, title: str):
    m_vals, g_vals, re_grid, im_grid, n_grid = _load_grid(df)
    m_mesh, g_mesh = np.meshgrid(m_vals, g_vals)

    fig, ax = plt.subplots(figsize=(8, 6))

    n_plot = np.log10(np.maximum(n_grid, 1e-16))
    bg = ax.imshow(
        n_plot,
        origin="lower",
        aspect="auto",
        extent=[m_vals.min(), m_vals.max(), g_vals.min(), g_vals.max()],
        cmap="Greys",
        alpha=0.35,
    )
    cbar = fig.colorbar(bg, ax=ax)
    cbar.set_label("log10 ||F|| (background)")

    cs_re = ax.contour(m_mesh, g_mesh, re_grid, levels=[0.0], colors=["crimson"], linewidths=2.0)
    cs_im = ax.contour(m_mesh, g_mesh, im_grid, levels=[0.0], colors=["royalblue"], linewidths=2.0)

    # Mark grid minimum of residual norm as a visual reference
    idx = np.unravel_index(np.argmin(n_grid), n_grid.shape)
    g0 = g_vals[idx[0]]
    m0 = m_vals[idx[1]]
    n0 = n_grid[idx]
    ax.plot([m0], [g0], marker="o", color="black", ms=5, label=f"grid min ||F||={n0:.3e}")

    ax.set_xlabel("mass m (fm^-1)")
    ax.set_ylabel("gamma (fm^-1)")
    ax.set_title(title)

    handles = [
        plt.Line2D([0], [0], color="crimson", lw=2, label="Re(F)=0 contour"),
        plt.Line2D([0], [0], color="royalblue", lw=2, label="Im(F)=0 contour"),
        plt.Line2D([0], [0], marker="o", color="black", lw=0, markersize=5, label=f"grid min ||F||={n0:.3e}"),
    ]
    ax.legend(handles=handles, loc="best")
    fig.tight_layout()
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Plot Re=0 and Im=0 contour overlay")
    parser.add_argument("--surface", required=True, help="Path to residual_surface.csv")
    parser.add_argument("--out", required=True, help="Output PNG path")
    parser.add_argument("--title", default="Zero-Contour Overlay", help="Figure title")
    args = parser.parse_args()

    df = pd.read_csv(args.surface)
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    _plot_overlay(df, out, args.title)
    print(f"Wrote overlay: {out}")


if __name__ == "__main__":
    main()
