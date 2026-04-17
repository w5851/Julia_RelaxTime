#!/usr/bin/env python3

import argparse
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _load_surface(path: pathlib.Path):
    df = pd.read_csv(path)
    m_vals = np.sort(df["m_fm"].unique())
    g_vals = np.sort(df["g_fm"].unique())
    n_grid = df.pivot(index="g_fm", columns="m_fm", values="residual_norm").values
    return m_vals, g_vals, n_grid


def _plot(surface_path, scatter_path, out_png, title):
    m_vals, g_vals, n_grid = _load_surface(pathlib.Path(surface_path))
    ds = pd.read_csv(scatter_path)

    fig, ax = plt.subplots(figsize=(8.4, 6.2))

    z = np.log10(np.maximum(n_grid, 1e-16))
    im = ax.imshow(
        z,
        origin="lower",
        aspect="auto",
        extent=[m_vals.min(), m_vals.max(), g_vals.min(), g_vals.max()],
        cmap="Greys",
        alpha=0.35,
    )
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("log10 ||F|| (surface background)")

    # Overlay direct-evaluation iso-classes
    n = ds["residual_norm"].to_numpy()
    m = ds["m_fm"].to_numpy()
    g = ds["g_fm"].to_numpy()

    mask_1e2 = n <= 1e-2
    mask_1e3 = n <= 1e-3
    mask_1e4 = n <= 1e-4

    if np.any(mask_1e2):
        ax.scatter(m[mask_1e2], g[mask_1e2], s=2, c="#f59e0b", alpha=0.35, label="direct ||F|| <= 1e-2")
    if np.any(mask_1e3):
        ax.scatter(m[mask_1e3], g[mask_1e3], s=3, c="#ef4444", alpha=0.45, label="direct ||F|| <= 1e-3")
    if np.any(mask_1e4):
        ax.scatter(m[mask_1e4], g[mask_1e4], s=5, c="#7c3aed", alpha=0.65, label="direct ||F|| <= 1e-4")

    ax.set_xlabel("mass m (fm^-1)")
    ax.set_ylabel("gamma (fm^-1)")
    ax.set_title(title)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description="Overlay direct-eval low-residual points on residual surface")
    p.add_argument("--surface", required=True)
    p.add_argument("--scatter", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--title", default="direct vs surface residual overlay")
    args = p.parse_args()

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    _plot(args.surface, args.scatter, out, args.title)
    print(f"Wrote direct-vs-surface overlay: {out}")


if __name__ == "__main__":
    main()
