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


def _load_trace(path: pathlib.Path):
    return pd.read_csv(path)


def _plot(m_vals, g_vals, n_grid, trace_df, out_png: pathlib.Path, title: str):
    fig, ax = plt.subplots(figsize=(8.2, 6.2))

    z = np.log10(np.maximum(n_grid, 1e-16))
    im = ax.imshow(
        z,
        origin="lower",
        aspect="auto",
        extent=[m_vals.min(), m_vals.max(), g_vals.min(), g_vals.max()],
        cmap="viridis",
    )
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("log10 ||F||")

    levels = np.linspace(np.nanmin(z), np.nanmax(z), 12)
    M, G = np.meshgrid(m_vals, g_vals)
    ax.contour(M, G, z, levels=levels, colors="white", linewidths=0.35, alpha=0.45)

    x = trace_df["m_fm"].to_numpy()
    y = trace_df["g_fm"].to_numpy()
    ax.plot(x, y, color="black", lw=1.2, alpha=0.9, label="solver trajectory")
    ax.plot([x[0]], [y[0]], marker="o", ms=5, color="lime", label="start")
    ax.plot([x[-1]], [y[-1]], marker="o", ms=5, color="red", label="end")

    ax.set_xlabel("mass m (fm^-1)")
    ax.set_ylabel("gamma (fm^-1)")
    ax.set_title(title)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description="Overlay NLsolve trajectory on residual surface")
    p.add_argument("--surface", required=True)
    p.add_argument("--trace", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--title", default="trajectory on residual surface")
    args = p.parse_args()

    m_vals, g_vals, n_grid = _load_surface(pathlib.Path(args.surface))
    trace = _load_trace(pathlib.Path(args.trace))
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    _plot(m_vals, g_vals, n_grid, trace, out, args.title)
    print(f"Wrote trajectory-overlay figure: {out}")


if __name__ == "__main__":
    main()
