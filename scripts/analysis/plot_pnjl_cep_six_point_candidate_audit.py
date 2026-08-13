#!/usr/bin/env python3
"""Render a temporary local rho-mu audit plot for the six-point replay."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--curve", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    points = sorted(read_csv(args.curve), key=lambda row: float(row["rho"]))
    candidates = [
        row for row in read_csv(args.candidates)
        if row["source"] == "full_fine_oracle"
        and abs(float(row["xi"]) + 0.34375) < 1e-8
        and abs(float(row["T_MeV"]) - 142.1875) < 1e-8
    ]
    x = [float(row["rho"]) for row in points]
    y = [float(row["muq_MeV"]) for row in points]

    # Discrete extrema are shown as observations from the sampled curve.  We
    # do not interpolate them into a production certificate.
    turns: list[int] = []
    slopes = [y[index + 1] - y[index] for index in range(len(y) - 1)]
    for index in range(1, len(slopes)):
        if slopes[index - 1] * slopes[index] < 0.0:
            turns.append(index)
    extrema = [(x[index], y[index], "spinodal max" if slopes[index - 1] > 0 else "spinodal min")
               for index in turns]

    colors = {1: "#007c91", 2: "#d1495b"}
    styles = {1: "-", 2: "--"}
    labels = {1: "candidate 1: sign-change root", 2: "candidate 2: tolerance grid-hit"}

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2), constrained_layout=True,
                             gridspec_kw={"width_ratios": [1.05, 1.0, 1.15]})
    panels = [
        (axes[0], (1.80, 2.50), (237.0868, 237.0900), "context around the weak S"),
        (axes[1], (2.10, 2.27), (237.08835, 237.08925), "local S and spinodals"),
        (axes[2], (2.135, 2.235), (237.08878, 237.08902), "candidate separation"),
    ]
    for axis, xlim, ylim, title in panels:
        mask = [xlim[0] <= value <= xlim[1] for value in x]
        axis.plot([value for value, keep in zip(x, mask) if keep],
                  [value for value, keep in zip(y, mask) if keep],
                  color="#264653", linewidth=1.35, label="full-fine rho-mu curve")
        for candidate in candidates:
            number = int(candidate["candidate"])
            mu = float(candidate["mu_MeV"])
            axis.axhline(mu, color=colors[number], linestyle=styles[number], linewidth=1.1,
                         label=labels[number])
            crossings = json.loads(candidate["crossings"])
            crossing_y = [mu] * len(crossings)
            axis.scatter(crossings, crossing_y, color=colors[number], s=28, zorder=4,
                         edgecolor="white", linewidth=0.45)
            if axis is axes[2]:
                for crossing in crossings:
                    axis.annotate(f"{crossing:.5f}", (crossing, mu),
                                  xytext=(0, 7 if number == 1 else -15),
                                  textcoords="offset points", ha="center", fontsize=7,
                                  color=colors[number], rotation=90)
        for rho, mu, kind in extrema:
            if xlim[0] <= rho <= xlim[1] and ylim[0] <= mu <= ylim[1]:
                axis.scatter([rho], [mu], marker="x", s=46, color="#6a4c93", zorder=5)
                axis.annotate(kind, (rho, mu), xytext=(4, 5), textcoords="offset points",
                              fontsize=7, color="#6a4c93")
        axis.set_xlim(*xlim)
        axis.set_ylim(*ylim)
        axis.set_xlabel(r"$\rho$")
        axis.set_ylabel(r"$\mu_q$ [MeV]")
        axis.set_title(title)
        axis.grid(alpha=0.22)

    axes[0].legend(fontsize=7, loc="best")
    fig.suptitle(r"$\xi=-0.34375$, $T=142.1875$ MeV: full-fine Maxwell candidate audit",
                 fontsize=13)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=260)
    plt.close(fig)
    print(json.dumps({"output": str(args.output), "candidate_count": len(candidates),
                      "extrema": extrema}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
