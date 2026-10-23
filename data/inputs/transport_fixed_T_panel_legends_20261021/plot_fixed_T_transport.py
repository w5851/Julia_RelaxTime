"""Fixed-temperature (mode B) transport comparison for the review revision.

Reads the publication_clean_v5 display layer of Julia_RelaxTime (v2 on-shell kernel,
prod_v2) without modifying it. Usage:
    python plot_fixed_T_transport.py <publication_clean_points.csv> <scan.csv> <out_dir>
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator, FuncFormatter, LogLocator, NullFormatter

points_csv, scan_csv, out_dir = map(Path, sys.argv[1:4])
out_dir.mkdir(parents=True, exist_ok=True)

matplotlib.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "TeX Gyre Termes", "Nimbus Roman"],
    "mathtext.fontset": "stix",
    "font.size": 13, "axes.labelsize": 13, "xtick.labelsize": 13,
    "ytick.labelsize": 13, "legend.fontsize": 12,
    "axes.linewidth": 0.7, "lines.linewidth": 1.2,
    "pdf.fonttype": 42,
})
COLORS = ["#0072B2", "#D55E00", "#009E73"]
STYLES = ["-", "--", "-."]
MUB = [0.0, 450.0, 900.0]
TEMPS = [120.0, 160.0, 200.0]
ROWS = [("eta_over_s", r"$\eta/s$", "log"),
        ("zeta_over_s", r"$\zeta/s$", "log"),
        ("sigma_over_T", r"$\sigma/T$", "linear")]

pts = pd.read_csv(points_csv)
pts = pts[pts.mode_key == "mode_b"]
scan = pd.read_csv(scan_csv, comment="#")


def branch_split(T, mu):
    # phase_curr also flips across crossovers; only a finite mass gap marks a first-order branch switch.
    s = scan[(scan.T_MeV == T) & (scan.muB_MeV == mu)].sort_values("xi")
    gap = np.abs(np.diff(s.m_u.values)) * 197.327
    first_order = s.phase_structure.values[1:] == "first_order_possible"
    jumps = np.flatnonzero((gap > 100.0) & first_order)
    return None if len(jumps) == 0 else 0.5 * (s.xi.values[jumps[0]] + s.xi.values[jumps[0] + 1])


fig, axes = plt.subplots(3, 3, figsize=(6.75, 6.9), sharex=True)
fig.subplots_adjust(left=0.115, right=0.985, bottom=0.08, top=0.875, hspace=0.24, wspace=0.36)
letters = iter("abcdefghi")
records = []
for r, (obs, label, scale) in enumerate(ROWS):
    for c, T in enumerate(TEMPS):
        ax = axes[r, c]
        for k, mu in enumerate(MUB):
            x = pts[(pts.plot_panel == f"T{T}") & (pts.plot_series == f"muB{mu}") & (pts.observable == obs)]
            x = x.sort_values("xi")
            xi, y = x.xi.values, x.clean_value.values
            cut = branch_split(T, mu)
            segs = [(xi, y)] if cut is None else [(xi[xi < cut], y[xi < cut]), (xi[xi > cut], y[xi > cut])]
            for xs, ys in segs:
                ax.plot(xs, ys, color=COLORS[k], ls=STYLES[k])
            if cut is not None:
                lo, hi = xi[xi < cut][-1], xi[xi > cut][0]
                ax.plot(lo, y[xi == lo][0], "o", mfc="white", mec=COLORS[k], ms=5.5, mew=1.0)
                ax.plot(hi, y[xi == hi][0], "s", mfc="white", mec=COLORS[k], ms=5.5, mew=1.0)
                records.append((obs, T, mu, lo, hi))
        ax.set_yscale(scale)
        if scale == "linear":
            ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        else:
            lo, hi = ax.get_ylim()
            subs = (1.0,) if np.log10(hi / lo) > 1.5 else (1.0, 2.0, 5.0)
            ax.yaxis.set_major_locator(LogLocator(base=10, subs=subs))
            ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
            ax.yaxis.set_minor_formatter(NullFormatter())
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.tick_params(which="both", direction="in", top=True, right=True)
        ax.set_xlim(-0.52, 0.52)
        ax.set_xticks([-0.5, 0, 0.5])
        ax.set_xticklabels(["-0.5", "0", "0.5"])
        title = f"({next(letters)})"
        if r == 0:
            title += rf"  $T={T:.0f}$ MeV"
        ax.text(0.0, 1.03, title, transform=ax.transAxes, ha="left", va="bottom", fontsize=12.5)
        if c == 0:
            ax.set_ylabel(label)
        if r == 2:
            ax.set_xlabel(r"$\xi$")

handles = [Line2D([], [], color=COLORS[k], ls=STYLES[k], label=rf"$\mu_B={mu:.0f}$ MeV") for k, mu in enumerate(MUB)]
handles += [Line2D([], [], ls="none", marker="o", mfc="white", mec="k", label="restored"),
            Line2D([], [], ls="none", marker="s", mfc="white", mec="k", label="broken")]
fig.legend(handles=handles, loc="upper center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 0.995),
           handlelength=1.8, columnspacing=1.0, handletextpad=0.4)

for ext in ("pdf", "png"):
    fig.savefig(out_dir / f"figure4_transport_fixed_T.{ext}", dpi=300)
pd.DataFrame(records, columns=["observable", "T_MeV", "muB_MeV", "xi_restored_end", "xi_broken_end"]).to_csv(
    out_dir / "first_order_endpoints.csv", index=False)
print("done", records)
