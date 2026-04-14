#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import subprocess
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt


AUX_BEGIN = "<!-- AUX_FIGURES:BEGIN -->"
AUX_END = "<!-- AUX_FIGURES:END -->"


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Generate PNJL Mott xi auxiliary figures")
    ap.add_argument("--derived-csv", type=Path, required=True)
    ap.add_argument("--scan-csv", type=Path, required=True)
    ap.add_argument("--gap-csv", type=Path, required=True)
    ap.add_argument("--mode-ab-dir", type=Path, required=True)
    ap.add_argument("--fig-dir", type=Path, required=True)
    ap.add_argument("--doc", type=Path, required=True)
    ap.add_argument("--skip-fig5", action="store_true")
    return ap.parse_args()


def require_file(path: Path, name: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"missing {name}: {path}")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as f:
        lines = [ln for ln in f if ln.strip() and not ln.lstrip().startswith("#")]
    return list(csv.DictReader(lines))


def require_columns(rows: list[dict[str, str]], required: list[str], name: str) -> None:
    if not rows:
        raise ValueError(f"empty csv: {name}")
    missing = [c for c in required if c not in rows[0]]
    if missing:
        raise ValueError(f"missing columns in {name}: {missing}")


def group_by_xi(rows: list[dict[str, str]]) -> dict[float, list[dict[str, str]]]:
    out: dict[float, list[dict[str, str]]] = {}
    for r in rows:
        xi = float(r["xi"])
        out.setdefault(xi, []).append(r)
    for xi in out:
        out[xi].sort(key=lambda rr: float(rr["T_MeV"]))
    return out


def estimate_crossing(ts: list[float], deltas: list[float]) -> float:
    for i in range(len(ts) - 1):
        d0 = deltas[i]
        d1 = deltas[i + 1]
        if d0 == 0.0:
            return ts[i]
        if d0 * d1 < 0.0:
            t0 = ts[i]
            t1 = ts[i + 1]
            return t0 + (0.0 - d0) * (t1 - t0) / (d1 - d0)
    raise ValueError("no crossing found")


def linear_fit(x: list[float], y: list[float]) -> tuple[float, float]:
    n = float(len(x))
    sx = sum(x)
    sy = sum(y)
    sxx = sum(v * v for v in x)
    sxy = sum(vx * vy for vx, vy in zip(x, y))
    den = n * sxx - sx * sx
    if abs(den) < 1e-14:
        return 0.0, sy / n
    a = (n * sxy - sx * sy) / den
    b = (sy - a * sx) / n
    return a, b


def solve_3x3(a: list[list[float]], b: list[float]) -> tuple[float, float, float]:
    m = [row[:] + [rhs] for row, rhs in zip(a, b)]
    n = 3
    for i in range(n):
        piv = i
        for r in range(i + 1, n):
            if abs(m[r][i]) > abs(m[piv][i]):
                piv = r
        if abs(m[piv][i]) < 1e-14:
            return 0.0, 0.0, 0.0
        if piv != i:
            m[i], m[piv] = m[piv], m[i]
        v = m[i][i]
        for c in range(i, n + 1):
            m[i][c] /= v
        for r in range(n):
            if r == i:
                continue
            f = m[r][i]
            for c in range(i, n + 1):
                m[r][c] -= f * m[i][c]
    return m[0][3], m[1][3], m[2][3]


def quadratic_fit(x: list[float], y: list[float]) -> tuple[float, float, float]:
    sx = sum(x)
    sx2 = sum(v * v for v in x)
    sx3 = sum(v * v * v for v in x)
    sx4 = sum(v * v * v * v for v in x)
    sy = sum(y)
    sxy = sum(vx * vy for vx, vy in zip(x, y))
    sx2y = sum((vx * vx) * vy for vx, vy in zip(x, y))
    a = [
        [sx4, sx3, sx2],
        [sx3, sx2, sx],
        [sx2, sx, float(len(x))],
    ]
    b = [sx2y, sxy, sy]
    return solve_3x3(a, b)


def lin_interp(x0: float, y0: float, x1: float, y1: float, x: float) -> float:
    if abs(x1 - x0) < 1e-14:
        return y0
    return y0 + (x - x0) * (y1 - y0) / (x1 - x0)


def interp_at_t(rows: list[dict[str, str]], col: str, t_mev: float) -> float:
    pts = sorted(((float(r["T_MeV"]), float(r[col])) for r in rows), key=lambda z: z[0])
    if t_mev <= pts[0][0]:
        return pts[0][1]
    if t_mev >= pts[-1][0]:
        return pts[-1][1]
    for i in range(len(pts) - 1):
        t0, y0 = pts[i]
        t1, y1 = pts[i + 1]
        if t0 <= t_mev <= t1:
            return lin_interp(t0, y0, t1, y1, t_mev)
    return pts[-1][1]


def interp_gap_value(gap_by_xi: dict[float, list[dict[str, str]]], xi: float, col: str, t_mev: float) -> float:
    if xi in gap_by_xi:
        return interp_at_t(gap_by_xi[xi], col, t_mev)

    xis = sorted(gap_by_xi.keys())
    if xi < xis[0] or xi > xis[-1]:
        raise ValueError(f"xi out of interpolation range for gap data: {xi}")

    lo = None
    hi = None
    for i in range(len(xis) - 1):
        x0 = xis[i]
        x1 = xis[i + 1]
        if x0 <= xi <= x1:
            lo = x0
            hi = x1
            break
    if lo is None or hi is None:
        raise ValueError(f"cannot bracket xi for gap interpolation: {xi}")

    y0 = interp_at_t(gap_by_xi[lo], col, t_mev)
    y1 = interp_at_t(gap_by_xi[hi], col, t_mev)
    return lin_interp(lo, y0, hi, y1, xi)


def plot_fig1_tmott_vs_xi(out_png: Path, xi: list[float], tmott_pi: list[float], tmott_k: list[float]) -> None:
    fig, ax = plt.subplots(figsize=(6.75, 4.6))
    ax.scatter(xi, tmott_pi, label="pi samples", marker="o", color="#4477AA")
    ax.scatter(xi, tmott_k, label="K samples", marker="s", color="#EE6677")

    p1_pi = linear_fit(xi, tmott_pi)
    p2_pi = quadratic_fit(xi, tmott_pi)
    p1_k = linear_fit(xi, tmott_k)
    p2_k = quadratic_fit(xi, tmott_k)

    x0 = min(xi)
    x1 = max(xi)
    xx = [x0 + (x1 - x0) * i / 200.0 for i in range(201)]
    ax.plot(xx, [p1_pi[0] * v + p1_pi[1] for v in xx], "--", color="#4477AA", label="pi linear")
    ax.plot(xx, [p2_pi[0] * v * v + p2_pi[1] * v + p2_pi[2] for v in xx], "-", color="#4477AA", label="pi quadratic")
    ax.plot(xx, [p1_k[0] * v + p1_k[1] for v in xx], "--", color="#EE6677", label="K linear")
    ax.plot(xx, [p2_k[0] * v * v + p2_k[1] * v + p2_k[2] for v in xx], "-", color="#EE6677", label="K quadratic")

    ax.set_xlabel("xi")
    ax.set_ylabel("T_mott [MeV]")
    ax.legend(frameon=False)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_fig4_orderparam(
    out_png: Path,
    xi: list[float],
    fixed_phi: list[float],
    traj_phi: list[float],
    fixed_mu: list[float],
    traj_mu: list[float],
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(6.75, 3.4), sharex=True)
    axes[0].plot(xi, fixed_phi, "o-", label="Phi @ T=200", color="#4477AA")
    axes[0].plot(xi, traj_phi, "s--", label="Phi @ T_mott(xi)", color="#EE6677")
    axes[0].set_ylabel("Phi")
    axes[0].legend(frameon=False)

    axes[1].plot(xi, fixed_mu, "o-", label="m_u @ T=200", color="#228833")
    axes[1].plot(xi, traj_mu, "s--", label="m_u @ T_mott(xi)", color="#CCBB44")
    axes[1].set_ylabel("m_u [fm^-1]")
    axes[1].legend(frameon=False)

    for ax in axes:
        ax.set_xlabel("xi")

    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_fig2_gamma_delta(out_png: Path, derived_by_xi: dict[float, list[dict[str, str]]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.8), sharex=True)
    channels = [
        ("pi", "Gamma_pi", "M_u_plus_M_d", "M_pi"),
        ("K", "Gamma_K", "M_u_plus_M_s", "M_K"),
    ]
    palette = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE"]

    for ax, (ch, gcol, thcol, mcol) in zip(axes, channels):
        ax2 = ax.twinx()
        for i, xi in enumerate(sorted(derived_by_xi.keys())):
            rows = derived_by_xi[xi]
            xs = [float(r["T_MeV"]) for r in rows]
            gamma = [float(r[gcol]) for r in rows]
            delta = [float(r[thcol]) - float(r[mcol]) for r in rows]
            color = palette[i % len(palette)]
            ax.plot(xs, gamma, color=color, linestyle="-", label=f"xi={xi:+.2f} Gamma")
            ax2.plot(xs, delta, color=color, linestyle="--", label=f"xi={xi:+.2f} Delta")

        ax.set_title(ch)
        ax.set_xlabel("T [MeV]")
        ax.set_ylabel("Gamma [fm^-1]")
        ax2.set_ylabel("Delta = M_thr - M_mes [fm^-1]")

        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, frameon=False, fontsize=7, loc="best")

    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def annotate_mode_ab(src: Path, dst: Path, text_lines: list[str]) -> None:
    img = mpimg.imread(src)
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.imshow(img)
    ax.axis("off")
    text = "\n".join(text_lines)
    ax.text(
        0.02,
        0.98,
        text,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox=dict(boxstyle="round", fc="white", ec="black", alpha=0.75),
    )
    fig.savefig(dst, dpi=300, bbox_inches="tight")
    plt.close(fig)


def can_build_fig5(scan_rows: list[dict[str, str]], gap_rows: list[dict[str, str]]) -> tuple[bool, str]:
    need_scan = {"xi", "T_MeV", "M_pi", "M_K", "m_u", "m_s"}
    need_gap = {"xi", "T_MeV", "Phi", "Phibar", "m_u", "m_s"}
    have_scan = set(scan_rows[0].keys()) if scan_rows else set()
    have_gap = set(gap_rows[0].keys()) if gap_rows else set()
    miss = [f"scan:{c}" for c in sorted(need_scan - have_scan)] + [f"gap:{c}" for c in sorted(need_gap - have_gap)]
    if miss:
        return False, f"missing required fields for fig5: {miss}"
    return True, "ok"


def _run_polarization_exporter(derived_csv: Path, scan_csv: Path, gap_csv: Path, out_csv: Path) -> None:
    script = Path("scripts") / "analysis" / "export_mott_polarization_decomposition.jl"
    cmd = [
        "julia",
        "--project=.",
        str(script),
        "--derived-csv",
        str(derived_csv),
        "--scan-csv",
        str(scan_csv),
        "--gap-csv",
        str(gap_csv),
        "--out-csv",
        str(out_csv),
    ]
    subprocess.run(cmd, check=True)


def _read_fig5_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def plot_fig5_from_decomposition_csv(out_png: Path, decomp_csv: Path) -> None:
    rows = _read_fig5_rows(decomp_csv)
    grouped: dict[str, list[dict[str, str]]] = {"pi": [], "K": []}
    for r in rows:
        ch = r["channel"]
        if ch in grouped:
            grouped[ch].append(r)

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.2), sharex=True)
    channels = [("pi", axes[0]), ("K", axes[1])]
    for ch, ax in channels:
        rs = sorted(grouped[ch], key=lambda rr: float(rr["xi"]))
        x = [float(r["xi"]) for r in rs]

        re_full = [float(r["re_full"]) for r in rs]
        re_ind = [float(r["re_indirect"]) for r in rs]
        re_dir = [float(r["re_direct"]) for r in rs]
        im_full = [float(r["im_full"]) for r in rs]
        im_ind = [float(r["im_indirect"]) for r in rs]
        im_dir = [float(r["im_direct"]) for r in rs]

        ax.plot(x, re_full, "o-", color="#4477AA", label="Re full")
        ax.plot(x, re_ind, "s--", color="#228833", label="Re indirect")
        ax.plot(x, re_dir, "d-.", color="#66CCEE", label="Re direct")
        ax.plot(x, im_full, "o-", color="#EE6677", label="Im full")
        ax.plot(x, im_ind, "s--", color="#AA3377", label="Im indirect")
        ax.plot(x, im_dir, "d-.", color="#CCBB44", label="Im direct")

        ax.set_title(ch)
        ax.set_xlabel("xi")
        ax.set_ylabel("polarization")
        ax.legend(frameon=False, fontsize=7)

    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def upsert_aux_section(doc: Path, *, fig5_generated: bool, fig5_msg: str) -> None:
    section = [
        "## 8. 辅助图集（新增）",
        "",
        "- 图1：`figures/fig1_tmott_vs_xi_fit.png`",
        "- 图4：`figures/fig4_orderparam_direct_indirect.png`",
        "- 图2：`figures/fig2_gamma_delta_dualaxis.png`",
        "- 图3（轻增强复用）：`mode_ab/mott_mode_ab__M_K__xi3_annotated.png` 与 `mode_ab/mott_mode_ab__M_pi__xi3_annotated.png`",
        f"- 图5：`figures/fig5_taylor_decomposition.png` ({'generated' if fig5_generated else 'skipped'})，定义为 ReΠ/ImΠ 的 xi 直接-间接分解示意（双面板：pi/K）。",
        "",
        "符号统一：`Delta = M_thr - M_mes`（若旧段落使用相反定义，请按符号翻转理解）。",
    ]
    if not fig5_generated and fig5_msg:
        section.append(f"fig5 skipped: {fig5_msg}")

    body = doc.read_text(encoding="utf-8")
    block = AUX_BEGIN + "\n" + "\n".join(section) + "\n" + AUX_END
    if AUX_BEGIN in body and AUX_END in body:
        pre = body.split(AUX_BEGIN)[0]
        post = body.split(AUX_END, 1)[1]
        out = pre.rstrip() + "\n\n" + block + "\n" + post.lstrip("\n")
    else:
        out = body.rstrip() + "\n\n" + block + "\n"
    doc.write_text(out, encoding="utf-8")


def main() -> None:
    args = parse_args()
    args.fig_dir.mkdir(parents=True, exist_ok=True)

    require_file(args.derived_csv, "derived-csv")
    require_file(args.scan_csv, "scan-csv")
    require_file(args.gap_csv, "gap-csv")
    require_file(args.doc, "doc")
    require_file(args.mode_ab_dir / "mott_mode_ab__M_K__xi3.png", "mode_ab K png")
    require_file(args.mode_ab_dir / "mott_mode_ab__M_pi__xi3.png", "mode_ab pi png")

    derived_rows = read_rows(args.derived_csv)
    scan_rows = read_rows(args.scan_csv)
    gap_rows = read_rows(args.gap_csv)

    require_columns(
        derived_rows,
        ["xi", "T_MeV", "M_pi", "M_K", "M_u_plus_M_d", "M_u_plus_M_s", "Gamma_pi", "Gamma_K"],
        "derived",
    )
    require_columns(scan_rows, ["xi", "T_MeV", "M_pi", "M_K", "Gamma_pi", "Gamma_K", "m_u", "m_s"], "scan")
    require_columns(gap_rows, ["xi", "T_MeV", "Phi", "m_u", "m_s"], "gap")

    derived_by_xi = group_by_xi(derived_rows)
    scan_by_xi = group_by_xi(scan_rows)
    gap_by_xi = group_by_xi(gap_rows)

    xis = sorted(derived_by_xi.keys())

    tmott_pi = []
    tmott_k = []
    for xi in xis:
        rows = derived_by_xi[xi]
        ts = [float(r["T_MeV"]) for r in rows]
        d_pi = [float(r["M_pi"]) - float(r["M_u_plus_M_d"]) for r in rows]
        d_k = [float(r["M_K"]) - float(r["M_u_plus_M_s"]) for r in rows]
        tmott_pi.append(estimate_crossing(ts, d_pi))
        tmott_k.append(estimate_crossing(ts, d_k))

    plot_fig1_tmott_vs_xi(args.fig_dir / "fig1_tmott_vs_xi_fit.png", xis, tmott_pi, tmott_k)

    phi_fixed = []
    phi_traj = []
    mu_fixed = []
    mu_traj = []
    for xi, tpi in zip(xis, tmott_pi):
        phi_fixed.append(interp_gap_value(gap_by_xi, xi, "Phi", 200.0))
        phi_traj.append(interp_gap_value(gap_by_xi, xi, "Phi", tpi))
        mu_fixed.append(interp_gap_value(gap_by_xi, xi, "m_u", 200.0))
        mu_traj.append(interp_gap_value(gap_by_xi, xi, "m_u", tpi))

    plot_fig4_orderparam(
        args.fig_dir / "fig4_orderparam_direct_indirect.png",
        xis,
        phi_fixed,
        phi_traj,
        mu_fixed,
        mu_traj,
    )

    plot_fig2_gamma_delta(args.fig_dir / "fig2_gamma_delta_dualaxis.png", derived_by_xi)

    dt_k = tmott_k[-1] - tmott_k[0]
    dt_pi = tmott_pi[-1] - tmott_pi[0]
    annotate_mode_ab(
        args.mode_ab_dir / "mott_mode_ab__M_K__xi3.png",
        args.mode_ab_dir / "mott_mode_ab__M_K__xi3_annotated.png",
        ["Crossing shifts right with xi", f"Delta T_mott^K ~= {dt_k:.2f} MeV"],
    )
    annotate_mode_ab(
        args.mode_ab_dir / "mott_mode_ab__M_pi__xi3.png",
        args.mode_ab_dir / "mott_mode_ab__M_pi__xi3_annotated.png",
        ["Crossing shifts right with xi", f"Delta T_mott^pi ~= {dt_pi:.2f} MeV"],
    )

    fig5_generated = False
    fig5_msg = ""
    fig5_ok, fig5_msg = can_build_fig5(scan_rows, gap_rows)
    if (not args.skip_fig5) and fig5_ok:
        fig5_data_csv = args.fig_dir / "fig5_polarization_decomposition.csv"
        _run_polarization_exporter(args.derived_csv, args.scan_csv, args.gap_csv, fig5_data_csv)
        plot_fig5_from_decomposition_csv(args.fig_dir / "fig5_taylor_decomposition.png", fig5_data_csv)
        fig5_generated = True

    upsert_aux_section(args.doc, fig5_generated=fig5_generated, fig5_msg=fig5_msg)

    print(f"wrote auxiliary figures under: {args.fig_dir}")


if __name__ == "__main__":
    main()
