"""Summarize finite-x resonance diagnosis and original-cache runtime costs."""

from __future__ import annotations

import argparse
import bisect
from datetime import datetime, timezone
from pathlib import Path
import shutil
import statistics

from analyze_threshold_cache_upgrade import audit_run, read_json, rows, sha256, write_json, write_rows


def main(root: Path, out: Path) -> None:
    root, out = root.resolve(), out.resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError(f"Use a new empty output directory: {out}")
    audits = [audit_run(root / name) for name in ("tail_feature_v1", "baseline_cost_v2")]
    feature = rows(root / "tail_feature_v1" / "tail_feature.csv")
    roots = rows(root / "tail_feature_v1" / "denominator_roots.csv")
    for before, after in zip(feature, feature[1:]):
        if float(after["no_s_mixed_scalar"]) > float(before["no_s_mixed_scalar"]) + 1e-12:
            raise ValueError("Removing the s mixed scalar did not remove the nonmonotone peak")
    maximum = max(feature, key=lambda r: float(r["sigma"]))
    minimum = min(feature[:feature.index(maximum)], key=lambda r: float(r["sigma"]))
    scalar_root = next(r for r in roots if r["family"] == "scalar_mixed")
    state = next(r for r in rows(root / "ablation_v1" / "inputs" / "review_points.csv")
                 if r["review_mode"] == "mode_a" and r["plot_panel"] == "muB0.0" and
                 r["plot_series"] == "alpha1.0" and float(r["xi"]) == 0.37)
    s0 = 4 * float(state["m_u"]) ** 2
    stock_x = sorted(float(r["s"]) - s0 for r in rows(root / "ablation_v1" / "knots" /
                     "xi0.37_uubar_to_uubar.csv") if r["stock_node"] == "true" and float(r["s"]) > s0)
    scales = rows(root / "candidates_v2" / "threshold_scales.csv")
    x_delta = float(next(r for r in scales if float(r["xi"]) == 0.37 and
                         r["channel"] == "uubar_to_uubar" and r["family"] == "mixed")["x_delta"])
    i = bisect.bisect_left(stock_x, x_delta)
    timings = rows(root / "baseline_cost_v2" / "timings.csv")
    cost = []
    for method in ("stock", "dense", "adaptive_x_sigma"):
        selected = [r for r in timings if r["method"] == method]
        if len(selected) != 3 or len({r["rate"] for r in selected}) != 1:
            raise ValueError(f"Missing or numerically inconsistent timing repetitions: {method}")
        expected = "false" if method == "stock" else "true"
        if any(r["meets_target"] != expected for r in selected):
            raise ValueError(f"Unexpected accuracy classification: {method}")
        fields = ("build_seconds", "build_cpu_seconds", "callable_rate_seconds",
                  "callable_rate_cpu_seconds", "original_rate_seconds", "original_rate_cpu_seconds",
                  "total_seconds", "total_cpu_seconds", "callable_total_seconds", "callable_total_cpu_seconds")
        record = {"method": method, "nodes": int(selected[0]["nodes"]),
                  "new_sigma_evaluations": int(selected[0]["new_sigma_evaluations"]),
                  "J_relative_error": float(selected[0]["J_relative_error"]),
                  "rate_relative_error": float(selected[0]["rate_relative_error"])}
        for field in fields:
            values = [float(r[field]) for r in selected if r[field]]
            for name, value in (("min", min(values) if values else None),
                                ("median", statistics.median(values) if values else None),
                                ("max", max(values) if values else None)):
                record[f"{field}_{name}"] = value
        cost.append(record)
    summary = {
        "status": "analysis_only", "numerical_defaults_promoted": False,
        "sampling_gap": {"x_delta_inv_fm2": x_delta, "stock_left_inv_fm2": stock_x[i-1],
                         "stock_right_inv_fm2": stock_x[i]},
        "finite_x_feature": {
            "local_maximum": maximum, "preceding_minimum": minimum,
            "real_axis_resonance": scalar_root,
            "removing_s_mixed_scalar_is_monotone": True,
            "verification": read_json(root / "tail_feature_v1" / "status.json"),
            "interpretation": "finite-width s-channel mixed-scalar resonance in this model; not a complex-pole mass determination",
        },
        "cost": cost,
        "cost_scope": "one incoming/outgoing scattering process at p128/angle20/phi36; includes cache construction and one thermal average",
        "original_consumer_note": "Stock native consumer is timed first after construction; callable stock timing separately checks equal-formula parity.",
    }
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / "audit.json", {"runs": audits,
        "excluded_failed_run": read_json(root / "baseline_cost_v1" / "status.json")})
    write_json(out / "summary.json", summary)
    write_rows(out / "cost_summary.csv", cost)
    make_plot(feature, scalar_root, out)
    (out / "source").mkdir(exist_ok=True)
    for script in (Path(__file__), Path(__file__).with_name("analyze_threshold_cache_upgrade.py")):
        shutil.copyfile(script, out / "source" / script.name)
    write_json(out / "manifest.json", {
        "status": "analysis_only", "created_utc": datetime.now(timezone.utc).isoformat(),
        "run_root": str(root), "script_sha256": sha256(Path(__file__)),
        "helper_sha256": sha256(Path(__file__).with_name("analyze_threshold_cache_upgrade.py")),
        "input_manifests": [{"run": r["run"], "sha256": r["manifest_sha256"]} for r in audits],
        "outputs": [{"path": p.name, "sha256": sha256(p)} for p in sorted(out.iterdir()) if p.is_file()],
        "numerical_defaults_promoted": False,
    })
    write_json(out / "status.json", {"complete": True, "runs_audited": len(audits),
        "source_files_verified": sum(r["source_files_verified"] for r in audits),
        "input_files_verified": sum(r["input_files_verified"] for r in audits)})
    print(f"COMPLETE {out}")
    print("Sampling gap:", summary["sampling_gap"])
    print("Finite-x maximum:", maximum["x"], "real-axis root:", scalar_root)
    for r in cost:
        print(r["method"], "build", r["build_seconds_median"],
              "native thermal", r["original_rate_seconds_median"],
              "callable thermal", r["callable_rate_seconds_median"],
              "total", r["total_seconds_median"], "cpu", r["total_cpu_seconds_median"])


def make_plot(data: list[dict], root: dict, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    x = [float(r["x"]) for r in data]
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.7), constrained_layout=True)
    for key, label, style in (("sigma", "Full cross section", "-"),
                               ("no_s_mixed_scalar", "s mixed scalar removed", "--"),
                               ("s_scalar", "s scalar squared contribution", ":")):
        ax[0].plot(x, [float(r[key]) for r in data], ls=style, lw=1.6, label=label)
    ax[0].set(xlabel=r"$x=s-4m_u^2$ [fm$^{-2}$]", ylabel=r"$\sigma$ or contribution [fm$^2$]")
    ax[0].legend(frameon=False, loc="lower right", fontsize=8)
    for key, label in (("scalar_mixed_re", r"$\mathrm{Re}\,H_S$"),
                       ("scalar_mixed_im", r"$\mathrm{Im}\,H_S$")):
        ax[1].plot(x, [float(r[key]) for r in data], label=label)
    ax[1].axhline(0, color="0.6", lw=0.8)
    ax[1].axvline(float(root["x"]), color="0.5", ls=":", label=r"$\mathrm{Re}\,H_S=0$")
    ax[1].set(xlabel=r"$x=s-4m_u^2$ [fm$^{-2}$]", ylabel=r"$H_S=\det M_S/\det K_-$")
    ax[1].legend(frameon=False, fontsize=9)
    for suffix in ("png", "pdf"):
        fig.savefig(out / f"finite_x_resonance.{suffix}", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_root", type=Path)
    parser.add_argument("new_output_directory", type=Path)
    args = parser.parse_args()
    main(args.run_root, args.new_output_directory)
