"""Audit and summarize the retained threshold algorithm experiments.

This reads completed runs and writes a new analysis directory. It never
recomputes a scattering kernel or changes a production result.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import statistics
import sys
from datetime import datetime, timezone


RUN_NAMES = (
    "ablation_v1", "endpoint_v1", "candidates_v1", "candidates_v2",
    "performance_v1", "performance_v2",
)


def sha256(path: Path) -> str:
    with path.open("rb") as stream:
        digest = hashlib.sha256()
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                    encoding="utf-8")


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        result = list(csv.DictReader(stream))
    if not result:
        raise ValueError(f"Empty result: {path}")
    return result


def write_rows(path: Path, data: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)


def audit_run(directory: Path, *, allow_interrupted: bool = False) -> dict:
    manifest = read_json(directory / "manifest.json")
    status = read_json(directory / "status.json")
    if status.get("complete") is not True and not (allow_interrupted and status.get("interrupted_for_continuation") is True):
        raise ValueError(f"Incomplete run: {directory}")
    if status.get("interrupted_for_continuation") is True:
        if sha256(directory / "rates.csv") != status.get("rates_sha256"):
            raise ValueError(f"Retained interrupted-run results changed: {directory}")
    for item in manifest["source_files"]:
        saved = directory / "source" / item["path"].replace("\\", "/")
        if sha256(saved) != item["sha256"]:
            raise ValueError(f"Archived source hash mismatch: {saved}")
    for item in manifest["input_files"]:
        saved = directory / item["snapshot"].replace("\\", "/")
        if sha256(saved) != item["sha256"]:
            raise ValueError(f"Archived input hash mismatch: {saved}")
    # Check historical snapshots against their own manifests, not current code.
    artifacts = [p for p in directory.iterdir() if p.is_file()]
    for name in ("knots", "tables"):
        if (directory / name).is_dir():
            artifacts.extend((directory / name).glob("*.csv"))
    return {
        "run": directory.name, "status": status,
        "manifest_sha256": sha256(directory / "manifest.json"),
        "source_files_verified": len(manifest["source_files"]),
        "input_files_verified": len(manifest["input_files"]),
        "input_hashes": [r["sha256"] for r in manifest["input_files"]],
        "branch": manifest["branch"], "head": manifest["head"],
        "artifacts": [{"path": str(p.relative_to(directory)), "sha256": sha256(p)}
                      for p in sorted(artifacts)],
    }


def candidate_key(row: dict) -> tuple:
    return (float(row["xi"]), row["channel"], row["loop"], row["representation"],
            float(row["rtol"]))


def summarize_performance(root: Path, name: str) -> list[dict]:
    reference = read_json(root / name / "reference.json")
    data = rows(root / name / "timings.csv")
    if abs(reference["rate384"] / reference["rate768"] - 1) >= 2e-5:
        raise ValueError("Dense thermal reference is insufficiently converged")
    result = []
    for method in ("dense", "adaptive_raw", "adaptive_x_sigma"):
        selected = [r for r in data if r["method"] == method]
        if len(selected) != 3 or len({r["rate"] for r in selected}) != 1:
            raise ValueError(f"Missing trials or numerical instability: {name}/{method}")
        for row in selected:
            rate_error = float(row["rate"]) / reference["rate768"] - 1
            if abs(rate_error) > 1e-4 or abs(float(row["J_relative_error"])) > 1e-4:
                raise ValueError(f"Matched-error target failed: {name}/{method}")
        entry = {"run": name, "method": method, "trials": len(selected),
                 "nodes": int(selected[0]["nodes"]),
                 "new_sigma_evaluations": int(selected[0]["new_sigma_evaluations"]),
                 "fallback_calls": int(selected[0]["fallback_calls"]),
                 "J_relative_error": float(selected[0]["J_relative_error"]),
                 "rate_relative_error": float(selected[0]["rate_relative_error"])}
        for field in ("build_seconds", "rate_seconds", "total_seconds", "query_seconds",
                      "build_cpu_seconds", "rate_cpu_seconds", "total_cpu_seconds",
                      "build_bytes", "rate_bytes"):
            if field not in selected[0]:
                continue
            values = [float(r[field]) for r in selected]
            for label, value in (("min", min(values)), ("median", statistics.median(values)),
                                 ("max", max(values))):
                entry[f"{field}_{label}"] = value
        result.append(entry)
    return result


def tau_comparison(ablation: list[dict], candidates: list[dict], kept: list[dict]) -> list[dict]:
    result = [{"xi": float(r["xi"]), "method": r["method"],
               "tau_fm": float(r["tau"]), "inverse_tau_inv_fm": float(r["inverse_tau"])}
              for r in kept]
    for xi in (0.36, 0.37, 0.38):
        source = {r["channel"]: r for r in ablation
                  if float(r["xi"]) == xi and r["method"] == "stock_on"}
        for loop, representation, label in (
            ("legacy", "raw_0.0001", "adaptive_raw"),
            ("legacy", "x_sigma_0.0001", "adaptive_x_sigma"),
            ("paired_pv64", "x_sigma_0.0001", "paired_pv64"),
        ):
            selected = [r for r in candidates if float(r["xi"]) == xi and
                        r["loop"] == loop and r["method"] == representation]
            if len(selected) != 2:
                raise ValueError("Two-channel partial tau reconstruction is incomplete")
            total = float(source[selected[0]["channel"]]["source_total"])
            for row in selected:
                old = source[row["channel"]]
                total += float(old["source_contribution"]) * (
                    float(row["rate"]) / float(old["source_rate"]) - 1)
            result.append({"xi": xi, "method": label, "tau_fm": 1 / total,
                           "inverse_tau_inv_fm": total})
    return sorted(result, key=lambda r: (r["method"], r["xi"]))


def make_figures(root: Path, out: Path, ablation: list[dict], tau: list[dict]) -> list[str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.dpi": 180})
    names = []
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.7), constrained_layout=True)
    methods = ("stock_on", "stock_off", "dense_on", "dense_off")
    for offset, channel, label in ((-0.18, "uubar_to_uubar", r"$u\bar u\to u\bar u$"),
                                   (0.18, "uubar_to_ddbar", r"$u\bar u\to d\bar d$")):
        values = [abs(float(next(r for r in ablation if float(r["xi"]) == 0.37 and
                                 r["method"] == m and r["channel"] == channel)
                            ["J_relative_error"])) for m in methods]
        axes[0].bar([i + offset for i in range(4)], values, width=0.34, label=label)
    axes[0].set(yscale="log", ylabel=r"Relative error in $J=\int_0^{0.6}x\sigma\,dx$",
                title=r"Local threshold integral, $\xi=0.37$")
    axes[0].set_xticks(range(4), ["Stock\non", "Stock\noff", "Dense\non", "Dense\noff"])
    axes[0].axhline(1e-4, color="0.4", ls=":", lw=1)
    axes[0].legend(loc="upper right", frameon=False)
    for method, label, marker in (("stock_on", "Stock + subtraction", "o"),
                                  ("stock_off", "Stock, subtraction off", "v"),
                                  ("dense_on", "Dense + subtraction", "s"),
                                  ("dense_off", "Dense, subtraction off", "^"),
                                  ("adaptive_x_sigma", r"Adaptive $x\sigma$", "x")):
        data = [r for r in tau if r["method"] == method]
        axes[1].plot([r["xi"] for r in data], [r["tau_fm"] for r in data],
                     marker=marker, ms=5, lw=1, label=label)
    axes[1].set(xlabel=r"$\xi$", ylabel=r"$\tau_u$ [fm]",
                title="Only two light channels replaced")
    axes[1].set_xticks([0.36, 0.37, 0.38])
    axes[1].legend(loc="lower right", fontsize=8, frameon=False)
    for extension in ("png", "pdf"):
        name = f"ablation_tau.{extension}"
        fig.savefig(out / name)
        names.append(name)
    plt.close(fig)

    inputs = rows(root / "ablation_v1" / "inputs" / "review_points.csv")
    state = next(r for r in inputs if r["review_mode"] == "mode_a" and
                 r["plot_panel"] == "muB0.0" and r["plot_series"] == "alpha1.0" and
                 float(r["xi"]) == 0.37)
    s0 = 4 * float(state["m_u"]) ** 2
    knots = rows(root / "ablation_v1" / "knots" / "xi0.37_uubar_to_uubar.csv")
    selected = [r for r in knots if 0 < float(r["s"]) - s0 <= 0.6 and
                math.isfinite(float(r["sigma"])) and float(r["sigma"]) > 0]
    selected.sort(key=lambda r: float(r["s"]))
    table = rows(root / "candidates_v2" / "tables" /
                 "xi0.37_uubar_to_uubar_x_sigma_0.0001.csv")
    paired = rows(root / "candidates_v2" / "tables" /
                  "xi0.37_uubar_to_uubar_paired_pv64.csv")
    paired_by_x = {float(r["x"]): float(r["sigma"]) for r in paired}
    endpoint_probe = [{"x_inv_fm2": float(r["x"]), "legacy_sigma_fm2": float(r["sigma"]),
                       "paired_sigma_fm2": paired_by_x[float(r["x"])]}
                      for r in table if float(r["x"]) < 1e-7 and float(r["x"]) in paired_by_x]
    write_rows(out / "endpoint_cross_section_probe.csv", endpoint_probe)
    scales = rows(root / "candidates_v2" / "threshold_scales.csv")
    scale = float(next(r for r in scales if float(r["xi"]) == 0.37 and
                       r["channel"] == "uubar_to_uubar" and r["family"] == "mixed")["x_delta"])
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.7), constrained_layout=True)
    for ax, weighted in zip(axes, (False, True)):
        x = [float(r["s"]) - s0 for r in selected]
        y = [float(r["sigma"]) * (xx if weighted else 1) for r, xx in zip(selected, x)]
        ax.loglog(x, y, lw=1, label="Legacy-loop direct samples")
        tx = [float(r["x"]) for r in table]
        ty = [float(r["sigma"]) * (xx if weighted else 1) for r, xx in zip(table, tx)]
        ax.loglog(tx, ty, ls="none", marker=".", ms=3, label=r"Adaptive $x\sigma$ knots")
        px = [float(r["x"]) for r in paired]
        py = [float(r["sigma"]) * (xx if weighted else 1) for r, xx in zip(paired, px)]
        ax.loglog(px, py, ls="--", lw=1, label="Paired PV, 64 nodes")
        ax.axvline(scale, color="0.5", ls=":", label=r"$x_\delta=\delta^2/c^2$")
        ax.set(xlabel=r"$x=s-4m_u^2$ [fm$^{-2}$]",
               ylabel=r"$x\sigma$" if weighted else r"$\sigma$ [fm$^2$]")
        ax.legend(loc="best", frameon=False, fontsize=8)
    for extension in ("png", "pdf"):
        name = f"weighted_threshold_shape.{extension}"
        fig.savefig(out / name)
        names.append(name)
    plt.close(fig)
    return names


def main(root: Path, out: Path) -> None:
    root, out = root.resolve(), out.resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError(f"Use a new empty output directory: {out}")
    audits = [audit_run(root / name) for name in RUN_NAMES]
    if len({tuple(r["input_hashes"]) for r in audits}) != 1:
        raise ValueError("Experiment inputs are not identical")
    ablation = rows(root / "ablation_v1" / "ablation.csv")
    direct = rows(root / "ablation_v1" / "direct_reference.csv")
    endpoints = rows(root / "endpoint_v1" / "endpoint_audit.csv")
    v1 = rows(root / "candidates_v1" / "adaptive_integrals.csv")
    v2 = rows(root / "candidates_v2" / "adaptive_integrals.csv")
    if any(abs(float(r["J_relative_error"])) > float(r["rtol"]) for r in v1 + v2):
        raise ValueError("An independently checked adaptive integral missed its target")
    previous = {candidate_key(r): r for r in v1}
    savings = []
    for row in v2:
        old = previous[candidate_key(row)]
        savings.append({"xi": row["xi"], "channel": row["channel"], "loop": row["loop"],
                        "representation": row["representation"], "rtol": row["rtol"],
                        "nodes": row["nodes"], "v1_evaluations": old["direct_evaluations"],
                        "v2_evaluations": row["direct_evaluations"],
                        "J_relative_error": row["J_relative_error"]})
    tight_names = sorted(p.name for p in (root / "candidates_v1" / "tables").glob("*.csv")
                         if p.name.endswith(("_0.0001.csv", "_paired_pv64.csv")))
    for name in tight_names:
        if sha256(root / "candidates_v1" / "tables" / name) != sha256(root / "candidates_v2" / "tables" / name):
            raise ValueError(f"Tight table changed; cannot reuse v1 full rates: {name}")
    parity = read_json(root / "candidates_v1" / "consumer_parity.json")
    if abs(parity["relative_error"]) > 1e-10:
        raise ValueError("Callable thermal consumer failed parity")
    reference = read_json(root / "performance_v2" / "reference.json")
    if reference["reused_sha256"] != sha256(root / "performance_v1" / "reference.json"):
        raise ValueError("Reused performance reference hash differs")
    rates = rows(root / "candidates_v1" / "candidate_rates.csv")
    sensitivity = []
    for row in rates:
        if row["loop"] != "paired_pv64":
            continue
        old = next(r for r in rates if (r["xi"], r["channel"], r["method"], r["loop"]) ==
                   (row["xi"], row["channel"], row["method"], "legacy"))
        sensitivity.append({"xi": row["xi"], "channel": row["channel"],
                            "legacy_rate_fm2": float(old["rate"]),
                            "paired_local_rate_fm2": float(row["rate"]),
                            "relative_change": float(row["rate"]) / float(old["rate"]) - 1})
    tau = tau_comparison(ablation, rates, rows(root / "ablation_v1" / "partial_tau.csv"))
    curvature = {}
    for method in sorted({r["method"] for r in tau}):
        t = {r["xi"]: r["tau_fm"] for r in tau if r["method"] == method}
        curvature[method] = t[0.37] / ((t[0.36] + t[0.38]) / 2) - 1
    finite_old = [r for r in endpoints if r["old_status"] == "finite" and
                  math.isfinite(float(r["old_absolute_error"]))]
    summary = {
        "status": "analysis_only", "numerical_defaults_promoted": False,
        "ablation_cases": len(direct), "ablation_comparisons": len(ablation),
        "ablation_center": [r for r in ablation if float(r["xi"]) == 0.37],
        "max_direct_order_relative_change": max(abs(float(r["relative_order_change"])) for r in direct),
        "max_direct_t_relative_change": max(abs(float(r["relative_t_change"])) for r in direct),
        "max_cache_quadrature_relative_change": max(abs(float(r["cache_quadrature_change"])) for r in ablation),
        "adaptive_comparisons_per_version": len(v2),
        "tight_tables_identical": len(tight_names),
        "consumer_parity": parity,
        "adaptive_max_abs_error_by_rtol": {str(tol): max(abs(float(r["J_relative_error"]))
            for r in v2 if float(r["rtol"]) == tol) for tol in (1e-3, 1e-4)},
        "endpoint": {
            "comparisons": len(endpoints),
            "old_failure_count": len(endpoints) - len(finite_old),
            "old_worst_finite": max(finite_old, key=lambda r: float(r["old_absolute_error"])),
            "max_candidate_error": {str(n): max(float(r[f"n{n}_error"]) for r in endpoints)
                                    for n in (32, 64, 128)},
            "max_oracle_estimated_error": max(float(r["oracle_error"]) for r in endpoints),
        },
        "max_paired_local_rate_relative_change": max(abs(r["relative_change"]) for r in sensitivity),
        "tau_center_relative_to_neighbor_mean": curvature,
        "performance_v1": summarize_performance(root, "performance_v1"),
        "performance_v2": summarize_performance(root, "performance_v2"),
    }
    out.mkdir(parents=True, exist_ok=True)
    write_json(out / "audit.json", {"runs": audits, "tight_table_names_identical": tight_names})
    write_json(out / "summary.json", summary)
    write_rows(out / "adaptive_probe_reuse.csv", savings)
    write_rows(out / "partial_tau_comparison.csv", tau)
    write_rows(out / "local_loop_sensitivity.csv", sensitivity)
    write_rows(out / "performance_summary.csv", summary["performance_v2"])
    figures = make_figures(root, out, ablation, tau)
    (out / "source").mkdir(exist_ok=True)
    shutil.copyfile(__file__, out / "source" / Path(__file__).name)
    outputs = [p for p in out.iterdir() if p.is_file()]
    write_json(out / "manifest.json", {
        "status": "analysis_only", "created_utc": datetime.now(timezone.utc).isoformat(),
        "run_root": str(root), "python": sys.version, "script_sha256": sha256(Path(__file__)),
        "input_manifests": [{"run": r["run"], "sha256": r["manifest_sha256"]} for r in audits],
        "outputs": [{"path": p.name, "sha256": sha256(p)} for p in sorted(outputs)],
        "figures": figures, "numerical_defaults_promoted": False,
    })
    write_json(out / "status.json", {"complete": True, "runs_audited": len(audits),
        "source_files_verified": sum(r["source_files_verified"] for r in audits),
        "input_files_verified": sum(r["input_files_verified"] for r in audits),
        "tight_tables_identical": len(tight_names), "all_matched_error_targets_passed": True})
    print(json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False))
    print(f"COMPLETE {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_root", type=Path)
    parser.add_argument("new_output_directory", type=Path)
    args = parser.parse_args()
    main(args.run_root, args.new_output_directory)
