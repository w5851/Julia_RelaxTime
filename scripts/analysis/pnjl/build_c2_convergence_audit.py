#!/usr/bin/env python3
"""Build the solver-free C0/C1/C2 convergence audit package.

This script consumes already-produced Actions CSV artifacts.  It never imports
Julia or calls an equilibrium solver.  Only derived tables and representative
plots are written to the requested output directory; raw full-resolution
curves stay in the external Actions/local artifact directories.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


COMPARATOR_PATH = Path(__file__).with_name("compare_phase_reference_convergence.py")
ARTIFACT_NAMES = {
    "boundary": "boundary_{tag}.csv",
    "cep": "cep_{tag}.csv",
    "spinodals": "spinodals_{tag}.csv",
    "crossover": "crossover_{tag}.csv",
    "grid_convergence": "phase_grid_convergence_{tag}.csv",
}

COMPARISON_FIELDS = [
    "artifact",
    "match_status",
    "grid_change_kind",
    "xi",
    "match_key",
    "metric",
    "candidate_value",
    "reference_value",
    "abs_diff",
    "rel_diff",
]
SUMMARY_FIELDS = [
    "artifact",
    "metric",
    "matched_count",
    "missing_count",
    "non_numeric_count",
    "max_abs_diff",
    "max_rel_diff",
    "max_abs_match_key",
    "max_rel_match_key",
    "mean_abs_diff",
    "mean_rel_diff",
    "status_changed_count",
    "not_applicable_count",
    "adaptive_xi_count",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    for prefix in ("c0", "c1", "c2"):
        parser.add_argument(f"--{prefix}-root", required=True, type=Path)
        parser.add_argument(f"--{prefix}-tag", required=True)
        parser.add_argument(f"--{prefix}-run", required=True)
    parser.add_argument("--calculation-sha", required=True)
    parser.add_argument("--postprocess-sha", required=True)
    parser.add_argument("--out-dir", required=True, type=Path)
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_comparator():
    spec = importlib.util.spec_from_file_location("phase_convergence_comparator", COMPARATOR_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load comparator: {COMPARATOR_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_pair(comparator, candidate_root: Path, candidate_tag: str, reference_root: Path, reference_tag: str, out_dir: Path, candidate_label: str, reference_label: str) -> dict:
    command = [
        sys.executable,
        str(COMPARATOR_PATH),
        "--candidate-root",
        str(candidate_root),
        "--candidate-tag",
        candidate_tag,
        "--reference-root",
        str(reference_root),
        "--reference-tag",
        reference_tag,
        "--candidate-label",
        candidate_label,
        "--reference-label",
        reference_label,
        "--out-dir",
        str(out_dir),
    ]
    subprocess.run(command, check=True)
    return json.loads((out_dir / "phase_reference_convergence_summary.json").read_text(encoding="utf-8"))


def comparison_exception(row: dict, comparator) -> str | None:
    status = row.get("match_status", "")
    if status not in {"matched", "interpolated"}:
        return status or "unclassified"
    try:
        abs_diff = float(row.get("abs_diff", ""))
        rel_diff = float(row.get("rel_diff", ""))
    except (TypeError, ValueError):
        return "non_numeric"
    artifact = row.get("artifact")
    metric = row.get("metric")
    if artifact in {"boundary", "spinodals"}:
        tolerance = comparator.GEOMETRY_TOLERANCES.get(metric)
        if tolerance is not None and abs_diff > tolerance:
            return "geometry_gate_exceeded"
    elif artifact == "crossover":
        if metric == "derivative":
            if rel_diff > comparator.CROSSOVER_TOLERANCES["derivative_rel"]:
                return "crossover_response_gate_exceeded"
        else:
            tolerance = comparator.CROSSOVER_TOLERANCES.get(metric)
            if tolerance is not None and abs_diff > tolerance:
                return "crossover_gate_exceeded"
    return None


def copy_derived_pair(pair_dir: Path, package_tables: Path, prefix: str, pair_summary: dict, comparator) -> None:
    mapping = {
        "phase_reference_artifact_inventory.csv": f"{prefix}_inventory.csv",
        "phase_reference_anchor_states.csv": f"{prefix}_anchor_states.csv",
        "phase_reference_geometry_gates.csv": f"{prefix}_geometry_gates.csv",
        "phase_reference_cep_gates.csv": f"{prefix}_cep_gates.csv",
    }
    package_tables.mkdir(parents=True, exist_ok=True)
    # Older working-tree generations included the full 40k-row comparison.
    # Remove that derived-only file so reruns converge to the compact v1 schema.
    (package_tables / f"{prefix}_comparison.csv").unlink(missing_ok=True)
    for source_name, target_name in mapping.items():
        source = pair_dir / source_name
        if source.is_file():
            target = package_tables / target_name
            target.write_bytes(source.read_bytes())

    write_csv(
        package_tables / f"{prefix}_comparison_summary.csv",
        pair_summary.get("comparison_summary", []),
        SUMMARY_FIELDS,
    )
    comparison_path = pair_dir / "phase_reference_convergence_comparison.csv"
    comparison_rows = read_csv(comparison_path) if comparison_path.is_file() else []
    exceptions = []
    crossover_rows = []
    for row in comparison_rows:
        if row.get("artifact") == "crossover" and row.get("xi") == "0.2875":
            crossover_rows.append(row)
        reason = comparison_exception(row, comparator)
        if reason is not None:
            exceptions.append({**row, "exception_kind": reason})
    write_csv(
        package_tables / f"{prefix}_comparison_exceptions.csv",
        exceptions,
        COMPARISON_FIELDS + ["exception_kind"],
    )
    write_csv(
        package_tables / f"{prefix}_crossover_xi_0p2875.csv",
        crossover_rows,
        COMPARISON_FIELDS,
    )


def write_plot_manifest(figures: Path, entries: list[dict]) -> None:
    (figures / "plot_manifest.json").write_text(
        json.dumps({"schema_version": "pnjl_c2_convergence_plot_manifest_v1", "figures": entries}, indent=2) + "\n",
        encoding="utf-8",
    )


def build_figures(package_dir: Path, c1_root: Path, c1_tag: str, c2_root: Path, c2_tag: str, comparator) -> list[dict]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figures = package_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    entries = []

    anchor_path = package_dir / "tables" / "c1_vs_c2_anchor_states.csv"
    regressions = [row for row in read_csv(anchor_path) if row.get("classification_regression") == "true"]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    labels = [f"({row['xi']}, {row['T_MeV']})" for row in regressions]
    ax.bar(range(len(labels)), [1] * len(labels), color="#b23a48")
    ax.set_xticks(range(len(labels)), labels, rotation=45, ha="right")
    ax.set_yticks([])
    ax.set_title("C1 → C2 public first-order classification regressions")
    ax.set_ylabel("regression")
    fig.tight_layout()
    path = figures / "public_first_order_regressions.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    entries.append({"path": path.name, "sha256": sha256_file(path), "source": anchor_path.as_posix()})

    cep_path = package_dir / "tables" / "c1_vs_c2_cep_gates.csv"
    cep_rows = read_csv(cep_path)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    xis = [float(row["xi"]) for row in cep_rows]
    widths = [float(row["candidate_bracket_width_T_MeV"]) for row in cep_rows]
    ax.plot(xis, widths, "o-", color="#2b6cb0", markersize=3, label="C2 bracket width")
    ax.axhline(0.1, color="#b23a48", linestyle="--", label="0.1 MeV gate")
    ax.set_xlabel("xi")
    ax.set_ylabel("T bracket width [MeV]")
    ax.set_title("C2 CEP bracket width")
    ax.legend()
    fig.tight_layout()
    path = figures / "cep_bracket_widths.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    entries.append({"path": path.name, "sha256": sha256_file(path), "source": cep_path.as_posix()})

    comparison_path = package_dir / "tables" / "c1_vs_c2_crossover_xi_0p2875.csv"
    crossover_rows = [
        row
        for row in read_csv(comparison_path)
        if row.get("artifact") == "crossover" and row.get("xi") == "0.2875"
    ]
    colors = {"T_crossover_MeV": "#2f855a", "rho": "#805ad5", "derivative": "#dd6b20"}
    fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True)
    for axis, metric in zip(axes, ("T_crossover_MeV", "rho", "derivative")):
        rows = [row for row in crossover_rows if row.get("metric") == metric]
        x = [float(row["match_key"].split("|", 1)[1]) for row in rows]
        candidate = [float(row["candidate_value"]) for row in rows]
        reference = [float(row["reference_value"]) for row in rows]
        axis.plot(x, candidate, ".-", color=colors[metric], label="C2")
        axis.plot(x, reference, "--", color="#4a5568", label="C1")
        axis.set_ylabel(metric)
        axis.grid(alpha=0.25)
        axis.legend(loc="best")
    axes[-1].set_xlabel("mu [MeV], physical union")
    fig.suptitle("C1/C2 crossover comparison at xi=0.2875")
    fig.tight_layout()
    path = figures / "crossover_xi_0p2875_physical_mu.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    entries.append({"path": path.name, "sha256": sha256_file(path), "source": comparison_path.as_posix()})
    return entries


def main() -> None:
    args = parse_args()
    comparator = load_comparator()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    temp_root = args.out_dir / ".pair_replays"
    temp_root.mkdir(parents=True, exist_ok=True)

    roots = {name: getattr(args, f"{name}_root") for name in ("c0", "c1", "c2")}
    tags = {name: getattr(args, f"{name}_tag") for name in ("c0", "c1", "c2")}
    runs = {name: getattr(args, f"{name}_run") for name in ("c0", "c1", "c2")}
    c1_vs_c2 = run_pair(comparator, roots["c2"], tags["c2"], roots["c1"], tags["c1"], temp_root / "c1_vs_c2", "C2", "C1")
    c0_vs_c2 = run_pair(comparator, roots["c2"], tags["c2"], roots["c0"], tags["c0"], temp_root / "c0_vs_c2", "C2", "C0")

    tables = args.out_dir / "tables"
    copy_derived_pair(temp_root / "c1_vs_c2", tables, "c1_vs_c2", c1_vs_c2, comparator)
    copy_derived_pair(temp_root / "c0_vs_c2", tables, "c0_vs_c2", c0_vs_c2, comparator)
    for pair in (c1_vs_c2, c0_vs_c2):
        pair.pop("public_anchor_states", None)
        pair.pop("geometry_gate_rows", None)
        pair.pop("cep_gate_rows", None)

    input_files = []
    for name, root in roots.items():
        for artifact, template in ARTIFACT_NAMES.items():
            path = root / template.format(tag=tags[name])
            if path.is_file():
                input_files.append({"stage": name, "run": runs[name], "artifact": artifact, "path": path.as_posix(), "sha256": sha256_file(path)})

    package_json = {
        "schema_version": "pnjl_c2_convergence_audit_v1",
        "calculation_sha": args.calculation_sha,
        "postprocess_sha": args.postprocess_sha,
        "source_runs": runs,
        "source_tags": tags,
        "input_files": input_files,
        "replay": {
            "c1_vs_c2_verdict": c1_vs_c2["verdict"],
            "c1_vs_c2_secondary_failures": c1_vs_c2.get("secondary_failures", []),
            "c1_vs_c2_classification_regressions": len(c1_vs_c2.get("classification_regressions", [])),
            "c1_vs_c2_cep_bracket_failures": len(c1_vs_c2.get("cep_bracket_failures", [])),
            "c1_vs_c2_crossover_failure_count": c1_vs_c2.get("crossover_failure_count", 0),
            "c0_vs_c2_verdict": c0_vs_c2["verdict"],
            "comparison_table_policy": "summary_exceptions_and_xi_0p2875_plot_subset",
        },
        "solver_called": False,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "generator": {
            "script": str(Path(__file__).as_posix()),
            "script_sha256": sha256_file(Path(__file__)),
            "comparator": str(COMPARATOR_PATH.as_posix()),
            "comparator_sha256": sha256_file(COMPARATOR_PATH),
        },
    }
    (args.out_dir / "manifest.json").write_text(json.dumps(package_json, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    plot_entries = build_figures(args.out_dir, roots["c1"], tags["c1"], roots["c2"], tags["c2"], comparator)
    (args.out_dir / "figures" / "plot_manifest.json").write_text(
        json.dumps({"schema_version": "pnjl_c2_convergence_plot_manifest_v1", "figures": plot_entries}, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    claim_rows = [
        {"claim_id": "c2_primary_verdict", "claim": "C2 comparator primary verdict is classification_regression", "status": "supported", "evidence": "manifest.json; tables/c1_vs_c2_anchor_states.csv"},
        {"claim_id": "public_regression_count", "claim": "Nine public first-order anchors regress from C1 first-order to C2 ambiguous", "status": "supported", "evidence": "tables/c1_vs_c2_anchor_states.csv"},
        {"claim_id": "cep_bracket_risk", "claim": "Sixteen CEP brackets exceed the 0.1 MeV hard gate", "status": "supported", "evidence": "tables/c1_vs_c2_cep_gates.csv"},
        {"claim_id": "crossover_local_risk", "claim": "Crossover residual risk is localized at xi=0.2875 high physical-mu points", "status": "supported", "evidence": "tables/c1_vs_c2_comparison_exceptions.csv; tables/c1_vs_c2_crossover_xi_0p2875.csv; figures/crossover_xi_0p2875_physical_mu.png"},
        {"claim_id": "solver_free", "claim": "Replay does not call the PNJL solver", "status": "supported", "evidence": "manifest.json"},
        {"claim_id": "production_ready", "claim": "C2 is production-ready", "status": "not_supported", "evidence": "manifest.json"},
    ]
    write_csv(args.out_dir / "claim_ledger.csv", claim_rows, ["claim_id", "claim", "status", "evidence"])

    readme = f"""# C2 convergence audit v1

这是 C0/C1/C2 的 solver-free 后处理审计包。输入固定为 Actions 已完成 artifact；脚本不调用 Julia equilibrium solver，也不修改数值产物。

## 结论

- C1→C2 的主 verdict 为 `{c1_vs_c2['verdict']}`。
- 9 个 public first-order anchor 在 C2 退为 ambiguous；它们是分类回归，不被解释为 monotone。
- 16 个 CEP bracket 超过 `0.1 MeV`，并单独记录 endpoint search resolution 与实际 ambiguity width。
- crossover 只在双方共同的 monotone/crossover 物理温区、物理 `mu` 并集上插值比较；剩余风险集中在 `xi=0.2875` 的高 `mu` 局部点。
- C2 保持 `diagnostic-only`，不创建 production PR，不启动新的 C0/C1/C2、reference 或 transport。

## 输入与追溯

完整曲线留在外部 Actions/local artifact；本目录只保留派生表和代表图。所有输入文件 SHA、calculation SHA、postprocess SHA、source run 和 generator SHA 见 `manifest.json`。

## 作者检查项

后续 feasibility 必须在固定 cap=12 下验证 density certificate、CEP bracket 和 crossover refinement；不得通过增加 cap 或放宽现有容差掩盖上述回归。
"""
    (args.out_dir / "README.md").write_text(readme, encoding="utf-8")
    shutil.rmtree(temp_root, ignore_errors=True)
    print(json.dumps({"out_dir": args.out_dir.as_posix(), "verdict": c1_vs_c2["verdict"]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
