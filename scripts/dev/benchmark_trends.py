#!/usr/bin/env python3
"""Combine PNJL reports and compare like-for-like runs. Standard library only."""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SNAPSHOT = "benchmark_snapshot.json"
SCHEMA = "pnjl_benchmark_snapshot_v1"
ENVIRONMENT_KEYS = (
    "julia_version", "julia_threads", "blas_threads", "blas_config", "os", "arch",
    "cpu_name", "cpu_model", "cpu_threads", "runner_os", "runner_arch", "runner_image",
    "pnjl_profile", "physics_profile",
    "project_sha256", "manifest_sha256", "config_sha256", "workload_sha256",
)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True,
                               allow_nan=False) + "\n", encoding="utf-8")


def timestamp(value: str) -> datetime:
    result = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if result.tzinfo is None:
        raise ValueError("benchmark timestamps must include a timezone")
    return result.astimezone(timezone.utc)


def number(value, name: str, *, positive: bool = False, maximum=None) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number")
    if value < 0 or (positive and value == 0) or (maximum is not None and value > maximum):
        raise ValueError(f"{name} is outside its valid range")
    return value


def build_snapshot(directory: Path, threshold_status: str) -> dict:
    single = read_json(directory / "single_point_benchmark.json")
    scan = read_json(directory / "scan_benchmark.json")
    for report in (single, scan):
        if report["schema_version"] != "pnjl_benchmark_raw_v1":
            raise ValueError("raw benchmark report has an unsupported schema")
        number(report["samples"], "requested samples", positive=True)
        if not report["benchmarks"]:
            raise ValueError("benchmark report is empty")
    if single["metadata"] != scan["metadata"]:
        raise ValueError("single-point and scan reports came from different environments or sources")
    metadata = single["metadata"]
    environment = metadata["environment"]
    if any(environment.get(key) in (None, "") for key in ENVIRONMENT_KEYS):
        raise ValueError("benchmark environment metadata is incomplete")
    metrics, samples = {}, {}

    def add(suite, label, metric, value, unit, direction="lower"):
        key = f"{suite} / {label} / {metric}"
        if key in metrics:
            raise ValueError(f"duplicate benchmark metric: {key}")
        metrics[key] = dict(value=number(value, key, positive=unit.startswith("ms")),
                            unit=unit, direction=direction)

    for suite, report in (("single", single), ("scan", scan)):
        for entry in report["benchmarks"]:
            label, stats = entry["label"], entry["stats"]
            samples[f"{suite} / {label}"] = number(stats["samples"], "actual samples", positive=True)
            if suite == "single":
                add(suite, label, "median", stats["median"], "ms")
                add(suite, label, "minimum", stats["min"], "ms")
                add(suite, label, "memory", stats["memory_bytes"], "bytes")
            else:
                number(entry["n_points"], "scan points", positive=True)
                add(suite, label, "median", stats["per_point_median_ms"], "ms/point")
                add(suite, label, "memory", stats["memory_mb"] * 1024**2, "bytes")
                rate = number(entry["convergence_rate"], "convergence", maximum=1)
                add(suite, label, "convergence", 100 * rate, "%", "higher")
            add(suite, label, "allocations", stats["allocs"], "count")
    provenance = metadata["provenance"]
    config = {
        "single": {"samples": single["samples"], "config": single["config"],
                   "methodology": single["methodology"],
                   "solvers": {entry["label"]: entry["params"] for entry in single["benchmarks"]}},
        "scan": {"samples": scan["samples"], "tmu": scan["tmu_config"], "trho": scan["trho_config"],
                 "methodology": scan["methodology"],
                 "points": {entry["label"]: entry["n_points"] for entry in scan["benchmarks"]}},
    }
    return dict(
        schema_version=SCHEMA,
        generated_at=max(timestamp(single["generated_at"]), timestamp(scan["generated_at"])).isoformat(),
        provenance=provenance, environment=environment, config=config, metrics=metrics, samples=samples,
        threshold_status=threshold_status,
        baseline_eligible=(threshold_status == "success" and not provenance["source_dirty"]
                           and provenance["branch"] == "main" and bool(provenance["run_id"])
                           and provenance["event"] in ("push", "schedule", "workflow_dispatch")),
    )


def compare(current: dict, baseline: dict, max_age_days: int = 35) -> dict:
    try:
        if baseline["schema_version"] != SCHEMA or baseline["threshold_status"] != "success":
            raise ValueError("baseline schema or absolute-threshold status is invalid")
        age = (timestamp(current["generated_at"]) - timestamp(baseline["generated_at"])).total_seconds() / 86400
        if age < 0:
            raise ValueError("baseline is newer than the current measurement")
        if age > max_age_days:
            return dict(status="expired_baseline", detail=f"基线已有 {age:.1f} 天，超过 {max_age_days} 天。", rows=[])
        mismatches = [key for key in ENVIRONMENT_KEYS
                      if baseline["environment"].get(key) != current["environment"][key]]
        if baseline["config"] != current["config"]:
            mismatches.append("benchmark_config")
        if set(baseline["metrics"]) != set(current["metrics"]):
            mismatches.append("metric_set")
        if mismatches:
            return dict(status="incompatible_baseline", detail="配置或环境不同：" + ", ".join(mismatches), rows=[])
        rows = []
        for key, metric in sorted(current["metrics"].items()):
            old = baseline["metrics"][key]
            before = number(old["value"], key, positive=metric["unit"].startswith("ms"),
                            maximum=100 if metric["unit"] == "%" else None)
            if old["unit"] != metric["unit"] or old["direction"] != metric["direction"]:
                raise ValueError(f"baseline metric contract differs: {key}")
            now = metric["value"]
            if metric["unit"] == "%":
                delta, delta_unit = now - before, "pp"
                warning = delta < -2
            else:
                delta = 100 * (now / before - 1) if before else (0.0 if now == 0 else None)
                delta_unit = "%"
                warning = (delta is None and now > 0) or (delta is not None and delta >= 25)
            rows.append(dict(metric=key, current=now, baseline=before, unit=metric["unit"],
                             delta=delta, delta_unit=delta_unit, warning=warning))
        return dict(status="compared", detail="同配置比较；耗时/内存/分配次数增长 ≥25% 或收敛率降低 >2 个百分点仅作预警。",
                    baseline_provenance=baseline["provenance"], age_days=age, rows=rows,
                    warning_count=sum(row["warning"] for row in rows))
    except (KeyError, TypeError, ValueError, AttributeError) as error:
        return dict(status="invalid_baseline", detail=str(error), rows=[])


def gh(*args: str) -> str:
    result = subprocess.run(["gh", *args], capture_output=True, text=True, encoding="utf-8", timeout=45)
    if result.returncode:
        raise RuntimeError(f"GitHub read failed (exit {result.returncode})")
    return result.stdout


def github_baseline(current: dict, repository: str, max_age_days: int) -> tuple[dict | None, dict]:
    """Only completed successful main runs can supply a baseline; PR artifacts never can."""
    rejected = []
    try:
        runs = json.loads(gh("api", f"repos/{repository}/actions/workflows/pnjl-benchmarks.yml/runs"
                             "?branch=main&status=success&per_page=20"))["workflow_runs"]
        for run in runs:
            if (run["head_branch"] != "main" or run["conclusion"] != "success"
                    or run["event"] not in ("push", "schedule", "workflow_dispatch")
                    or str(run["id"]) == current["provenance"]["run_id"]
                    or timestamp(run["created_at"]) > timestamp(current["generated_at"])):
                continue
            result = dict(status="missing_baseline", detail="该运行没有新版 benchmark snapshot。", rows=[])
            artifacts = json.loads(gh("api", f"repos/{repository}/actions/runs/{run['id']}/artifacts"
                                      "?per_page=100"))["artifacts"]
            candidates = [a for a in artifacts if a["name"] == "pnjl-benchmark"]
            if candidates and all(a["expired"] for a in candidates):
                result = dict(status="expired_baseline", detail="历史 benchmark artifact 已过期。", rows=[])
            elif candidates:
                with tempfile.TemporaryDirectory(prefix="pnjl-baseline-") as temp:
                    gh("run", "download", str(run["id"]), "--repo", repository,
                       "--name", "pnjl-benchmark", "--dir", temp)
                    path = Path(temp) / SNAPSHOT
                    if path.is_file():
                        try:
                            baseline = read_json(path)
                            if (not baseline.get("baseline_eligible")
                                    or baseline["provenance"]["commit"] != run["head_sha"]
                                    or baseline["provenance"]["run_id"] != str(run["id"])
                                    or baseline["provenance"]["repository"] != repository
                                    or baseline["provenance"]["branch"] != "main"):
                                raise ValueError("baseline provenance does not match the successful main run")
                            result = compare(current, baseline, max_age_days)
                            if result["status"] == "compared":
                                result.update(baseline_url=run["html_url"], rejected_candidates=rejected)
                                return baseline, result
                        except (KeyError, TypeError, ValueError) as error:
                            result = dict(status="invalid_baseline", detail=str(error), rows=[])
            rejected.append(dict(run_id=run["id"], status=result["status"], detail=result["detail"]))
        if rejected:
            first = rejected[0]
            return None, dict(status=first["status"], detail=first["detail"],
                              rejected_candidates=rejected, rows=[])
        return None, dict(status="missing_baseline", detail="尚无可用的 main 基线；本次结果可作为后续比较起点。", rows=[])
    except (OSError, RuntimeError, subprocess.TimeoutExpired, KeyError, ValueError) as error:
        return None, dict(status="baseline_unavailable", detail=f"无法读取 GitHub 历史基线：{error}",
                          rejected_candidates=rejected, rows=[])


def render_summary(current: dict, comparison: dict) -> str:
    provenance, environment = current["provenance"], current["environment"]
    lines = ["# PNJL benchmark", "",
             f"状态：**{comparison['status']}**。{comparison['detail']}",
             f"Commit：{provenance['commit']}；branch：{provenance['branch']}；source dirty：{provenance['source_dirty']}。",
             f"Julia {environment['julia_version']}；Julia/BLAS threads："
             f"{environment['julia_threads']}/{environment['blas_threads']}；"
             f"{environment['os']} {environment['arch']}；CPU：{environment['cpu_model']}。",
             f"绝对阈值：{current['threshold_status']}；可作为 main 基线：{current['baseline_eligible']}。"]
    if comparison.get("baseline_url"):
        lines.append(f"基线运行：{comparison['baseline_url']}")
    if comparison.get("baseline_provenance"):
        lines.append(f"基线 commit：{comparison['baseline_provenance']['commit']}。")
    lines += ["", "| 指标 | 当前 | 基线 | 变化 | 状态 |",
              "| --- | ---: | ---: | ---: | --- |"]
    rows = {row["metric"]: row for row in comparison["rows"]}
    for key, metric in sorted(current["metrics"].items()):
        label = key.replace("|", r"\|").replace("\n", " ")
        row = rows.get(key)
        old, delta, state = "—", "—", "仅记录"
        if row:
            old = f"{row['baseline']:.6g} {row['unit']}"
            delta = "从 0 增加" if row["delta"] is None else f"{row['delta']:+.2f}{row['delta_unit']}"
            state = "预警" if row["warning"] else "正常"
        lines.append(f"| {label} | {metric['value']:.6g} {metric['unit']} | {old} | {delta} | {state} |")
    lines += ["", "实际采样数：" + "；".join(f"{key}={value}" for key, value in sorted(current["samples"].items())),
              "", "原始数据、环境、配置和比较 JSON 随 pnjl-benchmark artifact 保留 90 天。"]
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=Path(os.environ.get(
        "PNJL_BENCHMARK_OUTPUT_DIR", ROOT / "tests/perf/results/pnjl")))
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--baseline", type=Path)
    group.add_argument("--github-baseline", action="store_true")
    parser.add_argument("--repository", default=os.environ.get("GITHUB_REPOSITORY", ""))
    parser.add_argument("--max-baseline-age-days", type=int, default=35)
    parser.add_argument("--threshold-status", choices=("success", "failure", "skipped", "not_checked"),
                        default="not_checked")
    args = parser.parse_args(argv)
    if args.max_baseline_age_days <= 0:
        parser.error("--max-baseline-age-days must be positive")
    if args.github_baseline and not args.repository:
        parser.error("--github-baseline requires --repository or GITHUB_REPOSITORY")
    args.results_dir.mkdir(parents=True, exist_ok=True)
    try:
        current = build_snapshot(args.results_dir, args.threshold_status)
        baseline = None
        if args.github_baseline:
            baseline, comparison = github_baseline(current, args.repository, args.max_baseline_age_days)
        elif args.baseline:
            try:
                baseline = read_json(args.baseline)
                comparison = compare(current, baseline, args.max_baseline_age_days)
            except (OSError, ValueError) as error:
                comparison = dict(status="invalid_baseline", detail=str(error), rows=[])
        else:
            comparison = dict(status="missing_baseline", detail="未指定基线，仅记录本次结果。", rows=[])
        write_json(args.results_dir / SNAPSHOT, current)
        write_json(args.results_dir / "environment.json", current["environment"])
        if baseline is not None and comparison["status"] == "compared":
            write_json(args.results_dir / "baseline.json", baseline)
        summary = render_summary(current, comparison)
        exit_code = 0
    except (OSError, KeyError, ValueError, TypeError, AttributeError) as error:
        comparison = dict(status="invalid_current", detail=str(error), rows=[])
        summary = f"# PNJL benchmark\n\n状态：**invalid_current**。\n\n{error}\n"
        exit_code = 1
    write_json(args.results_dir / "comparison.json", comparison)
    (args.results_dir / "comparison.md").write_text(summary, encoding="utf-8")
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with Path(os.environ["GITHUB_STEP_SUMMARY"]).open("a", encoding="utf-8") as stream:
            stream.write(summary)
    print(f"[benchmark-trends] {comparison['status']}: {comparison['detail']}")
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
