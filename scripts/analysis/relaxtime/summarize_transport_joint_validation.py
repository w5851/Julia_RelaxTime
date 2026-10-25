"""Collect completed joint qualification shards without hiding missing points."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import statistics
from pathlib import Path


OBSERVABLES = {f"tau_{name}" for name in ("u", "d", "s", "ubar", "dbar", "sbar")}
OBSERVABLES.update(("eta_over_s", "zeta_over_s", "sigma_over_T"))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def numerical_sources(manifest: dict) -> dict[str, str]:
    """Prevent combining measurements made with different numerical code."""
    return {item["path"].replace("\\", "/"): item["sha256"] for item in manifest["source_files"]
            if item["path"].replace("\\", "/").startswith(("src/", "config/models/", "config/physics/"))
            or item["path"] in ("Project.toml", "Manifest.toml")}


def summarize(runs: list[Path], outdir: Path, expected_points: int, benchmark: Path | None) -> dict:
    if expected_points < 1:
        raise ValueError("expected_points must be positive")
    if outdir.exists() and any(outdir.iterdir()):
        raise ValueError(f"Use a new empty output directory: {outdir}")
    comparisons: dict[tuple[str, str], dict] = {}
    sources = []
    failures = []
    configurations = None
    input_hash = None
    source_hashes = None
    for directory in runs:
        manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
        status = json.loads((directory / "status.json").read_text(encoding="utf-8"))
        if not status.get("complete"):
            raise ValueError(f"Unfinished joint run: {directory}")
        settings = manifest["parameters"]["configurations"]
        inputs = digest(directory / "inputs/review_points.csv")
        numerical = numerical_sources(manifest)
        if configurations is None:
            configurations, input_hash = settings, inputs
            source_hashes = numerical
        elif settings != configurations or inputs != input_hash:
            raise ValueError("Joint runs have different numerical configurations or inputs")
        elif numerical != source_hashes:
            raise ValueError("Joint runs have different numerical source versions")
        input_points = {row["point_id"] for row in read_rows(directory / "inputs/review_points.csv")}
        if not status.get("joint_target_passed"):
            failures.append({"run": str(directory.resolve()), "status": status})
        path = directory / "comparisons.csv"
        rows = read_rows(path) if path.exists() else []
        if len(rows) != status.get("comparisons"):
            raise ValueError(f"Comparison count differs from status: {directory}")
        for row in rows:
            key = row["point_id"], row["observable"]
            if key[0] not in input_points:
                raise ValueError(f"Joint point is absent from the recorded inputs: {key[0]}")
            if key in comparisons:
                raise ValueError(f"Duplicate comparison: {key}")
            candidate, reference = float(row["candidate"]), float(row["reference"])
            if not all(math.isfinite(v) and v > 0 for v in (candidate, reference)):
                raise ValueError(f"Invalid observable: {key}")
            change = candidate / reference - 1
            if not math.isclose(change, float(row["relative_change"]), rel_tol=1e-12, abs_tol=1e-14):
                raise ValueError(f"Recorded relative change is inconsistent: {key}")
            passed = abs(change) <= 0.01
            if passed != (row["passed"].lower() == "true"):
                raise ValueError(f"Recorded qualification flag is inconsistent: {key}")
            comparisons[key] = dict(point_id=key[0], observable=key[1], candidate=candidate,
                                    reference=reference, relative_change=change, passed=passed)
        sources.append(dict(directory=str(directory.resolve()), manifest_sha256=digest(directory / "manifest.json"),
                            status_sha256=digest(directory / "status.json"),
                            comparisons_sha256=digest(path) if path.exists() else None))
    points = sorted({point for point, _ in comparisons})
    missing = {point: sorted(OBSERVABLES - {obs for p, obs in comparisons if p == point}) for point in points}
    unexpected = {obs for _, obs in comparisons} - OBSERVABLES
    if unexpected or any(missing.values()):
        raise ValueError(f"Incomplete or unexpected observables: {missing}, {unexpected}")
    maxima = {}
    for name in sorted(OBSERVABLES):
        rows = [row for (_, obs), row in comparisons.items() if obs == name]
        if rows:
            worst = max(rows, key=lambda row: abs(row["relative_change"]))
            maxima[name] = dict(point_id=worst["point_id"], relative_change=worst["relative_change"])
    all_passed = (len(points) == expected_points and len(comparisons) == 9 * expected_points
                  and not failures and all(row["passed"] for row in comparisons.values()))
    point_maxima = []
    for point in points:
        worst = max((row for (key, _), row in comparisons.items() if key == point),
                    key=lambda row: abs(row["relative_change"]))
        point_maxima.append(dict(point_id=point, observable=worst["observable"],
                                 relative_change=worst["relative_change"]))
    result = dict(analysis_only=True, numerical_defaults_promoted=False,
                  expected_points=expected_points, completed_points=len(points), comparisons=len(comparisons),
                  joint_target_passed=all_passed, configurations=configurations,
                  input_sha256=input_hash, maximum_changes=maxima, point_maxima=point_maxima,
                  failed_runs=failures, sources=sources)
    if benchmark is not None:
        status = json.loads((benchmark / "status.json").read_text(encoding="utf-8"))
        timing_manifest = json.loads((benchmark / "manifest.json").read_text(encoding="utf-8"))
        if not status.get("complete") or timing_manifest["parameters"]["candidate"] != configurations["candidate"]:
            raise ValueError("Complete timing is unfinished or uses a different candidate")
        if (numerical_sources(timing_manifest) != source_hashes
                or digest(benchmark / "inputs/review_points.csv") != input_hash):
            raise ValueError("Complete timing uses different numerical sources or inputs")
        rows = read_rows(benchmark / "timings.csv")
        if len(rows) != status["measurements"]:
            raise ValueError("Complete timing count differs from status")
        expected_timing_points = set(timing_manifest["parameters"]["points"])
        if ({row["point_id"] for row in rows} != expected_timing_points
                or status["points"] != len(expected_timing_points)):
            raise ValueError("Complete timing is missing an advertised point")
        if {row["method"] for row in rows} != {"none", "grouped", "thermal_weighted"}:
            raise ValueError("Unexpected or missing complete timing method")
        for row in rows:
            for field in ("total_seconds", "cpu_seconds"):
                value = float(row[field])
                if not math.isfinite(value) or value <= 0:
                    raise ValueError(f"Invalid {field} in complete timing")
        timings = []
        for point in sorted(expected_timing_points):
            for method in ("none", "grouped", "thermal_weighted"):
                selected = [row for row in rows if row["point_id"] == point and row["method"] == method]
                repeats = timing_manifest["parameters"]["repeats"]
                if len(selected) != repeats or {int(row["repeat"]) for row in selected} != set(range(1, repeats + 1)):
                    raise ValueError(f"Missing complete timing repetitions: {point}, {method}")
                timings.append(dict(point_id=point, method=method, repeats=len(selected),
                                    median_wall_seconds=statistics.median(float(row["total_seconds"]) for row in selected),
                                    median_cpu_seconds=statistics.median(float(row["cpu_seconds"]) for row in selected),
                                    min_wall_seconds=min(float(row["total_seconds"]) for row in selected),
                                    max_wall_seconds=max(float(row["total_seconds"]) for row in selected)))
        output_rows = read_rows(benchmark / "observables.csv")
        output_keys = [(row["point_id"], row["repeat"], row["method"]) for row in output_rows]
        timing_keys = {(row["point_id"], row["repeat"], row["method"]) for row in rows}
        if len(set(output_keys)) != len(output_keys) or set(output_keys) != timing_keys:
            raise ValueError("Complete timing outputs do not match its measurements")
        replay_changes, reference_changes = [], []
        for row in output_rows:
            if row["method"] != "thermal_weighted":
                continue
            for name in OBSERVABLES:
                key = row["point_id"], name
                if key not in comparisons:
                    raise ValueError(f"Timed candidate point lacks joint qualification: {key[0]}")
                value = float(row[name])
                if not math.isfinite(value) or value <= 0:
                    raise ValueError(f"Invalid timed candidate observable: {key}")
                replay_changes.append(abs(value / comparisons[key]["candidate"] - 1))
                reference_changes.append(abs(value / comparisons[key]["reference"] - 1))
        if max(reference_changes) > 0.01:
            raise ValueError("Timed candidate fails the joint one-percent target")
        result["complete_timing"] = dict(directory=str(benchmark.resolve()),
            manifest_sha256=digest(benchmark / "manifest.json"), timings_sha256=digest(benchmark / "timings.csv"),
            observables_sha256=digest(benchmark / "observables.csv"),
            candidate_replay_max_relative=max(replay_changes),
            candidate_reference_max_relative=max(reference_changes),
            points=timings, interpretation="Two-point host timings; not a GitHub runner measurement or a scan-wide bound")
    outdir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(__file__, outdir / Path(__file__).name)
    (outdir / "summary.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if comparisons:
        with (outdir / "comparisons.csv").open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=next(iter(comparisons.values())).keys())
            writer.writeheader()
            writer.writerows(comparisons[key] for key in sorted(comparisons))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", nargs="+", type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--expected-points", required=True, type=int)
    parser.add_argument("--benchmark", type=Path)
    args = parser.parse_args()
    result = summarize(args.runs, args.outdir, args.expected_points, args.benchmark)
    print(json.dumps({key: result[key] for key in ("completed_points", "comparisons", "joint_target_passed")}, indent=2))
    raise SystemExit(0 if result["joint_target_passed"] else 1)
