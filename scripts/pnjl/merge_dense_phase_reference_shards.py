#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[2]

ARTIFACTS: dict[str, tuple[str, ...]] = {
    "boundary": ("xi", "T_MeV"),
    "cep": ("xi",),
    "spinodals": ("xi", "T_MeV"),
    "crossover": ("xi", "mu_MeV"),
    "grid_convergence": ("axis", "xi", "T_MeV", "level", "left", "right", "midpoint"),
}

CEP_LEGACY_COLUMNS = (
    "xi",
    "T_CEP_MeV",
    "muq_CEP_MeV",
    "muB_CEP_MeV",
    "uncertainty_T_MeV",
    "T_bracket_low_MeV",
    "T_bracket_high_MeV",
    "bracket_width_T_MeV",
)
CEP_COLUMNS = CEP_LEGACY_COLUMNS + (
    "result_status",
    "T_last_first_order_MeV",
    "muq_last_first_order_MeV",
    "muB_last_first_order_MeV",
    "T_first_monotone_MeV",
    "ambiguity_width_T_MeV",
    "temperature_resolution_target_MeV",
)


def fail(message: str) -> None:
    raise SystemExit(f"[dense-reference-merge] {message}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Deterministically merge PNJL dense-reference xi shards.")
    parser.add_argument("--shards-root", type=Path, required=True)
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument(
        "--xi-convergence-root",
        type=Path,
        help="optional root containing staged xi_grid_convergence_<tag>_level*.csv records",
    )
    parser.add_argument("--tag", required=True)
    parser.add_argument("--expected-xi-list", required=True)
    parser.add_argument(
        "--expected-calculation-git-commit",
        help="reject shards whose numerical calculation commit differs from this SHA",
    )
    parser.add_argument(
        "--postprocess-git-commit",
        help="commit providing merge/validation code; defaults to the calculation commit",
    )
    parser.add_argument(
        "--source-workflow-run-id",
        help="workflow run that produced the shard artifacts",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def parse_float_list(raw: str) -> list[float]:
    values = sorted({float(token.strip()) for token in raw.split(",") if token.strip()})
    if not values:
        fail("expected xi list is empty")
    if any(not math.isfinite(value) for value in values):
        fail("expected xi list contains a non-finite value")
    return values


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalized_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return resolved.as_posix()


def float_key(value: str) -> tuple[int, float | str]:
    text = value.strip()
    if text == "":
        return (1, "")
    try:
        number = float(text)
    except ValueError:
        return (2, text)
    return (0, number)


def row_key(row: dict[str, str], columns: tuple[str, ...]) -> tuple[tuple[int, float | str], ...]:
    return tuple(float_key(row.get(column, "")) for column in columns)


def load_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            fail(f"missing CSV header: {path}")
        fieldnames = list(reader.fieldnames)
        rows: list[dict[str, str]] = []
        for row in reader:
            extra = row.get(None) or []
            missing = [name for name in fieldnames if row.get(name) is None]
            if extra or missing:
                actual_count = len(fieldnames) + len(extra) - len(missing)
                fail(
                    f"malformed CSV row in {path} ending at physical line {reader.line_num}: "
                    f"expected {len(fieldnames)} fields, parsed {actual_count}; "
                    f"missing={missing}, extra_values={extra}"
                )
            rows.append(dict(row))
        return fieldnames, rows


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def normalize_config(config: dict[str, Any]) -> dict[str, Any]:
    ignored = {"xi_values", "requested_xi_values"}
    return {key: value for key, value in config.items() if key not in ignored}


def merge_csv_artifact(
    name: str,
    paths: list[Path],
    destination: Path,
) -> tuple[list[str], list[dict[str, str]]]:
    key_columns = ARTIFACTS[name]
    header: list[str] | None = None
    merged: dict[tuple[tuple[int, float | str], ...], dict[str, str]] = {}
    for path in sorted(paths, key=lambda item: item.as_posix()):
        current_header, rows = load_csv(path)
        if name == "cep":
            if tuple(current_header) == CEP_LEGACY_COLUMNS:
                converted: list[dict[str, str]] = []
                for row in rows:
                    has_point = bool(row.get("T_CEP_MeV", "").strip() and row.get("muq_CEP_MeV", "").strip())
                    converted.append({
                        **row,
                        "result_status": "resolved" if has_point else "not_found",
                        "T_last_first_order_MeV": row.get("T_bracket_low_MeV", ""),
                        "muq_last_first_order_MeV": row.get("muq_CEP_MeV", ""),
                        "muB_last_first_order_MeV": row.get("muB_CEP_MeV", ""),
                        "T_first_monotone_MeV": row.get("T_bracket_high_MeV", ""),
                        "ambiguity_width_T_MeV": row.get("bracket_width_T_MeV", ""),
                        "temperature_resolution_target_MeV": "",
                    })
                current_header = list(CEP_COLUMNS)
                rows = converted
            elif tuple(current_header) != CEP_COLUMNS:
                fail(f"unsupported CEP CSV schema in {path}: {current_header}")
        if header is None:
            header = current_header
        elif current_header != header:
            fail(f"CSV header mismatch for {name}: {path}")
        for row in rows:
            key = row_key(row, key_columns)
            previous = merged.get(key)
            if previous is not None and previous != row:
                fail(f"conflicting duplicate {name} row for key={key}: {path}")
            merged[key] = row
    if header is None:
        fail(f"no CSV shards found for {name}")
    ordered = [merged[key] for key in sorted(merged)]
    write_csv(destination, header, ordered)
    return header, ordered


def resolved_xi(rows_by_artifact: dict[str, list[dict[str, str]]]) -> list[float]:
    values: set[float] = set()
    for rows in rows_by_artifact.values():
        for row in rows:
            raw = row.get("xi", "").strip()
            if raw:
                values.add(float(raw))
    return sorted(values)


def update_dense_meaning(dense_meaning: dict[str, Any], xi_values: list[float]) -> dict[str, Any]:
    updated = json.loads(json.dumps(dense_meaning))
    sampling = updated.setdefault("xi_sampling", {})
    sampling["strategy"] = "adaptive_midpoint_grid" if len(xi_values) > 1 else "fixed_single_value"
    sampling["count"] = len(xi_values)
    sampling["values"] = xi_values
    diffs = [b - a for a, b in zip(xi_values, xi_values[1:])]
    sampling["min_step"] = min(diffs) if diffs else None
    sampling["max_step"] = max(diffs) if diffs else None
    sampling.pop("step", None)
    return updated


def main() -> None:
    args = parse_args()
    tag = args.tag
    manifests = sorted(args.shards_root.rglob(f"phase_reference_{tag}_manifest.json"))
    if not manifests:
        fail(f"no shard manifests found below {args.shards_root}")

    payloads = [json.loads(path.read_text(encoding="utf-8")) for path in manifests]
    commits = {payload.get("generator", {}).get("git_commit") for payload in payloads}
    if len(commits) != 1 or None in commits:
        fail(f"shards do not share one git commit: {sorted(str(value) for value in commits)}")
    calculation_git_commit = str(next(iter(commits)))
    if (
        args.expected_calculation_git_commit is not None
        and calculation_git_commit != args.expected_calculation_git_commit
    ):
        fail(
            "shard calculation commit mismatch: "
            f"expected {args.expected_calculation_git_commit}, got {calculation_git_commit}"
        )
    postprocess_git_commit = args.postprocess_git_commit or calculation_git_commit
    configs = [payload.get("config", {}) for payload in payloads]
    normalized = [normalize_config(config) for config in configs]
    if any(config != normalized[0] for config in normalized[1:]):
        fail("shard configuration mismatch outside xi sampling fields")

    args.reference_root.mkdir(parents=True, exist_ok=True)
    outputs = {
        "boundary": args.reference_root / f"boundary_{tag}.csv",
        "cep": args.reference_root / f"cep_{tag}.csv",
        "spinodals": args.reference_root / f"spinodals_{tag}.csv",
        "crossover": args.reference_root / f"crossover_{tag}.csv",
        "grid_convergence": args.reference_root / f"phase_grid_convergence_{tag}.csv",
    }
    meta_path = args.reference_root / f"crossover_{tag}.meta.json"
    manifest_path = args.reference_root / f"phase_reference_{tag}_manifest.json"
    for path in [*outputs.values(), meta_path, manifest_path]:
        if path.exists() and not args.overwrite:
            fail(f"output exists: {path}; rerun with --overwrite")

    rows_by_artifact: dict[str, list[dict[str, str]]] = {}
    for name, output in outputs.items():
        source_paths = sorted(path.parent / output.name for path in manifests if (path.parent / output.name).is_file())
        if not source_paths:
            if name in {"boundary", "cep", "spinodals"} and bool(configs[0].get("crossover_only")):
                continue
            fail(f"missing {name} CSV in all shards")
        _, rows = merge_csv_artifact(name, source_paths, output)
        rows_by_artifact[name] = rows

    xi_convergence_paths: list[Path] = []
    if args.xi_convergence_root is not None and args.xi_convergence_root.is_dir():
        xi_convergence_paths = sorted(
            args.xi_convergence_root.rglob(f"xi_grid_convergence_{tag}_level*.csv")
        )
    if xi_convergence_paths:
        grid_output = outputs["grid_convergence"]
        _, rows = merge_csv_artifact(
            "grid_convergence",
            [grid_output, *xi_convergence_paths],
            grid_output,
        )
        rows_by_artifact["grid_convergence"] = rows

    xis = resolved_xi(rows_by_artifact)
    expected = parse_float_list(args.expected_xi_list)
    missing = [value for value in expected if not any(abs(value - actual) <= 1e-12 for actual in xis)]
    if missing:
        fail(f"merged artifacts are missing requested xi anchors: {missing}")

    meta_sources = sorted(path.parent / meta_path.name for path in manifests if (path.parent / meta_path.name).is_file())
    if not meta_sources:
        fail("missing crossover metadata shards")
    meta_payloads = [json.loads(path.read_text(encoding="utf-8")) for path in meta_sources]
    column_definitions = meta_payloads[0].get("column_definitions", [])
    if any(payload.get("column_definitions", []) != column_definitions for payload in meta_payloads[1:]):
        fail("crossover metadata column definitions differ across shards")
    crossover_rows = rows_by_artifact.get("crossover", [])
    mu_values = sorted({float(row["mu_MeV"]) for row in crossover_rows if row.get("mu_MeV", "").strip()})
    generated_at = min(
        str(payload.get("generator", {}).get("generated_at", "")) for payload in payloads
    )
    crossover_meta = {
        "schema_version": meta_payloads[0].get("schema_version", "v1"),
        "artifact": {"path": normalized_path(outputs["crossover"]), "row_count": len(crossover_rows)},
        "generator": {
            "script": "scripts/pnjl/merge_dense_phase_reference_shards.py",
            "git_commit": calculation_git_commit,
            "generated_at": generated_at,
            "crossover_only": bool(configs[0].get("crossover_only")),
            "crossover_mu0_only": bool(configs[0].get("crossover_mu0_only")),
        },
        "provenance": {
            "calculation_git_commit": calculation_git_commit,
            "postprocess_git_commit": postprocess_git_commit,
            "source_workflow_run_id": args.source_workflow_run_id,
        },
        "xi_coverage": {
            "count": len(xis),
            "min": xis[0] if xis else None,
            "max": xis[-1] if xis else None,
            "values": xis,
        },
        "mu_q_coverage": {
            "count": len(mu_values),
            "min_MeV": mu_values[0] if mu_values else None,
            "max_MeV": mu_values[-1] if mu_values else None,
            "values_MeV": mu_values,
        },
        "column_definitions": column_definitions,
        "dense_meaning": update_dense_meaning(meta_payloads[0].get("dense_meaning", {}), xis),
    }
    meta_path.write_text(json.dumps(crossover_meta, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    runs: dict[float, dict[str, Any]] = {}
    for payload, source in zip(payloads, manifests):
        for run in payload.get("runs", []):
            xi = float(run["xi"])
            semantic = {key: value for key, value in run.items() if key not in {"run_dir"}}
            previous = runs.get(xi)
            if previous is not None:
                previous_semantic = {key: value for key, value in previous.items() if key not in {"run_dir", "source_manifest"}}
                if previous_semantic != semantic:
                    fail(f"conflicting duplicate run metadata for xi={xi}")
                if source.as_posix() >= str(previous["source_manifest"]):
                    continue
            runs[xi] = {**run, "source_manifest": source.as_posix()}

    final_config = dict(configs[0])
    final_config["requested_xi_values"] = expected
    final_config["xi_values"] = xis
    grid_rows = rows_by_artifact.get("grid_convergence", [])
    unconverged = sum(row.get("converged", "").strip().lower() != "true" for row in grid_rows)
    final_manifest = {
        "schema_version": payloads[0].get("schema_version", "v1"),
        "generator": {
            "script": "scripts/pnjl/merge_dense_phase_reference_shards.py",
            "git_commit": calculation_git_commit,
            "generated_at": generated_at,
            "crossover_only": bool(configs[0].get("crossover_only")),
            "crossover_mu0_only": bool(configs[0].get("crossover_mu0_only")),
        },
        "provenance": {
            "calculation_git_commit": calculation_git_commit,
            "postprocess_git_commit": postprocess_git_commit,
            "source_workflow_run_id": args.source_workflow_run_id,
        },
        "config": final_config,
        "output_root": payloads[0].get("output_root"),
        "reference_root": normalized_path(args.reference_root),
        "artifacts": {
            name: {"path": normalized_path(path), "row_count": len(rows_by_artifact.get(name, []))}
            for name, path in outputs.items()
            if name in rows_by_artifact
        },
        "runs": [{key: value for key, value in runs[xi].items() if key != "source_manifest"} for xi in sorted(runs)],
        "shards": [
            {"manifest": path.as_posix(), "sha256": sha256(path)} for path in manifests
        ],
        "xi_refinement_records": [
            {"path": path.as_posix(), "sha256": sha256(path)} for path in xi_convergence_paths
        ],
        "grid_convergence": {
            "record_count": len(grid_rows),
            "unconverged_count": unconverged,
        },
    }
    final_manifest["artifacts"]["crossover_meta"] = {"path": normalized_path(meta_path)}
    final_manifest["artifacts"]["manifest"] = {"path": normalized_path(manifest_path)}
    manifest_path.write_text(json.dumps(final_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"merged {len(manifests)} shards calculated at commit {calculation_git_commit}")
    print(f"postprocessed at commit {postprocess_git_commit}")
    print(f"resolved xi values: {xis}")
    print(f"manifest: {manifest_path}")


if __name__ == "__main__":
    main()
