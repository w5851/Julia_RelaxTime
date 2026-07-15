#!/usr/bin/env python3
"""Import and audit the Issue #130 phase-guided transport Action artifacts.

This importer is intentionally specific to the RS distribution/on-shell
transport-kernel semantic transition.  It merges sharded Action artifacts,
requires an already-passed p104-to-p128 convergence gate, verifies unchanged
thermodynamic and relaxation-time fields against the superseded dense-xi case,
records the expected transport drift, renders figures, and writes candidate
production directories without modifying the old case.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import struct
import subprocess
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path
from typing import Any

import merge_phase_guided_highxi_artifacts as legacy


PROJECT_ROOT = Path(__file__).resolve().parents[3]
RESULT_ROOT = legacy.RESULT_ROOT
FIGURE_ROOT = legacy.FIGURE_ROOT
MODE_A = legacy.MODE_A
MODE_B = legacy.MODE_B

TAU_FIELDS = [
    "tau_u",
    "tau_d",
    "tau_s",
    "tau_ubar",
    "tau_dbar",
    "tau_sbar",
    "tauinv_u",
    "tauinv_d",
    "tauinv_s",
    "tauinv_ubar",
    "tauinv_dbar",
    "tauinv_sbar",
]
TRANSPORT_FIELDS = [
    "eta",
    "sigma",
    "zeta",
    "eta_over_s",
    "zeta_over_s",
    "sigma_over_T",
]
UNCHANGED_FIELDS = [
    field for field in legacy.SCAN_CORE_COLUMNS if field not in TRANSPORT_FIELDS
]
NONNEGATIVE_SCAN_FIELDS = TAU_FIELDS + ["eta", "sigma", "zeta"]
PLOT_Y_COLUMNS = [
    "tau_u",
    "tau_d",
    "tau_s",
    "tau_ubar",
    "tau_dbar",
    "tau_sbar",
    "eta",
    "sigma",
    "zeta",
    "eta_over_s",
    "zeta_over_s",
    "sigma_over_T",
]
TEXT_OUTPUT_SUFFIXES = {".csv", ".json", ".md"}


def write_text_lf(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)


def write_json_lf(path: Path, payload: Any) -> None:
    write_text_lf(path, json.dumps(payload, indent=2, ensure_ascii=False) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download-root", required=True, type=Path)
    parser.add_argument("--run-manifest", required=True, type=Path)
    parser.add_argument("--convergence-dir", required=True, type=Path)
    parser.add_argument("--new-case", required=True)
    parser.add_argument("--source-candidate-case", required=True)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument(
        "--old-case",
        default="first_canonical_v1_p128_xi001_validated_anchored_prod_v1",
    )
    parser.add_argument("--invariance-rel-threshold", type=float, default=1e-3)
    parser.add_argument("--abs-floor", type=float, default=1e-10)
    parser.add_argument("--overlap-rel-threshold", type=float, default=1e-3)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--skip-plots", action="store_true")
    parser.add_argument("--validate-only", action="store_true")
    return parser.parse_args()


def compare_old_new(
    mode: str,
    old_case: str,
    new_rows: list[dict[str, str]],
    abs_floor: float,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    old_path = RESULT_ROOT / mode / old_case / "phase_guided_transport_scan.csv"
    old_rows = legacy.read_scan_csv(old_path).rows
    old_by_key = {legacy.old_anchor_key(mode, row): row for row in old_rows}
    new_keys = {legacy.old_anchor_key(mode, row) for row in new_rows}
    comparisons: list[dict[str, Any]] = []
    missing: list[str] = []

    for new_row in new_rows:
        key = legacy.old_anchor_key(mode, new_row)
        old_row = old_by_key.get(key)
        if old_row is None:
            missing.append("|".join(key))
            continue
        xi = float(new_row["xi"])
        for field in legacy.scan_validity_columns(mode):
            old = legacy.parse_float(old_row.get(field, ""))
            new = legacy.parse_float(new_row.get(field, ""))
            if old is None or new is None:
                continue
            if not math.isfinite(old) and not math.isfinite(new):
                continue
            if not (math.isfinite(old) and math.isfinite(new)):
                abs_delta = math.inf
                rel_delta = math.inf
                signed = math.inf
            else:
                abs_delta = abs(new - old)
                rel_delta = abs_delta / max(abs(new), abs(old), abs_floor)
                signed = (new - old) / max(abs(old), abs_floor)
            if field in TRANSPORT_FIELDS:
                role = "xi0_invariance" if abs(xi) <= 1e-12 else "transport_semantic_drift"
            else:
                role = "unchanged_contract"
            comparisons.append(
                {
                    "mode": mode,
                    "key": "|".join(key),
                    "xi": xi,
                    "field": field,
                    "role": role,
                    "old": old,
                    "new": new,
                    "abs_delta": abs_delta,
                    "rel_delta": rel_delta,
                    "signed_rel_delta": signed,
                }
            )

    return comparisons, {
        "old_rows": len(old_rows),
        "new_rows": len(new_rows),
        "matched_rows": len({row["key"] for row in comparisons}),
        "new_rows_missing_from_old": missing,
        "old_rows_missing_from_new": [
            "|".join(key) for key in sorted(set(old_by_key) - new_keys)
        ],
    }


def worst_row(rows: list[dict[str, Any]]) -> dict[str, Any] | None:
    return max(rows, key=lambda row: float(row["rel_delta"]), default=None)


def summarize_fields(rows: list[dict[str, Any]]) -> dict[str, Any]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[str(row["field"])].append(row)
    result: dict[str, Any] = {}
    for field, values in grouped.items():
        worst = worst_row(values)
        signed = [float(row["signed_rel_delta"]) for row in values]
        result[field] = {
            "count": len(values),
            "max_rel": float(worst["rel_delta"]) if worst else 0.0,
            "min_signed_rel": min(signed),
            "max_signed_rel": max(signed),
            "worst": worst,
        }
    return result


def write_comparison(path: Path, rows: list[dict[str, Any]]) -> None:
    header = [
        "mode",
        "key",
        "xi",
        "field",
        "role",
        "old",
        "new",
        "abs_delta",
        "rel_delta",
        "signed_rel_delta",
    ]
    legacy.write_csv(
        path,
        [],
        header,
        [{name: str(row.get(name, "")) for name in header} for row in rows],
    )


def validate_convergence_gate(path: Path, source_commit: str) -> dict[str, Any]:
    summary_path = path / "convergence_summary.json"
    if not summary_path.is_file():
        raise FileNotFoundError(f"convergence summary not found: {summary_path}")
    summary = json.loads(summary_path.read_text(encoding="utf-8-sig"))
    if summary.get("verdict") != "production-grade":
        raise ValueError(f"convergence verdict is not production-grade: {summary.get('verdict')}")
    if summary.get("source_commit") != source_commit:
        raise ValueError("convergence source commit does not match production source commit")
    return summary


def make_effective_config(
    mode: str,
    case_name: str,
    source_candidate_case: str,
    source_commit: str,
    rows: list[dict[str, str]],
) -> dict[str, Any]:
    config = legacy.effective_config(mode, case_name, rows)
    config["options"]["transport_kernel_energy_policy"] = "onshell_kinematic_energy"
    config["options"]["distribution_energy_policy"] = "rs_deformed"
    config["production_import"] = {
        "verdict": "production-grade",
        "manuscript_status": "current_candidate",
        "source_candidate_case": source_candidate_case,
        "formal_case_name": case_name,
        "source_commit": source_commit,
        "path_fields_rewritten_for_repository_import": True,
        "numerical_options": "p128_validated_anchored",
        "convergence_gate": "convergence/parameter_gate/convergence_summary.json",
        "old_new_gate": "convergence/production_import_summary.json",
        "transport_kernel_energy_policy": "onshell_kinematic_energy",
        "distribution_energy_policy": "rs_deformed",
        "propagator_xi_policy": "match_thermo",
    }
    return config


def write_readme_and_audit(
    case_dir: Path,
    figure_dir: Path,
    mode: str,
    case_name: str,
    source_candidate_case: str,
    source_commit: str,
    rows: list[dict[str, str]],
    source_runs: list[dict[str, Any]],
    summary: dict[str, Any],
) -> None:
    if mode == MODE_A:
        scope = "fixed mu_B with phase-scaled temperature bands"
    else:
        scope = "fixed temperature with sparse mu_B values"
    xis = legacy.collect_xi_values(rows)
    run_table = "\n".join(
        f"| {run['label']} | [{run['run_id']}]({run['url']}) |"
        for run in source_runs
    )
    figure_rel = figure_dir.relative_to(PROJECT_ROOT).as_posix()
    comparison = "convergence/old_vs_new_full_grid_comparison.csv"
    readme = f"""# {case_name} ({legacy.mode_title(mode)})

This is the candidate replacement for the superseded RS-shared-energy production case after Issue #130.

## Status

- numerical verdict: `production-grade`
- production registry status: `current_candidate`
- manuscript eligible: `false` until PR review approves promotion

## Physics policies

- transport kernel energy: `onshell_kinematic_energy`, E_kin=sqrt(p^2+M^2)
- distribution energy: `rs_deformed`, E_dist=sqrt(p^2+M^2+xi*(p*cos(theta))^2)
- propagator xi policy: `match_thermo`
- sigma cache policy: `validated_anchored`

## Scope

- mode: `{mode}` ({scope})
- xi grid: `{xis[0]:.2f}:0.01:{xis[-1]:.2f}`
- rows: `{len(rows)}`
- source commit: `{source_commit}`
- source candidate case: `{source_candidate_case}`

## Evidence

- p104-to-p128 gate: `convergence/parameter_gate/convergence_summary.json`
- full-grid old/new comparison: `{comparison}`
- production import summary: `convergence/production_import_summary.json`
- result manifest: `manifest.json`
- figure directory: `{figure_rel}`
- plot manifest: `{figure_rel}/plot_manifest.json`

The old production directories are retained byte-for-byte and are not modified by this case.
"""
    audit = f"""# Production Audit: {case_name} ({legacy.mode_title(mode)})

## Verdict

`production-grade` numerical artifact; registry status remains `current_candidate` pending human review.

## Production Case

- case: `{case_name}`
- mode: `{mode}`
- source commit: `{source_commit}`
- source candidate: `{source_candidate_case}`
- workflow: `.github/workflows/relaxtime-phase-guided-transport-production.yml`

## Physics Scope

- `{scope}`
- xi grid `{xis[0]:.2f}:0.01:{xis[-1]:.2f}`
- eta, sigma, zeta and their normalized ratios, plus six relaxation times and inverse times
- on-shell E_kin in explicit transport kernels; RS-deformed E_dist only in distributions

## Non-goals

- no change to s, v_rel, t limits, Pauli blocking, or relaxation-time semantics
- no longitudinal/transverse transport tensor decomposition
- this workflow does not output kappa_XY or lambda columns

## Selected Production Parameters

- tau nodes: `128/20/36/24`
- sigma grid: `560`
- sigma cache: `validated_anchored`
- propagator xi policy: `match_thermo`
- bulk viscosity: enabled
- channel diagnostics: enabled

## Convergence And Import Gates

- p104-to-p128 gate verdict: `production-grade`
- unchanged-contract worst relative drift: `{summary['unchanged_max_rel']:.8g}`
- tau old/new worst relative drift: `{summary['tau_max_rel']:.8g}`
- xi=0 transport worst relative drift: `{summary['xi0_transport_max_rel']:.8g}`
- source failed rows: `{summary['source_failed_rows']}`
- scan non-finite values: `{summary['scan_nonfinite']}`
- diagnostic non-finite values: `{summary['diagnostic_nonfinite']}`
- diagnostic negative rate/contribution values: `{summary['diagnostic_negative']}`

## Source Action Runs

| shard | run |
| --- | --- |
{run_table}

## Data Outputs

- `phase_guided_transport_scan.csv`
- `sampling_plan.csv`
- `channel_diagnostics.csv`
- `failed_points.csv`
- `effective_config.json`
- `convergence/`

## Figure Outputs

- `{figure_rel}/plot_manifest.json`
- `{figure_rel}/**/*.png`

## Known Limitations And Residual Risks

- Action outputs were `diagnostic-only`; this local audit is the explicit production gate.
- Dense-xi local structures still require neighboring-point and channel-diagnostic interpretation.
- Manuscript eligibility remains false until the independent production PR is reviewed.
"""
    write_text_lf(case_dir / "README.md", readme)
    write_text_lf(case_dir / "PRODUCTION_AUDIT.md", audit)


def make_manifest(
    case_dir: Path,
    figure_dir: Path,
    mode: str,
    case_name: str,
    source_candidate_case: str,
    source_commit: str,
    source_runs: list[dict[str, Any]],
    summary: dict[str, Any],
) -> dict[str, Any]:
    result_files = [
        "phase_guided_transport_scan.csv",
        "sampling_plan.csv",
        "channel_diagnostics.csv",
        "failed_points.csv",
        "effective_config.json",
        "README.md",
        "PRODUCTION_AUDIT.md",
        "convergence/production_download_manifest.csv",
        "convergence/old_vs_new_full_grid_comparison.csv",
        "convergence/production_import_summary.json",
        "convergence/parameter_gate/manifest.json",
    ]
    hashes = {
        name: legacy.sha256_file(case_dir / name)
        for name in result_files
        if (case_dir / name).is_file()
    }
    return {
        "format": "phase_guided_transport_energy_semantics_production_manifest_v1",
        "case_slug": case_name,
        "mode": mode,
        "date": date.today().isoformat(),
        "verdict": "production-grade",
        "registry_status": "current_candidate",
        "manuscript_eligible": False,
        "source_commit": source_commit,
        "remote_workflow": ".github/workflows/relaxtime-phase-guided-transport-production.yml",
        "source_candidate_case": source_candidate_case,
        "source_runs": {row["label"]: str(row["run_id"]) for row in source_runs},
        "production_parameters": {
            "tau_p_nodes": 128,
            "tau_angle_nodes": 20,
            "tau_phi_nodes": 36,
            "tau_n_sigma_points": 24,
            "sigma_grid_n": 560,
        },
        "policies": {
            "transport_kernel_energy_policy": "onshell_kinematic_energy",
            "distribution_energy_policy": "rs_deformed",
            "propagator_xi_policy": "match_thermo",
            "sigma_cache_policy": "validated_anchored",
            "channel_diagnostics": True,
            "compute_bulk": True,
        },
        "gates": summary,
        "result_files": result_files + ["manifest.json"],
        "hashes": hashes,
        "figure_files": [
            figure_dir.relative_to(PROJECT_ROOT).as_posix() + "/plot_manifest.json",
            figure_dir.relative_to(PROJECT_ROOT).as_posix() + "/**/*.png",
        ],
    }


def enrich_plot_manifest(case_dir: Path, figure_dir: Path) -> None:
    manifest_path = figure_dir / "plot_manifest.json"
    payload = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    source_csv = case_dir / "phase_guided_transport_scan.csv"
    figure_hashes = {
        rel.as_posix(): legacy.sha256_file(figure_dir / rel)
        for rel in sorted(
            (path.relative_to(figure_dir) for path in figure_dir.rglob("*.png")),
            key=lambda value: value.as_posix(),
        )
    }
    payload.update(
        {
            "schema_version": "v2",
            "source_csv": source_csv.relative_to(PROJECT_ROOT).as_posix(),
            "source_csv_sha256": legacy.sha256_file(source_csv),
            "plot_entrypoint": "scripts/relaxtime/run_phase_guided_transport_plots.jl",
            "x": "xi",
            "y_columns": PLOT_Y_COLUMNS,
            "format": "png",
            "dpi": 600,
            "figure_hashes": figure_hashes,
        }
    )
    write_json_lf(manifest_path, payload)


def refresh_manifest(case_dir: Path, figure_dir: Path) -> None:
    manifest_path = case_dir / "manifest.json"
    payload = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    payload["hashes"] = {
        name: legacy.sha256_file(case_dir / name)
        for name in payload["result_files"]
        if name != "manifest.json" and (case_dir / name).is_file()
    }
    payload["figure_summary"] = {
        "png_count": len(list(figure_dir.rglob("*.png"))),
        "plot_manifest_sha256": legacy.sha256_file(figure_dir / "plot_manifest.json"),
    }
    write_json_lf(manifest_path, payload)


def png_dpi(path: Path) -> tuple[float, float] | None:
    """Read PNG pHYs metadata without adding a Pillow dependency."""
    with path.open("rb") as stream:
        if stream.read(8) != b"\x89PNG\r\n\x1a\n":
            return None
        while True:
            raw_length = stream.read(4)
            if len(raw_length) != 4:
                return None
            length = struct.unpack(">I", raw_length)[0]
            chunk_type = stream.read(4)
            data = stream.read(length)
            stream.read(4)
            if chunk_type == b"pHYs" and len(data) == 9:
                pixels_x, pixels_y, unit = struct.unpack(">IIB", data)
                if unit == 1:
                    return pixels_x * 0.0254, pixels_y * 0.0254
                return None
            if chunk_type == b"IEND":
                return None


def finalize_audit(
    case_dir: Path,
    figure_dir: Path,
    mode: str,
    summary: dict[str, Any],
) -> None:
    audit_path = case_dir / "PRODUCTION_AUDIT.md"
    plot_manifest_path = figure_dir / "plot_manifest.json"
    payload = json.loads(plot_manifest_path.read_text(encoding="utf-8-sig"))
    figure_hashes = payload.get("figure_hashes", {})
    figure_hashes_ok = all(
        legacy.sha256_file(figure_dir / relative_path) == expected
        for relative_path, expected in figure_hashes.items()
    )
    unexpected_figure_files = [
        path.relative_to(figure_dir).as_posix()
        for path in figure_dir.rglob("*")
        if path.is_file()
        and path.suffix.lower() != ".png"
        and path.name != "plot_manifest.json"
    ]
    dpi_values = [
        dpi
        for relative_path in figure_hashes
        if (dpi := png_dpi(figure_dir / relative_path)) is not None
    ]
    dpi_min = min((min(value) for value in dpi_values), default=math.nan)
    dpi_max = max((max(value) for value in dpi_values), default=math.nan)
    if not figure_hashes or len(figure_hashes) != int(payload.get("count", -1)):
        raise ValueError(f"{mode}: figure manifest count/hash coverage mismatch")
    if not figure_hashes_ok:
        raise ValueError(f"{mode}: one or more figure SHA256 values do not match")
    if unexpected_figure_files:
        raise ValueError(f"{mode}: unexpected figure-side files: {unexpected_figure_files}")
    if payload.get("dpi") != 600 or not dpi_values or dpi_min < 599.0:
        raise ValueError(
            f"{mode}: plot DPI gate failed (target={payload.get('dpi')}, min={dpi_min})"
        )
    case_name = case_dir.name
    case_manifest = json.loads(
        (case_dir / "manifest.json").read_text(encoding="utf-8-sig")
    )
    source_commit = case_manifest["source_commit"]
    source_candidate_case = case_manifest["source_candidate_case"]
    marker = "## Final Local Validation And Command Log"
    text = audit_path.read_text(encoding="utf-8")
    if marker in text:
        text = text.split(marker, 1)[0].rstrip() + "\n\n"
    else:
        text = text.rstrip() + "\n\n"
    section = f"""{marker}

- mode: `{mode}`.
- Action collector: `collect_phase_guided_action_artifacts.py --case-name {source_candidate_case} --source-commit {source_commit} --expected-count 30 --download`.
- Import gate: `import_phase_guided_energy_semantics_artifacts.py ... --new-case {case_name} --validate-only` passed before repository writes.
- Formal import: the same importer without `--validate-only` wrote the new slug only.
- Figure command: `julia --project=. scripts/relaxtime/run_phase_guided_transport_plots.jl --case-dir {case_dir.relative_to(PROJECT_ROOT).as_posix()} --overwrite`.
- final scan rows: `{summary['scan_rows']}`; diagnostic rows: `{summary['diagnostic_rows']}`; converged rows: `{summary['converged_counts'].get('true', 0)}`.
- failed rows: `{summary['source_failed_rows']}`; scan NaN/Inf: `{summary['scan_nonfinite']}`; scan negative values: `{summary['scan_negative']}`.
- diagnostic NaN/Inf: `{summary['diagnostic_nonfinite']}`; diagnostic negative rate/contribution values: `{summary['diagnostic_negative']}`.
- figure manifest: `{len(figure_hashes)}` PNG SHA256 entries; hash verification: `passed`.
- plot DPI target: `{payload.get('dpi')}`; PNG pHYs metadata range: `{dpi_min:.4f}` to `{dpi_max:.4f}` dpi. The sub-600 display is PNG pixels-per-meter quantization of the 600 dpi target.
- figure-side hygiene: `passed`; only PNG files and `plot_manifest.json` remain after provenance sidecar cleanup.
- result and figure manifests are refreshed after this section is written.
- old-case byte-integrity is checked independently with pre/post tree SHA256 and recorded in the M2 task document and production PR.
"""
    write_text_lf(audit_path, text + section)


def normalize_text_tree_lf(root: Path) -> None:
    for path in root.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in TEXT_OUTPUT_SUFFIXES:
            continue
        raw = path.read_bytes()
        normalized = raw.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
        if normalized != raw:
            path.write_bytes(normalized)


def refresh_parameter_gate_manifest(case_dir: Path) -> None:
    manifest_path = case_dir / "convergence" / "parameter_gate" / "manifest.json"
    payload = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    payload["hashes"] = {
        name: legacy.sha256_file(manifest_path.parent / name)
        for name in payload.get("files", [])
    }
    write_json_lf(manifest_path, payload)


def finalize_repository_outputs(case_dir: Path, figure_dir: Path) -> None:
    normalize_text_tree_lf(case_dir)
    normalize_text_tree_lf(figure_dir)
    refresh_parameter_gate_manifest(case_dir)

    plot_manifest_path = figure_dir / "plot_manifest.json"
    plot_payload = json.loads(plot_manifest_path.read_text(encoding="utf-8-sig"))
    plot_payload["source_csv_sha256"] = legacy.sha256_file(
        case_dir / "phase_guided_transport_scan.csv"
    )
    write_json_lf(plot_manifest_path, plot_payload)

    refresh_manifest(case_dir, figure_dir)
    result_manifest = json.loads(
        (case_dir / "manifest.json").read_text(encoding="utf-8-sig")
    )
    mismatches = [
        name
        for name, expected in result_manifest.get("hashes", {}).items()
        if legacy.sha256_file(case_dir / name) != expected
    ]
    if mismatches:
        raise ValueError(f"result manifest hash mismatch after LF finalization: {mismatches}")


def main() -> int:
    args = parse_args()
    convergence_summary = validate_convergence_gate(
        args.convergence_dir.resolve(), args.source_commit
    )
    run_rows = json.loads(args.run_manifest.read_text(encoding="utf-8-sig"))
    discovered = legacy.discover_artifacts(args.download_root.resolve(), run_rows)
    scan_rows_raw, diag_rows_raw, source_manifest_rows = legacy.load_all_rows(discovered)

    legacy.SOURCE_CANDIDATE_CASE = args.source_candidate_case
    mode_dirs = {mode: RESULT_ROOT / mode / args.new_case for mode in (MODE_A, MODE_B)}
    figure_dirs = {mode: FIGURE_ROOT / mode / args.new_case for mode in (MODE_A, MODE_B)}

    merged_scans: dict[str, list[dict[str, str]]] = {}
    merged_diags: dict[str, list[dict[str, str]]] = {}
    summaries: dict[str, dict[str, Any]] = {}
    comparisons: dict[str, list[dict[str, Any]]] = {}
    gate_failures: list[str] = []

    if any(row["head_sha"] != args.source_commit for row in source_manifest_rows):
        gate_failures.append("one or more production runs use a different head SHA")
    if any(row["conclusion"] != "success" for row in source_manifest_rows):
        gate_failures.append("one or more production runs did not conclude successfully")

    for mode in (MODE_A, MODE_B):
        merged_scans[mode], overlap = legacy.merge_rows(
            scan_rows_raw[mode],
            legacy.scan_key,
            legacy.scan_sort_key,
            legacy.SCAN_CORE_COLUMNS,
            args.overlap_rel_threshold,
        )
        merged_diags[mode], diagnostic_overlap = legacy.merge_rows(
            diag_rows_raw[mode],
            legacy.diag_key,
            legacy.diag_sort_key,
            legacy.DIAG_COMPARE_COLUMNS,
            args.overlap_rel_threshold,
        )
        if not overlap["passed"]:
            gate_failures.append(f"{mode}: scan overlap gate failed")
        if not diagnostic_overlap["passed"]:
            gate_failures.append(f"{mode}: diagnostic overlap gate failed")

        mode_comparisons, coverage = compare_old_new(
            mode, args.old_case, merged_scans[mode], args.abs_floor
        )
        comparisons[mode] = mode_comparisons
        unchanged_rows = [row for row in mode_comparisons if row["role"] == "unchanged_contract"]
        tau_rows = [row for row in mode_comparisons if row["field"] in TAU_FIELDS]
        xi0_rows = [row for row in mode_comparisons if row["role"] == "xi0_invariance"]
        effective_unchanged_rows = [
            row for row in unchanged_rows if float(row["abs_delta"]) > args.abs_floor
        ]
        unchanged_worst = worst_row(effective_unchanged_rows)
        tau_worst = worst_row(tau_rows)
        xi0_worst = worst_row(xi0_rows)
        unchanged_max = float(unchanged_worst["rel_delta"]) if unchanged_worst else 0.0
        tau_max = float(tau_worst["rel_delta"]) if tau_worst else math.inf
        xi0_max = float(xi0_worst["rel_delta"]) if xi0_worst else math.inf
        stable_violations = [
            row
            for row in unchanged_rows
            if float(row["abs_delta"]) > args.abs_floor
            and float(row["rel_delta"]) > args.invariance_rel_threshold
        ]
        xi0_violations = [
            row
            for row in xi0_rows
            if float(row["abs_delta"]) > args.abs_floor
            and float(row["rel_delta"]) > args.invariance_rel_threshold
        ]
        if coverage["new_rows_missing_from_old"] or coverage["old_rows_missing_from_new"]:
            gate_failures.append(f"{mode}: old/new comparison is missing rows")
        if stable_violations:
            gate_failures.append(
                f"{mode}: unchanged-contract violations={len(stable_violations)}"
            )
        if xi0_violations:
            gate_failures.append(f"{mode}: xi=0 transport violations={len(xi0_violations)}")

        scan_nonfinite = 0
        scan_negative = 0
        for row in merged_scans[mode]:
            for field in legacy.scan_validity_columns(mode):
                value = legacy.parse_float(row.get(field, ""))
                if value is not None and not math.isfinite(value):
                    scan_nonfinite += 1
            for field in NONNEGATIVE_SCAN_FIELDS:
                value = legacy.parse_float(row.get(field, ""))
                if value is None or not math.isfinite(value):
                    continue
                if value < -1e-12:
                    scan_negative += 1
        diag_validity = legacy.numeric_validity(
            merged_diags[mode], legacy.DIAG_COMPARE_COLUMNS
        )
        source_runs = [row for row in source_manifest_rows if row["mode"] == mode]
        source_failed = sum(int(row["failed_rows"]) for row in source_runs)
        if source_failed:
            gate_failures.append(f"{mode}: source failed rows={source_failed}")
        if scan_nonfinite or scan_negative:
            gate_failures.append(
                f"{mode}: scan validity nonfinite={scan_nonfinite}, negative={scan_negative}"
            )
        if diag_validity["nonfinite"] or diag_validity["negative"]:
            gate_failures.append(f"{mode}: diagnostic validity failed")

        summaries[mode] = {
            "convergence_verdict": convergence_summary["verdict"],
            "source_runs": len(source_runs),
            "source_failed_rows": source_failed,
            "scan_rows": len(merged_scans[mode]),
            "diagnostic_rows": len(merged_diags[mode]),
            "scan_overlap": overlap,
            "diagnostic_overlap": diagnostic_overlap,
            "coverage": coverage,
            "unchanged_max_rel": unchanged_max,
            "unchanged_worst": unchanged_worst,
            "tau_max_rel": tau_max,
            "tau_worst": tau_worst,
            "xi0_transport_max_rel": xi0_max,
            "xi0_transport_worst": xi0_worst,
            "transport_drift_field_summary": summarize_fields(
                [row for row in mode_comparisons if row["role"] == "transport_semantic_drift"]
            ),
            "scan_nonfinite": scan_nonfinite,
            "scan_negative": scan_negative,
            "diagnostic_nonfinite": diag_validity["nonfinite"],
            "diagnostic_negative": diag_validity["negative"],
            "converged_counts": dict(Counter(row.get("converged", "") for row in merged_scans[mode])),
            "quality_reason_counts": dict(Counter(row.get("quality_reason", "") for row in merged_scans[mode])),
        }

    if gate_failures:
        raise SystemExit("Production import blocked: " + "; ".join(gate_failures))

    if args.validate_only:
        print(
            json.dumps(
                {
                    "ready": True,
                    "validate_only": True,
                    "new_case": args.new_case,
                    "source_commit": args.source_commit,
                    "summaries": summaries,
                },
                indent=2,
                ensure_ascii=False,
            )
        )
        return 0

    legacy.cleanup_outputs(
        [*mode_dirs.values(), *figure_dirs.values()], args.overwrite
    )

    scan_header = list(legacy.read_scan_csv(discovered[0]["scan"]).header)
    diag_header = list(legacy.read_scan_csv(discovered[0]["diagnostics"]).header)
    failed_header = list(legacy.read_scan_csv(discovered[0]["failed"]).header)

    for mode in (MODE_A, MODE_B):
        case_dir = mode_dirs[mode]
        figure_dir = figure_dirs[mode]
        convergence_dir = case_dir / "convergence"
        parameter_gate_dir = convergence_dir / "parameter_gate"
        case_dir.mkdir(parents=True, exist_ok=True)
        shutil.copytree(args.convergence_dir.resolve(), parameter_gate_dir)
        source_runs = [row for row in source_manifest_rows if row["mode"] == mode]
        legacy.write_download_manifest(
            convergence_dir / "production_download_manifest.csv", source_runs
        )
        write_comparison(
            convergence_dir / "old_vs_new_full_grid_comparison.csv",
            comparisons[mode],
        )
        legacy.write_json(
            convergence_dir / "production_import_summary.json", summaries[mode]
        )

        source_manifest_rel = "convergence/production_download_manifest.csv"
        scan_comments = legacy.production_comments(
            args.new_case,
            args.source_commit,
            source_manifest_rel,
            "phase_guided_transport_scan",
        )
        diag_comments = legacy.production_comments(
            args.new_case,
            args.source_commit,
            source_manifest_rel,
            "phase_guided_transport_scan_channel_diagnostics",
            mode_path=(case_dir / "phase_guided_transport_scan.csv").relative_to(PROJECT_ROOT).as_posix(),
        )
        failed_comments = [
            f"# production_case: {args.new_case}",
            f"# source_candidate_case: {args.source_candidate_case}",
            f"# source_sharded_action_runs: see {source_manifest_rel}",
        ]
        legacy.write_csv(
            case_dir / "phase_guided_transport_scan.csv",
            scan_comments,
            scan_header,
            merged_scans[mode],
        )
        legacy.write_csv(
            case_dir / "channel_diagnostics.csv",
            diag_comments,
            diag_header,
            merged_diags[mode],
        )
        legacy.write_csv(
            case_dir / "failed_points.csv", failed_comments, failed_header, []
        )
        legacy.write_csv(
            case_dir / "sampling_plan.csv",
            [],
            ["T_MeV", "muB_MeV", "xi", "mode", "plot_panel", "plot_series"],
            legacy.make_plan_rows(merged_scans[mode]),
        )
        legacy.write_json(
            case_dir / "effective_config.json",
            make_effective_config(
                mode,
                args.new_case,
                args.source_candidate_case,
                args.source_commit,
                merged_scans[mode],
            ),
        )
        write_readme_and_audit(
            case_dir,
            figure_dir,
            mode,
            args.new_case,
            args.source_candidate_case,
            args.source_commit,
            merged_scans[mode],
            source_runs,
            summaries[mode],
        )
        legacy.write_json(
            case_dir / "manifest.json",
            make_manifest(
                case_dir,
                figure_dir,
                mode,
                args.new_case,
                args.source_candidate_case,
                args.source_commit,
                source_runs,
                summaries[mode],
            ),
        )

        if not args.skip_plots:
            subprocess.run(
                [
                    "julia",
                    "--project=.",
                    "scripts/relaxtime/run_phase_guided_transport_plots.jl",
                    "--case-dir",
                    str(case_dir),
                    "--overwrite",
                ],
                cwd=PROJECT_ROOT,
                check=True,
            )
            enrich_plot_manifest(case_dir, figure_dir)
            finalize_audit(case_dir, figure_dir, mode, summaries[mode])
            finalize_repository_outputs(case_dir, figure_dir)

    print(
        json.dumps(
            {
                "new_case": args.new_case,
                "source_commit": args.source_commit,
                "verdict": "production-grade",
                "registry_status": "current_candidate",
                "result_dirs": {mode: str(path) for mode, path in mode_dirs.items()},
                "figure_dirs": {mode: str(path) for mode, path in figure_dirs.items()},
                "summaries": summaries,
            },
            indent=2,
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
