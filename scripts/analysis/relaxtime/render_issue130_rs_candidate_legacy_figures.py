#!/usr/bin/env python3
"""Render Issue #130 candidate/legacy RS scans with the production plot layout.

This is an audit-only renderer.  It reuses the existing Julia plotting entrypoint
and writes outside ``data/outputs/figures`` by default.  The numerical CSVs are
read as immutable inputs; no solver or transport calculation is executed here.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import tempfile
from pathlib import Path
from typing import Any


SCHEMA = "issue130_rs_candidate_legacy_figures_v1"
CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
CASE_NAME = "first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2"
OLD_CASE_NAME = "first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v1"
MODES = ("mode_a_fixed_muB_phase_scaled", "mode_b_fixed_T_sparse_muB")
REFERENCES = ("candidate_runtime", "legacy_runtime")
Y_COLUMNS = (
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
)
COMMON_COLUMNS = (
    "T_MeV",
    "muB_MeV",
    "xi",
    "mode",
    "plot_panel",
    "plot_series_label",
    *Y_COLUMNS,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    lines = [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not lines:
        raise ValueError(f"CSV has no data/header rows: {path}")
    reader = csv.DictReader(lines)
    if reader.fieldnames is None:
        raise ValueError(f"CSV has no header: {path}")
    return list(reader.fieldnames), list(reader)


def finite(value: str, *, allow_nan: bool = False) -> bool:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    if math.isfinite(number):
        return True
    return allow_nan and math.isnan(number)


def _key_fields(mode: str) -> tuple[str, ...]:
    if mode == "mode_a_fixed_muB_phase_scaled":
        return ("muB_MeV", "xi", "alpha_T")
    if mode == "mode_b_fixed_T_sparse_muB":
        return ("T_MeV", "muB_MeV", "xi")
    raise ValueError(f"unsupported phase-guided mode: {mode}")


def validate_scan(path: Path, expected_mode: str) -> dict[str, Any]:
    fields, rows = read_csv(path)
    missing = [field for field in COMMON_COLUMNS if field not in fields]
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(missing)}")
    if not rows:
        raise ValueError(f"{path} has no numerical rows")

    mode_values = {row.get("mode", "") for row in rows}
    if mode_values != {expected_mode}:
        raise ValueError(f"{path} mode values {sorted(mode_values)} != [{expected_mode!r}]")

    key_fields = _key_fields(expected_mode)
    seen: set[tuple[str, ...]] = set()
    duplicate_keys: list[tuple[str, ...]] = []
    nonfinite: list[tuple[int, str]] = []
    for index, row in enumerate(rows):
        key = tuple(row.get(field, "") for field in key_fields)
        if key in seen:
            duplicate_keys.append(key)
        seen.add(key)
        for field in ("T_MeV", "muB_MeV", "xi", *Y_COLUMNS):
            if not finite(row.get(field, "")):
                nonfinite.append((index, field))
        if expected_mode == "mode_a_fixed_muB_phase_scaled" and not finite(row.get("alpha_T", "")):
            nonfinite.append((index, "alpha_T"))
    if duplicate_keys:
        raise ValueError(f"{path} has {len(duplicate_keys)} duplicate numerical keys")
    if nonfinite:
        preview = ", ".join(f"row={row} field={field}" for row, field in nonfinite[:5])
        raise ValueError(f"{path} has non-finite plotted/input values: {preview}")

    return {
        "path": str(path),
        "sha256": sha256(path),
        "mode": expected_mode,
        "rows": len(rows),
        "key_fields": list(key_fields),
        "unique_keys": len(seen),
        "plot_panels": sorted({row["plot_panel"] for row in rows}),
        "series": sorted({row["plot_series_label"] for row in rows}),
        "solver_called": False,
    }


def _figure_hashes(fig_dir: Path) -> dict[str, str]:
    return {
        path.relative_to(fig_dir).as_posix(): sha256(path)
        for path in sorted(fig_dir.rglob("*.png"))
    }


def _write_plot_manifest(fig_dir: Path, source: dict[str, Any], command: list[str]) -> dict[str, Any]:
    hashes = _figure_hashes(fig_dir)
    payload: dict[str, Any] = {
        "schema_version": "v2",
        "figures": sorted(hashes),
        "count": len(hashes),
        "source_csv": source["path"],
        "source_csv_sha256": source["sha256"],
        "plot_entrypoint": "scripts/relaxtime/run_phase_guided_transport_plots.jl",
        "x": "xi",
        "y_columns": list(Y_COLUMNS),
        "format": "png",
        "dpi": 600,
        "figure_hashes": hashes,
        "solver_called": False,
        "command": command,
    }
    path = fig_dir / "plot_manifest.json"
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def _render_one(
    *,
    repo_root: Path,
    julia_cmd: str,
    csv_path: Path,
    mode: str,
    reference: str,
    output_root: Path,
    source: dict[str, Any],
) -> dict[str, Any]:
    case_name = CASE_NAME if reference == "candidate_runtime" else OLD_CASE_NAME
    mode_root = output_root / mode / reference / case_name
    mode_root.mkdir(parents=True, exist_ok=True)
    command = [
        julia_cmd,
        f"--project={repo_root}",
        str(repo_root / "scripts" / "relaxtime" / "run_phase_guided_transport_plots.jl"),
        "--case-dir",
        "<temporary-case-dir>",
        "--csv",
        str(csv_path),
        "--fig-dir",
        str(mode_root),
        "--overwrite",
    ]
    with tempfile.TemporaryDirectory(prefix="issue130_rs_plot_") as temp_dir:
        case_dir = Path(temp_dir) / "case"
        case_dir.mkdir()
        # The historical wrapper prefers effective_config.json for mode
        # detection.  This tiny temporary file is routing metadata only; it is
        # never copied to the audit artifact and contains no numerical state.
        (case_dir / "effective_config.json").write_text(
            json.dumps({"options": {"mode": mode}}, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        run_command = [
            julia_cmd,
            f"--project={repo_root}",
            str(repo_root / "scripts" / "relaxtime" / "run_phase_guided_transport_plots.jl"),
            "--case-dir",
            str(case_dir),
            "--csv",
            str(csv_path),
            "--fig-dir",
            str(mode_root),
            "--overwrite",
        ]
        subprocess.run(run_command, cwd=repo_root, check=True)
    manifest = _write_plot_manifest(mode_root, source, command)
    if manifest["count"] != 36:
        raise RuntimeError(f"expected 36 figures for {mode}/{reference}, found {manifest['count']}")
    return {
        "mode": mode,
        "reference": reference,
        "figure_root": str(mode_root),
        "plot_manifest": str(mode_root / "plot_manifest.json"),
        "figure_count": manifest["count"],
        "figure_hashes": manifest["figure_hashes"],
        "source": source,
        "solver_called": False,
    }


def render(args: argparse.Namespace) -> dict[str, Any]:
    repo_root = Path(args.repo_root).resolve()
    output_root = Path(args.output_root).resolve()
    if not output_root.is_absolute():
        raise ValueError("--output-root must be absolute")
    official_figures = (repo_root / "data" / "outputs" / "figures").resolve()
    if output_root == official_figures or official_figures in output_root.parents:
        raise ValueError("audit renderer refuses to write under data/outputs/figures")
    output_root.mkdir(parents=True, exist_ok=True)

    inputs = {
        ("mode_a_fixed_muB_phase_scaled", "candidate_runtime"): Path(args.candidate_mode_a).resolve(),
        ("mode_a_fixed_muB_phase_scaled", "legacy_runtime"): Path(args.legacy_mode_a).resolve(),
        ("mode_b_fixed_T_sparse_muB", "candidate_runtime"): Path(args.candidate_mode_b).resolve(),
        ("mode_b_fixed_T_sparse_muB", "legacy_runtime"): Path(args.legacy_mode_b).resolve(),
    }
    source_records: dict[str, Any] = {}
    for (mode, reference), path in inputs.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        source_records[f"{mode}/{reference}"] = validate_scan(path, mode)

    rendered = []
    for (mode, reference), path in inputs.items():
        rendered.append(
            _render_one(
                repo_root=repo_root,
                julia_cmd=args.julia,
                csv_path=path,
                mode=mode,
                reference=reference,
                output_root=output_root,
                source=source_records[f"{mode}/{reference}"],
            )
        )

    manifest = {
        "schema_version": SCHEMA,
        "generated_utc": __import__("datetime").datetime.now(__import__("datetime").timezone.utc).isoformat(),
        "calculation_sha": CALCULATION_SHA,
        "workflow_head_sha": args.workflow_head_sha,
        "source_aggregate": str(Path(args.aggregate_root).resolve()),
        "repo_root": str(repo_root),
        "plot_contract": {
            "entrypoint": "scripts/relaxtime/run_phase_guided_transport_plots.jl",
            "x": "xi",
            "y_columns": list(Y_COLUMNS),
            "split": "plot_panel",
            "group": "plot_series_label",
            "format": "png",
            "dpi": 600,
            "same_layout_as_legacy": True,
        },
        "solver_called": False,
        "production_write": False,
        "inputs": source_records,
        "rendered": rendered,
    }
    (output_root / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_root / "README.md").write_text(
        """# Issue #130 RS candidate/legacy audit figures

This directory contains solver-free diagnostic renders using the same
plot entrypoint, panel split, line grouping, variables, PNG format, and
600 dpi setting as the historical phase-guided production figures.

The four trees are `mode_a_fixed_muB_phase_scaled` and
`mode_b_fixed_T_sparse_muB`, each with `candidate_runtime` and
`legacy_runtime`. The images are for author review only; they do not
replace or modify formal production figures. See `manifest.json` and
the per-tree `plot_manifest.json` files for source and image hashes.
""",
        encoding="utf-8",
    )
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--aggregate-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--candidate-mode-a", required=True, type=Path)
    parser.add_argument("--legacy-mode-a", required=True, type=Path)
    parser.add_argument("--candidate-mode-b", required=True, type=Path)
    parser.add_argument("--legacy-mode-b", required=True, type=Path)
    parser.add_argument("--workflow-head-sha", required=True)
    parser.add_argument("--julia", default="julia")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    render(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
