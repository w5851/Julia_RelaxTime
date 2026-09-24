#!/usr/bin/env python3
"""Build the label-only ``publication_clean_v5`` display derivative.

v5 inherits the accepted publication_clean_v4 display values and audit tables.
It does not reapply a numerical replacement, call a solver, or modify raw
results.  The only rendered changes are publication labels: relaxation-time
y axes carry the ``fm`` unit, and first-order endpoint legend entries use
branch terminology instead of the potentially misleading ``hadron`` label.
"""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import importlib.util
import json
import shutil
import subprocess
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
V4_SCRIPT = ROOT / "scripts" / "analysis" / "relaxtime" / "build_phase_guided_publication_clean_v4.py"
V4_SPEC = importlib.util.spec_from_file_location("phase_guided_publication_clean_v4_v5", V4_SCRIPT)
if V4_SPEC is None or V4_SPEC.loader is None:  # pragma: no cover
    raise RuntimeError(f"unable to load v4 builder: {V4_SCRIPT}")
V4 = importlib.util.module_from_spec(V4_SPEC)
V4_SPEC.loader.exec_module(V4)


TRANSPORT_ANALYSIS_ROOT = ROOT / "docs" / "analysis" / "relaxtime" / "phase_guided_transport"
V4_ROOT = TRANSPORT_ANALYSIS_ROOT / "phase_guided_transport_publication_clean_v4"
V5_ROOT = TRANSPORT_ANALYSIS_ROOT / "phase_guided_transport_publication_clean_v5"
V4_TABLE_ROOT = V4_ROOT / "tables"
V5_TABLE_ROOT = V5_ROOT / "tables"
V4_FIGURE_ROOT = V4_ROOT / "figures"
V5_FIGURE_ROOT = V5_ROOT / "figures"
V4_MANIFEST = V4_ROOT / "manifest.json"
V4_PLOT_MANIFEST = V4_FIGURE_ROOT / "plot_manifest.json"
V4_FIGURE_LAYER_MANIFEST = TRANSPORT_ANALYSIS_ROOT / "phase_guided_transport_publication_clean_figure_layer_v4" / "figure_layer_manifest.json"
CURRENT_POINTER = TRANSPORT_ANALYSIS_ROOT / "publication_clean_current.json"
V4_ELIGIBILITY_RECORD = TRANSPORT_ANALYSIS_ROOT / "publication_clean_v4_manuscript_eligibility_v1.json"
V5_MANIFEST = V5_ROOT / "manifest.json"
V5_PLOT_MANIFEST = V5_FIGURE_ROOT / "plot_manifest.json"

TAU_FIELDS = ("tau_u", "tau_d", "tau_s", "tau_ubar", "tau_dbar", "tau_sbar")
V4_OBSERVABLE_LABELS = dict(V4.V3.V2.OBSERVABLE_LABELS)
V4_PHASE_LEGEND_LABEL = V4.V3.V2.phase_legend_label
V5_OBSERVABLE_LABELS = dict(V4_OBSERVABLE_LABELS)
for _field in TAU_FIELDS:
    _base = V5_OBSERVABLE_LABELS[_field].rstrip("$")
    V5_OBSERVABLE_LABELS[_field] = f"{_base}\\;[\\mathrm{{fm}}]$"

V5_ENDPOINT_LEGEND_LABELS = {
    "quark": "chirally restored branch endpoint",
    "hadron": "chirally broken branch endpoint",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relpath(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def git_head() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def phase_legend_label(phase_curr: str) -> str:
    return V5_ENDPOINT_LEGEND_LABELS.get(
        phase_curr,
        f"phase-unresolved ({phase_curr or 'unknown'}) endpoint",
    )


def label_map_rows() -> list[dict[str, str]]:
    rows = []
    for field in TAU_FIELDS:
        rows.append(
            {
                "scope": "relaxation_time_y_axis",
                "key": field,
                "old_label": V4_OBSERVABLE_LABELS[field],
                "new_label": V5_OBSERVABLE_LABELS[field],
                "unit_or_semantics": "fm",
                "canonical_data_modified": "False",
            }
        )
    for phase_curr, new_label in V5_ENDPOINT_LEGEND_LABELS.items():
        rows.append(
            {
                "scope": "first_order_endpoint_legend",
                "key": phase_curr,
                "old_label": V4_PHASE_LEGEND_LABEL(phase_curr),
                "new_label": new_label,
                "unit_or_semantics": "branch endpoint wording only",
                "canonical_data_modified": "False",
            }
        )
    return rows


def load_v4_inputs() -> tuple[list[dict[str, str]], dict[str, Any], dict[str, Any], list[dict[str, Any]], dict[str, Any]]:
    if not all(path.is_file() for path in (V4_MANIFEST, V4_PLOT_MANIFEST, V4_FIGURE_LAYER_MANIFEST, CURRENT_POINTER, V4_ELIGIBILITY_RECORD)):
        raise FileNotFoundError("publication_clean_v4 must exist before building v5")
    package_manifest = read_json(V4_MANIFEST)
    plot_manifest = read_json(V4_PLOT_MANIFEST)
    figure_layer_manifest = read_json(V4_FIGURE_LAYER_MANIFEST)
    current_pointer = read_json(CURRENT_POINTER)
    eligibility_record = read_json(V4_ELIGIBILITY_RECORD)
    if package_manifest.get("schema") != "phase_guided_transport_publication_clean_manifest_v4":
        raise ValueError("unexpected publication_clean_v4 package manifest schema")
    if plot_manifest.get("schema") != "phase_guided_transport_publication_clean_plot_manifest_v4":
        raise ValueError("unexpected publication_clean_v4 plot manifest schema")
    if package_manifest.get("manuscript_eligible") is not True:
        raise ValueError("v4 package manifest is not the accepted display layer")
    if plot_manifest.get("manuscript_eligible") is not True:
        raise ValueError("v4 plot manifest is not the accepted display layer")
    if figure_layer_manifest.get("manuscript_eligible") is not True:
        raise ValueError("v4 figure layer is not the accepted display layer")
    if current_pointer.get("current_analysis_package") != relpath(V4_ROOT):
        raise ValueError("current publication pointer does not select publication_clean_v4")
    if current_pointer.get("manuscript_eligible") is not True or eligibility_record.get("manuscript_eligible") is not True:
        raise ValueError("v4 eligibility record does not accept the display layer")
    if len(plot_manifest.get("figures", [])) != 72:
        raise ValueError("publication_clean_v4 must contain exactly 72 figures")
    points_path = V4_TABLE_ROOT / "publication_clean_points.csv"
    points = V4.read_csv(points_path)
    expected = int(package_manifest["derived_counts"]["publication_clean_point_rows"])
    if len(points) != expected:
        raise ValueError(f"v4 point rows {len(points)} != manifest count {expected}")

    # Reuse only the v4 renderer's audited gap construction.  No solver or
    # numerical transformation is called; the plotted values come from v4.
    _, v3_parent_manifest = V4.load_parent()
    _, _, gaps = V4.load_raw_context(v3_parent_manifest)
    return points, package_manifest, plot_manifest, gaps, {
        "figure_layer_manifest": figure_layer_manifest,
        "current_pointer": current_pointer,
        "eligibility_record": eligibility_record,
    }


def build_figure_assets(figure_specs: list[dict[str, Any]]) -> list[dict[str, Any]]:
    assets = []
    for spec in figure_specs:
        path = spec["path"]
        assets.append(
            {
                "path": relpath(path),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
                "mode_key": spec["mode_key"],
                "plot_panel": spec["plot_panel"],
                "observable": spec["observable"],
                "axis_scale": spec["axis_scale"],
                "axis_scale_reason": spec["axis_scale_reason"],
                "data_min": spec["data_min"],
                "data_max": spec["data_max"],
                "dynamic_range_ratio": spec["dynamic_range_ratio"],
                "first_order_gap_present": spec["first_order_gap_present"],
            }
        )
    return assets


def package_output_paths() -> list[Path]:
    return [
        V5_ROOT / "README.md",
        *sorted(V5_TABLE_ROOT.glob("*.csv")),
        V5_PLOT_MANIFEST,
        *sorted(V5_FIGURE_ROOT.rglob("*.png")),
        *sorted((V5_ROOT / "audit").glob("*.png")),
    ]


def render_readme(v4_package: dict[str, Any], figure_count: int) -> str:
    return f"""# Issue #130 RS `publication_clean_v5` figure-label review candidate

## Purpose and boundary

This package is a label-only display derivative of the author-accepted
`publication_clean_v4` figure layer.  It reuses the v4 point table, display
adjustments, phase gaps, and audit evidence byte-for-byte.  It changes no raw
CSV, production registry, solver output, or numerical display value.

The two requested rendering changes are:

1. The six relaxation-time y-axis labels now include the unit `[fm]`.
2. The first-order endpoint legend entries now read `chirally restored branch
   endpoint` and `chirally broken branch endpoint`.  The latter is wording for
   the quark transport branch and does not claim that hadronic transport was
   calculated.

The v4 package remains unchanged and remains the current formal publication
layer until this candidate is explicitly reviewed and accepted.  v5 is
`manuscript_eligible=false`, `solver_called=false`, and is not a convergence
or physical-model update.

## Provenance

- parent package: `{relpath(V4_ROOT)}`
- parent package manifest SHA256: `{sha256_file(V4_MANIFEST)}`
- parent plot manifest SHA256: `{sha256_file(V4_PLOT_MANIFEST)}`
- source point table and audit tables: copied from v4 without content changes
- publication figures: {figure_count}
- raw/production data modified: false
- solver called for this derivative: false
- canonical numerical data modified: false

The exact old/new strings are recorded in
`tables/v5_display_label_map.csv`.  The v5 plot and package manifests record
the updated PNG hashes and the unchanged v4 parent hashes.

## Reproduction

    python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v5.py
    python scripts/analysis/relaxtime/sync_phase_guided_publication_clean_figure_layer_v5.py
    python -m pytest tests/unit/python/test_phase_guided_publication_clean_v5.py
"""


def main() -> None:
    if V5_ROOT.exists():
        raise FileExistsError(f"refusing to overwrite completed v5 package: {V5_ROOT}")
    points, v4_package, v4_plot, gaps, v4_acceptance = load_v4_inputs()
    v4_table_hashes = {
        path.relative_to(V4_TABLE_ROOT).as_posix(): sha256_file(path)
        for path in V4_TABLE_ROOT.glob("*.csv")
    }

    shutil.copytree(V4_ROOT, V5_ROOT)
    shutil.rmtree(V5_FIGURE_ROOT)
    V5_FIGURE_ROOT.mkdir(parents=True, exist_ok=True)

    # V4's renderer is the authoritative figure implementation.  Patch only
    # its text maps and output root; points and gaps are the v4 values above.
    renderer = V4.V3.V2
    renderer.FIGURE_DIR = V5_FIGURE_ROOT
    renderer.OBSERVABLE_LABELS = V5_OBSERVABLE_LABELS
    renderer.phase_legend_label = phase_legend_label
    _, figure_specs = renderer.render_figures(points, gaps, V4.OBSERVABLES)
    if len(figure_specs) != 72:
        raise ValueError(f"v5 renderer produced {len(figure_specs)} figures, expected 72")

    for relative, expected_hash in v4_table_hashes.items():
        actual = sha256_file(V5_TABLE_ROOT / relative)
        if actual != expected_hash:
            raise ValueError(f"v5 changed inherited table bytes: {relative}")

    label_rows = label_map_rows()
    write_csv(
        V5_TABLE_ROOT / "v5_display_label_map.csv",
        label_rows,
        ["scope", "key", "old_label", "new_label", "unit_or_semantics", "canonical_data_modified"],
    )

    generated_at = dt.datetime.now(dt.timezone.utc).isoformat()
    generator_path = Path(__file__).resolve()
    assets = build_figure_assets(figure_specs)
    label_update = {
        "kind": "display_label_only",
        "tau_y_axis_unit": "fm",
        "tau_fields": list(TAU_FIELDS),
        "endpoint_legend_labels": V5_ENDPOINT_LEGEND_LABELS,
        "canonical_data_modified": False,
        "solver_called": False,
    }

    plot_manifest = dict(v4_plot)
    plot_manifest.update(
        {
            "schema": "phase_guided_transport_publication_clean_plot_manifest_v5",
            "generated_at": generated_at,
            "base_git_commit": git_head(),
            "generator": relpath(generator_path),
            "generator_sha256": sha256_file(generator_path),
            "source_parent_artifact": relpath(V4_ROOT),
            "source_parent_manifest_sha256": sha256_file(V4_MANIFEST),
            "source_parent_plot_manifest_sha256": sha256_file(V4_PLOT_MANIFEST),
            "source_parent_figure_layer_manifest_sha256": sha256_file(V4_FIGURE_LAYER_MANIFEST),
            "source_parent_eligibility_record": relpath(V4_ELIGIBILITY_RECORD),
            "source_parent_eligibility_record_sha256": sha256_file(V4_ELIGIBILITY_RECORD),
            "manuscript_eligible": False,
            "canonical_data_modified": False,
            "production_write": False,
            "solver_called": False,
            "display_label_update": label_update,
            "rendering_semantics": str(v4_plot.get("rendering_semantics", ""))
            + "; v5 adds fm units to tau y axes and branch endpoint wording only",
            "figures": assets,
        }
    )
    write_json(V5_PLOT_MANIFEST, plot_manifest)

    readme_path = V5_ROOT / "README.md"
    readme_path.write_text(render_readme(v4_package, len(assets)), encoding="utf-8")

    package_manifest = dict(v4_package)
    package_manifest.update(
        {
            "schema": "phase_guided_transport_publication_clean_manifest_v5",
            "generated_at": generated_at,
            "base_git_commit": git_head(),
            "generator": relpath(generator_path),
            "generator_sha256": sha256_file(generator_path),
            "source_parent_artifact": relpath(V4_ROOT),
            "source_parent_manifest_sha256": sha256_file(V4_MANIFEST),
            "source_parent_plot_manifest_sha256": sha256_file(V4_PLOT_MANIFEST),
            "source_parent_figure_layer_manifest_sha256": sha256_file(V4_FIGURE_LAYER_MANIFEST),
            "source_parent_eligibility_record": relpath(V4_ELIGIBILITY_RECORD),
            "source_parent_eligibility_record_sha256": sha256_file(V4_ELIGIBILITY_RECORD),
            "status": "derived_author_review_required",
            "manuscript_eligible": False,
            "canonical_data_modified": False,
            "production_write": False,
            "solver_called": False,
            "display_label_update": label_update,
            "derived_counts": {
                **dict(v4_package.get("derived_counts", {})),
                "publication_figure_count": len(assets),
            },
            "known_boundaries": [
                *list(v4_package.get("known_boundaries", [])),
                "v5 is a label-only review derivative of v4; v4 remains unchanged",
                "tau axis units and endpoint wording do not alter numerical values or phase semantics",
            ],
            "outputs": [],
        }
    )
    package_manifest["outputs"] = [
        {"path": relpath(path), "sha256": sha256_file(path), "bytes": path.stat().st_size}
        for path in package_output_paths()
    ]
    write_json(V5_MANIFEST, package_manifest)

    print(
        json.dumps(
            {
                "output": relpath(V5_ROOT),
                "manifest": relpath(V5_MANIFEST),
                "publication_figures": len(assets),
                "label_rows": len(label_rows),
                "solver_called": False,
                "canonical_data_modified": False,
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
