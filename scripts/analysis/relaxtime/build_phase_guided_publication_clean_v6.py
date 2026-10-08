#!/usr/bin/env python3
"""Build the publication_clean_v6 figure-style review layer.

v6 is a solver-free display derivative of the accepted publication_clean_v5
layer.  It keeps v5 point tables, audit tables, and provenance byte-for-byte
and regenerates only the PNG figures with publication-oriented tick and legend
styling.  No solver, convergence gate, or numerical production path is
called.
"""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import importlib.util
import json
import math
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[3]
V5_SCRIPT = ROOT / "scripts" / "analysis" / "relaxtime" / "build_phase_guided_publication_clean_v5.py"
V5_SPEC = importlib.util.spec_from_file_location("phase_guided_publication_clean_v5_v6", V5_SCRIPT)
if V5_SPEC is None or V5_SPEC.loader is None:  # pragma: no cover
    raise RuntimeError(f"unable to load v5 builder: {V5_SCRIPT}")
V5 = importlib.util.module_from_spec(V5_SPEC)
V5_SPEC.loader.exec_module(V5)

V4 = V5.V4
V3 = V4.V3
V2 = V3.V2

TRANSPORT_ANALYSIS_ROOT = ROOT / "docs" / "analysis" / "relaxtime" / "phase_guided_transport"
V5_ROOT = TRANSPORT_ANALYSIS_ROOT / "phase_guided_transport_publication_clean_v5"
V6_ROOT = TRANSPORT_ANALYSIS_ROOT / "phase_guided_transport_publication_clean_v6"
V5_TABLE_ROOT = V5_ROOT / "tables"
V6_TABLE_ROOT = V6_ROOT / "tables"
V5_FIGURE_ROOT = V5_ROOT / "figures"
V6_FIGURE_ROOT = V6_ROOT / "figures"
V5_MANIFEST = V5_ROOT / "manifest.json"
V5_PLOT_MANIFEST = V5_FIGURE_ROOT / "plot_manifest.json"
V5_FIGURE_LAYER_MANIFEST = (
    TRANSPORT_ANALYSIS_ROOT / "phase_guided_transport_publication_clean_figure_layer_v5" / "figure_layer_manifest.json"
)
CURRENT_POINTER = TRANSPORT_ANALYSIS_ROOT / "publication_clean_current.json"
V6_MANIFEST = V6_ROOT / "manifest.json"
V6_PLOT_MANIFEST = V6_FIGURE_ROOT / "plot_manifest.json"
V6_LABEL_MAP = V6_TABLE_ROOT / "v6_display_label_map.csv"
V6_STYLE_MANIFEST = V6_ROOT / "v6_display_style_manifest.json"

OBSERVABLES = list(V4.OBSERVABLES)
OBSERVABLE_LABELS = dict(V5.V5_OBSERVABLE_LABELS)

V6_ENDPOINT_LEGEND_LABELS = {
    "quark": "1st-order transition (chirally restored endpoint)",
    "hadron": "1st-order transition (chirally broken endpoint)",
}
V6_ENDPOINT_RENDER_LABELS = {
    phase: label.replace(" (", "\n(", 1)
    for phase, label in V6_ENDPOINT_LEGEND_LABELS.items()
}
LEGEND_FONTSIZE_PT = 20.0
LINEAR_MINOR_SUBDIVISIONS = 5
LOG_MINOR_SUBS = tuple(range(2, 10))


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


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        lines = [line for line in handle if line.strip() and not line.startswith("#")]
    return list(csv.DictReader(lines))


def write_csv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def phase_legend_label(phase_curr: str) -> str:
    return V6_ENDPOINT_RENDER_LABELS.get(
        phase_curr,
        f"phase-unresolved ({phase_curr or 'unknown'}) endpoint",
    )


def label_map_rows() -> list[dict[str, str]]:
    return [
        {
            "scope": "first_order_endpoint_legend",
            "key": phase,
            "old_label": V5.V5_ENDPOINT_LEGEND_LABELS[phase],
            "new_label": V6_ENDPOINT_LEGEND_LABELS[phase],
            "rendered_label": V6_ENDPOINT_RENDER_LABELS[phase],
            "font_size_pt": f"{LEGEND_FONTSIZE_PT:g}",
            "canonical_data_modified": "False",
        }
        for phase in ("quark", "hadron")
    ]


def load_v5_inputs() -> tuple[list[dict[str, str]], dict[str, Any], dict[str, Any], dict[str, Any], list[dict[str, Any]], dict[str, str]]:
    required = (V5_MANIFEST, V5_PLOT_MANIFEST, V5_FIGURE_LAYER_MANIFEST, CURRENT_POINTER)
    if not all(path.is_file() for path in required):
        raise FileNotFoundError("publication_clean_v5 must exist before building v6")

    package_manifest = read_json(V5_MANIFEST)
    plot_manifest = read_json(V5_PLOT_MANIFEST)
    figure_layer_manifest = read_json(V5_FIGURE_LAYER_MANIFEST)
    current_pointer = read_json(CURRENT_POINTER)
    if package_manifest.get("schema") != "phase_guided_transport_publication_clean_manifest_v5":
        raise ValueError("unexpected publication_clean_v5 package manifest schema")
    if plot_manifest.get("schema") != "phase_guided_transport_publication_clean_plot_manifest_v5":
        raise ValueError("unexpected publication_clean_v5 plot manifest schema")
    if figure_layer_manifest.get("schema") != "phase_guided_transport_publication_clean_figure_layer_manifest_v5":
        raise ValueError("unexpected publication_clean_v5 figure-layer manifest schema")
    if any(
        manifest.get("manuscript_eligible") is not True
        for manifest in (package_manifest, plot_manifest, figure_layer_manifest)
    ):
        raise ValueError("accepted v5 parent must remain manuscript_eligible=true")
    if current_pointer.get("current_analysis_package") != relpath(V5_ROOT):
        raise ValueError("v5 is not the current publication-clean analysis package")
    if current_pointer.get("manuscript_eligible") is not True:
        raise ValueError("current publication pointer no longer accepts v5")
    if len(plot_manifest.get("figures", [])) != 72:
        raise ValueError("publication_clean_v5 must contain exactly 72 figures")

    points_path = V5_TABLE_ROOT / "publication_clean_points.csv"
    points = read_csv(points_path)
    expected = int(package_manifest["derived_counts"]["publication_clean_point_rows"])
    if len(points) != expected:
        raise ValueError(f"v5 point rows {len(points)} != manifest count {expected}")
    for row in points:
        value = float(row["clean_value"])
        if not math.isfinite(value):
            raise ValueError(f"non-finite v5 clean value: {row}")

    # Reuse the already audited v2 phase gaps through the v4 loader.  This is
    # read-only provenance traversal; it does not invoke a solver.
    _, v3_parent_manifest = V4.load_parent()
    _, _, gaps = V4.load_raw_context(v3_parent_manifest)
    table_hashes = {
        path.name: sha256_file(path)
        for path in sorted(V5_TABLE_ROOT.glob("*.csv"))
    }
    return points, package_manifest, plot_manifest, figure_layer_manifest, gaps, table_hashes


def configure_axis_ticks(ax: Any, axis_scale: str) -> None:
    """Apply the v6 four-sided, inward major/minor tick contract."""
    from matplotlib.ticker import AutoMinorLocator, LogLocator, NullFormatter

    ax.xaxis.set_minor_locator(AutoMinorLocator(LINEAR_MINOR_SUBDIVISIONS))
    if axis_scale == "log":
        ax.yaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0,)))
        ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=LOG_MINOR_SUBS))
        ax.yaxis.set_minor_formatter(NullFormatter())
    elif axis_scale == "linear":
        ax.yaxis.set_minor_locator(AutoMinorLocator(LINEAR_MINOR_SUBDIVISIONS))
    else:
        raise ValueError(f"unsupported v6 axis scale: {axis_scale}")

    ax.tick_params(
        axis="both",
        which="major",
        direction="in",
        top=True,
        bottom=True,
        left=True,
        right=True,
        length=5.0,
        width=0.9,
    )
    ax.tick_params(
        axis="both",
        which="minor",
        direction="in",
        top=True,
        bottom=True,
        left=True,
        right=True,
        length=2.8,
        width=0.7,
    )


def render_v6_figures(
    points: list[dict[str, str]],
    boundary_gaps: list[dict[str, Any]],
    observables: Iterable[str],
) -> tuple[list[Path], list[dict[str, Any]]]:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("matplotlib is required to render publication-clean figures") from exc

    colors = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE"]
    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in points:
        grouped[(row["mode_key"], row["plot_panel"], row["plot_series"])].append(row)
    gaps_by_curve: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    for gap in boundary_gaps:
        gaps_by_curve[(gap["mode_key"], gap["plot_panel"], gap["plot_series"])].append(gap)

    paths: list[Path] = []
    figure_specs: list[dict[str, Any]] = []
    for mode_key in V2.MODE_CONFIG:
        panels = sorted({key[1] for key in grouped if key[0] == mode_key})
        for panel in panels:
            series_names = sorted({key[2] for key in grouped if key[:2] == (mode_key, panel)})
            for observable in observables:
                fig, ax = plt.subplots(figsize=(6.75, 4.6))
                phase_labels_seen: set[str] = set()
                for series_index, series in enumerate(series_names):
                    rows = [
                        row
                        for row in points
                        if row["mode_key"] == mode_key
                        and row["plot_panel"] == panel
                        and row["plot_series"] == series
                        and row["observable"] == observable
                    ]
                    rows.sort(key=lambda row: float(row["xi"]))
                    gaps = gaps_by_curve.get((mode_key, panel, series), [])
                    curve_color = colors[series_index % len(colors)]
                    segments = V2.split_curve_segments(rows, gaps)
                    for segment_index, segment in enumerate(segments):
                        ax.plot(
                            [float(row["xi"]) for row in segment],
                            [float(row["clean_value"]) for row in segment],
                            color=curve_color,
                            linewidth=1.5,
                            label=(
                                V2.display_series_label(mode_key, segment[0])
                                if segment_index == 0
                                else None
                            ),
                        )
                    point_index = {V2.canonical_xi(row["xi"]): row for row in rows}
                    for gap in gaps:
                        for endpoint_role, phase_marker, phase_edge in (
                            ("left", "o", "#3366CC"),
                            ("right", "s", "#CC3311"),
                        ):
                            endpoint = gap["endpoint_rows"][endpoint_role]
                            xi = V2.canonical_xi(endpoint["xi"])
                            point = point_index.get(xi)
                            if point is None:
                                raise ValueError(f"render endpoint missing from v5 point table: {gap['boundary_id']} {xi}")
                            label = V2.phase_label(endpoint["phase_curr"])
                            legend_label = phase_legend_label(endpoint["phase_curr"])
                            if label in phase_labels_seen:
                                legend_label = None
                            else:
                                phase_labels_seen.add(label)
                            ax.scatter(
                                [float(point["xi"])],
                                [float(point["clean_value"])],
                                marker=phase_marker,
                                s=48,
                                facecolor="white",
                                edgecolor=phase_edge,
                                linewidth=0.9,
                                zorder=6,
                                label=legend_label,
                            )

                figure_rows = [
                    row
                    for row in points
                    if row["mode_key"] == mode_key
                    and row["plot_panel"] == panel
                    and row["observable"] == observable
                ]
                panel_gaps = [
                    gap
                    for series in series_names
                    for gap in gaps_by_curve.get((mode_key, panel, series), [])
                ]
                axis_spec = V2.figure_axis_spec(figure_rows, panel_gaps)
                if axis_spec["axis_scale"] == "log":
                    ax.set_yscale("log")
                ax.set_xlabel(r"$\xi$")
                ax.set_ylabel(OBSERVABLE_LABELS[observable])
                ax.set_xlim(-0.52, 0.52)
                configure_axis_ticks(ax, axis_spec["axis_scale"])
                ax.legend(loc="best", fontsize=LEGEND_FONTSIZE_PT)
                fig.tight_layout()
                path = V6_FIGURE_ROOT / mode_key / f"plot_panel={panel}" / f"{observable}_vs_xi.png"
                path.parent.mkdir(parents=True, exist_ok=True)
                fig.savefig(path, dpi=600, bbox_inches="tight", pad_inches=0.08)
                plt.close(fig)
                paths.append(path)
                figure_specs.append(
                    {
                        "path": path,
                        "mode_key": mode_key,
                        "plot_panel": panel,
                        "observable": observable,
                        **axis_spec,
                    }
                )
    return paths, figure_specs


def build_figure_assets(figure_specs: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "path": relpath(spec["path"]),
            "sha256": sha256_file(spec["path"]),
            "bytes": spec["path"].stat().st_size,
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
        for spec in figure_specs
    ]


def package_output_paths() -> list[Path]:
    return [
        V6_ROOT / "README.md",
        V6_STYLE_MANIFEST,
        *sorted(V6_TABLE_ROOT.glob("*.csv")),
        V6_PLOT_MANIFEST,
        *sorted(V6_FIGURE_ROOT.rglob("*.png")),
        *sorted((V6_ROOT / "audit").glob("*.png")),
    ]


def style_manifest_payload(v5_package: dict[str, Any], generator_path: Path) -> dict[str, Any]:
    return {
        "schema": "phase_guided_transport_publication_clean_v6_display_style_manifest_v1",
        "status": "derived_author_review_required",
        "manuscript_eligible": False,
        "canonical_data_modified": False,
        "solver_called": False,
        "source_parent_artifact": relpath(V5_ROOT),
        "source_parent_manifest_sha256": sha256_file(V5_MANIFEST),
        "source_parent_plot_manifest_sha256": sha256_file(V5_PLOT_MANIFEST),
        "source_parent_figure_layer_manifest": relpath(V5_FIGURE_LAYER_MANIFEST),
        "source_parent_figure_layer_manifest_sha256": sha256_file(V5_FIGURE_LAYER_MANIFEST),
        "source_calculation_sha": v5_package.get("calculation_sha"),
        "source_workflow_head_sha": v5_package.get("workflow_head_sha"),
        "generator": relpath(generator_path),
        "generator_sha256": sha256_file(generator_path),
        "legend": {
            "font_size_pt": LEGEND_FONTSIZE_PT,
            "old_endpoint_labels": V5.V5_ENDPOINT_LEGEND_LABELS,
            "new_endpoint_labels": V6_ENDPOINT_LEGEND_LABELS,
            "rendered_endpoint_labels": V6_ENDPOINT_RENDER_LABELS,
            "line_break_policy": "insert newline before endpoint parentheses for source single-panel figures",
        },
        "ticks": {
            "major_ticks": True,
            "minor_ticks": True,
            "sides": ["top", "bottom", "left", "right"],
            "direction": "in",
            "linear_minor_locator": "matplotlib.ticker.AutoMinorLocator(5)",
            "log_major_locator": "matplotlib.ticker.LogLocator(base=10, subs=(1,))",
            "log_minor_locator": "matplotlib.ticker.LogLocator(base=10, subs=(2,3,4,5,6,7,8,9))",
            "x_axis": "linear",
            "y_axis": "inherited panel scale; linear or log",
        },
        "figure_count": 72,
        "mode_counts": {"mode_a": 36, "mode_b": 36},
    }


def render_readme(parent_manifest: dict[str, Any], figure_count: int) -> str:
    return f"""# Issue #130 RS `publication_clean_v6` figure-style review candidate

## Purpose and boundary

This package is a solver-free display derivative of the accepted
`publication_clean_v5` layer.  It preserves the v5 point table, audit tables,
smooth-display records, and raw provenance byte-for-byte.  v6 changes only
PNG rendering style and endpoint legend wording.

The v6 figure contract is:

- major and minor ticks on every axis;
- ticks on all four sides and directed inward;
- `AutoMinorLocator(5)` on linear axes and explicit `LogLocator` minor ticks
  on log axes;
- 20 pt legends in source single-panel figures;
- first-order endpoint labels rendered as `1st-order transition` followed by
  a line-broken `(chirally restored endpoint)` or `(chirally broken endpoint)`.

The v5 package and `publication_clean_current.json` remain unchanged and
continue to define the current formal publication layer.  v6 is a review
candidate with `manuscript_eligible=false`.

## Provenance

- parent package: `{relpath(V5_ROOT)}`
- parent package manifest SHA256: `{sha256_file(V5_MANIFEST)}`
- parent plot manifest SHA256: `{sha256_file(V5_PLOT_MANIFEST)}`
- parent figure-layer manifest SHA256: `{sha256_file(V5_FIGURE_LAYER_MANIFEST)}`
- publication figures: {figure_count} (mode A: 36; mode B: 36)
- solver called for this derivative: false
- canonical numerical data modified: false

The exact label mapping is recorded in `tables/v6_display_label_map.csv` and
the complete style contract plus parent hashes is recorded in
`v6_display_style_manifest.json`.

## Reproduction

```powershell
python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v6.py
python scripts/analysis/relaxtime/sync_phase_guided_publication_clean_figure_layer_v6.py
python -m pytest tests/unit/python/test_phase_guided_publication_clean_v6.py
```
"""


def main() -> None:
    if V6_ROOT.exists():
        raise FileExistsError(f"refusing to overwrite completed v6 package: {V6_ROOT}")

    points, v5_package, v5_plot, v5_layer, gaps, v5_table_hashes = load_v5_inputs()
    shutil.copytree(V5_ROOT, V6_ROOT)
    shutil.rmtree(V6_FIGURE_ROOT)
    V6_FIGURE_ROOT.mkdir(parents=True, exist_ok=True)

    for name, expected_hash in v5_table_hashes.items():
        actual_hash = sha256_file(V6_TABLE_ROOT / name)
        if actual_hash != expected_hash:
            raise ValueError(f"v6 changed inherited v5 table bytes: {name}")

    figure_paths, figure_specs = render_v6_figures(points, gaps, OBSERVABLES)
    if len(figure_specs) != 72:
        raise ValueError(f"v6 renderer produced {len(figure_specs)} figures, expected 72")
    mode_counts = {"mode_a": 0, "mode_b": 0}
    for spec in figure_specs:
        mode_counts[spec["mode_key"]] += 1
    if mode_counts != {"mode_a": 36, "mode_b": 36}:
        raise ValueError(f"unexpected v6 mode counts: {mode_counts}")

    generated_at = dt.datetime.now(dt.timezone.utc).isoformat()
    generator_path = Path(__file__).resolve()
    style_manifest = style_manifest_payload(v5_package, generator_path)
    style_manifest["generated_at"] = generated_at
    write_json(V6_STYLE_MANIFEST, style_manifest)
    write_csv(
        V6_LABEL_MAP,
        label_map_rows(),
        ["scope", "key", "old_label", "new_label", "rendered_label", "font_size_pt", "canonical_data_modified"],
    )

    assets = build_figure_assets(figure_specs)
    style_update = {
        "kind": "display_style_and_label_only",
        "legend_fontsize_pt": LEGEND_FONTSIZE_PT,
        "endpoint_legend_labels": V6_ENDPOINT_LEGEND_LABELS,
        "rendered_endpoint_legend_labels": V6_ENDPOINT_RENDER_LABELS,
        "tick_contract": style_manifest["ticks"],
        "canonical_data_modified": False,
        "solver_called": False,
    }

    plot_manifest = dict(v5_plot)
    plot_manifest.update(
        {
            "schema": "phase_guided_transport_publication_clean_plot_manifest_v6",
            "generated_at": generated_at,
            "base_git_commit": git_head(),
            "generator": relpath(generator_path),
            "generator_sha256": sha256_file(generator_path),
            "source_parent_artifact": relpath(V5_ROOT),
            "source_parent_manifest_sha256": sha256_file(V5_MANIFEST),
            "source_parent_plot_manifest_sha256": sha256_file(V5_PLOT_MANIFEST),
            "source_parent_figure_layer_manifest_sha256": sha256_file(V5_FIGURE_LAYER_MANIFEST),
            "manuscript_eligible": False,
            "canonical_data_modified": False,
            "production_write": False,
            "solver_called": False,
            "display_style_update": style_update,
            "display_style_manifest": relpath(V6_STYLE_MANIFEST),
            "display_style_manifest_sha256": sha256_file(V6_STYLE_MANIFEST),
            "rendering_semantics": str(v5_plot.get("rendering_semantics", ""))
            + "; v6 changes only tick, legend-size, and endpoint legend rendering",
            "figures": assets,
        }
    )
    write_json(V6_PLOT_MANIFEST, plot_manifest)

    (V6_ROOT / "README.md").write_text(render_readme(v5_package, len(assets)), encoding="utf-8")
    package_manifest = dict(v5_package)
    package_manifest.update(
        {
            "schema": "phase_guided_transport_publication_clean_manifest_v6",
            "generated_at": generated_at,
            "base_git_commit": git_head(),
            "generator": relpath(generator_path),
            "generator_sha256": sha256_file(generator_path),
            "source_parent_artifact": relpath(V5_ROOT),
            "source_parent_manifest_sha256": sha256_file(V5_MANIFEST),
            "source_parent_plot_manifest_sha256": sha256_file(V5_PLOT_MANIFEST),
            "source_parent_figure_layer_manifest_sha256": sha256_file(V5_FIGURE_LAYER_MANIFEST),
            "status": "derived_author_review_required",
            "manuscript_eligible": False,
            "canonical_data_modified": False,
            "production_write": False,
            "solver_called": False,
            "display_style_update": style_update,
            "display_style_manifest": relpath(V6_STYLE_MANIFEST),
            "display_style_manifest_sha256": sha256_file(V6_STYLE_MANIFEST),
            "derived_counts": {
                **dict(v5_package.get("derived_counts", {})),
                "publication_figure_count": len(assets),
            },
            "known_boundaries": [
                *list(v5_package.get("known_boundaries", [])),
                "v6 is a figure-style and endpoint-label review derivative of v5; v5 remains current",
                "v6 does not assert a new numerical convergence result or physical-model change",
            ],
            "outputs": [],
        }
    )
    package_manifest["outputs"] = [
        {"path": relpath(path), "sha256": sha256_file(path), "bytes": path.stat().st_size}
        for path in package_output_paths()
    ]
    write_json(V6_MANIFEST, package_manifest)

    print(
        json.dumps(
            {
                "output": relpath(V6_ROOT),
                "manifest": relpath(V6_MANIFEST),
                "publication_figures": len(assets),
                "mode_counts": mode_counts,
                "solver_called": False,
                "canonical_data_modified": False,
                "manuscript_eligible": False,
                "v5_parent_manifest_sha256": sha256_file(V5_MANIFEST),
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
