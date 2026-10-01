#!/usr/bin/env python3
"""Render diagnostic heatmaps from charged-GBU screening artifacts.

This script is deliberately solver-free.  It consumes the immutable CSV/JSON
files emitted by ``run_charged_gbu_contour_scan.jl`` and writes PNG review
figures plus a provenance manifest.  Failed points remain masked; no
interpolation, smoothing, or zero-filling is performed.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
import tomllib
from typing import Any, Iterable


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
REQUIRED_COLUMNS = {
    "T_MeV",
    "muB_MeV",
    "status",
    "pi_plus_density_inv_fm3",
    "K_plus_density_inv_fm3",
    "Kplus_over_pi_plus",
    "pi_minus_density_inv_fm3",
    "K_minus_density_inv_fm3",
    "Kminus_over_pi_minus",
    "pi_plus_passed",
    "K_plus_passed",
    "pi_minus_passed",
    "K_minus_passed",
}
VALUE_COLUMNS = (
    "pi_plus_density_inv_fm3",
    "K_plus_density_inv_fm3",
    "Kplus_over_pi_plus",
    "pi_minus_density_inv_fm3",
    "K_minus_density_inv_fm3",
    "Kminus_over_pi_minus",
)
SCHEMA_VERSION = "charged_gbu_contour_screening_plot_manifest_v1"
DEFAULT_FREEZEOUT_PROFILE = REPOSITORY_ROOT / "config" / "physics" / "freezeout" / "default.toml"
DEFAULT_PHASE_REFERENCE_ROOT = (
    REPOSITORY_ROOT
    / "data"
    / "reference"
    / "pnjl"
    / "issue130_phase_reference_v2"
    / "accepted"
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _finite_float(value: str, *, field: str, path: Path) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{path}: {field} is not numeric: {value!r}") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{path}: {field} is not finite: {value!r}")
    return parsed


def _optional_float(value: str, *, field: str, path: Path) -> float | None:
    if value is None or str(value).strip().lower() in {"", "nan", "null", "none"}:
        return None
    return _finite_float(value, field=field, path=path)


def _bool(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _canonical(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def _repo_head() -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPOSITORY_ROOT, text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _discover(input_root: Path) -> tuple[list[Path], list[Path]]:
    if not input_root.is_dir():
        raise ValueError(f"input root does not exist: {input_root}")
    csv_paths = sorted(input_root.rglob("contour_points.csv"))
    manifest_paths = sorted(input_root.rglob("manifest.json"))
    if not csv_paths:
        raise ValueError(f"no contour_points.csv found below {input_root}")
    if len(csv_paths) != len(manifest_paths):
        raise ValueError(
            f"expected one manifest per CSV, found {len(csv_paths)} CSV and "
            f"{len(manifest_paths)} manifest files"
        )
    return csv_paths, manifest_paths


def _read_rows(csv_paths: Iterable[Path]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in csv_paths:
        with path.open("r", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            fields = set(reader.fieldnames or ())
            missing = sorted(REQUIRED_COLUMNS - fields)
            if missing:
                raise ValueError(f"{path}: missing CSV columns: {', '.join(missing)}")
            for row_number, raw in enumerate(reader, start=2):
                T = _finite_float(raw["T_MeV"], field="T_MeV", path=path)
                muB = _finite_float(raw["muB_MeV"], field="muB_MeV", path=path)
                row: dict[str, Any] = dict(raw)
                row["T_MeV"] = T
                row["muB_MeV"] = muB
                row["_source"] = str(path)
                row["_row_number"] = row_number
                row["status"] = str(raw["status"] or "").strip()
                for field in VALUE_COLUMNS:
                    row[field] = _optional_float(raw.get(field), field=field, path=path)
                for field in ("pi_plus_passed", "K_plus_passed", "pi_minus_passed", "K_minus_passed"):
                    row[field] = _bool(raw.get(field, "false"))
                rows.append(row)
    return rows


def load_dataset(input_root: str | Path) -> dict[str, Any]:
    """Validate and load a complete regular screening grid."""

    root = Path(input_root).resolve()
    csv_paths, manifest_paths = _discover(root)
    manifests: list[dict[str, Any]] = []
    for path in manifest_paths:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("schema") != "charged_gbu_contour_scan_v2":
            raise ValueError(f"unsupported scan manifest schema in {path}")
        if not isinstance(payload.get("source_hashes"), dict):
            raise ValueError(f"{path}: source_hashes must be an object")
        manifests.append(payload)

    first = manifests[0]
    invariant_fields = ("T_grid", "muB_grid", "channels", "settings", "config", "git_head")
    for field in invariant_fields:
        expected = _canonical(first.get(field))
        if any(_canonical(manifest.get(field)) != expected for manifest in manifests[1:]):
            raise ValueError(f"scan manifests disagree on {field}")
    expected_source_hashes = _canonical(first["source_hashes"])
    if any(_canonical(manifest["source_hashes"]) != expected_source_hashes for manifest in manifests[1:]):
        raise ValueError("scan manifests disagree on source hashes")

    rows = _read_rows(csv_paths)
    key_to_row: dict[tuple[float, float], dict[str, Any]] = {}
    for row in rows:
        key = (row["T_MeV"], row["muB_MeV"])
        if key in key_to_row:
            raise ValueError(f"duplicate screening key: {key}")
        key_to_row[key] = row

    T_grid = [float(value) for value in first["T_grid"]]
    muB_grid = [float(value) for value in first["muB_grid"]]
    expected_keys = {(T, muB) for T in T_grid for muB in muB_grid}
    actual_keys = set(key_to_row)
    if actual_keys != expected_keys:
        missing = sorted(expected_keys - actual_keys)
        extra = sorted(actual_keys - expected_keys)
        raise ValueError(f"grid keys differ; missing={missing[:5]} extra={extra[:5]}")
    expected_count = sum(int(manifest.get("point_count", -1)) for manifest in manifests)
    if expected_count != len(rows):
        raise ValueError(f"manifest point count {expected_count} != CSV rows {len(rows)}")

    status_counts: dict[str, int] = {}
    for row in rows:
        status_counts[row["status"]] = status_counts.get(row["status"], 0) + 1
        if row["status"] == "screened":
            for field in VALUE_COLUMNS:
                if row[field] is None:
                    raise ValueError(f"screened row has missing {field}: {row}")

    return {
        "input_root": root,
        "csv_paths": csv_paths,
        "manifest_paths": manifest_paths,
        "manifests": manifests,
        "rows": sorted(rows, key=lambda row: (row["T_MeV"], row["muB_MeV"])),
        "T_grid": sorted(T_grid),
        "muB_grid": sorted(muB_grid),
        "channels": list(first["channels"]),
        "settings": first["settings"],
        "config": first["config"],
        "git_head": first["git_head"],
        "source_hashes": first["source_hashes"],
        "status_counts": status_counts,
    }


def matrix(dataset: dict[str, Any], field: str) -> list[list[float | None]]:
    """Return an exact grid matrix, retaining failed values as None."""

    values = {(row["T_MeV"], row["muB_MeV"]): row for row in dataset["rows"]}
    return [
        [
            values[(T, muB)][field] if values[(T, muB)]["status"] == "screened" else None
            for muB in dataset["muB_grid"]
        ]
        for T in dataset["T_grid"]
    ]


def mask_matrix(dataset: dict[str, Any], side: str | None = None) -> list[list[int]]:
    """Return categorical masks without converting failed points to zero."""

    values = {(row["T_MeV"], row["muB_MeV"]): row for row in dataset["rows"]}
    result: list[list[int]] = []
    for T in dataset["T_grid"]:
        current: list[int] = []
        for muB in dataset["muB_grid"]:
            row = values[(T, muB)]
            if row["status"] != "screened":
                current.append(2)
                continue
            if side is None:
                current.append(0)
            elif side == "plus":
                current.append(0 if row["pi_plus_passed"] and row["K_plus_passed"] else 1)
            elif side == "minus":
                current.append(0 if row["pi_minus_passed"] and row["K_minus_passed"] else 1)
            else:
                raise ValueError(f"unknown mask side: {side}")
        result.append(current)
    return result


def _edges(values: list[float]) -> list[float]:
    if len(values) == 1:
        return [values[0] - 0.5, values[0] + 0.5]
    mids = [(left + right) / 2.0 for left, right in zip(values[:-1], values[1:])]
    return [2 * values[0] - mids[0], *mids, 2 * values[-1] - mids[-1]]


def _reference_points(reference_lines: dict[str, Any]) -> list[tuple[float, float]]:
    return [
        (float(point["muB_MeV"]), float(point["T_MeV"]))
        for line in reference_lines.values()
        for point in line["points"]
    ]


def load_reference_lines(
    *,
    freezeout_profile: Path = DEFAULT_FREEZEOUT_PROFILE,
    phase_reference_root: Path = DEFAULT_PHASE_REFERENCE_ROOT,
    xi: float = 0.0,
    sqrt_s_min_GeV: float = 3.0,
    sqrt_s_max_GeV: float = 200.0,
    freezeout_points: int = 256,
) -> dict[str, Any]:
    """Load solver-free reference curves with explicit coordinate provenance.

    The phase-reference tables use ``mu_MeV`` as the quark chemical potential;
    the contour axes use ``muB_MeV`` and therefore apply ``muB=3*mu_q``.  The
    accepted phase tables are an equal-flavor historical reference, not a
    solved BQS phase boundary; that distinction is retained in the manifest and
    plot labels.
    """

    if freezeout_points < 2:
        raise ValueError("freezeout_points must be at least 2")
    if not freezeout_profile.is_file():
        raise ValueError(f"freeze-out profile does not exist: {freezeout_profile}")
    cfg = tomllib.loads(freezeout_profile.read_text(encoding="utf-8"))
    coeff = cfg.get("coefficients", {})
    required = ("a_GeV", "b_GeV_inv1", "c_GeV_inv3", "d_GeV", "e_GeV_inv1")
    if any(key not in coeff for key in required):
        raise ValueError(f"freeze-out profile missing coefficients: {freezeout_profile}")
    a, b, c = (float(coeff[key]) for key in required[:3])
    d, e = (float(coeff[key]) for key in required[3:])
    if not all(math.isfinite(value) for value in (a, b, c, d, e)):
        raise ValueError(f"freeze-out profile has non-finite coefficients: {freezeout_profile}")
    if sqrt_s_min_GeV <= 0 or sqrt_s_max_GeV < sqrt_s_min_GeV:
        raise ValueError("invalid freeze-out sqrt(s) range")
    freezeout = []
    for index in range(freezeout_points):
        sqrt_s = sqrt_s_min_GeV + (sqrt_s_max_GeV - sqrt_s_min_GeV) * index / (freezeout_points - 1)
        muB_GeV = d / (1.0 + e * sqrt_s)
        T_GeV = a - b * muB_GeV**2 - c * muB_GeV**4
        if math.isfinite(muB_GeV) and math.isfinite(T_GeV) and T_GeV > 0:
            freezeout.append({"sqrt_s_NN_GeV": sqrt_s, "muB_MeV": 1000.0 * muB_GeV, "T_MeV": 1000.0 * T_GeV})
    lines: dict[str, Any] = {
        "freezeout": {
            "label": "chemical freeze-out (default profile; quark-only BQS target)",
            "style": {"color": "white", "linestyle": "-", "linewidth": 1.35},
            "points": freezeout,
            "source": str(freezeout_profile),
            "source_sha256": sha256_file(freezeout_profile),
            "coordinate_convention": "muB_MeV; Cleymans-like profile evaluated from config",
        }
    }
    phase_files = {
        "crossover": "crossover_surface_accepted_phase_map_v1.csv",
        "first_order_maxwell": "maxwell_surface_accepted_phase_map_v1.csv",
    }
    phase_styles = {
        "crossover": {"color": "cyan", "linestyle": "--", "linewidth": 1.0},
        "first_order_maxwell": {"color": "yellow", "linestyle": ":", "linewidth": 1.15},
    }
    for key, filename in phase_files.items():
        path = phase_reference_root / "tables" / filename
        if not path.is_file():
            raise ValueError(f"phase reference table does not exist: {path}")
        points = []
        with path.open("r", newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                try:
                    row_xi = float(row["xi"])
                    temperature = float(row["T_MeV"])
                    mu_q = float(row["mu_MeV"])
                except (KeyError, TypeError, ValueError) as exc:
                    raise ValueError(f"invalid phase reference row in {path}: {row}") from exc
                if abs(row_xi - xi) > 1.0e-9:
                    continue
                if not all(math.isfinite(value) for value in (temperature, mu_q)):
                    raise ValueError(f"non-finite phase reference row in {path}: {row}")
                points.append({"muB_MeV": 3.0 * mu_q, "T_MeV": temperature})
        points.sort(key=lambda point: (point["muB_MeV"], point["T_MeV"]))
        lines[key] = {
            "label": f"{key.replace('_', ' ')} reference (equal-flavor; muB=3 muq)",
            "style": phase_styles[key],
            "points": points,
            "source": str(path),
            "source_sha256": sha256_file(path),
            "xi": xi,
            "coordinate_convention": "input mu_MeV interpreted as mu_q; plotted muB_MeV=3*mu_q",
            "background_scope": "equal-flavor historical reference; not current BQS phase boundary",
        }
    return lines


def _plot_reference_lines(ax: Any, reference_lines: dict[str, Any]) -> None:
    import matplotlib.patheffects as path_effects

    for line in reference_lines.values():
        points = line["points"]
        if not points:
            continue
        plotted = ax.plot(
            [point["muB_MeV"] for point in points],
            [point["T_MeV"] for point in points],
            label=line["label"],
            **line["style"],
            zorder=4,
        )
        # Keep the white freeze-out guide visible both over the heatmap and in
        # the blank region where a reference line extends beyond the scan.
        if line.get("style", {}).get("color") == "white":
            plotted[0].set_path_effects([
                path_effects.Stroke(linewidth=float(line["style"].get("linewidth", 1.0)) + 1.6, foreground="black"),
                path_effects.Normal(),
            ])


def _axis_limits(dataset: dict[str, Any], reference_lines: dict[str, Any]) -> tuple[tuple[float, float], tuple[float, float]]:
    x_edges = _edges(dataset["muB_grid"])
    y_edges = _edges(dataset["T_grid"])
    points = _reference_points(reference_lines)
    x_values = [x_edges[0], x_edges[-1], *[point[0] for point in points]]
    y_values = [y_edges[0], y_edges[-1], *[point[1] for point in points]]
    return (min(x_values), max(x_values)), (min(y_values), max(y_values))


def _render_figures(dataset: dict[str, Any], output_dir: Path, reference_lines: dict[str, Any]) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    import numpy as np

    matplotlib.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 9.0,
            "axes.titlesize": 9.0,
            "axes.labelsize": 9.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.08,
        }
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    x_edges = _edges(dataset["muB_grid"])
    y_edges = _edges(dataset["T_grid"])
    x_limits, y_limits = _axis_limits(dataset, reference_lines)
    figures: list[Path] = []

    def heatmap(field: str, filename: str, title: str, colorbar: str, cmap: str) -> None:
        values = np.array(
            [[np.nan if value is None else value for value in row] for row in matrix(dataset, field)],
            dtype=float,
        )
        fig, ax = plt.subplots(figsize=(6.75, 4.6))
        color_map = plt.get_cmap(cmap).copy()
        color_map.set_bad("#d9d9d9")
        image = ax.pcolormesh(
            x_edges,
            y_edges,
            np.ma.masked_invalid(values),
            shading="flat",
            cmap=color_map,
        )
        ax.set_xlabel(r"$\mu_B$ (MeV)")
        ax.set_ylabel(r"$T$ (MeV)")
        ax.set_title(title + "\nfixed quark-only BQS; diagnostic screening")
        _plot_reference_lines(ax, reference_lines)
        ax.axvline(x_edges[-1], color="white", linestyle="-.", linewidth=0.8, alpha=0.8,
                   label="screening domain boundary" if field == "Kplus_over_pi_plus" else None)
        ax.set_xlim(*x_limits)
        ax.set_ylim(*y_limits)
        if field == "Kplus_over_pi_plus":
            ax.legend(loc="best", fontsize=6.5, framealpha=0.75)
        fig.colorbar(image, ax=ax, label=colorbar)
        path = output_dir / filename
        fig.savefig(path)
        plt.close(fig)
        figures.append(path)

    heatmap("pi_plus_density_inv_fm3", "n_pi_plus_screening.png", r"Charged GBU screening: $n_{\pi^+}$", r"$n_{\pi^+}$ (fm$^{-3}$)", "viridis")
    heatmap("K_plus_density_inv_fm3", "n_K_plus_screening.png", "Charged GBU screening: $n_{K^+}$", r"$n_{K^+}$ (fm$^{-3}$)", "plasma")
    heatmap("Kplus_over_pi_plus", "Kplus_over_pi_plus_screening.png", r"Charged GBU screening: $K^+/\pi^+$", r"$K^+/\pi^+$", "magma")

    masks = [
        ("status", mask_matrix(dataset), ["screened", "gate failed", "other failure"]),
        ("plus", mask_matrix(dataset, "plus"), ["plus channels pass", "plus channel failure", "point failure"]),
        ("minus", mask_matrix(dataset, "minus"), ["minus channels pass", "minus channel failure", "point failure"]),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(9.8, 3.6), sharex=True, sharey=True)
    mask_cmap = ListedColormap(["#2ca25f", "#de2d26", "#756bb1"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], mask_cmap.N)
    for ax, (name, values, labels) in zip(axes, masks):
        image = ax.pcolormesh(
            x_edges,
            y_edges,
            np.asarray(values, dtype=float),
            shading="flat",
            cmap=mask_cmap,
            norm=norm,
        )
        ax.set_title(name)
        ax.set_xlabel(r"$\mu_B$ (MeV)")
        _plot_reference_lines(ax, reference_lines)
        ax.set_xlim(*x_limits)
        ax.set_ylim(*y_limits)
    axes[0].set_ylabel(r"$T$ (MeV)")
    cbar = fig.colorbar(image, ax=axes.ravel().tolist(), ticks=[0, 1, 2], shrink=0.88)
    cbar.ax.set_yticklabels(["pass", "channel fail", "point fail"])
    fig.suptitle("Charged GBU screening masks (failed points are not zero-filled)\nfixed quark-only BQS; reference phase lines are not BQS solves")
    mask_path = output_dir / "screening_failure_masks.png"
    fig.savefig(mask_path)
    plt.close(fig)
    figures.append(mask_path)
    return figures


def _write_merged_csv(dataset: dict[str, Any], path: Path) -> None:
    fieldnames = [
        "T_MeV",
        "muB_MeV",
        "status",
        *VALUE_COLUMNS,
        "pi_plus_passed",
        "K_plus_passed",
        "pi_minus_passed",
        "K_minus_passed",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in dataset["rows"]:
            output = {field: row.get(field, "") for field in fieldnames}
            writer.writerow(output)


def _manifest(dataset: dict[str, Any], output_dir: Path, figures: list[Path], merged_csv: Path, source_run_id: str | None, git_sha: str | None, reference_lines: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "diagnostic_only": True,
        "production_default": False,
        "delivery_stage": "png_review",
        "manuscript_eligible": False,
        "current_publication_layer": False,
        "vector_delivery_pending": True,
        "source_run_id": source_run_id,
        "source_git_sha": git_sha,
        "repository_head_at_plot": _repo_head(),
        "renderer": {"python": platform.python_version(), "matplotlib": _matplotlib_version()},
        "input_root": str(dataset["input_root"]),
        "input_manifests": [
            {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256_file(path)}
            for path in dataset["manifest_paths"]
        ],
        "input_csvs": [
            {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256_file(path)}
            for path in dataset["csv_paths"]
        ],
        "scan_source_hashes": dataset["source_hashes"],
        "config": dataset["config"],
        "channels": dataset["channels"],
        "settings": dataset["settings"],
        "grid": {"T_MeV": dataset["T_grid"], "muB_MeV": dataset["muB_grid"]},
        "row_count": len(dataset["rows"]),
        "status_counts": dataset["status_counts"],
        "selection_rule": "exact screened rows from scan CSV; failed rows remain masked",
        "interpolation_policy": "none",
        "missing_value_policy": "mask; never zero-fill",
        "units": {"T": "MeV", "muB": "MeV", "density": "fm^-3"},
        "background_contract": {
            "scope": "fixed quark-only BQS GBU partial yield",
            "rho_Q_over_rho_B": 0.4,
            "rho_S_fm3": 0.0,
            "meson_feedback": False,
            "diagnostic_only": True,
        },
        "reference_lines": {
            key: {field: value for field, value in line.items() if field != "points"}
            | {"point_count": len(line["points"])}
            for key, line in reference_lines.items()
        },
        "reference_coordinate_note": "phase tables use mu_q=mu_MeV; plots convert to muB_MeV=3*mu_q; phase references are equal-flavor historical guides, not current BQS boundaries",
        "screening_domain": {
            "T_MeV": [min(dataset["T_grid"]), max(dataset["T_grid"])],
            "muB_MeV": [min(dataset["muB_grid"]), max(dataset["muB_grid"])],
            "reference_lines_may_extend_beyond_domain": True,
        },
        "merged_csv": {"path": str(merged_csv), "sha256": sha256_file(merged_csv)},
        "figures": [
            {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256_file(path), "format": "png"}
            for path in figures
        ],
    }


def _matplotlib_version() -> str:
    try:
        import matplotlib

        return str(matplotlib.__version__)
    except ImportError:
        return "unavailable"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--source-run-id", default=None)
    parser.add_argument("--git-sha", default=None)
    parser.add_argument("--freezeout-profile", type=Path, default=DEFAULT_FREEZEOUT_PROFILE)
    parser.add_argument("--phase-reference-root", type=Path, default=DEFAULT_PHASE_REFERENCE_ROOT)
    parser.add_argument("--reference-xi", type=float, default=0.0)
    parser.add_argument("--no-reference-lines", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    input_root = args.input_root.resolve()
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"refusing to overwrite non-empty output directory: {output_dir}")
    dataset = load_dataset(input_root)
    output_dir.mkdir(parents=True, exist_ok=True)
    reference_lines = {} if args.no_reference_lines else load_reference_lines(
        freezeout_profile=args.freezeout_profile.resolve(),
        phase_reference_root=args.phase_reference_root.resolve(),
        xi=args.reference_xi,
    )
    merged_csv = output_dir / "contour_points_merged.csv"
    _write_merged_csv(dataset, merged_csv)
    figures = _render_figures(dataset, output_dir, reference_lines)
    manifest = _manifest(dataset, output_dir, figures, merged_csv, args.source_run_id, args.git_sha, reference_lines)
    manifest_path = output_dir / "plot_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "complete", "row_count": len(dataset["rows"]), "figures": [path.name for path in figures], "manifest": str(manifest_path)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
