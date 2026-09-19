#!/usr/bin/env python3
"""Mirror publication_clean_v4 publication figures to the public figure layer."""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import shutil
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
ANALYSIS_ROOT = (
    ROOT
    / "docs"
    / "analysis"
    / "relaxtime"
    / "phase_guided_transport"
    / "phase_guided_transport_publication_clean_v4"
)
ANALYSIS_FIGURE_ROOT = ANALYSIS_ROOT / "figures"
ANALYSIS_PLOT_MANIFEST = ANALYSIS_FIGURE_ROOT / "plot_manifest.json"
ANALYSIS_PACKAGE_MANIFEST = ANALYSIS_ROOT / "manifest.json"
PUBLIC_ROOT = (
    ROOT
    / "data"
    / "outputs"
    / "figures"
    / "relaxtime"
    / "transport"
    / "phase_guided"
    / "publication_clean_v4"
)
PUBLIC_PLOT_MANIFEST = PUBLIC_ROOT / "plot_manifest.json"
FIGURE_LAYER_ROOT = (
    ROOT
    / "docs"
    / "analysis"
    / "relaxtime"
    / "phase_guided_transport"
    / "phase_guided_transport_publication_clean_figure_layer_v4"
)
FIGURE_LAYER_MANIFEST = FIGURE_LAYER_ROOT / "figure_layer_manifest.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relpath(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def load_manifests() -> tuple[dict[str, Any], dict[str, Any]]:
    if not ANALYSIS_PLOT_MANIFEST.is_file() or not ANALYSIS_PACKAGE_MANIFEST.is_file():
        raise FileNotFoundError("build publication_clean_v4 before syncing its figure layer")
    plot_manifest = json.loads(ANALYSIS_PLOT_MANIFEST.read_text(encoding="utf-8"))
    package_manifest = json.loads(ANALYSIS_PACKAGE_MANIFEST.read_text(encoding="utf-8"))
    if plot_manifest.get("schema") != "phase_guided_transport_publication_clean_plot_manifest_v4":
        raise ValueError("analysis plot manifest is not publication_clean_v4")
    if plot_manifest.get("marker_contract") != "phase_endpoint_gap_v2":
        raise ValueError("analysis plot manifest has an unexpected marker contract")
    if plot_manifest.get("manuscript_eligible") is not False:
        raise ValueError("v4 analysis layer must remain manuscript_eligible=false")
    figures = plot_manifest.get("figures", [])
    if len(figures) != 72:
        raise ValueError(f"v4 publication figure count is {len(figures)}, expected 72")
    if package_manifest.get("schema") != "phase_guided_transport_publication_clean_manifest_v4":
        raise ValueError("package manifest is not publication_clean_v4")
    if package_manifest.get("manuscript_eligible") is not False:
        raise ValueError("v4 package manifest must remain manuscript_eligible=false")
    return plot_manifest, package_manifest


def public_path_for(source_path: Path) -> Path:
    relative = source_path.resolve().relative_to(ANALYSIS_FIGURE_ROOT.resolve())
    mode_alias = {
        "mode_a": "mode_a_fixed_muB_phase_scaled",
        "mode_b": "mode_b_fixed_T_sparse_muB",
    }
    if not relative.parts or relative.parts[0] not in mode_alias:
        raise ValueError(f"unexpected v4 analysis figure layout: {relative}")
    return PUBLIC_ROOT / mode_alias[relative.parts[0]] / Path(*relative.parts[1:])


def sync_figures(plot_manifest: dict[str, Any]) -> list[dict[str, Any]]:
    assets: list[dict[str, Any]] = []
    expected: set[Path] = set()
    for asset in plot_manifest["figures"]:
        source = ROOT / Path(asset["path"])
        if not source.resolve().is_relative_to(ANALYSIS_FIGURE_ROOT.resolve()):
            raise ValueError(f"analysis asset is outside v4 figure root: {asset['path']}")
        if source.suffix.lower() != ".png" or not source.is_file():
            raise ValueError(f"missing or non-PNG v4 publication asset: {source}")
        if sha256_file(source) != asset.get("sha256"):
            raise ValueError(f"v4 analysis manifest hash mismatch: {asset['path']}")
        target = public_path_for(source)
        expected.add(target.resolve())
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        if sha256_file(target) != asset["sha256"]:
            raise ValueError(f"byte-preserving v4 copy failed: {target}")
        assets.append(
            {
                **{key: value for key, value in asset.items() if key not in {"path", "bytes", "sha256"}},
                "path": relpath(target),
                "bytes": target.stat().st_size,
                "sha256": sha256_file(target),
            }
        )
    existing = {path.resolve() for path in PUBLIC_ROOT.rglob("*.png") if path.is_file()}
    stale = sorted(existing - expected)
    if stale:
        raise ValueError("v4 public figure layer contains stale PNGs: " + ", ".join(str(path) for path in stale))
    return assets


def write_layer_readme() -> None:
    FIGURE_LAYER_ROOT.mkdir(parents=True, exist_ok=True)
    (FIGURE_LAYER_ROOT / "README.md").write_text(
        """# RS publication_clean_v4 figure-layer review candidate

This directory records the byte-preserving publication figure mirror for
publication_clean_v4.

The v4 analysis package inherits v3, adjusts the three mode-B composite curves
at T=200, muB=0, xi=0.36, and applies one one-sided left-branch display
adjustment to mode-A tau_sbar at the audited first-order left endpoint.
The first-order gap remains split. Raw results, production registries, and
solver outputs are unchanged.

The layer is for author review only and remains manuscript_eligible=false.
The mechanism review is retained as diagnostic context; it is not a
propagator regularization or a production convergence certificate.

Reproduction:

    python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v4.py
    python scripts/analysis/relaxtime/sync_phase_guided_publication_clean_figure_layer_v4.py
""",
        encoding="utf-8",
    )


def main() -> None:
    if PUBLIC_ROOT.exists() or FIGURE_LAYER_ROOT.exists():
        raise FileExistsError("refusing to overwrite an existing publication_clean_v4 figure layer")
    analysis_manifest, package_manifest = load_manifests()
    assets = sync_figures(analysis_manifest)
    mode_counts = {"mode_a": 0, "mode_b": 0}
    for asset in assets:
        if "/mode_a_fixed_muB_phase_scaled/" in asset["path"]:
            mode_counts["mode_a"] += 1
        elif "/mode_b_fixed_T_sparse_muB/" in asset["path"]:
            mode_counts["mode_b"] += 1
        else:
            raise ValueError(f"unable to classify v4 asset: {asset['path']}")
    if mode_counts != {"mode_a": 36, "mode_b": 36}:
        raise ValueError(f"unexpected v4 mode counts: {mode_counts}")

    public_manifest = dict(analysis_manifest)
    public_manifest["figure_layer_role"] = "publication_clean_display_mirror"
    public_manifest["source_analysis_plot_manifest"] = relpath(ANALYSIS_PLOT_MANIFEST)
    public_manifest["source_analysis_plot_manifest_sha256"] = sha256_file(ANALYSIS_PLOT_MANIFEST)
    public_manifest["figures"] = assets
    public_manifest["rendering_semantics"] = (
        str(public_manifest.get("rendering_semantics", ""))
        + "; public figure layer is a byte-preserving mirror of the v4 analysis figures"
    )
    write_json(PUBLIC_PLOT_MANIFEST, public_manifest)

    layer_manifest = {
        "schema": "phase_guided_transport_publication_clean_figure_layer_manifest_v4",
        "figure_layer_role": "publication_clean_display_mirror",
        "marker_contract": "phase_endpoint_gap_v2",
        "status": "derived_author_review_required",
        "manuscript_eligible": False,
        "solver_called": False,
        "source_analysis_package": relpath(ANALYSIS_ROOT),
        "source_analysis_package_manifest": relpath(ANALYSIS_PACKAGE_MANIFEST),
        "source_analysis_package_manifest_sha256": sha256_file(ANALYSIS_PACKAGE_MANIFEST),
        "source_plot_manifest": relpath(ANALYSIS_PLOT_MANIFEST),
        "source_plot_manifest_sha256": sha256_file(ANALYSIS_PLOT_MANIFEST),
        "target_plot_manifest": relpath(PUBLIC_PLOT_MANIFEST),
        "target_plot_manifest_sha256": sha256_file(PUBLIC_PLOT_MANIFEST),
        "sync_generator": relpath(Path(__file__).resolve()),
        "sync_generator_sha256": sha256_file(Path(__file__).resolve()),
        "source_parent_artifact": package_manifest["source_parent_artifact"],
        "source_parent_manifest_sha256": package_manifest["source_parent_manifest_sha256"],
        "adjustment_recipe": package_manifest["adjustment_recipe"],
        "adjustment_recipe_sha256": package_manifest["adjustment_recipe_sha256"],
        "mechanism_review_artifact": package_manifest["mechanism_review_artifact"],
        "calculation_sha": package_manifest["calculation_sha"],
        "workflow_head_sha": package_manifest["workflow_head_sha"],
        "source_case": package_manifest["case"],
        "mode_counts": mode_counts,
        "figure_count": len(assets),
        "content_policy": "byte_preserve",
        "source_package_unchanged": True,
        "raw_results_unchanged": True,
        "production_registry_unchanged": True,
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "assets": assets,
    }
    write_json(FIGURE_LAYER_MANIFEST, layer_manifest)
    write_layer_readme()
    print(
        json.dumps(
            {
                "public_root": relpath(PUBLIC_ROOT),
                "figure_layer_manifest": relpath(FIGURE_LAYER_MANIFEST),
                "figures": len(assets),
                "mode_counts": mode_counts,
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
