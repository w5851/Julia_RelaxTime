#!/usr/bin/env python3
"""Mirror the solver-free publication-clean analysis figures to the public layer.

The build script owns the analysis package under ``docs/analysis``.  This
small, explicit sync step owns only the byte-preserving figure mirror under
``data/outputs/figures`` and its two manifests.  It refuses to copy anything
whose source hash disagrees with the analysis plot manifest and never reads or
writes a raw result, registry, or production figure.
"""

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
    / "phase_guided_transport_publication_clean_v1"
)
ANALYSIS_FIGURE_ROOT = ANALYSIS_ROOT / "figures"
PUBLIC_ROOT = (
    ROOT
    / "data"
    / "outputs"
    / "figures"
    / "relaxtime"
    / "transport"
    / "phase_guided"
    / "publication_clean_v1"
)
FIGURE_LAYER_ROOT = (
    ROOT
    / "docs"
    / "analysis"
    / "relaxtime"
    / "phase_guided_transport"
    / "phase_guided_transport_publication_clean_figure_layer_v1"
)
ANALYSIS_PLOT_MANIFEST = ANALYSIS_FIGURE_ROOT / "plot_manifest.json"
ANALYSIS_PACKAGE_MANIFEST = ANALYSIS_ROOT / "manifest.json"
PUBLIC_PLOT_MANIFEST = PUBLIC_ROOT / "plot_manifest.json"
FIGURE_LAYER_MANIFEST = FIGURE_LAYER_ROOT / "figure_layer_manifest.json"
ANALYSIS_FIGURE_REL = ANALYSIS_FIGURE_ROOT.relative_to(ROOT).as_posix()
PUBLIC_FIGURE_REL = PUBLIC_ROOT.relative_to(ROOT).as_posix()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def relpath(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def load_manifest() -> dict[str, Any]:
    if not ANALYSIS_PLOT_MANIFEST.exists():
        raise FileNotFoundError(f"missing analysis plot manifest: {ANALYSIS_PLOT_MANIFEST}")
    payload = json.loads(ANALYSIS_PLOT_MANIFEST.read_text(encoding="utf-8"))
    figures = payload.get("figures", [])
    if not figures:
        raise ValueError("analysis plot manifest has no figures")
    if payload.get("manuscript_eligible") is not False:
        raise ValueError("analysis layer must remain manuscript_eligible=false")
    return payload


def load_package_manifest() -> dict[str, Any]:
    if not ANALYSIS_PACKAGE_MANIFEST.exists():
        raise FileNotFoundError(
            f"missing analysis package manifest: {ANALYSIS_PACKAGE_MANIFEST}"
        )
    return json.loads(ANALYSIS_PACKAGE_MANIFEST.read_text(encoding="utf-8"))


def public_path_for(source_path: Path) -> Path:
    relative = source_path.resolve().relative_to(ANALYSIS_FIGURE_ROOT.resolve())
    mode_alias = {
        "mode_a": "mode_a_fixed_muB_phase_scaled",
        "mode_b": "mode_b_fixed_T_sparse_muB",
    }
    if not relative.parts or relative.parts[0] not in mode_alias:
        raise ValueError(f"unexpected analysis figure layout: {relative}")
    return PUBLIC_ROOT / mode_alias[relative.parts[0]] / Path(*relative.parts[1:])


def sync_figures(analysis_manifest: dict[str, Any]) -> list[dict[str, Any]]:
    assets: list[dict[str, Any]] = []
    expected_public_paths: set[Path] = set()
    for asset in analysis_manifest["figures"]:
        source = ROOT / Path(asset["path"])
        if source.parent != ANALYSIS_FIGURE_ROOT and not source.is_relative_to(ANALYSIS_FIGURE_ROOT):
            raise ValueError(f"analysis asset is outside figure root: {asset['path']}")
        if source.suffix.lower() != ".png":
            raise ValueError(f"publication-clean mirror only accepts PNG assets: {asset['path']}")
        if not source.exists():
            raise FileNotFoundError(f"missing analysis figure: {source}")
        source_hash = sha256_file(source)
        if source_hash != asset.get("sha256"):
            raise ValueError(f"analysis manifest hash mismatch: {asset['path']}")
        target = public_path_for(source)
        expected_public_paths.add(target.resolve())
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        target_hash = sha256_file(target)
        if target_hash != source_hash:
            raise ValueError(f"byte-preserving copy failed: {target}")
        assets.append(
            {
                "path": relpath(target),
                "bytes": target.stat().st_size,
                "sha256": target_hash,
            }
        )

    existing_pngs = {
        path.resolve()
        for path in PUBLIC_ROOT.rglob("*.png")
        if path.is_file()
    }
    stale = sorted(existing_pngs - expected_public_paths)
    if stale:
        raise ValueError(
            "public figure layer contains PNGs not represented by analysis manifest: "
            + ", ".join(str(path) for path in stale)
        )
    return assets


def build_public_plot_manifest(analysis_manifest: dict[str, Any], assets: list[dict[str, Any]]) -> dict[str, Any]:
    payload = dict(analysis_manifest)
    payload["figure_layer_role"] = "publication_clean_display_mirror"
    payload["source_analysis_plot_manifest"] = relpath(ANALYSIS_PLOT_MANIFEST)
    payload["source_analysis_plot_manifest_sha256"] = sha256_file(ANALYSIS_PLOT_MANIFEST)
    payload["figures"] = assets
    payload["rendering_semantics"] = (
        str(payload.get("rendering_semantics", ""))
        + "; public figure layer is a byte-preserving mirror of the analysis figures"
    )
    return payload


def build_figure_layer_manifest(
    analysis_package_manifest: dict[str, Any],
    analysis_manifest: dict[str, Any],
    public_plot_manifest: dict[str, Any],
    assets: list[dict[str, Any]],
) -> dict[str, Any]:
    mode_counts = {"mode_a": 0, "mode_b": 0}
    for asset in assets:
        if "/mode_a_fixed_muB_phase_scaled/" in asset["path"]:
            mode_counts["mode_a"] += 1
        elif "/mode_b_fixed_T_sparse_muB/" in asset["path"]:
            mode_counts["mode_b"] += 1
    if sum(mode_counts.values()) != len(assets):
        raise ValueError("could not classify every public figure by mode")
    return {
        "schema": "phase_guided_transport_publication_clean_figure_layer_manifest_v1",
        "figure_layer_role": "publication_clean_display_mirror",
        "status": "derived_author_review_required",
        "manuscript_eligible": False,
        "solver_called": False,
        "source_analysis_package": relpath(ANALYSIS_ROOT),
        "source_plot_manifest": relpath(ANALYSIS_PLOT_MANIFEST),
        "source_plot_manifest_sha256": sha256_file(ANALYSIS_PLOT_MANIFEST),
        "target_plot_manifest": relpath(PUBLIC_PLOT_MANIFEST),
        "target_plot_manifest_sha256": sha256_file(PUBLIC_PLOT_MANIFEST),
        "sync_generator": relpath(Path(__file__).resolve()),
        "sync_generator_sha256": sha256_file(Path(__file__).resolve()),
        "calculation_sha": analysis_package_manifest["calculation_sha"],
        "workflow_head_sha": analysis_package_manifest["workflow_head_sha"],
        "source_case": analysis_package_manifest["source_case"],
        "mode_counts": mode_counts,
        "figure_count": len(assets),
        "content_policy": "byte_preserve",
        "source_package_unchanged": True,
        "raw_results_unchanged": True,
        "production_registry_unchanged": True,
        "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "assets": assets,
    }


def main() -> None:
    analysis_manifest = load_manifest()
    analysis_package_manifest = load_package_manifest()
    assets = sync_figures(analysis_manifest)
    public_plot_manifest = build_public_plot_manifest(analysis_manifest, assets)
    write_json(PUBLIC_PLOT_MANIFEST, public_plot_manifest)
    layer_manifest = build_figure_layer_manifest(
        analysis_package_manifest, analysis_manifest, public_plot_manifest, assets
    )
    write_json(FIGURE_LAYER_MANIFEST, layer_manifest)
    print(
        json.dumps(
            {
                "public_root": relpath(PUBLIC_ROOT),
                "figure_layer_manifest": relpath(FIGURE_LAYER_MANIFEST),
                "figures": len(assets),
                "mode_counts": layer_manifest["mode_counts"],
            },
            ensure_ascii=False,
        )
    )


if __name__ == "__main__":
    main()
