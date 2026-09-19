#!/usr/bin/env python3
"""Record author acceptance of the v4 publication-clean display layer.

This is a metadata-only promotion.  It does not call a solver, rewrite raw
results, or replace the immutable v3 parent snapshot.  The accepted figure
layer remains numerically diagnostic-only because its local adjustments are
display transformations rather than new production calculations.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
TRANSPORT_ROOT = ROOT / "docs" / "analysis" / "relaxtime" / "phase_guided_transport"
V3_ANALYSIS_ROOT = TRANSPORT_ROOT / "phase_guided_transport_publication_clean_v3"
V3_PUBLIC_ROOT = (
    ROOT
    / "data"
    / "outputs"
    / "figures"
    / "relaxtime"
    / "transport"
    / "phase_guided"
    / "publication_clean_v3"
)
V4_ANALYSIS_ROOT = TRANSPORT_ROOT / "phase_guided_transport_publication_clean_v4"
V4_PUBLIC_ROOT = (
    ROOT
    / "data"
    / "outputs"
    / "figures"
    / "relaxtime"
    / "transport"
    / "phase_guided"
    / "publication_clean_v4"
)
V4_LAYER_ROOT = TRANSPORT_ROOT / "phase_guided_transport_publication_clean_figure_layer_v4"
V4_PACKAGE_MANIFEST = V4_ANALYSIS_ROOT / "manifest.json"
V4_PLOT_MANIFEST = V4_ANALYSIS_ROOT / "figures" / "plot_manifest.json"
V4_PUBLIC_PLOT_MANIFEST = V4_PUBLIC_ROOT / "plot_manifest.json"
V4_LAYER_MANIFEST = V4_LAYER_ROOT / "figure_layer_manifest.json"
V3_PACKAGE_MANIFEST = V3_ANALYSIS_ROOT / "manifest.json"
V3_PLOT_MANIFEST = V3_ANALYSIS_ROOT / "figures" / "plot_manifest.json"
V3_PUBLIC_PLOT_MANIFEST = V3_PUBLIC_ROOT / "plot_manifest.json"
FORMALIZATION_RECORD = TRANSPORT_ROOT / "publication_clean_formalization_v1.json"
CURRENT_POINTER = TRANSPORT_ROOT / "publication_clean_current.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relpath(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def require_files() -> None:
    required = (
        V3_PACKAGE_MANIFEST,
        V3_PLOT_MANIFEST,
        V3_PUBLIC_PLOT_MANIFEST,
        V4_PACKAGE_MANIFEST,
        V4_PLOT_MANIFEST,
        V4_PUBLIC_PLOT_MANIFEST,
        V4_LAYER_MANIFEST,
    )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing formalization inputs: " + ", ".join(missing))


def validate_plot_manifest(manifest: dict[str, Any], *, expected_schema: str, root_prefix: str) -> None:
    if manifest.get("schema") != expected_schema:
        raise ValueError(f"unexpected plot manifest schema: {manifest.get('schema')}")
    figures = manifest.get("figures")
    if not isinstance(figures, list) or len(figures) != 72:
        raise ValueError("publication-clean plot manifest must contain 72 figures")
    for item in figures:
        path = ROOT / Path(str(item["path"]))
        if not str(item["path"]).startswith(root_prefix):
            raise ValueError(f"figure escapes expected root: {item['path']}")
        if not path.is_file():
            raise FileNotFoundError(path)
        if int(item["bytes"]) != path.stat().st_size or item["sha256"] != sha256_file(path):
            raise ValueError(f"figure hash mismatch: {path}")


def validate_inputs() -> dict[str, Any]:
    require_files()
    v3_package = read_json(V3_PACKAGE_MANIFEST)
    v3_plot = read_json(V3_PLOT_MANIFEST)
    v3_public = read_json(V3_PUBLIC_PLOT_MANIFEST)
    v4_package = read_json(V4_PACKAGE_MANIFEST)
    v4_plot = read_json(V4_PLOT_MANIFEST)
    v4_public = read_json(V4_PUBLIC_PLOT_MANIFEST)
    v4_layer = read_json(V4_LAYER_MANIFEST)

    if v3_package.get("schema") != "phase_guided_transport_publication_clean_manifest_v3":
        raise ValueError("v3 parent package is not the expected immutable v3 package")
    if v3_plot.get("schema") != "phase_guided_transport_publication_clean_plot_manifest_v3":
        raise ValueError("v3 parent plot manifest is not the expected immutable v3 package")
    validate_plot_manifest(
        v3_plot,
        expected_schema="phase_guided_transport_publication_clean_plot_manifest_v3",
        root_prefix="docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v3/",
    )
    validate_plot_manifest(
        v4_plot,
        expected_schema="phase_guided_transport_publication_clean_plot_manifest_v4",
        root_prefix="docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v4/",
    )
    validate_plot_manifest(
        v4_public,
        expected_schema="phase_guided_transport_publication_clean_plot_manifest_v4",
        root_prefix="data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v4/",
    )
    if v4_package.get("schema") != "phase_guided_transport_publication_clean_manifest_v4":
        raise ValueError("v4 package is not the expected v4 source")
    if v4_layer.get("schema") != "phase_guided_transport_publication_clean_figure_layer_manifest_v4":
        raise ValueError("v4 figure layer is not the expected v4 source")
    if v4_package.get("solver_called") is not False or v4_package.get("production_write") is not False:
        raise ValueError("v4 source is not solver-free/display-only")
    if v4_layer.get("figure_count") != 72 or v4_layer.get("mode_counts") != {"mode_a": 36, "mode_b": 36}:
        raise ValueError("v4 figure layer counts are invalid")
    recorded_parent_sha = v4_package.get("source_parent_manifest_sha256")
    recorded_parent_plot_sha = v4_package.get("source_parent_plot_manifest_sha256")
    if recorded_parent_sha != sha256_file(V3_PACKAGE_MANIFEST):
        raise ValueError("v3 parent manifest hash drifted; refusing promotion")
    if recorded_parent_plot_sha != sha256_file(V3_PLOT_MANIFEST):
        raise ValueError("v3 parent plot manifest hash drifted; refusing promotion")
    if v4_public.get("source_analysis_plot_manifest_sha256") != sha256_file(V4_PLOT_MANIFEST):
        raise ValueError("v4 public plot manifest does not match its analysis source")
    if v4_layer.get("target_plot_manifest_sha256") != sha256_file(V4_PUBLIC_PLOT_MANIFEST):
        raise ValueError("v4 layer manifest does not match its public plot manifest")
    return {
        "v3_package": v3_package,
        "v3_plot": v3_plot,
        "v3_public": v3_public,
        "v4_package": v4_package,
        "v4_plot": v4_plot,
        "v4_public": v4_public,
        "v4_layer": v4_layer,
        "v3_package_sha256": sha256_file(V3_PACKAGE_MANIFEST),
        "v3_plot_sha256": sha256_file(V3_PLOT_MANIFEST),
        "v3_public_plot_sha256": sha256_file(V3_PUBLIC_PLOT_MANIFEST),
        "v4_package_sha256_before": sha256_file(V4_PACKAGE_MANIFEST),
        "v4_plot_sha256_before": sha256_file(V4_PLOT_MANIFEST),
        "v4_public_plot_sha256_before": sha256_file(V4_PUBLIC_PLOT_MANIFEST),
        "v4_layer_sha256_before": sha256_file(V4_LAYER_MANIFEST),
    }


def formal_readme(recorded_at: str) -> str:
    return f"""# Issue #130 RS `publication_clean_v4` formal publication-clean layer

## Acceptance and scope

This is the author-accepted current publication-clean figure layer, promoted
from the reviewed v4 candidate on `{recorded_at}`.  It is a solver-free,
display-only derivative of the approved raw transport case.  The raw CSVs,
production registry, canonical data, and solver outputs are unchanged.

The numerical layer remains `diagnostic_only`: local display interpolation and
one-sided endpoint extrapolation are not new equilibrium or transport
solutions, and they are not a convergence certificate.  The figure layer is
accepted as the formal publication layout.  `manuscript_eligible=false` on
the numerical manifest therefore refers to raw numerical claims, not to the
author acceptance of these display figures.

## Accepted display rules

| scope | rule |
| --- | --- |
| mode-B composite curves | Linear interpolation between xi=0.35 and xi=0.37 at T=200 MeV, muB=0 MeV for zeta, sigma/T, and sigma. |
| mode-A tau_sbar endpoint | Log-linear extrapolation from the left-branch anchors xi=-0.02 and xi=-0.01 to xi=-0.003. |
| phase gate | Crossover-only interpolation; the first-order gap [-0.003, +0.003] remains a hard split. |

The endpoint raw value remains in `tables/publication_clean_points.csv` as
`raw_value`; the display value is retained separately as `clean_value`.  Do
not use the display endpoint for exact branch derivatives or jump amplitudes.

## Provenance

- promoted source: `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v4/`
- v4 source manifest SHA256: recorded in `publication_clean_formalization_v1.json`
- immutable intermediate parent: `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v3/`
- solver called for this layer: `false`
- raw/production data modified: `false`
- formal figure status: `author_accepted_formal_layout`
- numerical status: `diagnostic_only`

The v3 package is retained unchanged as the parent and intermediate evidence
snapshot.  The v4 package is the current publication layer; no lower-version
directory is overwritten or silently reinterpreted.

## Reproduction and audit

The v4 source build remains reproducible with:

    python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v4.py
    python scripts/analysis/relaxtime/sync_phase_guided_publication_clean_figure_layer_v4.py

The explicit author-acceptance record is
`docs/analysis/relaxtime/phase_guided_transport/publication_clean_formalization_v1.json`.
"""


def formal_layer_readme() -> str:
    return """# RS `publication_clean_v4` formal publication figure layer

This directory records the byte-preserving publication figure mirror for the
author-accepted current publication-clean layer.

The v4 analysis package inherits v3, adjusts the three mode-B composite curves
at T=200, muB=0, xi=0.36, and applies one one-sided left-branch display
adjustment to mode-A tau_sbar at the audited first-order left endpoint. The
first-order gap remains split. Raw results, production registries, and solver
outputs are unchanged.

`figure_status=author_accepted_formal_layout` and
`numerical_status=diagnostic_only`. The latter means that display values are
not new solver results or a production convergence certificate.

Acceptance record:

    docs/analysis/relaxtime/phase_guided_transport/publication_clean_formalization_v1.json

Reproduction of the source package:

    python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v4.py
    python scripts/analysis/relaxtime/sync_phase_guided_publication_clean_figure_layer_v4.py
"""


def update_claim_ledger() -> None:
    path = V4_ANALYSIS_ROOT / "tables" / "claim_ledger.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        fields = list(rows[0].keys()) if rows else []
    for row in rows:
        if row.get("claim_id") == "PC-V4-001":
            row["scope_limit"] = "不修改 raw/production 数据，不重新求解；正式化只接受 figure layout，不晋升 raw numerical production。"
        elif row.get("claim_id") == "PC-V4-005":
            row["status"] = "supported_with_scope_limit"
            row["claim_zh"] = "作者已审核并接受 v4 图层作为当前 publication-clean formal figure layer；raw/value provenance 和 display-only 边界保持不变。"
            row["scope_limit"] = "该接受不把显示值变成新的 branch solution、传播子正则化或 production-grade 收敛证明。"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def accepted_fields(recorded_at: str, source: dict[str, Any]) -> dict[str, Any]:
    return {
        "status": "derived_author_accepted_display_only",
        "figure_status": "author_accepted_formal_layout",
        "numerical_status": "diagnostic_only",
        "manuscript_eligible": False,
        "current_publication_layer": True,
        "supersedes": [
            "docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v3",
            "data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v3",
        ],
        "author_acceptance": {
            "status": "accepted",
            "scope": "formal_publication_figure_layer",
            "recorded_at_utc": recorded_at,
            "runtime_default_unchanged": True,
            "raw_results_unchanged": True,
            "production_registry_unchanged": True,
        },
        "formalization_record": relpath(FORMALIZATION_RECORD),
        "promoted_from_artifact": relpath(V4_ANALYSIS_ROOT),
        "promoted_from_manifest_sha256": source["v4_package_sha256_before"],
        "promoted_from_plot_manifest_sha256": source["v4_plot_sha256_before"],
        "promoted_from_public_plot_manifest_sha256": source["v4_public_plot_sha256_before"],
        "promoted_from_figure_layer_manifest_sha256": source["v4_layer_sha256_before"],
    }


def update_manifests(recorded_at: str, source: dict[str, Any]) -> dict[str, str]:
    metadata = accepted_fields(recorded_at, source)
    package = source["v4_package"]
    package.update(metadata)
    package["known_boundaries"] = [
        "publication_clean_v3 remains the immutable intermediate parent snapshot",
        "publication_clean_v2 and publication_clean_v2_residual_smoothed remain unchanged",
        "this is an author-accepted display-only publication layer, not a solver or production rerun",
        "raw_value, v2_clean_value and inherited v3 display values remain available",
        "the mode-A first-order gap remains a hard split and is not filled",
        "the tau_sbar endpoint display value is not a replacement equilibrium branch solution",
        "denominator-chain evidence is mechanism context, not a finite-width regularization",
        "the local high-rate convergence gate was not run",
    ]
    V4_ANALYSIS_ROOT.joinpath("README.md").write_text(
        formal_readme(recorded_at), encoding="utf-8", newline="\n"
    )
    V4_LAYER_ROOT.joinpath("README.md").write_text(
        formal_layer_readme(), encoding="utf-8", newline="\n"
    )
    update_claim_ledger()
    write_json(V4_PACKAGE_MANIFEST, package)

    plot = source["v4_plot"]
    plot.update(metadata)
    plot["formalization_record"] = relpath(FORMALIZATION_RECORD)
    write_json(V4_PLOT_MANIFEST, plot)
    analysis_plot_sha = sha256_file(V4_PLOT_MANIFEST)

    public = source["v4_public"]
    public.update(metadata)
    public["source_analysis_plot_manifest"] = relpath(V4_PLOT_MANIFEST)
    public["source_analysis_plot_manifest_sha256"] = analysis_plot_sha
    public["formalization_record"] = relpath(FORMALIZATION_RECORD)
    write_json(V4_PUBLIC_PLOT_MANIFEST, public)
    public_plot_sha = sha256_file(V4_PUBLIC_PLOT_MANIFEST)

    layer = source["v4_layer"]
    layer.update(metadata)
    layer["source_analysis_package_manifest_sha256"] = sha256_file(V4_PACKAGE_MANIFEST)
    layer["source_plot_manifest_sha256"] = analysis_plot_sha
    layer["target_plot_manifest_sha256"] = public_plot_sha
    layer["formalization_record"] = relpath(FORMALIZATION_RECORD)
    layer["generated_at_utc"] = recorded_at
    write_json(V4_LAYER_MANIFEST, layer)

    # The package manifest records all package files except itself.  Refresh
    # those records after the README, claim ledger and plot manifest changed.
    package = read_json(V4_PACKAGE_MANIFEST)
    records = []
    for item in package.get("outputs", []):
        path = ROOT / Path(str(item["path"]))
        if not path.is_file():
            raise FileNotFoundError(path)
        records.append({"path": relpath(path), "sha256": sha256_file(path), "bytes": path.stat().st_size})
    package["outputs"] = records
    write_json(V4_PACKAGE_MANIFEST, package)

    # The layer manifest points at the final package hash, so refresh it once
    # more after the package output records were updated.
    layer = read_json(V4_LAYER_MANIFEST)
    layer["source_analysis_package_manifest_sha256"] = sha256_file(V4_PACKAGE_MANIFEST)
    write_json(V4_LAYER_MANIFEST, layer)
    return {
        "package_manifest_sha256": sha256_file(V4_PACKAGE_MANIFEST),
        "plot_manifest_sha256": sha256_file(V4_PLOT_MANIFEST),
        "public_plot_manifest_sha256": sha256_file(V4_PUBLIC_PLOT_MANIFEST),
        "figure_layer_manifest_sha256": sha256_file(V4_LAYER_MANIFEST),
        "claim_ledger_sha256": sha256_file(V4_ANALYSIS_ROOT / "tables" / "claim_ledger.csv"),
    }


def write_records(recorded_at: str, source: dict[str, Any], final_hashes: dict[str, str]) -> None:
    record = {
        "schema": "phase_guided_transport_publication_clean_formalization_v1",
        "status": "author_accepted_formal_publication_figure_layer",
        "recorded_at_utc": recorded_at,
        "current_publication_layer": relpath(V4_PUBLIC_ROOT),
        "current_analysis_package": relpath(V4_ANALYSIS_ROOT),
        "current_figure_layer_manifest": relpath(V4_LAYER_MANIFEST),
        "source_v4_manifest_sha256": source["v4_package_sha256_before"],
        "source_v4_plot_manifest_sha256": source["v4_plot_sha256_before"],
        "source_v4_public_plot_manifest_sha256": source["v4_public_plot_sha256_before"],
        "source_v4_figure_layer_manifest_sha256": source["v4_layer_sha256_before"],
        "immutable_v3_parent_manifest_sha256": source["v3_package_sha256"],
        "immutable_v3_parent_plot_manifest_sha256": source["v3_plot_sha256"],
        "immutable_v3_parent_public_plot_manifest_sha256": source["v3_public_plot_sha256"],
        "final_hashes": final_hashes,
        "figure_count": 72,
        "mode_counts": {"mode_a": 36, "mode_b": 36},
        "solver_called": False,
        "production_write": False,
        "canonical_data_modified": False,
        "raw_results_unchanged": True,
        "production_registry_unchanged": True,
        "numerical_status": "diagnostic_only",
        "figure_status": "author_accepted_formal_layout",
        "superseded_intermediate": {
            "analysis": relpath(V3_ANALYSIS_ROOT),
            "public": relpath(V3_PUBLIC_ROOT),
            "status": "retained_unchanged_as_parent_evidence",
        },
        "known_boundaries": [
            "display values are not new solver results",
            "the tau_sbar endpoint is not a replacement branch solution",
            "mechanism evidence does not constitute a high-rate convergence gate",
        ],
    }
    write_json(FORMALIZATION_RECORD, record)
    write_json(
        CURRENT_POINTER,
        {
            "schema": "phase_guided_transport_publication_clean_current_v1",
            "status": "author_accepted_formal_publication_figure_layer",
            "current_analysis_package": relpath(V4_ANALYSIS_ROOT),
            "current_public_figure_root": relpath(V4_PUBLIC_ROOT),
            "current_figure_layer": relpath(V4_LAYER_ROOT),
            "formalization_record": relpath(FORMALIZATION_RECORD),
            "manifest_sha256": final_hashes["package_manifest_sha256"],
            "plot_manifest_sha256": final_hashes["plot_manifest_sha256"],
            "public_plot_manifest_sha256": final_hashes["public_plot_manifest_sha256"],
            "superseded_intermediate": relpath(V3_ANALYSIS_ROOT),
            "manuscript_eligible": False,
            "numerical_status": "diagnostic_only",
            "figure_status": "author_accepted_formal_layout",
        },
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true", help="write the acceptance metadata")
    parser.add_argument("--recorded-at", help="UTC ISO-8601 acceptance timestamp")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    source = validate_inputs()
    recorded_at = args.recorded_at or dt.datetime.now(dt.timezone.utc).isoformat()
    if not args.apply:
        print(json.dumps({"status": "dry_run", "v3_parent_unchanged_hash": source["v3_package_sha256"]}, ensure_ascii=False))
        return 0
    if FORMALIZATION_RECORD.exists():
        raise FileExistsError(
            f"formalization already recorded at {FORMALIZATION_RECORD}; refusing a second promotion"
        )
    final_hashes = update_manifests(recorded_at, source)
    write_records(recorded_at, source, final_hashes)
    print(json.dumps({"status": "accepted", "recorded_at_utc": recorded_at, **final_hashes}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
