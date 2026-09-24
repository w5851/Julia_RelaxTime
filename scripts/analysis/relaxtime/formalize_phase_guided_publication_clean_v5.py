#!/usr/bin/env python3
"""Record formal manuscript adoption of the v5 display layer.

This is a metadata-only promotion of the label-only v5 derivative.  It does
not call a solver, modify raw numerical data, change the production registry,
or claim that the local high-rate convergence gate passed.
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

V4_ROOT = TRANSPORT_ROOT / "phase_guided_transport_publication_clean_v4"
V5_ROOT = TRANSPORT_ROOT / "phase_guided_transport_publication_clean_v5"
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
V5_PUBLIC_ROOT = (
    ROOT
    / "data"
    / "outputs"
    / "figures"
    / "relaxtime"
    / "transport"
    / "phase_guided"
    / "publication_clean_v5"
)
V4_LAYER_ROOT = TRANSPORT_ROOT / "phase_guided_transport_publication_clean_figure_layer_v4"
V5_LAYER_ROOT = TRANSPORT_ROOT / "phase_guided_transport_publication_clean_figure_layer_v5"

V4_PACKAGE_MANIFEST = V4_ROOT / "manifest.json"
V4_PLOT_MANIFEST = V4_ROOT / "figures" / "plot_manifest.json"
V4_PUBLIC_PLOT_MANIFEST = V4_PUBLIC_ROOT / "plot_manifest.json"
V4_LAYER_MANIFEST = V4_LAYER_ROOT / "figure_layer_manifest.json"
V4_ELIGIBILITY_RECORD = TRANSPORT_ROOT / "publication_clean_v4_manuscript_eligibility_v1.json"

V5_PACKAGE_MANIFEST = V5_ROOT / "manifest.json"
V5_PLOT_MANIFEST = V5_ROOT / "figures" / "plot_manifest.json"
V5_PUBLIC_PLOT_MANIFEST = V5_PUBLIC_ROOT / "plot_manifest.json"
V5_LAYER_MANIFEST = V5_LAYER_ROOT / "figure_layer_manifest.json"
V5_FORMALIZATION_RECORD = TRANSPORT_ROOT / "publication_clean_v5_formalization_v1.json"
V5_ELIGIBILITY_RECORD = TRANSPORT_ROOT / "publication_clean_v5_manuscript_eligibility_v1.json"
CURRENT_POINTER = TRANSPORT_ROOT / "publication_clean_current.json"

CRLF_JSON_PATHS = {
    path.resolve()
    for path in (
        V5_PACKAGE_MANIFEST,
        V5_PLOT_MANIFEST,
        V5_PUBLIC_PLOT_MANIFEST,
        V5_LAYER_MANIFEST,
    )
}

TAU_FIELDS = ("tau_u", "tau_d", "tau_s", "tau_ubar", "tau_dbar", "tau_sbar")
EXPECTED_LABELS = {
    ("relaxation_time_y_axis", "tau_u"): (r"$\tau_u$", r"$\tau_u\;[\mathrm{fm}]$"),
    ("relaxation_time_y_axis", "tau_d"): (r"$\tau_d$", r"$\tau_d\;[\mathrm{fm}]$"),
    ("relaxation_time_y_axis", "tau_s"): (r"$\tau_s$", r"$\tau_s\;[\mathrm{fm}]$"),
    ("relaxation_time_y_axis", "tau_ubar"): (r"$\tau_{\bar u}$", r"$\tau_{\bar u}\;[\mathrm{fm}]$"),
    ("relaxation_time_y_axis", "tau_dbar"): (r"$\tau_{\bar d}$", r"$\tau_{\bar d}\;[\mathrm{fm}]$"),
    ("relaxation_time_y_axis", "tau_sbar"): (r"$\tau_{\bar s}$", r"$\tau_{\bar s}\;[\mathrm{fm}]$"),
    ("first_order_endpoint_legend", "quark"): (
        "chiral-restored (quark) endpoint",
        "chirally restored branch endpoint",
    ),
    ("first_order_endpoint_legend", "hadron"): (
        "chiral-broken (hadron) endpoint",
        "chirally broken branch endpoint",
    ),
}


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
    newline = "\r\n" if path.resolve() in CRLF_JSON_PATHS else "\n"
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline=newline,
    )


def validate_plot_manifest(
    manifest: dict[str, Any], *, expected_schema: str, root_prefix: str
) -> None:
    if manifest.get("schema") != expected_schema:
        raise ValueError(f"unexpected plot manifest schema: {manifest.get('schema')}")
    figures = manifest.get("figures")
    if not isinstance(figures, list) or len(figures) != 72:
        raise ValueError("publication-clean plot manifest must contain 72 figures")
    for item in figures:
        path_text = str(item["path"])
        path = ROOT / Path(path_text)
        if not path_text.startswith(root_prefix):
            raise ValueError(f"figure escapes expected root: {path_text}")
        if not path.is_file():
            raise FileNotFoundError(path)
        if int(item["bytes"]) != path.stat().st_size or item["sha256"] != sha256_file(path):
            raise ValueError(f"figure hash mismatch: {path}")


def validate_label_map() -> None:
    path = V5_ROOT / "tables" / "v5_display_label_map.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    actual = {
        (row.get("scope", ""), row.get("key", "")): (
            row.get("old_label", ""),
            row.get("new_label", ""),
        )
        for row in rows
    }
    if actual != EXPECTED_LABELS:
        raise ValueError("v5 display label map does not match the reviewed label changes")
    if any(row.get("canonical_data_modified") != "False" for row in rows):
        raise ValueError("v5 label map claims a canonical-data modification")


def validate_inherited_tables() -> None:
    v4_tables = V4_ROOT / "tables"
    v5_tables = V5_ROOT / "tables"
    for path in v4_tables.glob("*.csv"):
        # The claim ledger is extended when v5 adoption is recorded; all
        # numerical/audit tables remain byte-identical to the v4 parent.
        if path.name == "claim_ledger.csv":
            continue
        target = v5_tables / path.name
        if not target.is_file() or sha256_file(path) != sha256_file(target):
            raise ValueError(f"v5 changed inherited table bytes: {path.name}")


def validate_inputs() -> dict[str, Any]:
    required = (
        V4_PACKAGE_MANIFEST,
        V4_PLOT_MANIFEST,
        V4_PUBLIC_PLOT_MANIFEST,
        V4_LAYER_MANIFEST,
        V4_ELIGIBILITY_RECORD,
        V5_PACKAGE_MANIFEST,
        V5_PLOT_MANIFEST,
        V5_PUBLIC_PLOT_MANIFEST,
        V5_LAYER_MANIFEST,
        CURRENT_POINTER,
    )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing v5 formalization inputs: " + ", ".join(missing))

    v4_package = read_json(V4_PACKAGE_MANIFEST)
    v4_plot = read_json(V4_PLOT_MANIFEST)
    v4_public = read_json(V4_PUBLIC_PLOT_MANIFEST)
    v4_layer = read_json(V4_LAYER_MANIFEST)
    v4_eligibility = read_json(V4_ELIGIBILITY_RECORD)
    v5_package = read_json(V5_PACKAGE_MANIFEST)
    v5_plot = read_json(V5_PLOT_MANIFEST)
    v5_public = read_json(V5_PUBLIC_PLOT_MANIFEST)
    v5_layer = read_json(V5_LAYER_MANIFEST)
    pointer = read_json(CURRENT_POINTER)

    if v4_package.get("schema") != "phase_guided_transport_publication_clean_manifest_v4":
        raise ValueError("v4 parent package is not the expected accepted package")
    if v4_package.get("manuscript_eligible") is not True:
        raise ValueError("v4 parent package is not manuscript eligible")
    if v4_eligibility.get("manuscript_eligible") is not True:
        raise ValueError("v4 parent eligibility record is not accepted")
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
    if v4_layer.get("schema") != "phase_guided_transport_publication_clean_figure_layer_manifest_v4":
        raise ValueError("v4 parent figure layer is not the expected accepted layer")
    if v4_layer.get("manuscript_eligible") is not True:
        raise ValueError("v4 parent figure layer is not manuscript eligible")

    if v5_package.get("schema") != "phase_guided_transport_publication_clean_manifest_v5":
        raise ValueError("v5 package is not the expected label-only derivative")
    if v5_plot.get("schema") != "phase_guided_transport_publication_clean_plot_manifest_v5":
        raise ValueError("v5 plot manifest is not the expected label-only derivative")
    if v5_layer.get("schema") != "phase_guided_transport_publication_clean_figure_layer_manifest_v5":
        raise ValueError("v5 figure layer is not the expected label-only derivative")
    validate_plot_manifest(
        v5_plot,
        expected_schema="phase_guided_transport_publication_clean_plot_manifest_v5",
        root_prefix="docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v5/",
    )
    validate_plot_manifest(
        v5_public,
        expected_schema="phase_guided_transport_publication_clean_plot_manifest_v5",
        root_prefix="data/outputs/figures/relaxtime/transport/phase_guided/publication_clean_v5/",
    )
    if v5_package.get("manuscript_eligible") is not False:
        raise ValueError("v5 package is not at the expected pre-formalization state")
    if v5_plot.get("manuscript_eligible") is not False:
        raise ValueError("v5 plot manifest is not at the expected pre-formalization state")
    if v5_public.get("manuscript_eligible") is not False:
        raise ValueError("v5 public plot manifest is not at the expected pre-formalization state")
    if v5_layer.get("manuscript_eligible") is not False:
        raise ValueError("v5 figure layer is not at the expected pre-formalization state")
    if v5_package.get("solver_called") is not False or v5_package.get("canonical_data_modified") is not False:
        raise ValueError("v5 is not solver-free and canonical-data preserving")
    if v5_package.get("source_parent_manifest_sha256") != sha256_file(V4_PACKAGE_MANIFEST):
        raise ValueError("v5 parent package hash drifted")
    if v5_package.get("source_parent_plot_manifest_sha256") != sha256_file(V4_PLOT_MANIFEST):
        raise ValueError("v5 parent plot hash drifted")
    if v5_package.get("source_parent_eligibility_record_sha256") != sha256_file(V4_ELIGIBILITY_RECORD):
        raise ValueError("v5 parent eligibility hash drifted")
    if v5_public.get("source_analysis_plot_manifest_sha256") != sha256_file(V5_PLOT_MANIFEST):
        raise ValueError("v5 public plot manifest does not match its analysis source")
    if v5_layer.get("target_plot_manifest_sha256") != sha256_file(V5_PUBLIC_PLOT_MANIFEST):
        raise ValueError("v5 layer manifest does not match its public plot manifest")
    validate_label_map()
    validate_inherited_tables()

    return {
        "v4_package": v4_package,
        "v4_plot": v4_plot,
        "v4_public": v4_public,
        "v4_layer": v4_layer,
        "v4_eligibility": v4_eligibility,
        "v5_package": v5_package,
        "v5_plot": v5_plot,
        "v5_public": v5_public,
        "v5_layer": v5_layer,
        "pointer": pointer,
        "v4_package_sha256": sha256_file(V4_PACKAGE_MANIFEST),
        "v4_plot_sha256": sha256_file(V4_PLOT_MANIFEST),
        "v4_public_plot_sha256": sha256_file(V4_PUBLIC_PLOT_MANIFEST),
        "v4_layer_sha256": sha256_file(V4_LAYER_MANIFEST),
        "v4_eligibility_sha256": sha256_file(V4_ELIGIBILITY_RECORD),
        "v5_package_sha256_before": sha256_file(V5_PACKAGE_MANIFEST),
        "v5_plot_sha256_before": sha256_file(V5_PLOT_MANIFEST),
        "v5_public_plot_sha256_before": sha256_file(V5_PUBLIC_PLOT_MANIFEST),
        "v5_layer_sha256_before": sha256_file(V5_LAYER_MANIFEST),
    }


def accepted_fields(recorded_at: str) -> dict[str, Any]:
    return {
        "status": "derived_author_accepted_display_only",
        "figure_status": "author_accepted_formal_layout",
        "numerical_status": "author_accepted_display_only",
        "manuscript_eligible": True,
        "manuscript_eligibility_scope": "publication_clean_v5_display_layer_only",
        "raw_numerical_status": "diagnostic_only",
        "raw_manuscript_eligible": False,
        "convergence_gate_status": "local_high_rate_gate_not_run",
        "convergence_claim": False,
        "current_publication_layer": True,
        "eligibility_decision_basis": "explicit_author_confirmation_that_formal_manuscript_uses_v5",
        "author_acceptance": {
            "status": "accepted",
            "scope": "formal_publication_figure_layer",
            "recorded_at_utc": recorded_at,
            "source": "formal_manuscript_uses_publication_clean_v5",
            "runtime_default_unchanged": True,
            "raw_results_unchanged": True,
            "production_registry_unchanged": True,
        },
        "formalization_record": relpath(V5_FORMALIZATION_RECORD),
        "manuscript_eligibility_record": relpath(V5_ELIGIBILITY_RECORD),
        "parent_accepted_display_layer": relpath(V4_ROOT),
        "parent_eligibility_record": relpath(V4_ELIGIBILITY_RECORD),
    }


def update_readmes(recorded_at: str) -> None:
    (V5_ROOT / "README.md").write_text(
        f"""# Issue #130 RS `publication_clean_v5` formal publication-clean layer

## Acceptance and scope

This is the author-accepted current publication-clean figure layer, adopted
for the formal manuscript on `{recorded_at}`.  It is a solver-free, label-only
display derivative of the accepted `publication_clean_v4` layer.  The raw CSVs,
production registry, canonical data, display values, and solver outputs are
unchanged.

The accepted rendering changes are:

1. The six relaxation-time y-axis labels include the unit `[fm]`.
2. The first-order endpoint legend entries read `chirally restored branch
   endpoint` and `chirally broken branch endpoint`.  The latter describes the
   displayed quark transport branch and does not claim that hadronic transport
   was calculated.

`manuscript_eligible=true` applies only to this disclosed v5 display layer.
The source raw numerical results remain `diagnostic_only` and
`raw_manuscript_eligible=false`.  The local high-rate convergence gate was not
run and is not represented as passed.

## Provenance

- parent display layer: `docs/analysis/relaxtime/phase_guided_transport/phase_guided_transport_publication_clean_v4`
- parent manuscript eligibility record: `docs/analysis/relaxtime/phase_guided_transport/publication_clean_v4_manuscript_eligibility_v1.json`
- publication figures: 72 (36 mode A and 36 mode B)
- raw/production data modified: false
- solver called for this derivative: false
- canonical numerical data modified: false
- label audit: `tables/v5_display_label_map.csv`
- raw/display provenance inherited from v4: `tables/publication_clean_points.csv` and `tables/v4_display_adjustment_map.csv`

The v4 parent remains retained unchanged as provenance.  No first-order gap is
filled, no branch solution is replaced, and no display label change is a
physical-model or convergence claim.

## Records and reproduction

The formalization record is
`../publication_clean_v5_formalization_v1.json`; the manuscript display
eligibility record is `../publication_clean_v5_manuscript_eligibility_v1.json`.

The candidate layer was generated with:

    python scripts/analysis/relaxtime/build_phase_guided_publication_clean_v5.py
    python scripts/analysis/relaxtime/sync_phase_guided_publication_clean_figure_layer_v5.py
""",
        encoding="utf-8",
        newline="\n",
    )
    (V5_LAYER_ROOT / "README.md").write_text(
        f"""# RS `publication_clean_v5` formal publication figure layer

This directory records the byte-preserving public mirror for the
author-accepted current publication-clean layer, adopted on `{recorded_at}`.
It is the label-only derivative of `publication_clean_v4`: relaxation-time
axes carry `[fm]`, and the endpoint legend uses chirally restored/broken branch
wording.  Numerical values, raw data, audit tables, phase gaps, and solver
outputs are unchanged.

`figure_status=author_accepted_formal_layout` and
`manuscript_eligible=true` apply to the disclosed display layer only.
`numerical_status=author_accepted_display_only`, while the source raw numerical
status remains `diagnostic_only` and `raw_manuscript_eligible=false`.  No local
high-rate convergence gate was run.

Acceptance records:

    docs/analysis/relaxtime/phase_guided_transport/publication_clean_v5_formalization_v1.json
    docs/analysis/relaxtime/phase_guided_transport/publication_clean_v5_manuscript_eligibility_v1.json
""",
        encoding="utf-8",
        newline="\n",
    )


def update_claim_ledger() -> None:
    path = V5_ROOT / "tables" / "claim_ledger.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    existing = {row.get("claim_id") for row in rows}
    additions = [
        {
            "claim_id": "PC-V5-001",
            "status": "supported_with_scope_limit",
            "claim_zh": "publication_clean_v5 仅更新六个弛豫时间纵轴单位和一阶端点图例措辞；数值、相别语义和 raw/display provenance 不变。",
            "evidence": "manifest.json; figures/plot_manifest.json; tables/v5_display_label_map.csv",
            "scope_limit": "这是显示标签更新，不是 solver、传播子正则化或 production numerical update。",
        },
        {
            "claim_id": "PC-V5-002",
            "status": "author_accepted_with_scope_limit",
            "claim_zh": "作者已将 publication_clean_v5 的 72 张图用于正式文稿；v5 只在 v4 已接受显示层上改变标签。",
            "evidence": "publication_clean_v5_manuscript_eligibility_v1.json; README.md; tables/v5_display_label_map.csv",
            "scope_limit": "资格仅覆盖 v5 display layer；source raw numerical status 仍为 diagnostic_only，local high-rate gate 未运行。",
        },
    ]
    for row in additions:
        if row["claim_id"] not in existing:
            rows.append(row)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def refresh_package_outputs(package: dict[str, Any]) -> None:
    records = []
    for item in package.get("outputs", []):
        path = ROOT / Path(str(item["path"]))
        if not path.is_file():
            raise FileNotFoundError(path)
        records.append({"path": relpath(path), "sha256": sha256_file(path), "bytes": path.stat().st_size})
    package["outputs"] = records


def update_manifests(source: dict[str, Any], recorded_at: str) -> dict[str, str]:
    metadata = accepted_fields(recorded_at)
    package = source["v5_package"]
    plot = source["v5_plot"]
    public = source["v5_public"]
    layer = source["v5_layer"]
    for manifest in (package, plot, public, layer):
        manifest.update(metadata)

    inherited_boundaries = [
        item
        for item in package.get("known_boundaries", [])
        if "manuscript eligibility applies to the author-accepted v4 display layer only" not in item
    ]
    package["known_boundaries"] = list(
        dict.fromkeys(
            inherited_boundaries
            + [
                "publication_clean_v4 remains the immutable accepted parent display snapshot",
                "this is an author-accepted display-only publication layer, not a solver or production rerun",
                "v5 changes only tau-axis units and first-order endpoint legend wording",
                "raw numerical status remains diagnostic_only and raw_manuscript_eligible=false",
                "the local high-rate convergence gate was not run and is not claimed as passed",
                "the first-order gap and all inherited display adjustment boundaries remain unchanged",
            ]
        )
    )
    update_readmes(recorded_at)
    update_claim_ledger()

    write_json(V5_PLOT_MANIFEST, plot)
    analysis_plot_sha = sha256_file(V5_PLOT_MANIFEST)
    public["source_analysis_plot_manifest"] = relpath(V5_PLOT_MANIFEST)
    public["source_analysis_plot_manifest_sha256"] = analysis_plot_sha
    write_json(V5_PUBLIC_PLOT_MANIFEST, public)
    public_plot_sha = sha256_file(V5_PUBLIC_PLOT_MANIFEST)

    refresh_package_outputs(package)
    write_json(V5_PACKAGE_MANIFEST, package)
    package_sha = sha256_file(V5_PACKAGE_MANIFEST)

    layer["source_analysis_package_manifest_sha256"] = package_sha
    layer["source_plot_manifest_sha256"] = analysis_plot_sha
    layer["target_plot_manifest_sha256"] = public_plot_sha
    layer["formalization_record"] = relpath(V5_FORMALIZATION_RECORD)
    layer["manuscript_eligibility_record"] = relpath(V5_ELIGIBILITY_RECORD)
    write_json(V5_LAYER_MANIFEST, layer)
    layer_sha = sha256_file(V5_LAYER_MANIFEST)

    final_hashes = {
        "package_manifest_sha256": package_sha,
        "plot_manifest_sha256": analysis_plot_sha,
        "public_plot_manifest_sha256": public_plot_sha,
        "figure_layer_manifest_sha256": layer_sha,
        "claim_ledger_sha256": sha256_file(V5_ROOT / "tables" / "claim_ledger.csv"),
    }
    formalization = {
        "schema": "phase_guided_transport_publication_clean_v5_formalization_v1",
        "status": "author_accepted_formal_publication_figure_layer",
        "recorded_at_utc": recorded_at,
        "decision_basis": "explicit_author_confirmation_that_formal_manuscript_uses_v5",
        "current_publication_layer": relpath(V5_PUBLIC_ROOT),
        "current_analysis_package": relpath(V5_ROOT),
        "current_figure_layer_manifest": relpath(V5_LAYER_MANIFEST),
        "parent_v4_eligibility_record": relpath(V4_ELIGIBILITY_RECORD),
        "parent_v4_eligibility_record_sha256": source["v4_eligibility_sha256"],
        "source_v5_hashes_before": {
            "package_manifest_sha256": source["v5_package_sha256_before"],
            "plot_manifest_sha256": source["v5_plot_sha256_before"],
            "public_plot_manifest_sha256": source["v5_public_plot_sha256_before"],
            "figure_layer_manifest_sha256": source["v5_layer_sha256_before"],
        },
        "parent_v4_hashes": {
            "package_manifest_sha256": source["v4_package_sha256"],
            "plot_manifest_sha256": source["v4_plot_sha256"],
            "public_plot_manifest_sha256": source["v4_public_plot_sha256"],
            "figure_layer_manifest_sha256": source["v4_layer_sha256"],
        },
        "final_hashes": final_hashes,
        "figure_count": 72,
        "mode_counts": {"mode_a": 36, "mode_b": 36},
        "solver_called": False,
        "production_write": False,
        "canonical_data_modified": False,
        "raw_results_unchanged": True,
        "production_registry_unchanged": True,
        "numerical_status": "author_accepted_display_only",
        "manuscript_eligible": True,
        "manuscript_eligibility_scope": "publication_clean_v5_display_layer_only",
        "raw_numerical_status": "diagnostic_only",
        "raw_manuscript_eligible": False,
        "convergence_gate_status": "local_high_rate_gate_not_run",
        "convergence_claim": False,
        "display_label_update": {
            "kind": "display_label_only",
            "tau_y_axis_unit": "fm",
            "tau_fields": list(TAU_FIELDS),
            "endpoint_legend_labels": {
                "quark": "chirally restored branch endpoint",
                "hadron": "chirally broken branch endpoint",
            },
            "canonical_data_modified": False,
            "solver_called": False,
        },
        "known_boundaries": [
            "v4 remains retained unchanged as the accepted parent display layer",
            "raw numerical results remain diagnostic_only",
            "no high-rate convergence claim is made",
            "v5 label changes do not alter values, phase gaps, or branch semantics",
        ],
    }
    write_json(V5_FORMALIZATION_RECORD, formalization)
    formalization_sha = sha256_file(V5_FORMALIZATION_RECORD)

    eligibility = {
        "schema": "phase_guided_transport_publication_clean_v5_manuscript_eligibility_v1",
        "status": "author_accepted_display_only_manuscript_eligible",
        "recorded_at_utc": recorded_at,
        "decision_basis": "explicit_author_confirmation_that_formal_manuscript_uses_v5",
        "current_analysis_package": relpath(V5_ROOT),
        "current_public_figure_root": relpath(V5_PUBLIC_ROOT),
        "current_figure_layer_manifest": relpath(V5_LAYER_MANIFEST),
        "formalization_record": relpath(V5_FORMALIZATION_RECORD),
        "formalization_record_sha256": formalization_sha,
        "parent_v4_eligibility_record": relpath(V4_ELIGIBILITY_RECORD),
        "parent_v4_eligibility_record_sha256": source["v4_eligibility_sha256"],
        "final_hashes": final_hashes,
        "figure_count": 72,
        "mode_counts": {"mode_a": 36, "mode_b": 36},
        "numerical_status": "author_accepted_display_only",
        "manuscript_eligible": True,
        "manuscript_eligibility_scope": "publication_clean_v5_display_layer_only",
        "raw_numerical_status": "diagnostic_only",
        "raw_manuscript_eligible": False,
        "source_production_status": "diagnostic_only",
        "source_solver_called": True,
        "display_derivative_solver_called": False,
        "production_write": False,
        "raw_results_modified": False,
        "production_registry_modified": False,
        "convergence_gate_status": "local_high_rate_gate_not_run",
        "convergence_claim": False,
        "display_label_map": relpath(V5_ROOT / "tables" / "v5_display_label_map.csv"),
        "raw_display_provenance": {
            "point_table": relpath(V5_ROOT / "tables" / "publication_clean_points.csv"),
            "adjustment_map": relpath(V5_ROOT / "tables" / "v4_display_adjustment_map.csv"),
            "inherited_parent": relpath(V4_ROOT),
        },
        "accepted_scope": [
            "six tau y-axis units are labeled in fm",
            "first-order endpoint legend wording uses chirally restored/broken branch terminology",
        ],
        "boundaries": [
            "eligibility covers the v5 display derivative only",
            "raw numerical results remain diagnostic_only and are not promoted to production_grade",
            "the local high-rate convergence gate was not run and is not claimed as passed",
            "v5 does not change numerical values, phase gaps, branch solutions, or transport kernels",
        ],
    }
    write_json(V5_ELIGIBILITY_RECORD, eligibility)
    eligibility_sha = sha256_file(V5_ELIGIBILITY_RECORD)

    pointer = {
        "schema": "phase_guided_transport_publication_clean_current_v1",
        "status": "author_accepted_formal_publication_figure_layer",
        "current_analysis_package": relpath(V5_ROOT),
        "current_public_figure_root": relpath(V5_PUBLIC_ROOT),
        "current_figure_layer": relpath(V5_LAYER_ROOT),
        "formalization_record": relpath(V5_FORMALIZATION_RECORD),
        "manifest_sha256": final_hashes["package_manifest_sha256"],
        "plot_manifest_sha256": final_hashes["plot_manifest_sha256"],
        "public_plot_manifest_sha256": final_hashes["public_plot_manifest_sha256"],
        "figure_layer_manifest_sha256": final_hashes["figure_layer_manifest_sha256"],
        "superseded_intermediate": relpath(V4_ROOT),
        "superseded_publication_layer": relpath(V4_PUBLIC_ROOT),
        "manuscript_eligible": True,
        "numerical_status": "author_accepted_display_only",
        "figure_status": "author_accepted_formal_layout",
        "manuscript_eligibility_scope": "publication_clean_v5_display_layer_only",
        "raw_numerical_status": "diagnostic_only",
        "raw_manuscript_eligible": False,
        "convergence_gate_status": "local_high_rate_gate_not_run",
        "manuscript_eligibility_record": relpath(V5_ELIGIBILITY_RECORD),
        "manuscript_eligibility_record_sha256": eligibility_sha,
        "eligibility_decision_basis": "explicit_author_confirmation_that_formal_manuscript_uses_v5",
        "convergence_claim": False,
    }
    write_json(CURRENT_POINTER, pointer)
    return {**final_hashes, "formalization_record_sha256": formalization_sha, "eligibility_record_sha256": eligibility_sha}


def refresh_accepted_state() -> dict[str, str]:
    """Refresh cross-manifest hashes after a representation-only repair.

    This is intentionally limited to an already accepted v5 state.  It exists
    so line-ending or equivalent hash-bound representation repairs do not
    require re-running the figure builder or editing provenance hashes by hand.
    """
    package = read_json(V5_PACKAGE_MANIFEST)
    plot = read_json(V5_PLOT_MANIFEST)
    public = read_json(V5_PUBLIC_PLOT_MANIFEST)
    layer = read_json(V5_LAYER_MANIFEST)
    formalization = read_json(V5_FORMALIZATION_RECORD)
    eligibility = read_json(V5_ELIGIBILITY_RECORD)
    pointer = read_json(CURRENT_POINTER)
    if any(manifest.get("manuscript_eligible") is not True for manifest in (package, plot, public, layer)):
        raise ValueError("v5 is not already in the accepted state")
    if formalization.get("schema") != "phase_guided_transport_publication_clean_v5_formalization_v1":
        raise ValueError("v5 formalization record is missing")
    validate_label_map()

    write_json(V5_PLOT_MANIFEST, plot)
    analysis_plot_sha = sha256_file(V5_PLOT_MANIFEST)
    public["source_analysis_plot_manifest_sha256"] = analysis_plot_sha
    write_json(V5_PUBLIC_PLOT_MANIFEST, public)
    public_plot_sha = sha256_file(V5_PUBLIC_PLOT_MANIFEST)

    refresh_package_outputs(package)
    write_json(V5_PACKAGE_MANIFEST, package)
    package_sha = sha256_file(V5_PACKAGE_MANIFEST)

    layer["source_analysis_package_manifest_sha256"] = package_sha
    layer["source_plot_manifest_sha256"] = analysis_plot_sha
    layer["target_plot_manifest_sha256"] = public_plot_sha
    write_json(V5_LAYER_MANIFEST, layer)
    layer_sha = sha256_file(V5_LAYER_MANIFEST)

    final_hashes = {
        "package_manifest_sha256": package_sha,
        "plot_manifest_sha256": analysis_plot_sha,
        "public_plot_manifest_sha256": public_plot_sha,
        "figure_layer_manifest_sha256": layer_sha,
        "claim_ledger_sha256": sha256_file(V5_ROOT / "tables" / "claim_ledger.csv"),
    }
    formalization["final_hashes"] = final_hashes
    write_json(V5_FORMALIZATION_RECORD, formalization)
    formalization_sha = sha256_file(V5_FORMALIZATION_RECORD)
    eligibility["final_hashes"] = final_hashes
    eligibility["formalization_record_sha256"] = formalization_sha
    write_json(V5_ELIGIBILITY_RECORD, eligibility)
    eligibility_sha = sha256_file(V5_ELIGIBILITY_RECORD)
    pointer["manifest_sha256"] = package_sha
    pointer["plot_manifest_sha256"] = analysis_plot_sha
    pointer["public_plot_manifest_sha256"] = public_plot_sha
    pointer["figure_layer_manifest_sha256"] = layer_sha
    pointer["manuscript_eligibility_record_sha256"] = eligibility_sha
    write_json(CURRENT_POINTER, pointer)
    return {**final_hashes, "formalization_record_sha256": formalization_sha, "eligibility_record_sha256": eligibility_sha}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true", help="write v5 formalization and eligibility metadata")
    parser.add_argument("--refresh", action="store_true", help="refresh hashes in an already accepted v5 state")
    parser.add_argument("--recorded-at", help="UTC ISO-8601 author adoption timestamp")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.refresh:
        result = refresh_accepted_state()
        print(json.dumps({"status": "refreshed", **result}, ensure_ascii=False))
        return 0
    source = validate_inputs()
    if not args.apply:
        print(
            json.dumps(
                {
                    "status": "dry_run",
                    "manuscript_eligible": True,
                    "manuscript_eligibility_scope": "publication_clean_v5_display_layer_only",
                    "raw_numerical_status": "diagnostic_only",
                    "raw_manuscript_eligible": False,
                    "convergence_gate_status": "local_high_rate_gate_not_run",
                    "figure_count": 72,
                },
                ensure_ascii=False,
            )
        )
        return 0
    if V5_FORMALIZATION_RECORD.exists() or V5_ELIGIBILITY_RECORD.exists():
        raise FileExistsError("v5 formalization or eligibility record already exists; refusing a second promotion")
    recorded_at = args.recorded_at or dt.datetime.now(dt.timezone.utc).isoformat()
    result = update_manifests(source, recorded_at)
    print(json.dumps({"status": "accepted", "recorded_at_utc": recorded_at, **result}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
