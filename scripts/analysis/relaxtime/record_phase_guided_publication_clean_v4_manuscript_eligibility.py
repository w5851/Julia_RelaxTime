#!/usr/bin/env python3
"""Record author acceptance of publication_clean_v4 for display-layer manuscript use.

This does not promote source raw numerical data, run a solver, or run a
high-rate convergence gate.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
from pathlib import Path
from typing import Any

import formalize_phase_guided_publication_clean_v4 as formalization


ROOT = formalization.ROOT
TRANSPORT_ROOT = formalization.TRANSPORT_ROOT
V4_ROOT = formalization.V4_ANALYSIS_ROOT
V4_PUBLIC_ROOT = formalization.V4_PUBLIC_ROOT
V4_LAYER_ROOT = formalization.V4_LAYER_ROOT
V4_PACKAGE_MANIFEST = formalization.V4_PACKAGE_MANIFEST
V4_PLOT_MANIFEST = formalization.V4_PLOT_MANIFEST
V4_PUBLIC_PLOT_MANIFEST = formalization.V4_PUBLIC_PLOT_MANIFEST
V4_LAYER_MANIFEST = formalization.V4_LAYER_MANIFEST
FORMALIZATION_RECORD = formalization.FORMALIZATION_RECORD
CURRENT_POINTER = formalization.CURRENT_POINTER
ELIGIBILITY_RECORD = TRANSPORT_ROOT / "publication_clean_v4_manuscript_eligibility_v1.json"

ELIGIBILITY_FIELDS = {
    "numerical_status": "author_accepted_display_only",
    "manuscript_eligible": True,
    "manuscript_eligibility_scope": "publication_clean_v4_display_layer_only",
    "raw_numerical_status": "diagnostic_only",
    "raw_manuscript_eligible": False,
    "convergence_gate_status": "local_high_rate_gate_not_run",
}


def sha256_file(path: Path) -> str:
    import hashlib

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relpath(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def verify_source_inputs(package: dict[str, Any]) -> None:
    filenames = {
        "scan_sha256": "phase_guided_transport_scan.csv",
        "diagnostics_sha256": "channel_diagnostics.csv",
        "failed_sha256": "failed_points.csv",
        "manifest_sha256": "manifest.json",
        "effective_config_sha256": "effective_config.json",
    }
    for item in package.get("source_inputs", []):
        case_root = (
            ROOT
            / "data"
            / "outputs"
            / "results"
            / "relaxtime"
            / "transport"
            / "phase_guided"
            / str(item["mode"])
            / str(package["case"])
        )
        for field, filename in filenames.items():
            source = case_root / filename
            if not source.is_file() or sha256_file(source) != item.get(field):
                raise ValueError(f"v4 source input hash mismatch: {relpath(source)}")


def verify_package_outputs(package: dict[str, Any]) -> None:
    for item in package.get("outputs", []):
        path = ROOT / Path(str(item["path"]))
        if not path.is_file():
            raise FileNotFoundError(path)
        if int(item["bytes"]) != path.stat().st_size or sha256_file(path) != item["sha256"]:
            raise ValueError(f"v4 package output hash mismatch: {relpath(path)}")


def verify_adjustment_scope() -> list[dict[str, str]]:
    path = V4_ROOT / "tables" / "v4_display_adjustment_map.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 4:
        raise ValueError(f"expected 4 disclosed v4 display adjustments, found {len(rows)}")
    expected = {
        "v4_mode_b_T200p0_muB0p0_xip0p36_composite": 3,
        "v4_mode_a_muB900p0_alpha1p0_xim0p003_tau_sbar": 1,
    }
    actual: dict[str, int] = {}
    for row in rows:
        key = row["adjustment_id"]
        actual[key] = actual.get(key, 0) + 1
        for field in ("raw_value", "derived_v4_display_value", "method", "phase_policy"):
            if not row.get(field):
                raise ValueError(f"v4 adjustment {key} has no {field}")
        if not math.isfinite(float(row["raw_value"])) or not math.isfinite(float(row["derived_v4_display_value"])):
            raise ValueError(f"v4 adjustment {key} contains a non-finite value")
    if actual != expected:
        raise ValueError(f"unexpected v4 adjustment inventory: {actual}")
    return rows


def verify_gate_boundary() -> None:
    path = (
        TRANSPORT_ROOT
        / "phase_guided_transport_publication_clean_v4_mechanism_review"
        / "tables"
        / "mechanism_manifest.json"
    )
    manifest = formalization.read_json(path)
    if manifest.get("convergence_gate_enabled") is not False:
        raise ValueError("the recorded local high-rate gate state changed")
    gate_path = path.with_name("local_rate_convergence_gate.csv")
    with gate_path.open(newline="", encoding="utf-8") as handle:
        if list(csv.DictReader(handle)):
            raise ValueError("local high-rate gate has rows; review before recording this decision")


def replace_readme_block(path: Path, old: str, new: str) -> None:
    raw = path.read_bytes()
    text = raw.decode("utf-8")
    newline = "\r\n" if b"\r\n" in raw else "\n"
    old = old.replace("\n", newline)
    new = new.replace("\n", newline)
    if text.count(old) != 1:
        raise ValueError(f"expected one status block in {relpath(path)}")
    path.write_bytes(text.replace(old, new, 1).encode("utf-8"))


def update_readmes() -> None:
    replace_readme_block(
        V4_ROOT / "README.md",
        """The numerical layer remains `diagnostic_only`: local display interpolation and
one-sided endpoint extrapolation are not new equilibrium or transport
solutions, and they are not a convergence certificate.  The figure layer is
accepted as the formal publication layout.  `manuscript_eligible=false` on
the numerical manifest therefore refers to raw numerical claims, not to the
author acceptance of these display figures.""",
        """The v4 display derivative is `author_accepted_display_only` and is
`manuscript_eligible=true` for the disclosed display layer.  Source raw
numerical data remain `diagnostic_only` and `raw_manuscript_eligible=false`;
this decision does not certify a new solver result or convergence.  The local
high-rate gate was not run and is not represented as passed.

Manuscript eligibility covers the 72 v4 figures and the explicit raw/display
provenance in `tables/publication_clean_points.csv` and
`tables/v4_display_adjustment_map.csv`.  The mode-B composite mechanism
remains unassessed.  The mode-A `tau_sbar` endpoint remains a display
extrapolation, not a replacement branch solution or a basis for an exact jump
amplitude.  The first-order gap remains unfilled.

The eligibility decision is recorded in
`../publication_clean_v4_manuscript_eligibility_v1.json`.""",
    )
    replace_readme_block(
        V4_LAYER_ROOT / "README.md",
        """The layer is author-accepted as the current formal publication layout.
`manuscript_eligible=false` still applies to raw numerical claims: the
underlying layer remains solver-free and numerically diagnostic-only. The
mechanism review is retained as diagnostic context; it is not a propagator
regularization or a production convergence certificate.""",
        """The layer is author-accepted and `manuscript_eligible=true` for display
use with raw/display values and local adjustments disclosed.  Source raw
numerical data remain `diagnostic_only` and are not promoted by this decision.
No local high-rate convergence gate was run; the figure status does not imply
that such a gate passed.

The mode-B composite adjustment has no separate mechanism verdict in the v4
record.  The `tau_sbar` endpoint remains a display extrapolation and cannot be
used as a replacement branch solution or an exact jump amplitude.

Eligibility record:

    docs/analysis/relaxtime/phase_guided_transport/publication_clean_v4_manuscript_eligibility_v1.json""",
    )


def update_claim_ledger() -> None:
    path = V4_ROOT / "tables" / "claim_ledger.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = list(reader.fieldnames or [])
        rows = list(reader)
    for row in rows:
        if row.get("claim_id") == "PC-V4-001":
            row["scope_limit"] = (
                "manuscript eligibility applies only to the author-accepted v4 display layer; "
                "source raw numerics remain diagnostic_only and are not promoted."
            )
        elif row.get("claim_id") == "PC-V4-005":
            row["claim_zh"] = (
                "作者已接受 v4 显示层用于稿件，并要求显式提供 raw/display 数值及局部调整记录；"
                "该资格仅覆盖显示派生层。"
            )
            row["evidence"] = (
                "publication_clean_v4_manuscript_eligibility_v1.json; README.md; "
                "tables/publication_clean_points.csv; tables/v4_display_adjustment_map.csv"
            )
            row["scope_limit"] = (
                "source raw numerical status remains diagnostic_only; no new solver run or "
                "local high-rate convergence gate is claimed."
            )
    if not any(row.get("claim_id") == "PC-V4-006" for row in rows):
        rows.append(
            {
                "claim_id": "PC-V4-006",
                "status": "author_accepted_with_scope_limit",
                "claim_zh": (
                    "publication_clean_v4 的 72 张图及其明确标注的显示派生值可用于稿件；"
                    "四个局部调整均保留原始值、显示值、锚点和方法。"
                ),
                "evidence": (
                    "publication_clean_v4_manuscript_eligibility_v1.json; "
                    "tables/publication_clean_points.csv; "
                    "tables/v4_display_adjustment_map.csv"
                ),
                "scope_limit": (
                    "raw numerical results remain diagnostic_only; mode-B composite mechanism "
                    "is unassessed; the local high-rate gate was not run."
                ),
            }
        )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def apply_eligibility(payload: dict[str, Any], record_path: str) -> None:
    payload.update(ELIGIBILITY_FIELDS)
    payload["manuscript_eligibility_record"] = record_path
    payload["eligibility_decision_basis"] = "explicit_author_acceptance_of_v4_display_layer"
    payload["convergence_claim"] = False


def refresh_package_outputs(package: dict[str, Any]) -> None:
    records = []
    for item in package.get("outputs", []):
        path = ROOT / Path(str(item["path"]))
        if not path.is_file():
            raise FileNotFoundError(path)
        records.append({"path": relpath(path), "sha256": sha256_file(path), "bytes": path.stat().st_size})
    package["outputs"] = records


def load_and_validate() -> dict[str, Any]:
    if ELIGIBILITY_RECORD.exists():
        raise FileExistsError(f"eligibility already recorded: {ELIGIBILITY_RECORD}")
    source = formalization.validate_inputs()
    pointer = formalization.read_json(CURRENT_POINTER)
    initial_record = formalization.read_json(FORMALIZATION_RECORD)
    package = source["v4_package"]
    plot = source["v4_plot"]
    public = source["v4_public"]
    layer = source["v4_layer"]
    if pointer.get("current_analysis_package") != relpath(V4_ROOT):
        raise ValueError("v4 is not the current publication package")
    if pointer.get("manuscript_eligible") is not False or pointer.get("numerical_status") != "diagnostic_only":
        raise ValueError("current pointer is not at the expected pre-eligibility state")
    if initial_record.get("status") != "author_accepted_formal_publication_figure_layer":
        raise ValueError("initial v4 figure acceptance record is missing")
    for manifest in (package, plot, public, layer):
        if manifest.get("manuscript_eligible") is not False:
            raise ValueError("v4 is not at the expected pre-eligibility state")
        if manifest.get("numerical_status") != "diagnostic_only":
            raise ValueError("v4 source numerical status changed")
    final = initial_record.get("final_hashes", {})
    expected_hashes = {
        "package_manifest_sha256": source["v4_package_sha256_before"],
        "plot_manifest_sha256": source["v4_plot_sha256_before"],
        "public_plot_manifest_sha256": source["v4_public_plot_sha256_before"],
        "figure_layer_manifest_sha256": source["v4_layer_sha256_before"],
    }
    for key, value in expected_hashes.items():
        if final.get(key) != value:
            raise ValueError(f"initial formalization hash mismatch: {key}")
    verify_package_outputs(package)
    verify_source_inputs(package)
    adjustments = verify_adjustment_scope()
    verify_gate_boundary()
    return {
        "source": source,
        "pointer": pointer,
        "initial_record": initial_record,
        "package": package,
        "plot": plot,
        "public": public,
        "layer": layer,
        "adjustments": adjustments,
    }


def record_document(
    source: dict[str, Any],
    initial_record: dict[str, Any],
    adjustments: list[dict[str, str]],
    final_hashes: dict[str, str],
    recorded_at: str,
) -> dict[str, Any]:
    return {
        "schema": "phase_guided_transport_publication_clean_v4_manuscript_eligibility_v1",
        "status": "author_accepted_display_only_manuscript_eligible",
        "recorded_at_utc": recorded_at,
        "decision_basis": "explicit_author_acceptance_in_computational_project",
        "current_analysis_package": relpath(V4_ROOT),
        "current_public_figure_root": relpath(V4_PUBLIC_ROOT),
        "current_figure_layer_manifest": relpath(V4_LAYER_MANIFEST),
        "initial_figure_acceptance_record": relpath(FORMALIZATION_RECORD),
        "initial_figure_acceptance_record_sha256": sha256_file(FORMALIZATION_RECORD),
        "initial_v4_hashes": {
            "package_manifest_sha256": source["v4_package_sha256_before"],
            "plot_manifest_sha256": source["v4_plot_sha256_before"],
            "public_plot_manifest_sha256": source["v4_public_plot_sha256_before"],
            "figure_layer_manifest_sha256": source["v4_layer_sha256_before"],
        },
        "final_hashes": final_hashes,
        "numerical_status": ELIGIBILITY_FIELDS["numerical_status"],
        "manuscript_eligible": True,
        "manuscript_eligibility_scope": ELIGIBILITY_FIELDS["manuscript_eligibility_scope"],
        "raw_numerical_status": "diagnostic_only",
        "raw_manuscript_eligible": False,
        "source_case": "first_canonical_v2_p128_xi001_onshellkernel_validated_anchored_prod_v2",
        "source_production_status": "diagnostic_only",
        "source_solver_called": True,
        "display_derivative_solver_called": False,
        "production_write": False,
        "raw_results_modified": False,
        "production_registry_modified": False,
        "adjustment_count": len(adjustments),
        "adjustments": adjustments,
        "data_disclosure": {
            "point_table": relpath(V4_ROOT / "tables" / "publication_clean_points.csv"),
            "adjustment_map": relpath(V4_ROOT / "tables" / "v4_display_adjustment_map.csv"),
            "adjustment_recipe": (
                "docs/analysis/relaxtime/phase_guided_transport/"
                "phase_guided_transport_v4_residual_smoothing/tables/"
                "publication_clean_v4_display_adjustments.csv"
            ),
            "raw_value_and_display_value_retained": True,
        },
        "accepted_display_rules": [
            "mode-B T=200 MeV, muB=0 MeV, xi=0.36: linear interpolation for zeta, sigma/T, and sigma using xi=0.35 and 0.37",
            "mode-A muB=900 MeV, alpha_T=1.0, xi=-0.003: left-branch log-linear display extrapolation using xi=-0.02 and -0.01",
            "first-order gap [-0.003,+0.003] remains split; no cross-phase interpolation",
        ],
        "numerical_boundaries": [
            "eligibility covers the author-accepted v4 display derivative and its disclosed value table only",
            "source raw numerical data remain diagnostic_only and are not promoted to production_grade",
            "the local high-rate convergence gate was not run; this decision is not a gate pass",
            "accepted-primary p104-to-p128 evidence uses a different accepted/runtime reference contract and is not direct validation of the v4 source case",
            "the three mode-B composite adjustments have mechanism_status=inherited_parent_scope and mechanism_evidence=not_assessed",
            "the tau_sbar endpoint display value is not a new branch solution or an exact first-order jump value",
        ],
        "initial_figure_acceptance_status": initial_record.get("status"),
    }


def update_manifests(state: dict[str, Any], recorded_at: str) -> dict[str, Any]:
    source = state["source"]
    package, plot = state["package"], state["plot"]
    public, layer = state["public"], state["layer"]
    record_path = relpath(ELIGIBILITY_RECORD)
    for manifest in (package, plot, public, layer):
        apply_eligibility(manifest, record_path)
    package["known_boundaries"] = list(dict.fromkeys([
        *package.get("known_boundaries", []),
        "manuscript eligibility applies to the author-accepted v4 display layer only",
        "source raw numerical status remains diagnostic_only",
        "the local high-rate convergence gate was not run and is not claimed as passed",
        "the mode-B composite mechanism remains unassessed",
    ]))
    package["author_manuscript_acceptance"] = {
        "status": "accepted",
        "scope": ELIGIBILITY_FIELDS["manuscript_eligibility_scope"],
        "recorded_at_utc": recorded_at,
        "raw_results_unchanged": True,
        "local_high_rate_gate_run": False,
    }
    plot["formalization_record"] = relpath(FORMALIZATION_RECORD)
    public["formalization_record"] = relpath(FORMALIZATION_RECORD)
    update_readmes()
    update_claim_ledger()

    write_json(V4_PLOT_MANIFEST, plot)
    analysis_plot_sha = sha256_file(V4_PLOT_MANIFEST)
    public["source_analysis_plot_manifest"] = relpath(V4_PLOT_MANIFEST)
    public["source_analysis_plot_manifest_sha256"] = analysis_plot_sha
    write_json(V4_PUBLIC_PLOT_MANIFEST, public)
    public_plot_sha = sha256_file(V4_PUBLIC_PLOT_MANIFEST)

    refresh_package_outputs(package)
    write_json(V4_PACKAGE_MANIFEST, package)
    package_sha = sha256_file(V4_PACKAGE_MANIFEST)

    layer["source_analysis_package_manifest_sha256"] = package_sha
    layer["source_plot_manifest_sha256"] = analysis_plot_sha
    layer["target_plot_manifest_sha256"] = public_plot_sha
    layer["generated_at_utc"] = recorded_at
    write_json(V4_LAYER_MANIFEST, layer)
    layer_sha = sha256_file(V4_LAYER_MANIFEST)

    final_hashes = {
        "package_manifest_sha256": package_sha,
        "plot_manifest_sha256": analysis_plot_sha,
        "public_plot_manifest_sha256": public_plot_sha,
        "figure_layer_manifest_sha256": layer_sha,
        "claim_ledger_sha256": sha256_file(V4_ROOT / "tables" / "claim_ledger.csv"),
    }
    decision = record_document(
        source, state["initial_record"], state["adjustments"], final_hashes, recorded_at
    )
    write_json(ELIGIBILITY_RECORD, decision)
    record_sha = sha256_file(ELIGIBILITY_RECORD)

    pointer = state["pointer"]
    apply_eligibility(pointer, relpath(ELIGIBILITY_RECORD))
    pointer["manuscript_eligibility_record_sha256"] = record_sha
    pointer["eligibility_decision_basis"] = "explicit_author_acceptance_of_v4_display_layer"
    pointer["convergence_claim"] = False
    pointer["manifest_sha256"] = package_sha
    pointer["plot_manifest_sha256"] = analysis_plot_sha
    pointer["public_plot_manifest_sha256"] = public_plot_sha
    pointer["figure_layer_manifest_sha256"] = layer_sha
    write_json(CURRENT_POINTER, pointer)
    return {"final_hashes": final_hashes, "eligibility_record_sha256": record_sha}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true", help="write display-layer eligibility metadata")
    parser.add_argument("--recorded-at", help="UTC ISO-8601 author acceptance timestamp")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    state = load_and_validate()
    if not args.apply:
        print(json.dumps({
            "status": "dry_run",
            **ELIGIBILITY_FIELDS,
            "adjustment_count": len(state["adjustments"]),
            "local_high_rate_gate": "not_run_and_not_claimed",
        }, ensure_ascii=False))
        return 0
    recorded_at = args.recorded_at or dt.datetime.now(dt.timezone.utc).isoformat()
    result = update_manifests(state, recorded_at)
    print(json.dumps({"status": "accepted", **result}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
