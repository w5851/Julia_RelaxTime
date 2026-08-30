#!/usr/bin/env python3
"""Build the Issue #130 three-layer phase-reference candidate.

The public v2 package is intentionally organized as ``strict -> render ->
accepted``.  The existing v1 derived tables are an internal build input only:
they are retained byte-for-byte as provenance, while the render layer is made
complete by adding its missing spinodal table.  The accepted layer is a
deterministic downstream view that preserves native/unresolved/interpolated
status; it does not turn interpolation into strict certification.

This command is solver-free.  It reads CSV/JSON files, validates their keys and
numeric fields, copies the strict tables, completes the structured render
tables and writes an accepted candidate package.  It never calls Julia,
equilibrium, Maxwell or an oracle, and it refuses to overwrite an existing
output directory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


SCHEMA = "pnjl_issue130_phase_reference_v2"
SOURCE_SCHEMA = "pnjl_issue130_phase_reference_import_v1"
LAYER_SCHEMA = "pnjl_issue130_phase_reference_layer_v2"
CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
TABLES = ("boundary", "crossover", "cep", "spinodals")

TABLE_FILES: dict[str, dict[str, str]] = {
    "boundary": {
        "strict": "maxwell_surface_strict_reference_v1.csv",
        "derived": "maxwell_surface_derived_reference_v1.csv",
        "render": "maxwell_surface_render.csv",
        "accepted": "maxwell_surface_accepted_phase_map_v1.csv",
    },
    "crossover": {
        "strict": "crossover_surface_strict_reference_v1.csv",
        "derived": "crossover_surface_derived_reference_v1.csv",
        "render": "crossover_surface_render.csv",
        "accepted": "crossover_surface_accepted_phase_map_v1.csv",
    },
    "cep": {
        "strict": "cep_boundary_strict_reference_v1.csv",
        "derived": "cep_boundary_derived_reference_v1.csv",
        "render": "cep_boundary_render.csv",
        "accepted": "cep_boundary_accepted_phase_map_v1.csv",
    },
    "spinodals": {
        "strict": "spinodal_surface_strict_reference_v1.csv",
        "derived": "spinodal_surface_derived_reference_v1.csv",
        "render": "spinodal_surface_render.csv",
        "accepted": "spinodal_surface_accepted_phase_map_v1.csv",
    },
}

KEYS: dict[str, tuple[str, ...]] = {
    "boundary": ("xi", "T_MeV"),
    "crossover": ("xi", "mu_MeV"),
    "cep": ("xi",),
    "spinodals": ("xi", "T_MeV"),
}

NUMERIC_FIELDS: dict[str, tuple[str, ...]] = {
    "boundary": ("xi", "T_MeV", "mu_MeV", "rho_hadron", "rho_quark", "area_residual"),
    "crossover": ("xi", "mu_MeV", "T_MeV", "rho", "mu_CEP_proxy_MeV"),
    "cep": ("xi", "mu_CEP_proxy_MeV", "T_low_MeV", "T_high_MeV", "T_midpoint_MeV"),
    "spinodals": ("xi", "T_MeV", "mu_spinodal_hadron_MeV", "mu_spinodal_quark_MeV"),
}
COVERAGE_NUMERIC_FIELDS = ("xi", "xi_left", "xi_right", "axis_low", "axis_high", "source_rows")


class PhaseReferenceV2Error(ValueError):
    """Raised when the source cannot satisfy the v2 materialization contract."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_record(path: Path, *, root: Path | None = None, path_label: str | None = None) -> dict[str, Any]:
    if path_label is not None:
        display_path = path_label
    elif root is not None:
        display_path = path.relative_to(root).as_posix()
    else:
        display_path = str(path)
    record: dict[str, Any] = {
        "path": display_path,
        "sha256": sha256(path),
        "bytes": path.stat().st_size,
    }
    return record


def read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PhaseReferenceV2Error(f"invalid JSON: {path}") from exc
    if not isinstance(value, dict):
        raise PhaseReferenceV2Error(f"expected JSON object: {path}")
    return value


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    try:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            fields = list(reader.fieldnames or [])
            rows = list(reader)
    except OSError as exc:
        raise PhaseReferenceV2Error(f"cannot read CSV: {path}") from exc
    if not fields:
        raise PhaseReferenceV2Error(f"CSV has no header: {path}")
    if any(value is None for row in rows for value in row.values()):
        raise PhaseReferenceV2Error(f"CSV has a column-width mismatch: {path}")
    return fields, rows


def write_csv(path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def finite(value: str, *, path: Path, row_number: int, field: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise PhaseReferenceV2Error(
            f"non-numeric {field} at {path}:{row_number}"
        ) from exc
    if not math.isfinite(parsed):
        raise PhaseReferenceV2Error(f"non-finite {field} at {path}:{row_number}")
    return parsed


def validate_table(path: Path, table: str) -> tuple[list[str], list[dict[str, str]]]:
    fields, rows = read_csv(path)
    missing = [field for field in NUMERIC_FIELDS[table] if field not in fields]
    if missing:
        raise PhaseReferenceV2Error(
            f"{path} is missing required {table} fields: {', '.join(missing)}"
        )
    seen: set[tuple[str, ...]] = set()
    for row_number, row in enumerate(rows, start=2):
        key = tuple(row.get(field, "") for field in KEYS[table])
        if key in seen:
            raise PhaseReferenceV2Error(f"duplicate {table} key at {path}:{row_number}: {key}")
        seen.add(key)
        for field in NUMERIC_FIELDS[table]:
            finite(row.get(field, ""), path=path, row_number=row_number, field=field)
    return fields, rows


def validate_coverage(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Validate the mesh/coverage index without requiring a surface schema."""

    fields, rows = read_csv(path)
    missing = [field for field in ("surface", "xi", "coverage_status") if field not in fields]
    if missing:
        raise PhaseReferenceV2Error(
            f"{path} is missing required coverage fields: {', '.join(missing)}"
        )
    seen: set[tuple[str, str]] = set()
    for row_number, row in enumerate(rows, start=2):
        key = (row.get("surface", ""), row.get("xi", ""))
        if key in seen:
            raise PhaseReferenceV2Error(f"duplicate coverage key at {path}:{row_number}: {key}")
        seen.add(key)
        for field in COVERAGE_NUMERIC_FIELDS:
            value = row.get(field, "")
            if value not in (None, ""):
                finite(value, path=path, row_number=row_number, field=field)
    return fields, rows


def required_source_file(source_root: Path, layer: str, table: str) -> Path:
    path = source_root / layer / "tables" / TABLE_FILES[table][layer]
    if not path.is_file():
        raise PhaseReferenceV2Error(f"missing {layer} {table} table: {path}")
    return path


def copy_file(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)


def source_status(row: Mapping[str, str]) -> str:
    layer = row.get("layer", "").lower()
    status = row.get("status", "").lower()
    method = row.get("interpolation_method", "").lower()
    if "interpolat" in layer or "interpolat" in status or "interpolat" in method:
        return "interpolated_noncertified"
    if (
        "unresolved" in status
        or "not_converged" in status
        or row.get("grid_unresolved", "").strip().lower() in {"true", "1", "yes"}
        or row.get("geometry_converged", "").strip().lower() in {"false", "0", "no"}
        or row.get("finite_and_converged", "").strip().lower() in {"false", "0", "no"}
    ):
        return "unresolved"
    return "strict_certified"


def coverage_status(row: Mapping[str, str]) -> str:
    status = source_status(row)
    if status == "unresolved":
        return "native_unresolved_support"
    if status == "interpolated_noncertified":
        return "interpolated_common_support"
    return "native_support"


def accepted_rows(rows: Sequence[Mapping[str, str]]) -> tuple[list[str], list[dict[str, Any]]]:
    if not rows:
        raise PhaseReferenceV2Error("cannot build an accepted table from zero rows")
    fields = list(rows[0].keys())
    additions = [
        "source_status",
        "acceptance_status",
        "extrapolation",
        "coverage_status",
        "acceptance_scope",
    ]
    output_fields = fields + [field for field in additions if field not in fields]
    output: list[dict[str, Any]] = []
    for row in rows:
        item = dict(row)
        item.update(
            {
                "source_status": source_status(row),
                "acceptance_status": "candidate_pending_author_review",
                "extrapolation": "False",
                "coverage_status": coverage_status(row),
                "acceptance_scope": "downstream_phase_map_candidate",
            }
        )
        output.append(item)
    return output_fields, output


def layer_manifest(
    *,
    layer: str,
    semantics: str,
    source_root: Path,
    output_root: Path,
    table_records: Mapping[str, Any],
    generated_at: str,
    common: Mapping[str, Any],
    source_root_label: str,
    **extra: Any,
) -> dict[str, Any]:
    return {
        **common,
        "schema_version": LAYER_SCHEMA,
        "layer": layer,
        "semantics": semantics,
        "promotion_status": "awaiting_author_review",
        "reference_write": False,
        "runtime_consumption": False,
        "generated_at_utc": generated_at,
        "source_root": source_root_label,
        "tables": dict(table_records),
        **extra,
    }


def build_package(
    source_root: Path,
    output_root: Path,
    *,
    figure: Path | None = None,
    refresh_existing: bool = False,
) -> dict[str, Any]:
    """Materialize a v2 package from an immutable v1 candidate."""

    source_root = source_root.resolve()
    output_root = output_root.resolve()
    if output_root.exists() and not refresh_existing:
        raise FileExistsError(f"refusing to overwrite phase-reference v2 output: {output_root}")
    if output_root.exists() and refresh_existing:
        # This option is intentionally explicit and limited to a generated v2
        # candidate directory; it is not a general cleanup command.
        shutil.rmtree(output_root)
    source_manifest_path = source_root / "manifest.json"
    if not source_manifest_path.is_file():
        raise PhaseReferenceV2Error(f"missing source manifest: {source_manifest_path}")
    source_manifest = read_json(source_manifest_path)
    if source_manifest.get("schema_version") != SOURCE_SCHEMA:
        raise PhaseReferenceV2Error(
            f"source schema must be {SOURCE_SCHEMA}, got {source_manifest.get('schema_version')!r}"
        )
    if source_manifest.get("runtime_consumption") is not False:
        raise PhaseReferenceV2Error("source candidate must keep runtime_consumption=false")
    if source_manifest.get("calculation_sha") != CALCULATION_SHA:
        raise PhaseReferenceV2Error("source calculation SHA mismatch")

    source_files: dict[str, dict[str, Path]] = {layer: {} for layer in ("strict", "render", "derived")}
    source_data: dict[str, dict[str, list[dict[str, str]]]] = {layer: {} for layer in source_files}
    source_headers: dict[str, dict[str, list[str]]] = {layer: {} for layer in source_files}
    for table in TABLES:
        strict_path = required_source_file(source_root, "strict", table)
        strict_fields, strict_rows = validate_table(strict_path, table)
        source_files["strict"][table] = strict_path
        source_headers["strict"][table] = strict_fields
        source_data["strict"][table] = strict_rows

        render_path = required_source_file(source_root, "render", table) if table != "spinodals" else None
        if render_path is not None:
            render_fields, render_rows = validate_table(render_path, table)
            source_files["render"][table] = render_path
            source_headers["render"][table] = render_fields
            source_data["render"][table] = render_rows

        derived_path = required_source_file(source_root, "derived", table)
        derived_fields, derived_rows = validate_table(derived_path, table)
        source_files["derived"][table] = derived_path
        source_headers["derived"][table] = derived_fields
        source_data["derived"][table] = derived_rows

    mesh_path = source_root / "render" / "tables" / "mesh_coverage.csv"
    if not mesh_path.is_file():
        raise PhaseReferenceV2Error(f"missing render coverage table: {mesh_path}")
    mesh_fields, mesh_rows = validate_coverage(mesh_path)

    generated_at = datetime.now(timezone.utc).isoformat()
    common = {
        "calculation_sha": CALCULATION_SHA,
        "source_manifest": file_record(
            source_manifest_path,
            root=source_root,
            path_label="data/reference/pnjl/issue130_phase_reference_v1/manifest.json",
        ),
        "source_manifest_sha256": sha256(source_manifest_path),
        "source_run_id": source_manifest.get("source_run_id"),
        "replay_run_id": source_manifest.get("replay_run_id"),
        "uniform_xi_step": source_manifest.get("uniform_xi_step"),
        "uniform_xi_count": source_manifest.get("uniform_xi_count"),
        "solver_called": False,
        "oracle_labels_consumed": False,
    }
    output_root.mkdir(parents=True, exist_ok=True)

    figure_record: dict[str, Any] | None = None
    if figure is not None:
        figure = figure.resolve()
        if not figure.is_file():
            raise PhaseReferenceV2Error(f"figure path does not exist: {figure}")
        try:
            # Keep the provenance relocatable across worktrees.  The figure is
            # a repository asset, so its manifest path must never encode the
            # temporary checkout used to materialize this candidate.
            figure_record = file_record(
                figure,
                path_label=(
                    "data/outputs/figures/pnjl/phase_reference/"
                    "issue130_phase_reference_v1/phase_surface_render_mu_xi_T.png"
                ),
            )
        except ValueError:
            figure_record = file_record(figure)

    # The strict layer is a byte-preserving copy of v1 strict tables, with a
    # new v2 manifest that points back to every source file.
    strict_root = output_root / "strict"
    strict_table_records: dict[str, Any] = {}
    for table in TABLES:
        destination = strict_root / "tables" / TABLE_FILES[table]["strict"]
        copy_file(source_files["strict"][table], destination)
        strict_table_records[table] = {
            "source": file_record(source_files["strict"][table], root=source_root),
            "output": file_record(destination, root=output_root),
            "byte_preserved": sha256(destination) == sha256(source_files["strict"][table]),
            "rows": len(source_data["strict"][table]),
        }
    write_json(strict_root / "manifest.json", layer_manifest(
        layer="strict",
        semantics="native C2/v6 plus reviewed Maxwell additions; strict certificate and unresolved status preserved",
        source_root=source_root,
        output_root=output_root,
        table_records=strict_table_records,
        generated_at=generated_at,
        common=common,
        source_root_label="data/reference/pnjl/issue130_phase_reference_v1",
        source_layer="strict_reference_v1",
        strict_unresolved_preserved=True,
        rows={table: len(source_data["strict"][table]) for table in TABLES},
    ))
    (strict_root / "README.md").write_text(
        "# strict\n\n"
        "This v2 strict layer is a byte-preserving copy of the v1 strict tables. "
        "Its certificate and unresolved statuses are unchanged; author review of "
        "render/accepted data does not certify additional rows.\n",
        encoding="utf-8",
    )

    # Render keeps the existing three tables byte-for-byte and adds the missing
    # spinodal table from derived.  The source hash is explicit in the manifest.
    render_root = output_root / "render"
    render_table_records: dict[str, Any] = {}
    for table in ("boundary", "crossover", "cep"):
        destination = render_root / "tables" / TABLE_FILES[table]["render"]
        copy_file(source_files["render"][table], destination)
        render_table_records[table] = {
            "source": file_record(source_files["render"][table], root=source_root),
            "output": file_record(destination, root=output_root),
            "byte_preserved": sha256(destination) == sha256(source_files["render"][table]),
            "rows": len(source_data["render"][table]),
        }
    render_spinodal = render_root / "tables" / TABLE_FILES["spinodals"]["render"]
    copy_file(source_files["derived"]["spinodals"], render_spinodal)
    render_table_records["spinodals"] = {
        "source": file_record(source_files["derived"]["spinodals"], root=source_root),
        "output": file_record(render_spinodal, root=output_root),
        "byte_preserved": sha256(render_spinodal) == sha256(source_files["derived"]["spinodals"]),
        "rows": len(source_data["derived"]["spinodals"]),
        "source_layer": "derived_reference_v1",
    }
    render_mesh = render_root / "tables" / "mesh_coverage.csv"
    copy_file(mesh_path, render_mesh)
    render_table_records["coverage"] = {
        "source": file_record(mesh_path, root=source_root),
        "output": file_record(render_mesh, root=output_root),
        "byte_preserved": sha256(render_mesh) == sha256(mesh_path),
        "rows": len(mesh_rows),
    }
    render_manifest = layer_manifest(
        layer="render",
        semantics="complete structured display tables and author-reviewed figure provenance; no triangulation",
        source_root=source_root,
        output_root=output_root,
        table_records=render_table_records,
        generated_at=generated_at,
        common=common,
        source_root_label="data/reference/pnjl/issue130_phase_reference_v1",
        source_layer="phase_surface_render_v1_plus_derived_spinodal",
        surfaces=TABLES,
        complete_structured_tables=True,
        figure_review_status="existing_public_figure_reviewed; spinodal_table_completion_pending_review",
        runtime_input=False,
        figure=figure_record,
        rows={
            table: len(source_data["render"].get(table, source_data["derived"][table]))
            for table in TABLES
        },
    )
    write_json(render_root / "manifest.json", render_manifest)
    (render_root / "README.md").write_text(
        "# render\n\n"
        "This is the complete structured render layer. It preserves the existing "
        "crossover, Maxwell/first-order and CEP tables and adds the spinodal "
        "coordinate table from the immutable derived build input. The existing "
        "public figure remains referenced by provenance; no triangulation or "
        "strict-certificate promotion occurs here.\n",
        encoding="utf-8",
    )

    # Accepted is a data-level view, not a PNG-derived table.  It preserves all
    # source columns and appends explicit downstream acceptance metadata.
    accepted_root = output_root / "accepted"
    accepted_table_records: dict[str, Any] = {}
    accepted_counts: dict[str, int] = {}
    for table in TABLES:
        render_source = render_spinodal if table == "spinodals" else render_root / "tables" / TABLE_FILES[table]["render"]
        fields, rows = validate_table(render_source, table)
        accepted_fields, accepted = accepted_rows(rows)
        destination = accepted_root / "tables" / TABLE_FILES[table]["accepted"]
        write_csv(destination, accepted_fields, accepted)
        accepted_table_records[table] = {
            "source": file_record(render_source, root=output_root),
            "output": file_record(destination, root=output_root),
            "rows": len(accepted),
            "source_status_counts": {
                state: sum(1 for row in accepted if row["source_status"] == state)
                for state in ("strict_certified", "unresolved", "interpolated_noncertified")
            },
        }
        accepted_counts[table] = len(accepted)
    accepted_coverage = accepted_root / "tables" / "surface_coverage_mask.csv"
    copy_file(render_mesh, accepted_coverage)
    accepted_table_records["coverage"] = {
        "source": file_record(render_mesh, root=output_root),
        "output": file_record(accepted_coverage, root=output_root),
        "byte_preserved": sha256(accepted_coverage) == sha256(render_mesh),
        "rows": len(mesh_rows),
    }
    accepted_manifest = layer_manifest(
        layer="accepted",
        semantics="author-review candidate downstream phase map derived from complete structured render tables",
        source_root=source_root,
        output_root=output_root,
        table_records=accepted_table_records,
        generated_at=generated_at,
        common=common,
        source_root_label="data/reference/pnjl/issue130_phase_reference_v1",
        source_layer="render",
        promotion_status="awaiting_author_review",
        acceptance_status="candidate_pending_author_review",
        acceptance_scope="downstream_phase_map_only",
        interpolation_is_noncertified=True,
        strict_certificate_preserved=True,
        extrapolation=False,
        rows=accepted_counts,
    )
    write_json(accepted_root / "manifest.json", accepted_manifest)
    (accepted_root / "README.md").write_text(
        "# accepted\n\n"
        "This is a candidate downstream phase-map view generated from structured "
        "render tables. Native, unresolved and interpolated rows remain explicitly "
        "labelled; `candidate_pending_author_review` is not strict certification. "
        "The layer is opt-in and keeps `runtime_consumption=false` until a separate "
        "consumer audit and author promotion decision.\n",
        encoding="utf-8",
    )

    claims = [
        {
            "claim_id": "strict_byte_preserved",
            "claim": "v2 strict tables are byte-preserving copies of v1 strict tables",
            "status": "supported",
            "evidence": "strict/manifest.json",
            "boundary": "certificate semantics unchanged",
        },
        {
            "claim_id": "render_spinodal_complete",
            "claim": "render contains structured crossover, Maxwell, CEP and spinodal tables",
            "status": "supported",
            "evidence": "render/manifest.json; render/tables/spinodal_surface_render.csv",
            "boundary": "spinodal values are copied from derived build input; no solver call",
        },
        {
            "claim_id": "accepted_status_preserved",
            "claim": "accepted preserves native/unresolved/interpolated source status",
            "status": "supported",
            "evidence": "accepted/manifest.json; accepted/tables/*.csv",
            "boundary": "interpolated rows remain non-certified and author review is pending",
        },
        {
            "claim_id": "runtime_promotion",
            "claim": "v2 is an automatic replacement for strict runtime or legacy fallback",
            "status": "blocked",
            "evidence": "manifest.json",
            "boundary": "requires explicit consumer audit and author promotion",
        },
    ]
    claim_path = output_root / "claim_ledger.json"
    write_json(claim_path, {"schema_version": "pnjl_issue130_phase_reference_claim_ledger_v2", "claims": claims})

    root_manifest = {
        "schema_version": SCHEMA,
        "generated_at_utc": generated_at,
        "reference_status": "candidate",
        "promotion_status": "awaiting_author_review",
        "public_layers": ["strict", "render", "accepted"],
        "internal_build_inputs": ["data/reference/pnjl/issue130_phase_reference_v1/derived"],
        "reference_write": False,
        "runtime_consumption": False,
        **common,
        "layers": {
            "strict": file_record(strict_root / "manifest.json", root=output_root),
            "render": file_record(render_root / "manifest.json", root=output_root),
            "accepted": file_record(accepted_root / "manifest.json", root=output_root),
        },
        "build_input": {
            "derived_layer": {
                "path": "data/reference/pnjl/issue130_phase_reference_v1/derived",
                "manifest": file_record(
                    source_root / "derived" / "manifest.json",
                    root=source_root,
                ),
                "role": "internal_build_input_only",
            },
            "source_v1_manifest": file_record(source_manifest_path, root=source_root),
        },
        "figure": figure_record,
        "constraints": {
            "v1_unchanged": True,
            "derived_public_layer": False,
            "strict_certificate_preserved": True,
            "interpolated_rows_noncertified": True,
            "render_has_spinodal_table": True,
            "solver_called": False,
            "oracle_labels_consumed": False,
            "runtime_switch": False,
            "legacy_fallback_preserved": True,
            "extrapolation": False,
        },
        "claim_ledger": file_record(claim_path, root=output_root),
    }
    write_json(output_root / "manifest.json", root_manifest)
    (output_root / "README.md").write_text(
        "# Issue #130 phase-reference v2 candidate\n\n"
        "The public semantic layers are `strict`, `render` and `accepted`. The "
        "former v1 `derived` layer is retained only as an internal build input and "
        "provenance record. `render` is a complete structured display package; "
        "`accepted` is an opt-in downstream phase-map candidate. All layers are "
        "solver-free and keep `runtime_consumption=false` pending author review.\n",
        encoding="utf-8",
    )
    strict_counts_text = ", ".join(
        f"{table}={len(source_data['strict'][table])}" for table in TABLES
    )
    render_counts_text = ", ".join(
        f"{table}={len(source_data['render'].get(table, source_data['derived'][table]))}"
        for table in TABLES
    )
    accepted_counts_text = ", ".join(
        f"{table}={accepted_counts[table]}" for table in TABLES
    )
    (output_root / "AUDIT.md").write_text(
        "# Issue #130 phase-reference v2 solver-free materialization audit\n\n"
        "- source candidate: `data/reference/pnjl/issue130_phase_reference_v1`\n"
        f"- calculation SHA: `{CALCULATION_SHA}`\n"
        f"- strict rows: {strict_counts_text}\n"
        f"- render rows: {render_counts_text}\n"
        f"- accepted rows: {accepted_counts_text}\n"
        "- render spinodal is copied from the immutable derived build input because v1 render had no spinodal coordinate table.\n"
        "- accepted rows retain `strict_certified`, `unresolved` and `interpolated_noncertified`; no status is upgraded.\n"
        "- `solver_called=false`, `oracle_labels_consumed=false`, `reference_write=false`, `runtime_consumption=false`.\n\n"
        "## Boundary\n\n"
        "This package is a candidate for author-reviewed downstream research. It is not a strict convergence certificate and does not authorize runtime switching or PNJL legacy physical deletion.\n",
        encoding="utf-8",
    )
    return root_manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-root",
        type=Path,
        default=Path("data/reference/pnjl/issue130_phase_reference_v1"),
        help="immutable v1 candidate root",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("data/reference/pnjl/issue130_phase_reference_v2"),
        help="new v2 candidate root; must not already exist",
    )
    parser.add_argument(
        "--figure",
        type=Path,
        default=None,
        help="optional author-reviewed figure to record by path/hash",
    )
    parser.add_argument(
        "--refresh-existing",
        action="store_true",
        help="refresh only the generated v2 output directory",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = build_package(
        args.source_root,
        args.output_root,
        figure=args.figure,
        refresh_existing=args.refresh_existing,
    )
    print(json.dumps({
        "output_root": str(Path(args.output_root).resolve()),
        "schema_version": result["schema_version"],
        "layers": list(result["layers"]),
        "solver_called": result["solver_called"],
        "runtime_consumption": result["runtime_consumption"],
    }, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
