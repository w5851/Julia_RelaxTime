"""Explicit adapter for the versioned Issue #130 phase-reference candidate.

The candidate tables are deliberately not shaped like the historical files in
``data/reference/pnjl``.  This module is the small, solver-free boundary between
those contracts.  It never writes the legacy files and it never selects the
candidate implicitly.  Callers must provide the candidate root and layer.

``load_phase_reference`` returns normalized in-memory rows.  ``to_legacy_views``
is an opt-in, strict conversion for consumers that still need the historical
column names; it refuses unresolved/non-certified rows by default so a caller
cannot silently turn diagnostic evidence into runtime input.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping


SCHEMA_VERSION = "pnjl_issue130_phase_reference_adapter_v1"
IMPORT_SCHEMA = "pnjl_issue130_phase_reference_import_v1"
LAYERS = ("strict", "derived", "render", "accepted")
IMPORT_SCHEMA_V2 = "pnjl_issue130_phase_reference_v2"

_TABLES: dict[str, dict[str, Any]] = {
    "boundary": {
        "strict": "maxwell_surface_strict_reference_v1.csv",
        "derived": "maxwell_surface_derived_reference_v1.csv",
        "render": "maxwell_surface_render.csv",
        "keys": ("xi", "T_MeV"),
        "numeric": ("xi", "T_MeV", "mu_MeV", "rho_hadron", "rho_quark", "area_residual"),
    },
    "crossover": {
        "strict": "crossover_surface_strict_reference_v1.csv",
        "derived": "crossover_surface_derived_reference_v1.csv",
        "render": "crossover_surface_render.csv",
        "keys": ("xi", "mu_MeV"),
        "numeric": ("xi", "mu_MeV", "T_MeV", "rho", "mu_CEP_proxy_MeV"),
    },
    "cep": {
        "strict": "cep_boundary_strict_reference_v1.csv",
        "derived": "cep_boundary_derived_reference_v1.csv",
        "render": "cep_boundary_render.csv",
        "keys": ("xi",),
        "numeric": ("xi", "mu_CEP_proxy_MeV", "T_low_MeV", "T_high_MeV", "T_midpoint_MeV"),
    },
    "spinodals": {
        "strict": "spinodal_surface_strict_reference_v1.csv",
        "derived": "spinodal_surface_derived_reference_v1.csv",
        "render": None,
        "keys": ("xi", "T_MeV"),
        "numeric": ("xi", "T_MeV", "mu_spinodal_hadron_MeV", "mu_spinodal_quark_MeV"),
    },
}

# v2 intentionally has no public ``derived`` layer.  Its render layer is a
# complete structured display package, and accepted is an explicit downstream
# candidate.  Keys and required numeric fields remain identical to v1 so the
# normalizer and legacy-view conversion can be reused without changing units.
_TABLES_V2: dict[str, dict[str, Any]] = {
    "boundary": {
        **_TABLES["boundary"],
        "render": "maxwell_surface_render.csv",
        "accepted": "maxwell_surface_accepted_phase_map_v1.csv",
    },
    "crossover": {
        **_TABLES["crossover"],
        "render": "crossover_surface_render.csv",
        "accepted": "crossover_surface_accepted_phase_map_v1.csv",
    },
    "cep": {
        **_TABLES["cep"],
        "render": "cep_boundary_render.csv",
        "accepted": "cep_boundary_accepted_phase_map_v1.csv",
    },
    "spinodals": {
        **_TABLES["spinodals"],
        "render": "spinodal_surface_render.csv",
        "accepted": "spinodal_surface_accepted_phase_map_v1.csv",
    },
}


class PhaseReferenceContractError(ValueError):
    """Raised when candidate data cannot satisfy the explicit adapter contract."""


@dataclass(frozen=True)
class PhaseReferenceBundle:
    """Validated, normalized candidate data for one explicitly selected layer."""

    root: Path
    layer: str
    manifest: Mapping[str, Any]
    layer_manifest: Mapping[str, Any]
    tables: Mapping[str, tuple[Mapping[str, Any], ...]]
    diagnostics: Mapping[str, Any]

    @property
    def runtime_consumption(self) -> bool:
        return bool(self.manifest.get("runtime_consumption", False))


@dataclass(frozen=True)
class PhaseReferenceRuntimeView:
    """Certified candidate rows merged with an explicit legacy fallback view."""

    layer: str
    source: str
    tables: Mapping[str, tuple[Mapping[str, Any], ...]]
    diagnostics: Mapping[str, Any]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise PhaseReferenceContractError(f"invalid JSON: {path}") from exc
    if not isinstance(value, dict):
        raise PhaseReferenceContractError(f"expected JSON object: {path}")
    return value


def default_downstream_layer(reference_root: str | Path) -> str:
    """Select the default layer for downstream display/analysis consumers.

    The v2 package has an author-accepted ``accepted`` view for downstream
    phase maps.  The v1 import remains compatible with ``strict``.  Runtime
    callers must continue to select strict explicitly through the runtime
    adapter; this helper is intentionally not used by solver paths.
    """

    manifest = _read_json(Path(reference_root).resolve() / "manifest.json")
    schema_version = manifest.get("schema_version")
    if schema_version == IMPORT_SCHEMA_V2:
        if (
            manifest.get("promotion_status") == "accepted_for_downstream"
            and manifest.get("downstream_default_layer") == "accepted"
        ):
            return "accepted"
        # A v2 package can exist before the author promotion step.  Keep the
        # safe strict view as the implicit analysis default until that step is
        # recorded in the root manifest.
        return "strict"
    if schema_version == IMPORT_SCHEMA:
        return "strict"
    raise PhaseReferenceContractError(
        f"unsupported candidate schema for downstream default: {schema_version!r}"
    )


def _float(value: str | None, *, path: Path, row_number: int, field: str) -> float:
    try:
        parsed = float(value if value is not None else "")
    except (TypeError, ValueError) as exc:
        raise PhaseReferenceContractError(
            f"non-numeric {field} at {path}:{row_number}"
        ) from exc
    if not math.isfinite(parsed):
        raise PhaseReferenceContractError(f"non-finite {field} at {path}:{row_number}")
    return parsed


def _bool(value: str | None, *, default: bool = False) -> bool:
    if value is None or value == "":
        return default
    normalized = value.strip().lower()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no"}:
        return False
    raise PhaseReferenceContractError(f"invalid boolean value: {value!r}")


def _optional_float(value: str | None, *, path: Path, row_number: int, field: str) -> float | None:
    if value is None or value == "":
        return None
    return _float(value, path=path, row_number=row_number, field=field)


def _read_table(path: Path, spec: Mapping[str, Any]) -> tuple[list[str], list[dict[str, str]]]:
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            header = list(reader.fieldnames or [])
            rows = list(reader)
    except OSError as exc:
        raise PhaseReferenceContractError(f"cannot read candidate table: {path}") from exc
    missing = [field for field in spec["numeric"] if field not in header]
    if missing:
        raise PhaseReferenceContractError(
            f"{path} is missing required fields: {', '.join(missing)}"
        )
    keys = tuple(spec["keys"])
    seen: set[tuple[str, ...]] = set()
    duplicates: list[tuple[str, ...]] = []
    for row_number, row in enumerate(rows, start=2):
        key = tuple(row.get(field, "") for field in keys)
        if key in seen:
            duplicates.append(key)
        seen.add(key)
        for field in spec["numeric"]:
            _float(row.get(field), path=path, row_number=row_number, field=field)
    if duplicates:
        raise PhaseReferenceContractError(
            f"duplicate {keys} in {path}: {duplicates[0]}"
        )
    return header, rows


def _is_interpolated(layer: str, row: Mapping[str, str]) -> bool:
    del layer  # The row's explicit layer/status is authoritative for all views.
    source_status = row.get("source_status", "").lower()
    return (
        row.get("layer") == "interpolated_noncertified"
        or "interpolat" in row.get("status", "").lower()
        or "interpolat" in row.get("interpolation_method", "").lower()
        or "interpolat" in source_status
        or "unresolved" in source_status
    )


def _maxwell_certified(layer: str, row: Mapping[str, str]) -> bool:
    if _is_interpolated(layer, row):
        return False
    if "grid_unresolved" in row and _bool(row.get("grid_unresolved")):
        return False
    if "geometry_converged" in row and not _bool(row.get("geometry_converged"), default=False):
        return False
    if "finite_and_converged" in row and not _bool(row.get("finite_and_converged"), default=False):
        return False
    status = row.get("status", "").lower()
    return not any(token in status for token in ("unresolved", "not_converged", "ambiguous"))


def _crossover_certified(layer: str, row: Mapping[str, str]) -> bool:
    if _is_interpolated(layer, row):
        return False
    if row.get("physical_region", "").lower() not in {"", "crossover_below_cep"}:
        return False
    status = row.get("status", "").lower()
    return "ambiguous" not in status and "unresolved" not in status


def _cep_certified(row: Mapping[str, str]) -> bool:
    if _is_interpolated("", row):
        return False
    width = _float(row.get("T_high_MeV"), path=Path("<row>"), row_number=0, field="T_high_MeV") - _float(
        row.get("T_low_MeV"), path=Path("<row>"), row_number=0, field="T_low_MeV"
    )
    return width >= 0.0 and width <= 0.1 and "ambiguous" not in row.get("status", "").lower()


def _normalize(table: str, layer: str, rows: Iterable[Mapping[str, str]]) -> tuple[dict[str, Any], ...]:
    normalized: list[dict[str, Any]] = []
    for row in rows:
        if table == "boundary":
            muq = float(row["mu_MeV"])
            item = {
                "surface": "maxwell",
                "xi": float(row["xi"]),
                "T_MeV": float(row["T_MeV"]),
                "muq_MeV": muq,
                "muB_MeV": 3.0 * muq,
                "mu_transition_MeV": muq,
                "rho_hadron": float(row["rho_hadron"]),
                "rho_quark": float(row["rho_quark"]),
                "area_residual": float(row["area_residual"]),
                "layer": row.get("layer", layer),
                "status": row.get("status", ""),
                "certified": _maxwell_certified(layer, row),
                "source_layer": row.get("source_layer", ""),
            }
        elif table == "crossover":
            muq = float(row["mu_MeV"])
            item = {
                "surface": "crossover",
                "xi": float(row["xi"]),
                "muq_MeV": muq,
                "muB_MeV": 3.0 * muq,
                "mu_MeV": muq,
                "T_MeV": float(row["T_MeV"]),
                "rho": float(row["rho"]),
                "mu_CEP_proxy_MeV": float(row["mu_CEP_proxy_MeV"]),
                "physical_region": row.get("physical_region", ""),
                "layer": row.get("layer", layer),
                "status": row.get("status", ""),
                "certified": _crossover_certified(layer, row),
                "source_layer": row.get("source_layer", ""),
            }
        elif table == "cep":
            muq = float(row["mu_CEP_proxy_MeV"])
            low = float(row["T_low_MeV"])
            high = float(row["T_high_MeV"])
            item = {
                "surface": "cep_boundary",
                "xi": float(row["xi"]),
                "muq_CEP_proxy_MeV": muq,
                "muB_CEP_proxy_MeV": 3.0 * muq,
                "T_low_MeV": low,
                "T_high_MeV": high,
                "T_midpoint_MeV": float(row["T_midpoint_MeV"]),
                "bracket_width_MeV": high - low,
                "layer": row.get("layer", layer),
                "status": row.get("status", ""),
                "boundary_mode": row.get("boundary_mode", "estimated_midpoint"),
                "certified": _cep_certified(row),
                "source_layer": row.get("source_layer", ""),
            }
        elif table == "spinodals":
            item = {
                "surface": "spinodal",
                "xi": float(row["xi"]),
                "T_MeV": float(row["T_MeV"]),
                "muq_spinodal_hadron_MeV": float(row["mu_spinodal_hadron_MeV"]),
                "muq_spinodal_quark_MeV": float(row["mu_spinodal_quark_MeV"]),
                "muB_spinodal_hadron_MeV": 3.0 * float(row["mu_spinodal_hadron_MeV"]),
                "muB_spinodal_quark_MeV": 3.0 * float(row["mu_spinodal_quark_MeV"]),
                "rho_spinodal_hadron": _optional_float(
                    row.get("rho_spinodal_hadron"), path=Path("<row>"), row_number=0, field="rho_spinodal_hadron"
                ),
                "rho_spinodal_quark": _optional_float(
                    row.get("rho_spinodal_quark"), path=Path("<row>"), row_number=0, field="rho_spinodal_quark"
                ),
                "layer": row.get("layer", layer),
                "status": row.get("status", ""),
                "certified": not _is_interpolated(layer, row)
                and "unresolved" not in row.get("status", "").lower(),
                "source_layer": row.get("source_layer", ""),
            }
        else:
            raise AssertionError(f"unknown phase-reference table: {table}")
        normalized.append(item)
    return tuple(normalized)


def load_phase_reference(
    reference_root: str | Path,
    *,
    layer: str = "strict",
    allow_runtime: bool = False,
) -> PhaseReferenceBundle:
    """Load and normalize one explicitly selected candidate layer.

    ``allow_runtime`` is intentionally opt-in and only verifies the contract;
    it does not switch any repository default.  It rejects unresolved or
    interpolated rows, and therefore cannot turn diagnostic evidence into a
    runtime reference accidentally.
    """

    root = Path(reference_root).resolve()
    if layer not in LAYERS:
        raise PhaseReferenceContractError(f"layer must be one of {LAYERS}, got {layer!r}")
    manifest_path = root / "manifest.json"
    if not manifest_path.is_file():
        raise PhaseReferenceContractError(f"missing candidate manifest: {manifest_path}")
    manifest = _read_json(manifest_path)
    schema_version = manifest.get("schema_version")
    if schema_version not in {IMPORT_SCHEMA, IMPORT_SCHEMA_V2}:
        raise PhaseReferenceContractError("candidate import schema mismatch")
    if schema_version == IMPORT_SCHEMA:
        table_specs = _TABLES
        allowed_layers = {"strict", "derived", "render"}
    else:
        table_specs = _TABLES_V2
        allowed_layers = {"strict", "render", "accepted"}
    if layer not in allowed_layers:
        raise PhaseReferenceContractError(
            f"layer {layer!r} is not available for schema {schema_version!r}"
        )
    layer_root = root / layer
    layer_manifest_path = layer_root / "manifest.json"
    if not layer_manifest_path.is_file():
        raise PhaseReferenceContractError(f"missing candidate manifest: {layer_manifest_path}")
    layer_manifest = _read_json(layer_manifest_path)
    if manifest.get("runtime_consumption") is not False:
        raise PhaseReferenceContractError("candidate manifest must keep runtime_consumption=false")
    if manifest.get("reference_status") not in {"candidate", "imported_candidate"}:
        raise PhaseReferenceContractError("candidate reference_status is not a candidate state")

    table_rows: dict[str, tuple[Mapping[str, Any], ...]] = {}
    row_counts: dict[str, int] = {}
    for table, spec in table_specs.items():
        filename = spec[layer]
        if filename is None:
            continue
        path = layer_root / "tables" / filename
        if not path.is_file():
            raise PhaseReferenceContractError(f"missing {layer} table: {path}")
        _, rows = _read_table(path, spec)
        table_rows[table] = _normalize(table, layer, rows)
        row_counts[table] = len(rows)

    uncertified = sum(
        1 for rows in table_rows.values() for row in rows if not bool(row.get("certified", False))
    )
    if allow_runtime and uncertified:
        raise PhaseReferenceContractError(
            f"runtime view rejected {uncertified} unresolved/non-certified candidate rows"
        )
    if allow_runtime and layer == "render":
        raise PhaseReferenceContractError("render layer is visualization-only and cannot be runtime input")
    if allow_runtime and layer == "accepted":
        raise PhaseReferenceContractError(
            "accepted layer is a downstream candidate and requires explicit promotion before runtime input"
        )

    diagnostics = {
        "adapter_schema": SCHEMA_VERSION,
        "reference_status": manifest.get("reference_status"),
        "runtime_consumption": False,
        "layer": layer,
        "schema_version": schema_version,
        "row_counts": row_counts,
        "uncertified_rows": uncertified,
        "allow_runtime": allow_runtime,
        "manifest_sha256": sha256(manifest_path),
        "layer_manifest_sha256": sha256(layer_manifest_path),
    }
    return PhaseReferenceBundle(root, layer, manifest, layer_manifest, table_rows, diagnostics)


def _require_certified(rows: Iterable[Mapping[str, Any]], table: str) -> None:
    bad = next((row for row in rows if not row.get("certified", False)), None)
    if bad is not None:
        raise PhaseReferenceContractError(
            f"{table} contains unresolved/non-certified row; choose diagnostic mode explicitly"
        )


def to_legacy_views(
    bundle: PhaseReferenceBundle,
    *,
    require_certified: bool = True,
    cep_mode: str = "estimated_midpoint",
) -> dict[str, tuple[dict[str, Any], ...]]:
    """Convert a bundle to historical table names without writing files.

    ``cep_mode='bracket'`` keeps the bracket fields and omits the midpoint
    aliases; ``estimated_midpoint`` explicitly exposes the midpoint plus the
    original bracket for consumers that need a scalar endpoint.
    """

    if cep_mode not in {"bracket", "estimated_midpoint"}:
        raise PhaseReferenceContractError("cep_mode must be 'bracket' or 'estimated_midpoint'")
    if require_certified:
        for table, rows in bundle.tables.items():
            _require_certified(rows, table)

    boundary = tuple(
        {
            "xi": row["xi"],
            "T_MeV": row["T_MeV"],
            "mu_transition_MeV": row["muq_MeV"],
            "rho_hadron": row["rho_hadron"],
            "rho_quark": row["rho_quark"],
            "area_residual": row["area_residual"],
            "converged": row["certified"],
        }
        for row in bundle.tables.get("boundary", ())
    )
    crossover = tuple(
        {
            "xi": row["xi"],
            "mu_MeV": row["muq_MeV"],
            "T_crossover_MeV": row["T_MeV"],
            "rho": row["rho"],
            "method": "issue130_candidate",
            "converged": row["certified"],
            "derivative": None,
            "variable": "mu_q",
        }
        for row in bundle.tables.get("crossover", ())
    )
    cep = []
    for row in bundle.tables.get("cep", ()):
        item = {
            "xi": row["xi"],
            "T_CEP_MeV": row["T_midpoint_MeV"] if cep_mode == "estimated_midpoint" else None,
            "muq_CEP_MeV": row["muq_CEP_proxy_MeV"] if cep_mode == "estimated_midpoint" else None,
            "muB_CEP_MeV": row["muB_CEP_proxy_MeV"] if cep_mode == "estimated_midpoint" else None,
            "T_low_MeV": row["T_low_MeV"],
            "T_high_MeV": row["T_high_MeV"],
            "T_midpoint_MeV": row["T_midpoint_MeV"],
            "bracket_width_MeV": row["bracket_width_MeV"],
            "boundary_mode": row["boundary_mode"],
        }
        cep.append(item)
    spinodals = tuple(
        {
            "xi": row["xi"],
            "T_MeV": row["T_MeV"],
            "mu_spinodal_hadron_MeV": row["muq_spinodal_hadron_MeV"],
            "mu_spinodal_quark_MeV": row["muq_spinodal_quark_MeV"],
            "rho_spinodal_hadron": row["rho_spinodal_hadron"],
            "rho_spinodal_quark": row["rho_spinodal_quark"],
        }
        for row in bundle.tables.get("spinodals", ())
    )
    return {
        "boundary": boundary,
        "crossover": crossover,
        "cep": tuple(cep),
        "spinodals": spinodals,
    }


def build_runtime_view(
    candidate: PhaseReferenceBundle,
    *,
    legacy_tables: Mapping[str, Iterable[Mapping[str, Any]]] | None = None,
) -> PhaseReferenceRuntimeView:
    """Build a solver-free certified-only view with per-key legacy fallback.

    This mirrors the Julia runtime switch without writing either source.  The
    candidate remains the preferred source; unresolved/interpolated rows are
    omitted and only keys absent from that certified view are filled from the
    caller-provided legacy tables.
    """

    if candidate.layer == "render":
        raise PhaseReferenceContractError(
            "render layer is visualization-only and cannot be a runtime view"
        )
    if candidate.layer == "accepted":
        raise PhaseReferenceContractError(
            "accepted layer requires an explicit downstream consumer mode"
        )
    legacy_tables = legacy_tables or {}
    merged: dict[str, tuple[Mapping[str, Any], ...]] = {}
    candidate_counts: dict[str, int] = {}
    fallback_counts: dict[str, int] = {}

    def key(table: str, row: Mapping[str, Any]) -> tuple[Any, ...]:
        spec = _TABLES[table]
        return tuple(row.get(field) for field in spec["keys"])

    for table in _TABLES:
        certified = [row for row in candidate.tables.get(table, ()) if row.get("certified", False)]
        rows = list(certified)
        seen = {key(table, row) for row in rows}
        n_fallback = 0
        for row in legacy_tables.get(table, ()):
            row_key = key(table, row)
            if row_key in seen:
                continue
            fallback = dict(row)
            fallback.update({"source_layer": "legacy_fallback", "status": "legacy_fallback", "certified": True})
            rows.append(fallback)
            seen.add(row_key)
            n_fallback += 1
        merged[table] = tuple(rows)
        candidate_counts[table] = len(certified)
        fallback_counts[table] = n_fallback

    return PhaseReferenceRuntimeView(
        layer=candidate.layer,
        source="candidate",
        tables=merged,
        diagnostics={
            "runtime_view": "certified_candidate_with_legacy_fallback",
            "candidate_manifest_sha256": candidate.diagnostics.get("manifest_sha256", ""),
            "candidate_row_counts": candidate_counts,
            "fallback_row_counts": fallback_counts,
            "fallback_enabled": True,
            "fallback_reason": "candidate_key_absent_or_uncertified",
        },
    )


__all__ = [
    "LAYERS",
    "PhaseReferenceBundle",
    "PhaseReferenceContractError",
    "PhaseReferenceRuntimeView",
    "SCHEMA_VERSION",
    "default_downstream_layer",
    "build_runtime_view",
    "load_phase_reference",
    "sha256",
    "to_legacy_views",
]
