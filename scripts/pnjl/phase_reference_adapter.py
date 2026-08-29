"""Explicit adapter for the versioned Issue #130 phase-reference candidate.

The candidate tables are deliberately not shaped like the historical files in
``data/reference/pnjl``.  This module is the small, solver-free boundary between
those contracts.  It never writes the legacy files.  Runtime callers must
provide the candidate root and layer explicitly; display/analysis callers may
use ``default_downstream_layer`` for the author-promoted v2 package.

``load_phase_reference`` returns normalized in-memory rows.  ``to_legacy_views``
is an opt-in, strict conversion for consumers that still need the historical
column names; it refuses unresolved/non-certified rows by default so a caller
cannot silently turn diagnostic evidence into runtime input.  The default
runtime consumer is the author-accepted v2 layer.  The strict layer remains an
explicit certified-only mode; the historical legacy snapshot is not a runtime
fallback or rollback path.
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
    """Rows exposed to a runtime consumer with an auditable source policy."""

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
    phase maps and runtime consumers.  The v1 import remains compatible with
    ``strict``.  ``strict`` is still available as an explicit runtime choice.
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


def _bool(value: str | bool | None, *, default: bool = False) -> bool:
    if value is None or value == "":
        return default
    if isinstance(value, bool):
        return value
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
        # Keep the acceptance and support metadata available to the runtime
        # merger.  ``certified`` remains the strict numerical certificate;
        # ``runtime_eligible`` is only enabled for a row after an explicit
        # fallback policy admits it.
        item.update(
            {
                "source_status": row.get("source_status", ""),
                "acceptance_status": row.get("acceptance_status", ""),
                "interpolation_method": row.get("interpolation_method", ""),
                "extrapolation": row.get("extrapolation", ""),
                "coverage_status": row.get("coverage_status", ""),
                "acceptance_scope": row.get("acceptance_scope", ""),
                "runtime_eligible": bool(item["certified"]),
            }
        )
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
    if allow_runtime and layer != "accepted" and uncertified:
        raise PhaseReferenceContractError(
            f"runtime view rejected {uncertified} unresolved/non-certified candidate rows"
        )
    if allow_runtime and layer == "render":
        raise PhaseReferenceContractError("render layer is visualization-only and cannot be runtime input")
    if allow_runtime and layer == "accepted":
        if manifest.get("promotion_status") != "accepted_for_downstream":
            raise PhaseReferenceContractError(
                "accepted layer is not author-promoted for runtime consumption"
            )
        ineligible = [
            (table, row)
            for table, rows in table_rows.items()
            for row in rows
            if not _accepted_fallback_eligible(table, row)
        ]
        if ineligible:
            table, row = ineligible[0]
            raise PhaseReferenceContractError(
                f"accepted runtime rejected ineligible {table} row with key "
                f"{tuple(row.get(field) for field in table_specs[table]['keys'])}"
            )
        table_rows = {
            table: tuple(
                {
                    **row,
                    "runtime_eligible": True,
                    "runtime_source_layer": "accepted_primary",
                    "accepted_source_layer": row.get("source_layer", ""),
                }
                for row in rows
            )
            for table, rows in table_rows.items()
        }

    runtime_row_counts = {
        table: len(rows) for table, rows in table_rows.items()
    }
    runtime_noncertified_counts = {
        table: sum(1 for row in rows if not row.get("certified", False))
        for table, rows in table_rows.items()
    }

    diagnostics = {
        "adapter_schema": SCHEMA_VERSION,
        "reference_status": manifest.get("reference_status"),
        "runtime_consumption": bool(allow_runtime),
        "layer": layer,
        "schema_version": schema_version,
        "row_counts": row_counts,
        "uncertified_rows": uncertified,
        "allow_runtime": allow_runtime,
        "runtime_view": (
            "accepted_primary" if allow_runtime and layer == "accepted"
            else "strict_certified_only" if allow_runtime and layer == "strict"
            else "diagnostic_all_rows"
        ),
        "primary_layer": (
            "accepted" if allow_runtime and layer == "accepted"
            else "strict" if allow_runtime and layer == "strict"
            else ""
        ),
        "fallback_enabled": False,
        "fallback_order": (
            "accepted_primary" if allow_runtime and layer == "accepted"
            else "strict_primary" if allow_runtime and layer == "strict"
            else ""
        ),
        "accepted_primary_row_counts": runtime_row_counts if allow_runtime and layer == "accepted" else {},
        "accepted_primary_noncertified_row_counts": (
            runtime_noncertified_counts if allow_runtime and layer == "accepted" else {}
        ),
        "accepted_manifest_sha256": (
            sha256(manifest_path) if allow_runtime and layer == "accepted" else ""
        ),
        "accepted_layer_manifest_sha256": (
            sha256(layer_manifest_path) if allow_runtime and layer == "accepted" else ""
        ),
        "legacy_fallback_row_counts": {table: 0 for table in table_rows},
        "legacy_excluded_row_counts": {table: 0 for table in table_rows},
        "manifest_sha256": sha256(manifest_path),
        "layer_manifest_sha256": sha256(layer_manifest_path),
    }
    return PhaseReferenceBundle(root, layer, manifest, layer_manifest, table_rows, diagnostics)


def _strict_runtime_bundle(bundle: PhaseReferenceBundle) -> PhaseReferenceBundle:
    """Return the certified-only view used by an explicit strict consumer."""

    if bundle.layer != "strict":
        raise PhaseReferenceContractError("strict runtime requires the strict layer")
    tables = {
        table: tuple(
            {
                **row,
                "runtime_eligible": True,
                "runtime_source_layer": "strict_primary",
            }
            for row in rows
            if row.get("certified", False)
        )
        for table, rows in bundle.tables.items()
    }
    counts = {table: len(rows) for table, rows in tables.items()}
    diagnostics = dict(bundle.diagnostics)
    diagnostics.update(
        {
            "runtime_consumption": True,
            "runtime_view": "strict_certified_only",
            "primary_layer": "strict",
            "runtime_row_counts": counts,
            "runtime_ineligible_row_counts": {
                table: len(bundle.tables.get(table, ())) - len(rows)
                for table, rows in tables.items()
            },
            "fallback_enabled": False,
        }
    )
    return PhaseReferenceBundle(
        bundle.root,
        bundle.layer,
        bundle.manifest,
        bundle.layer_manifest,
        tables,
        diagnostics,
    )


def load_phase_reference_runtime(
    reference_root: str | Path,
    *,
    layer: str = "accepted",
) -> PhaseReferenceBundle:
    """Load the explicit runtime layer.

    ``accepted`` is the default runtime layer and admits only rows that carry
    the author-acceptance, support and phase-exclusive metadata.  ``strict``
    is an explicit opt-in certified-only view; uncertified strict rows are
    omitted rather than silently delegated to a legacy snapshot.
    """

    if layer == "accepted":
        return load_phase_reference(reference_root, layer=layer, allow_runtime=True)
    if layer == "strict":
        return _strict_runtime_bundle(load_phase_reference(reference_root, layer=layer))
    raise PhaseReferenceContractError(
        "runtime layer must be 'accepted' (default) or 'strict' (explicit)"
    )


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


def _accepted_fallback_eligible(table: str, row: Mapping[str, Any]) -> bool:
    """Return whether an accepted row may fill a runtime candidate gap.

    This deliberately does not change ``certified``.  The accepted layer is
    admitted only as an explicitly marked, non-certified fallback when its
    author-acceptance and common-support metadata are intact.
    """

    if row.get("acceptance_status") != "author_accepted_for_downstream":
        return False
    if _bool(row.get("extrapolation"), default=False):
        return False
    if row.get("coverage_status") not in {"native_support", "interpolated_common_support"}:
        return False
    if table == "crossover":
        if row.get("physical_region", "").lower() not in {"", "crossover_below_cep"}:
            return False
        endpoint = row.get("mu_CEP_proxy_MeV")
        if endpoint not in (None, "") and float(row["muq_MeV"]) > float(endpoint) + 1e-9:
            return False
    status = f"{row.get('status', '')} {row.get('source_status', '')}".lower()
    return not any(token in status for token in ("unresolved", "ambiguous", "not_converged"))


def build_accepted_primary_runtime_view(
    accepted: PhaseReferenceBundle,
) -> PhaseReferenceRuntimeView:
    """Build the accepted-primary view without any legacy fallback.

    The row-level ``certified`` field is intentionally preserved.  A row may
    therefore be runtime-eligible while remaining ``certified=False`` when it
    came from an author-reviewed interpolation.  The additional provenance
    fields make that distinction visible to downstream diagnostics.
    """

    if accepted.layer != "accepted":
        raise PhaseReferenceContractError("accepted primary view requires the accepted layer")
    if accepted.manifest.get("promotion_status") != "accepted_for_downstream":
        raise PhaseReferenceContractError(
            "accepted layer is not author-promoted for runtime consumption"
        )

    tables: dict[str, tuple[Mapping[str, Any], ...]] = {}
    row_counts: dict[str, int] = {}
    noncertified_counts: dict[str, int] = {}
    for table, rows in accepted.tables.items():
        ineligible = next(
            (row for row in rows if not _accepted_fallback_eligible(table, row)),
            None,
        )
        if ineligible is not None:
            raise PhaseReferenceContractError(
                f"accepted primary view rejected ineligible {table} row"
            )
        normalized = []
        for row in rows:
            item = dict(row)
            item.update(
                {
                    "runtime_eligible": True,
                    "runtime_source_layer": "accepted_primary",
                    "accepted_source_layer": row.get("source_layer", ""),
                }
            )
            normalized.append(item)
        tables[table] = tuple(normalized)
        row_counts[table] = len(normalized)
        noncertified_counts[table] = sum(
            1 for row in normalized if not row.get("certified", False)
        )

    diagnostics = {
        "runtime_view": "accepted_primary",
        "primary_layer": "accepted",
        "fallback_enabled": False,
        "fallback_order": "accepted_primary",
        "fallback_reason": "",
        "candidate_manifest_sha256": accepted.diagnostics.get("manifest_sha256", ""),
        "accepted_manifest_sha256": accepted.diagnostics.get("manifest_sha256", ""),
        "accepted_layer_manifest_sha256": accepted.diagnostics.get("layer_manifest_sha256", ""),
        "accepted_primary_row_counts": row_counts,
        "accepted_primary_noncertified_row_counts": noncertified_counts,
        "legacy_fallback_row_counts": {table: 0 for table in tables},
        "legacy_excluded_row_counts": {table: 0 for table in tables},
    }
    return PhaseReferenceRuntimeView(
        layer="accepted",
        source="accepted",
        tables=tables,
        diagnostics=diagnostics,
    )


def build_runtime_view(
    candidate: PhaseReferenceBundle,
    *,
    accepted_bundle: PhaseReferenceBundle | None = None,
    accepted_tables: Mapping[str, Iterable[Mapping[str, Any]]] | None = None,
    legacy_tables: Mapping[str, Iterable[Mapping[str, Any]]] | None = None,
) -> PhaseReferenceRuntimeView:
    """Build a runtime view under the current accepted/strict contract.

    ``accepted_bundle`` is the preferred path and produces an accepted-primary
    view.  Without it, the candidate is reduced to an explicit strict
    certified-only view.  Legacy rows are rejected here; historical snapshot
    comparison belongs to the dedicated retirement audit and is not a runtime
    fallback.
    """

    if accepted_bundle is not None and accepted_tables is not None:
        raise PhaseReferenceContractError("pass accepted_bundle or accepted_tables, not both")
    if legacy_tables:
        raise PhaseReferenceContractError(
            "legacy fallback is retired from runtime; use the historical retirement audit"
        )
    if accepted_bundle is not None:
        return build_accepted_primary_runtime_view(accepted_bundle)
    if accepted_tables is not None:
        raise PhaseReferenceContractError(
            "accepted_tables must be loaded as an author-promoted accepted bundle"
        )
    if candidate.layer != "strict":
        raise PhaseReferenceContractError("strict runtime view requires the strict layer")
    strict_rows = {
        table: tuple(
            {
                **row,
                "runtime_eligible": True,
                "runtime_source_layer": "strict_primary",
            }
            for row in rows
            if row.get("certified", False)
        )
        for table, rows in candidate.tables.items()
    }
    return PhaseReferenceRuntimeView(
        layer="strict",
        source="strict",
        tables=strict_rows,
        diagnostics={
            "runtime_view": "strict_certified_only",
            "primary_layer": "strict",
            "fallback_enabled": False,
            "fallback_order": "strict_primary",
            "candidate_manifest_sha256": candidate.diagnostics.get("manifest_sha256", ""),
            "strict_certified_row_counts": {
                table: len(rows) for table, rows in strict_rows.items()
            },
        },
    )


__all__ = [
    "LAYERS",
    "PhaseReferenceBundle",
    "PhaseReferenceContractError",
    "PhaseReferenceRuntimeView",
    "SCHEMA_VERSION",
    "default_downstream_layer",
    "build_accepted_primary_runtime_view",
    "build_runtime_view",
    "load_phase_reference",
    "load_phase_reference_runtime",
    "sha256",
    "to_legacy_views",
]
