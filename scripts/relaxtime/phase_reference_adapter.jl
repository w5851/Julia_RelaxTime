"""Explicit, solver-free adapter for phase-reference consumers.

The Issue #130 candidate uses versioned tables whose schema is intentionally
different from the historical CSV files.  This module is the only boundary
between those contracts.  Runtime selection is explicit and auditable:
certified candidate rows are the preferred source, author-accepted
common-support rows may be used as a marked non-certified fallback, and the
immutable legacy source is the final fallback and explicit rollback path.  This
module never deletes either source.
"""
module PhaseReferenceAdapter

using CSV
using JSON3
using SHA: sha256

export PhaseReferenceAdapterError
export PhaseReferenceSource
export load_phase_reference
export load_phase_reference_runtime
export load_phase_reference_runtime_with_fallback
export load_default_phase_reference_runtime
export load_legacy_phase_reference
export source_kind
export source_layer
export source_runtime_enabled
export source_summary
export available_xi
export boundary_data
export crossover_rows

struct PhaseReferenceAdapterError <: Exception
    message::String
end

Base.showerror(io::IO, err::PhaseReferenceAdapterError) = print(io, err.message)

"""Normalized source consumed by phase/transport reference queries."""
struct PhaseReferenceSource
    kind::Symbol
    layer::Symbol
    root::String
    runtime_enabled::Bool
    tables::Dict{Symbol, Vector{NamedTuple}}
    diagnostics::NamedTuple
end

const _TABLE_FILES = Dict{Symbol, Dict{Symbol, Union{Nothing,String}}}(
    :strict => Dict{Symbol,Union{Nothing,String}}(
        :boundary => "maxwell_surface_strict_reference_v1.csv",
        :crossover => "crossover_surface_strict_reference_v1.csv",
        :cep => "cep_boundary_strict_reference_v1.csv",
        :spinodals => "spinodal_surface_strict_reference_v1.csv",
    ),
    :derived => Dict{Symbol,Union{Nothing,String}}(
        :boundary => "maxwell_surface_derived_reference_v1.csv",
        :crossover => "crossover_surface_derived_reference_v1.csv",
        :cep => "cep_boundary_derived_reference_v1.csv",
        :spinodals => "spinodal_surface_derived_reference_v1.csv",
    ),
    :render => Dict{Symbol,Union{Nothing,String}}(
        :boundary => "maxwell_surface_render.csv",
        :crossover => "crossover_surface_render.csv",
        :cep => "cep_boundary_render.csv",
        :spinodals => nothing,
    ),
)

# v2 exposes exactly three public layers: strict, render and accepted.  The
# former v1 ``derived`` tables remain an internal build input and are not
# addressable through this map.  ``accepted`` is not a primary runtime source,
# but an explicit fallback policy may admit its author-accepted common-support
# rows while preserving their non-certified status.
const _TABLE_FILES_V2 = Dict{Symbol, Dict{Symbol, Union{Nothing,String}}}(
    :strict => Dict{Symbol,Union{Nothing,String}}(
        :boundary => "maxwell_surface_strict_reference_v1.csv",
        :crossover => "crossover_surface_strict_reference_v1.csv",
        :cep => "cep_boundary_strict_reference_v1.csv",
        :spinodals => "spinodal_surface_strict_reference_v1.csv",
    ),
    :render => Dict{Symbol,Union{Nothing,String}}(
        :boundary => "maxwell_surface_render.csv",
        :crossover => "crossover_surface_render.csv",
        :cep => "cep_boundary_render.csv",
        :spinodals => "spinodal_surface_render.csv",
    ),
    :accepted => Dict{Symbol,Union{Nothing,String}}(
        :boundary => "maxwell_surface_accepted_phase_map_v1.csv",
        :crossover => "crossover_surface_accepted_phase_map_v1.csv",
        :cep => "cep_boundary_accepted_phase_map_v1.csv",
        :spinodals => "spinodal_surface_accepted_phase_map_v1.csv",
    ),
)

const _IMPORT_SCHEMA_V1 = "pnjl_issue130_phase_reference_import_v1"
const _IMPORT_SCHEMA_V2 = "pnjl_issue130_phase_reference_v2"

const _REQUIRED_COLUMNS = Dict{Symbol, Tuple}(
    :boundary => (:xi, :T_MeV, :mu_MeV, :rho_hadron, :rho_quark, :area_residual),
    :crossover => (:xi, :mu_MeV, :T_MeV, :rho, :mu_CEP_proxy_MeV),
    :cep => (:xi, :mu_CEP_proxy_MeV, :T_low_MeV, :T_high_MeV, :T_midpoint_MeV),
    :spinodals => (:xi, :T_MeV, :mu_spinodal_hadron_MeV, :mu_spinodal_quark_MeV),
)
const _TABLE_ORDER = (:boundary, :crossover, :cep, :spinodals)

@inline function _error(message)
    throw(PhaseReferenceAdapterError(String(message)))
end

@inline function _string(value)
    value === missing && return ""
    return strip(string(value))
end

function _float(value, table::Symbol, row_number::Int, field::Symbol)
    raw = _string(value)
    parsed = tryparse(Float64, raw)
    parsed === nothing && _error("non-numeric $(field) at $(table) row $(row_number)")
    isfinite(parsed) || _error("non-finite $(field) at $(table) row $(row_number)")
    return parsed
end

function _bool(value, table::Symbol, row_number::Int, field::Symbol; default::Bool=false)
    raw = lowercase(_string(value))
    isempty(raw) && return default
    raw in ("true", "1", "yes") && return true
    raw in ("false", "0", "no") && return false
    _error("invalid boolean $(field) at $(table) row $(row_number): $(raw)")
end

function _row_value(row, field::Symbol)
    hasproperty(row, field) || return missing
    return getproperty(row, field)
end

function _read_raw_rows(path::AbstractString)
    isfile(path) || _error("missing phase-reference table: $(path)")
    rows = try
        collect(CSV.File(path; normalizenames=false))
    catch err
        _error("cannot read phase-reference table $(path): $(sprint(showerror, err))")
    end
    isempty(rows) && _error("empty phase-reference table: $(path)")
    return rows
end

function _read_rows(path::AbstractString, table::Symbol)
    rows = _read_raw_rows(path)
    first_row = first(rows)
    missing_columns = Symbol[field for field in _REQUIRED_COLUMNS[table] if !hasproperty(first_row, field)]
    isempty(missing_columns) || _error("$(path) missing columns: $(join(string.(missing_columns), ", "))")
    return rows
end

function _status_certified(table::Symbol, row, layer::Symbol, row_number::Int)
    status = lowercase(_string(_row_value(row, :status)))
    source_status = lowercase(_string(_row_value(row, :source_status)))
    interpolation = occursin("interpolat", status) ||
        occursin("interpolat", source_status) ||
        occursin("unresolved", source_status) ||
        lowercase(_string(_row_value(row, :layer))) == "interpolated_noncertified" ||
        occursin("interpolat", lowercase(_string(_row_value(row, :interpolation_method))))
    interpolation && return false
    occursin("unresolved", status) && return false
    occursin("not_converged", status) && return false
    occursin("ambiguous", status) && return false

    if table === :boundary
        _bool(_row_value(row, :grid_unresolved), table, row_number, :grid_unresolved; default=false) && return false
        hasproperty(row, :geometry_converged) && !_bool(_row_value(row, :geometry_converged), table, row_number, :geometry_converged; default=false) && return false
        hasproperty(row, :finite_and_converged) && !_bool(_row_value(row, :finite_and_converged), table, row_number, :finite_and_converged; default=false) && return false
    elseif table === :crossover
        region = lowercase(_string(_row_value(row, :physical_region)))
        !isempty(region) && region != "crossover_below_cep" && return false
    elseif table === :cep
        low = _float(_row_value(row, :T_low_MeV), table, row_number, :T_low_MeV)
        high = _float(_row_value(row, :T_high_MeV), table, row_number, :T_high_MeV)
        (high >= low && high - low <= 0.1) || return false
    end
    return true
end

function _duplicate_key(table::Symbol, row, row_number::Int)
    if table === :boundary || table === :spinodals
        return (_float(_row_value(row, :xi), table, row_number, :xi), _float(_row_value(row, :T_MeV), table, row_number, :T_MeV))
    elseif table === :crossover
        return (_float(_row_value(row, :xi), table, row_number, :xi), _float(_row_value(row, :mu_MeV), table, row_number, :mu_MeV))
    else
        return (_float(_row_value(row, :xi), table, row_number, :xi),)
    end
end

function _row_metadata(row, certified::Bool)
    return (
        source_status=_string(_row_value(row, :source_status)),
        acceptance_status=_string(_row_value(row, :acceptance_status)),
        interpolation_method=_string(_row_value(row, :interpolation_method)),
        extrapolation=_string(_row_value(row, :extrapolation)),
        coverage_status=_string(_row_value(row, :coverage_status)),
        acceptance_scope=_string(_row_value(row, :acceptance_scope)),
        runtime_eligible=certified,
    )
end

function _normalize_row(table::Symbol, row, layer::Symbol, row_number::Int)
    certified = _status_certified(table, row, layer, row_number)
    status = _string(_row_value(row, :status))
    source_layer = _string(_row_value(row, :source_layer))
    isempty(source_layer) && (source_layer = String(layer))
    if table === :boundary
        muq = _float(_row_value(row, :mu_MeV), table, row_number, :mu_MeV)
        return merge((
            xi=_float(_row_value(row, :xi), table, row_number, :xi),
            T_MeV=_float(_row_value(row, :T_MeV), table, row_number, :T_MeV),
            muq_MeV=muq,
            muB_MeV=3.0 * muq,
            rho_hadron=_float(_row_value(row, :rho_hadron), table, row_number, :rho_hadron),
            rho_quark=_float(_row_value(row, :rho_quark), table, row_number, :rho_quark),
            area_residual=_float(_row_value(row, :area_residual), table, row_number, :area_residual),
            certified=certified,
            status=status,
            source_layer=source_layer,
        ), _row_metadata(row, certified))
    elseif table === :crossover
        muq = _float(_row_value(row, :mu_MeV), table, row_number, :mu_MeV)
        return merge((
            xi=_float(_row_value(row, :xi), table, row_number, :xi),
            muq_MeV=muq,
            muB_MeV=3.0 * muq,
            T_MeV=_float(_row_value(row, :T_MeV), table, row_number, :T_MeV),
            rho=_float(_row_value(row, :rho), table, row_number, :rho),
            mu_CEP_proxy_MeV=_float(_row_value(row, :mu_CEP_proxy_MeV), table, row_number, :mu_CEP_proxy_MeV),
            physical_region=_string(_row_value(row, :physical_region)),
            certified=certified,
            status=status,
            source_layer=source_layer,
        ), _row_metadata(row, certified))
    elseif table === :cep
        muq = _float(_row_value(row, :mu_CEP_proxy_MeV), table, row_number, :mu_CEP_proxy_MeV)
        low = _float(_row_value(row, :T_low_MeV), table, row_number, :T_low_MeV)
        high = _float(_row_value(row, :T_high_MeV), table, row_number, :T_high_MeV)
        midpoint = _float(_row_value(row, :T_midpoint_MeV), table, row_number, :T_midpoint_MeV)
        return merge((
            xi=_float(_row_value(row, :xi), table, row_number, :xi),
            muq_CEP_MeV=muq,
            muB_CEP_MeV=3.0 * muq,
            T_low_MeV=low,
            T_high_MeV=high,
            T_midpoint_MeV=midpoint,
            bracket_width_MeV=high - low,
            certified=certified,
            status=status,
            boundary_mode=_string(_row_value(row, :boundary_mode)),
            source_layer=source_layer,
        ), _row_metadata(row, certified))
    else
        hadron = _float(_row_value(row, :mu_spinodal_hadron_MeV), table, row_number, :mu_spinodal_hadron_MeV)
        quark = _float(_row_value(row, :mu_spinodal_quark_MeV), table, row_number, :mu_spinodal_quark_MeV)
        return merge((
            xi=_float(_row_value(row, :xi), table, row_number, :xi),
            T_MeV=_float(_row_value(row, :T_MeV), table, row_number, :T_MeV),
            muq_spinodal_hadron_MeV=hadron,
            muq_spinodal_quark_MeV=quark,
            muB_spinodal_hadron_MeV=3.0 * hadron,
            muB_spinodal_quark_MeV=3.0 * quark,
            certified=certified,
            status=status,
            source_layer=source_layer,
        ), _row_metadata(row, certified))
    end
end

function _normalize_tables(paths::Dict{Symbol,String}, layer::Symbol; require_certified::Bool=false)
    tables = Dict{Symbol,Vector{NamedTuple}}()
    row_counts = Dict{String,Int}()
    uncertified = 0
    for (table, path) in paths
        rows = _read_rows(path, table)
        seen = Set{Any}()
        normalized = NamedTuple[]
        for (index, row) in enumerate(rows)
            row_number = index + 1
            key = _duplicate_key(table, row, row_number)
            key in seen && _error("duplicate $(table) key at row $(row_number): $(key)")
            push!(seen, key)
            item = _normalize_row(table, row, layer, row_number)
            item.certified || (uncertified += 1)
            push!(normalized, item)
        end
        tables[table] = normalized
        row_counts[String(table)] = length(normalized)
    end
    require_certified && uncertified > 0 && _error("runtime candidate rejected $(uncertified) uncertified rows")
    return tables, row_counts, uncertified
end

function _certified_runtime_source(source::PhaseReferenceSource)
    source.runtime_enabled && return source
    tables = Dict{Symbol,Vector{NamedTuple}}()
    row_counts = Dict{String,Int}()
    certified_counts = Dict{String,Int}()
    for (table, rows) in source.tables
        certified = filter(row -> row.certified, rows)
        tables[table] = [merge(row, (runtime_eligible=true,)) for row in certified]
        row_counts[String(table)] = length(rows)
        certified_counts[String(table)] = length(certified)
    end
    diagnostics = merge(source.diagnostics, (
        runtime_view="certified_only",
        certified_row_counts=certified_counts,
        fallback_enabled=false,
    ))
    return PhaseReferenceSource(source.kind, source.layer, source.root, true, tables, diagnostics)
end

function _table_key(table::Symbol, row)
    if table === :boundary || table === :spinodals
        return (row.xi, row.T_MeV)
    elseif table === :crossover
        return (row.xi, row.muq_MeV)
    else
        return (row.xi,)
    end
end

function _accepted_fallback_eligible(table::Symbol, row)
    hasproperty(row, :acceptance_status) || return false
    lowercase(_string(row.acceptance_status)) == "author_accepted_for_downstream" || return false
    extrapolation = hasproperty(row, :extrapolation) ?
        _bool(row.extrapolation, table, 0, :extrapolation; default=false) : false
    extrapolation && return false
    coverage = hasproperty(row, :coverage_status) ? lowercase(_string(row.coverage_status)) : ""
    coverage in ("native_support", "interpolated_common_support") || return false
    if table === :crossover && hasproperty(row, :physical_region)
        region = lowercase(_string(row.physical_region))
        region in ("", "crossover_below_cep") || return false
        if hasproperty(row, :mu_CEP_proxy_MeV) && isfinite(row.mu_CEP_proxy_MeV) &&
            row.muq_MeV > row.mu_CEP_proxy_MeV + 1e-9
            return false
        end
    end
    status = lowercase(join((_string(row.status), _string(row.source_status)), " "))
    any(token -> occursin(token, status), ("unresolved", "ambiguous", "not_converged")) && return false
    return true
end

function _merge_candidate_with_accepted_and_legacy(
    candidate::PhaseReferenceSource,
    accepted::Union{Nothing,PhaseReferenceSource},
    legacy::PhaseReferenceSource,
)
    candidate_runtime = _certified_runtime_source(candidate)
    tables = Dict{Symbol,Vector{NamedTuple}}()
    candidate_counts = Dict{String,Int}()
    accepted_counts = Dict{String,Int}()
    legacy_counts = Dict{String,Int}()
    legacy_excluded_counts = Dict{String,Int}()
    cep_limits = NamedTuple[]
    for row in get(candidate_runtime.tables, :cep, NamedTuple[])
        _runtime_eligible(row) && push!(cep_limits, row)
    end
    if accepted !== nothing
        for row in get(accepted.tables, :cep, NamedTuple[])
            _accepted_fallback_eligible(:cep, row) && push!(cep_limits, row)
        end
    end
    function cep_limit(xi)
        exact = findfirst(row -> isapprox(row.xi, xi; atol=1e-6, rtol=0.0), cep_limits)
        exact === nothing ? nothing : cep_limits[exact].muq_CEP_MeV
    end
    for table in _TABLE_ORDER
        rows = NamedTuple[]
        seen = Set{Any}()
        for row in get(candidate_runtime.tables, table, NamedTuple[])
            push!(rows, row)
            push!(seen, _table_key(table, row))
        end
        candidate_counts[String(table)] = length(rows)
        n_accepted = 0
        if accepted !== nothing
            for row in get(accepted.tables, table, NamedTuple[])
                _accepted_fallback_eligible(table, row) || continue
                key = _table_key(table, row)
                key in seen && continue
                accepted_row = merge(row, (
                    source_layer="accepted_fallback",
                    status="accepted_fallback",
                    accepted_source_status=hasproperty(row, :source_status) ? row.source_status : "",
                    runtime_eligible=true,
                    certified=false,
                ))
                push!(rows, accepted_row)
                push!(seen, key)
                n_accepted += 1
            end
        end
        accepted_counts[String(table)] = n_accepted
        n_legacy = 0
        for row in get(legacy.tables, table, NamedTuple[])
            key = _table_key(table, row)
            key in seen && continue
            fallback_row = merge(row, (
                source_layer="legacy_fallback",
                status="legacy_fallback",
                runtime_eligible=true,
                certified=true,
            ))
            push!(rows, fallback_row)
            push!(seen, key)
            n_legacy += 1
        end
        n_legacy_excluded = 0
        # The historical crossover file contains derivative peaks on both
        # sides of the endpoint.  Only the below-CEP part has crossover
        # semantics in a runtime phase view.
        if table === :crossover
            filtered_rows = NamedTuple[]
            for row in rows
                if row.source_layer == "legacy_fallback"
                    limit = cep_limit(row.xi)
                    if limit !== nothing && row.muq_MeV > limit + 1e-9
                        n_legacy_excluded += 1
                        continue
                    end
                end
                push!(filtered_rows, row)
            end
            rows = filtered_rows
        end
        # Count only rows that remain available to this runtime view.  Rows
        # above a known CEP are retained in the byte-preserving legacy
        # snapshot, but are not runtime crossover fallback rows.
        n_legacy -= n_legacy_excluded
        legacy_counts[String(table)] = n_legacy
        legacy_excluded_counts[String(table)] = n_legacy_excluded
        tables[table] = rows
    end
    accepted_enabled = accepted !== nothing
    diagnostics = merge(candidate_runtime.diagnostics, (
        runtime_view=accepted_enabled ?
            "certified_candidate_with_accepted_then_legacy_fallback" :
            "certified_candidate_with_legacy_fallback",
        fallback_enabled=true,
        fallback_order=accepted_enabled ?
            "strict_candidate>accepted_downstream>legacy_snapshot" :
            "strict_candidate>legacy_snapshot",
        fallback_reason="candidate_key_absent_or_uncertified",
        candidate_row_counts=candidate_counts,
        accepted_fallback_row_counts=accepted_counts,
        legacy_fallback_row_counts=legacy_counts,
        legacy_excluded_row_counts=legacy_excluded_counts,
        # Preserve the historical diagnostic key for existing consumers.
        fallback_row_counts=legacy_counts,
        legacy_row_counts=legacy.diagnostics.row_counts,
        accepted_manifest_sha256=accepted === nothing ? "" :
            get(accepted.diagnostics, :candidate_manifest_sha256, ""),
        accepted_layer_manifest_sha256=accepted === nothing ? "" :
            get(accepted.diagnostics, :candidate_layer_manifest_sha256, ""),
        legacy_manifest_reference_status="legacy",
    ))
    return PhaseReferenceSource(:candidate, candidate.layer, candidate.root, true, tables, diagnostics)
end

function _merge_certified_candidate_with_legacy(candidate::PhaseReferenceSource, legacy::PhaseReferenceSource)
    return _merge_candidate_with_accepted_and_legacy(candidate, nothing, legacy)
end

function _manifest(path::AbstractString)
    isfile(path) || _error("missing phase-reference manifest: $(path)")
    try
        return JSON3.read(Base.read(path, String), Dict{String,Any})
    catch err
        _error("invalid phase-reference manifest $(path): $(sprint(showerror, err))")
    end
end

function _file_sha256(path::AbstractString)
    isfile(path) || return ""
    return bytes2hex(sha256(read(path)))
end

function load_phase_reference(root::AbstractString; layer::Symbol=:strict, allow_runtime::Bool=false)
    root_abs = normpath(abspath(root))
    manifest = _manifest(joinpath(root_abs, "manifest.json"))
    schema_version = String(get(manifest, "schema_version", ""))
    table_files = if schema_version == _IMPORT_SCHEMA_V1
        layer in (:strict, :derived, :render) ||
            _error("phase-reference layer must be strict, derived, or render for v1")
        _TABLE_FILES
    elseif schema_version == _IMPORT_SCHEMA_V2
        layer in (:strict, :render, :accepted) ||
            _error("phase-reference layer must be strict, render, or accepted for v2")
        _TABLE_FILES_V2
    else
        _error("candidate import schema mismatch")
    end
    get(manifest, "runtime_consumption", true) === false || _error("candidate runtime_consumption must remain false")
    get(manifest, "reference_status", "") in ("candidate", "imported_candidate") ||
        _error("candidate reference_status is not a candidate state")
    layer_root = joinpath(root_abs, String(layer))
    layer_manifest = _manifest(joinpath(layer_root, "manifest.json"))
    get(layer_manifest, "layer", "") isa String || _error("candidate layer manifest has no layer identifier")
    if allow_runtime && layer in (:render, :accepted)
        layer === :render && _error("render layer is visualization-only and cannot be runtime input")
        _error("accepted layer cannot be the primary runtime source; pass it through the accepted fallback policy")
    end
    paths = Dict{Symbol,String}()
    for (table, filename) in table_files[layer]
        filename === nothing && continue
        paths[table] = joinpath(layer_root, "tables", filename)
    end
    tables, counts, uncertified = _normalize_tables(paths, layer; require_certified=allow_runtime)
    return PhaseReferenceSource(
        :candidate,
        layer,
        root_abs,
        allow_runtime,
        tables,
        (schema_version="pnjl_issue130_phase_reference_adapter_v1", candidate_schema_version=schema_version,
         row_counts=counts,
         uncertified_rows=uncertified,
         manifest_reference_status=String(get(manifest, "reference_status", "")),
         promotion_status=String(get(manifest, "promotion_status", "")),
         downstream_default_layer=String(get(manifest, "downstream_default_layer", "")),
         candidate_manifest_sha256=_file_sha256(joinpath(root_abs, "manifest.json")),
         candidate_layer_manifest_sha256=_file_sha256(joinpath(layer_root, "manifest.json")),
         source_root=root_abs,
         runtime_view=allow_runtime ? "strict_all_rows" : "diagnostic_all_rows",
         fallback_enabled=false),
    )
end

function load_phase_reference_runtime(root::AbstractString; layer::Symbol=:strict)
    return load_phase_reference(root; layer=layer, allow_runtime=true)
end

function load_phase_reference_runtime_with_fallback(candidate_root::AbstractString;
    layer::Symbol=:strict,
    boundary_path::AbstractString,
    cep_path::AbstractString,
    crossover_path::AbstractString,
    spinodals_path::AbstractString,
    accepted_root::Union{Nothing,AbstractString}=nothing,
)
    layer in (:render, :accepted) &&
        _error("$(layer) layer is downstream-only and cannot be used for runtime fallback")
    legacy = load_legacy_phase_reference(
        boundary_path=boundary_path,
        cep_path=cep_path,
        crossover_path=crossover_path,
        spinodals_path=spinodals_path,
    )
    accepted = nothing
    accepted_error = ""
    if accepted_root !== nothing && !isempty(String(accepted_root))
        accepted = try
            source = load_phase_reference(String(accepted_root); layer=:accepted)
            get(source.diagnostics, :promotion_status, "") == "accepted_for_downstream" ||
                _error("accepted package is not author-promoted for downstream fallback")
            source
        catch err
            accepted_error = sprint(showerror, err)
            nothing
        end
    end
    candidate = try
        load_phase_reference(candidate_root; layer=layer)
    catch err
        if accepted !== nothing
            empty_candidate = PhaseReferenceSource(
                :candidate,
                layer,
                String(candidate_root),
                false,
                Dict{Symbol,Vector{NamedTuple}}(),
                (
                    schema_version="candidate_load_failed",
                    candidate_schema_version="",
                    row_counts=Dict{String,Int}(),
                    uncertified_rows=0,
                    manifest_reference_status="candidate",
                    promotion_status="",
                    downstream_default_layer="",
                    candidate_manifest_sha256="",
                    candidate_layer_manifest_sha256="",
                    source_root=String(candidate_root),
                    runtime_view="candidate_load_failed",
                    fallback_enabled=true,
                    candidate_load_error=sprint(showerror, err),
                ),
            )
            merged = _merge_candidate_with_accepted_and_legacy(empty_candidate, accepted, legacy)
            return PhaseReferenceSource(
                merged.kind,
                merged.layer,
                merged.root,
                merged.runtime_enabled,
                merged.tables,
                merge(merged.diagnostics, (accepted_load_error=accepted_error,)),
            )
        end
        diagnostics = merge(legacy.diagnostics, (
            runtime_view="legacy_fallback",
            fallback_enabled=true,
            fallback_reason="candidate_load_failed: " * sprint(showerror, err),
            accepted_load_error=accepted_error,
        ))
        return PhaseReferenceSource(:legacy, :legacy, "", true, legacy.tables, diagnostics)
    end
    merged = accepted === nothing ?
        _merge_certified_candidate_with_legacy(candidate, legacy) :
        _merge_candidate_with_accepted_and_legacy(candidate, accepted, legacy)
    accepted_error == "" && return merged
    return PhaseReferenceSource(merged.kind, merged.layer, merged.root, merged.runtime_enabled, merged.tables,
        merge(merged.diagnostics, (accepted_load_error=accepted_error,)))
end

function load_default_phase_reference_runtime(; project_root::AbstractString,
    layer::Symbol=:strict,
    source::Symbol=:candidate,
)
    source in (:candidate, :legacy) || _error("phase-reference source must be candidate or legacy")
    reference_root = joinpath(project_root, "data", "reference", "pnjl")
    legacy_root = joinpath(reference_root, "legacy_phase_reference_v1")
    legacy_paths = (
        boundary_path=joinpath(legacy_root, "boundary.csv"),
        cep_path=joinpath(legacy_root, "cep.csv"),
        crossover_path=joinpath(legacy_root, "crossover_dense.csv"),
        spinodals_path=joinpath(legacy_root, "spinodals.csv"),
    )
    source === :legacy && return load_legacy_phase_reference(; legacy_paths...)
    candidate_root = joinpath(reference_root, "issue130_phase_reference_v1")
    accepted_root = joinpath(reference_root, "issue130_phase_reference_v2")
    return load_phase_reference_runtime_with_fallback(
        candidate_root;
        layer=layer,
        accepted_root=accepted_root,
        legacy_paths...,
    )
end

function _legacy_rows(path::AbstractString, table::Symbol)
    # The legacy files have a stable schema but may contain a small number of
    # optional diagnostic columns.  Reuse the same normalized query contract.
    return _read_raw_rows(path)
end

function load_legacy_phase_reference(; boundary_path::AbstractString, cep_path::AbstractString,
    crossover_path::AbstractString, spinodals_path::AbstractString)
    paths = Dict{Symbol,String}(
        :boundary => String(boundary_path), :cep => String(cep_path),
        :crossover => String(crossover_path), :spinodals => String(spinodals_path),
    )
    # Convert legacy headers to the candidate-shaped rows before normalization.
    tables = Dict{Symbol,Vector{NamedTuple}}()
    row_counts = Dict{String,Int}()
    for (table, path) in paths
        raw_rows = _legacy_rows(path, table)
        normalized = NamedTuple[]
        for (index, row) in enumerate(raw_rows)
            row_number = index + 1
            item = if table === :boundary
                (
                    xi=_float(_row_value(row, :xi), table, row_number, :xi),
                    T_MeV=_float(_row_value(row, :T_MeV), table, row_number, :T_MeV),
                    muq_MeV=_float(_row_value(row, :mu_transition_MeV), table, row_number, :mu_transition_MeV),
                    muB_MeV=3.0 * _float(_row_value(row, :mu_transition_MeV), table, row_number, :mu_transition_MeV),
                    rho_hadron=_float(_row_value(row, :rho_hadron), table, row_number, :rho_hadron),
                    rho_quark=_float(_row_value(row, :rho_quark), table, row_number, :rho_quark),
                    area_residual=hasproperty(row, :area_residual) ? _float(_row_value(row, :area_residual), table, row_number, :area_residual) : NaN,
                    certified=true, status="legacy", source_layer="legacy",
                )
            elseif table === :crossover
                muq = _float(_row_value(row, :mu_MeV), table, row_number, :mu_MeV)
                (
                    xi=_float(_row_value(row, :xi), table, row_number, :xi),
                    muq_MeV=muq, muB_MeV=3.0 * muq,
                    T_MeV=_float(_row_value(row, :T_crossover_MeV), table, row_number, :T_crossover_MeV),
                    rho=hasproperty(row, :rho) ? _float(_row_value(row, :rho), table, row_number, :rho) : NaN,
                    mu_CEP_proxy_MeV=NaN, physical_region="crossover_below_cep",
                    certified=true, status="legacy", source_layer="legacy",
                )
            elseif table === :cep
                muq = if hasproperty(row, :muq_CEP_MeV)
                    _float(_row_value(row, :muq_CEP_MeV), table, row_number, :muq_CEP_MeV)
                else
                    _float(_row_value(row, :mu_CEP_MeV), table, row_number, :mu_CEP_MeV)
                end
                t = _float(_row_value(row, :T_CEP_MeV), table, row_number, :T_CEP_MeV)
                (
                    xi=_float(_row_value(row, :xi), table, row_number, :xi),
                    muq_CEP_MeV=muq, muB_CEP_MeV=3.0 * muq,
                    T_low_MeV=t, T_high_MeV=t, T_midpoint_MeV=t,
                    bracket_width_MeV=0.0, certified=true, status="legacy",
                    boundary_mode="estimated_midpoint", source_layer="legacy",
                )
            else
                hadron = _float(_row_value(row, :mu_spinodal_hadron_MeV), table, row_number, :mu_spinodal_hadron_MeV)
                quark = _float(_row_value(row, :mu_spinodal_quark_MeV), table, row_number, :mu_spinodal_quark_MeV)
                (
                    xi=_float(_row_value(row, :xi), table, row_number, :xi),
                    T_MeV=_float(_row_value(row, :T_MeV), table, row_number, :T_MeV),
                    muq_spinodal_hadron_MeV=hadron, muq_spinodal_quark_MeV=quark,
                    muB_spinodal_hadron_MeV=3.0 * hadron, muB_spinodal_quark_MeV=3.0 * quark,
                    certified=true, status="legacy", source_layer="legacy",
                )
            end
            push!(normalized, item)
        end
        tables[table] = normalized
        row_counts[String(table)] = length(normalized)
    end
    return PhaseReferenceSource(:legacy, :legacy, "", true, tables,
        (schema_version="legacy", row_counts=row_counts, uncertified_rows=0,
         manifest_reference_status="legacy", promotion_status="", downstream_default_layer="",
         candidate_manifest_sha256="", candidate_layer_manifest_sha256="", source_root="",
         runtime_view="legacy", fallback_enabled=false))
end

source_kind(source::PhaseReferenceSource) = source.kind
source_layer(source::PhaseReferenceSource) = source.layer
source_runtime_enabled(source::PhaseReferenceSource) = source.runtime_enabled
source_summary(source::PhaseReferenceSource) = source.diagnostics

@inline function _runtime_eligible(row)
    hasproperty(row, :runtime_eligible) && Bool(row.runtime_eligible) && return true
    return hasproperty(row, :certified) && Bool(row.certified)
end

function available_xi(source::PhaseReferenceSource, table::Symbol)
    haskey(source.tables, table) || return Float64[]
    return unique(sort(Float64[row.xi for row in source.tables[table]]))
end

function _nearest_xi(source::PhaseReferenceSource, table::Symbol, xi::Float64)
    values = available_xi(source, table)
    isempty(values) && return nothing
    return values[argmin(abs.(values .- xi))]
end

function _rows_at_xi(source::PhaseReferenceSource, table::Symbol, xi::Float64)
    haskey(source.tables, table) || return NamedTuple[]
    exact = filter(row -> isapprox(row.xi, xi; atol=1e-6, rtol=0.0), source.tables[table])
    isempty(exact) || return exact
    nearest = _nearest_xi(source, table, xi)
    nearest === nothing && return NamedTuple[]
    return filter(row -> isapprox(row.xi, nearest; atol=1e-6, rtol=0.0), source.tables[table])
end

function boundary_data(source::PhaseReferenceSource, xi::Float64; require_certified::Bool=true)
    require_certified && !source.runtime_enabled && _error("phase-reference source is not runtime-enabled; pass require_certified=false for an explicit diagnostic view")
    rows = sort(_rows_at_xi(source, :boundary, xi), by=row -> row.T_MeV)
    require_certified && any(!_runtime_eligible(row) for row in rows) && _error("phase-reference boundary contains runtime-ineligible rows")
    cep_rows = _rows_at_xi(source, :cep, xi)
    cep = isempty(cep_rows) ? nothing : first(sort(cep_rows, by=row -> row.xi))
    require_certified && cep !== nothing && !_runtime_eligible(cep) && _error("phase-reference CEP contains a runtime-ineligible row")
    return (
        T_values=Float64[row.T_MeV for row in rows],
        mu_values=Float64[row.muq_MeV for row in rows],
        T_CEP=cep === nothing ? NaN : cep.T_midpoint_MeV,
        mu_CEP=cep === nothing ? NaN : cep.muq_CEP_MeV,
        muq_CEP=cep === nothing ? NaN : cep.muq_CEP_MeV,
        muB_CEP=cep === nothing ? NaN : cep.muB_CEP_MeV,
        xi=isempty(rows) ? Float64(xi) : first(rows).xi,
    )
end

function crossover_rows(source::PhaseReferenceSource, xi::Float64; require_certified::Bool=true)
    require_certified && !source.runtime_enabled && _error("phase-reference source is not runtime-enabled; pass require_certified=false for an explicit diagnostic view")
    rows = _rows_at_xi(source, :crossover, xi)
    require_certified && any(!_runtime_eligible(row) for row in rows) && _error("phase-reference crossover contains runtime-ineligible rows")
    return sort(rows, by=row -> row.muq_MeV)
end

end # module PhaseReferenceAdapter
