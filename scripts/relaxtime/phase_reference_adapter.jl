"""Explicit, solver-free adapter for phase-reference consumers.

The Issue #130 candidate uses versioned tables whose schema is intentionally
different from the historical CSV files.  This module is the only boundary
between those contracts.  Runtime selection is explicit and auditable: the
candidate strict layer is the preferred source, only certified candidate rows
are exposed, and missing keys may fall back to the legacy source.  The legacy
source remains an explicit rollback path and is never deleted by this module.
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
    interpolation = occursin("interpolat", status) ||
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

function _normalize_row(table::Symbol, row, layer::Symbol, row_number::Int)
    certified = _status_certified(table, row, layer, row_number)
    status = _string(_row_value(row, :status))
    source_layer = _string(_row_value(row, :source_layer))
    isempty(source_layer) && (source_layer = String(layer))
    if table === :boundary
        muq = _float(_row_value(row, :mu_MeV), table, row_number, :mu_MeV)
        return (
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
        )
    elseif table === :crossover
        muq = _float(_row_value(row, :mu_MeV), table, row_number, :mu_MeV)
        return (
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
        )
    elseif table === :cep
        muq = _float(_row_value(row, :mu_CEP_proxy_MeV), table, row_number, :mu_CEP_proxy_MeV)
        low = _float(_row_value(row, :T_low_MeV), table, row_number, :T_low_MeV)
        high = _float(_row_value(row, :T_high_MeV), table, row_number, :T_high_MeV)
        midpoint = _float(_row_value(row, :T_midpoint_MeV), table, row_number, :T_midpoint_MeV)
        return (
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
        )
    else
        hadron = _float(_row_value(row, :mu_spinodal_hadron_MeV), table, row_number, :mu_spinodal_hadron_MeV)
        quark = _float(_row_value(row, :mu_spinodal_quark_MeV), table, row_number, :mu_spinodal_quark_MeV)
        return (
            xi=_float(_row_value(row, :xi), table, row_number, :xi),
            T_MeV=_float(_row_value(row, :T_MeV), table, row_number, :T_MeV),
            muq_spinodal_hadron_MeV=hadron,
            muq_spinodal_quark_MeV=quark,
            muB_spinodal_hadron_MeV=3.0 * hadron,
            muB_spinodal_quark_MeV=3.0 * quark,
            certified=certified,
            status=status,
            source_layer=source_layer,
        )
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
        tables[table] = certified
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

function _merge_certified_candidate_with_legacy(candidate::PhaseReferenceSource, legacy::PhaseReferenceSource)
    candidate_runtime = _certified_runtime_source(candidate)
    tables = Dict{Symbol,Vector{NamedTuple}}()
    candidate_counts = Dict{String,Int}()
    fallback_counts = Dict{String,Int}()
    for table in _TABLE_ORDER
        rows = NamedTuple[]
        seen = Set{Any}()
        for row in get(candidate_runtime.tables, table, NamedTuple[])
            push!(rows, row)
            push!(seen, _table_key(table, row))
        end
        candidate_counts[String(table)] = length(rows)
        n_fallback = 0
        for row in get(legacy.tables, table, NamedTuple[])
            key = _table_key(table, row)
            key in seen && continue
            fallback_row = merge(row, (source_layer="legacy_fallback", status="legacy_fallback", certified=true))
            push!(rows, fallback_row)
            push!(seen, key)
            n_fallback += 1
        end
        fallback_counts[String(table)] = n_fallback
        tables[table] = rows
    end
    diagnostics = merge(candidate_runtime.diagnostics, (
        runtime_view="certified_candidate_with_legacy_fallback",
        fallback_enabled=true,
        fallback_reason="candidate_key_absent_or_uncertified",
        candidate_row_counts=candidate_counts,
        fallback_row_counts=fallback_counts,
        legacy_row_counts=legacy.diagnostics.row_counts,
        legacy_manifest_reference_status="legacy",
    ))
    return PhaseReferenceSource(:candidate, candidate.layer, candidate.root, true, tables, diagnostics)
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
    layer in (:strict, :derived, :render) || _error("phase-reference layer must be strict, derived, or render")
    root_abs = normpath(abspath(root))
    manifest = _manifest(joinpath(root_abs, "manifest.json"))
    get(manifest, "schema_version", "") == "pnjl_issue130_phase_reference_import_v1" ||
        _error("candidate import schema mismatch")
    get(manifest, "runtime_consumption", true) === false || _error("candidate runtime_consumption must remain false")
    get(manifest, "reference_status", "") in ("candidate", "imported_candidate") ||
        _error("candidate reference_status is not a candidate state")
    layer_root = joinpath(root_abs, String(layer))
    layer_manifest = _manifest(joinpath(layer_root, "manifest.json"))
    get(layer_manifest, "layer", "") isa String || _error("candidate layer manifest has no layer identifier")
    if layer === :render && allow_runtime
        _error("render layer is visualization-only and cannot be runtime input")
    end
    paths = Dict{Symbol,String}()
    for (table, filename) in _TABLE_FILES[layer]
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
        (schema_version="pnjl_issue130_phase_reference_adapter_v1", row_counts=counts,
         uncertified_rows=uncertified,
         manifest_reference_status=String(get(manifest, "reference_status", "")),
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
)
    legacy = load_legacy_phase_reference(
        boundary_path=boundary_path,
        cep_path=cep_path,
        crossover_path=crossover_path,
        spinodals_path=spinodals_path,
    )
    candidate = try
        load_phase_reference(candidate_root; layer=layer)
    catch err
        diagnostics = merge(legacy.diagnostics, (
            runtime_view="legacy_fallback",
            fallback_enabled=true,
            fallback_reason="candidate_load_failed: " * sprint(showerror, err),
        ))
        return PhaseReferenceSource(:legacy, :legacy, "", true, legacy.tables, diagnostics)
    end
    return _merge_certified_candidate_with_legacy(candidate, legacy)
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
    return load_phase_reference_runtime_with_fallback(candidate_root; layer=layer, legacy_paths...)
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
         manifest_reference_status="legacy", candidate_manifest_sha256="",
         candidate_layer_manifest_sha256="", source_root="", runtime_view="legacy", fallback_enabled=false))
end

source_kind(source::PhaseReferenceSource) = source.kind
source_layer(source::PhaseReferenceSource) = source.layer
source_runtime_enabled(source::PhaseReferenceSource) = source.runtime_enabled
source_summary(source::PhaseReferenceSource) = source.diagnostics

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
    require_certified && any(!row.certified for row in rows) && _error("phase-reference boundary contains uncertified rows")
    cep_rows = _rows_at_xi(source, :cep, xi)
    cep = isempty(cep_rows) ? nothing : first(sort(cep_rows, by=row -> row.xi))
    require_certified && cep !== nothing && !cep.certified && _error("phase-reference CEP contains an uncertified row")
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
    require_certified && any(!row.certified for row in rows) && _error("phase-reference crossover contains uncertified rows")
    return sort(rows, by=row -> row.muq_MeV)
end

end # module PhaseReferenceAdapter
