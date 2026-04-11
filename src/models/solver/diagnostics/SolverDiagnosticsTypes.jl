"""SolverDiagnosticsTypes

诊断数据内部类型化表示，以及对旧 NamedTuple 视图的兼容映射。
"""

struct SolverDiagnosticSummary
    diagnostic_version::Symbol
    attempt_origin::Symbol
    seed_source::Union{Symbol,Nothing}
    hard_constraint_ok::Union{Bool,Nothing}
    failed_constraints::Vector{Symbol}
    error_kind::Symbol
    error_msg::String
    selection_reason::Symbol
    endpoint_cause::Union{Symbol,Nothing}
    continuity_distance::Union{Float64,Nothing}
    selection_reason_source::Symbol
end

struct SolverDiagnosticCandidate
    summary::SolverDiagnosticSummary
end

struct SolverDiagnosticFull
    summary::SolverDiagnosticSummary
    candidates::Vector{SolverDiagnosticCandidate}
end

const SOLVER_DIAGNOSTIC_VERSION_V1 = :v1
const SOLVER_DIAGNOSTIC_PUBLIC_FIELDS = (
    :diagnostic_version,
    :attempt_origin,
    :seed_source,
    :hard_constraint_ok,
    :failed_constraints,
    :error_kind,
    :error_msg,
    :selection_reason,
    :endpoint_cause,
    :continuity_distance,
)
const SOLVER_DIAGNOSTIC_INTERNAL_FIELDS = (
    :selection_reason_source,
)

@inline solver_diagnostic_public_fields() = SOLVER_DIAGNOSTIC_PUBLIC_FIELDS
@inline solver_diagnostic_internal_fields() = SOLVER_DIAGNOSTIC_INTERNAL_FIELDS

@inline function _normalize_solver_diagnostic_version(version)::Symbol
    if version isa Symbol
        version == SOLVER_DIAGNOSTIC_VERSION_V1 || throw(ArgumentError("unsupported solver diagnostic_version: $(version), expected $(SOLVER_DIAGNOSTIC_VERSION_V1)"))
        return version
    end
    version isa AbstractString || throw(ArgumentError("diagnostic_version must be Symbol or AbstractString, got $(typeof(version))"))
    normalized = Symbol(version)
    normalized == SOLVER_DIAGNOSTIC_VERSION_V1 || throw(ArgumentError("unsupported solver diagnostic_version: $(normalized), expected $(SOLVER_DIAGNOSTIC_VERSION_V1)"))
    return normalized
end

function SolverDiagnosticSummary(
    attempt_origin::Symbol,
    seed_source::Union{Symbol,Nothing},
    hard_constraint_ok::Union{Bool,Nothing},
    failed_constraints::Vector{Symbol},
    error_kind::Symbol,
    error_msg::String,
    selection_reason::Symbol,
    endpoint_cause::Union{Symbol,Nothing},
    continuity_distance::Union{Float64,Nothing},
    selection_reason_source::Symbol,
)
    return SolverDiagnosticSummary(
        SOLVER_DIAGNOSTIC_VERSION_V1,
        attempt_origin,
        seed_source,
        hard_constraint_ok,
        failed_constraints,
        error_kind,
        error_msg,
        selection_reason,
        endpoint_cause,
        continuity_distance,
        selection_reason_source,
    )
end

function SolverDiagnosticSummary(
    diagnostic_version,
    attempt_origin::Symbol,
    seed_source::Union{Symbol,Nothing},
    hard_constraint_ok::Union{Bool,Nothing},
    failed_constraints::Vector{Symbol},
    error_kind::Symbol,
    error_msg::String,
    selection_reason::Symbol,
    endpoint_cause::Union{Symbol,Nothing},
    continuity_distance::Union{Float64,Nothing},
    selection_reason_source::Symbol,
)
    return SolverDiagnosticSummary(
        _normalize_solver_diagnostic_version(diagnostic_version),
        attempt_origin,
        seed_source,
        hard_constraint_ok,
        failed_constraints,
        error_kind,
        error_msg,
        selection_reason,
        endpoint_cause,
        continuity_distance,
        selection_reason_source,
    )
end

@inline solver_diagnostic_version(summary::SolverDiagnosticSummary) = summary.diagnostic_version

@inline function to_namedtuple(summary::SolverDiagnosticSummary)
    return (
        diagnostic_version=summary.diagnostic_version,
        attempt_origin=summary.attempt_origin,
        seed_source=summary.seed_source,
        hard_constraint_ok=summary.hard_constraint_ok,
        failed_constraints=summary.failed_constraints,
        error_kind=summary.error_kind,
        error_msg=summary.error_msg,
        selection_reason=summary.selection_reason,
        endpoint_cause=summary.endpoint_cause,
        continuity_distance=summary.continuity_distance,
        selection_reason_source=summary.selection_reason_source,
    )
end

@inline function to_namedtuple(candidate::SolverDiagnosticCandidate)
    return to_namedtuple(candidate.summary)
end

@inline function to_public_namedtuple(summary::SolverDiagnosticSummary)
    base = to_namedtuple(summary)
    return NamedTuple{SOLVER_DIAGNOSTIC_PUBLIC_FIELDS}(base)
end

@inline function to_public_namedtuple(candidate::SolverDiagnosticCandidate)
    return to_public_namedtuple(candidate.summary)
end

@inline function to_public_namedtuple(full::SolverDiagnosticFull)
    summary_public = to_public_namedtuple(full.summary)
    candidates_public = [to_public_namedtuple(c) for c in full.candidates]
    return (; summary_public..., candidates=candidates_public)
end

@inline function to_namedtuple(full::SolverDiagnosticFull)
    summary = to_namedtuple(full.summary)
    candidates = [to_namedtuple(c) for c in full.candidates]
    return (; summary..., candidates=candidates)
end

@inline function _summary_from_candidate(
    candidate;
    seed_source::Union{Symbol,Nothing}=nothing,
    selection_reason::Union{Nothing,Symbol}=nothing,
    selection_reason_source::Symbol=:problem_spec_selector,
)
    hard_ok = hasproperty(candidate, :hard_constraint_ok) ? Bool(getproperty(candidate, :hard_constraint_ok)) : nothing
    failed = hasproperty(candidate, :failed_constraints) ? Symbol.(getproperty(candidate, :failed_constraints)) : Symbol[]
    continuity = if hasproperty(candidate, :continuity_distance)
        value = getproperty(candidate, :continuity_distance)
        value === nothing ? nothing : Float64(value)
    else
        nothing
    end
    error_kind = hasproperty(candidate, :error_kind) ? Symbol(getproperty(candidate, :error_kind)) : :none
    error_msg = hasproperty(candidate, :error_msg) ? String(getproperty(candidate, :error_msg)) : ""
    reason = selection_reason === nothing ?
             (hasproperty(candidate, :selection_reason) ? Symbol(getproperty(candidate, :selection_reason)) : :none) :
             selection_reason
    endpoint = _diagnostic_endpoint_cause(candidate)
    return SolverDiagnosticSummary(
        _diagnostic_attempt_origin(candidate),
        seed_source,
        hard_ok,
        failed,
        error_kind,
        error_msg,
        reason,
        endpoint,
        continuity,
        selection_reason_source,
    )
end

@inline function coerce_solver_diagnostic_summary(diag)
    if diag isa SolverDiagnosticSummary
        return diag
    end

    diag_nt = diag isa NamedTuple ? diag : (applicable(to_namedtuple, diag) ? to_namedtuple(diag) : nothing)
    diag_nt === nothing && throw(ArgumentError("diag must be SolverDiagnosticSummary-compatible, got $(typeof(diag))"))

    return SolverDiagnosticSummary(
        get(diag_nt, :diagnostic_version, SOLVER_DIAGNOSTIC_VERSION_V1),
        get(diag_nt, :attempt_origin, :fallback) |> Symbol,
        get(diag_nt, :seed_source, nothing),
        get(diag_nt, :hard_constraint_ok, nothing),
        Symbol.(get(diag_nt, :failed_constraints, Symbol[])),
        get(diag_nt, :error_kind, :none) |> Symbol,
        String(get(diag_nt, :error_msg, "")),
        get(diag_nt, :selection_reason, :none) |> Symbol,
        begin
            endpoint = get(diag_nt, :endpoint_cause, nothing)
            endpoint === nothing ? nothing : Symbol(endpoint)
        end,
        begin
            continuity = get(diag_nt, :continuity_distance, nothing)
            continuity === nothing ? nothing : Float64(continuity)
        end,
        get(diag_nt, :selection_reason_source, :problem_spec_selector) |> Symbol,
    )
end

@inline coerce_solver_diagnostic_public_view(diag) = to_public_namedtuple(coerce_solver_diagnostic_summary(diag))
