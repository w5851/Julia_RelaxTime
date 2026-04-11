using Dates: DateTime

abstract type AbstractPipelineIOContract end

struct PipelineIOContract <: AbstractPipelineIOContract
    contract_version::Symbol
    required_inputs::Vector{Symbol}
    required_outputs::Vector{Symbol}
    artifact_schema_version::Symbol
    manifest_schema_version::Symbol

    function PipelineIOContract(
        contract_version,
        required_inputs::AbstractVector,
        required_outputs::AbstractVector,
        artifact_schema_version,
        manifest_schema_version,
    )
        contract_version isa Symbol || throw(ArgumentError("contract_version must be Symbol, got $(typeof(contract_version))"))
        artifact_schema_version isa Symbol || throw(ArgumentError("artifact_schema_version must be Symbol, got $(typeof(artifact_schema_version))"))
        manifest_schema_version isa Symbol || throw(ArgumentError("manifest_schema_version must be Symbol, got $(typeof(manifest_schema_version))"))
        return new(
            contract_version,
            _coerce_symbol_vector(required_inputs; field_name="required_inputs"),
            _coerce_symbol_vector(required_outputs; field_name="required_outputs"),
            artifact_schema_version,
            manifest_schema_version,
        )
    end
end

function _coerce_symbol_vector(values::AbstractVector; field_name::AbstractString)
    syms = Symbol[]
    for value in values
        value isa Symbol || throw(ArgumentError("$(field_name) entries must be Symbol, got $(typeof(value))"))
        push!(syms, value)
    end
    return syms
end

struct PipelineProvenance
    git_commit::String
    config_hash::String
    run_id::String
    timestamp::DateTime

    function PipelineProvenance(
        git_commit::AbstractString,
        config_hash::AbstractString,
        run_id::AbstractString,
        timestamp::DateTime,
    )
        isempty(strip(git_commit)) && throw(ArgumentError("git_commit must be non-empty"))
        isempty(strip(config_hash)) && throw(ArgumentError("config_hash must be non-empty"))
        isempty(strip(run_id)) && throw(ArgumentError("run_id must be non-empty"))
        return new(String(git_commit), String(config_hash), String(run_id), timestamp)
    end
end

struct PipelineSpec{P<:NamedTuple, I}
    name::String
    version::String
    model_kind::Symbol
    stages::Vector{Symbol}
    params::P
    io_contract::I

    function PipelineSpec(
        name::AbstractString,
        version::AbstractString,
        model_kind::Symbol,
        stages::AbstractVector,
        params::P,
        io_contract::I,
    ) where {P<:NamedTuple, I}
        io_contract isa AbstractPipelineIOContract || io_contract isa NamedTuple || throw(ArgumentError("io_contract must be AbstractPipelineIOContract or NamedTuple, got $(typeof(io_contract))"))
        isempty(strip(name)) && throw(ArgumentError("name must be non-empty"))
        isempty(strip(version)) && throw(ArgumentError("version must be non-empty"))

        stage_syms = _coerce_symbol_vector(stages; field_name="stages")
        isempty(stage_syms) && throw(ArgumentError("stages must be non-empty"))

        return new{P, I}(String(name), String(version), model_kind, stage_syms, params, io_contract)
    end
end

struct PipelineStage{F}
    id::Symbol
    requires::Vector{Symbol}
    provides::Vector{Symbol}
    run!::F

    function PipelineStage(id, requires::AbstractVector, provides::AbstractVector, run!::F) where {F}
        id isa Symbol || throw(ArgumentError("stage id must be Symbol, got $(typeof(id))"))
        is_callable = (run! isa Function) || !isempty(methods(run!))
        is_callable || throw(ArgumentError("stage run! must be callable, got $(typeof(run!))"))
        return new{F}(
            id,
            _coerce_symbol_vector(requires; field_name="requires"),
            _coerce_symbol_vector(provides; field_name="provides"),
            run!,
        )
    end
end

struct PipelineContext
    state::Dict{Symbol, Any}
    provenance::PipelineProvenance
end

function PipelineContext(state::AbstractDict, provenance)
    keys_are_symbol = all(k -> k isa Symbol, keys(state))
    keys_are_symbol || throw(ArgumentError("state keys must be Symbol"))

    typed_state = Dict{Symbol, Any}()
    for (k, v) in pairs(state)
        typed_state[k] = v
    end

    return PipelineContext(typed_state, _coerce_pipeline_provenance(provenance))
end

@inline function _coerce_pipeline_provenance(provenance::PipelineProvenance)
    return provenance
end

function _coerce_pipeline_provenance(provenance::NamedTuple)
    hasproperty(provenance, :git_commit) || throw(ArgumentError("provenance missing :git_commit"))
    hasproperty(provenance, :config_hash) || throw(ArgumentError("provenance missing :config_hash"))
    hasproperty(provenance, :run_id) || throw(ArgumentError("provenance missing :run_id"))
    hasproperty(provenance, :timestamp) || throw(ArgumentError("provenance missing :timestamp"))

    timestamp = provenance.timestamp
    timestamp isa DateTime || throw(ArgumentError("provenance.timestamp must be DateTime, got $(typeof(timestamp))"))
    return PipelineProvenance(
        provenance.git_commit,
        provenance.config_hash,
        provenance.run_id,
        timestamp,
    )
end

_coerce_pipeline_provenance(provenance) =
    throw(ArgumentError("provenance must be PipelineProvenance or NamedTuple, got $(typeof(provenance))"))

struct PipelineArtifact
    path::String
    hash::String
    schema_version::String

    function PipelineArtifact(path::AbstractString, hash::AbstractString, schema_version::AbstractString)
        isempty(strip(path)) && throw(ArgumentError("path must be non-empty"))
        isempty(strip(hash)) && throw(ArgumentError("hash must be non-empty"))
        isempty(strip(schema_version)) && throw(ArgumentError("schema_version must be non-empty"))
        return new(String(path), String(hash), String(schema_version))
    end
end

struct StageResult
    produced::Dict{Symbol, Any}
    artifacts::Vector{PipelineArtifact}
    metrics::Dict{Symbol, Float64}
end

function StageResult(
    produced::AbstractDict,
    artifacts::AbstractVector,
    metrics::AbstractDict,
)
    produced_keys_ok = all(k -> k isa Symbol, keys(produced))
    produced_keys_ok || throw(ArgumentError("produced keys must be Symbol"))

    metrics_keys_ok = all(k -> k isa Symbol, keys(metrics))
    metrics_keys_ok || throw(ArgumentError("metrics keys must be Symbol"))

    produced_typed = Dict{Symbol, Any}()
    for (k, v) in pairs(produced)
        produced_typed[k] = v
    end

    metrics_typed = Dict{Symbol, Float64}()
    for (k, v) in pairs(metrics)
        v isa Real || throw(ArgumentError("metric $(k) must be Real, got $(typeof(v))"))
        metrics_typed[k] = Float64(v)
    end

    artifacts_typed = PipelineArtifact[]
    for a in artifacts
        a isa PipelineArtifact || throw(ArgumentError("artifacts entries must be PipelineArtifact, got $(typeof(a))"))
        push!(artifacts_typed, a)
    end

    return StageResult(produced_typed, artifacts_typed, metrics_typed)
end

function persisted_symbol_to_string(value::Symbol)::String
    return _normalize_persisted_key(String(value))
end

function persisted_string_to_symbol(value::AbstractString)::Symbol
    return Symbol(_normalize_persisted_key(value))
end

function _normalize_persisted_key(raw::AbstractString)::String
    stripped = strip(raw)
    isempty(stripped) && throw(ArgumentError("persisted key must be non-empty"))

    for ch in stripped
        isascii(ch) || throw(ArgumentError("persisted key must be ASCII"))
        if !(isletter(ch) || isnumeric(ch) || ch in ('_', '-', ' '))
            throw(ArgumentError("persisted key contains invalid character: $(repr(ch))"))
        end
    end

    lowered = lowercase(stripped)
    replaced = replace(lowered, '-' => '_', ' ' => '_')
    collapsed = replace(replaced, r"_+" => "_")
    normalized = strip(collapsed, '_')
    isempty(normalized) && throw(ArgumentError("persisted key must contain at least one alnum character"))

    return normalized
end
