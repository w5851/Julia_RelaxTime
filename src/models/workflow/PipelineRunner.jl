using Dates
using JSON3
using SHA

struct PipelineStageRecord
    id::Symbol
    status::Symbol
    started_at::Union{Nothing, String}
    ended_at::Union{Nothing, String}
    error_kind::Union{Nothing, Symbol}
    error_msg::Union{Nothing, String}
end

struct PipelineRunResult
    success::Bool
    failed_stage::Union{Nothing, Symbol}
    error_kind::Union{Nothing, Symbol}
    error_msg::Union{Nothing, String}
    completed_stages::Vector{Symbol}
    stage_records::Vector{PipelineStageRecord}
    manifest_path::String
    config_hash::String
    artifact_hash::String
end

@inline function _utc_now_iso8601()
    return Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sssZ")
end

@inline _error_kind_symbol(err) = Symbol(nameof(typeof(err)))

function compute_pipeline_config_hash(spec::PipelineSpec, ctx::PipelineContext, stages::Vector{<:PipelineStage})
    stage_tokens = String[]
    for stage in sort(stages; by=s -> String(s.id))
        token = String(stage.id) * ":" * join(sort(String.(stage.requires)), ",") * ":" * join(sort(String.(stage.provides)), ",")
        push!(stage_tokens, token)
    end
    payload = join([
            String(spec.model_kind),
            spec.name,
            spec.version,
            join(sort(String.(spec.stages)), ","),
            String(ctx.provenance.git_commit),
            join(stage_tokens, ";"),
        ], "|")
    return bytes2hex(sha2_256(payload))
end

function compute_pipeline_artifact_hash(state::Dict{Symbol, Any}, stage_records::Vector{PipelineStageRecord})
    state_tokens = String[]
    for key in sort(collect(keys(state)); by=String)
        value_repr = _stable_hash_repr(state[key])
        push!(state_tokens, String(key) * "=" * value_repr)
    end
    record_tokens = String[]
    for rec in sort(stage_records; by=r -> String(r.id))
        push!(record_tokens,
            String(rec.id) * ":" * String(rec.status) * ":" *
            (isnothing(rec.error_kind) ? "none" : String(rec.error_kind)) * ":" *
            (isnothing(rec.error_msg) ? "" : String(rec.error_msg)))
    end
    payload = join([join(state_tokens, ";"), join(record_tokens, ";")], "|")
    return bytes2hex(sha2_256(payload))
end

function _stable_hash_repr(value)
    if value === nothing
        return "null"
    elseif value isa Bool
        return value ? "true" : "false"
    elseif value isa Real
        return isnan(value) ? "NaN" : string(value)
    elseif value isa AbstractString
        return repr(String(value))
    elseif value isa Symbol
        return "symbol:" * String(value)
    elseif value isa DateTime
        return Dates.format(value, dateformat"yyyy-mm-ddTHH:MM:SS.sss")
    elseif value isa AbstractVector
        return "[" * join((_stable_hash_repr(v) for v in value), ",") * "]"
    elseif value isa AbstractSet
        normalized = sort([_stable_hash_repr(v) for v in value])
        return "set{" * join(normalized, ",") * "}"
    elseif value isa AbstractDict
        tokens = String[]
        for key in sort(collect(keys(value)); by=k -> _stable_hash_repr(k))
            push!(tokens, _stable_hash_repr(key) * ":" * _stable_hash_repr(value[key]))
        end
        return "{" * join(tokens, ",") * "}"
    elseif value isa NamedTuple
        names = sort!(collect(keys(value)); by=String)
        tokens = String[]
        for name in names
            push!(tokens, String(name) * ":" * _stable_hash_repr(getproperty(value, name)))
        end
        return "(" * join(tokens, ",") * ")"
    end
    return repr(value)
end

function _validate_stage_id_uniqueness(stages::Vector{<:PipelineStage})
    seen = Set{Symbol}()
    for stage in stages
        if stage.id in seen
            throw(ArgumentError("duplicate stage id: $(stage.id)"))
        end
        push!(seen, stage.id)
    end
    return nothing
end

function _validate_duplicate_provides(stages::Vector{<:PipelineStage})
    provider_by_symbol = Dict{Symbol, Symbol}()
    for stage in stages
        for sym in stage.provides
            if haskey(provider_by_symbol, sym)
                throw(ArgumentError("duplicate provides symbol without merge policy: $(sym)"))
            end
            provider_by_symbol[sym] = stage.id
        end
    end
    return provider_by_symbol
end

function _validate_missing_dependencies(stages::Vector{<:PipelineStage}, ctx::PipelineContext)
    known = Set{Symbol}(keys(ctx.state))
    for stage in stages
        union!(known, stage.provides)
    end
    for stage in stages
        for dep in stage.requires
            dep in known || throw(ArgumentError("missing required dependency $(dep) for stage $(stage.id)"))
        end
    end
    return nothing
end

function _build_stage_graph(stages::Vector{<:PipelineStage}, provider_by_symbol::Dict{Symbol, Symbol})
    outgoing = Dict{Symbol, Set{Symbol}}(stage.id => Set{Symbol}() for stage in stages)
    indegree = Dict{Symbol, Int}(stage.id => 0 for stage in stages)

    for stage in stages
        for dep in stage.requires
            if haskey(provider_by_symbol, dep)
                provider = provider_by_symbol[dep]
                if provider != stage.id && !(stage.id in outgoing[provider])
                    push!(outgoing[provider], stage.id)
                    indegree[stage.id] += 1
                end
            end
        end
    end
    return outgoing, indegree
end

function _validate_no_cycles(stages::Vector{<:PipelineStage}, provider_by_symbol::Dict{Symbol, Symbol})
    outgoing, indegree = _build_stage_graph(stages, provider_by_symbol)
    pending = sort([id for (id, deg) in indegree if deg == 0]; by=String)
    visited = 0
    while !isempty(pending)
        id = popfirst!(pending)
        visited += 1
        for nxt in sort(collect(outgoing[id]); by=String)
            indegree[nxt] -= 1
            if indegree[nxt] == 0
                push!(pending, nxt)
                sort!(pending; by=String)
            end
        end
    end
    visited == length(stages) || throw(ArgumentError("cyclic dependencies detected in pipeline stages"))
    return nothing
end

function _validate_stage_coverage(spec::PipelineSpec, stages::Vector{<:PipelineStage})
    stage_ids = Set{Symbol}(stage.id for stage in stages)
    spec_ids = Set{Symbol}(spec.stages)
    for id in spec.stages
        id in stage_ids || throw(ArgumentError("spec stage not provided: $(id)"))
    end
    for id in stage_ids
        id in spec_ids || throw(ArgumentError("stage provided but not listed in spec: $(id)"))
    end
    return nothing
end

function _resolve_manifest_schema_version(io_contract)::String
    if io_contract isa AbstractPipelineIOContract
        return String(getproperty(io_contract, :manifest_schema_version))
    end
    io_contract isa NamedTuple || throw(ArgumentError("io_contract must be AbstractPipelineIOContract or NamedTuple, got $(typeof(io_contract))"))
    hasproperty(io_contract, :manifest_schema_version) || throw(ArgumentError("io_contract missing :manifest_schema_version"))
    value = getproperty(io_contract, :manifest_schema_version)
    value isa Symbol || throw(ArgumentError("io_contract.manifest_schema_version must be Symbol, got $(typeof(value))"))
    return String(value)
end

function _resolve_required_outputs(io_contract)::Vector{Symbol}
    raw = if io_contract isa AbstractPipelineIOContract
        getproperty(io_contract, :required_outputs)
    elseif io_contract isa NamedTuple
        hasproperty(io_contract, :required_outputs) || throw(ArgumentError("io_contract missing :required_outputs"))
        getproperty(io_contract, :required_outputs)
    else
        throw(ArgumentError("io_contract must be AbstractPipelineIOContract or NamedTuple, got $(typeof(io_contract))"))
    end
    raw isa AbstractVector || throw(ArgumentError("io_contract.required_outputs must be AbstractVector, got $(typeof(raw))"))
    return Symbol.(raw)
end

function _validate_stage_runtime_inputs(stage::PipelineStage, ctx::PipelineContext)
    missing = Symbol[]
    for dep in stage.requires
        haskey(ctx.state, dep) || push!(missing, dep)
    end
    isempty(missing) || throw(ArgumentError("missing runtime dependencies for stage $(stage.id): $(join(String.(missing), ", "))"))
    return nothing
end

function _record_to_manifest(record::PipelineStageRecord)
    return Dict(
        "id" => String(record.id),
        "status" => String(record.status),
        "started_at" => record.started_at,
        "ended_at" => record.ended_at,
        "error_kind" => isnothing(record.error_kind) ? nothing : String(record.error_kind),
        "error_msg" => record.error_msg,
    )
end

function _write_manifest(
        manifest_path::String,
        spec::PipelineSpec,
        ctx::PipelineContext,
        result::PipelineRunResult)
    pipeline_meta = Dict(
        "name" => spec.name,
        "version" => spec.version,
        "model_kind" => String(spec.model_kind),
        "run_id" => ctx.provenance.run_id,
        "git_commit" => ctx.provenance.git_commit,
        "manifest_schema_version" => _resolve_manifest_schema_version(spec.io_contract),
        "timestamp" => _utc_now_iso8601(),
        "success" => result.success,
        "failed_stage" => isnothing(result.failed_stage) ? nothing : String(result.failed_stage),
        "error_kind" => isnothing(result.error_kind) ? nothing : String(result.error_kind),
        "error_msg" => result.error_msg,
        "config_hash" => result.config_hash,
        "artifact_hash" => result.artifact_hash,
    )
    if haskey(ctx.state, :manifest_extensions)
        extensions = ctx.state[:manifest_extensions]
        extensions isa AbstractDict || throw(ArgumentError("ctx.state[:manifest_extensions] must be AbstractDict, got $(typeof(extensions))"))
        merge!(pipeline_meta, normalize_manifest_extensions(extensions))
    end
    payload = Dict(
        "pipeline" => pipeline_meta,
        "completed_stages" => String.(result.completed_stages),
        "stage_records" => [_record_to_manifest(rec) for rec in result.stage_records],
    )
    mkpath(dirname(manifest_path))
    open(manifest_path, "w") do io
        write(io, JSON3.write(payload))
    end
    return nothing
end

function run_pipeline(
        spec::PipelineSpec,
        stages::AbstractVector,
        ctx::PipelineContext;
        manifest_path::String)
    typed_stages = PipelineStage[stage for stage in stages]
    _validate_stage_id_uniqueness(typed_stages)
    _validate_stage_coverage(spec, typed_stages)
    provider_by_symbol = _validate_duplicate_provides(typed_stages)
    _validate_missing_dependencies(typed_stages, ctx)
    _validate_no_cycles(typed_stages, provider_by_symbol)

    stage_by_id = Dict(stage.id => stage for stage in typed_stages)
    ordered_ids = copy(spec.stages)
    records = Dict{Symbol, PipelineStageRecord}(
        id => PipelineStageRecord(id, :pending, nothing, nothing, nothing, nothing)
        for id in ordered_ids
    )

    completed = Symbol[]
    failed_stage = nothing
    error_kind = nothing
    error_msg = nothing

    for id in ordered_ids
        stage = stage_by_id[id]
        started_at = _utc_now_iso8601()
        try
            _validate_stage_runtime_inputs(stage, ctx)
            stage_result = stage.run!(ctx)
            stage_result isa StageResult || throw(ArgumentError("stage $(id) must return StageResult"))
            merge!(ctx.state, stage_result.produced)
            records[id] = PipelineStageRecord(id, :completed, started_at, _utc_now_iso8601(), nothing, nothing)
            push!(completed, id)
        catch err
            failed_stage = id
            error_kind = _error_kind_symbol(err)
            error_msg = sprint(showerror, err)
            records[id] = PipelineStageRecord(id, :failed, started_at, _utc_now_iso8601(), error_kind, error_msg)
            break
        end
    end

    if failed_stage !== nothing
        seen_failure = false
        for id in ordered_ids
            if id == failed_stage
                seen_failure = true
                continue
            end
            seen_failure || continue
            rec = records[id]
            if rec.status == :pending
                records[id] = PipelineStageRecord(id, :skipped, nothing, nothing, nothing, nothing)
            end
        end
    else
        required_outputs = _resolve_required_outputs(spec.io_contract)
        missing = Symbol[]
        for sym in required_outputs
            haskey(ctx.state, sym) || push!(missing, sym)
        end
        if !isempty(missing)
            error_kind = :missing_required_outputs
            error_msg = "missing required outputs: " * join(String.(missing), ", ")
        end
    end

    ordered_records = [records[id] for id in ordered_ids]
    success = isnothing(failed_stage) && isnothing(error_kind)
    artifact_hash = compute_pipeline_artifact_hash(ctx.state, ordered_records)
    config_hash = compute_pipeline_config_hash(spec, ctx, typed_stages)

    result = PipelineRunResult(
        success,
        failed_stage,
        error_kind,
        error_msg,
        completed,
        ordered_records,
        manifest_path,
        config_hash,
        artifact_hash,
    )
    _write_manifest(manifest_path, spec, ctx, result)
    return result
end
