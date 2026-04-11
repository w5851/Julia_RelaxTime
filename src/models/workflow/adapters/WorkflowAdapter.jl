using Dates
using SHA

const _WORKFLOW_ADAPTER_STAGE_IDS = [
    :prepare_inputs,
    :solve_core,
    :postprocess,
    :export_artifacts,
    :emit_repro_manifest,
]

@inline function _workflow_adapter_output_dir(output_dir)
    output_dir === nothing && return mktempdir()
    return String(output_dir)
end

function _workflow_adapter_git_commit()
    try
        return readchomp(`git -C $(joinpath(@__DIR__, "..", "..", "..")) rev-parse HEAD`)
    catch
        return "unknown"
    end
end

function _workflow_adapter_config_hash(kind::Symbol, prepared::NamedTuple)
    payload = String[]
    push!(payload, String(kind))
    for key in sort(collect(keys(prepared)); by=String)
        push!(payload, String(key) * "=" * sprint(show, getproperty(prepared, key)))
    end
    return bytes2hex(SHA.sha2_256(join(payload, "|")))
end

function _resolve_transport_inputs(kwargs::NamedTuple)
    hasproperty(kwargs, :T_fm) || throw(ArgumentError("run_workflow_pipeline(:transport) requires keyword T_fm"))
    hasproperty(kwargs, :mu_fm) || throw(ArgumentError("run_workflow_pipeline(:transport) requires keyword mu_fm"))

    T_fm = Float64(kwargs.T_fm)
    mu_fm = Float64(kwargs.mu_fm)
    xi = hasproperty(kwargs, :xi) ? Float64(kwargs.xi) : 0.0
    output_dir = _workflow_adapter_output_dir(get(kwargs, :output_dir, nothing))

    passthrough_pairs = Pair{Symbol, Any}[]
    for (k, v) in pairs(kwargs)
        if k in (:T_fm, :mu_fm, :xi, :output_dir)
            continue
        end
        push!(passthrough_pairs, k => v)
    end

    return (
        T_fm=T_fm,
        mu_fm=mu_fm,
        xi=xi,
        output_dir=output_dir,
        passthrough=(; passthrough_pairs...),
    )
end

function _build_workflow_stages()
    stages = PipelineStage[]

    push!(stages,
        PipelineStage(
            :prepare_inputs,
            [:workflow_inputs],
            [:prepared_inputs],
            (ctx) -> StageResult(
                Dict{Symbol, Any}(:prepared_inputs => ctx.state[:workflow_inputs]),
                PipelineArtifact[],
                Dict{Symbol, Float64}(),
            ),
        ),
    )

    push!(stages,
        PipelineStage(
            :solve_core,
            [:prepared_inputs],
            [:workflow_result, :transport],
            (ctx) -> begin
                prepared = ctx.state[:prepared_inputs]
                passthrough = prepared.passthrough
                result = solve_gap_and_transport(
                    prepared.T_fm,
                    prepared.mu_fm;
                    xi=prepared.xi,
                    passthrough...
                )
                StageResult(
                    Dict{Symbol, Any}(
                        :workflow_result => result,
                        :transport => result.transport,
                    ),
                    PipelineArtifact[],
                    Dict{Symbol, Float64}(),
                )
            end,
        ),
    )

    push!(stages,
        PipelineStage(
            :postprocess,
            [:workflow_result],
            [:postprocessed],
            (ctx) -> StageResult(
                Dict{Symbol, Any}(:postprocessed => ctx.state[:workflow_result]),
                PipelineArtifact[],
                Dict{Symbol, Float64}(),
            ),
        ),
    )

    push!(stages,
        PipelineStage(
            :export_artifacts,
            [:postprocessed, :output_dir],
            [:artifact_paths],
            (ctx) -> begin
                manifest_path = joinpath(ctx.state[:output_dir], "run_manifest.json")
                StageResult(
                    Dict{Symbol, Any}(:artifact_paths => Dict{String, String}("run_manifest" => manifest_path)),
                    PipelineArtifact[],
                    Dict{Symbol, Float64}(),
                )
            end,
        ),
    )

    push!(stages,
        PipelineStage(
            :emit_repro_manifest,
            [:transport, :artifact_paths],
            [:workflow_pipeline_output],
            (ctx) -> begin
                run_manifest = ctx.state[:artifact_paths]["run_manifest"]
                output = (
                    transport=ctx.state[:transport],
                    run_manifest=run_manifest,
                )
                StageResult(
                    Dict{Symbol, Any}(:workflow_pipeline_output => output),
                    PipelineArtifact[],
                    Dict{Symbol, Float64}(),
                )
            end,
        ),
    )

    return stages
end

function _run_transport_workflow_pipeline(; kwargs...)
    prepared = _resolve_transport_inputs((; kwargs...))
    manifest_path = joinpath(prepared.output_dir, "run_manifest.json")
    run_id = basename(normpath(prepared.output_dir))

    ctx = PipelineContext(
        Dict{Symbol, Any}(
            :workflow_inputs => prepared,
            :output_dir => prepared.output_dir,
            :manifest_extensions => build_manifest_extensions((
                pipeline_family="workflow",
                baseline_suite="smoke",
                physics_profile=get(ENV, "PHYSICS_PARAM_PROFILE", "default"),
                adapter_version="v1",
            )),
        ),
        PipelineProvenance(
            _workflow_adapter_git_commit(),
            _workflow_adapter_config_hash(:transport, prepared),
            run_id,
            Dates.now(Dates.UTC),
        ),
    )

    spec = PipelineSpec(
        "workflow_adapter_transport",
        "v1",
        :PNJL,
        _WORKFLOW_ADAPTER_STAGE_IDS,
        (; kind=:transport),
        PipelineIOContract(
            :v1,
            [:workflow_inputs],
            [:workflow_pipeline_output],
            :workflow_artifact_v1,
            :workflow_manifest_v1,
        ),
    )

    run_result = run_pipeline(spec, _build_workflow_stages(), ctx; manifest_path=manifest_path)
    if !run_result.success
        run_result.error_kind == :ArgumentError && throw(ArgumentError(run_result.error_msg))
        throw(ErrorException("workflow adapter pipeline failed at $(run_result.failed_stage): $(run_result.error_msg)"))
    end

    return ctx.state[:workflow_pipeline_output]
end

function run_workflow_pipeline_adapter(kind::Symbol; kwargs...)
    kind === :transport || throw(ArgumentError("unsupported workflow kind: $(kind)"))
    return _run_transport_workflow_pipeline(; kwargs...)
end
