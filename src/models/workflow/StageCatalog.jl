using Dates
using SHA

const PHASE_PIPELINE_STAGE_IDS = (
    :build_model,
    :prepare_grid,
    :solve_points,
    :collect_diagnostics,
    :analyze_phase,
    :export_artifacts,
    :emit_repro_manifest,
)

function _resolve_phase_runner_output_dir(
        model_kind::Symbol;
        output_dir::Union{Nothing, String},
        profile,
        run_id::Union{Nothing, String},
        policy::Symbol)
    target = resolve_phase_output_target(model_kind; profile=profile, run_id=run_id, policy=policy)
    return isnothing(output_dir) ? target.run_dir : output_dir
end

function _phase_runner_run_id(run_id::Union{Nothing, String}, run_dir::String)
    if run_id === nothing || isempty(strip(run_id))
        return basename(normpath(run_dir))
    end
    return String(run_id)
end

function _phase_runner_git_commit()
    try
        return readchomp(`git -C $(joinpath(@__DIR__, "..", "..", "..")) rev-parse HEAD`)
    catch
        return "unknown"
    end
end

function _phase_runner_config_hash(model_kind::Symbol, phase_kwargs)
    payload = String[]
    push!(payload, String(model_kind))
    for key in sort(collect(keys(phase_kwargs)); by=String)
        push!(payload, String(key) * "=" * sprint(show, getproperty(phase_kwargs, key)))
    end
    return bytes2hex(SHA.sha2_256(join(payload, "|")))
end

function _build_phase_pipeline_stages(core_run_phase_pipeline)
    build_model_stage = PipelineStage(
        :build_model,
        [:model_kind],
        [:phase_model_kind],
        (ctx) -> StageResult(Dict{Symbol, Any}(:phase_model_kind => ctx.state[:model_kind]), PipelineArtifact[], Dict{Symbol, Float64}()),
    )

    prepare_grid_stage = PipelineStage(
        :prepare_grid,
        [:phase_kwargs],
        [:phase_kwargs_prepared],
        (ctx) -> StageResult(Dict{Symbol, Any}(:phase_kwargs_prepared => ctx.state[:phase_kwargs]), PipelineArtifact[], Dict{Symbol, Float64}()),
    )

    solve_points_stage = PipelineStage(
        :solve_points,
        [:phase_model_kind, :phase_kwargs_prepared],
        [:phase_result],
        (ctx) -> begin
            phase_kwargs = ctx.state[:phase_kwargs_prepared]
            phase_result = core_run_phase_pipeline(ctx.state[:phase_model_kind]; phase_kwargs...)
            StageResult(Dict{Symbol, Any}(:phase_result => phase_result), PipelineArtifact[], Dict{Symbol, Float64}())
        end,
    )

    collect_diagnostics_stage = PipelineStage(
        :collect_diagnostics,
        [:phase_result],
        [:phase_diagnostics],
        (ctx) -> StageResult(Dict{Symbol, Any}(:phase_diagnostics => ctx.state[:phase_result].diagnostics), PipelineArtifact[], Dict{Symbol, Float64}()),
    )

    analyze_phase_stage = PipelineStage(
        :analyze_phase,
        [:phase_result],
        [:phase_analysis],
        (ctx) -> begin
            analysis = Dict{String, Any}(
                "boundary_count" => length(ctx.state[:phase_result].first_order_boundary),
                "spinodal_count" => length(ctx.state[:phase_result].spinodal),
                "crossover_count" => length(ctx.state[:phase_result].crossover_line),
            )
            StageResult(Dict{Symbol, Any}(:phase_analysis => analysis), PipelineArtifact[], Dict{Symbol, Float64}())
        end,
    )

    export_artifacts_stage = PipelineStage(
        :export_artifacts,
        [:phase_result],
        [:phase_artifact_paths],
        (ctx) -> StageResult(Dict{Symbol, Any}(:phase_artifact_paths => ctx.state[:phase_result].artifact_paths), PipelineArtifact[], Dict{Symbol, Float64}()),
    )

    emit_repro_manifest_stage = PipelineStage(
        :emit_repro_manifest,
        [:phase_result],
        [:phase_pipeline_output],
        (ctx) -> StageResult(Dict{Symbol, Any}(:phase_pipeline_output => ctx.state[:phase_result]), PipelineArtifact[], Dict{Symbol, Float64}()),
    )

    return PipelineStage[
        build_model_stage,
        prepare_grid_stage,
        solve_points_stage,
        collect_diagnostics_stage,
        analyze_phase_stage,
        export_artifacts_stage,
        emit_repro_manifest_stage,
    ]
end

function run_phase_pipeline_via_runner(core_run_phase_pipeline, model_kind::Symbol=:PNJL; kwargs...)
    profile = get(kwargs, :profile, :default)
    run_id = get(kwargs, :run_id, nothing)
    policy = get(kwargs, :policy, :processed_first)
    output_dir = get(kwargs, :output_dir, nothing)

    run_dir = _resolve_phase_runner_output_dir(model_kind;
        output_dir=output_dir,
        profile=profile,
        run_id=run_id,
        policy=policy)
    manifest_path = joinpath(run_dir, "run_manifest.json")

    core_kwargs = merge((; output_dir=run_dir), (; kwargs...))
    effective_run_id = _phase_runner_run_id(run_id, run_dir)
    provenance_config_hash = _phase_runner_config_hash(model_kind, core_kwargs)

    ctx = PipelineContext(
        Dict{Symbol, Any}(
            :model_kind => model_kind,
            :phase_kwargs => core_kwargs,
        ),
        PipelineProvenance(
            _phase_runner_git_commit(),
            provenance_config_hash,
            effective_run_id,
            Dates.now(Dates.UTC),
        ),
    )

    spec = PipelineSpec(
        "phase_pipeline_runner",
        "v1",
        model_kind,
        collect(PHASE_PIPELINE_STAGE_IDS),
        (;),
        PipelineIOContract(
            :v1,
            [:model_kind, :phase_kwargs],
            [:phase_pipeline_output],
            :phase_artifact_v1,
            :phase_manifest_v1,
        ),
    )

    stages = _build_phase_pipeline_stages(core_run_phase_pipeline)
    run_result = run_pipeline(spec, stages, ctx; manifest_path=manifest_path)
    if !run_result.success
        throw(ErrorException("phase runner stage failed at $(run_result.failed_stage): $(run_result.error_msg)"))
    end

    return ctx.state[:phase_pipeline_output]
end
