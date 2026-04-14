using Dates
using SHA

function _scan_adapter_git_commit()
    try
        return readchomp(`git -C $(joinpath(@__DIR__, "..", "..", "..")) rev-parse HEAD`)
    catch
        return "unknown"
    end
end

function _scan_adapter_config_hash(kind::Symbol, kwargs::NamedTuple)
    payload = String[]
    push!(payload, String(kind))
    for key in sort(collect(keys(kwargs)); by=String)
        push!(payload, String(key) * "=" * sprint(show, getproperty(kwargs, key)))
    end
    return bytes2hex(SHA.sha2_256(join(payload, "|")))
end

function _scan_adapter_output_dir(scan_kwargs::NamedTuple)
    if hasproperty(scan_kwargs, :output_dir)
        return String(getproperty(scan_kwargs, :output_dir))
    end
    if hasproperty(scan_kwargs, :output_path)
        return dirname(String(getproperty(scan_kwargs, :output_path)))
    end
    return mktempdir()
end

function _run_scan(kind::Symbol, scan_kwargs::NamedTuple)
    if kind === :tmu
        return run_tmu_scan(; scan_kwargs...)
    elseif kind === :trho
        return run_trho_scan(; scan_kwargs...)
    end
    throw(ArgumentError("unsupported scan kind: $(kind)"))
end

function _build_scan_pipeline_stages()
    stages = PipelineStage[]

    push!(stages,
        PipelineStage(
            :prepare_inputs,
            [:scan_inputs],
            [:prepared_scan_inputs],
            (ctx) -> StageResult(
                Dict{Symbol, Any}(:prepared_scan_inputs => ctx.state[:scan_inputs]),
                PipelineArtifact[],
                Dict{Symbol, Float64}(),
            ),
        ),
    )

    push!(stages,
        PipelineStage(
            :solve_core,
            [:prepared_scan_inputs],
            [:scan_stats],
            (ctx) -> begin
                prepared = ctx.state[:prepared_scan_inputs]
                stats = _run_scan(prepared.kind, prepared.kwargs)
                StageResult(
                    Dict{Symbol, Any}(:scan_stats => stats),
                    PipelineArtifact[],
                    Dict{Symbol, Float64}(),
                )
            end,
        ),
    )

    push!(stages,
        PipelineStage(
            :postprocess,
            [:scan_stats],
            [:scan_summary],
            (ctx) -> StageResult(
                Dict{Symbol, Any}(:scan_summary => ctx.state[:scan_stats]),
                PipelineArtifact[],
                Dict{Symbol, Float64}(),
            ),
        ),
    )

    push!(stages,
        PipelineStage(
            :export_artifacts,
            [:scan_summary, :output_dir],
            [:scan_artifact_paths],
            (ctx) -> begin
                scan_output = hasproperty(ctx.state[:scan_summary], :output) ? String(ctx.state[:scan_summary].output) : ""
                artifact_paths = Dict{String, String}("run_manifest" => joinpath(ctx.state[:output_dir], "run_manifest.json"))
                if !isempty(scan_output)
                    artifact_paths["scan_output"] = scan_output
                end
                StageResult(
                    Dict{Symbol, Any}(:scan_artifact_paths => artifact_paths),
                    PipelineArtifact[],
                    Dict{Symbol, Float64}(),
                )
            end,
        ),
    )

    push!(stages,
        PipelineStage(
            :emit_diagnostics,
            [:scan_inputs, :scan_summary, :scan_artifact_paths],
            [:scan_diagnostics],
            (ctx) -> begin
                scan_inputs = ctx.state[:scan_inputs]
                summary = ctx.state[:scan_summary]
                diagnostics = (
                    status=:success,
                    mode=:builtin,
                    kind=scan_inputs.kind,
                    total=getproperty(summary, :total),
                    success=getproperty(summary, :success),
                    failure=getproperty(summary, :failure),
                    skipped=getproperty(summary, :skipped),
                )

                extensions = ctx.state[:manifest_extensions]
                extensions["diagnostics_mode"] = String(diagnostics.mode)
                extensions["diagnostics_status"] = String(diagnostics.status)

                StageResult(
                    Dict{Symbol, Any}(:scan_diagnostics => diagnostics),
                    PipelineArtifact[],
                    Dict{Symbol, Float64}(),
                )
            end,
        ),
    )

    push!(stages,
        PipelineStage(
            :emit_repro_manifest,
            [:scan_summary, :scan_artifact_paths, :scan_diagnostics],
            [:scan_pipeline_output],
            (ctx) -> begin
                summary = ctx.state[:scan_summary]
                output = merge(
                    (; summary...),
                    (
                        manifest_path=ctx.state[:scan_artifact_paths]["run_manifest"],
                        diagnostics=ctx.state[:scan_diagnostics],
                    ),
                )
                StageResult(
                    Dict{Symbol, Any}(:scan_pipeline_output => output),
                    PipelineArtifact[],
                    Dict{Symbol, Float64}(),
                )
            end,
        ),
    )

    return stages
end

function run_scan_pipeline_adapter(kind::Symbol; kwargs...)
    (kind === :tmu || kind === :trho) || throw(ArgumentError("unsupported scan kind: $(kind)"))
    scan_kwargs = (; kwargs...)
    output_dir = _scan_adapter_output_dir(scan_kwargs)
    manifest_path = joinpath(output_dir, "run_manifest.json")
    run_id = basename(normpath(output_dir))
    model_kind = hasproperty(scan_kwargs, :model_kind) ? Symbol(getproperty(scan_kwargs, :model_kind)) : :PNJL

    ctx = PipelineContext(
        Dict{Symbol, Any}(
            :scan_inputs => (kind=kind, kwargs=scan_kwargs),
            :output_dir => output_dir,
            :manifest_extensions => build_manifest_extensions((
                pipeline_family="scan",
                baseline_suite="smoke",
                physics_profile=get(ENV, "PHYSICS_PARAM_PROFILE", "default"),
                adapter_version="v1",
            )),
        ),
        PipelineProvenance(
            _scan_adapter_git_commit(),
            _scan_adapter_config_hash(kind, scan_kwargs),
            run_id,
            Dates.now(Dates.UTC),
        ),
    )

    spec = PipelineSpec(
        "scan_pipeline_runner",
        "v1",
        model_kind,
        collect(scan_pipeline_stage_ids()),
        (; kind=kind),
        PipelineIOContract(
            :v1,
            [:scan_inputs],
            [:scan_pipeline_output],
            :scan_artifact_v1,
            :scan_manifest_v1,
        ),
    )

    run_result = run_pipeline(spec, _build_scan_pipeline_stages(), ctx; manifest_path=manifest_path)
    if !run_result.success
        run_result.error_kind == :ArgumentError && throw(ArgumentError(run_result.error_msg))
        throw(ErrorException("scan runner stage failed at $(run_result.failed_stage): $(run_result.error_msg)"))
    end

    return ctx.state[:scan_pipeline_output]
end
