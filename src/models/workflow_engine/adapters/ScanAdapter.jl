using Dates


function _run_scan(kind::Symbol, scan_kwargs::NamedTuple)
    if kind === :tmu
        return run_tmu_scan(; scan_kwargs...)
    elseif kind === :trho
        return run_trho_scan(; scan_kwargs...)
    elseif kind === :magnetic
        return run_magnetic_scan(; scan_kwargs...)
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
                if hasproperty(ctx.state[:scan_summary], :selected_output)
                    selected_output = String(ctx.state[:scan_summary].selected_output)
                    artifact_paths["scan_selected_output"] = selected_output
                end
                if hasproperty(ctx.state[:scan_summary], :candidates_output)
                    candidates_output = String(ctx.state[:scan_summary].candidates_output)
                    artifact_paths["scan_candidates_output"] = candidates_output
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
                summary = ctx.state[:scan_summary]
                hasproperty(summary, :selected_output) &&
                    (extensions["scan_output_path"] = String(summary.selected_output))
                hasproperty(summary, :candidates_output) &&
                    (extensions["scan_candidates_output_path"] = String(summary.candidates_output))

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
    (kind === :tmu || kind === :trho || kind === :magnetic) || throw(ArgumentError("unsupported scan kind: $(kind)"))
    scan_kwargs = (; kwargs...)
    output_dir = resolve_adapter_output_dir(scan_kwargs, (:output_dir, :output_path); path_keys=(:output_path,))
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
            adapter_git_commit(),
            adapter_config_hash(kind, scan_kwargs),
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
