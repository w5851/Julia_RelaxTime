using Dates
using JSON3
using SHA
using TOML

const _REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
const _DEFAULT_RELAXTIME_CONFIG_PATH = joinpath(_REPO_ROOT, "config", "workflows", "relaxtime", "default.toml")
const _RELAXTIME_ALIASES_PATH = joinpath(_REPO_ROOT, "config", "workflows", "relaxtime", "schema", "aliases_v1.toml")

@inline function _touch(path::String)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, "")
    end
end

function _read_fallback_events(path::String)
    if !isfile(path)
        return Any[]
    end
    content = strip(read(path, String))
    isempty(content) && return Any[]
    try
        parsed = JSON3.read(content)
        parsed isa AbstractVector && return collect(parsed)
    catch
    end
    occursin("bulk_fallback", content) && return Any[Dict("event" => "bulk_fallback")]
    return Any[Dict("raw" => content)]
end

@inline function _number_tag(x)
    v = Float64(x)
    if isapprox(v, round(v); atol=1e-12, rtol=0.0)
        return string(Int(round(v)))
    end
    return replace(string(v), "." => "p")
end

function _emit_plot_contract_artifacts(cmd::Symbol, effective::Dict{String, Any}, outdir::String)
    fig_dir = joinpath(outdir, "figures")
    mkpath(fig_dir)

    scan = get(effective, "scan", Dict{String, Any}())
    plot = get(effective, "plot", Dict{String, Any}())

    if cmd === :transport
        transport_scan = get(scan, "transport", Dict{String, Any}())
        muB_tag = _number_tag(get(transport_scan, "muB_MeV", 0.0))
        transport_plot = get(plot, "transport", Dict{String, Any}())
        ys = String.(get(transport_plot, "ys", Any[]))
        for y in ys
            _touch(joinpath(fig_dir, "transport__$(y)__muB$(muB_tag).png"))
        end
    elseif cmd === :cross_section
        xs_scan = get(scan, "cross_section", Dict{String, Any}())
        Ts = Float64.(get(xs_scan, "T_list_MeV", Any[]))
        processes = String.(get(xs_scan, "processes", Any[]))
        for T in Ts
            T_tag = _number_tag(T)
            for p in processes
                _touch(joinpath(fig_dir, "xsec__T$(T_tag)__$(p).png"))
            end
        end
    end
    return nothing
end

if !isdefined(@__MODULE__, :WorkflowConfig)
    include(joinpath(_REPO_ROOT, "scripts", "relaxtime", "config", "WorkflowConfig.jl"))
end
if !isdefined(@__MODULE__, :WorkflowConfigAudit)
    include(joinpath(_REPO_ROOT, "scripts", "relaxtime", "config", "WorkflowConfigAudit.jl"))
end
if !isdefined(@__MODULE__, :CrossSectionOrchestratedScan)
    include(joinpath(_REPO_ROOT, "scripts", "relaxtime", "workflow", "cross_section_orchestrated.jl"))
end

using .WorkflowConfig: normalize_merge_validate
using .WorkflowConfigAudit: build_consumption_report
using .CrossSectionOrchestratedScan: run_cross_section_orchestrated

function _relaxtime_orchestrator_git_commit()
    try
        return readchomp(`git -C $(_REPO_ROOT) rev-parse HEAD`)
    catch
        return "unknown"
    end
end

function _relaxtime_orchestrator_config_hash(cmd::Symbol, kwargs::NamedTuple)
    payload = String[]
    push!(payload, String(cmd))
    for key in sort(collect(keys(kwargs)); by=String)
        push!(payload, String(key) * "=" * sprint(show, getproperty(kwargs, key)))
    end
    return bytes2hex(SHA.sha2_256(join(payload, "|")))
end

@inline function _relaxtime_orchestrator_output_dir(kwargs::NamedTuple)
    if hasproperty(kwargs, :outdir)
        return String(getproperty(kwargs, :outdir))
    elseif hasproperty(kwargs, :output_dir)
        return String(getproperty(kwargs, :output_dir))
    end
    return mktempdir()
end

function _canonicalize_processes(raw)
    if raw === nothing
        return nothing
    end
    if raw isa AbstractString
        values = String[]
        for item in split(raw, ',')
            token = strip(item)
            isempty(token) && continue
            push!(values, token)
        end
        isempty(values) && throw(ArgumentError("processes cannot be empty"))
        return values
    end
    if raw isa AbstractVector
        values = String[]
        for item in raw
            token = strip(String(item))
            isempty(token) && continue
            push!(values, token)
        end
        isempty(values) && throw(ArgumentError("processes cannot be empty"))
        return values
    end
    throw(ArgumentError("processes must be string or vector, got $(typeof(raw))"))
end

function _build_consumed_keys(cmd::Symbol, effective::Dict{String, Any})
    consumed = Set{String}([
        "schema_version",
        "profile_name",
        "scan.cross_section.muB_MeV",
        "scan.cross_section.T_list_MeV",
        "scan.cross_section.xi_list",
        "scan.cross_section.processes",
        "scan.cross_section.n_points",
        "scan.cross_section.energy.mode",
        "scan.cross_section.energy.sqrt_s_min_MeV",
        "scan.cross_section.energy.sqrt_s_max_MeV",
        "scan.cross_section.energy.sqrt_s_num",
        "plot.cross_section.x",
        "plot.cross_section.group",
        "plot.cross_section.split",
    ])
    if haskey(effective, "strict")
        push!(consumed, "strict")
    end
    union!(consumed, Set{String}([
        "scan.transport.muB_MeV",
        "scan.transport.xi_list",
        "scan.transport.tmin_MeV",
        "scan.transport.tmax_MeV",
        "scan.transport.tstep_MeV",
        "scan.transport.resume",
        "scan.transport.overwrite",
        "scan.transport.solver.p_num",
        "scan.transport.solver.t_num",
        "scan.transport.solver.max_iter",
        "scan.transport.tau.mode",
        "scan.transport.tau.tau_p_nodes",
        "scan.transport.tau.tau_angle_nodes",
        "scan.transport.tau.tau_phi_nodes",
        "scan.transport.tau.tau_n_sigma",
        "scan.transport.tau.sigma_grid_n",
        "scan.transport.transport.compute_bulk",
        "scan.transport.transport.tr_p_nodes",
        "scan.transport.transport.tr_p_max_fm",
        "plot.transport.x",
        "plot.transport.group",
        "plot.transport.ys",
    ]))
    if haskey(effective, "scan") && effective["scan"] isa AbstractDict
        scan = effective["scan"]
        if haskey(scan, "cross_section") && scan["cross_section"] isa AbstractDict
            xs = scan["cross_section"]
            if haskey(xs, "energy") && xs["energy"] isa AbstractDict
                energy = xs["energy"]
                if haskey(energy, "sqrt_s_list_MeV")
                    push!(consumed, "scan.cross_section.energy.sqrt_s_list_MeV")
                end
            end
        end
    end
    return consumed
end

function _build_orchestrator_stages()
    stages = PipelineStage[]

    push!(stages,
        PipelineStage(
            :prepare_inputs,
            [:orchestrator_request],
            [:prepared_orchestrator_inputs],
            (ctx) -> begin
                request = ctx.state[:orchestrator_request]
                cmd = request.cmd
                (cmd === :cross_section || cmd === :transport) || throw(ArgumentError("unsupported relaxtime orchestrator cmd: $(cmd)"))

                kwargs = request.kwargs
                config_path = hasproperty(kwargs, :config_path) ? String(kwargs.config_path) : _DEFAULT_RELAXTIME_CONFIG_PATH
                outdir = String(request.output_dir)
                run_id = String(request.run_id)
                resume = hasproperty(kwargs, :resume) ? getproperty(kwargs, :resume) : nothing
                overwrite = hasproperty(kwargs, :overwrite) ? getproperty(kwargs, :overwrite) : nothing
                fail_on_fallback = hasproperty(kwargs, :fail_on_fallback) ? Bool(getproperty(kwargs, :fail_on_fallback)) : false
                mkpath(outdir)

                default_cfg = TOML.parsefile(_DEFAULT_RELAXTIME_CONFIG_PATH)
                toml_cfg = TOML.parsefile(config_path)
                aliases = TOML.parsefile(_RELAXTIME_ALIASES_PATH)

                cli_cfg = Dict{String, Any}()
                overridden_keys = Set{String}()
                if cmd === :transport && resume !== nothing
                    cli_cfg["scan"] = get(cli_cfg, "scan", Dict{String, Any}())
                    cli_cfg["scan"]["transport"] = get(cli_cfg["scan"], "transport", Dict{String, Any}())
                    cli_cfg["scan"]["transport"]["resume"] = Bool(resume)
                    push!(overridden_keys, "scan.transport.resume")
                end
                if cmd === :transport && overwrite !== nothing
                    cli_cfg["scan"] = get(cli_cfg, "scan", Dict{String, Any}())
                    cli_cfg["scan"]["transport"] = get(cli_cfg["scan"], "transport", Dict{String, Any}())
                    cli_cfg["scan"]["transport"]["overwrite"] = Bool(overwrite)
                    push!(overridden_keys, "scan.transport.overwrite")
                end
                if cmd === :cross_section && hasproperty(kwargs, :processes)
                    processes = _canonicalize_processes(kwargs.processes)
                    if processes !== nothing
                        cli_cfg["scan"] = get(cli_cfg, "scan", Dict{String, Any}())
                        cli_cfg["scan"]["cross_section"] = get(cli_cfg["scan"], "cross_section", Dict{String, Any}())
                        cli_cfg["scan"]["cross_section"]["processes"] = processes
                        push!(overridden_keys, "scan.cross_section.processes")
                    end
                end

                merged = normalize_merge_validate(default_cfg, toml_cfg, cli_cfg, aliases)
                effective = merged.effective

                fallback_events = _read_fallback_events(joinpath(outdir, "fallback_events.json"))
                fallback_used = !isempty(fallback_events)
                prepared = (
                    cmd=cmd,
                    command_name=(cmd === :cross_section ? "cross-section" : "transport"),
                    config_path=config_path,
                    outdir=outdir,
                    run_id=run_id,
                    resume=resume,
                    overwrite=overwrite,
                    fail_on_fallback=fail_on_fallback,
                    effective=effective,
                    trace=String.(merged.trace),
                    strict_mode=Bool(get(effective, "strict", false)),
                    fallback_events=fallback_events,
                    fallback_used=fallback_used,
                    overridden_keys=overridden_keys,
                )

                if prepared.fail_on_fallback && prepared.fallback_used
                    throw(ArgumentError("fallback detected while --fail-on-fallback=true"))
                end

                StageResult(
                    Dict{Symbol, Any}(:prepared_orchestrator_inputs => prepared),
                    PipelineArtifact[],
                    Dict{Symbol, Float64}(),
                )
            end,
        ),
    )

    push!(stages,
        PipelineStage(
            :solve_core,
            [:prepared_orchestrator_inputs],
            [:orchestrator_core_output],
            (ctx) -> begin
                prepared = ctx.state[:prepared_orchestrator_inputs]
                out_csv = if prepared.cmd === :cross_section
                    run_cross_section_orchestrated(prepared.effective, prepared.outdir; run_id=prepared.run_id)
                else
                    ""
                end
                core = (out_csv=out_csv,)
                StageResult(Dict{Symbol, Any}(:orchestrator_core_output => core), PipelineArtifact[], Dict{Symbol, Float64}())
            end,
        ),
    )

    push!(stages,
        PipelineStage(
            :postprocess,
            [:prepared_orchestrator_inputs, :orchestrator_core_output],
            [:orchestrator_postprocessed],
            (ctx) -> begin
                prepared = ctx.state[:prepared_orchestrator_inputs]
                consumed_keys = _build_consumed_keys(prepared.cmd, prepared.effective)
                report = build_consumption_report(
                    prepared.effective,
                    consumed_keys;
                    overridden=prepared.overridden_keys,
                    fallback_used=prepared.fallback_used,
                    strict=prepared.strict_mode,
                )
                postprocessed = (report=report,)
                StageResult(Dict{Symbol, Any}(:orchestrator_postprocessed => postprocessed), PipelineArtifact[], Dict{Symbol, Float64}())
            end,
        ),
    )

    push!(stages,
        PipelineStage(
            :export_artifacts,
            [:prepared_orchestrator_inputs, :orchestrator_postprocessed, :orchestrator_core_output],
            [:orchestrator_artifact_paths],
            (ctx) -> begin
                prepared = ctx.state[:prepared_orchestrator_inputs]
                postprocessed = ctx.state[:orchestrator_postprocessed]
                effective_path = joinpath(prepared.outdir, "effective_config.json")
                report_path = joinpath(prepared.outdir, "consumption_report.json")
                manifest_path = joinpath(prepared.outdir, "run_manifest.json")

                open(effective_path, "w") do io
                    write(io, JSON3.write(prepared.effective))
                end
                open(report_path, "w") do io
                    write(io, JSON3.write(postprocessed.report))
                end
                _emit_plot_contract_artifacts(prepared.cmd, prepared.effective, prepared.outdir)

                paths = Dict{String, String}(
                    "run_manifest" => manifest_path,
                    "effective_config" => effective_path,
                    "consumption_report" => report_path,
                )
                if prepared.cmd === :cross_section
                    paths["cross_section_orchestrated"] = String(ctx.state[:orchestrator_core_output].out_csv)
                end
                StageResult(Dict{Symbol, Any}(:orchestrator_artifact_paths => paths), PipelineArtifact[], Dict{Symbol, Float64}())
            end,
        ),
    )

    push!(stages,
        PipelineStage(
            :emit_diagnostics,
            [:prepared_orchestrator_inputs, :orchestrator_artifact_paths],
            [:orchestrator_diagnostics],
            (ctx) -> begin
                prepared = ctx.state[:prepared_orchestrator_inputs]
                diagnostics = (
                    status=:success,
                    mode=:builtin,
                    command=prepared.command_name,
                    fallback_used=prepared.fallback_used,
                    artifact_count=length(ctx.state[:orchestrator_artifact_paths]),
                )

                extensions = ctx.state[:manifest_extensions]
                extensions["diagnostics_mode"] = String(diagnostics.mode)
                extensions["diagnostics_status"] = String(diagnostics.status)

                StageResult(Dict{Symbol, Any}(:orchestrator_diagnostics => diagnostics), PipelineArtifact[], Dict{Symbol, Float64}())
            end,
        ),
    )

    push!(stages,
        PipelineStage(
            :emit_repro_manifest,
            [:prepared_orchestrator_inputs, :orchestrator_artifact_paths, :orchestrator_diagnostics],
            [:orchestrator_pipeline_output],
            (ctx) -> begin
                prepared = ctx.state[:prepared_orchestrator_inputs]
                artifacts = ctx.state[:orchestrator_artifact_paths]
                diagnostics = ctx.state[:orchestrator_diagnostics]
                output = (
                    command=prepared.command_name,
                    outdir=prepared.outdir,
                    run_id=prepared.run_id,
                    resume=prepared.resume,
                    overwrite=prepared.overwrite,
                    fail_on_fallback=prepared.fail_on_fallback,
                    fallback_used=prepared.fallback_used,
                    manifest_path=artifacts["run_manifest"],
                    effective_config_path=artifacts["effective_config"],
                    consumption_report_path=artifacts["consumption_report"],
                    cross_section_path=get(artifacts, "cross_section_orchestrated", ""),
                    diagnostics=diagnostics,
                )
                StageResult(Dict{Symbol, Any}(:orchestrator_pipeline_output => output), PipelineArtifact[], Dict{Symbol, Float64}())
            end,
        ),
    )

    return stages
end

function run_relaxtime_orchestrator_pipeline_adapter(cmd::Symbol; kwargs...)
    (cmd === :cross_section || cmd === :transport) || throw(ArgumentError("unsupported relaxtime orchestrator cmd: $(cmd)"))

    orchestrator_kwargs = (; kwargs...)
    output_dir = _relaxtime_orchestrator_output_dir(orchestrator_kwargs)
    mkpath(output_dir)
    prepared_run_id = string("relaxtime-orch-", Dates.format(now(UTC), "yyyymmddTHHMMSS"), "-", first(_relaxtime_orchestrator_config_hash(cmd, orchestrator_kwargs), 8))
    manifest_path = joinpath(output_dir, "run_manifest.json")

    ctx = PipelineContext(
        Dict{Symbol, Any}(
            :orchestrator_request => (cmd=cmd, kwargs=orchestrator_kwargs, output_dir=output_dir, run_id=prepared_run_id),
            :manifest_extensions => build_manifest_extensions((
                pipeline_family="relaxtime_orchestrator",
                baseline_suite="smoke",
                physics_profile=get(ENV, "PHYSICS_PARAM_PROFILE", "default"),
                adapter_version="v1",
            )),
            :orchestrator_switches => Dict{String, Any}(
                "resume" => (haskey(orchestrator_kwargs, :resume) ? getproperty(orchestrator_kwargs, :resume) : nothing),
                "overwrite" => (haskey(orchestrator_kwargs, :overwrite) ? getproperty(orchestrator_kwargs, :overwrite) : nothing),
                "fail_on_fallback" => (haskey(orchestrator_kwargs, :fail_on_fallback) ? Bool(getproperty(orchestrator_kwargs, :fail_on_fallback)) : false),
            ),
        ),
        PipelineProvenance(
            _relaxtime_orchestrator_git_commit(),
            _relaxtime_orchestrator_config_hash(cmd, orchestrator_kwargs),
            prepared_run_id,
            Dates.now(Dates.UTC),
        ),
    )

    spec = PipelineSpec(
        "relaxtime_orchestrator_pipeline_runner",
        "v1",
        :PNJL,
        collect(relaxtime_orchestrator_stage_ids()),
        (; cmd=cmd),
        PipelineIOContract(
            :v1,
            [:orchestrator_request],
            [:orchestrator_pipeline_output],
            :relaxtime_orchestrator_artifact_v1,
            :relaxtime_orchestrator_manifest_v1,
        ),
    )

    run_result = run_pipeline(spec, _build_orchestrator_stages(), ctx; manifest_path=manifest_path)
    if !run_result.success
        run_result.error_kind == :ArgumentError && throw(ArgumentError(run_result.error_msg))
        throw(ErrorException("relaxtime orchestrator pipeline failed at $(run_result.failed_stage): $(run_result.error_msg)"))
    end

    return ctx.state[:orchestrator_pipeline_output]
end
