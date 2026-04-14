using Dates
using JSON3
using SHA

@inline function _workflow_repo_root()
    return normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
end

function _trapz(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    n == length(y) || throw(ArgumentError("trapz length mismatch"))
    n <= 1 && return 0.0
    acc = 0.0
    @inbounds for i in 1:(n - 1)
        acc += 0.5 * (y[i] + y[i + 1]) * (x[i + 1] - x[i])
    end
    return acc
end

function _ensure_t190_chain_lib_loaded(repo_root::String)
    lib_path = joinpath(repo_root, "scripts", "analysis", "relaxtime", "t190_sigma_chain_decomposition_lib.jl")
    isfile(lib_path) || throw(ArgumentError("t190 diagnostics library not found: $(lib_path)"))
    if !isdefined(Main, :PROJECT_ROOT)
        Core.eval(Main, :(PROJECT_ROOT = $repo_root))
    end
    required_symbols = (
        :build_state_point,
        :process_threshold_info,
        :decompose_mixed_p_propagator_chain,
    )
    if any(sym -> !isdefined(Main, sym), required_symbols)
        Base.include(Main, lib_path)
    end
    return nothing
end

function _write_t190_mixed_p_chain_artifacts(diagnostics_dir::String)
    repo_root = _workflow_repo_root()
    _ensure_t190_chain_lib_loaded(repo_root)

    mkpath(diagnostics_dir)
    out_csv = joinpath(diagnostics_dir, "t190_mixed_p_chain.csv")
    out_summary = joinpath(diagnostics_dir, "t190_mixed_p_chain_summary.csv")

    T_MeV = 190.0
    muB_MeV = 0.0
    xi_A = -0.10
    xi_B = -0.08
    processes = Symbol[:uubar_to_ddbar, :uubar_to_uubar]
    ds_vals = collect(range(1e-6, 2.0; length=16))

    stA = Base.invokelatest((T, muB, xi) -> getfield(Main, :build_state_point)(T, muB, xi), T_MeV, muB_MeV, xi_A)
    stB = Base.invokelatest((T, muB, xi) -> getfield(Main, :build_state_point)(T, muB, xi), T_MeV, muB_MeV, xi_B)

    summary = Dict{String, Dict{String, Float64}}()

    io = open(out_csv, "w")
    io_sum = open(out_summary, "w")
    try
        println(io, "process,xi_state,s_minus_sth,abs_detM_sq,abs_D_mixed_P_sq,abs_JMJ_sq,detK_plus")
        println(io_sum, "process,ratio_detM_area_B_over_A,ratio_Dmixed_area_B_over_A,ratio_detM_point_B_over_A,ratio_Dmixed_point_B_over_A,ratio_JMJ_area_B_over_A,ratio_detK_area_B_over_A")

        for process in processes
            thA = Base.invokelatest((proc, qp) -> getfield(Main, :process_threshold_info)(proc, qp), process, stA.quark_params)
            thB = Base.invokelatest((proc, qp) -> getfield(Main, :process_threshold_info)(proc, qp), process, stB.quark_params)

            detM_A = Float64[]
            detM_B = Float64[]
            Dmix_A = Float64[]
            Dmix_B = Float64[]
            JMJ_A = Float64[]
            JMJ_B = Float64[]
            detK_A = Float64[]
            detK_B = Float64[]

            detM_point_ratio = NaN
            Dmix_point_ratio = NaN

            for (idx, ds) in enumerate(ds_vals)
                sA = thA.s_th + ds
                sB = thB.s_th + ds
                tbA = Main.TotalCrossSection.calculate_t_bounds(sA, thA.mi, thA.mj, thA.mc, thA.md)
                tbB = Main.TotalCrossSection.calculate_t_bounds(sB, thB.mi, thB.mj, thB.mc, thB.md)
                tA = 0.5 * (tbA.t_min + tbA.t_max)
                tB = 0.5 * (tbB.t_min + tbB.t_max)

                pA = Base.invokelatest(
                    (proc, s, t, qp, tp, kc) -> getfield(Main, :decompose_mixed_p_propagator_chain)(proc, s, t, qp, tp, kc),
                    process,
                    sA,
                    tA,
                    stA.quark_params,
                    stA.thermo_params,
                    stA.K_coeffs,
                )
                pB = Base.invokelatest(
                    (proc, s, t, qp, tp, kc) -> getfield(Main, :decompose_mixed_p_propagator_chain)(proc, s, t, qp, tp, kc),
                    process,
                    sB,
                    tB,
                    stB.quark_params,
                    stB.thermo_params,
                    stB.K_coeffs,
                )

                push!(detM_A, pA.abs_detM_sq)
                push!(detM_B, pB.abs_detM_sq)
                push!(Dmix_A, pA.abs_D_mixed_P_sq)
                push!(Dmix_B, pB.abs_D_mixed_P_sq)
                push!(JMJ_A, pA.abs_JMJ_sq)
                push!(JMJ_B, pB.abs_JMJ_sq)
                push!(detK_A, pA.detK_plus)
                push!(detK_B, pB.detK_plus)

                if idx == 1
                    detM_point_ratio = pB.abs_detM_sq / pA.abs_detM_sq
                    Dmix_point_ratio = pB.abs_D_mixed_P_sq / pA.abs_D_mixed_P_sq
                end

                println(io, join((string(process), string(xi_A), ds, pA.abs_detM_sq, pA.abs_D_mixed_P_sq, pA.abs_JMJ_sq, pA.detK_plus), ','))
                println(io, join((string(process), string(xi_B), ds, pB.abs_detM_sq, pB.abs_D_mixed_P_sq, pB.abs_JMJ_sq, pB.detK_plus), ','))
            end

            ratio_detM = _trapz(ds_vals, detM_B) / _trapz(ds_vals, detM_A)
            ratio_Dmix = _trapz(ds_vals, Dmix_B) / _trapz(ds_vals, Dmix_A)
            ratio_JMJ = _trapz(ds_vals, JMJ_B) / _trapz(ds_vals, JMJ_A)
            ratio_detK = _trapz(ds_vals, detK_B) / _trapz(ds_vals, detK_A)

            println(io_sum, join((string(process), ratio_detM, ratio_Dmix, detM_point_ratio, Dmix_point_ratio, ratio_JMJ, ratio_detK), ','))

            summary[string(process)] = Dict(
                "ratio_detM_area_B_over_A" => ratio_detM,
                "ratio_Dmixed_area_B_over_A" => ratio_Dmix,
                "ratio_detM_point_B_over_A" => detM_point_ratio,
                "ratio_Dmixed_point_B_over_A" => Dmix_point_ratio,
            )
        end
    finally
        close(io)
        close(io_sum)
    end

    return (csv_path=out_csv, summary_path=out_summary, summary=summary)
end

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

    diagnostics_mode = hasproperty(kwargs, :diagnostics_mode) ? Symbol(getproperty(kwargs, :diagnostics_mode)) : :none
    diagnostics_mode in (:none, :t190_chain) || throw(ArgumentError("unsupported diagnostics_mode: $(diagnostics_mode)"))

    diagnostics_strict = hasproperty(kwargs, :diagnostics_strict) ? Bool(getproperty(kwargs, :diagnostics_strict)) : false

    diagnostics_output_dir = if hasproperty(kwargs, :diagnostics_output_dir)
        String(getproperty(kwargs, :diagnostics_output_dir))
    else
        joinpath(output_dir, "diagnostics")
    end

    passthrough_pairs = Pair{Symbol, Any}[]
    for (k, v) in pairs(kwargs)
        if k in (:T_fm, :mu_fm, :xi, :output_dir, :diagnostics_mode, :diagnostics_output_dir, :diagnostics_strict)
            continue
        end
        push!(passthrough_pairs, k => v)
    end

    return (
        T_fm=T_fm,
        mu_fm=mu_fm,
        xi=xi,
        output_dir=output_dir,
        diagnostics_mode=diagnostics_mode,
        diagnostics_output_dir=diagnostics_output_dir,
        diagnostics_strict=diagnostics_strict,
        passthrough=(; passthrough_pairs...),
    )
end

function _emit_workflow_diagnostics_index(prepared, output)::NamedTuple
    mode = prepared.diagnostics_mode
    mode === :none && return (mode=:none, status=:disabled, index_path="", artifacts=Dict{String, String}())

    mkpath(prepared.diagnostics_output_dir)
    index_path = joinpath(prepared.diagnostics_output_dir, "diagnostics_index.json")

    run_context = Dict(
        "T_fm" => prepared.T_fm,
        "mu_fm" => prepared.mu_fm,
        "xi" => prepared.xi,
        "run_manifest" => output.run_manifest,
    )

    artifacts = Dict{String, String}(
        "run_manifest" => output.run_manifest,
    )
    status = :success
    error_payload = nothing
    chain_summary = Dict{String, Any}()

    if mode === :t190_chain
        try
            chain_artifacts = _write_t190_mixed_p_chain_artifacts(prepared.diagnostics_output_dir)
            artifacts["t190_mixed_p_chain_csv"] = chain_artifacts.csv_path
            artifacts["t190_mixed_p_chain_summary_csv"] = chain_artifacts.summary_path
            chain_summary = chain_artifacts.summary
        catch err
            status = err isa ArgumentError ? :unavailable : :failed
            error_payload = Dict(
                "type" => string(typeof(err)),
                "message" => sprint(showerror, err),
            )
        end
    end

    eta_val = isfinite(output.transport.eta) ? Float64(output.transport.eta) : nothing
    sigma_val = isfinite(output.transport.sigma) ? Float64(output.transport.sigma) : nothing
    zeta_val = isfinite(output.transport.zeta) ? Float64(output.transport.zeta) : nothing

    payload = Dict(
        "mode" => String(mode),
        "status" => String(status),
        "run_context" => run_context,
        "artifacts" => artifacts,
        "summary" => Dict(
            "eta" => eta_val,
            "sigma" => sigma_val,
            "zeta" => zeta_val,
        ),
        "error" => error_payload,
        "t190_chain_summary" => chain_summary,
    )

    open(index_path, "w") do io
        write(io, JSON3.write(payload))
    end

    t190_summary = if isempty(chain_summary)
        Dict{String, Float64}()
    else
        get(chain_summary, "uubar_to_ddbar", Dict{String, Float64}())
    end

    return (
        mode=mode,
        status=status,
        index_path=index_path,
        artifacts=artifacts,
        t190_summary=t190_summary,
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
            :emit_diagnostics,
            [:prepared_inputs, :artifact_paths],
            [:workflow_diagnostics],
            (ctx) -> begin
                prepared = ctx.state[:prepared_inputs]
                run_manifest = ctx.state[:artifact_paths]["run_manifest"]
                diagnostics = _emit_workflow_diagnostics_index(prepared, (transport=ctx.state[:transport], run_manifest=run_manifest))

                extensions = ctx.state[:manifest_extensions]
                extensions["diagnostics_mode"] = String(prepared.diagnostics_mode)
                extensions["diagnostics_status"] = String(diagnostics.status)
                extensions["diagnostics_index_path"] = String(diagnostics.index_path)
                if prepared.diagnostics_mode == :t190_chain && !isempty(diagnostics.t190_summary)
                    extensions["diagnostics_t190_ratio_detM_area_B_over_A"] = diagnostics.t190_summary["ratio_detM_area_B_over_A"]
                    extensions["diagnostics_t190_ratio_Dmixed_area_B_over_A"] = diagnostics.t190_summary["ratio_Dmixed_area_B_over_A"]
                    extensions["diagnostics_t190_ratio_detM_point_B_over_A"] = diagnostics.t190_summary["ratio_detM_point_B_over_A"]
                    extensions["diagnostics_t190_ratio_Dmixed_point_B_over_A"] = diagnostics.t190_summary["ratio_Dmixed_point_B_over_A"]
                end

                if prepared.diagnostics_mode != :none && prepared.diagnostics_strict && diagnostics.status != :success
                    throw(ArgumentError("workflow diagnostics failed in strict mode"))
                end

                StageResult(
                    Dict{Symbol, Any}(:workflow_diagnostics => diagnostics),
                    PipelineArtifact[],
                    Dict{Symbol, Float64}(),
                )
            end,
        ),
    )

    push!(stages,
        PipelineStage(
            :emit_repro_manifest,
            [:transport, :artifact_paths, :workflow_diagnostics],
            [:workflow_pipeline_output],
            (ctx) -> begin
                run_manifest = ctx.state[:artifact_paths]["run_manifest"]
                output_with_diag = (
                    transport=ctx.state[:transport],
                    run_manifest=run_manifest,
                    diagnostics=ctx.state[:workflow_diagnostics],
                )
                StageResult(
                    Dict{Symbol, Any}(:workflow_pipeline_output => output_with_diag),
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
        collect(workflow_pipeline_stage_ids()),
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
