module GapTransportScanProvenance

export build_effective_config
export build_summary
export collect_artifacts
export write_scan_sidecars

function _phase_reference_summary_dict(source)
    summary = Main.PhaseReferenceAdapter.source_summary(source)
    return Dict{String,Any}(String(key) => value for (key, value) in pairs(summary))
end

function build_effective_config(opts; phase_reference=nothing,muB_values=nothing)
    resolved_root = if phase_reference !== nothing
        get(Main.PhaseReferenceAdapter.source_summary(phase_reference), :source_root, opts.phase_reference_root)
    elseif opts.phase_reference_mode === :runtime && isdefined(Main, :DEFAULT_PHASE_REFERENCE_ROOT)
        Main.DEFAULT_PHASE_REFERENCE_ROOT
    else
        opts.phase_reference_root
    end
    cfg = Dict{String,Any}(
        "output" => opts.output,
        "channel_diagnostics_output" => opts.channel_diagnostics_output,
        "failed_points_output" => opts.failed_points_output,
        "xi_values" => opts.xi_values,
        "tmin_mev" => opts.tmin_mev,
        "tmax_mev" => opts.tmax_mev,
        "tstep_mev" => opts.tstep_mev,
        "mubmin_mev" => opts.mubmin_mev,
        "mubmax_mev" => opts.mubmax_mev,
        "mubstep_mev" => opts.mubstep_mev,
        "overwrite" => opts.overwrite,
        "resume" => opts.resume,
        "compute_bulk" => opts.compute_bulk,
        "p_num" => opts.p_num,
        "t_num" => opts.t_num,
        "max_iter" => opts.max_iter,
        "tau_p_nodes" => opts.tau_p_nodes,
        "tau_angle_nodes" => opts.tau_angle_nodes,
        "tau_phi_nodes" => opts.tau_phi_nodes,
        "tau_n_sigma_points" => opts.tau_n_sigma_points,
        "propagator_xi_policy" => String(opts.propagator_xi_policy),
        "sigma_cache_policy" => String(opts.sigma_cache_policy),
        "integration_mode" => String(opts.integration_mode),
        "gc_every_n" => opts.gc_every_n,
        "tr_p_nodes" => opts.tr_p_nodes,
        "tr_p_max_fm" => opts.tr_p_max_fm,
        "phase_reference_root" => resolved_root,
        "phase_reference_layer" => String(opts.phase_reference_layer),
        "phase_reference_mode" => String(opts.phase_reference_mode),
        "phase_reference_source" => phase_reference === nothing ? "none" : String(Main.PhaseReferenceAdapter.source_kind(phase_reference)),
    )
    phase_reference !== nothing && (cfg["phase_reference_summary"] = _phase_reference_summary_dict(phase_reference))
    muB_values===nothing || (cfg["muB_values_MeV"]=muB_values)
    Main.TransportSamplingScanConfig.record_profile!(cfg,opts)
    return cfg
end

function build_summary(stats_success::Integer, stats_error::Integer, stats_skipped::Integer)
    return Dict{String,Any}(
        "points_total" => stats_success + stats_error,
        "success_count" => stats_success,
        "error_count" => stats_error,
        "skipped_count" => stats_skipped,
    )
end

function collect_artifacts(opts)
    artifacts = String[opts.output]
    if hasproperty(opts,:sampling_profile) && opts.sampling_profile !== nothing
        path = Main.TransportSamplingScanConfig.sampling_diagnostics_path(opts)
        isfile(path) && push!(artifacts,path)
        contract = Main.TransportSamplingScanConfig.sampling_contract_path(opts)
        isfile(contract) && push!(artifacts,contract)
    end
    if opts.channel_diagnostics_output !== nothing
        push!(artifacts, opts.channel_diagnostics_output)
    end
    if opts.failed_points_output !== nothing
        push!(artifacts, opts.failed_points_output)
    end
    return artifacts
end

function write_scan_sidecars(provenance_dir::AbstractString, ctx, opts, stats_success::Integer, stats_error::Integer, stats_skipped::Integer; phase_reference=nothing,muB_values=nothing)
    effective_config = build_effective_config(opts; phase_reference=phase_reference,muB_values=muB_values)
    summary = build_summary(stats_success, stats_error, stats_skipped)
    artifacts = collect_artifacts(opts)
    Main.ProvenanceMetadata.write_run_sidecars(
        provenance_dir;
        ctx=ctx,
        effective_config=effective_config,
        artifacts=artifacts,
        summary=summary,
    )
    return (; effective_config, summary, artifacts)
end

end # module GapTransportScanProvenance
