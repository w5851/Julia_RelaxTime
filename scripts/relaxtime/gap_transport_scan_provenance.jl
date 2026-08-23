module GapTransportScanProvenance

export build_effective_config
export build_summary
export collect_artifacts
export write_scan_sidecars

function build_effective_config(opts)
    return Dict{String,Any}(
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
        "tau_threshold_subtraction" => opts.tau_threshold_subtraction,
        "tau_asym_window" => opts.tau_asym_window,
        "tau_asym_fit_min_points" => opts.tau_asym_fit_min_points,
        "tau_asym_extra_points" => opts.tau_asym_extra_points,
        "tau_interpolation_mode" => String(opts.tau_interpolation_mode),
        "propagator_xi_policy" => String(opts.propagator_xi_policy),
        "sigma_cache_policy" => String(opts.sigma_cache_policy),
        "sigma_grid_n" => opts.sigma_grid_n,
        "integration_mode" => String(opts.integration_mode),
        "gc_every_n" => opts.gc_every_n,
        "tr_p_nodes" => opts.tr_p_nodes,
        "tr_p_max_fm" => opts.tr_p_max_fm,
        "phase_reference_root" => opts.phase_reference_root,
        "phase_reference_layer" => String(opts.phase_reference_layer),
        "phase_reference_mode" => String(opts.phase_reference_mode),
    )
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
    if opts.channel_diagnostics_output !== nothing
        push!(artifacts, opts.channel_diagnostics_output)
    end
    if opts.failed_points_output !== nothing
        push!(artifacts, opts.failed_points_output)
    end
    return artifacts
end

function write_scan_sidecars(provenance_dir::AbstractString, ctx, opts, stats_success::Integer, stats_error::Integer, stats_skipped::Integer)
    effective_config = build_effective_config(opts)
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
