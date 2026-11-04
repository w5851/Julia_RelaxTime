module GapTransportScanOrchestration

@inline function _csv_quote(text::AbstractString)
    return "\"" * replace(String(text), "\"" => "\"\"") * "\""
end

function build_scan_runtime(opts)
    profile = Main.TransportSamplingScanConfig.profile_for_options(opts)
    profile === nothing || return Main.Models.TransportSamplingProfiles.build_runtime(profile)
    p_grid, p_w, sigma_cutoff = Main.integration_grids(opts)
    cos_grid, cos_w = Main.gauleg(-1.0, 1.0, opts.tau_angle_nodes)
    phi_grid, phi_w = Main.gauleg(0.0, 2 * pi, opts.tau_phi_nodes)
    transport_config = Main.TransportWorkflow.TransportIntegrationConfig(
        p_nodes=opts.tr_p_nodes,
        p_max=opts.tr_p_max_fm,
    )
    return (
        profile=nothing,
        model=Main.PNJL_MODEL,
        p_grid=p_grid,
        p_w=p_w,
        sigma_cutoff=sigma_cutoff,
        cos_grid=cos_grid,
        cos_w=cos_w,
        phi_grid=phi_grid,
        phi_w=phi_w,
        tau_kwargs=(
            p_nodes=opts.tau_p_nodes,
            angle_nodes=opts.tau_angle_nodes,
            phi_nodes=opts.tau_phi_nodes,
            n_sigma_points=opts.tau_n_sigma_points,
            p_grid=p_grid,
            p_w=p_w,
            cos_grid=cos_grid,
            cos_w=cos_w,
            phi_grid=phi_grid,
            phi_w=phi_w,
            sigma_cutoff=sigma_cutoff,
            propagator_xi_policy=opts.propagator_xi_policy,
            sigma_cache_policy=opts.sigma_cache_policy,
        ),
        transport_config=transport_config,
    )
end

function _solve_transport(eq, T_fm::Float64, muq_fm::Float64, xi::Float64, opts, runtime, K_coeffs; run_id="")
    if runtime.profile !== nothing
        callback = row -> _write_sampling_diagnostic(opts,T_fm,muq_fm,xi,row;run_id=run_id)
        result = Main.Models.TransportSamplingProfiles.solve_profile_transport(runtime,eq,T_fm,muq_fm;
            xi=xi,compute_bulk=opts.compute_bulk,propagator_xi_policy=opts.propagator_xi_policy,on_channel=callback)
        return result, opts.compute_bulk
    end
    result = Main.Models.solve_transport_from_equilibrium(
            eq,
            T_fm,
            muq_fm;
            xi=xi,
            compute_tau=true,
            K_coeffs=K_coeffs,
            compute_bulk=opts.compute_bulk,
            p_num=opts.p_num,
            t_num=opts.t_num,
            tau_kwargs=runtime.tau_kwargs,
            transport_config=runtime.transport_config,
        )
    return result, opts.compute_bulk
end

function _build_scan_row(eq, diag, res, T_mev::Float64, muB_mev::Float64, xi::Float64, opts, ctx, point_meta;
    model=Main.PNJL_MODEL)
    T_fm = T_mev / Main.ħc_MeV_fm
    muq_mev = muB_mev / 3.0
    muq_fm = muq_mev / Main.ħc_MeV_fm

    Φ = Float64(eq.x_state[4])
    Φbar = Float64(eq.x_state[5])
    masses = (u=Float64(eq.masses[1]), d=Float64(eq.masses[2]), s=Float64(eq.masses[3]))

    dens = res.densities
    tau = res.tau
    tauinv = res.tau_inv
    tr = res.transport
    rates = res.rates

    P_fm4inv, _, s_fm3inv, epsilon_fm4inv = Main.Models.model_thermo(
        model,
        eq.x_state,
        eq.mu_vec,
        T_fm,
        p_num=opts.p_num,
        t_num=opts.t_num,
        xi=xi,
    )

    rho_quark_net = (dens.u - dens.ubar) + (dens.d - dens.dbar) + (dens.s - dens.sbar)
    rho_baryon = rho_quark_net / 3.0
    rho_norm = rho_baryon / Main.ρ0_inv_fm3

    omega_fm4inv = Main.Models.omega(
        model,
        eq.x_state,
        T_fm,
        eq.mu_vec;
        p_num=opts.p_num,
        t_num=opts.t_num,
        xi=xi,
    )
    omega_MeV_fm3 = omega_fm4inv * Main.ħc_MeV_fm
    P_MeV_fm3 = P_fm4inv * Main.ħc_MeV_fm
    epsilon_MeV_fm3 = epsilon_fm4inv * Main.ħc_MeV_fm
    eps_minus_3P_over_T4 = (isfinite(epsilon_fm4inv) && isfinite(P_fm4inv) && isfinite(T_fm) && T_fm != 0.0) ? ((epsilon_fm4inv - 3.0 * P_fm4inv) / T_fm^4) : NaN
    eta_over_s = (isfinite(tr.eta) && isfinite(s_fm3inv) && s_fm3inv != 0.0) ? (tr.eta / s_fm3inv) : NaN
    zeta_over_s = (isfinite(tr.zeta) && isfinite(s_fm3inv) && s_fm3inv != 0.0) ? (tr.zeta / s_fm3inv) : NaN
    sigma_over_T = (isfinite(tr.sigma) && isfinite(T_fm) && T_fm != 0.0) ? (tr.sigma / T_fm) : NaN
    sigma_over_T_over_eta_over_s = (isfinite(sigma_over_T) && isfinite(eta_over_s) && eta_over_s != 0.0) ? (sigma_over_T / eta_over_s) : NaN
    zeta_over_s_over_eta_over_s = (isfinite(zeta_over_s) && isfinite(eta_over_s) && eta_over_s != 0.0) ? (zeta_over_s / eta_over_s) : NaN
    quality_flag, quality_reason, quality_metric = Main.assess_point_quality(tau, eta_over_s, sigma_over_T)

    row = join([
        string(T_mev), string(muq_mev), string(muB_mev), string(xi),
        string(get(point_meta, :mode, :grid)),
        string(get(point_meta, :phase_reference_kind, :regular_grid)),
        string(get(point_meta, :scan_group, "")),
        _csv_quote(string(get(point_meta, :group_label, ""))),
        string(get(point_meta, :plot_panel, "")),
        _csv_quote(string(get(point_meta, :plot_panel_label, ""))),
        string(get(point_meta, :plot_series, "")),
        _csv_quote(string(get(point_meta, :plot_series_label, ""))),
        string(get(point_meta, :T_phase_base_MeV, NaN)),
        string(get(point_meta, :alpha_T, NaN)),
        string(T_fm), string(muq_fm),
        Main.csv_bool(eq.converged), string(eq.iterations), string(eq.residual_norm), string(eq.solver_backend), diag.seed_source, string(diag.phase_prev), string(diag.phase_curr), string(diag.phase_structure), string(diag.phase_boundary_xi_used),
        string(Φ), string(Φbar),
        string(masses.u), string(masses.d), string(masses.s),
        string(rho_baryon), string(rho_norm),
        string(omega_fm4inv), string(P_fm4inv), string(epsilon_fm4inv), string(s_fm3inv),
        string(omega_MeV_fm3), string(P_MeV_fm3), string(epsilon_MeV_fm3),
        string(eps_minus_3P_over_T4),
        string(dens.u), string(dens.d), string(dens.s), string(dens.ubar), string(dens.dbar), string(dens.sbar),
        string(tau.u), string(tau.d), string(tau.s), string(tau.ubar), string(tau.dbar), string(tau.sbar),
        string(tauinv.u), string(tauinv.d), string(tauinv.s), string(tauinv.ubar), string(tauinv.dbar), string(tauinv.sbar),
        string(tr.eta), string(tr.sigma), string(tr.zeta), string(eta_over_s), string(zeta_over_s),
        string(sigma_over_T), string(sigma_over_T_over_eta_over_s), string(zeta_over_s_over_eta_over_s),
        Main.csv_bool(quality_flag), quality_reason, string(quality_metric), string(ctx.run_id),
    ], ',')

    return (
        row=row,
        muq_mev=muq_mev,
        densities=dens,
        tau_inv=tauinv,
        rates=rates,
    )
end

function execute_gap_transport_scan_point!(io, channel_io, failed_io,
    T_mev::Float64, muB_mev::Float64, xi::Float64, opts, ctx, runtime;
    previous_solution=nothing,
    previous_phase::Symbol=:unknown,
    phase_reference=nothing,
    phase_reference_mode::Symbol=:runtime,
    point_meta=(; mode=:grid, phase_reference_kind=:regular_grid, scan_group="", group_label="", plot_panel="", plot_panel_label="", plot_series="", plot_series_label="", T_phase_base_MeV=NaN, alpha_T=NaN),
)
    next_solution = previous_solution
    next_phase = previous_phase
    diag_hint = (seed_source="unknown", phase_prev=previous_phase, phase_curr=:unknown)

    try
        eq, next_solution, next_phase, diag = Main.solve_equilibrium_with_diagnostics(
            T_mev,
            muB_mev,
            xi,
            opts;
            previous_solution=previous_solution,
            previous_phase=previous_phase,
            phase_reference=phase_reference,
            phase_reference_mode=phase_reference_mode,
        )
        diag_hint = diag

        T_fm = T_mev / Main.ħc_MeV_fm
        muq_mev = muB_mev / 3.0
        muq_fm = muq_mev / Main.ħc_MeV_fm
        masses = (u=Float64(eq.masses[1]), d=Float64(eq.masses[2]), s=Float64(eq.masses[3]))
        ktmp = runtime.profile === nothing ?
            Main.build_K_data(T_fm, muq_fm, masses, Float64(eq.x_state[4]), Float64(eq.x_state[5])) :
            (K_coeffs=nothing,)
        res, _ = _solve_transport(eq, T_fm, muq_fm, xi, opts, runtime, ktmp.K_coeffs;run_id=ctx.run_id)
        payload = _build_scan_row(eq, diag, res, T_mev, muB_mev, xi, opts, ctx, point_meta;model=runtime.model)

        println(io, payload.row)
        flush(io)

        if channel_io !== nothing && payload.rates !== nothing
            Main.write_channel_diagnostics_rows!(
                channel_io,
                T_mev,
                payload.muq_mev,
                muB_mev,
                xi,
                payload.densities,
                payload.rates,
                payload.tau_inv,
                eq.solver_backend,
                diag,
            )
        end

        return (
            success=true,
            next_solution=next_solution,
            next_phase=next_phase,
        )
    catch point_err
        point_err isa InterruptException && rethrow()
        @warn "SCAN POINT FAILED — skipped" T_mev=T_mev muB_mev=muB_mev xi=xi err=point_err
        if failed_io !== nothing
            Main.write_failed_point_row!(failed_io, T_mev, muB_mev, xi, diag_hint, point_err)
        end
        return (
            success=false,
            next_solution=next_solution,
            next_phase=next_phase,
        )
    end
end

export build_scan_runtime
export execute_gap_transport_scan_point!

function _write_sampling_diagnostic(opts,T_fm,muq_fm,xi,row;run_id="")
    path = Main.TransportSamplingScanConfig.sampling_diagnostics_path(opts)
    Main.ensure_parent_dir(path)
    needs_header = !isfile(path) || filesize(path)==0
    open(path,"a") do io
        needs_header && println(io,"T_MeV,muB_MeV,xi,channel,seconds,rate,estimated_error,nodes,evaluations,events,converged,reason,guard_phase,feature_count,run_id")
        println(io,join((T_fm*Main.ħc_MeV_fm,3muq_fm*Main.ħc_MeV_fm,xi,row.channel,row.seconds,
            row.rate,row.estimated_error,row.nodes,row.evaluations,row.events,true,row.reason,row.guard_phase,row.feature_count,run_id),','))
    end
    return nothing
end

end # module GapTransportScanOrchestration
