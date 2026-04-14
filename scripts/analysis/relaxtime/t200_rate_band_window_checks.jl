#!/usr/bin/env julia

using CSV

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(@__DIR__, "t190_sigma_chain_decomposition_lib.jl"))

const BAND_EDGES = [0.0, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, Inf]

function ensure_parent_dir(path::AbstractString)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
end

function read_scan_metadata(path::AbstractString)
    meta = Dict{String, String}()
    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            startswith(s, "#") || break
            body = strip(s[2:end])
            isempty(body) && continue
            k_v = split(body, ':'; limit=2)
            length(k_v) == 2 || continue
            key = strip(k_v[1])
            val = strip(k_v[2])
            meta[key] = val
        end
    end
    return meta
end

function tau_settings_from_scan_metadata(path::AbstractString)
    meta = read_scan_metadata(path)
    p_nodes = parse(Int, get(meta, "tau_p_nodes", "20"))
    angle_nodes = parse(Int, get(meta, "tau_angle_nodes", "4"))
    phi_nodes = parse(Int, get(meta, "tau_phi_nodes", "8"))
    n_sigma_points = parse(Int, get(meta, "tau_n_sigma_points", "6"))
    interpolation_mode = Symbol(get(meta, "tau_interpolation_mode", "linear"))
    integration_mode = Symbol(get(meta, "integration_mode", "finite_15"))
    return (
        p_nodes=p_nodes,
        angle_nodes=angle_nodes,
        phi_nodes=phi_nodes,
        n_sigma_points=n_sigma_points,
        interpolation_mode=interpolation_mode,
        integration_mode=integration_mode,
    )
end

function build_state_from_main_row(row)
    T_fm = Float64(row.T_fm)
    muq_fm = Float64(row.muq_fm)
    xi = Float64(row.xi)
    masses = (u=Float64(row.m_u), d=Float64(row.m_d), s=Float64(row.m_s))
    Phi = Float64(row.Phi)
    Phibar = Float64(row.Phibar)

    qp_basic = (m=masses, μ=(u=muq_fm, d=muq_fm, s=muq_fm))
    thermo_params = (T=T_fm, Φ=Phi, Φbar=Phibar, ξ=xi)
    A_vals = Main.AFieldBuilder.build_A_triplet(qp_basic, thermo_params)

    A_u_iso = Main.A(masses.u, muq_fm, T_fm, Phi, Phibar, Main.DEFAULT_MOMENTUM_NODES_LIB, Main.DEFAULT_MOMENTUM_WEIGHTS_LIB)
    A_s_iso = Main.A(masses.s, muq_fm, T_fm, Phi, Phibar, Main.DEFAULT_MOMENTUM_NODES_LIB, Main.DEFAULT_MOMENTUM_WEIGHTS_LIB)
    G_u = Main.calculate_G_from_A(A_u_iso, masses.u)
    G_s = Main.calculate_G_from_A(A_s_iso, masses.s)
    K_coeffs = Main.calculate_effective_couplings(Main.G_fm2, Main.K_fm5, G_u, G_s)

    quark_params = (m=masses, μ=(u=muq_fm, d=muq_fm, s=muq_fm), A=(u=A_vals.u, d=A_vals.d, s=A_vals.s))
    return (quark_params=quark_params, thermo_params=thermo_params, K_coeffs=K_coeffs)
end

function compute_rate_band_stats(process::Symbol, state;
    p_nodes::Int=20, angle_nodes::Int=4, phi_nodes::Int=8,
    n_sigma_points::Int=6, interpolation_mode::Symbol=:linear,
    integration_mode::Symbol=:finite_15)

    p_grid, p_w = if integration_mode == :finite_15
        Main.GaussLegendre.gauleg(0.0, 15.0, p_nodes)
    elseif integration_mode == :finite_lambda
        Main.GaussLegendre.gauleg(0.0, Main.Constants_PNJL.Λ_inv_fm, p_nodes)
    else
        (nothing, nothing)
    end
    cos_grid, cos_w = Main.GaussLegendre.gauleg(-1.0, 1.0, angle_nodes)
    phi_grid, phi_w = Main.GaussLegendre.gauleg(0.0, 2 * pi, phi_nodes)

    band_omega_ref = Ref(Vector{Float64}())
    band_omega_sigma_ref = Ref(Vector{Float64}())

    rate = Main.AverageScatteringRate.average_scattering_rate(
        process,
        state.quark_params,
        state.thermo_params,
        state.K_coeffs;
        p_nodes=p_nodes,
        angle_nodes=angle_nodes,
        phi_nodes=phi_nodes,
        p_grid=p_grid,
        p_w=p_w,
        cos_grid=cos_grid,
        cos_w=cos_w,
        phi_grid=phi_grid,
        phi_w=phi_w,
        n_sigma_points=n_sigma_points,
        sigma_cutoff=Main.Constants_PNJL.Λ_inv_fm,
        threshold_subtraction=true,
        asym_window=0.6,
        asym_fit_min_points=8,
        asym_extra_points=10,
        interpolation_mode=interpolation_mode,
        band_edges=BAND_EDGES,
        band_omega_out=band_omega_ref,
        band_omega_sigma_out=band_omega_sigma_ref,
    )

    pi_sym, pj_sym, pc_sym, pd_sym = Main.TotalCrossSection.parse_particles_from_process(process)
    mi = Main.AverageScatteringRate.get_mass(pi_sym, state.quark_params)
    mj = Main.AverageScatteringRate.get_mass(pj_sym, state.quark_params)
    mc = Main.AverageScatteringRate.get_mass(pc_sym, state.quark_params)
    md = Main.AverageScatteringRate.get_mass(pd_sym, state.quark_params)
    μi = Main.AverageScatteringRate.get_mu(pi_sym, state.quark_params)
    μj = Main.AverageScatteringRate.get_mu(pj_sym, state.quark_params)

    ρ_i = Main.AverageScatteringRate.number_density(
        pi_sym,
        mi,
        μi,
        state.thermo_params.T,
        state.thermo_params.Φ,
        state.thermo_params.Φbar,
        state.thermo_params.ξ;
        p_grid=p_grid,
        p_w=p_w,
        angle_nodes=angle_nodes,
    )
    ρ_j = Main.AverageScatteringRate.number_density(
        pj_sym,
        mj,
        μj,
        state.thermo_params.T,
        state.thermo_params.Φ,
        state.thermo_params.Φbar,
        state.thermo_params.ξ;
        p_grid=p_grid,
        p_w=p_w,
        angle_nodes=angle_nodes,
    )

    prefactor = (Main.AverageScatteringRate.DQ^2) / (32.0 * pi^5 * ρ_i * ρ_j)

    return (
        rate=rate,
        prefactor=prefactor,
        rho_i=ρ_i,
        rho_j=ρ_j,
        omega_by_bin=band_omega_ref[],
        omega_sigma_by_bin=band_omega_sigma_ref[],
        omega_total=sum(band_omega_sigma_ref[]),
    )
end

function main()
    scenarios = [
        (
            name="tauu_neg_udbar",
            process=:udbar_to_udbar,
            species="u",
            xis=[-0.34, -0.32, -0.30],
            main_csv=raw"D:\Desktop\Temp\relaxtime_t200_window\t200_dual_window_main.csv",
            channel_csv=raw"D:\Desktop\Temp\relaxtime_t200_window\t200_dual_window_channel_diag.csv",
        ),
        (
            name="taus_usbar",
            process=:usbar_to_usbar,
            species="s",
            xis=[-0.22, -0.20, -0.18],
            main_csv=raw"D:\Desktop\Temp\relaxtime_t200_window\t200_taus_window_main.csv",
            channel_csv=raw"D:\Desktop\Temp\relaxtime_t200_window\t200_taus_window_channel_diag.csv",
        ),
        (
            name="tauu_pos_uubarddbar",
            process=:uubar_to_ddbar,
            species="u",
            xis=[0.34, 0.36, 0.38],
            main_csv=raw"D:\Desktop\Temp\relaxtime_t200_window\t200_dual_window_main.csv",
            channel_csv=raw"D:\Desktop\Temp\relaxtime_t200_window\t200_dual_window_channel_diag.csv",
        ),
        (
            name="tauu_pos_uubaruubar",
            process=:uubar_to_uubar,
            species="u",
            xis=[0.34, 0.36, 0.38],
            main_csv=raw"D:\Desktop\Temp\relaxtime_t200_window\t200_dual_window_main.csv",
            channel_csv=raw"D:\Desktop\Temp\relaxtime_t200_window\t200_dual_window_channel_diag.csv",
        ),
    ]

    out_detail = raw"D:\Desktop\Temp\relaxtime_t200_window\t200_rate_band_window_checks_detail.csv"
    out_pair = raw"D:\Desktop\Temp\relaxtime_t200_window\t200_rate_band_window_checks_pair_delta.csv"
    ensure_parent_dir(out_detail)
    ensure_parent_dir(out_pair)

    io_d = open(out_detail, "w")
    io_p = open(out_pair, "w")
    try
        println(io_d, "scenario,process,xi,rate_func,rate_direct,contrib_direct,bin,ds_left,ds_right,omega_bin,omega_sigma_bin,omega_sigma_frac,sigma_eff_bin")
        println(io_p, "scenario,process,xi_left,xi_right,delta_rate_func,delta_rate_direct,delta_contrib_direct,bin,delta_omega_sigma_bin,delta_rate_bin")

        cache_stats = Dict{Tuple{String, Float64}, Any}()

        for sc in scenarios
            main_rows = collect(CSV.File(sc.main_csv; comment="#"))
            channel_rows = collect(CSV.File(sc.channel_csv; comment="#"))

            for xi in sc.xis
                row = only(filter(r -> Float64(r.T_MeV) == 200.0 && Float64(r.muB_MeV) == 0.0 && Float64(r.xi) == xi, main_rows))
                st = build_state_from_main_row(row)
                tau_cfg = tau_settings_from_scan_metadata(sc.main_csv)
                stats = compute_rate_band_stats(
                    sc.process,
                    st;
                    p_nodes=tau_cfg.p_nodes,
                    angle_nodes=tau_cfg.angle_nodes,
                    phi_nodes=tau_cfg.phi_nodes,
                    n_sigma_points=tau_cfg.n_sigma_points,
                    interpolation_mode=tau_cfg.interpolation_mode,
                    integration_mode=tau_cfg.integration_mode,
                )
                cache_stats[(sc.name, xi)] = stats

                ch = filter(r -> Float64(r.T_MeV) == 200.0 && Float64(r.muB_MeV) == 0.0 && Float64(r.xi) == xi && String(r.species) == sc.species && String(r.channel) == string(sc.process), channel_rows)
                rate_direct = isempty(ch) ? NaN : Float64(ch[1].rate)
                contrib_direct = isempty(ch) ? NaN : Float64(ch[1].contribution)

                total = stats.omega_total
                for b in 1:(length(BAND_EDGES) - 1)
                    omega_bin = stats.omega_by_bin[b]
                    omega_sigma_bin = stats.omega_sigma_by_bin[b]
                    frac = total == 0.0 ? NaN : (omega_sigma_bin / total)
                    sigma_eff = omega_bin == 0.0 ? NaN : (omega_sigma_bin / omega_bin)
                    println(io_d, join((
                        sc.name,
                        sc.process,
                        xi,
                        stats.rate,
                        rate_direct,
                        contrib_direct,
                        b,
                        BAND_EDGES[b],
                        BAND_EDGES[b + 1],
                        omega_bin,
                        omega_sigma_bin,
                        frac,
                        sigma_eff,
                    ), ','))
                end
            end

            for i in 1:(length(sc.xis) - 1)
                xl = sc.xis[i]
                xr = sc.xis[i + 1]
                sl = cache_stats[(sc.name, xl)]
                sr = cache_stats[(sc.name, xr)]

                ch_l = filter(r -> Float64(r.T_MeV) == 200.0 && Float64(r.muB_MeV) == 0.0 && Float64(r.xi) == xl && String(r.species) == sc.species && String(r.channel) == string(sc.process), channel_rows)
                ch_r = filter(r -> Float64(r.T_MeV) == 200.0 && Float64(r.muB_MeV) == 0.0 && Float64(r.xi) == xr && String(r.species) == sc.species && String(r.channel) == string(sc.process), channel_rows)
                rate_direct_l = isempty(ch_l) ? NaN : Float64(ch_l[1].rate)
                rate_direct_r = isempty(ch_r) ? NaN : Float64(ch_r[1].rate)
                contrib_direct_l = isempty(ch_l) ? NaN : Float64(ch_l[1].contribution)
                contrib_direct_r = isempty(ch_r) ? NaN : Float64(ch_r[1].contribution)

                for b in 1:(length(BAND_EDGES) - 1)
                    dωσ = sr.omega_sigma_by_bin[b] - sl.omega_sigma_by_bin[b]
                    dr_bin = sr.prefactor * sr.omega_sigma_by_bin[b] - sl.prefactor * sl.omega_sigma_by_bin[b]
                    println(io_p, join((
                        sc.name,
                        sc.process,
                        xl,
                        xr,
                        sr.rate - sl.rate,
                        rate_direct_r - rate_direct_l,
                        contrib_direct_r - contrib_direct_l,
                        b,
                        dωσ,
                        dr_bin,
                    ), ','))
                end
            end
        end
    finally
        close(io_d)
        close(io_p)
    end
end

main()
