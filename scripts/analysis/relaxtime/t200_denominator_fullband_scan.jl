#!/usr/bin/env julia

using CSV

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(@__DIR__, "t190_sigma_chain_decomposition_lib.jl"))

const BAND_EDGES = [0.0, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, Inf]
const SAMPLE_DS_POINTS = [1.0e-3, 5.0e-3, 1.0e-2, 2.0e-2, 5.0e-2, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0]

function ensure_parent_dir(path::AbstractString)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
end

function trapz(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    n == length(y) || error("trapz length mismatch")
    n <= 1 && return 0.0
    acc = 0.0
    @inbounds for i in 1:(n - 1)
        acc += 0.5 * (y[i] + y[i + 1]) * (x[i + 1] - x[i])
    end
    return acc
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

function p_grid_from_mode(integration_mode::Symbol, p_nodes::Int)
    if integration_mode == :finite_15
        return Main.GaussLegendre.gauleg(0.0, 15.0, p_nodes)
    elseif integration_mode == :finite_lambda
        return Main.GaussLegendre.gauleg(0.0, Main.Constants_PNJL.Λ_inv_fm, p_nodes)
    end
    return (nothing, nothing)
end

@inline function meson_K(meson::Symbol, K_coeffs)
    if meson === :pi
        return K_coeffs.K123_plus
    elseif meson === :K
        return K_coeffs.K4567_plus
    elseif meson === :sigma_pi
        return K_coeffs.K123_minus
    elseif meson === :sigma_K
        return K_coeffs.K4567_minus
    end
    error("unsupported meson: $meson")
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

function simple_p_branch_for_channel(process::Symbol, channel::Symbol, k0::Float64, k_norm::Float64,
    quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple)
    process_map = Main.Constants_PNJL.SCATTERING_MESON_MAP[process]
    channel_info = process_map[:channels][channel]
    simple_mesons = channel_info[:simple]
    T1, T2 = Main.TotalPropagator.get_flavor_factors_for_channel(process, channel)

    D_simple = 0.0 + 0.0im
    den_ref = NaN + NaN * im

    for meson in simple_mesons
        (meson === :pi || meson === :K) || continue
        pol = Main.TotalPropagator.get_polarization_params(meson, quark_params)
        Π_re, Π_im = Main.PolarizationCache.polarization_aniso_cached(
            pol.channel, k0, k_norm,
            pol.m1, pol.m2,
            pol.μ1, pol.μ2,
            thermo_params.T, thermo_params.Φ, thermo_params.Φbar, thermo_params.ξ,
            pol.A1, pol.A2,
            pol.num_s_quark,
        )
        Π = ComplexF64(Π_re, Π_im)
        K = meson_K(meson, K_coeffs)
        den = 1.0 - 4.0 * K * Π
        D = Main.MesonPropagator.meson_propagator_simple(meson, K_coeffs, Π)
        D_simple += T1 * D * T2
        if isnan(real(den_ref))
            den_ref = den
        end
    end

    return (D_simple=D_simple, den_ref=den_ref)
end

function mixed_p_branch_for_channel(process::Symbol, channel::Symbol, k0::Float64, k_norm::Float64,
    quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple)
    process_map = Main.Constants_PNJL.SCATTERING_MESON_MAP[process]
    channel_info = process_map[:channels][channel]
    if !Bool(channel_info[:mixed_P])
        return (has_mixed=false, D_mixed=0.0 + 0.0im, detM=NaN + NaN * im)
    end

    D_mixed = Main.TotalPropagator.total_propagator_mixed(
        process, channel, :P, k0, k_norm, quark_params, thermo_params, K_coeffs,
    )

    Π_uu_re, Π_uu_im = Main.PolarizationCache.polarization_aniso_cached(
        :P, k0, k_norm,
        quark_params.m.u, quark_params.m.u,
        quark_params.μ.u, quark_params.μ.u,
        thermo_params.T, thermo_params.Φ, thermo_params.Φbar, thermo_params.ξ,
        quark_params.A.u, quark_params.A.u,
        0,
    )
    Π_ss_re, Π_ss_im = Main.PolarizationCache.polarization_aniso_cached(
        :P, k0, k_norm,
        quark_params.m.s, quark_params.m.s,
        quark_params.μ.s, quark_params.μ.s,
        thermo_params.T, thermo_params.Φ, thermo_params.Φbar, thermo_params.ξ,
        quark_params.A.s, quark_params.A.s,
        2,
    )

    Π_uu = ComplexF64(Π_uu_re, Π_uu_im)
    Π_ss = ComplexF64(Π_ss_re, Π_ss_im)
    M00, M08, M88 = Main.MesonPropagator.calculate_coupling_elements(Π_uu, Π_ss, K_coeffs, :P)
    detM = M00 * M88 - M08 * M08

    return (has_mixed=true, D_mixed=D_mixed, detM=detM)
end

function compute_rate_band_stats(process::Symbol, state, tau_cfg)
    p_grid, p_w = p_grid_from_mode(tau_cfg.integration_mode, tau_cfg.p_nodes)
    cos_grid, cos_w = Main.GaussLegendre.gauleg(-1.0, 1.0, tau_cfg.angle_nodes)
    phi_grid, phi_w = Main.GaussLegendre.gauleg(0.0, 2 * pi, tau_cfg.phi_nodes)

    band_omega_ref = Ref(Vector{Float64}())
    band_omega_sigma_ref = Ref(Vector{Float64}())

    rate = Main.AverageScatteringRate.average_scattering_rate(
        process,
        state.quark_params,
        state.thermo_params,
        state.K_coeffs;
        p_nodes=tau_cfg.p_nodes,
        angle_nodes=tau_cfg.angle_nodes,
        phi_nodes=tau_cfg.phi_nodes,
        p_grid=p_grid,
        p_w=p_w,
        cos_grid=cos_grid,
        cos_w=cos_w,
        phi_grid=phi_grid,
        phi_w=phi_w,
        n_sigma_points=tau_cfg.n_sigma_points,
        sigma_cutoff=Main.Constants_PNJL.Λ_inv_fm,
        threshold_subtraction=true,
        asym_window=0.6,
        asym_fit_min_points=8,
        asym_extra_points=10,
        interpolation_mode=tau_cfg.interpolation_mode,
        band_edges=BAND_EDGES,
        band_omega_out=band_omega_ref,
        band_omega_sigma_out=band_omega_sigma_ref,
    )

    pi_sym, pj_sym, _, _ = Main.TotalCrossSection.parse_particles_from_process(process)
    mi = Main.AverageScatteringRate.get_mass(pi_sym, state.quark_params)
    mj = Main.AverageScatteringRate.get_mass(pj_sym, state.quark_params)
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
        angle_nodes=tau_cfg.angle_nodes,
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
        angle_nodes=tau_cfg.angle_nodes,
    )

    prefactor = (Main.AverageScatteringRate.DQ^2) / (32.0 * pi^5 * ρ_i * ρ_j)
    omega_sigma = band_omega_sigma_ref[]
    omega = band_omega_ref[]

    return (
        rate=rate,
        prefactor=prefactor,
        omega_by_bin=omega,
        omega_sigma_by_bin=omega_sigma,
        omega_sigma_total=sum(omega_sigma),
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

    ds_grid = collect(range(1.0e-4, 20.0; length=260))

    out_detail = raw"D:\Desktop\Temp\relaxtime_t200_window\t200_denominator_fullband_detail.csv"
    out_pair = raw"D:\Desktop\Temp\relaxtime_t200_window\t200_denominator_fullband_pair_delta.csv"
    out_chain = raw"D:\Desktop\Temp\relaxtime_t200_window\t200_fullchain_band_table.csv"
    out_sample = raw"D:\Desktop\Temp\relaxtime_t200_window\t200_denominator_ds_samples.csv"
    ensure_parent_dir(out_detail)
    ensure_parent_dir(out_pair)
    ensure_parent_dir(out_chain)
    ensure_parent_dir(out_sample)

    io_d = open(out_detail, "w")
    io_p = open(out_pair, "w")
    io_c = open(out_chain, "w")
    io_s = open(out_sample, "w")
    try
        println(io_d, "scenario,process,xi,band,ds_left,ds_right,area_invabs_den_simple,area_invabs_detM,area_abs_D_simple_sq,area_abs_D_mixed_sq,area_abs_D_total_sq,area_sigma")
        println(io_p, "scenario,process,xi_left,xi_right,band,delta_area_invabs_den_simple,delta_area_invabs_detM,delta_area_abs_D_simple_sq,delta_area_abs_D_mixed_sq,delta_area_abs_D_total_sq,delta_area_sigma")
        println(io_c, "scenario,process,species,xi,rate_func,rate_direct,contrib_direct,band,ds_left,ds_right,area_invabs_den_simple,area_invabs_detM,area_abs_D_total_sq,area_sigma,omega_bin,omega_sigma_bin,sigma_eff_bin,rate_bin")
        println(io_s, "scenario,process,xi,ds,s_value,den_simple_re,den_simple_im,detM_re,detM_im,invabs_den_simple,invabs_detM,D_simple_re,D_simple_im,D_mixed_re,D_mixed_im,D_total_re,D_total_im,abs_D_simple_sq,abs_D_mixed_sq,abs_D_total_sq,sigma")

        cache_area = Dict{Tuple{String, Float64, Int}, NamedTuple}()
        cache_rate = Dict{Tuple{String, Float64}, NamedTuple}()

        for sc in scenarios
            main_rows = collect(CSV.File(sc.main_csv; comment="#"))
            channel_rows = collect(CSV.File(sc.channel_csv; comment="#"))
            tau_cfg = tau_settings_from_scan_metadata(sc.main_csv)

            states = Dict{Float64, Any}()
            for xi in sc.xis
                row = only(filter(r -> Float64(r.T_MeV) == 200.0 && Float64(r.muB_MeV) == 0.0 && Float64(r.xi) == xi, main_rows))
                states[xi] = build_state_from_main_row(row)
                cache_rate[(sc.name, xi)] = compute_rate_band_stats(sc.process, states[xi], tau_cfg)
            end

            for xi in sc.xis
                st = states[xi]
                th = process_threshold_info(sc.process, st.quark_params)

                den_simple_invabs = Float64[]
                detM_invabs = Float64[]
                abs_D_simple_sq = Float64[]
                abs_D_mixed_sq = Float64[]
                abs_D_total_sq = Float64[]
                sigma_vals = Float64[]

                for ds in ds_grid
                    s = th.s_th + ds
                    tb = Main.TotalCrossSection.calculate_t_bounds(s, th.mi, th.mj, th.mc, th.md)
                    t = 0.5 * (tb.t_min + tb.t_max)

                    cms = Main.TotalPropagator.calculate_cms_momentum(sc.process, s, t, :s, st.quark_params)
                    simple = simple_p_branch_for_channel(sc.process, :s, cms.k0, cms.k, st.quark_params, st.thermo_params, st.K_coeffs)
                    mixed = mixed_p_branch_for_channel(sc.process, :s, cms.k0, cms.k, st.quark_params, st.thermo_params, st.K_coeffs)

                    D_total = simple.D_simple + mixed.D_mixed
                    sigma = Main.TotalCrossSection.total_cross_section(sc.process, s, st.quark_params, st.thermo_params, st.K_coeffs; n_points=tau_cfg.n_sigma_points)

                    push!(den_simple_invabs, 1.0 / max(abs2(simple.den_ref), 1.0e-30))
                    push!(detM_invabs, mixed.has_mixed ? 1.0 / max(abs2(mixed.detM), 1.0e-30) : NaN)
                    push!(abs_D_simple_sq, abs2(simple.D_simple))
                    push!(abs_D_mixed_sq, abs2(mixed.D_mixed))
                    push!(abs_D_total_sq, abs2(D_total))
                    push!(sigma_vals, sigma)
                end

                for b in 1:(length(BAND_EDGES) - 1)
                    l = BAND_EDGES[b]
                    r = BAND_EDGES[b + 1]
                    idx = findall(d -> d > l && d <= r, ds_grid)
                    length(idx) < 2 && continue

                    x = ds_grid[idx]
                    area_den = trapz(x, den_simple_invabs[idx])
                    v_det = detM_invabs[idx]
                    area_det = all(isnan, v_det) ? NaN : trapz(x, replace(v_det, NaN => 0.0))
                    area_ds = trapz(x, abs_D_simple_sq[idx])
                    area_dm = trapz(x, abs_D_mixed_sq[idx])
                    area_dt = trapz(x, abs_D_total_sq[idx])
                    area_sigma = trapz(x, sigma_vals[idx])

                    cache_area[(sc.name, xi, b)] = (
                        area_den=area_den,
                        area_det=area_det,
                        area_ds=area_ds,
                        area_dm=area_dm,
                        area_dt=area_dt,
                        area_sigma=area_sigma,
                    )

                    println(io_d, join((
                        sc.name,
                        sc.process,
                        xi,
                        b,
                        l,
                        r,
                        area_den,
                        area_det,
                        area_ds,
                        area_dm,
                        area_dt,
                        area_sigma,
                    ), ','))
                end

                ch = filter(r -> Float64(r.T_MeV) == 200.0 && Float64(r.muB_MeV) == 0.0 && Float64(r.xi) == xi && String(r.species) == sc.species && String(r.channel) == string(sc.process), channel_rows)
                rate_direct = isempty(ch) ? NaN : Float64(ch[1].rate)
                contrib_direct = isempty(ch) ? NaN : Float64(ch[1].contribution)
                rate_stats = cache_rate[(sc.name, xi)]

                for b in 1:(length(BAND_EDGES) - 1)
                    haskey(cache_area, (sc.name, xi, b)) || continue
                    a = cache_area[(sc.name, xi, b)]
                    omega_bin = rate_stats.omega_by_bin[b]
                    omega_sigma_bin = rate_stats.omega_sigma_by_bin[b]
                    sigma_eff = omega_bin == 0.0 ? NaN : (omega_sigma_bin / omega_bin)
                    rate_bin = rate_stats.prefactor * omega_sigma_bin
                    println(io_c, join((
                        sc.name,
                        sc.process,
                        sc.species,
                        xi,
                        rate_stats.rate,
                        rate_direct,
                        contrib_direct,
                        b,
                        BAND_EDGES[b],
                        BAND_EDGES[b + 1],
                        a.area_den,
                        a.area_det,
                        a.area_dt,
                        a.area_sigma,
                        omega_bin,
                        omega_sigma_bin,
                        sigma_eff,
                        rate_bin,
                    ), ','))
                end

                for ds in SAMPLE_DS_POINTS
                    s = th.s_th + ds
                    tb = Main.TotalCrossSection.calculate_t_bounds(s, th.mi, th.mj, th.mc, th.md)
                    t = 0.5 * (tb.t_min + tb.t_max)
                    cms = Main.TotalPropagator.calculate_cms_momentum(sc.process, s, t, :s, st.quark_params)
                    simple = simple_p_branch_for_channel(sc.process, :s, cms.k0, cms.k, st.quark_params, st.thermo_params, st.K_coeffs)
                    mixed = mixed_p_branch_for_channel(sc.process, :s, cms.k0, cms.k, st.quark_params, st.thermo_params, st.K_coeffs)
                    D_total = simple.D_simple + mixed.D_mixed
                    sigma = Main.TotalCrossSection.total_cross_section(sc.process, s, st.quark_params, st.thermo_params, st.K_coeffs; n_points=tau_cfg.n_sigma_points)

                    invabs_den = 1.0 / max(abs2(simple.den_ref), 1.0e-30)
                    invabs_det = mixed.has_mixed ? 1.0 / max(abs2(mixed.detM), 1.0e-30) : NaN
                    det_re = mixed.has_mixed ? real(mixed.detM) : NaN
                    det_im = mixed.has_mixed ? imag(mixed.detM) : NaN

                    println(io_s, join((
                        sc.name,
                        sc.process,
                        xi,
                        ds,
                        s,
                        real(simple.den_ref),
                        imag(simple.den_ref),
                        det_re,
                        det_im,
                        invabs_den,
                        invabs_det,
                        real(simple.D_simple),
                        imag(simple.D_simple),
                        real(mixed.D_mixed),
                        imag(mixed.D_mixed),
                        real(D_total),
                        imag(D_total),
                        abs2(simple.D_simple),
                        abs2(mixed.D_mixed),
                        abs2(D_total),
                        sigma,
                    ), ','))
                end
            end

            for i in 1:(length(sc.xis) - 1)
                xl = sc.xis[i]
                xr = sc.xis[i + 1]
                for b in 1:(length(BAND_EDGES) - 1)
                    haskey(cache_area, (sc.name, xl, b)) || continue
                    haskey(cache_area, (sc.name, xr, b)) || continue
                    l = cache_area[(sc.name, xl, b)]
                    r = cache_area[(sc.name, xr, b)]

                    println(io_p, join((
                        sc.name,
                        sc.process,
                        xl,
                        xr,
                        b,
                        r.area_den - l.area_den,
                        r.area_det - l.area_det,
                        r.area_ds - l.area_ds,
                        r.area_dm - l.area_dm,
                        r.area_dt - l.area_dt,
                        r.area_sigma - l.area_sigma,
                    ), ','))
                end
            end
        end
    finally
        close(io_d)
        close(io_p)
        close(io_c)
        close(io_s)
    end
end

main()
