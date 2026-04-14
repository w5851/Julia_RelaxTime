#!/usr/bin/env julia

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(@__DIR__, "t190_sigma_chain_decomposition_lib.jl"))

function ensure_parent_dir(path::AbstractString)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
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

function simple_p_for(process::Symbol, channel::Symbol, k0::Float64, k_norm::Float64,
    quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple)
    info = Main.Constants_PNJL.SCATTERING_MESON_MAP[process][:channels][channel]
    simple_mesons = info[:simple]
    T1, T2 = Main.TotalPropagator.get_flavor_factors_for_channel(process, channel)

    D_simple = 0.0 + 0.0im
    den_ref = NaN + NaN * im
    meson_ref = Symbol("")
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
        if meson_ref == Symbol("")
            meson_ref = meson
            den_ref = den
        end
    end
    return (D_simple=D_simple, den_ref=den_ref, meson_ref=meson_ref)
end

function mixed_p_for(process::Symbol, channel::Symbol, k0::Float64, k_norm::Float64,
    quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple)
    info = Main.Constants_PNJL.SCATTERING_MESON_MAP[process][:channels][channel]
    if !Bool(info[:mixed_P])
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

function main()
    T_MeV = 200.0
    muB_MeV = 0.0
    process = :usbar_to_usbar
    xi_values = collect(range(-0.3, -0.1; length=11))
    ds = 1.0e-8

    out = raw"D:\Desktop\Temp\relaxtime_t200_window\t200_taus_usbar_denominator_trace.csv"
    ensure_parent_dir(out)

    io = open(out, "w")
    try
        println(io, "xi,channel,s_minus_sth,k0,k_norm,meson_simple_ref,re_den_simple,im_den_simple,has_mixed,re_detM,im_detM,abs_D_simple_sq,abs_D_mixed_sq,abs_D_total_sq,M2_total,M_s_sq,M_t_sq,M_interf,sigma_total")
        for xi in xi_values
            st = build_state_point(T_MeV, muB_MeV, xi)
            th = process_threshold_info(process, st.quark_params)
            s = th.s_th + ds
            tb = Main.TotalCrossSection.calculate_t_bounds(s, th.mi, th.mj, th.mc, th.md)
            t = 0.5 * (tb.t_min + tb.t_max)
            amp = decompose_qqbar_amplitude_terms(process, s, t, st.quark_params, st.thermo_params, st.K_coeffs)
            sigma_val = Main.TotalCrossSection.total_cross_section(process, s, st.quark_params, st.thermo_params, st.K_coeffs; n_points=24)

            for channel in (:s, :t)
                cms = Main.TotalPropagator.calculate_cms_momentum(process, s, t, channel, st.quark_params)
                simple = simple_p_for(process, channel, cms.k0, cms.k, st.quark_params, st.thermo_params, st.K_coeffs)
                mixed = mixed_p_for(process, channel, cms.k0, cms.k, st.quark_params, st.thermo_params, st.K_coeffs)
                D_total = simple.D_simple + mixed.D_mixed
                println(io, join((
                    xi,
                    string(channel),
                    ds,
                    cms.k0,
                    cms.k,
                    string(simple.meson_ref),
                    real(simple.den_ref),
                    imag(simple.den_ref),
                    mixed.has_mixed,
                    real(mixed.detM),
                    imag(mixed.detM),
                    abs2(simple.D_simple),
                    abs2(mixed.D_mixed),
                    abs2(D_total),
                    amp.M2_total,
                    amp.M_s_sq,
                    amp.M_t_sq,
                    amp.M_interf,
                    sigma_val,
                ), ','))
            end
        end
    finally
        close(io)
    end
end

main()
