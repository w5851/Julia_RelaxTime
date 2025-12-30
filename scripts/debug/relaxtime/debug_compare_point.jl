#!/usr/bin/env julia

using StaticArrays

function find_project_root(start_dir::AbstractString)
    dir = abspath(start_dir)
    while true
        if isfile(joinpath(dir, "Project.toml"))
            return dir
        end
        parent = dirname(dir)
        parent == dir && error("Could not find Project.toml from: $start_dir")
        dir = parent
    end
end

const PROJECT_ROOT = find_project_root(@__DIR__)

push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src"))
push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src", "relaxtime"))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "ScatteringAmplitude.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "DifferentialCrossSection.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalCrossSection.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalPropagator.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "PolarizationCache.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using .OneLoopIntegrals: A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .ScatteringAmplitude: scattering_amplitude_squared
using .DifferentialCrossSection: differential_cross_section
using .TotalCrossSection: total_cross_section
using .TotalPropagator: calculate_all_propagators_by_channel, calculate_cms_momentum
using .PolarizationCache: polarization_aniso_cached

function main()
    T_MeV = 150.0
    muB_MeV = 800.0
    muq_MeV = muB_MeV / 3

    T = T_MeV / ħc_MeV_fm
    muq = muq_MeV / ħc_MeV_fm

    res = PNJL.solve(PNJL.FixedMu(), T, muq; xi=0.0)
    res.converged || error("gap solver did not converge")

    x = res.solution
    phi = SVector{3, Float64}(x[1], x[2], x[3])
    Phi = Float64(x[4])
    Phibar = Float64(x[5])

    m_u, m_d, m_s = res.masses

    println("gap: m_u(MeV)=", m_u * ħc_MeV_fm, " m_d(MeV)=", m_d * ħc_MeV_fm, " Phi=", Phi, " Phibar=", Phibar)

    nodes_p = DEFAULT_MOMENTUM_NODES
    weights_p = DEFAULT_MOMENTUM_WEIGHTS

    A_u = A(m_u, muq, T, Phi, Phibar, nodes_p, weights_p)
    A_s = A(m_s, muq, T, Phi, Phibar, nodes_p, weights_p)
    G_u = calculate_G_from_A(A_u, m_u)
    G_s = calculate_G_from_A(A_s, m_s)
    Kc = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

    println("Kc: K123_plus=", Kc.K123_plus, " K123_minus=", Kc.K123_minus,
        "  K0_plus=", Kc.K0_plus, " K0_minus=", Kc.K0_minus,
        "  K8_plus=", Kc.K8_plus, " K8_minus=", Kc.K8_minus,
        "  K08_plus=", Kc.K08_plus, " K08_minus=", Kc.K08_minus,
        "  detK_plus=", Kc.det_K_plus, " detK_minus=", Kc.det_K_minus)

    quark_params = (m=(u=m_u, d=m_d, s=m_s), μ=(u=muq, d=muq, s=muq), A=(u=A_u, d=A_u, s=A_s))
    thermo_params = (T=T, Φ=Phi, Φbar=Phibar, ξ=0.0)

    process = :udbar_to_udbar
    s = 13.7258

    println("process=", process, " s=", s)

    # include the same representative t as C++ debug_amp_vs_t (midpoint of t range)
    for t in (-13.0321, -6.52157, -0.0110493)
        Msq = scattering_amplitude_squared(process, s, t, quark_params, thermo_params, Kc)
        splus = s - (m_u + m_d)^2
        sminus = s - (m_u - m_d)^2
        dsdt = differential_cross_section(splus, sminus, Msq)

        cms_s = calculate_cms_momentum(process, s, t, :s, quark_params)
        cms_t = calculate_cms_momentum(process, s, t, :t, quark_params)

        # Polarization cross-checks (should match C++ header PiP/PiS)
        PiP_s_re, PiP_s_im = polarization_aniso_cached(:P, cms_s.k0, cms_s.k, m_u, m_d, muq, muq, T, Phi, Phibar, 0.0, A_u, A_u, 0)
        PiS_s_re, PiS_s_im = polarization_aniso_cached(:S, cms_s.k0, cms_s.k, m_u, m_d, muq, muq, T, Phi, Phibar, 0.0, A_u, A_u, 0)
        PiP_t_re, PiP_t_im = polarization_aniso_cached(:P, cms_t.k0, cms_t.k, m_u, m_u, muq, muq, T, Phi, Phibar, 0.0, A_u, A_u, 0)
        PiS_t_re, PiS_t_im = polarization_aniso_cached(:S, cms_t.k0, cms_t.k, m_u, m_u, muq, muq, T, Phi, Phibar, 0.0, A_u, A_u, 0)
        props_s = calculate_all_propagators_by_channel(process, cms_s.k0, cms_s.k, quark_params, thermo_params, Kc)
        props_t = calculate_all_propagators_by_channel(process, cms_t.k0, cms_t.k, quark_params, thermo_params, Kc)

        println("t=", t)
        println("  Pi_s: P(re,im)=", PiP_s_re, ",", PiP_s_im, "  S(re,im)=", PiS_s_re, ",", PiS_s_im)
        println("  Pi_t: P(re,im)=", PiP_t_re, ",", PiP_t_im, "  S(re,im)=", PiS_t_re, ",", PiS_t_im)
        println("  cms_s: k0=", cms_s.k0, " k=", cms_s.k, "  D_s_S=", props_s.s_S, "  D_s_P=", props_s.s_P)
        println("  cms_t: k0=", cms_t.k0, " k=", cms_t.k, "  D_t_S=", props_t.t_S, "  D_t_P=", props_t.t_P)
        println("  Msq=", Msq, " dsdt=", dsdt)
    end

    println("sigma_total(n6)=", total_cross_section(process, s, quark_params, thermo_params, Kc; n_points=6))
    println("sigma_total(n64)=", total_cross_section(process, s, quark_params, thermo_params, Kc; n_points=64))
end

main()
