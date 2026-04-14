using StaticArrays

if !isdefined(Main, :Constants_PNJL)
    include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
end
if !isdefined(Main, :GaussLegendre)
    include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
end
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Main.Models.pnjl_module()
if !isdefined(Main, :OneLoopIntegrals)
    include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))
end
if !isdefined(Main, :EffectiveCouplings)
    include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
end
if !isdefined(Main, :ScatteringAmplitude)
    include(joinpath(PROJECT_ROOT, "src", "relaxtime", "ScatteringAmplitude.jl"))
end
if !isdefined(Main, :TotalCrossSection)
    include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalCrossSection.jl"))
end
if !isdefined(Main, :AverageScatteringRate)
    include(joinpath(PROJECT_ROOT, "src", "relaxtime", "AverageScatteringRate.jl"))
end
if !isdefined(Main, :TotalPropagator)
    include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalPropagator.jl"))
end

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5
using .OneLoopIntegrals: A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .ScatteringAmplitude: calculate_mandelstam_variables

const PNJL_LIB = Main.Models.pnjl_module()
const solve_lib = getproperty(PNJL_LIB, :solve)
const FixedMu_lib = getproperty(PNJL_LIB, :FixedMu)
const Integrals_lib = getproperty(PNJL_LIB, :Integrals)
const DEFAULT_MOMENTUM_NODES_LIB = getproperty(Integrals_lib, :DEFAULT_MOMENTUM_NODES)
const DEFAULT_MOMENTUM_WEIGHTS_LIB = getproperty(Integrals_lib, :DEFAULT_MOMENTUM_WEIGHTS)

struct StatePoint
    xi::Float64
    quark_params::NamedTuple
    thermo_params::NamedTuple
    K_coeffs::NamedTuple
end

@inline function get_mass(flavor::Symbol, quark_params::NamedTuple)
    if flavor === :u || flavor === :ubar
        return quark_params.m.u
    elseif flavor === :d || flavor === :dbar
        return quark_params.m.d
    elseif flavor === :s || flavor === :sbar
        return quark_params.m.s
    end
    error("unknown flavor: $flavor")
end

@inline function get_mu(flavor::Symbol, quark_params::NamedTuple)
    if flavor === :u || flavor === :ubar
        return quark_params.μ.u
    elseif flavor === :d || flavor === :dbar
        return quark_params.μ.d
    elseif flavor === :s || flavor === :sbar
        return quark_params.μ.s
    end
    error("unknown flavor: $flavor")
end

function build_state_point(T_MeV::Float64, muB_MeV::Float64, xi::Float64)
    T_fm = T_MeV / ħc_MeV_fm
    muq_fm = (muB_MeV / 3.0) / ħc_MeV_fm

    base = solve_lib(FixedMu_lib(), T_fm, muq_fm; xi=xi, p_num=12, t_num=6)
    Bool(base.converged) || error("equilibrium not converged for xi=$xi")

    Phi = Float64(base.x_state[4])
    Phibar = Float64(base.x_state[5])
    masses = (u=Float64(base.masses[1]), d=Float64(base.masses[2]), s=Float64(base.masses[3]))

    A_u = A(masses.u, muq_fm, T_fm, Phi, Phibar, DEFAULT_MOMENTUM_NODES_LIB, DEFAULT_MOMENTUM_WEIGHTS_LIB)
    A_s = A(masses.s, muq_fm, T_fm, Phi, Phibar, DEFAULT_MOMENTUM_NODES_LIB, DEFAULT_MOMENTUM_WEIGHTS_LIB)
    G_u = calculate_G_from_A(A_u, masses.u)
    G_s = calculate_G_from_A(A_s, masses.s)
    K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

    quark_params = (m=masses, μ=(u=muq_fm, d=muq_fm, s=muq_fm), A=(u=A_u, d=A_u, s=A_s))
    thermo_params = (T=T_fm, Φ=Phi, Φbar=Phibar, ξ=xi)
    return StatePoint(xi, quark_params, thermo_params, K_coeffs)
end

function process_threshold_info(process::Symbol, quark_params::NamedTuple)
    i, j, c, d = Main.AverageScatteringRate.parse_particles_from_process(process)
    mi = get_mass(i, quark_params)
    mj = get_mass(j, quark_params)
    mc = get_mass(c, quark_params)
    md = get_mass(d, quark_params)
    s_th = max((mi + mj)^2, (mc + md)^2)
    return (i=i, j=j, c=c, d=d, mi=mi, mj=mj, mc=mc, md=md, s_th=s_th)
end

function decompose_qqbar_amplitude_terms(process::Symbol, s::Float64, t::Float64,
    quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple)
    m1, m2, m3, m4 = Main.ParticleSymbols.get_quark_masses_for_process(process, quark_params)
    u = m1^2 + m2^2 + m3^2 + m4^2 - s - t

    cms_s = Main.TotalPropagator.calculate_cms_momentum(process, s, t, :s, quark_params; u=u)
    cms_t = Main.TotalPropagator.calculate_cms_momentum(process, s, t, :t, quark_params; u=u)
    D_s_S, D_s_P = Main.TotalPropagator.calculate_propagator_by_channel(process, :s, cms_s.k0, cms_s.k, quark_params, thermo_params, K_coeffs)
    D_t_S, D_t_P = Main.TotalPropagator.calculate_propagator_by_channel(process, :t, cms_t.k0, cms_t.k, quark_params, thermo_params, K_coeffs)

    vars = calculate_mandelstam_variables(s, t, u, m1, m2, m3, m4)

    M_s_S = abs2(D_s_S) * vars.s_12_plus * vars.s_34_plus
    M_s_P = abs2(D_s_P) * vars.s_12_minus * vars.s_34_minus
    M_t_S = abs2(D_t_S) * vars.t_13_plus * vars.t_24_plus
    M_t_P = abs2(D_t_P) * vars.t_13_minus * vars.t_24_minus
    M_s_sq = M_s_S + M_s_P
    M_t_sq = M_t_S + M_t_P

    cross_factor = 1.0 / (4.0 * Main.Constants_PNJL.N_color)
    D_t_S_conj = conj(D_t_S)
    D_t_P_conj = conj(D_t_P)

    term1 = D_s_S * D_t_S_conj * (vars.s_12_plus * vars.s_34_plus - vars.u_14_plus * vars.u_23_plus + vars.t_13_plus * vars.t_24_plus)
    term2 = D_s_S * D_t_P_conj * (vars.s_12_plus * vars.s_34_plus - vars.u_14_minus * vars.u_23_minus + vars.t_13_minus * vars.t_24_minus)
    term3 = D_s_P * D_t_S_conj * (vars.s_12_minus * vars.s_34_minus - vars.u_14_minus * vars.u_23_minus + vars.t_13_plus * vars.t_24_plus)
    term4 = D_s_P * D_t_P_conj * (vars.s_12_minus * vars.s_34_minus - vars.u_14_plus * vars.u_23_plus + vars.t_13_minus * vars.t_24_minus)
    cross_term = cross_factor * (term1 - term2 - term3 + term4)
    M_interf = -2.0 * real(cross_term)

    M2_total = M_s_sq + M_t_sq + M_interf
    return (
        M2_total=M2_total,
        M_s_sq=M_s_sq,
        M_t_sq=M_t_sq,
        M_interf=M_interf,
        M_s_S=M_s_S,
        M_s_P=M_s_P,
        M_t_S=M_t_S,
        M_t_P=M_t_P,
    )
end

function decompose_qqbar_s_channel_factors(process::Symbol, s::Float64, t::Float64,
    quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple)
    m1, m2, m3, m4 = Main.ParticleSymbols.get_quark_masses_for_process(process, quark_params)
    u = m1^2 + m2^2 + m3^2 + m4^2 - s - t
    cms_s = Main.TotalPropagator.calculate_cms_momentum(process, s, t, :s, quark_params; u=u)
    D_s_S, D_s_P = Main.TotalPropagator.calculate_propagator_by_channel(process, :s, cms_s.k0, cms_s.k, quark_params, thermo_params, K_coeffs)
    vars = calculate_mandelstam_variables(s, t, u, m1, m2, m3, m4)

    abs_D_s_S_sq = abs2(D_s_S)
    abs_D_s_P_sq = abs2(D_s_P)
    kin_s_S = vars.s_12_plus * vars.s_34_plus
    kin_s_P = vars.s_12_minus * vars.s_34_minus
    M_s_S = abs_D_s_S_sq * kin_s_S
    M_s_P = abs_D_s_P_sq * kin_s_P

    return (
        abs_D_s_S_sq=abs_D_s_S_sq,
        abs_D_s_P_sq=abs_D_s_P_sq,
        kin_s_S=kin_s_S,
        kin_s_P=kin_s_P,
        M_s_S=M_s_S,
        M_s_P=M_s_P,
    )
end

function decompose_p_channel_propagator_strength(process::Symbol, s::Float64, t::Float64,
    quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple)
    m1, m2, m3, m4 = Main.ParticleSymbols.get_quark_masses_for_process(process, quark_params)
    u = m1^2 + m2^2 + m3^2 + m4^2 - s - t
    vars = calculate_mandelstam_variables(s, t, u, m1, m2, m3, m4)
    kin_s_P = vars.s_12_minus * vars.s_34_minus

    cms_s = Main.TotalPropagator.calculate_cms_momentum(process, s, t, :s, quark_params; u=u)
    k0 = cms_s.k0
    k_norm = cms_s.k

    D_s_P_simple = 0.0 + 0.0im
    D_s_P_mixed = 0.0 + 0.0im

    process_map = Main.Constants_PNJL.SCATTERING_MESON_MAP[process]
    s_channel_info = process_map[:channels][:s]
    simple_mesons = s_channel_info[:simple]

    T1, T2 = Main.TotalPropagator.get_flavor_factors_for_channel(process, :s)
    for meson in simple_mesons
        if !(meson === :pi || meson === :K)
            continue
        end
        pol_params = Main.TotalPropagator.get_polarization_params(meson, quark_params)
        Π_real, Π_imag = Main.PolarizationCache.polarization_aniso_cached(
            pol_params.channel, k0, k_norm,
            pol_params.m1, pol_params.m2,
            pol_params.μ1, pol_params.μ2,
            thermo_params.T, thermo_params.Φ, thermo_params.Φbar, thermo_params.ξ,
            pol_params.A1, pol_params.A2,
            pol_params.num_s_quark,
        )
        Π = ComplexF64(Π_real, Π_imag)
        D_meson = Main.MesonPropagator.meson_propagator_simple(meson, K_coeffs, Π)
        D_s_P_simple += T1 * D_meson * T2
    end

    if Bool(s_channel_info[:mixed_P])
        D_s_P_mixed = Main.TotalPropagator.total_propagator_mixed(
            process, :s, :P, k0, k_norm, quark_params, thermo_params, K_coeffs,
        )
    end

    D_s_P_total = D_s_P_simple + D_s_P_mixed
    abs_D_s_P_sq_simple = abs2(D_s_P_simple)
    abs_D_s_P_sq_mixed = abs2(D_s_P_mixed)
    abs_D_s_P_sq_total = abs2(D_s_P_total)
    return (
        k0=k0,
        k_norm=k_norm,
        kin_s_P=kin_s_P,
        D_s_P_simple=D_s_P_simple,
        D_s_P_mixed=D_s_P_mixed,
        D_s_P_total=D_s_P_total,
        abs_D_s_P_sq_simple=abs_D_s_P_sq_simple,
        abs_D_s_P_sq_mixed=abs_D_s_P_sq_mixed,
        abs_D_s_P_sq_total=abs_D_s_P_sq_total,
        M_s_P_total=abs_D_s_P_sq_total * kin_s_P,
    )
end

function decompose_mixed_p_propagator_chain(process::Symbol, s::Float64, t::Float64,
    quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple)
    m1, m2, m3, m4 = Main.ParticleSymbols.get_quark_masses_for_process(process, quark_params)
    u = m1^2 + m2^2 + m3^2 + m4^2 - s - t
    cms_s = Main.TotalPropagator.calculate_cms_momentum(process, s, t, :s, quark_params; u=u)
    k0 = cms_s.k0
    k_norm = cms_s.k

    Π_uu_real, Π_uu_imag = Main.PolarizationCache.polarization_aniso_cached(
        :P, k0, k_norm,
        quark_params.m.u, quark_params.m.u,
        quark_params.μ.u, quark_params.μ.u,
        thermo_params.T, thermo_params.Φ, thermo_params.Φbar, thermo_params.ξ,
        quark_params.A.u, quark_params.A.u,
        0,
    )
    Π_ss_real, Π_ss_imag = Main.PolarizationCache.polarization_aniso_cached(
        :P, k0, k_norm,
        quark_params.m.s, quark_params.m.s,
        quark_params.μ.s, quark_params.μ.s,
        thermo_params.T, thermo_params.Φ, thermo_params.Φbar, thermo_params.ξ,
        quark_params.A.s, quark_params.A.s,
        2,
    )

    Π_uu = ComplexF64(Π_uu_real, Π_uu_imag)
    Π_ss = ComplexF64(Π_ss_real, Π_ss_imag)

    detK = K_coeffs.det_K_plus
    c_diag = (4.0 / 3.0) * detK
    c_offdiag = (4.0 / 3.0) * sqrt(2.0) * detK

    M00_from_K0 = K_coeffs.K0_plus
    M00_from_Piuu = -c_diag * Π_uu
    M00_from_Piss = -2.0 * c_diag * Π_ss

    M08_from_K08 = K_coeffs.K08_plus
    M08_from_Piuu = c_offdiag * Π_uu
    M08_from_Piss = -c_offdiag * Π_ss

    M88_from_K8 = K_coeffs.K8_plus
    M88_from_Piuu = -2.0 * c_diag * Π_uu
    M88_from_Piss = -c_diag * Π_ss

    M00, M08, M88 = Main.MesonPropagator.calculate_coupling_elements(Π_uu, Π_ss, K_coeffs, :P)
    detM_term_M00M88 = M00 * M88
    detM_term_M08sq = M08 * M08
    detM = detM_term_M00M88 - detM_term_M08sq
    detM_cross_term = -2.0 * real(detM_term_M00M88 * conj(detM_term_M08sq))

    q1, q2, q3, q4 = Main.ParticleSymbols.parse_scattering_process(process)
    J = Main.MesonPropagator.calculate_current_vector(q1, q2, :s)
    Jp = Main.MesonPropagator.calculate_current_vector(q3, q4, :s)
    J0, J8 = J[1], J[2]
    Jp0, Jp8 = Jp[1], Jp[2]

    JMJ = J0 * (M00 * Jp0 + M08 * Jp8) + J8 * (M08 * Jp0 + M88 * Jp8)
    prefactor = 2.0 * detK / detM
    D_mixed_P = prefactor * JMJ
    abs_detM_sq = abs2(detM)
    abs_JMJ_sq = abs2(JMJ)
    abs_D_mixed_P_sq = abs2(D_mixed_P)

    return (
        k0=k0,
        k_norm=k_norm,
        Π_uu=Π_uu,
        Π_ss=Π_ss,
        M00=M00,
        M08=M08,
        M88=M88,
        M00_from_K0=M00_from_K0,
        M00_from_Piuu=M00_from_Piuu,
        M00_from_Piss=M00_from_Piss,
        M08_from_K08=M08_from_K08,
        M08_from_Piuu=M08_from_Piuu,
        M08_from_Piss=M08_from_Piss,
        M88_from_K8=M88_from_K8,
        M88_from_Piuu=M88_from_Piuu,
        M88_from_Piss=M88_from_Piss,
        detM_term_M00M88=detM_term_M00M88,
        detM_term_M08sq=detM_term_M08sq,
        abs_detM_term_M00M88_sq=abs2(detM_term_M00M88),
        abs_detM_term_M08sq_sq=abs2(detM_term_M08sq),
        detM_cross_term=detM_cross_term,
        detM=detM,
        abs_detM_sq=abs_detM_sq,
        J0=J0,
        J8=J8,
        Jp0=Jp0,
        Jp8=Jp8,
        JMJ=JMJ,
        abs_JMJ_sq=abs_JMJ_sq,
        detK_plus=detK,
        prefactor=prefactor,
        D_mixed_P=D_mixed_P,
        abs_D_mixed_P_sq=abs_D_mixed_P_sq,
    )
end
