#!/usr/bin/env julia

using Printf

println("Starting integrand_vs_s_average_rate.jl ...")
flush(stdout)

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src"))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "AverageScatteringRate.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalCrossSection.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5, Λ_inv_fm
using .GaussLegendre: gauleg, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using .OneLoopIntegrals: A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .TotalCrossSection: DEFAULT_T_INTEGRAL_POINTS

const ASR = AverageScatteringRate

"""
Build quark/thermo parameters for inspecting the ω integrand.
Matches the defaults used in the cross-section scan helpers.
"""
function build_params(; T_MeV=150.0, mu_u_MeV=0.0, mu_d_MeV=0.0, mu_s_MeV=0.0,
    m_u_MeV=300.0, m_d_MeV=300.0, m_s_MeV=500.0, phi=0.5, phibar=0.5, xi=0.0)
    T = T_MeV / ħc_MeV_fm
    μ_u = mu_u_MeV / ħc_MeV_fm
    μ_d = mu_d_MeV / ħc_MeV_fm
    μ_s = mu_s_MeV / ħc_MeV_fm
    m_u = m_u_MeV / ħc_MeV_fm
    m_d = m_d_MeV / ħc_MeV_fm
    m_s = m_s_MeV / ħc_MeV_fm

    nodes_p = DEFAULT_MOMENTUM_NODES
    weights_p = DEFAULT_MOMENTUM_WEIGHTS

    A_u = A(m_u, μ_u, T, phi, phibar, nodes_p, weights_p)
    A_d = A(m_d, μ_d, T, phi, phibar, nodes_p, weights_p)
    A_s = A(m_s, μ_s, T, phi, phibar, nodes_p, weights_p)

    G_u = calculate_G_from_A(A_u)
    G_s = calculate_G_from_A(A_s)

    quark_params = (m=(u=m_u, d=m_d, s=m_s), μ=(u=μ_u, d=μ_d, s=μ_s), A=(u=A_u, d=A_d, s=A_s))
    thermo_params = (T=T, Φ=phi, Φbar=phibar, ξ=xi)
    K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
    return (quark_params=quark_params, thermo_params=thermo_params, K_coeffs=K_coeffs)
end

"""
Compute the raw ω increment (integrand term) for the given phase-space point.
Returns `nothing` when the kinematic cuts null the contribution.
"""
function omega_term_at(process::Symbol, p_i::Float64, p_j::Float64, cθi::Float64, cθj::Float64, φ::Float64;
    params, cache=ASR.CrossSectionCache(process), n_sigma_points::Int=DEFAULT_T_INTEGRAL_POINTS)
    pi_sym, pj_sym, _, _ = ASR.parse_particles_from_process(process)
    mi = ASR.get_mass(pi_sym, params.quark_params)
    mj = ASR.get_mass(pj_sym, params.quark_params)
    μi = ASR.get_mu(pi_sym, params.quark_params)
    μj = ASR.get_mu(pj_sym, params.quark_params)

    T = params.thermo_params.T
    Φ = params.thermo_params.Φ
    Φbar = params.thermo_params.Φbar
    ξ = hasproperty(params.thermo_params, :ξ) ? params.thermo_params.ξ : 0.0

    Ei = ASR.energy_from_p(p_i, mi)
    Ej = ASR.energy_from_p(p_j, mj)
    sθi = sqrt(max(1.0 - cθi * cθi, 0.0))
    sθj = sqrt(max(1.0 - cθj * cθj, 0.0))

    f_i = ASR.distribution_with_anisotropy(pi_sym, p_i, mi, μi, T, Φ, Φbar, ξ, cθi)
    f_j = ASR.distribution_with_anisotropy(pj_sym, p_j, mj, μj, T, Φ, Φbar, ξ, cθj)
    if f_i == 0.0 || f_j == 0.0
        return nothing
    end

    cosΘ = cθi * cθj + sθi * sθj * cos(φ)
    p_dot = Ei * Ej - p_i * p_j * cosΘ
    s = mi^2 + mj^2 + 2.0 * p_dot
    if s <= (mi + mj)^2
        return nothing
    end

    v_rel_num = p_dot^2 - mi^2 * mj^2
    v_rel = v_rel_num > 0 ? sqrt(v_rel_num) / (Ei * Ej) : 0.0
    if v_rel == 0.0
        return nothing
    end

    σ = ASR.get_sigma(cache, s, params.quark_params, params.thermo_params, params.K_coeffs; n_points=n_sigma_points)
    term = (p_i ^ 2) * (p_j ^ 2) * f_i * f_j * v_rel * σ
    return (s=s, term=term, v_rel=v_rel, σ=σ, p_i=p_i, p_j=p_j, cθi=cθi, cθj=cθj, φ=φ)
end

"""
Sweep quadrature nodes and collect ω integrand contributions versus s.
"""
function sample_omega_terms(process::Symbol; params, p_nodes::Int=6, angle_nodes::Int=4, phi_nodes::Int=8,
    n_sigma_points::Int=DEFAULT_T_INTEGRAL_POINTS)
    p_grid, p_w = gauleg(0.0, Λ_inv_fm, p_nodes)
    cos_grid, cos_w = gauleg(0.0, 1.0, angle_nodes)
    phi_grid, phi_w = gauleg(0.0, Float64(π), phi_nodes)

    cache = ASR.CrossSectionCache(process)
    rows = NamedTuple[]

    for (p_i, w_pi) in zip(p_grid, p_w)
        for (p_j, w_pj) in zip(p_grid, p_w)
            for (cθi, w_cθi) in zip(cos_grid, cos_w)
                for (cθj, w_cθj) in zip(cos_grid, cos_w)
                    for (φ, wφ) in zip(phi_grid, phi_w)
                        rec = omega_term_at(process, p_i, p_j, cθi, cθj, φ; params=params, cache=cache, n_sigma_points=n_sigma_points)
                        if rec === nothing
                            continue
                        end
                        weighted = rec.term * w_pi * w_pj * w_cθi * w_cθj * wφ
                        push!(rows, merge(rec, (; weighted=weighted)))
                    end
                end
            end
        end
    end

    sort!(rows, by = r -> r.s)
    return rows
end

function print_table(rows; max_print=50)
    n = length(rows)
    @printf("Collected %d contributing phase-space points. Showing up to %d by ascending s.\n", n, min(max_print, n))
    println("s\tterm\tweighted\tv_rel\tσ\tp_i\tp_j\tcθi\tcθj\tφ")
    for (idx, r) in enumerate(rows)
        if idx > max_print
            break
        end
        @printf("%.8e\t%.8e\t%.8e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3f\t%.3f\t%.3f\n",
            r.s, r.term, r.weighted, r.v_rel, r.σ, r.p_i, r.p_j, r.cθi, r.cθj, r.φ)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    params = build_params()
    rows = sample_omega_terms(:ssbar_to_uubar; params=params, p_nodes=6, angle_nodes=4, phi_nodes=8, n_sigma_points=16)
    print_table(rows; max_print=60)
end
