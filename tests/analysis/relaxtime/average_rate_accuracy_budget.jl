#!/usr/bin/env julia

using Printf

println("Starting average_rate_accuracy_budget.jl ...")
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
using .TotalCrossSection: total_cross_section

const ASR = AverageScatteringRate

@inline function relerr(a::Float64, b::Float64)
    return abs(a - b) / max(1e-12, abs(b))
end

"""Lightweight parameter builder.

Provides fields required by `TotalCrossSection.total_cross_section`:
- `quark_params.m`, `quark_params.μ`, and `quark_params.A`
- `thermo_params` and `K_coeffs`

This is meant for *numerical integration diagnostics* (not the PNJL gap solver).
"""
function build_params(; T_MeV::Float64, muB_MeV::Float64, phi::Float64, phibar::Float64,
    m_u_MeV::Float64=300.0, m_d_MeV::Float64=300.0, m_s_MeV::Float64=500.0,
    mu_s_MeV::Float64=0.0, xi::Float64=0.0)

    T = T_MeV / ħc_MeV_fm
    μ_q = (muB_MeV / 3.0) / ħc_MeV_fm
    μ_u = μ_q
    μ_d = μ_q
    μ_s = mu_s_MeV / ħc_MeV_fm

    m_u = m_u_MeV / ħc_MeV_fm
    m_d = m_d_MeV / ħc_MeV_fm
    m_s = m_s_MeV / ħc_MeV_fm

    A_u = A(m_u, μ_u, T, phi, phibar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    A_d = A(m_d, μ_d, T, phi, phibar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    A_s = A(m_s, μ_s, T, phi, phibar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)

    G_u = calculate_G_from_A(A_u, m_u)
    G_s = calculate_G_from_A(A_s, m_s)

    quark_params = (m=(u=m_u, d=m_d, s=m_s), μ=(u=μ_u, d=μ_d, s=μ_s), A=(u=A_u, d=A_d, s=A_s))
    thermo_params = (T=T, Φ=phi, Φbar=phibar, ξ=xi)
    K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

    return (quark_params=quark_params, thermo_params=thermo_params, K_coeffs=K_coeffs)
end

function omega_tensor_direct_sigma(process::Symbol; params, p_nodes::Int, angle_nodes::Int, phi_nodes::Int, n_sigma_points::Int)
    quark_params = params.quark_params
    thermo_params = params.thermo_params
    K_coeffs = params.K_coeffs

    p_grid, p_w = gauleg(0.0, Λ_inv_fm, p_nodes)
    cos_grid, cos_w = gauleg(0.0, 1.0, angle_nodes)
    phi_grid, phi_w = gauleg(0.0, Float64(π), phi_nodes)

    pi_sym, pj_sym, _, _ = ASR.parse_particles_from_process(process)
    mi = ASR.get_mass(pi_sym, quark_params)
    mj = ASR.get_mass(pj_sym, quark_params)
    μi = ASR.get_mu(pi_sym, quark_params)
    μj = ASR.get_mu(pj_sym, quark_params)

    T = thermo_params.T
    Φ = thermo_params.Φ
    Φbar = thermo_params.Φbar
    ξ = hasproperty(thermo_params, :ξ) ? thermo_params.ξ : 0.0

    # Use the same density definition as the production function, with the *same* grids.
    ρ_i = ASR.number_density(pi_sym, mi, μi, T, Φ, Φbar, ξ; p_nodes=p_nodes, angle_nodes=angle_nodes, p_grid=p_grid, p_w=p_w, cos_grid=cos_grid, cos_w=cos_w)
    ρ_j = ASR.number_density(pj_sym, mj, μj, T, Φ, Φbar, ξ; p_nodes=p_nodes, angle_nodes=angle_nodes, p_grid=p_grid, p_w=p_w, cos_grid=cos_grid, cos_w=cos_w)
    if ρ_i == 0.0 || ρ_j == 0.0
        return 0.0
    end

    prefactor = (ASR.DQ^2) / (4.0 * π^5 * ρ_i * ρ_j)

    acc = 0.0
    for (p_i, w_pi) in zip(p_grid, p_w)
        Ei = ASR.energy_from_p(p_i, mi)
        for (p_j, w_pj) in zip(p_grid, p_w)
            Ej = ASR.energy_from_p(p_j, mj)
            for (cθi, w_cθi) in zip(cos_grid, cos_w)
                sθi = sqrt(max(1.0 - cθi * cθi, 0.0))
                for (cθj, w_cθj) in zip(cos_grid, cos_w)
                    sθj = sqrt(max(1.0 - cθj * cθj, 0.0))
                    f_i = ASR.distribution_with_anisotropy(pi_sym, p_i, mi, μi, T, Φ, Φbar, ξ, cθi)
                    f_j = ASR.distribution_with_anisotropy(pj_sym, p_j, mj, μj, T, Φ, Φbar, ξ, cθj)
                    if f_i == 0.0 || f_j == 0.0
                        continue
                    end
                    for (φ, wφ) in zip(phi_grid, phi_w)
                        cosΘ = cθi * cθj + sθi * sθj * cos(φ)
                        p_dot = Ei * Ej - p_i * p_j * cosΘ
                        s = mi^2 + mj^2 + 2.0 * p_dot
                        if s <= (mi + mj)^2
                            continue
                        end
                        v_rel_num = p_dot^2 - mi^2 * mj^2
                        v_rel = v_rel_num > 0 ? sqrt(v_rel_num) / (Ei * Ej) : 0.0
                        if v_rel == 0.0
                            continue
                        end
                        σ = total_cross_section(process, s, quark_params, thermo_params, K_coeffs; n_points=n_sigma_points)
                        acc += w_pi * w_pj * w_cθi * w_cθj * wφ * (p_i^2) * (p_j^2) * f_i * f_j * v_rel * σ
                    end
                end
            end
        end
    end

    return prefactor * acc
end

function omega_tensor_cached(process::Symbol; params, p_nodes::Int, angle_nodes::Int, phi_nodes::Int, n_sigma_points::Int,
    N_sigma_grid::Int, design_p_nodes::Int, design_angle_nodes::Int, design_phi_nodes::Int)

    quark_params = params.quark_params
    thermo_params = params.thermo_params
    K_coeffs = params.K_coeffs

    cache = ASR.build_w0cdf_pchip_cache(
        process,
        quark_params,
        thermo_params,
        K_coeffs;
        N=N_sigma_grid,
        design_p_nodes=design_p_nodes,
        design_angle_nodes=design_angle_nodes,
        design_phi_nodes=design_phi_nodes,
        n_sigma_points=n_sigma_points,
    )

    ω = ASR.average_scattering_rate(process, quark_params, thermo_params, K_coeffs;
        p_nodes=p_nodes,
        angle_nodes=angle_nodes,
        phi_nodes=phi_nodes,
        cs_cache=cache,
        n_sigma_points=n_sigma_points,
    )

    return (value=ω, cache_points=length(cache.s_vals))
end

function main()
    process = Symbol(get(ENV, "PROCESS", "udbar_to_udbar"))

    params = build_params(
        T_MeV=parse(Float64, get(ENV, "T_MEV", "150.0")),
        muB_MeV=parse(Float64, get(ENV, "MUB_MEV", "800.0")),
        phi=parse(Float64, get(ENV, "PHI", "0.5")),
        phibar=parse(Float64, get(ENV, "PHIBAR", "0.5")),
        xi=parse(Float64, get(ENV, "XI", "0.0")),
    )

    # Reference budget (direct σ)
    ref_p = parse(Int, get(ENV, "REF_P", "8"))
    ref_a = parse(Int, get(ENV, "REF_ANGLE", "4"))
    ref_phi = parse(Int, get(ENV, "REF_PHI", "6"))
    ref_nσ = parse(Int, get(ENV, "REF_N_SIGMA", "8"))

    @printf("process=%s\n", string(process))
    @printf("ref (direct σ): p=%d angle=%d phi=%d n_sigma_points=%d\n", ref_p, ref_a, ref_phi, ref_nσ)

    t_ref = @elapsed ω_ref = omega_tensor_direct_sigma(process; params=params, p_nodes=ref_p, angle_nodes=ref_a, phi_nodes=ref_phi, n_sigma_points=ref_nσ)
    @printf("ref ω=%.6e  time=%.3fs\n", ω_ref, t_ref)

    # Sweep: which dimension buys the most accuracy (relative to the above ref)
    p_list = [3, 4, 5, 6, 8, 10]
    a_list = [2, 3, 4, 6]
    phi_list = [2, 3, 4, 6, 8]
    nσ_list = [4, 6, 8, 12, 16]

    anchor_p = parse(Int, get(ENV, "ANCHOR_P", "6"))
    anchor_a = parse(Int, get(ENV, "ANCHOR_ANGLE", "2"))
    anchor_phi = parse(Int, get(ENV, "ANCHOR_PHI", "4"))
    anchor_nσ = parse(Int, get(ENV, "ANCHOR_N_SIGMA", "6"))

    @printf("\nAnchor for sweeps: p=%d angle=%d phi=%d n_sigma_points=%d\n", anchor_p, anchor_a, anchor_phi, anchor_nσ)

    function run_case(label::AbstractString; p::Int, a::Int, ph::Int, nσ::Int)
        t = @elapsed ω = omega_tensor_direct_sigma(process; params=params, p_nodes=p, angle_nodes=a, phi_nodes=ph, n_sigma_points=nσ)
        e = relerr(ω, ω_ref)
        @printf("%-18s p=%-2d a=%-2d phi=%-2d nσ=%-2d | ω=%.6e relerr_vs_ref=%.3e time=%.3fs\n", label, p, a, ph, nσ, ω, e, t)
        return nothing
    end

    println("\nSweep p_nodes (direct σ):")
    for p in p_list
        run_case("p_sweep"; p=p, a=anchor_a, ph=anchor_phi, nσ=anchor_nσ)
    end

    println("\nSweep angle_nodes (direct σ):")
    for a in a_list
        run_case("angle_sweep"; p=anchor_p, a=a, ph=anchor_phi, nσ=anchor_nσ)
    end

    println("\nSweep phi_nodes (direct σ):")
    for ph in phi_list
        run_case("phi_sweep"; p=anchor_p, a=anchor_a, ph=ph, nσ=anchor_nσ)
    end

    println("\nSweep n_sigma_points (direct σ):")
    for nσ in nσ_list
        run_case("nσ_sweep"; p=anchor_p, a=anchor_a, ph=anchor_phi, nσ=nσ)
    end

    # Compare cache strategies against direct σ for a fixed node budget
    test_p = parse(Int, get(ENV, "TEST_P", string(anchor_p)))
    test_a = parse(Int, get(ENV, "TEST_ANGLE", string(anchor_a)))
    test_phi = parse(Int, get(ENV, "TEST_PHI", string(anchor_phi)))
    test_nσ = parse(Int, get(ENV, "TEST_N_SIGMA", string(anchor_nσ)))

    @printf("\nCompare σ strategy at fixed budget: p=%d angle=%d phi=%d n_sigma_points=%d\n", test_p, test_a, test_phi, test_nσ)
    t_dir = @elapsed ω_dir = omega_tensor_direct_sigma(process; params=params, p_nodes=test_p, angle_nodes=test_a, phi_nodes=test_phi, n_sigma_points=test_nσ)
    @printf("direct: ω=%.6e  time=%.3fs\n", ω_dir, t_dir)

    N_sigma_grid = parse(Int, get(ENV, "N_COARSE", "60"))
    design_p_nodes = parse(Int, get(ENV, "DESIGN_P", "14"))
    design_angle_nodes = parse(Int, get(ENV, "DESIGN_ANGLE", "4"))
    design_phi_nodes = parse(Int, get(ENV, "DESIGN_PHI", "8"))

    t = @elapsed r_interp = omega_tensor_cached(process; params=params, p_nodes=test_p, angle_nodes=test_a, phi_nodes=test_phi, n_sigma_points=test_nσ,
        N_sigma_grid=N_sigma_grid, design_p_nodes=design_p_nodes, design_angle_nodes=design_angle_nodes, design_phi_nodes=design_phi_nodes)
    @printf("w0cdf+pchip: ω=%.6e relerr_vs_direct=%.3e time=%.3fs cache_points=%d\n", r_interp.value, relerr(r_interp.value, ω_dir), t, r_interp.cache_points)

    println("\nDone.")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
