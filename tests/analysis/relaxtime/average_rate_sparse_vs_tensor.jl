#!/usr/bin/env julia

using Printf

println("Starting average_rate_sparse_vs_tensor.jl ...")
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

# --------------------------- parameter helpers ---------------------------

function build_params(; T_MeV=150.0, mu_u_MeV=0.0, mu_d_MeV=0.0, mu_s_MeV=0.0,
    m_u_MeV=300.0, m_d_MeV=300.0, m_s_MeV=500.0, phi=0.5, phibar=0.5, xi=0.0)
    T = T_MeV / ħc_MeV_fm
    μ_u = mu_u_MeV / ħc_MeV_fm
    μ_d = mu_d_MeV / ħc_MeV_fm
    μ_s = mu_s_MeV / ħc_MeV_fm
    m_u = m_u_MeV / ħc_MeV_fm
    m_d = m_d_MeV / ħc_MeV_fm
    m_s = m_s_MeV / ħc_MeV_fm

    # reuse default Gauss-Legendre nodes for A integrals
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

# --------------------------- integrand ---------------------------

function omega_core(process::Symbol, params, cache::ASR.CrossSectionCache,
    p_i::Float64, p_j::Float64, cθi::Float64, cθj::Float64, φ::Float64; n_sigma_points::Int)
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
        return 0.0
    end

    cosΘ = cθi * cθj + sθi * sθj * cos(φ)
    p_dot = Ei * Ej - p_i * p_j * cosΘ
    s = mi^2 + mj^2 + 2.0 * p_dot
    if s <= (mi + mj)^2
        return 0.0
    end

    v_rel_num = p_dot^2 - mi^2 * mj^2
    v_rel = v_rel_num > 0 ? sqrt(v_rel_num) / (Ei * Ej) : 0.0
    if v_rel == 0.0
        return 0.0
    end

    σ = ASR.get_sigma(cache, s, params.quark_params, params.thermo_params, params.K_coeffs; n_points=n_sigma_points)
    return (p_i ^ 2) * (p_j ^ 2) * f_i * f_j * v_rel * σ
end

# --------------------------- tensor quadrature ---------------------------

function tensor_average_rate(process::Symbol; params, cache, p_nodes::Int, angle_nodes::Int, phi_nodes::Int, n_sigma_points::Int)
    p_grid, p_w = gauleg(0.0, Λ_inv_fm, p_nodes)
    cos_grid, cos_w = gauleg(0.0, 1.0, angle_nodes)
    phi_grid, phi_w = gauleg(0.0, Float64(π), phi_nodes)

    # densities use the module defaults
    pi_sym, pj_sym, _, _ = ASR.parse_particles_from_process(process)
    ρ_i = ASR.number_density(pi_sym, ASR.get_mass(pi_sym, params.quark_params), ASR.get_mu(pi_sym, params.quark_params), params.thermo_params.T, params.thermo_params.Φ, params.thermo_params.Φbar, hasproperty(params.thermo_params, :ξ) ? params.thermo_params.ξ : 0.0)
    ρ_j = ASR.number_density(pj_sym, ASR.get_mass(pj_sym, params.quark_params), ASR.get_mu(pj_sym, params.quark_params), params.thermo_params.T, params.thermo_params.Φ, params.thermo_params.Φbar, hasproperty(params.thermo_params, :ξ) ? params.thermo_params.ξ : 0.0)
    prefactor = (ASR.DQ ^ 2) / (4.0 * π ^ 5 * ρ_i * ρ_j)

    acc = 0.0
    evals = 0
    for (p_i, w_pi) in zip(p_grid, p_w)
        for (p_j, w_pj) in zip(p_grid, p_w)
            for (cθi, w_cθi) in zip(cos_grid, cos_w)
                for (cθj, w_cθj) in zip(cos_grid, cos_w)
                    for (φ, wφ) in zip(phi_grid, phi_w)
                        val = omega_core(process, params, cache, p_i, p_j, cθi, cθj, φ; n_sigma_points=n_sigma_points)
                        if val != 0.0
                            acc += w_pi * w_pj * w_cθi * w_cθj * wφ * val
                        end
                        evals += 1
                    end
                end
            end
        end
    end
    return (value = prefactor * acc, evals = evals)
end

# --------------------------- sparse grid (Smolyak) ---------------------------

@inline function level_to_order(level::Int)
    return 2 * level - 1
end

function smolyak_multi_indices(dim::Int, q::Int)
    idx = Vector{Int}(undef, dim)
    results = Vector{Vector{Int}}()
    function rec(pos::Int, current_sum::Int)
        if pos > dim
            push!(results, copy(idx))
            return
        end
        max_i = q - current_sum - (dim - pos)
        for i in 1:max_i
            idx[pos] = i
            rec(pos + 1, current_sum + i)
        end
    end
    rec(1, 0)
    return results
end

function smolyak_average_rate(process::Symbol; params, cache, level::Int, n_sigma_points::Int)
    dim = 5
    q = level + dim - 1
    ranges = ((0.0, Λ_inv_fm), (0.0, Λ_inv_fm), (0.0, 1.0), (0.0, 1.0), (0.0, Float64(π)))

    # densities with defaults to keep consistent with tensor version
    pi_sym, pj_sym, _, _ = ASR.parse_particles_from_process(process)
    ρ_i = ASR.number_density(pi_sym, ASR.get_mass(pi_sym, params.quark_params), ASR.get_mu(pi_sym, params.quark_params), params.thermo_params.T, params.thermo_params.Φ, params.thermo_params.Φbar, hasproperty(params.thermo_params, :ξ) ? params.thermo_params.ξ : 0.0)
    ρ_j = ASR.number_density(pj_sym, ASR.get_mass(pj_sym, params.quark_params), ASR.get_mu(pj_sym, params.quark_params), params.thermo_params.T, params.thermo_params.Φ, params.thermo_params.Φbar, hasproperty(params.thermo_params, :ξ) ? params.thermo_params.ξ : 0.0)
    prefactor = (ASR.DQ ^ 2) / (4.0 * π ^ 5 * ρ_i * ρ_j)

    # cache 1D rules per (dim_index, level)
    rule_cache = Dict{Tuple{Int,Int},Tuple{Vector{Float64},Vector{Float64}}}()
    function get_rule(d::Int, l::Int)
        key = (d, l)
        if haskey(rule_cache, key)
            return rule_cache[key]
        end
        a, b = ranges[d]
        n = level_to_order(l)
        nodes, weights = gauleg(a, b, n)
        rule_cache[key] = (nodes, weights)
        return rule_cache[key]
    end

    multi_indices = smolyak_multi_indices(dim, q)
    acc = 0.0
    evals = 0

    for idx in multi_indices
        s = sum(idx)
        coeff = (-1) ^ (q - s) * binomial(dim - 1, q - s)
        rules = [get_rule(d, idx[d]) for d in 1:dim]
        for (p_i, w1) in zip(rules[1][1], rules[1][2])
            for (p_j, w2) in zip(rules[2][1], rules[2][2])
                for (cθi, w3) in zip(rules[3][1], rules[3][2])
                    for (cθj, w4) in zip(rules[4][1], rules[4][2])
                        for (φ, w5) in zip(rules[5][1], rules[5][2])
                            w = coeff * w1 * w2 * w3 * w4 * w5
                            val = omega_core(process, params, cache, p_i, p_j, cθi, cθj, φ; n_sigma_points=n_sigma_points)
                            if val != 0.0
                                acc += w * val
                            end
                            evals += 1
                        end
                    end
                end
            end
        end
    end

    return (value = prefactor * acc, evals = evals)
end

# --------------------------- driver ---------------------------

function compare_methods()
    params = build_params()
    process = :ssbar_to_uubar
    cache = ASR.CrossSectionCache(process)
    n_sigma_points = 16

    # reference: denser tensor grid
    ref_nodes = (p=10, angle=6, phi=8)
    @printf("Reference (tensor) nodes p=%d angle=%d phi=%d ...\n", ref_nodes.p, ref_nodes.angle, ref_nodes.phi)
    t_ref = @elapsed begin
        global ref = tensor_average_rate(process; params=params, cache=cache, p_nodes=ref_nodes.p, angle_nodes=ref_nodes.angle, phi_nodes=ref_nodes.phi, n_sigma_points=n_sigma_points)
    end
    per_eval_ref = t_ref / ref.evals
    @printf("Reference ω = %.6e (evals=%d), time=%.3fs, per-eval=%.3es\n", ref.value, ref.evals, t_ref, per_eval_ref)

    # matched-evals tensor (close to sparse eval count)
    matched_nodes = (p=4, angle=3, phi=3) # evals = p^2 * angle^2 * phi = 144
    t_matched = @elapsed begin
        matched = tensor_average_rate(process; params=params, cache=cache, p_nodes=matched_nodes.p, angle_nodes=matched_nodes.angle, phi_nodes=matched_nodes.phi, n_sigma_points=n_sigma_points)
    end
    rel_err_matched = abs(matched.value - ref.value) / abs(ref.value)
    per_eval_matched = t_matched / matched.evals
    @printf("Matched tensor ω = %.6e (evals=%d), rel_err vs ref = %.3e, time=%.3fs, per-eval=%.3es\n", matched.value, matched.evals, rel_err_matched, t_matched, per_eval_matched)

    # sparse grid (Smolyak)
    level = 3
    t_sparse = @elapsed begin
        sparse = smolyak_average_rate(process; params=params, cache=cache, level=level, n_sigma_points=n_sigma_points)
    end
    rel_err_sparse = abs(sparse.value - ref.value) / abs(ref.value)
    per_eval_sparse = t_sparse / sparse.evals
    @printf("Sparse grid (level=%d) ω = %.6e (evals=%d), rel_err vs ref = %.3e, time=%.3fs, per-eval=%.3es\n", level, sparse.value, sparse.evals, rel_err_sparse, t_sparse, per_eval_sparse)

    # crude comparison of evals
    @printf("Eval counts: tensor matched %d vs sparse %d (level=%d)\n", matched.evals, sparse.evals, level)

    # ---------------- sweeps: which dimension needs more nodes? ----------------
    p_anchor, angle_anchor, phi_anchor = 6, 4, 6
    p_list = [3, 4, 5, 6, 8, 10]
    angle_list = [2, 3, 4, 6, 8]
    phi_list = [2, 3, 4, 6, 8]

    @printf("\nSweep p-nodes (angle=%d, phi=%d):\n", angle_anchor, phi_anchor)
    for p_n in p_list
        res = tensor_average_rate(process; params=params, cache=cache, p_nodes=p_n, angle_nodes=angle_anchor, phi_nodes=phi_anchor, n_sigma_points=n_sigma_points)
        rel_err = abs(res.value - ref.value) / abs(ref.value)
        @printf("p=%2d angle=%2d phi=%2d | ω=%.6e rel_err=%.3e evals=%5d\n", p_n, angle_anchor, phi_anchor, res.value, rel_err, res.evals)
    end

    @printf("\nSweep angle-nodes (p=%d, phi=%d):\n", p_anchor, phi_anchor)
    for a_n in angle_list
        res = tensor_average_rate(process; params=params, cache=cache, p_nodes=p_anchor, angle_nodes=a_n, phi_nodes=phi_anchor, n_sigma_points=n_sigma_points)
        rel_err = abs(res.value - ref.value) / abs(ref.value)
        @printf("p=%2d angle=%2d phi=%2d | ω=%.6e rel_err=%.3e evals=%5d\n", p_anchor, a_n, phi_anchor, res.value, rel_err, res.evals)
    end

    @printf("\nSweep phi-nodes (p=%d, angle=%d):\n", p_anchor, angle_anchor)
    for ph_n in phi_list
        res = tensor_average_rate(process; params=params, cache=cache, p_nodes=p_anchor, angle_nodes=angle_anchor, phi_nodes=ph_n, n_sigma_points=n_sigma_points)
        rel_err = abs(res.value - ref.value) / abs(ref.value)
        @printf("p=%2d angle=%2d phi=%2d | ω=%.6e rel_err=%.3e evals=%5d\n", p_anchor, angle_anchor, ph_n, res.value, rel_err, res.evals)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    compare_methods()
end
