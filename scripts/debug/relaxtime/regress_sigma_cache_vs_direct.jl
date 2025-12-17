#!/usr/bin/env julia

using Printf
using Random

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
push!(LOAD_PATH, joinpath(PROJECT_ROOT, "src"))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "AverageScatteringRate.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalCrossSection.jl"))

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5
using .GaussLegendre: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using .OneLoopIntegrals: A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .TotalCrossSection: total_cross_section

const ASR = AverageScatteringRate

"""Build quark/thermo parameters for σ(s) benchmarking.

Note: this is *not* the PNJL gap-solver output; it's a lightweight parameter builder
that provides the fields required by `total_cross_section` (including `A`).
"""
function build_params(; T_MeV=150.0, muB_MeV=800.0,
    mu_s_MeV=0.0, m_u_MeV=300.0, m_d_MeV=300.0, m_s_MeV=500.0,
    phi=0.5, phibar=0.5, xi=0.0)

    T = T_MeV / ħc_MeV_fm
    μ_q = (muB_MeV / 3.0) / ħc_MeV_fm
    μ_u = μ_q
    μ_d = μ_q
    μ_s = mu_s_MeV / ħc_MeV_fm

    m_u = m_u_MeV / ħc_MeV_fm
    m_d = m_d_MeV / ħc_MeV_fm
    m_s = m_s_MeV / ħc_MeV_fm

    nodes_p = DEFAULT_MOMENTUM_NODES
    weights_p = DEFAULT_MOMENTUM_WEIGHTS

    A_u = A(m_u, μ_u, T, phi, phibar, nodes_p, weights_p)
    A_d = A(m_d, μ_d, T, phi, phibar, nodes_p, weights_p)
    A_s = A(m_s, μ_s, T, phi, phibar, nodes_p, weights_p)

    G_u = calculate_G_from_A(A_u, m_u)
    G_s = calculate_G_from_A(A_s, m_s)

    quark_params = (m=(u=m_u, d=m_d, s=m_s), μ=(u=μ_u, d=μ_d, s=μ_s), A=(u=A_u, d=A_d, s=A_s))
    thermo_params = (T=T, Φ=phi, Φbar=phibar, ξ=xi)
    K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

    return (quark_params=quark_params, thermo_params=thermo_params, K_coeffs=K_coeffs)
end

@inline function relerr(a::Float64, b::Float64)
    return abs(a - b) / max(1e-12, abs(b))
end

function main()
    process = Symbol(get(ENV, "PROCESS", "udbar_to_udbar"))
    rng_seed = parse(Int, get(ENV, "SEED", "1"))

    # Random sampling range is configured in √s (MeV) for ergonomics.
    sqrt_s_min_MeV = parse(Float64, get(ENV, "SQRT_S_MIN_MEV", "250.0"))
    sqrt_s_max_MeV = parse(Float64, get(ENV, "SQRT_S_MAX_MEV", "1100.0"))

    n_samples = parse(Int, get(ENV, "N_SAMPLES", "24"))
    n_coarse = parse(Int, get(ENV, "N_COARSE", "12"))
    n_points = parse(Int, get(ENV, "N_POINTS", "6"))
    rtol = parse(Float64, get(ENV, "RTOL", "1e-2"))
    max_refine = parse(Int, get(ENV, "MAX_REFINE", "12"))

    params = build_params(
        T_MeV=parse(Float64, get(ENV, "T_MEV", "150.0")),
        muB_MeV=parse(Float64, get(ENV, "MUB_MEV", "800.0")),
        mu_s_MeV=parse(Float64, get(ENV, "MU_S_MEV", "0.0")),
        m_u_MeV=parse(Float64, get(ENV, "M_U_MEV", "300.0")),
        m_d_MeV=parse(Float64, get(ENV, "M_D_MEV", "300.0")),
        m_s_MeV=parse(Float64, get(ENV, "M_S_MEV", "500.0")),
        phi=parse(Float64, get(ENV, "PHI", "0.5")),
        phibar=parse(Float64, get(ENV, "PHIBAR", "0.5")),
        xi=parse(Float64, get(ENV, "XI", "0.0")),
    )

    s_min = (sqrt_s_min_MeV / ħc_MeV_fm)^2
    s_max = (sqrt_s_max_MeV / ħc_MeV_fm)^2
    s_min < s_max || error("require SQRT_S_MIN_MEV < SQRT_S_MAX_MEV")

    rng = MersenneTwister(rng_seed)
    s_samples = [s_min + rand(rng) * (s_max - s_min) for _ in 1:n_samples]

    # Shared coarse grid (precompute endpoints for both cache strategies)
    s_grid = collect(range(s_min, s_max; length=n_coarse))

    function build_interp_cache()
        cache = ASR.CrossSectionCache(process; compute_missing=false)
        ASR.precompute_cross_section!(cache, s_grid, params.quark_params, params.thermo_params, params.K_coeffs; n_points=n_points)
        return cache
    end

    function build_adaptive_cache()
        cache = ASR.CrossSectionCache(process; compute_missing=true, rtol=rtol, max_refine=max_refine)
        ASR.precompute_cross_section!(cache, s_grid, params.quark_params, params.thermo_params, params.K_coeffs; n_points=n_points)
        return cache
    end

    function summarize_errors(label::AbstractString, s_vals::Vector{Float64}, σ_pred::Vector{Float64}, σ_true::Vector{Float64})
        max_rel = 0.0
        max_abs = 0.0
        worst = (idx=0, s=0.0, σp=0.0, σt=0.0, rel=0.0, abs=0.0)
        for i in eachindex(s_vals)
            σp = σ_pred[i]
            σt = σ_true[i]
            re = relerr(σp, σt)
            ae = abs(σp - σt)
            max_rel = max(max_rel, re)
            max_abs = max(max_abs, ae)
            if re >= worst.rel
                worst = (idx=i, s=s_vals[i], σp=σp, σt=σt, rel=re, abs=ae)
            end
        end

        @printf("error[%s]: max_rel=%.3e  max_abs=%.3e\n", label, max_rel, max_abs)
        @printf("worst[%s](idx=%d): s=%.6e  sqrt_s=%.3f MeV  σ_pred=%.6e  σ_true=%.6e  rel=%.3e  abs=%.3e\n",
            label,
            worst.idx,
            worst.s,
            sqrt(worst.s) * ħc_MeV_fm,
            worst.σp,
            worst.σt,
            worst.rel,
            worst.abs,
        )
        return nothing
    end

    # 0) Direct (reference)
    t0 = time_ns()
    σ_direct = [total_cross_section(process, s, params.quark_params, params.thermo_params, params.K_coeffs; n_points=n_points) for s in s_samples]
    t_direct_s = (time_ns() - t0) / 1e9

    # 1) Interpolation-only cache (no refinement, no missing compute)
    cache_interp = build_interp_cache()
    t0 = time_ns()
    σ_interp = [ASR.get_sigma(cache_interp, s, params.quark_params, params.thermo_params, params.K_coeffs; n_points=n_points) for s in s_samples]
    t_interp_s = (time_ns() - t0) / 1e9

    # 2) Adaptive cache (local midpoint refinement + missing compute)
    cache_adapt = build_adaptive_cache()
    n0 = length(cache_adapt.s_vals)
    t0 = time_ns()
    σ_adapt = [ASR.get_sigma(cache_adapt, s, params.quark_params, params.thermo_params, params.K_coeffs; n_points=n_points) for s in s_samples]
    t_adapt_s = (time_ns() - t0) / 1e9
    n1 = length(cache_adapt.s_vals)
    added = n1 - n0

    @printf("process=%s\n", string(process))
    @printf("range: sqrt_s=[%.1f, %.1f] MeV  (s=[%.6e, %.6e])\n", sqrt_s_min_MeV, sqrt_s_max_MeV, s_min, s_max)
    @printf("config: N_SAMPLES=%d  N_COARSE=%d  N_POINTS=%d  SEED=%d\n", n_samples, n_coarse, n_points, rng_seed)
    @printf("adaptive: RTOL=%.3e  MAX_REFINE=%d\n", rtol, max_refine)
    @printf("cache(interp): points=%d (precomputed)  compute_missing=%s\n", length(cache_interp.s_vals), string(cache_interp.compute_missing))
    @printf("cache(adapt):  points=%d (precomputed) -> %d (after queries)  added=%d\n", n0, n1, added)

    @printf("timing: direct=%.3fs  interp=%.3fs  adapt=%.3fs\n", t_direct_s, t_interp_s, t_adapt_s)
    @printf("speedup vs direct: interp=%.2fx  adapt=%.2fx\n", t_direct_s / max(1e-12, t_interp_s), t_direct_s / max(1e-12, t_adapt_s))

    summarize_errors("interp", s_samples, σ_interp, σ_direct)
    summarize_errors("adapt", s_samples, σ_adapt, σ_direct)

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
