#!/usr/bin/env julia

raw"""Benchmark *one full scan point* (T, μB) as done in the real scan script.

It measures a realistic per-point cost breakdown:
- PNJL equilibrium solve (gap + Polyakov loop)
- Build A/K coeffs + densities
- Precompute σ(s) caches for all REQUIRED_PROCESSES
- Compute relaxation times τ via RelaxationTime.relaxation_times

Configure via environment variables (all optional):
  T_MEV, MUB_MEV, XI
  P_NUM, T_NUM, MAX_ITER
  TAU_P_NODES, TAU_ANGLE_NODES, TAU_PHI_NODES, TAU_N_SIGMA
  SIGMA_GRID_N, P_CUT_FACTOR
  SIGMA_COMPUTE_MISSING (0/1), SIGMA_CACHE_RTOL, SIGMA_CACHE_MAX_REFINE
  N_REPEAT (default 3)

Example (PowerShell):
  $env:T_MEV=150; $env:MUB_MEV=800; $env:N_REPEAT=3
  julia --project=. scripts/debug/relaxtime/bench_scan_point.jl
"""

using Printf
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

# Reuse the exact implementations used by the production scan script
include(joinpath(PROJECT_ROOT, "scripts", "relaxtime", "scan_relaxation_times_vs_T.jl"))

# Symbols are defined by the included scan script
# - solve_equilibrium_mu
# - cached_nodes, calculate_mass_vec, calculate_number_densities, DEFAULT_MU_GUESS
# - relaxation_times, REQUIRED_PROCESSES
# - build_K_coeffs, densities_from_equilibrium, build_sigma_caches

function getenv_float(name::AbstractString, default::Float64)
    v = get(ENV, name, "")
    return isempty(v) ? default : parse(Float64, v)
end

function getenv_int(name::AbstractString, default::Int)
    v = get(ENV, name, "")
    return isempty(v) ? default : parse(Int, v)
end

function getenv_bool01(name::AbstractString, default::Bool)
    v = get(ENV, name, "")
    if isempty(v)
        return default
    end
    x = parse(Int, v)
    return x != 0
end

function timed_block(label::AbstractString, f::Function)
    GC.gc()
    t = @timed f()
    @printf("%-22s  time=%8.3f s  alloc=%8.1f MB  gctime=%6.3f s\n",
        label,
        t.time,
        t.bytes / 1024^2,
        t.gctime,
    )
    return t.value
end

function run_one_point(; T_mev::Float64, muB_mev::Float64, xi::Float64,
    p_num::Int, t_num::Int, max_iter::Int,
    tau_p_nodes::Int, tau_angle_nodes::Int, tau_phi_nodes::Int, tau_n_sigma_points::Int,
    sigma_grid_n::Int, p_cut_factor::Float64,
    sigma_cache_rtol::Float64, sigma_cache_max_refine::Int, sigma_compute_missing::Bool)

    muq_mev = muB_mev / 3.0
    muq_fm = muq_mev / ħc_MeV_fm
    T_fm = T_mev / ħc_MeV_fm

    base = timed_block("equilibrium", () -> solve_equilibrium_mu(
        T_fm,
        muq_fm;
        xi=xi,
        seed_state=DEFAULT_MU_GUESS,
        p_num=p_num,
        t_num=t_num,
        iterations=max_iter,
    ))

    if !Bool(base.converged)
        @printf("equilibrium NOT converged (iters=%d, residual=%.3e)\n", Int(base.iterations), Float64(base.residual_norm))
        return nothing
    end

    Φ = Float64(base.x_state[4])
    Φbar = Float64(base.x_state[5])
    masses_vec = calculate_mass_vec(SVector{3}(base.x_state[1], base.x_state[2], base.x_state[3]))
    masses = (u=Float64(masses_vec[1]), d=Float64(masses_vec[2]), s=Float64(masses_vec[3]))

    thermo_params = (T=Float64(T_fm), Φ=Φ, Φbar=Φbar, ξ=Float64(xi))
    quark_params_basic = (m=masses, μ=(u=Float64(muq_fm), d=Float64(muq_fm), s=Float64(muq_fm)))

    ktmp = timed_block("A+K_coeffs", () -> build_K_coeffs(Float64(T_fm), Float64(muq_fm), masses, Φ, Φbar))
    quark_params = (m=quark_params_basic.m, μ=quark_params_basic.μ, A=ktmp.A_vals)
    K_coeffs = ktmp.K_coeffs

    thermal_nodes = cached_nodes(p_num, t_num)
    densities = timed_block("densities", () -> densities_from_equilibrium(base.x_state, base.mu_vec, T_fm, thermal_nodes, Float64(xi)))

    cs_caches = timed_block("sigma_caches", () -> build_sigma_caches(
        REQUIRED_PROCESSES,
        quark_params,
        thermo_params,
        K_coeffs;
        n_sigma_points=tau_n_sigma_points,
        sigma_grid_n=sigma_grid_n,
        p_cut_factor=p_cut_factor,
        sigma_cache_rtol=sigma_cache_rtol,
        sigma_cache_max_refine=sigma_cache_max_refine,
        compute_missing=sigma_compute_missing,
    ))

    tau_res = timed_block("relaxation_times", () -> relaxation_times(
        quark_params,
        thermo_params,
        K_coeffs;
        densities=densities,
        cs_caches=cs_caches,
        p_nodes=tau_p_nodes,
        angle_nodes=tau_angle_nodes,
        phi_nodes=tau_phi_nodes,
        n_sigma_points=tau_n_sigma_points,
    ))

    tau = tau_res.tau
    @printf("tau: u=%.3e  s=%.3e  ubar=%.3e  sbar=%.3e\n", tau.u, tau.s, tau.ubar, tau.sbar)

    # lightweight cache size summary (helps diagnose memory blow-ups)
    total_pts = 0
    max_pts = 0
    for (proc, cache) in cs_caches
        n = length(cache.s_vals)
        total_pts += n
        max_pts = max(max_pts, n)
    end
    @printf("cache points: processes=%d  total_s_points=%d  max_per_process=%d  compute_missing=%s\n",
        length(REQUIRED_PROCESSES), total_pts, max_pts, string(sigma_compute_missing))

    return nothing
end

function main()
    T_mev = getenv_float("T_MEV", 150.0)
    muB_mev = getenv_float("MUB_MEV", 800.0)
    xi = getenv_float("XI", 0.0)

    p_num = getenv_int("P_NUM", 12)
    t_num = getenv_int("T_NUM", 6)
    max_iter = getenv_int("MAX_ITER", 40)

    tau_p_nodes = getenv_int("TAU_P_NODES", 6)
    tau_angle_nodes = getenv_int("TAU_ANGLE_NODES", 2)
    tau_phi_nodes = getenv_int("TAU_PHI_NODES", 6)
    tau_n_sigma_points = getenv_int("TAU_N_SIGMA", 6)

    sigma_grid_n = getenv_int("SIGMA_GRID_N", 18)
    p_cut_factor = getenv_float("P_CUT_FACTOR", 8.0)

    sigma_compute_missing = getenv_bool01("SIGMA_COMPUTE_MISSING", false)
    sigma_cache_rtol = getenv_float("SIGMA_CACHE_RTOL", 5e-2)
    sigma_cache_max_refine = getenv_int("SIGMA_CACHE_MAX_REFINE", 6)

    n_repeat = getenv_int("N_REPEAT", 3)

    @printf("point: T=%.1f MeV  muB=%.1f MeV  xi=%.3g\n", T_mev, muB_mev, xi)
    @printf("equilibrium: p_num=%d  t_num=%d  max_iter=%d\n", p_num, t_num, max_iter)
    @printf("tau nodes: p=%d  angle=%d  phi=%d  n_sigma=%d\n", tau_p_nodes, tau_angle_nodes, tau_phi_nodes, tau_n_sigma_points)
    @printf("sigma grid: n=%d  p_cut_factor=%.2f  compute_missing=%s  rtol=%.2e  max_refine=%d\n",
        sigma_grid_n, p_cut_factor, string(sigma_compute_missing), sigma_cache_rtol, sigma_cache_max_refine)

    # Warmup for JIT (do not trust first run)
    @printf("\n-- warmup (JIT) --\n")
    run_one_point(
        T_mev=T_mev,
        muB_mev=muB_mev,
        xi=xi,
        p_num=min(p_num, 8),
        t_num=min(t_num, 4),
        max_iter=max_iter,
        tau_p_nodes=min(tau_p_nodes, 4),
        tau_angle_nodes=tau_angle_nodes,
        tau_phi_nodes=min(tau_phi_nodes, 4),
        tau_n_sigma_points=min(tau_n_sigma_points, 4),
        sigma_grid_n=min(sigma_grid_n, 12),
        p_cut_factor=p_cut_factor,
        sigma_cache_rtol=sigma_cache_rtol,
        sigma_cache_max_refine=sigma_cache_max_refine,
        sigma_compute_missing=sigma_compute_missing,
    )

    @printf("\n-- measured runs --\n")
    for i in 1:n_repeat
        @printf("\nrun %d/%d\n", i, n_repeat)
        run_one_point(
            T_mev=T_mev,
            muB_mev=muB_mev,
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            max_iter=max_iter,
            tau_p_nodes=tau_p_nodes,
            tau_angle_nodes=tau_angle_nodes,
            tau_phi_nodes=tau_phi_nodes,
            tau_n_sigma_points=tau_n_sigma_points,
            sigma_grid_n=sigma_grid_n,
            p_cut_factor=p_cut_factor,
            sigma_cache_rtol=sigma_cache_rtol,
            sigma_cache_max_refine=sigma_cache_max_refine,
            sigma_compute_missing=sigma_compute_missing,
        )
    end

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
