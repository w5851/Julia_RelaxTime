"""
单点诊断：检查某个 (T, μB, ξ) 下的密度、平均散射率、截面缓存覆盖范围与最终 τ。

用法：
  julia --project=. scripts/debug/debug_tau_point.jl --T 150 --muB 800 --xi 0

输出：
- PNJL 求解状态、Φ/Φbar、质量
- 数密度（来自 PNJL 的 cached_nodes 积分）
- 每个散射过程的 w0cdf σ(s) cache: [s_min, s_max], [σ_min, σ_max]
- 估计的 s_max(Λ)（由 p_cut=Λ 和 cosΘ=-1 给出）与 cache 覆盖对比
- relaxation_times 的 rates / τ^{-1} / τ
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "RelaxationTime.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "OneLoopIntegrals.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))

using Printf
using StaticArrays

using .Constants_PNJL: ħc_MeV_fm, G_fm2, K_fm5, Λ_inv_fm
using .PNJL: solve, FixedMu, cached_nodes, calculate_number_densities
using .PNJL.Integrals: DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
using .RelaxationTime: relaxation_times, REQUIRED_PROCESSES
using .OneLoopIntegrals: A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings

const RT_ASR = RelaxationTime.AverageScatteringRate

struct Opt
    T_mev::Float64
    muB_mev::Float64
    xi::Float64
    p_num::Int
    t_num::Int
    tau_p_nodes::Int
    tau_angle_nodes::Int
    tau_phi_nodes::Int
    tau_n_sigma_points::Int
    sigma_grid_n::Int
end

function parse_args(args)
    # defaults match scan script typical debug-ish values
    T = 150.0
    muB = 800.0
    xi = 0.0
    p_num = 12
    t_num = 6
    tau_p_nodes = 6
    tau_angle_nodes = 2
    tau_phi_nodes = 6
    tau_n_sigma = 6
    sigma_grid_n = 18

    i = 1
    while i <= length(args)
        a = args[i]
        function needval()
            i == length(args) && error("missing value for $a")
            v = args[i+1]
            i += 1
            return v
        end

        if a == "--T"
            T = parse(Float64, needval())
        elseif a == "--muB"
            muB = parse(Float64, needval())
        elseif a == "--xi"
            xi = parse(Float64, needval())
        elseif a == "--p-num"
            p_num = parse(Int, needval())
        elseif a == "--t-num"
            t_num = parse(Int, needval())
        elseif a == "--tau-p-nodes"
            tau_p_nodes = parse(Int, needval())
        elseif a == "--tau-angle-nodes"
            tau_angle_nodes = parse(Int, needval())
        elseif a == "--tau-phi-nodes"
            tau_phi_nodes = parse(Int, needval())
        elseif a == "--tau-n-sigma"
            tau_n_sigma = parse(Int, needval())
        elseif a == "--sigma-grid-n"
            sigma_grid_n = parse(Int, needval())
        elseif a in ("-h", "--help")
            println("Usage: julia --project=. scripts/debug/debug_tau_point.jl [--T <MeV>] [--muB <MeV>] [--xi <val>] ...")
            exit(0)
        else
            error("unknown arg: $a")
        end
        i += 1
    end

    return Opt(T, muB, xi, p_num, t_num, tau_p_nodes, tau_angle_nodes, tau_phi_nodes, tau_n_sigma, sigma_grid_n)
end

function build_K_coeffs(T_fm::Float64, muq_fm::Float64, masses::NamedTuple, Φ::Float64, Φbar::Float64)
    A_u = A(masses.u, muq_fm, T_fm, Φ, Φbar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    A_s = A(masses.s, muq_fm, T_fm, Φ, Φbar, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS)
    G_u = calculate_G_from_A(A_u, masses.u)
    G_s = calculate_G_from_A(A_s, masses.s)
    return (K_coeffs=calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s), A_vals=(u=A_u, d=A_u, s=A_s))
end

function densities_from_equilibrium(x_state, mu_vec, T_fm, thermal_nodes, xi)
    nd = calculate_number_densities(x_state, mu_vec, T_fm, thermal_nodes, xi)
    return (
        u=Float64(nd.quark[1]),
        d=Float64(nd.quark[2]),
        s=Float64(nd.quark[3]),
        ubar=Float64(nd.antiquark[1]),
        dbar=Float64(nd.antiquark[2]),
        sbar=Float64(nd.antiquark[3]),
    )
end

function safe_total_cross_section(process::Symbol, s::Float64, quark_params, thermo_params, K_coeffs; n_points::Int)
    # Keep it simple here: if it errors, mark NaN.
    try
        σ = RelaxationTime.TotalCrossSection.total_cross_section(process, s, quark_params, thermo_params, K_coeffs; n_points=n_points)
        return isfinite(σ) ? σ : NaN
    catch
        return NaN
    end
end

function build_sigma_caches(quark_params, thermo_params, K_coeffs; n_sigma_points::Int, sigma_grid_n::Int)
    cs_caches = Dict{Symbol,RT_ASR.CrossSectionCache}()
    for process in REQUIRED_PROCESSES
        s_grid = RT_ASR.design_w0cdf_s_grid(process, quark_params, thermo_params; N=sigma_grid_n)
        cache = RT_ASR.CrossSectionCache(process)
        for s in s_grid
            σ = safe_total_cross_section(process, s, quark_params, thermo_params, K_coeffs; n_points=n_sigma_points)
            isfinite(σ) || continue
            RT_ASR.insert_sigma!(cache, s, σ)
        end
        length(cache.s_vals) >= 2 || error("sigma cache too small for $process")
        cs_caches[process] = cache
    end
    return cs_caches
end

@inline function fmt(x)
    return @sprintf("%.6e", x)
end

function main()
    opt = parse_args(ARGS)

    muq_mev = opt.muB_mev / 3.0
    T_fm = opt.T_mev / ħc_MeV_fm
    muq_fm = muq_mev / ħc_MeV_fm

    @printf("Point: T=%.3f MeV (%.6f fm^-1), muB=%.3f MeV, muq=%.3f MeV (%.6f fm^-1), xi=%.3f\n", opt.T_mev, T_fm, opt.muB_mev, muq_mev, muq_fm, opt.xi)

    base = solve(FixedMu(), T_fm, muq_fm; xi=opt.xi, p_num=opt.p_num, t_num=opt.t_num)
    @printf("PNJL: converged=%s, iterations=%d, residual=%.3e\n", string(Bool(base.converged)), Int(base.iterations), Float64(base.residual_norm))

    Φ = Float64(base.x_state[4])
    Φbar = Float64(base.x_state[5])
    masses = (u=Float64(base.masses[1]), d=Float64(base.masses[2]), s=Float64(base.masses[3]))
    @printf("Phi=%.6f, Phibar=%.6f\n", Φ, Φbar)
    @printf("masses (fm^-1): mu=%s, md=%s, ms=%s\n", fmt(masses.u), fmt(masses.d), fmt(masses.s))

    thermo_params = (T=Float64(T_fm), Φ=Φ, Φbar=Φbar, ξ=Float64(opt.xi))
    quark_params_basic = (m=masses, μ=(u=Float64(muq_fm), d=Float64(muq_fm), s=Float64(muq_fm)))

    ktmp = build_K_coeffs(Float64(T_fm), Float64(muq_fm), masses, Φ, Φbar)
    quark_params = (m=quark_params_basic.m, μ=quark_params_basic.μ, A=ktmp.A_vals)
    K_coeffs = ktmp.K_coeffs

    thermal_nodes = cached_nodes(opt.p_num, opt.t_num)
    densities = densities_from_equilibrium(base.x_state, base.mu_vec, T_fm, thermal_nodes, Float64(opt.xi))

    println("densities (fm^-3):")
    @printf("  n_u=%s  n_d=%s  n_s=%s\n", fmt(densities.u), fmt(densities.d), fmt(densities.s))
    @printf("  n_ubar=%s  n_dbar=%s  n_sbar=%s\n", fmt(densities.ubar), fmt(densities.dbar), fmt(densities.sbar))

    println("building sigma caches...")
    cs_caches = build_sigma_caches(quark_params, thermo_params, K_coeffs; n_sigma_points=opt.tau_n_sigma_points, sigma_grid_n=opt.sigma_grid_n)

    # `Λ_inv_fm` is Λ (fm⁻¹) despite its historical name.
    Λ = Λ_inv_fm
    println("sigma cache summary (per process):")
    for process in REQUIRED_PROCESSES
        cache = cs_caches[process]
        smin = cache.s_vals[1]
        smax = cache.s_vals[end]
        σmin = minimum(cache.sigma_vals)
        σmax = maximum(cache.sigma_vals)

        # rough upper bound of s reachable inside numerator when using p_cut = Λ and cosΘ = -1
        pi_sym, pj_sym, _, _ = RT_ASR.parse_particles_from_process(process)
        mi = RT_ASR.get_mass(pi_sym, quark_params)
        mj = RT_ASR.get_mass(pj_sym, quark_params)
        Ei = sqrt(Λ^2 + mi^2)
        Ej = sqrt(Λ^2 + mj^2)
        s_up_est = mi^2 + mj^2 + 2.0 * (Ei * Ej + Λ^2)

        @printf("  %-20s s:[%s, %s]  σ:[%s, %s]  s_up_est(Λ)=%s\n",
            string(process), fmt(smin), fmt(smax), fmt(σmin), fmt(σmax), fmt(s_up_est))
    end

    println("computing relaxation times...")
    res = relaxation_times(
        quark_params,
        thermo_params,
        K_coeffs;
        densities=densities,
        cs_caches=cs_caches,
        p_nodes=opt.tau_p_nodes,
        angle_nodes=opt.tau_angle_nodes,
        phi_nodes=opt.tau_phi_nodes,
        n_sigma_points=opt.tau_n_sigma_points,
    )

    tau = res.tau
    tauinv = res.tau_inv

    println("tau_inv (1/fm):")
    @printf("  u=%s  s=%s  ubar=%s  sbar=%s\n", fmt(tauinv.u), fmt(tauinv.s), fmt(tauinv.ubar), fmt(tauinv.sbar))
    println("tau (fm):")
    @printf("  u=%s  s=%s  ubar=%s  sbar=%s\n", fmt(tau.u), fmt(tau.s), fmt(tau.ubar), fmt(tau.sbar))

    println("done")
end

main()
