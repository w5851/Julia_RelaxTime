#!/usr/bin/env julia
# 诊断单个 (T, muB, xi) 点的 PNJL 能隙方程求解。
#
# 目标：通过“多初值 + 多求解器方法(NLsolve)”找到收敛且物理合理的解，
# 并按热力学势/压强排序输出。
#
# 用法示例：
#   julia --project=. scripts/pnjl/diagnose_gap_point.jl --T=349 --muB=521 --xi=0.577 --p_num=24 --t_num=8
#
# 可选：
#   --methods=newton,trust_region,anderson
#   --verbose

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Printf

include(joinpath(@__DIR__, "..", "..", "src", "Constants_PNJL.jl"))
include(joinpath(@__DIR__, "..", "..", "src", "pnjl", "PNJL.jl"))

using .PNJL
using StaticArrays
using NLsolve

const HBARC = 197.327 # MeV·fm

struct Config
    T_MeV::Float64
    muB_MeV::Float64
    xi::Float64
    p_num::Int
    t_num::Int
    methods::Vector{Symbol}
    verbose::Bool
end

function parse_args(args)
    T = 349.0
    muB = 521.0
    xi = 0.0
    p_num = 24
    t_num = 8
    methods = Symbol[:newton, :trust_region]
    verbose = false

    for a in args
        if startswith(a, "--T=")
            T = parse(Float64, a[5:end])
        elseif startswith(a, "--muB=")
            muB = parse(Float64, a[7:end])
        elseif startswith(a, "--xi=")
            xi = parse(Float64, a[6:end])
        elseif startswith(a, "--p_num=")
            p_num = parse(Int, a[9:end])
        elseif startswith(a, "--t_num=")
            t_num = parse(Int, a[9:end])
        elseif startswith(a, "--methods=")
            ms = split(a[11:end], ',')
            methods = Symbol.(ms)
        elseif a == "--verbose"
            verbose = true
        elseif a == "--help"
            println("Usage: julia --project=. scripts/pnjl/diagnose_gap_point.jl --T=<MeV> --muB=<MeV> --xi=<...> [--p_num=24 --t_num=8] [--methods=newton,trust_region] [--verbose]")
            exit(0)
        end
    end

    return Config(T, muB, xi, p_num, t_num, methods, verbose)
end

"""物理性判据（与求解器默认口径一致）：Polyakov loop 合理、有效质量为正。"""
function is_physical(x_state::SVector{5,Float64}, masses::SVector{3,Float64})
    φu, φd, φs, Φ, Φbar = x_state
    if !(0.0 <= Φ <= 1.0 && 0.0 <= Φbar <= 1.0)
        return false
    end
    if any(!isfinite, masses) || any(m -> m <= 0.0, masses)
        return false
    end
    return true
end

function build_seed_bank()
    # 内置种子（来自 SeedStrategies.jl）
    bank = Vector{Pair{String,Vector{Float64}}}()
    push!(bank, "HADRON_SEED_5" => Float64.(PNJL.HADRON_SEED_5))
    push!(bank, "MEDIUM_SEED_5" => Float64.(PNJL.MEDIUM_SEED_5))
    push!(bank, "HIGH_DENSITY_SEED_5" => Float64.(PNJL.HIGH_DENSITY_SEED_5))
    push!(bank, "HIGH_TEMP_SEED_5" => Float64.(PNJL.HIGH_TEMP_SEED_5))

    # 人工高温候选：凝聚幅度更小、Polyakov loop 更高
    push!(bank, "HT_GUESS_0p8" => [-0.50, -0.50, -1.20, 0.80, 0.80])
    push!(bank, "HT_GUESS_0p9" => [-0.30, -0.30, -0.90, 0.90, 0.90])
    push!(bank, "HT_GUESS_0p95" => [-0.20, -0.20, -0.70, 0.95, 0.95])

    # 更偏禁闭：Φ 很小，但凝聚较弱（用于避免卡到坏分支）
    push!(bank, "WEAK_CHIRAL_CONF" => [-0.50, -0.50, -1.20, 1e-3, 1e-3])

    return bank
end

function solve_once(T_fm::Float64, mu_fm::Float64, xi::Float64, p_num::Int, t_num::Int, seed::Vector{Float64}, method::Symbol; verbose::Bool=false)
    thermal_nodes = PNJL.cached_nodes(p_num, t_num)
    params = PNJL.GapParams(T_fm, thermal_nodes, xi)
    mu_vec = SVector{3}(mu_fm, mu_fm, mu_fm)
    residual_fn! = PNJL.build_residual!(PNJL.FixedMu(), mu_vec, params)

    # 一些 NLsolve 选项：trust_region 往往需要更宽松的 initial step；anderson 用来兜底
    t0 = time()
    res = nlsolve(residual_fn!, seed;
        autodiff = :forward,
        method = method,
        xtol = 1e-10,
        ftol = 1e-10,
        iterations = 400,
    )
    elapsed = time() - t0

    x_state = SVector{5}(Tuple(res.zero))
    pressure, rho_norm, entropy, energy = PNJL.calculate_thermo(x_state, mu_vec, T_fm, thermal_nodes, xi)
    omega = -pressure
    masses = PNJL.calculate_mass_vec(x_state)

    if verbose
        @printf("    converged=%s iter=%d |res|=%.3e  omega=%.6e  time=%.3fs\n", string(res.f_converged), res.iterations, res.residual_norm, omega, elapsed)
        @printf("    x=[% .4f % .4f % .4f % .4f % .4f] masses=[%.4f %.4f %.4f]\n", x_state[1], x_state[2], x_state[3], x_state[4], x_state[5], masses[1], masses[2], masses[3])
    end

    return (converged = res.f_converged,
            iterations = res.iterations,
            residual_norm = res.residual_norm,
            elapsed = elapsed,
            x_state = x_state,
            masses = masses,
            omega = omega,
            pressure = pressure,
            rho_norm = rho_norm,
            entropy = entropy,
            energy = energy)
end

function main()
    cfg = parse_args(ARGS)

    T_fm = cfg.T_MeV / HBARC
    mu_fm = (cfg.muB_MeV / 3.0) / HBARC

    println("diagnose_gap_point")
    @printf("  T=%.3f MeV  muB=%.3f MeV (muq=%.3f)  xi=%.3f\n", cfg.T_MeV, cfg.muB_MeV, cfg.muB_MeV/3.0, cfg.xi)
    @printf("  T_fm=%.6f  mu_fm=%.6f  grid(p_num=%d,t_num=%d)\n", T_fm, mu_fm, cfg.p_num, cfg.t_num)
    println("  methods=", join(string.(cfg.methods), ","))

    seeds = build_seed_bank()
    results = []

    for m in cfg.methods
        for (name, seed) in seeds
            println("- try method=", m, " seed=", name)
            local out
            try
                out = solve_once(T_fm, mu_fm, cfg.xi, cfg.p_num, cfg.t_num, seed, m; verbose=cfg.verbose)
            catch err
                println("    ERROR: ", err)
                continue
            end
            phys = is_physical(out.x_state, out.masses)
            push!(results, (method=m, seed=name, phys=phys, out...))
            @printf("    ok=%s phys=%s iter=%d |res|=%.3e  Phi=%.4f Phibar=%.4f  masses=[%.4f %.4f %.4f]  omega=%.6e\n",
                string(out.converged), string(phys), out.iterations, out.residual_norm,
                out.x_state[4], out.x_state[5], out.masses[1], out.masses[2], out.masses[3], out.omega)
        end
    end

    good = filter(r -> r.converged && r.phys, results)
    if isempty(good)
        println("\nNo converged physical solution found with current seed/method bank.")
        best = sort(results; by=r->(r.converged ? 0 : 1, r.residual_norm))
        println("Top 5 (by converged then residual_norm):")
        for r in Iterators.take(best, min(5, length(best)))
            @printf("  method=%s seed=%s ok=%s phys=%s |res|=%.3e Phi=%.4f Phibar=%.4f masses=[%.4f %.4f %.4f]\n",
                string(r.method), r.seed, string(r.converged), string(r.phys), r.residual_norm,
                r.x_state[4], r.x_state[5], r.masses[1], r.masses[2], r.masses[3])
        end
        exit(2)
    end

    best_phys = sort(good; by=r->r.omega)  # omega 越小越好（P 越大）
    println("\nConverged physical solutions (sorted by omega):")
    for r in best_phys
        @printf("  omega=%.6e  method=%s seed=%s |res|=%.3e Phi=%.4f Phibar=%.4f masses=[%.4f %.4f %.4f] x=[% .4f % .4f % .4f]\n",
            r.omega, string(r.method), r.seed, r.residual_norm, r.x_state[4], r.x_state[5],
            r.masses[1], r.masses[2], r.masses[3], r.x_state[1], r.x_state[2], r.x_state[3])
    end

    top = first(best_phys)
    println("\nSelected best physical solution:")
    @printf("  method=%s seed=%s omega=%.6e |res|=%.3e\n", string(top.method), top.seed, top.omega, top.residual_norm)
    @printf("  x_state=[% .6f, % .6f, % .6f, %.6f, %.6f] (phi_u,phi_d,phi_s,Phi,Phibar)\n",
        top.x_state[1], top.x_state[2], top.x_state[3], top.x_state[4], top.x_state[5])
    @printf("  masses=[%.6f, %.6f, %.6f] (fm^-1)\n", top.masses[1], top.masses[2], top.masses[3])
    @printf("  masses(MeV)=[%.3f, %.3f, %.3f]\n", top.masses[1]*HBARC, top.masses[2]*HBARC, top.masses[3]*HBARC)
end

main()
