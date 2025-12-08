# Sensitivity scan for momentum/angle quadrature nodes in solve_fixed_mu
using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))

module PNJLLocal
    using ..Constants_PNJL
    include(joinpath(Main.PROJECT_ROOT, "src", "pnjl", "solvers", "AnisoGapSolver.jl"))
    export AnisoGapSolver
end

using .PNJLLocal: AnisoGapSolver
using Main.GaussLegendre: DEFAULT_momentum_POINTS, DEFAULT_theta_POINTS

const PRESSURE_TOL = 5e-7
const RHO_TOL = 5e-8
const ENERGY_TOL = 5e-7

const P_SWEEP = [16, 24, 32, 48, 64]
const T_SWEEP = [4, 8, 12, 16]
const REF_P_NUM = 96
const REF_T_NUM = 32

const CASES = [
    (name = "mid_xi", T_mev = 120.0, mu_mev = 30.0, xi = 0.3),
    (name = "high_xi", T_mev = 150.0, mu_mev = 60.0, xi = 1.0),
]

function solve_case(case_cfg; p_num = DEFAULT_momentum_POINTS, t_num = DEFAULT_theta_POINTS)
    mu_fm = case_cfg.mu_mev / Constants_PNJL.ħc_MeV_fm
    return AnisoGapSolver.solve_fixed_mu(case_cfg.T_mev, mu_fm; xi = case_cfg.xi, p_num = p_num, t_num = t_num)
end

within_tol(dp, dr, de) = dp <= PRESSURE_TOL && dr <= RHO_TOL && de <= ENERGY_TOL

function scan_case(case_cfg)
    println("\n== Case $(case_cfg.name): T=$(case_cfg.T_mev) MeV, mu=$(case_cfg.mu_mev) MeV, xi=$(case_cfg.xi) ==")
    ref = solve_case(case_cfg; p_num = REF_P_NUM, t_num = REF_T_NUM)
    println("reference pressure=$(ref.pressure), rho=$(ref.rho), energy=$(ref.energy)")

    default_res = solve_case(case_cfg)
    dp_def = abs(default_res.pressure - ref.pressure)
    dr_def = abs(default_res.rho - ref.rho)
    de_def = abs(default_res.energy - ref.energy)
    println("default nodes (p=$(DEFAULT_momentum_POINTS), t=$(DEFAULT_theta_POINTS)): Δp=$(dp_def) Δrho=$(dr_def) Δe=$(de_def)")
    default_ok = within_tol(dp_def, dr_def, de_def)

    println("\n-- momentum sweep (t=$(DEFAULT_theta_POINTS)) --")
    @printf("%6s %10s %14s %14s %14s\n", "p_num", "iters", "|Δp|", "|Δrho|", "|Δenergy|")
    min_p = nothing
    for p_num in P_SWEEP
        res = solve_case(case_cfg; p_num = p_num, t_num = DEFAULT_theta_POINTS)
        dp = abs(res.pressure - ref.pressure)
        dr = abs(res.rho - ref.rho)
        de = abs(res.energy - ref.energy)
        @printf("%6d %10d %14.6e %14.6e %14.6e\n", p_num, res.iterations, dp, dr, de)
        within_tol(dp, dr, de) && isnothing(min_p) && (min_p = p_num)
    end

    println("\n-- theta sweep (p=$(DEFAULT_momentum_POINTS)) --")
    @printf("%6s %10s %14s %14s %14s\n", "t_num", "iters", "|Δp|", "|Δrho|", "|Δenergy|")
    min_t = nothing
    for t_num in T_SWEEP
        res = solve_case(case_cfg; p_num = DEFAULT_momentum_POINTS, t_num = t_num)
        dp = abs(res.pressure - ref.pressure)
        dr = abs(res.rho - ref.rho)
        de = abs(res.energy - ref.energy)
        @printf("%6d %10d %14.6e %14.6e %14.6e\n", t_num, res.iterations, dp, dr, de)
        within_tol(dp, dr, de) && isnothing(min_t) && (min_t = t_num)
    end

    if default_ok
        println("Default nodes are within tolerance for this case.")
    else
        println("Default nodes exceed tolerance for this case.")
        min_p !== nothing && println("  Suggested p_num >= $(min_p)")
        min_t !== nothing && println("  Suggested t_num >= $(min_t)")
    end
end

function main()
    for case_cfg in CASES
        scan_case(case_cfg)
    end
end

main()
