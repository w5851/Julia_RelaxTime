# Quantify sensitivity of angle quadrature nodes on PNJL single-point observables under multiple regimes
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
using Main.GaussLegendre: DEFAULT_theta_POINTS

const THETA_SWEEP = [4, 8, 12, 16, 24, 32, 40, 48, 64]
const DEFAULT_T_NUM = DEFAULT_theta_POINTS

const PRESSURE_TOL = 1e-8
const RHO_TOL = 1e-9
const ENERGY_TOL = 1e-8

const CASES = [
    (name = "baseline", T_mev = 120.0, mu_mev = 30.0, xi = 0.3, p_num = 32),
    (name = "high_xi", T_mev = 120.0, mu_mev = 30.0, xi = 1.2, p_num = 48),
    (name = "dense_p", T_mev = 160.0, mu_mev = 60.0, xi = 0.8, p_num = 64),
]

function solve_with_theta(case_cfg, t_num)
    mu_fm = case_cfg.mu_mev / Constants_PNJL.ħc_MeV_fm
    return AnisoGapSolver.solve_fixed_mu(case_cfg.T_mev, mu_fm;
        xi = case_cfg.xi,
        p_num = case_cfg.p_num,
        t_num = t_num,
    )
end

within_tol(errors) = errors.dp <= PRESSURE_TOL && errors.dr <= RHO_TOL && errors.de <= ENERGY_TOL

function evaluate_default(errors_by_theta)
    default_errors = get(errors_by_theta, DEFAULT_T_NUM, nothing)
    default_errors === nothing && return (true, maximum(THETA_SWEEP))
    within_tol(default_errors) && return (false, DEFAULT_T_NUM)
    for t_num in THETA_SWEEP
        errors = errors_by_theta[t_num]
        if within_tol(errors)
            return (true, t_num)
        end
    end
    return (true, maximum(THETA_SWEEP))
end

function analyze_case(case_cfg)
    println("\n== Case $(case_cfg.name): T=$(case_cfg.T_mev) MeV, μ=$(case_cfg.mu_mev) MeV, ξ=$(case_cfg.xi), p_num=$(case_cfg.p_num) ==")
    reference = solve_with_theta(case_cfg, maximum(THETA_SWEEP))
    println("reference pressure=$(reference.pressure)")
    errors_by_theta = Dict{Int, NamedTuple{(:dp, :dr, :de), NTuple{3, Float64}}}()
    @printf("%8s %10s %14s %14s %14s\n", "t_num", "iters", "|Δp|", "|Δρ|", "|Δenergy|")
    for t_num in THETA_SWEEP
        res = solve_with_theta(case_cfg, t_num)
        dp = abs(res.pressure - reference.pressure)
        dr = abs(res.rho - reference.rho)
        de = abs(res.energy - reference.energy)
        errors_by_theta[t_num] = (dp = dp, dr = dr, de = de)
        @printf("%8d %10d %14.6e %14.6e %14.6e\n", t_num, res.iterations, dp, dr, de)
    end
    violates, suggestion = evaluate_default(errors_by_theta)
    if violates
        println("!! Default t_num=$(DEFAULT_T_NUM) exceeds tolerance for this case. Suggested t_num=$(suggestion) or higher.")
    else
        println("Default t_num=$(DEFAULT_T_NUM) satisfies tolerance for this case.")
    end
    return violates
end

function main()
    violations = false
    for case_cfg in CASES
        violations |= analyze_case(case_cfg)
    end
    if violations
        println("\nSummary: At least one case requires more than t_num=$(DEFAULT_T_NUM). Consider increasing θ nodes for demanding runs.")
    else
        println("\nSummary: All evaluated cases are within tolerance using t_num=$(DEFAULT_T_NUM).")
    end
end

main()
