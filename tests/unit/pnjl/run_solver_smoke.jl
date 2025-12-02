# Smoke tests for PNJL solvers (solve_fixed_rho / solve_fixed_mu)
using Printf
using NLsolve: LineSearches

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))

module PNJLLocal
    using ..Constants_PNJL
    include(joinpath(Main.PROJECT_ROOT, "src", "pnjl", "solvers", "AnisoGapSolver.jl"))
    export AnisoGapSolver
end

using .PNJLLocal: AnisoGapSolver

function run_case_mu(T_mev, mu_mev; xi=0.0, p_num=32, t_num=16)
    @printf("\nCase solve_fixed_mu: T=%.2f MeV, mu=%.2f MeV, xi=%.3f\n", T_mev, mu_mev, xi)
    try
        res = AnisoGapSolver.solve_fixed_mu(T_mev, mu_mev/Constants_PNJL.ħc_MeV_fm; xi=xi, p_num=p_num, t_num=t_num)
        @printf("  converged=%s iterations=%d residual=%.3e pressure=%.6f rho=%.6f energy=%.6f\n",
            string(res.converged), res.iterations, res.residual_norm, res.pressure, res.rho, res.energy)
    catch e
        println("  ERROR: ", e)
        @show stacktrace(catch_backtrace())
    end
end

function run_case_rho(T_mev, rho_target; xi=0.0, p_num=32, t_num=16)
    @printf("\nCase solve_fixed_rho: T=%.2f MeV, rho=%.4f, xi=%.3f\n", T_mev, rho_target, xi)
    try
        res = AnisoGapSolver.solve_fixed_rho(T_mev, rho_target; xi=xi, p_num=p_num, t_num=t_num, linesearch=LineSearches.BackTracking())
        @printf("  converged=%s iterations=%d residual=%.3e pressure=%.6f rho=%.6f energy=%.6f\n",
            string(res.converged), res.iterations, res.residual_norm, res.pressure, res.rho, res.energy)
    catch e
        println("  ERROR: ", e)
        @show stacktrace(catch_backtrace())
    end
end

# Run a few smoke cases
run_case_mu(50.0, 10.0; xi=0.0, p_num=24, t_num=12)
run_case_mu(100.0, 0.0; xi=0.5, p_num=24, t_num=12)
run_case_rho(50.0, 0.05; xi=0.0, p_num=24, t_num=12)
run_case_rho(100.0, 0.1; xi=0.5, p_num=24, t_num=12)
