# Smoke tests for PNJL solvers (using new architecture)
using Printf
using NLsolve: LineSearches

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

using .Constants_PNJL: ħc_MeV_fm
using .PNJL: solve, FixedMu, FixedRho

function run_case_mu(T_mev, mu_mev; xi=0.0, p_num=32, t_num=16)
    @printf("\nCase solve(FixedMu()): T=%.2f MeV, mu=%.2f MeV, xi=%.3f\n", T_mev, mu_mev, xi)
    try
        T_fm = T_mev / ħc_MeV_fm
        mu_fm = mu_mev / ħc_MeV_fm
        res = solve(FixedMu(), T_fm, mu_fm; xi=xi, p_num=p_num, t_num=t_num)
        @printf("  converged=%s iterations=%d residual=%.3e pressure=%.6f rho=%.6f energy=%.6f\n",
            string(res.converged), res.iterations, res.residual_norm, res.pressure, res.rho_norm, res.energy)
    catch e
        println("  ERROR: ", e)
        @show stacktrace(catch_backtrace())
    end
end

function run_case_rho(T_mev, rho_target; xi=0.0, p_num=32, t_num=16)
    @printf("\nCase solve(FixedRho()): T=%.2f MeV, rho=%.4f, xi=%.3f\n", T_mev, rho_target, xi)
    try
        T_fm = T_mev / ħc_MeV_fm
        res = solve(FixedRho(rho_target), T_fm; xi=xi, p_num=p_num, t_num=t_num, linesearch=LineSearches.BackTracking())
        @printf("  converged=%s iterations=%d residual=%.3e pressure=%.6f rho=%.6f energy=%.6f\n",
            string(res.converged), res.iterations, res.residual_norm, res.pressure, res.rho_norm, res.energy)
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
