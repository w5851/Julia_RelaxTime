"""
Status: deprecate-candidate
Purpose: PNJL single-point smoke self-check script.
Replacement: prefer unit smoke profile and `scripts/pnjl/run_tmu_scan.jl` for maintained scan workflows.
"""

# Smoke tests for PNJL solvers (using new architecture)
using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
Models.pnjl_module()

const DEFAULT_SCAN_NUMERIC = Models.default_scan_numeric_options()

function run_case_mu(T_mev, mu_mev; xi=0.0, p_num=DEFAULT_SCAN_NUMERIC.p_num, t_num=DEFAULT_SCAN_NUMERIC.t_num)
    @printf("\nCase solve(FixedMu()): T=%.2f MeV, mu=%.2f MeV, xi=%.3f\n", T_mev, mu_mev, xi)
    try
        res = Models.solve_pnjl_point(T_mev=T_mev, mu_mev=mu_mev; xi=xi, p_num=p_num, t_num=t_num)
        @printf("  converged=%s iterations=%d residual=%.3e pressure=%.6f rho=%.6f energy=%.6f\n",
            string(res.converged), res.iterations, res.residual_norm, res.pressure, res.rho_norm, res.energy)
    catch e
        println("  ERROR: ", e)
        @show stacktrace(catch_backtrace())
    end
end

function run_case_rho(T_mev, rho_target; xi=0.0, p_num=DEFAULT_SCAN_NUMERIC.p_num, t_num=DEFAULT_SCAN_NUMERIC.t_num)
    @printf("\nCase solve(FixedRho()): T=%.2f MeV, rho=%.4f, xi=%.3f\n", T_mev, rho_target, xi)
    try
        res = Models.solve_pnjl_point(T_mev=T_mev, rho_target=rho_target; xi=xi, p_num=p_num, t_num=t_num)
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
