# Sweep PNJL solver over a T-μ grid to validate stability (using new architecture)
using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
Models.pnjl_module()

const DEFAULT_SCAN_NUMERIC = Models.default_scan_numeric_options()

const T_CASES = [50.0, 100.0, 150.0]
const MU_CASES = [-50.0, 0.0, 50.0, 100.0]
const XI_CASES = [0.0, 0.5]

function run_case_tmu(T_mev, mu_mev; xi=0.0, p_num=DEFAULT_SCAN_NUMERIC.p_num, t_num=DEFAULT_SCAN_NUMERIC.t_num)
    @printf("Case solve(FixedMu()) (T-μ grid): T=%.2f MeV, μ=%.2f MeV, xi=%.3f\n", T_mev, mu_mev, xi)
    res = Models.solve_pnjl_point(T_mev=T_mev, mu_mev=mu_mev; xi=xi, p_num=p_num, t_num=t_num)
    @assert res.converged "solve(FixedMu()) did not converge for T=$(T_mev), μ=$(mu_mev), xi=$(xi)"
    @assert isfinite(res.pressure) && isfinite(res.energy) && isfinite(res.rho_norm) "Non-finite thermodynamic results"
    @printf("  pressure=%.6f rho=%.6f energy=%.6f residual=%.3e iterations=%d\n",
        res.pressure,
        res.rho_norm,
        res.energy,
        res.residual_norm,
        res.iterations,
    )
end

function run_grid()
    for xi in XI_CASES, T_mev in T_CASES, mu_mev in MU_CASES
        run_case_tmu(T_mev, mu_mev; xi=xi)
    end
end

run_grid()
