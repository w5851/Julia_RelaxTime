# PNJL solver robustness smoke test (deterministic random sampling)
#
# This is intentionally *deterministic* (fixed RNG seed) so it can live under
# unit tests without flakiness.
#
# It validates:
# - solve(FixedMu()) returns a converged solution
# - masses are positive and finite
# - Phi/Phibar are within a tolerant physical range

using Test
using Random

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

using .Constants_PNJL: ħc_MeV_fm
using .PNJL: solve
using .PNJL.ConstraintModes: FixedMu

@inline function _phys_ok(res; phi_tol::Float64=1e-8)
    Φ = res.x_state[4]
    Φbar = res.x_state[5]
    phi_ok = isfinite(Φ) && isfinite(Φbar) && (-phi_tol <= Φ <= 1 + phi_tol) && (-phi_tol <= Φbar <= 1 + phi_tol)
    mass_ok = all(isfinite, res.masses) && all(>(0.0), res.masses)
    return phi_ok && mass_ok
end

@testset "PNJL solve(FixedMu) random physical smoke" begin
    # Two fixed regression points found by the random sampler/diagnostics.
    @testset "Regression points" begin
        # 1) Previously: residual small but negative masses from some seeds.
        T_MeV = 48.269077
        muB_MeV = 972.640900
        xi = 0.479315
        p_num = 12
        t_num = 4
        res = solve(FixedMu(), T_MeV / ħc_MeV_fm, (muB_MeV / 3.0) / ħc_MeV_fm; xi=xi, p_num=p_num, t_num=t_num)
        @test res.converged
        @test _phys_ok(res)

        # 2) Previously: default seeds could converge to Phi/Phibar out-of-range.
        T_MeV = 288.453662
        muB_MeV = 563.275131
        xi = -0.652469
        res = solve(FixedMu(), T_MeV / ħc_MeV_fm, (muB_MeV / 3.0) / ħc_MeV_fm; xi=xi, p_num=p_num, t_num=t_num)
        @test res.converged
        @test _phys_ok(res)
    end

    # Deterministic random stress: keep modest by default.
    @testset "Deterministic random sampling" begin
        # Keep this small by default so it remains unit-test friendly.
        # For heavier stress runs, override via ENV (e.g. PNJL_RANDOM_SMOKE_N=2000).
        n_samples = parse(Int, get(ENV, "PNJL_RANDOM_SMOKE_N", "50"))
        seed = parse(Int, get(ENV, "PNJL_RANDOM_SMOKE_SEED", "2"))
        p_num = parse(Int, get(ENV, "PNJL_RANDOM_SMOKE_P_NUM", "12"))
        t_num = parse(Int, get(ENV, "PNJL_RANDOM_SMOKE_T_NUM", "4"))

        rng = MersenneTwister(seed)

        # Allow disabling the random loop (while keeping regression points) with N=0.
        if n_samples <= 0
            @test true
            return
        end

        # Sampling ranges (MeV)
        T_min, T_max = 0.0, 350.0
        muB_min, muB_max = 0.0, 1800.0
        xi_min, xi_max = -0.8, 0.8

        for _ in 1:n_samples
            T_MeV = T_min + rand(rng) * (T_max - T_min)
            muB_MeV = muB_min + rand(rng) * (muB_max - muB_min)
            xi = xi_min + rand(rng) * (xi_max - xi_min)

            # Avoid exact T=0.
            T_fm = max(T_MeV, 1e-6) / ħc_MeV_fm
            muq_fm = (muB_MeV / 3.0) / ħc_MeV_fm

            res = solve(FixedMu(), T_fm, muq_fm; xi=xi, p_num=p_num, t_num=t_num)
            @test res.converged
            @test _phys_ok(res)
            @test isfinite(res.omega) && isfinite(res.pressure) && isfinite(res.rho_norm)
        end
    end
end
