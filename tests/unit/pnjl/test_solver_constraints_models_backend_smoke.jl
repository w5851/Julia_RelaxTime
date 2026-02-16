# PNJL constraints under thermo_backend=:models (smoke)
#
# Goal:
# - Ensure FixedRho / FixedEntropy / FixedSigma run end-to-end with the models thermo backend.
# - Catch accidental legacy fallbacks inside Conditions/constraints.

using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

P = PNJL

const ħc = 197.327  # MeV·fm

@testset "constraints with models thermo backend (smoke)" begin
    T_fm = 100.0 / ħc

    # keep these small-ish for smoke runtime
    p_num = 16
    t_num = 6

    @testset "FixedRho" begin
        target_rho = 1.0
        result = P.solve(P.FixedRho(target_rho), T_fm; thermo_backend=:models, p_num=p_num, t_num=t_num)
        @test result.converged
        @test isfinite(result.pressure)
        @test isfinite(result.rho_norm)
        @test all(isfinite.(result.masses))
        @test isapprox(result.rho_norm, target_rho; rtol=0.05)
    end

    @testset "FixedEntropy" begin
        target_s = 0.5
        result = P.solve(P.FixedEntropy(target_s), T_fm; thermo_backend=:models, p_num=p_num, t_num=t_num)
        @test result.converged
        @test isfinite(result.entropy)
        @test all(isfinite.(result.masses))
        @test isapprox(result.entropy, target_s; rtol=0.05)
    end

    @testset "FixedSigma" begin
        target_sigma = 10.0
        result = P.solve(P.FixedSigma(target_sigma), T_fm; thermo_backend=:models, p_num=p_num, t_num=t_num)
        @test result.converged
        @test isfinite(result.entropy)
        @test isfinite(result.rho_norm)
        @test all(isfinite.(result.masses))

        n_B = result.rho_norm * P.ρ0
        sigma = n_B > 1e-12 ? result.entropy / n_B : 0.0
        @test isapprox(sigma, target_sigma; rtol=0.05)
    end

    @testset "FixedAsymmetricRho (PNJL)" begin
        mode = P.FixedAsymmetricRho(0.05, 1.0, 0.0)
        result = P.solve(mode, T_fm;
            thermo_backend=:models,
            model_kind=:PNJL,
            p_num=p_num,
            t_num=t_num,
            iterations=120,
        )

        @test result.converged
        @test isfinite(result.pressure)
        @test isfinite(result.rho_norm)
        @test all(isfinite.(result.masses))
    end

    @testset "FixedAsymmetricRho (RPNJL)" begin
        mode = P.FixedAsymmetricRho(0.05, 1.0, 0.0)
        result = P.solve(mode, T_fm;
            thermo_backend=:models,
            model_kind=:RPNJL,
            p_num=p_num,
            t_num=t_num,
            iterations=120,
        )

        @test result.converged
        @test isfinite(result.pressure)
        @test isfinite(result.rho_norm)
        @test all(isfinite.(result.masses))
    end
end
