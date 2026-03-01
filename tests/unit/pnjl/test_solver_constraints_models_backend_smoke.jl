# PNJL constraints under thermo_backend=:models (smoke)
#
# Goal:
# - Ensure FixedRho / FixedEntropy / FixedSigma run end-to-end with the models thermo backend.
# - Catch accidental legacy fallbacks inside Conditions/constraints.

using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Constants_PNJL)
    include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
end
if !isdefined(Main, :GaussLegendre)
    include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
end
if !isdefined(Main, :PNJL)
    if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
    Models.pnjl_module()
end

P = PNJL

const ħc = 197.327  # MeV·fm

@testset "constraints with models thermo backend (smoke)" begin
    T_fm = 100.0 / ħc
    μ_fm = 20.0 / ħc

    # keep these small-ish for smoke runtime
    p_num = 16
    t_num = 6

    @testset "FixedMu" begin
        result = P.solve(P.FixedMu(), T_fm, μ_fm; thermo_backend=:models, p_num=p_num, t_num=t_num)
        @test result.converged
        @test isfinite(result.pressure)
        @test isfinite(result.rho_norm)
        @test all(isfinite.(result.masses))
        @test result.residual_norm <= 1e-2
    end

    @testset "FixedMu MultiSeed (models backend)" begin
        result = P.solve(P.FixedMu(), T_fm, μ_fm;
            thermo_backend=:models,
            seed_strategy=P.MultiSeed(),
            p_num=p_num,
            t_num=t_num,
        )
        @test result.converged
        @test isfinite(result.pressure)
        @test isfinite(result.rho_norm)
        @test all(isfinite.(result.masses))
        @test result.residual_norm <= 1e-2
    end

    @testset "FixedMu PhaseAware bootstrap MultiSeed (models backend)" begin
        tracker = P.PhaseAwareContinuitySeed(0.0; bootstrap_multiseed=true)
        result = P.solve(P.FixedMu(), T_fm, μ_fm;
            thermo_backend=:models,
            seed_strategy=tracker,
            p_num=p_num,
            t_num=t_num,
        )
        @test result.converged
        @test isfinite(result.pressure)
        @test isfinite(result.rho_norm)
        @test all(isfinite.(result.masses))
        @test result.residual_norm <= 1e-2
    end

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

    @testset "Implicit solver (models backend)" begin
        solver = P.create_implicit_solver(thermo_backend=:models, p_num=p_num, t_num=t_num)
        θ = [T_fm, μ_fm]
        x, _ = solver(θ)
        @test length(x) == 5
        @test all(isfinite.(x))

        d = P.solve_with_derivatives(T_fm, μ_fm;
            order=1,
            thermo_backend=:models,
            model_kind=:PNJL,
            p_num=p_num,
            t_num=t_num,
        )
        @test length(d.dx_dT) == 5
        @test length(d.dx_dμ) == 5
        @test all(isfinite.(d.dx_dT))
        @test all(isfinite.(d.dx_dμ))
    end
end
