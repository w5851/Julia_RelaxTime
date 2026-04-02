using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const P = Models.pnjl_module()

@testset "dimension-agnostic scan result adapters smoke" begin
    @testset "TmuScan adapter keeps dynamic state/mu dims" begin
        raw = (
            converged = true,
            solution = [1.0, 2.0, 3.0, 4.0],
            x_state = [0.1, 0.2, 0.3, 0.4],
            mu_vec = [0.01, 0.02],
            omega = -1.0,
            pressure = 1.0,
            rho_norm = 0.0,
            entropy = 0.0,
            energy = 0.0,
            masses = [0.31, 0.32, 0.33],
            iterations = 1,
            residual_norm = 1e-8,
        )

        adapted = P.TmuScan._to_solver_result(P.FixedMu(), raw, 0.0)
        @test length(adapted.x_state) == 4
        @test length(adapted.mu_vec) == 2
    end

    @testset "TrhoScan adapter keeps dynamic state/mu dims" begin
        raw = (
            converged = true,
            solution = [1.0, 2.0, 3.0, 4.0, 5.0],
            x_state = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
            mu_vec = [0.01, 0.02, 0.03, 0.04],
            omega = -1.0,
            pressure = 1.0,
            rho_norm = 0.0,
            entropy = 0.0,
            energy = 0.0,
            masses = [0.31, 0.32, 0.33],
            iterations = 1,
            residual_norm = 1e-8,
        )

        adapted = P.TrhoScan._to_solver_result(P.FixedRho(1.0), raw, 0.0)
        @test length(adapted.x_state) == 6
        @test length(adapted.mu_vec) == 4
    end
end
