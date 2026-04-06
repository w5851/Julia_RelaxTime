using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "solver named/vec parity" begin
    model = Models.create_model(:PNJL)

    @testset "FixedMu parity" begin
        mode = Models.FixedMu()
        theta_vec = [100.0 / 197.327, 0.0]
        theta_named = (T_fm=theta_vec[1], μ_fm=theta_vec[2])

        r_vec = Models.solve_vec(model, mode, theta_vec;
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )
        r_named = Models.solve_named(model, mode, theta_named;
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )

        @test r_vec isa Models.SolverResult
        @test r_named isa Models.SolverResult
        @test r_vec.converged == r_named.converged
        @test isapprox(r_vec.pressure, r_named.pressure; rtol=1e-3, atol=1e-5)
        @test isapprox(r_vec.residual_norm, r_named.residual_norm; rtol=1e-3, atol=1e-7)
    end

    @testset "FixedEntropy parity" begin
        mode = Models.FixedEntropy(0.5)
        theta_vec = [100.0 / 197.327]
        theta_named = (T_fm=theta_vec[1],)

        seed = copy(Models.pnjl_module().HADRON_SEED_8)

        r_vec = Models.solve_vec(model, mode, theta_vec;
            seed_guess=seed,
            seed_candidates=(seed,),
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )
        r_named = Models.solve_named(model, mode, theta_named;
            seed_guess=seed,
            seed_candidates=(seed,),
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )

        @test r_vec isa Models.SolverResult
        @test r_named isa Models.SolverResult
        @test r_vec.converged == r_named.converged
        @test isapprox(r_vec.pressure, r_named.pressure; rtol=1e-3, atol=1e-5)
        @test isapprox(r_vec.residual_norm, r_named.residual_norm; rtol=1e-3, atol=1e-7)
    end
end
