using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

using .Models: run_trho_scan

@testset "trho scan semantic modes smoke" begin
    tmp_dir = mktempdir()

    ground_output = joinpath(tmp_dir, "trho_semantic_ground.csv")
    manifold_output = joinpath(tmp_dir, "trho_semantic_manifold.csv")

    ground = run_trho_scan(
        T_values=[150.0],
        rho_values=[0.2],
        xi_values=[0.0],
        output_path=ground_output,
        overwrite=true,
        resume=false,
        reverse_rho=false,
        solver_backend=:models,
        semantic_mode=:ground_state,
        p_num=12,
        t_num=4,
        iterations=80,
    )

    manifold = run_trho_scan(
        T_values=[150.0],
        rho_values=[0.2],
        xi_values=[0.0],
        output_path=manifold_output,
        overwrite=true,
        resume=false,
        reverse_rho=false,
        solver_backend=:models,
        semantic_mode=:constrained_manifold,
        p_num=12,
        t_num=4,
        iterations=80,
    )

    @test isfile(ground_output)
    @test isfile(manifold_output)
    @test ground.total == 1
    @test manifold.total == 1
    @test ground.success == 1
    @test manifold.success == 1
end
