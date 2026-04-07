using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "solver/implicit residual parity smoke" begin
    model = Models.create_model(:PNJL)
    p_num = 12
    t_num = 4
    xi = 0.0

    @testset "fixed-mu parity" begin
        adapters = Models.build_pnjl_fixedmu_adapters(model; p_num=p_num, t_num=t_num, xi=xi)
        theta = [0.5, 0.12]
        x, meta = adapters.forward_solve(theta)
        @test meta === nothing

        thermal_nodes = Models.cached_nodes(p_num, t_num)
        params = Models.GapParams(theta[1], thermal_nodes, xi; p_num=p_num, t_num=t_num, model_kind=:PNJL)

        residual_solver! = Models.build_residual!(Models.FixedMu(), SVector{3}(theta[2], theta[2], theta[2]), params)
        F_solver = zeros(5)
        residual_solver!(F_solver, x)

        F_implicit = adapters.conditions(theta, x, meta)
        @test length(F_implicit) == 5
        @test all(isfinite, F_implicit)

        for i in 1:5
            @test isapprox(F_solver[i], F_implicit[i]; rtol=1e-8, atol=1e-9)
        end
    end

    @testset "flavor-mu parity" begin
        adapters = Models.build_pnjl_flavor_mu_adapters(model; p_num=p_num, t_num=t_num, xi=xi)
        theta = [0.5, 0.18, 0.12, 0.06]
        x, meta = adapters.forward_solve(theta)
        @test meta === nothing

        thermal_nodes = Models.cached_nodes(p_num, t_num)
        params = Models.GapParams(theta[1], thermal_nodes, xi; p_num=p_num, t_num=t_num, model_kind=:PNJL)
        F_solver = zeros(5)
        mu_vec = SVector{3}(theta[2], theta[3], theta[4])
        Models.gap_core_residual!(F_solver, SVector{5}(Tuple(x)), mu_vec, params)

        F_implicit = adapters.conditions(theta, x, meta)
        @test length(F_implicit) == 5
        @test all(isfinite, F_implicit)

        for i in 1:5
            @test isapprox(F_solver[i], F_implicit[i]; rtol=1e-8, atol=1e-9)
        end
    end
end
