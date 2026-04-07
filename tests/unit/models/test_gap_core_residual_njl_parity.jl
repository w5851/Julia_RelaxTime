using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "gap_core_residual njl parity" begin
    p_num = 8
    t_num = 4
    xi = 0.0
    thermal_nodes = Models.cached_nodes(p_num, t_num)

    @testset "NJL2 dim=2" begin
        model = Models.create_model(:NJL2)
        x_state = SVector{2}(-1.84, -1.82)
        mu_vec = SVector{3}(0.0, 0.0, 0.0)
        params = Models.GapParams(0.5, thermal_nodes, xi; p_num=p_num, t_num=t_num, model_kind=:NJL2)

        F_core = zeros(2)
        Models.gap_core_residual!(F_core, x_state, mu_vec, params)
        F_ref = Models.gap_residual(model, x_state, params.T_fm, mu_vec; p_num=p_num, t_num=t_num, xi=xi)

        @test length(F_core) == 2
        for i in 1:2
            @test isapprox(F_core[i], F_ref[i]; rtol=1e-9, atol=1e-9)
        end
    end

    @testset "NJL dim=3" begin
        model = Models.create_model(:NJL)
        x_state = SVector{3}(-1.84, -1.84, -2.23)
        mu_vec = SVector{3}(0.0, 0.0, 0.0)
        params = Models.GapParams(0.5, thermal_nodes, xi; p_num=p_num, t_num=t_num, model_kind=:NJL)

        F_core = zeros(3)
        Models.gap_core_residual!(F_core, x_state, mu_vec, params)
        F_ref = Models.gap_residual(model, x_state, params.T_fm, mu_vec; p_num=p_num, t_num=t_num, xi=xi)

        @test length(F_core) == 3
        for i in 1:3
            @test isapprox(F_core[i], F_ref[i]; rtol=1e-9, atol=1e-9)
        end
    end
end
