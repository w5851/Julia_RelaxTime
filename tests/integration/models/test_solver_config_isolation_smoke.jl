using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "solver config isolation smoke" begin
    @testset "pnjl solve_gap calls keep independent config" begin
        T_fm = 0.52
        μ_fm = 0.16
        μ_vec = SVector(μ_fm, μ_fm, μ_fm)
        model = Models.create_model(:PNJL)

        state_a = Models.solve_gap(model, T_fm, μ_vec; xi=0.0, p_num=12, t_num=4)
        state_b = Models.solve_gap(model, T_fm, μ_vec; xi=0.35, p_num=24, t_num=8)

        x_a = collect(Models.state_vector(state_a))
        x_b = collect(Models.state_vector(state_b))

        @test length(x_a) == 5
        @test length(x_b) == 5
        @test all(isfinite, x_a)
        @test all(isfinite, x_b)

        @test !all(isapprox.(x_a, x_b; rtol=1e-10, atol=1e-12))
    end
end
