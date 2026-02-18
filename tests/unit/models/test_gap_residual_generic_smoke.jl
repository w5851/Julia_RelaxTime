using Test

# New models entry
_models_entry = joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    include(_models_entry)
end
using .Models

using StaticArrays

@testset "Models gap_residual generic smoke" begin
    T = 0.8
    mu_vec = @SVector [0.1, 0.1, 0.1]
    xi = 0.0

    p_num = 32
    t_num = 8

    @testset "NJL (3D)" begin
        m = Models.create_model(:NJL; profile="default")

        x3 = @SVector [-1.84329, -1.84329, -2.22701]
        r3 = Models.gap_residual(m, x3, T, mu_vec; p_num=p_num, t_num=t_num, xi=xi)

        @test length(r3) == 3
        @test all(isfinite, r3)

        st3 = Models.MeanFieldState(@SVector [x3[1], x3[2], x3[3]]; Phi=0.5, PhiBar=0.5)
        r3b = Models.gap_residual(m, st3, T, mu_vec; p_num=p_num, t_num=t_num, xi=xi)
        @test r3b == r3
    end

    @testset "PNJL (5D)" begin
        m = Models.create_model(:PNJL; profile="default", physics_profile="default")

        x5 = @SVector [0.02, 0.02, 0.03, 0.5, 0.5]
        r5 = Models.gap_residual(m, x5, T, mu_vec; p_num=p_num, t_num=t_num, xi=xi)

        @test length(r5) == 5
        @test all(isfinite, r5)

        st5 = Models.MeanFieldState(@SVector [x5[1], x5[2], x5[3]]; Phi=x5[4], PhiBar=x5[5])
        r5b = Models.gap_residual(m, st5, T, mu_vec; p_num=p_num, t_num=t_num, xi=xi)
        @test r5b == r5
    end
end
