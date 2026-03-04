# EllipsoidCalculation 模块单元测试

using Test
using LinearAlgebra

const _ELLIPSOID_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "simulation", "EllipsoidCalculation.jl"))
if !isdefined(Main, :EllipsoidCalculation)
    Base.include(Main, _ELLIPSOID_PATH)
end

using Main.EllipsoidCalculation: EllipsoidParams, calculate_ellipsoid_parameters, sample_ellipsoid_surface

@testset "EllipsoidCalculation" begin
    @testset "EllipsoidParams 构建" begin
        # 位置参数构造: EllipsoidParams(center, axes_directions, half_lengths)
        params = EllipsoidParams(
            [0.0, 0.0, 0.0],
            Matrix(1.0I, 3, 3),
            [1.0, 2.0, 3.0]
        )
        @test params.center == [0.0, 0.0, 0.0]
        @test length(params.half_lengths) == 3
    end

    @testset "calculate_ellipsoid_parameters 正定矩阵" begin
        A = [2.0 0.0 0.0; 0.0 3.0 0.0; 0.0 0.0 4.0]
        p_star = 1.0
        params = calculate_ellipsoid_parameters(A, p_star)
        @test params isa EllipsoidParams
        @test all(h -> h > 0, params.half_lengths)
    end

    @testset "sample_ellipsoid_surface 生成采样点" begin
        A = Matrix(1.0I, 3, 3)
        b = [0.0, 0.0, 0.0]
        points = sample_ellipsoid_surface(10, 1.0, A, b)
        @test length(points) == 10
        for pt in points
            @test length(pt) == 3
            @test all(isfinite, pt)
        end
    end
end
