# PhaseSpaceSampling 模块单元测试

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
using Test

const _GL_PATH_PSS = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "integration", "GaussLegendre.jl"))
if !isdefined(Main, :GaussLegendre)
    Base.include(Main, _GL_PATH_PSS)
end
const _PSS_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "integration", "PhaseSpaceSampling.jl"))
if !isdefined(Main, :PhaseSpaceSampling)
    Base.include(Main, _PSS_PATH)
end

using Main.PhaseSpaceSampling: p_nodes_weights, cos_nodes_weights, integrate_p, integrate_p_cos

@testset "PhaseSpaceSampling" begin
    @testset "p_nodes_weights 返回正确维度" begin
        nodes, weights = p_nodes_weights(16, 5.0, nothing, nothing)
        @test length(nodes) == 16
        @test length(weights) == 16
        @test all(n -> n ≥ 0, nodes)
        @test all(w -> w > 0, weights)
    end

    @testset "cos_nodes_weights 返回正确范围" begin
        nodes, weights = cos_nodes_weights(8, nothing, nothing)
        @test length(nodes) == 8
        @test all(n -> -1.0 ≤ n ≤ 1.0, nodes)
    end

    @testset "integrate_p 简单积分" begin
        nodes, weights = p_nodes_weights(32, 5.0, nothing, nothing)
        # 积分 p^2 从 0 到 p_max
        result = integrate_p(p -> p^2, nodes, weights)
        @test isfinite(result)
        @test result > 0
    end

    @testset "integrate_p_cos 二维积分" begin
        p_n, p_w = p_nodes_weights(16, 3.0, nothing, nothing)
        c_n, c_w = cos_nodes_weights(8, nothing, nothing)
        result = integrate_p_cos((p, c) -> p^2, p_n, p_w, c_n, c_w)
        @test isfinite(result)
        @test result > 0
    end
end
