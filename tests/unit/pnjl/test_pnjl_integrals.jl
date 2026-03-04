# PNJLIntegrals 模块单元测试
#
# 测试内容：
# 1. cached_nodes 节点缓存
# 2. calculate_log_sum 热项积分

using Test
using StaticArrays

const PROJECT_ROOT_PI = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT_PI, "src", "models", "Models.jl"))
end
Models.pnjl_module()

const PNJL_PI = Models.pnjl_module()
const cached_nodes_pi = Models.cached_nodes  # Models 导出的包装函数
const calculate_log_sum_pi = Models.PNJLIntegrals.calculate_log_sum

@testset "PNJLIntegrals" begin
    @testset "cached_nodes 返回正确维度" begin
        nodes = cached_nodes_pi(24, 8)
        @test nodes isa NTuple{3, Matrix{Float64}}
        p_mesh, cosθ_mesh, coeff_mesh = nodes
        @test size(p_mesh) == (24, 8)
        @test size(cosθ_mesh) == (24, 8)
        @test size(coeff_mesh) == (24, 8)
    end

    @testset "cached_nodes 缓存一致性" begin
        nodes1 = cached_nodes_pi(24, 8)
        nodes2 = cached_nodes_pi(24, 8)
        @test nodes1[1] === nodes2[1]  # 同一对象（缓存命中）
    end

    @testset "calculate_log_sum 返回有限值" begin
        masses = SVector{3, Float64}(0.3, 0.31, 0.5)
        nodes = cached_nodes_pi(24, 8)
        p_mesh, cosθ_mesh, coeff_mesh = nodes
        Φ, Φbar = 0.5, 0.5
        mu_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
        T_fm = 0.8
        xi = 0.0

        result = calculate_log_sum_pi(masses, p_mesh, cosθ_mesh, coeff_mesh, Φ, Φbar, mu_vec, T_fm, xi)
        @test isfinite(result)
    end
end
