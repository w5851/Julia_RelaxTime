# PNJLIntegrals 模块单元测试
#
# 测试内容：
# 1. cached_nodes 节点缓存
# 2. calculate_log_sum 热项积分

using Test
using ForwardDiff
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

    @testset "RS 角约化只作用于分布自变量和测度" begin
        masses = SVector(0.3, 0.31, 0.5)
        mu_vec = SVector(0.4, 0.35, 0.2)
        Φ, Φbar = 0.3, 0.25
        T_fm = 0.2
        xi = 0.2

        adaptive = Models.PNJLIntegrals.calculate_log_sum_rs_reduced_adaptive(
            masses, Φ, Φbar, mu_vec, T_fm, xi;
            rtol=1e-10, atol=1e-12,
        )
        isotropic = Models.PNJLIntegrals.calculate_log_sum_rs_reduced_adaptive(
            masses, Φ, Φbar, mu_vec, T_fm, 0.0;
            rtol=1e-10, atol=1e-12,
        )
        tensor_nodes = cached_nodes_pi(256, 24)
        tensor = calculate_log_sum_pi(
            masses, tensor_nodes..., Φ, Φbar, mu_vec, T_fm, xi,
        )

        @test adaptive ≈ isotropic / sqrt(1 + xi) rtol=5e-10
        @test adaptive ≈ tensor rtol=2e-8
        estimated = Models.PNJLIntegrals.calculate_log_sum_rs_reduced_adaptive_with_error(
            masses, Φ, Φbar, mu_vec, T_fm, xi;
            rtol=1e-10, atol=1e-12,
        )
        @test estimated.value == adaptive
        @test isfinite(estimated.error)
        @test estimated.error >= 0.0
    end

    @testset "策略和 RS 定义域校验" begin
        @test Models.PNJLIntegrals.validate_thermal_quadrature_policy(:tensor_gauss) === :tensor_gauss
        @test Models.PNJLIntegrals.validate_thermal_quadrature_policy(:rs_reduced_adaptive) === :rs_reduced_adaptive
        @test_throws ArgumentError Models.PNJLIntegrals.validate_thermal_quadrature_policy(:unknown)
        @test_throws ArgumentError Models.PNJLIntegrals.validate_rs_anisotropy(-1.0)
        @test_throws ArgumentError Models.PNJLIntegrals.validate_rs_anisotropy(-1.1)
        @test_throws ArgumentError Models.PNJLIntegrals.validate_rs_anisotropy(Inf)
        @test_throws ArgumentError Models.PNJLIntegrals.validate_thermal_quadrature_controls(Inf, 1e-10, 100)
        @test_throws ArgumentError Models.PNJLIntegrals.validate_thermal_quadrature_controls(1e-8, Inf, 100)
        @test Models.PNJLIntegrals.rs_anisotropy_measure_factor(-0.5) ≈ sqrt(2.0)
    end

    @testset "低温和严格零温固定态极限" begin
        masses = SVector(0.3, 0.31, 0.5)
        mu_vec = SVector(0.4, 0.35, 0.2)
        args = (masses, 0.3, 0.25, mu_vec)
        low_T = Models.PNJLIntegrals.calculate_log_sum_rs_reduced_adaptive(
            args..., 1e-4, 0.2; rtol=1e-8, atol=1e-10,
        )
        zero_T = Models.PNJLIntegrals.calculate_log_sum_rs_reduced_adaptive(
            args..., 0.0, 0.2; rtol=1e-8, atol=1e-10,
        )
        @test isfinite(low_T)
        @test isfinite(zero_T)
        @test low_T ≈ zero_T rtol=2e-5 atol=1e-10

        derivative = ForwardDiff.derivative(0.2) do xi
            Models.PNJLIntegrals.calculate_log_sum_rs_reduced_adaptive(
                args..., 0.1, xi; rtol=1e-8, atol=1e-10,
            )
        end
        @test isfinite(derivative)
    end
end
