"""
测试 OneLoopIntegrals 模块

运行方式：
```julia
using Pkg
Pkg.activate(".")
include("test/test_oneloopintegrals.jl")
```
"""


using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Test
using Base: time_ns
using QuadGK: quadgk

include("../src/relaxtime/OneLoopIntegrals.jl")
using .OneLoopIntegrals: B0, A

const Λ_INV_FM = OneLoopIntegrals.Λ_inv_fm

const DEFAULT_P_MAX = 20.0
const DEFAULT_GAUSS_POINTS = 128

function build_gauss_legendre_nodes()
    return OneLoopIntegrals.gauleg(0.0, DEFAULT_P_MAX, DEFAULT_GAUSS_POINTS)
end

const TEST_PARAMS = (
    λ = 0.45,
    m1 = 0.24,
    μ1 = 0.12,
    m2 = 0.38,
    μ2 = -0.05,
    T = 0.17,
    Φ = 0.15,
    Φbar = 0.15,
)

# 辅助函数：粗略基准测试，输出单次调用的平均耗时（毫秒）
function benchmark_B0_call(; iterations::Int=8)
    λ = 0.62
    k = 0.35
    m1 = 0.28
    μ1 = 0.08
    m2 = 0.40
    μ2 = -0.04
    T = 0.18
    Φ = 0.20
    Φbar = 0.20

    GC.gc() # 减少基准测试的波动
    accumulator = 0.0
    start = time_ns()
    for _ in 1:iterations
        real_part, imag_part = B0(λ, k, m1, μ1, m2, μ2, T; Φ=Φ, Φbar=Φbar)
        accumulator += real_part + imag_part
    end
    elapsed_ms = (time_ns() - start) / iterations / 1.0e6
    return elapsed_ms, accumulator
end

@testset "OneLoopIntegrals.B0" begin
    λ = TEST_PARAMS.λ
    m1 = TEST_PARAMS.m1
    μ1 = TEST_PARAMS.μ1
    m2 = TEST_PARAMS.m2
    μ2 = TEST_PARAMS.μ2
    T = TEST_PARAMS.T

    @testset "k = 0 branch" begin
        result_default = B0(λ, 0.0, m1, μ1, m2, μ2, T)
        result_tight = B0(λ, 0.0, m1, μ1, m2, μ2, T; rtol=1e-6)
        mirrored = B0(-λ, 0.0, m2, μ2, m1, μ1, T)
        sample_complex = B0(-0.05, 0.0, m1, μ1, m2, μ2, T; Φ=TEST_PARAMS.Φ, Φbar=TEST_PARAMS.Φbar)

        @info "B0(k=0) 不同 rtol 对比" rtol_default=result_default rtol_1e6=result_tight mirrored_swap=mirrored
        @info "B0(k=0) 虚部示例" sample_complex

        @test result_default isa NTuple{2, Float64}
        @test all(isfinite, result_default)
        @test isapprox(result_default[1], result_tight[1]; rtol=5e-3, atol=1e-8)
        @test isapprox(result_default[2], result_tight[2]; rtol=5e-3, atol=1e-8)
        @test isapprox(result_default[1], mirrored[1]; rtol=1e-8, atol=1e-10)
        @test isapprox(result_default[2], mirrored[2]; rtol=1e-8, atol=1e-10)
        @test abs(sample_complex[2]) > 0
    end

    @testset "k > 0 branch" begin
        k = 0.27
        result_default = B0(λ, k, m1, μ1, m2, μ2, T; Φ=TEST_PARAMS.Φ, Φbar=TEST_PARAMS.Φbar)
        result_tight = B0(λ, k, m1, μ1, m2, μ2, T; Φ=TEST_PARAMS.Φ, Φbar=TEST_PARAMS.Φbar, rtol=1e-6)
        sample_complex = B0(-1.0, 0.1, m1, μ1, m2, μ2, T; Φ=TEST_PARAMS.Φ, Φbar=TEST_PARAMS.Φbar)

        @info "B0(k>0) 不同 rtol 对比" rtol_default=result_default rtol_1e6=result_tight
        @info "B0(k>0) 虚部示例" sample_complex

        k_small = 1.0e-4
        near_zero = B0(λ, k_small, m1, μ1, m2, μ2, T)
        zero_branch = B0(λ, 0.0, m1, μ1, m2, μ2, T)

        @test all(isfinite, result_default)
        @test isapprox(result_default[1], result_tight[1]; rtol=1e-3, atol=1e-7)
        @test isapprox(result_default[2], result_tight[2]; rtol=1e-3, atol=1e-7)
        @test isapprox(near_zero[1], zero_branch[1]; rtol=5e-3, atol=1e-4)
        @test isapprox(near_zero[2], zero_branch[2]; rtol=5e-3, atol=1e-4)
        @test abs(sample_complex[2]) > 0
    end

    @testset "Performance sanity check" begin
        avg_ms, sum_acc = benchmark_B0_call(iterations=50000)
        @info "B0 单次调用平均耗时 (毫秒)" avg_ms
        @test isfinite(sum_acc)
        @test avg_ms < 200.0
    end
end

@testset "OneLoopIntegrals.const_integral_term_A" begin
    m = TEST_PARAMS.m1
    integrand(p) = p^2 / sqrt(p^2 + m^2)
    numeric, _ = quadgk(integrand, 0.0, Λ_INV_FM; rtol=1e-10, atol=1e-12)
    analytic = OneLoopIntegrals.const_integral_term_A(m)
    @test isapprox(analytic, numeric; rtol=1e-9, atol=1e-11)
end

@testset "OneLoopIntegrals.A" begin
    nodes, weights = build_gauss_legendre_nodes()

    m = TEST_PARAMS.m1
    μ = TEST_PARAMS.μ1
    T = TEST_PARAMS.T
    Φ = TEST_PARAMS.Φ
    Φbar = TEST_PARAMS.Φbar

    result = A(m, μ, T, Φ, Φbar, nodes, weights)

    integrand(p) = begin
        E = sqrt(p^2 + m^2)
        dist = OneLoopIntegrals.quark_distribution(E, μ, T, Φ, Φbar) +
            OneLoopIntegrals.antiquark_distribution(E, μ, T, Φ, Φbar)
        return p^2 / E * dist
    end

    dist_integral, _ = quadgk(integrand, 0.0, DEFAULT_P_MAX; rtol=1e-8, atol=1e-10)
    expected = 4.0 * (-OneLoopIntegrals.const_integral_term_A(m) + dist_integral)

    @test isapprox(result, expected; rtol=5e-5, atol=1e-6)

    nodes_fine, weights_fine = OneLoopIntegrals.gauleg(0.0, DEFAULT_P_MAX, DEFAULT_GAUSS_POINTS * 2)
    result_fine = A(m, μ, T, Φ, Φbar, nodes_fine, weights_fine)
    @test isapprox(result, result_fine; rtol=5e-5, atol=1e-6)
end

println("\n" * "="^70)
println("OneLoopIntegrals 模块测试完成！")
println("="^70)
