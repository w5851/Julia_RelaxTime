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
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))
using Test
using Base: time_ns
using FastGaussQuadrature: gausslegendre

const _CONSTANTS_PNJL_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "constants", "Constants_PNJL.jl"))
if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, _CONSTANTS_PNJL_PATH)
end
const _GAUSS_LEGENDRE_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "integration", "GaussLegendre.jl"))
if !isdefined(Main, :GaussLegendre)
    Base.include(Main, _GAUSS_LEGENDRE_PATH)
end
const _QUARK_DISTRIBUTION_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "models", "pnjl_physics", "QuarkDistribution.jl"))
if !isdefined(Main, :PNJLQuarkDistributions)
    Base.include(Main, _QUARK_DISTRIBUTION_PATH)
end
const _ONE_LOOP_INTEGRALS_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "OneLoopIntegrals.jl"))
if !isdefined(Main, :OneLoopIntegrals)
    Base.include(Main, _ONE_LOOP_INTEGRALS_PATH)
end
using Main.OneLoopIntegrals: B0, B0_retarded, A

const Λ_INV_FM = Main.OneLoopIntegrals.Λ_inv_fm

const DEFAULT_P_MAX = 20.0
const DEFAULT_GAUSS_POINTS = 128

function fixed_gl_oracle(f, a::Float64, b::Float64; n::Int=512)
    nodes, weights = gausslegendre(n)
    half = (b - a) / 2
    center = (a + b) / 2
    total = 0.0
    @inbounds for i in eachindex(nodes)
        total += weights[i] * f(center + half * nodes[i])
    end
    return half * total
end

function build_gauss_legendre_nodes()
    return OneLoopIntegrals.gauleg(0.0, DEFAULT_P_MAX, DEFAULT_GAUSS_POINTS)
end

TEST_PARAMS = (
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
        # 注意：当前实现不再接受 `rtol` keyword；这里重复调用用于一致性检查。
        result_repeat = B0(λ, 0.0, m1, μ1, m2, μ2, T)
        mirrored = B0(-λ, 0.0, m2, μ2, m1, μ1, T)
        sample_complex = B0(-0.05, 0.0, m1, μ1, m2, μ2, T; Φ=TEST_PARAMS.Φ, Φbar=TEST_PARAMS.Φbar)

        @info "B0(k=0) 对比" default=result_default repeat=result_repeat mirrored_swap=mirrored
        @info "B0(k=0) 虚部示例" sample_complex

        @test result_default isa NTuple{2, Float64}
        @test all(isfinite, result_default)
        @test isapprox(result_default[1], result_repeat[1]; rtol=1e-12, atol=0.0)
        @test isapprox(result_default[2], result_repeat[2]; rtol=1e-12, atol=0.0)
        @test isapprox(result_default[1], mirrored[1]; rtol=1e-8, atol=1e-10)
        @test isapprox(result_default[2], mirrored[2]; rtol=1e-8, atol=1e-10)
        @test abs(sample_complex[2]) > 0
    end

    @testset "k > 0 branch" begin
        k = 0.27
        result_default = B0(λ, k, m1, μ1, m2, μ2, T; Φ=TEST_PARAMS.Φ, Φbar=TEST_PARAMS.Φbar)
        # 注意：当前实现不再接受 `rtol` keyword；这里重复调用用于一致性检查。
        result_repeat = B0(λ, k, m1, μ1, m2, μ2, T; Φ=TEST_PARAMS.Φ, Φbar=TEST_PARAMS.Φbar)
        sample_complex = B0(-1.0, 0.1, m1, μ1, m2, μ2, T; Φ=TEST_PARAMS.Φ, Φbar=TEST_PARAMS.Φbar)

        @info "B0(k>0) 对比" default=result_default repeat=result_repeat
        @info "B0(k>0) 虚部示例" sample_complex

        k_small = 1.0e-4
        near_zero = B0(λ, k_small, m1, μ1, m2, μ2, T)
        zero_branch = B0(λ, 0.0, m1, μ1, m2, μ2, T)

        @test all(isfinite, result_default)
        @test isapprox(result_default[1], result_repeat[1]; rtol=1e-12, atol=0.0)
        @test isapprox(result_default[2], result_repeat[2]; rtol=1e-12, atol=0.0)
        @test isapprox(near_zero[1], zero_branch[1]; rtol=5e-3, atol=1e-4)
        @test isapprox(near_zero[2], zero_branch[2]; rtol=5e-3, atol=1e-4)
        @test abs(sample_complex[2]) > 0
    end

    @testset "Performance sanity check" begin
        if get(ENV, "RUN_UNIT_PERF", "0") == "1"
            iterations = parse(Int, get(ENV, "UNIT_B0_BENCH_ITERS", "2000"))
            avg_ms, sum_acc = benchmark_B0_call(iterations=iterations)
            @info "B0 单次调用平均耗时 (毫秒)" avg_ms iterations
            @test isfinite(sum_acc)
            @test avg_ms < 200.0
        else
            @test true
        end
    end
end

@testset "OneLoopIntegrals.B0_retarded" begin
    params = TEST_PARAMS

    @testset "upper-half-plane ordered continuation" begin
        for q in (0.0, 0.27)
            forward = B0_retarded(
                params.λ, q, params.m1, params.μ1, params.m2, params.μ2, params.T;
                Φ=params.Φ,
                Φbar=params.Φbar,
                eta_inv_fm=0.01,
                energy_nodes=128,
            )
            reverse = B0_retarded(
                -params.λ, q, params.m2, params.μ2, params.m1, params.μ1, params.T;
                Φ=params.Φ,
                Φbar=params.Φbar,
                eta_inv_fm=0.01,
                energy_nodes=128,
            )

            @test forward isa ComplexF64
            @test isfinite(real(forward)) && isfinite(imag(forward))
            @test forward ≈ conj(reverse) rtol=1e-12 atol=1e-12
        end
    end

    @testset "energy-node convergence away from a sharp cut" begin
        coarse = B0_retarded(
            params.λ, 0.27, params.m1, params.μ1, params.m2, params.μ2, params.T;
            Φ=params.Φ,
            Φbar=params.Φbar,
            eta_inv_fm=0.01,
            energy_nodes=64,
        )
        fine = B0_retarded(
            params.λ, 0.27, params.m1, params.μ1, params.m2, params.μ2, params.T;
            Φ=params.Φ,
            Φbar=params.Φbar,
            eta_inv_fm=0.01,
            energy_nodes=128,
        )
        @test coarse ≈ fine rtol=1e-5 atol=3e-6
    end

    @test_throws ArgumentError B0_retarded(params.λ, -0.1, params.m1, params.μ1, params.m2, params.μ2, params.T)
    @test_throws ArgumentError B0_retarded(params.λ, 0.1, -params.m1, params.μ1, params.m2, params.μ2, params.T)
    @test_throws ArgumentError B0_retarded(params.λ, 0.1, params.m1, params.μ1, params.m2, params.μ2, params.T; eta_inv_fm=0.0)
    @test_throws ArgumentError B0_retarded(params.λ, 0.1, params.m1, params.μ1, params.m2, params.μ2, params.T; energy_nodes=3)
end

@testset "OneLoopIntegrals.const_integral_term_A" begin
    m = TEST_PARAMS.m1
    integrand(p) = p^2 / sqrt(p^2 + m^2)
    numeric = fixed_gl_oracle(integrand, 0.0, Λ_INV_FM)
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

    dist_integral = fixed_gl_oracle(integrand, 0.0, DEFAULT_P_MAX)
    expected = 4.0 * (-OneLoopIntegrals.const_integral_term_A(m) + dist_integral)

    @test isapprox(result, expected; rtol=5e-5, atol=1e-6)

    nodes_fine, weights_fine = OneLoopIntegrals.gauleg(0.0, DEFAULT_P_MAX, DEFAULT_GAUSS_POINTS * 2)
    result_fine = A(m, μ, T, Φ, Φbar, nodes_fine, weights_fine)
    @test isapprox(result, result_fine; rtol=5e-5, atol=1e-6)
end

println("\n" * "="^70)
println("OneLoopIntegrals 模块测试完成！")
println("="^70)


