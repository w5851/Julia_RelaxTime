"""
测试 Polarization 模块

运行方式：
```julia
using Pkg
Pkg.activate(".")
include("test/test_polarization.jl")
```
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using Test
using Base: time_ns

include("../src/Constants_PNJL.jl")
include("../src/relaxtime/OneLoopIntegrals.jl")
include("../src/relaxtime/Polarization.jl")
using .Constants_PNJL: N_color
using .OneLoopIntegrals: B0, A, gauleg
using .Polarization: polarization


const POL_TEST_PARAMS = (
    k0 = 0.45,
    k_norm = 0.3,
    m1 = 0.25,
    μ1 = 0.12,
    m2 = 0.36,
    μ2 = -0.05,
    T = 0.17,
    Φ = 0.15,
    Φbar = 0.15,
)

const GAUSS_POINTS = 32
const P_MAX = 20.0
const NUM_S_QUARK_DEFAULT = 1

function prepare_A_values(m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64)
    nodes, weights = gauleg(0.0, P_MAX, GAUSS_POINTS)
    return A(m, μ, T, Φ, Φbar, nodes, weights)
end

function benchmark_polarization_call(channel::Symbol; iterations::Int=4)
    params = POL_TEST_PARAMS
    nodes, weights = gauleg(0.0, P_MAX, GAUSS_POINTS)
    A1 = A(params.m1, params.μ1, params.T, params.Φ, params.Φbar, nodes, weights)
    A2 = A(params.m2, params.μ2, params.T, params.Φ, params.Φbar, nodes, weights)

    GC.gc()
    real_acc = 0.0
    imag_acc = 0.0
    start = time_ns()
    for _ in 1:iterations
        value = polarization(channel, params.k0, params.k_norm, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, A1, A2, NUM_S_QUARK_DEFAULT)
        real_acc += value[1]
        imag_acc += value[2]
    end
    elapsed_ms = (time_ns() - start) / iterations / 1.0e6
    return elapsed_ms, (real_acc, imag_acc)
end

@testset "Polarization" begin
    params = POL_TEST_PARAMS
    nodes, weights = gauleg(0.0, P_MAX, GAUSS_POINTS)
    nodes_fine, weights_fine = gauleg(0.0, P_MAX, GAUSS_POINTS * 2)

    A1 = A(params.m1, params.μ1, params.T, params.Φ, params.Φbar, nodes, weights)
    A2 = A(params.m2, params.μ2, params.T, params.Φ, params.Φbar, nodes, weights)
    A1_fine = A(params.m1, params.μ1, params.T, params.Φ, params.Φbar, nodes_fine, weights_fine)
    A2_fine = A(params.m2, params.μ2, params.T, params.Φ, params.Φbar, nodes_fine, weights_fine)

    ΠP = polarization(:P, params.k0, params.k_norm, params.m1, params.m2,
        params.μ1, params.μ2, params.T, params.Φ, params.Φbar, A1, A2, NUM_S_QUARK_DEFAULT)
    ΠS = polarization(:S, params.k0, params.k_norm, params.m1, params.m2,
        params.μ1, params.μ2, params.T, params.Φ, params.Φbar, A1, A2, NUM_S_QUARK_DEFAULT)

    ΠP_ref = polarization(:P, params.k0, params.k_norm, params.m1, params.m2,
        params.μ1, params.μ2, params.T, params.Φ, params.Φbar, A1_fine, A2_fine, NUM_S_QUARK_DEFAULT)
    ΠS_ref = polarization(:S, params.k0, params.k_norm, params.m1, params.m2,
        params.μ1, params.μ2, params.T, params.Φ, params.Φbar, A1_fine, A2_fine, NUM_S_QUARK_DEFAULT)

    @info "Polarization(:P) 对比不同 A 精度" coarse=ΠP fine=ΠP_ref
    @info "Polarization(:S) 对比不同 A 精度" coarse=ΠS fine=ΠS_ref

    @test all(isfinite, ΠP)
    @test all(isfinite, ΠS)
    @test isapprox(ΠP[1], ΠP_ref[1]; rtol=5e-4, atol=1e-6)
    @test isapprox(ΠP[2], ΠP_ref[2]; rtol=5e-4, atol=1e-6)
    @test isapprox(ΠS[1], ΠS_ref[1]; rtol=5e-4, atol=1e-6)
    @test isapprox(ΠS[2], ΠS_ref[2]; rtol=5e-4, atol=1e-6)

    @info "Polarization(:P) 示例" value=ΠP
    @info "Polarization(:S) 示例" value=ΠS

    λ_imag = -1.0
    k_norm_imag = 0.1
    k0_imag = λ_imag - (params.μ1 - params.μ2)
    ΠP_imag = polarization(:P, k0_imag, k_norm_imag, params.m1, params.m2,
        params.μ1, params.μ2, params.T, params.Φ, params.Φbar, A1, A2, NUM_S_QUARK_DEFAULT)
    ΠS_imag = polarization(:S, k0_imag, k_norm_imag, params.m1, params.m2,
        params.μ1, params.μ2, params.T, params.Φ, params.Φbar, A1, A2, NUM_S_QUARK_DEFAULT)

    @info "Polarization(:P) 虚部示例" params=(λ=λ_imag, k=k_norm_imag, k0=k0_imag) value=ΠP_imag
    @info "Polarization(:S) 虚部示例" params=(λ=λ_imag, k=k_norm_imag, k0=k0_imag) value=ΠS_imag

    @test abs(ΠP_imag[2]) > 0
    @test abs(ΠS_imag[2]) > 0

    @testset "num_s_quark = 1 对称平均" begin
        num_s = 1
        ΠP_sym = polarization(:P, params.k0, params.k_norm, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, A1, A2, num_s)
        ΠS_sym = polarization(:S, params.k0, params.k_norm, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, A1, A2, num_s)

        λ = params.k0 + params.μ1 - params.μ2
        λ_extra = -params.k0 + params.μ1 - params.μ2
        B0_main = B0(λ, params.k_norm, params.m1, params.μ1, params.m2, params.μ2, params.T; Φ=params.Φ, Φbar=params.Φbar)
        B0_extra = B0(λ_extra, params.k_norm, params.m1, params.μ1, params.m2, params.μ2, params.T; Φ=params.Φ, Φbar=params.Φbar)
        B0_avg_real = 0.5 * (B0_main[1] + B0_extra[1])
        B0_avg_imag = 0.5 * (B0_main[2] + B0_extra[2])

        prefactor_P = params.k_norm^2 - λ^2 + (params.m1 - params.m2)^2
        prefactor_S = params.k_norm^2 - λ^2 + (params.m1 + params.m2)^2
        A_sum = A1 + A2
    factor = -N_color / (8 * pi^2)

        expected_P = (factor * (A_sum + prefactor_P * B0_avg_real), factor * (prefactor_P * B0_avg_imag))
        expected_S = (factor * (A_sum + prefactor_S * B0_avg_real), factor * (prefactor_S * B0_avg_imag))

        @info "Polarization(:P, num_s=1)" average_B0=(B0_avg_real, B0_avg_imag) value=ΠP_sym
        @info "Polarization(:S, num_s=1)" average_B0=(B0_avg_real, B0_avg_imag) value=ΠS_sym

        @test isapprox(ΠP_sym[1], expected_P[1]; rtol=1e-8, atol=1e-10)
        @test isapprox(ΠP_sym[2], expected_P[2]; rtol=1e-8, atol=1e-10)
        @test isapprox(ΠS_sym[1], expected_S[1]; rtol=1e-8, atol=1e-10)
        @test isapprox(ΠS_sym[2], expected_S[2]; rtol=1e-8, atol=1e-10)
    end

    @testset "Invalid channel" begin
        @test_throws ArgumentError polarization(:X, params.k0, params.k_norm, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, A1, A2, NUM_S_QUARK_DEFAULT)
    end

    @testset "Performance sanity" begin
        avg_ms_P, acc_P = benchmark_polarization_call(:P; iterations=2000)
        avg_ms_S, acc_S = benchmark_polarization_call(:S; iterations=2000)
        @info "Polarization(:P) 单次调用平均耗时 (毫秒)" avg_ms_P
        @info "Polarization(:S) 单次调用平均耗时 (毫秒)" avg_ms_S
        @test isfinite(acc_P[1]) && isfinite(acc_P[2])
        @test isfinite(acc_S[1]) && isfinite(acc_S[2])
        @test avg_ms_P < 0.5
        @test avg_ms_S < 0.5
    end
end

println("\n" * "="^70)
println("Polarization 模块测试完成！")
println("="^70)