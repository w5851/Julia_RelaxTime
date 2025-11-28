"""
测试 PolarizationAniso 模块

运行方式：
```julia
using Pkg
Pkg.activate(".")
include("test/test_polarization_aniso.jl")
```
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
using Test
using Base: time_ns
using Printf

include("../../src/Constants_PNJL.jl")
include("../../src/relaxtime/OneLoopIntegrals.jl")
include("../../src/relaxtime/OneLoopIntegralsAniso.jl")
include("../../src/relaxtime/PolarizationAniso.jl")
using .Constants_PNJL: N_color
using .OneLoopIntegrals: B0, A, gauleg, EPS_SEGMENT
using .OneLoopIntegralsCorrection: B0_correction
using .PolarizationAniso: polarization_aniso


const POL_ANISO_TEST_PARAMS = (
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

# 与 test_b0_correction.jl 相同的参数，用于验证一致性
const B0_CORRECTION_TEST_PARAMS = (
    λ = 0.45,
    k = 0.30,
    m1 = 0.24,
    m2 = 0.38,
    μ1 = 0.12,
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

function benchmark_polarization_aniso_call(channel::Symbol, ξ::Float64, A1::Float64, A2::Float64; iterations::Int=4, warmup::Int=100)
    params = POL_ANISO_TEST_PARAMS

    # 预热阶段
    for _ in 1:warmup
        polarization_aniso(channel, params.k0, params.k_norm, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, ξ, A1, A2, NUM_S_QUARK_DEFAULT)
    end

    GC.gc()
    real_acc = 0.0
    imag_acc = 0.0
    start = time_ns()
    for _ in 1:iterations
        value = polarization_aniso(channel, params.k0, params.k_norm, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, ξ, A1, A2, NUM_S_QUARK_DEFAULT)
        real_acc += value[1]
        imag_acc += value[2]
    end
    elapsed_ms = (time_ns() - start) / iterations / 1.0e6
    return elapsed_ms, (real_acc, imag_acc)
end

@testset "PolarizationAniso" begin
    params = POL_ANISO_TEST_PARAMS
    nodes, weights = gauleg(0.0, P_MAX, GAUSS_POINTS)
    nodes_fine, weights_fine = gauleg(0.0, P_MAX, GAUSS_POINTS * 2)

    A1 = A(params.m1, params.μ1, params.T, params.Φ, params.Φbar, nodes, weights)
    A2 = A(params.m2, params.μ2, params.T, params.Φ, params.Φbar, nodes, weights)
    A1_fine = A(params.m1, params.μ1, params.T, params.Φ, params.Φbar, nodes_fine, weights_fine)
    A2_fine = A(params.m2, params.μ2, params.T, params.Φ, params.Φbar, nodes_fine, weights_fine)

    @testset "基本功能测试" begin
        # 测试各向同性情况 (ξ = 0)
        ΠP_iso = polarization_aniso(:P, params.k0, params.k_norm, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, 0.0, A1, A2, NUM_S_QUARK_DEFAULT)
        ΠS_iso = polarization_aniso(:S, params.k0, params.k_norm, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, 0.0, A1, A2, NUM_S_QUARK_DEFAULT)

        @info "PolarizationAniso(:P, ξ=0) 各向同性" value=ΠP_iso
        @info "PolarizationAniso(:S, ξ=0) 各向同性" value=ΠS_iso

        @test all(isfinite, ΠP_iso)
        @test all(isfinite, ΠS_iso)

        # 测试各向异性情况 (ξ = 0.5)
        ΠP_aniso = polarization_aniso(:P, params.k0, params.k_norm, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, 0.5, A1, A2, NUM_S_QUARK_DEFAULT)
        ΠS_aniso = polarization_aniso(:S, params.k0, params.k_norm, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, 0.5, A1, A2, NUM_S_QUARK_DEFAULT)

        @info "PolarizationAniso(:P, ξ=0.5) 各向异性" value=ΠP_aniso
        @info "PolarizationAniso(:S, ξ=0.5) 各向异性" value=ΠS_aniso

        @test all(isfinite, ΠP_aniso)
        @test all(isfinite, ΠS_aniso)
    end

    @testset "ξ 参数影响测试" begin
        # 测试 ξ 从 -1 到 1 的影响
        ξ_values = [-1.0, -0.8, -0.5, -0.2, 0.0, 0.2, 0.5, 0.8, 1.0]
        
        println("\n" * "="^70)
        println("测试各向异性参数 ξ 对极化函数的影响")
        println("="^70)
        println("ξ 取值范围: -1 到 1，0 对应各向同性")
        println("-"^70)
        
        results_P = Float64[]
        results_S = Float64[]
        
        for ξ in ξ_values
            ΠP = polarization_aniso(:P, params.k0, params.k_norm, params.m1, params.m2,
                params.μ1, params.μ2, params.T, params.Φ, params.Φbar, ξ, A1, A2, NUM_S_QUARK_DEFAULT)
            ΠS = polarization_aniso(:S, params.k0, params.k_norm, params.m1, params.m2,
                params.μ1, params.μ2, params.T, params.Φ, params.Φbar, ξ, A1, A2, NUM_S_QUARK_DEFAULT)
            
            push!(results_P, ΠP[1])
            push!(results_S, ΠS[1])
            
            @test all(isfinite, ΠP)
            @test all(isfinite, ΠS)
            
            println(@sprintf("ξ = %5.2f | Π_P: Re=%12.6e, Im=%12.6e | Π_S: Re=%12.6e, Im=%12.6e", 
                    ξ, ΠP[1], ΠP[2], ΠS[1], ΠS[2]))
        end
        
        println("-"^70)
        
        # 检查 ξ 的影响是否显著
        max_P = maximum(results_P)
        min_P = minimum(results_P)
        max_S = maximum(results_S)
        min_S = minimum(results_S)
        
        variation_P = (max_P - min_P) / abs(min_P)
        variation_S = (max_S - min_S) / abs(min_S)
        
        @info "Π_P 实部随 ξ 的相对变化" min=min_P max=max_P variation=variation_P
        @info "Π_S 实部随 ξ 的相对变化" min=min_S max=max_S variation=variation_S
        
        # 各向异性应该产生可观测的效应（至少万分之一的变化）
        @test variation_P > 0.0001 || variation_S > 0.0001
        
        println("="^70)
    end

    @testset "ξ = 0 与各向同性一致性" begin
        # 当 ξ 非常小时，各向异性修正应该可以忽略
        ξ_small = EPS_SEGMENT / 2.0
        ΠP_small_xi = polarization_aniso(:P, params.k0, params.k_norm, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, ξ_small, A1, A2, NUM_S_QUARK_DEFAULT)
        ΠS_small_xi = polarization_aniso(:S, params.k0, params.k_norm, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, ξ_small, A1, A2, NUM_S_QUARK_DEFAULT)
        
        # 与 ξ = 0 的结果对比
        ΠP_zero = polarization_aniso(:P, params.k0, params.k_norm, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, 0.0, A1, A2, NUM_S_QUARK_DEFAULT)
        ΠS_zero = polarization_aniso(:S, params.k0, params.k_norm, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, 0.0, A1, A2, NUM_S_QUARK_DEFAULT)
        
        @info "PolarizationAniso(:P) 小 ξ 与 ξ=0 对比" small_xi=ΠP_small_xi zero_xi=ΠP_zero
        @info "PolarizationAniso(:S) 小 ξ 与 ξ=0 对比" small_xi=ΠS_small_xi zero_xi=ΠS_zero
        
        # 小 ξ 时结果应该与 ξ=0 非常接近
        @test isapprox(ΠP_small_xi[1], ΠP_zero[1]; rtol=1e-10, atol=1e-12)
        @test isapprox(ΠP_small_xi[2], ΠP_zero[2]; rtol=1e-10, atol=1e-12)
        @test isapprox(ΠS_small_xi[1], ΠS_zero[1]; rtol=1e-10, atol=1e-12)
        @test isapprox(ΠS_small_xi[2], ΠS_zero[2]; rtol=1e-10, atol=1e-12)
    end

    @testset "精度对比测试" begin
        ΠP_coarse = polarization_aniso(:P, params.k0, params.k_norm, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, 0.3, A1, A2, NUM_S_QUARK_DEFAULT)
        ΠS_coarse = polarization_aniso(:S, params.k0, params.k_norm, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, 0.3, A1, A2, NUM_S_QUARK_DEFAULT)

        ΠP_fine = polarization_aniso(:P, params.k0, params.k_norm, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, 0.3, A1_fine, A2_fine, NUM_S_QUARK_DEFAULT)
        ΠS_fine = polarization_aniso(:S, params.k0, params.k_norm, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, 0.3, A1_fine, A2_fine, NUM_S_QUARK_DEFAULT)

        @info "PolarizationAniso(:P, ξ=0.3) 对比不同 A 精度" coarse=ΠP_coarse fine=ΠP_fine
        @info "PolarizationAniso(:S, ξ=0.3) 对比不同 A 精度" coarse=ΠS_coarse fine=ΠS_fine

        @test isapprox(ΠP_coarse[1], ΠP_fine[1]; rtol=5e-4, atol=1e-6)
        @test isapprox(ΠP_coarse[2], ΠP_fine[2]; rtol=5e-4, atol=1e-6)
        @test isapprox(ΠS_coarse[1], ΠS_fine[1]; rtol=5e-4, atol=1e-6)
        @test isapprox(ΠS_coarse[2], ΠS_fine[2]; rtol=5e-4, atol=1e-6)
    end

    @testset "虚部测试" begin
        λ_imag = -1.0
        k_norm_imag = 0.1
        k0_imag = λ_imag - (params.μ1 - params.μ2)
        
        ΠP_imag_iso = polarization_aniso(:P, k0_imag, k_norm_imag, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, 0.0, A1, A2, NUM_S_QUARK_DEFAULT)
        ΠS_imag_iso = polarization_aniso(:S, k0_imag, k_norm_imag, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, 0.0, A1, A2, NUM_S_QUARK_DEFAULT)
        
        ΠP_imag_aniso = polarization_aniso(:P, k0_imag, k_norm_imag, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, 0.6, A1, A2, NUM_S_QUARK_DEFAULT)
        ΠS_imag_aniso = polarization_aniso(:S, k0_imag, k_norm_imag, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, 0.6, A1, A2, NUM_S_QUARK_DEFAULT)

        @info "PolarizationAniso(:P) 虚部示例 (ξ=0)" params=(λ=λ_imag, k=k_norm_imag, k0=k0_imag) value=ΠP_imag_iso
        @info "PolarizationAniso(:S) 虚部示例 (ξ=0)" params=(λ=λ_imag, k=k_norm_imag, k0=k0_imag) value=ΠS_imag_iso
        @info "PolarizationAniso(:P) 虚部示例 (ξ=0.6)" params=(λ=λ_imag, k=k_norm_imag, k0=k0_imag) value=ΠP_imag_aniso
        @info "PolarizationAniso(:S) 虚部示例 (ξ=0.6)" params=(λ=λ_imag, k=k_norm_imag, k0=k0_imag) value=ΠS_imag_aniso

        @test abs(ΠP_imag_iso[2]) > 0
        @test abs(ΠS_imag_iso[2]) > 0
        @test abs(ΠP_imag_aniso[2]) > 0
        @test abs(ΠS_imag_aniso[2]) > 0
    end

    @testset "num_s_quark = 1 对称平均" begin
        num_s = 1
        ξ_test = 0.4
        
        ΠP_sym = polarization_aniso(:P, params.k0, params.k_norm, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, ξ_test, A1, A2, num_s)
        ΠS_sym = polarization_aniso(:S, params.k0, params.k_norm, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, ξ_test, A1, A2, num_s)

        λ = params.k0 + params.μ1 - params.μ2
        λ_extra = -params.k0 + params.μ1 - params.μ2
        
        # 计算 B0 及其修正
        B0_main = B0(λ, params.k_norm, params.m1, params.μ1, params.m2, params.μ2, params.T; Φ=params.Φ, Φbar=params.Φbar)
        B0_extra = B0(λ_extra, params.k_norm, params.m1, params.μ1, params.m2, params.μ2, params.T; Φ=params.Φ, Φbar=params.Φbar)
        
        B0_corr_main = B0_correction(λ, params.k_norm, params.m1, params.m2, params.μ1, params.μ2, params.T, params.Φ, params.Φbar, ξ_test)
        B0_corr_extra = B0_correction(λ_extra, params.k_norm, params.m1, params.m2, params.μ1, params.μ2, params.T, params.Φ, params.Φbar, ξ_test)
        
        B0_avg_real = 0.5 * ((B0_main[1] + B0_corr_main[1]) + (B0_extra[1] + B0_corr_extra[1]))
        B0_avg_imag = 0.5 * ((B0_main[2] + B0_corr_main[2]) + (B0_extra[2] + B0_corr_extra[2]))

        prefactor_P = params.k_norm^2 - λ^2 + (params.m1 - params.m2)^2
        prefactor_S = params.k_norm^2 - λ^2 + (params.m1 + params.m2)^2
        A_sum = A1 + A2
        factor = -N_color / (8 * pi^2)

        expected_P = (factor * (A_sum + prefactor_P * B0_avg_real), factor * (prefactor_P * B0_avg_imag))
        expected_S = (factor * (A_sum + prefactor_S * B0_avg_real), factor * (prefactor_S * B0_avg_imag))

        @info "PolarizationAniso(:P, num_s=1, ξ=$ξ_test)" average_B0=(B0_avg_real, B0_avg_imag) value=ΠP_sym
        @info "PolarizationAniso(:S, num_s=1, ξ=$ξ_test)" average_B0=(B0_avg_real, B0_avg_imag) value=ΠS_sym

        @test isapprox(ΠP_sym[1], expected_P[1]; rtol=1e-8, atol=1e-10)
        @test isapprox(ΠP_sym[2], expected_P[2]; rtol=1e-8, atol=1e-10)
        @test isapprox(ΠS_sym[1], expected_S[1]; rtol=1e-8, atol=1e-10)
        @test isapprox(ΠS_sym[2], expected_S[2]; rtol=1e-8, atol=1e-10)
    end

    @testset "Invalid channel" begin
        @test_throws ArgumentError polarization_aniso(:X, params.k0, params.k_norm, params.m1, params.m2,
            params.μ1, params.μ2, params.T, params.Φ, params.Φbar, 0.5, A1, A2, NUM_S_QUARK_DEFAULT)
    end

    @testset "Performance sanity" begin
        # 预计算 A 值（避免在性能测试中重复计算）
        nodes, weights = gauleg(0.0, P_MAX, GAUSS_POINTS)
        A1_perf = A(params.m1, params.μ1, params.T, params.Φ, params.Φbar, nodes, weights)
        A2_perf = A(params.m2, params.μ2, params.T, params.Φ, params.Φbar, nodes, weights)
        
        # 使用更多迭代次数以获得更准确的性能数据
        avg_ms_P_iso, acc_P_iso = benchmark_polarization_aniso_call(:P, 0.0, A1_perf, A2_perf; iterations=10000, warmup=500)
        avg_ms_S_iso, acc_S_iso = benchmark_polarization_aniso_call(:S, 0.0, A1_perf, A2_perf; iterations=10000, warmup=500)
        avg_ms_P_aniso, acc_P_aniso = benchmark_polarization_aniso_call(:P, 0.5, A1_perf, A2_perf; iterations=10000, warmup=500)
        avg_ms_S_aniso, acc_S_aniso = benchmark_polarization_aniso_call(:S, 0.5, A1_perf, A2_perf; iterations=10000, warmup=500)
        
        @info "PolarizationAniso(:P, ξ=0) 单次调用平均耗时 (毫秒)" avg_ms_P_iso
        @info "PolarizationAniso(:S, ξ=0) 单次调用平均耗时 (毫秒)" avg_ms_S_iso
        @info "PolarizationAniso(:P, ξ=0.5) 单次调用平均耗时 (毫秒)" avg_ms_P_aniso
        @info "PolarizationAniso(:S, ξ=0.5) 单次调用平均耗时 (毫秒)" avg_ms_S_aniso
        
        @test isfinite(acc_P_iso[1]) && isfinite(acc_P_iso[2])
        @test isfinite(acc_S_iso[1]) && isfinite(acc_S_iso[2])
        @test isfinite(acc_P_aniso[1]) && isfinite(acc_P_aniso[2])
        @test isfinite(acc_S_aniso[1]) && isfinite(acc_S_aniso[2])
        
        # 各向异性计算应该比各向同性慢一些（因为有额外的修正项）
        @test avg_ms_P_aniso >= avg_ms_P_iso
        @test avg_ms_S_aniso >= avg_ms_S_iso
        
        # 合理的性能范围（基于详细分析）
        @test avg_ms_P_iso < 0.05  # ξ=0 应该约 0.026ms
        @test avg_ms_S_iso < 0.05
        @test avg_ms_P_aniso < 0.15  # ξ=0.5 应该约 0.077ms
        @test avg_ms_S_aniso < 0.15
    end

    @testset "使用 B0_correction 测试参数验证一致性" begin
        # 使用与 test_b0_correction.jl 相同的参数
        # 这些参数已经在 B0_correction 文档中记录了详细的行为
        test_params = B0_CORRECTION_TEST_PARAMS
        
        # 从 λ = k0 + μ1 - μ2 反推 k0
        k0 = test_params.λ - (test_params.μ1 - test_params.μ2)
        
        nodes, weights = gauleg(0.0, P_MAX, GAUSS_POINTS)
        A1 = A(test_params.m1, test_params.μ1, test_params.T, test_params.Φ, test_params.Φbar, nodes, weights)
        A2 = A(test_params.m2, test_params.μ2, test_params.T, test_params.Φ, test_params.Φbar, nodes, weights)
        
        println("\n" * "="^80)
        println("使用 B0_correction 测试参数 (已记录在文档中)")
        println("="^80)
        println("参数设置：")
        println(@sprintf("  λ = %.3f, k = %.3f, k0 = %.3f", test_params.λ, test_params.k, k0))
        println(@sprintf("  m1 = %.3f, m2 = %.3f", test_params.m1, test_params.m2))
        println(@sprintf("  μ1 = %.3f, μ2 = %.3f, T = %.3f", test_params.μ1, test_params.μ2, test_params.T))
        println("-"^80)
        
        # 测试与文档中记录的 ξ 值范围一致
        ξ_values = [0.0, 0.05, 0.10, 0.20, 0.30, 0.50]
        
        println(@sprintf("%-10s %-22s %-22s %-22s %-22s", 
                        "ξ", "Π_P 实部", "Π_P 虚部", "Π_S 实部", "Π_S 虚部"))
        println("-"^80)
        
        results_P = []
        results_S = []
        
        for ξ in ξ_values
            ΠP = polarization_aniso(:P, k0, test_params.k, test_params.m1, test_params.m2,
                test_params.μ1, test_params.μ2, test_params.T, test_params.Φ, test_params.Φbar, 
                ξ, A1, A2, 0)
            ΠS = polarization_aniso(:S, k0, test_params.k, test_params.m1, test_params.m2,
                test_params.μ1, test_params.μ2, test_params.T, test_params.Φ, test_params.Φbar, 
                ξ, A1, A2, 0)
            
            push!(results_P, (ξ=ξ, real=ΠP[1], imag=ΠP[2]))
            push!(results_S, (ξ=ξ, real=ΠS[1], imag=ΠS[2]))
            
            println(@sprintf("%-10.2f %-22.10e %-22.10e %-22.10e %-22.10e", 
                            ξ, ΠP[1], ΠP[2], ΠS[1], ΠS[2]))
            
            @test all(isfinite, ΠP)
            @test all(isfinite, ΠS)
        end
        
        println("="^80)
        
        # 验证线性关系（与 B0_correction 文档中记录的一致）
        println("\n线性度验证（小 ξ 下）：")
        ξ1_idx = findfirst(r -> r.ξ == 0.05, results_P)
        ξ2_idx = findfirst(r -> r.ξ == 0.10, results_P)
        
        if !isnothing(ξ1_idx) && !isnothing(ξ2_idx)
            ratio_P = (results_P[ξ2_idx].real - results_P[1].real) / (results_P[ξ1_idx].real - results_P[1].real)
            ratio_S = (results_S[ξ2_idx].real - results_S[1].real) / (results_S[ξ1_idx].real - results_S[1].real)
            expected_ratio = 0.10 / 0.05
            
            println(@sprintf("  Π_P 实部: ξ=0.10 相对 ξ=0.05 的变化比值 = %.4f (预期 %.2f)", ratio_P, expected_ratio))
            println(@sprintf("  Π_S 实部: ξ=0.10 相对 ξ=0.05 的变化比值 = %.4f (预期 %.2f)", ratio_S, expected_ratio))
            
            # B0_correction 文档显示完美的 2.0 比值，polarization 应该也保持这个比例
            @test isapprox(ratio_P, expected_ratio; rtol=0.1)
            @test isapprox(ratio_S, expected_ratio; rtol=0.1)
        end
        
        # 计算相对变化
        println("\n相对于各向同性的变化百分比：")
        println(@sprintf("%-10s %-20s %-20s", "ξ", "Π_P 变化 (%)", "Π_S 变化 (%)"))
        println("-"^80)
        
        Π0_P = results_P[1].real
        Π0_S = results_S[1].real
        
        for i in 2:length(results_P)
            change_P = abs((results_P[i].real - Π0_P) / Π0_P) * 100
            change_S = abs((results_S[i].real - Π0_S) / Π0_S) * 100
            println(@sprintf("%-10.2f %-20.6f %-20.6f", results_P[i].ξ, change_P, change_S))
        end
        println("="^80 * "\n")
    end

    @testset "不同动量 k 下的 ξ 影响（对应 B0_correction 文档）" begin
        test_params = B0_CORRECTION_TEST_PARAMS
        k0 = test_params.λ - (test_params.μ1 - test_params.μ2)
        ξ_test = 0.2  # 与 B0_correction 文档中使用的 ξ 一致
        
        nodes, weights = gauleg(0.0, P_MAX, GAUSS_POINTS)
        A1 = A(test_params.m1, test_params.μ1, test_params.T, test_params.Φ, test_params.Φbar, nodes, weights)
        A2 = A(test_params.m2, test_params.μ2, test_params.T, test_params.Φ, test_params.Φbar, nodes, weights)
        
        k_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
        
        println("\n" * "="^80)
        println("不同 k 值下的 Polarization 行为 (ξ = 0.2, 对应 B0_correction 文档)")
        println("="^80)
        println(@sprintf("%-8s %-20s %-20s %-20s %-20s", 
                        "k", "Π_P 实部", "Π_P 虚部", "Π_S 实部", "Π_S 虚部"))
        println("-"^80)
        
        for k in k_values
            ΠP = polarization_aniso(:P, k0, k, test_params.m1, test_params.m2,
                test_params.μ1, test_params.μ2, test_params.T, test_params.Φ, test_params.Φbar, 
                ξ_test, A1, A2, 0)
            ΠS = polarization_aniso(:S, k0, k, test_params.m1, test_params.m2,
                test_params.μ1, test_params.μ2, test_params.T, test_params.Φ, test_params.Φbar, 
                ξ_test, A1, A2, 0)
            
            println(@sprintf("%-8.2f %-20.10e %-20.10e %-20.10e %-20.10e", 
                            k, ΠP[1], ΠP[2], ΠS[1], ΠS[2]))
            
            @test all(isfinite, ΠP)
            @test all(isfinite, ΠS)
            
            # 根据 B0_correction 文档，k=0.5 时虚部应该非零（λ² < k²）
            if k == 0.5
                @info "k=0.5 (λ² < k²) 虚部行为" Π_P_imag=ΠP[2] Π_S_imag=ΠS[2]
            end
        end
        println("="^80 * "\n")
    end

    @testset "多组物理参数的 ξ 依赖性测试" begin
        # 使用与 B0_correction 测试类似的多组参数
        test_cases = [
            (λ=0.3, k=0.2, m1=0.2, m2=0.3, μ1=0.1, μ2=0.0, T=0.15, name="参数组 1"),
            (λ=0.5, k=0.4, m1=0.25, m2=0.35, μ1=0.15, μ2=-0.1, T=0.18, name="参数组 2"),
            (λ=0.6, k=0.5, m1=0.3, m2=0.4, μ1=0.2, μ2=-0.05, T=0.20, name="参数组 3"),
        ]
        
        Φ = POL_ANISO_TEST_PARAMS.Φ
        Φbar = POL_ANISO_TEST_PARAMS.Φbar
        ξ = 0.2
        
        nodes, weights = gauleg(0.0, P_MAX, GAUSS_POINTS)
        
        println("\n" * "="^80)
        println("多组物理参数下的 Polarization (ξ = 0.2)")
        println("="^80)
        println(@sprintf("%-12s %-20s %-20s %-20s %-20s", 
                        "参数组", "Π_P 实部", "Π_P 虚部", "Π_S 实部", "Π_S 虚部"))
        println("-"^80)
        
        for case in test_cases
            k0 = case.λ - (case.μ1 - case.μ2)
            A1 = A(case.m1, case.μ1, case.T, Φ, Φbar, nodes, weights)
            A2 = A(case.m2, case.μ2, case.T, Φ, Φbar, nodes, weights)
            
            ΠP = polarization_aniso(:P, k0, case.k, case.m1, case.m2,
                case.μ1, case.μ2, case.T, Φ, Φbar, ξ, A1, A2, 0)
            ΠS = polarization_aniso(:S, k0, case.k, case.m1, case.m2,
                case.μ1, case.μ2, case.T, Φ, Φbar, ξ, A1, A2, 0)
            
            println(@sprintf("%-12s %-20.10e %-20.10e %-20.10e %-20.10e", 
                            case.name, ΠP[1], ΠP[2], ΠS[1], ΠS[2]))
            
            @test all(isfinite, ΠP)
            @test all(isfinite, ΠS)
        end
        println("="^80 * "\n")
    end

    @testset "A值在极化函数中的权重分析" begin
        # 分析A函数项相对于B0项的贡献
        test_params = B0_CORRECTION_TEST_PARAMS
        k0 = test_params.λ - (test_params.μ1 - test_params.μ2)
        
        nodes, weights = gauleg(0.0, P_MAX, GAUSS_POINTS)
        A1 = A(test_params.m1, test_params.μ1, test_params.T, test_params.Φ, test_params.Φbar, nodes, weights)
        A2 = A(test_params.m2, test_params.μ2, test_params.T, test_params.Φ, test_params.Φbar, nodes, weights)
        
        println("\n" * "="^90)
        println("A值在极化函数中的权重分析 (解释为何 ξ 对实部影响较小)")
        println("="^90)
        println("参数设置：k0 = $(k0), k = $(test_params.k), m1 = $(test_params.m1), m2 = $(test_params.m2)")
        println("         μ1 = $(test_params.μ1), μ2 = $(test_params.μ2), T = $(test_params.T)")
        println("-"^90)
        
        # 计算各个组成部分
        λ = test_params.λ
        
        # 计算 A 项贡献
        A_sum = A1 + A2
        factor = -N_color / (8 * pi^2)
        A_contribution = factor * A_sum
        
        println("\n1. A 函数项的贡献：")
        println(@sprintf("   A1 = %.10e", A1))
        println(@sprintf("   A2 = %.10e", A2))
        println(@sprintf("   A1 + A2 = %.10e", A_sum))
        println(@sprintf("   factor × (A1 + A2) = %.10e", A_contribution))
        
        # 测试不同 ξ 值下的 B0 及其修正项
        ξ_values = [0.0, 0.1, 0.2, 0.3, 0.5]
        
        println("\n2. B0 及其修正项随 ξ 的变化：")
        println(@sprintf("%-8s %-18s %-18s %-18s %-18s", "ξ", "B0 实部", "B0_corr 实部", "B0_total 实部", "修正/零阶比值(%)"))
        println("-"^90)
        
        for ξ in ξ_values
            # 计算 B0（零阶项）
            B0_result = B0(λ, test_params.k, test_params.m1, test_params.μ1, 
                          test_params.m2, test_params.μ2, test_params.T; 
                          Φ=test_params.Φ, Φbar=test_params.Φbar)
            
            # 计算 B0_correction（一阶修正）
            if ξ > EPS_SEGMENT
                B0_corr = B0_correction(λ, test_params.k, test_params.m1, test_params.m2,
                                       test_params.μ1, test_params.μ2, test_params.T,
                                       test_params.Φ, test_params.Φbar, ξ)
            else
                B0_corr = (0.0, 0.0)
            end
            
            B0_total_real = B0_result[1] + B0_corr[1]
            correction_ratio = abs(B0_result[1]) > 1e-15 ? abs(B0_corr[1] / B0_result[1]) * 100 : 0.0
            
            println(@sprintf("%-8.2f %-18.10e %-18.10e %-18.10e %-18.4f", 
                            ξ, B0_result[1], B0_corr[1], B0_total_real, correction_ratio))
        end
        
        # 计算 prefactor × B0 项的贡献
        println("\n3. B0 项在 Polarization 中的实际贡献：")
        prefactor_P = test_params.k^2 - λ^2 + (test_params.m1 - test_params.m2)^2
        prefactor_S = test_params.k^2 - λ^2 + (test_params.m1 + test_params.m2)^2
        
        println(@sprintf("   赝标量通道 prefactor_P = %.10e", prefactor_P))
        println(@sprintf("   标量通道 prefactor_S = %.10e", prefactor_S))
        
        println("\n4. 各项对 Polarization 实部的权重分析：")
        println(@sprintf("%-8s %-20s %-20s %-20s %-20s", "ξ", "A项贡献", "B0项贡献(P)", "总贡献(P)", "A项占比(%)"))
        println("-"^90)
        
        for ξ in ξ_values
            B0_result = B0(λ, test_params.k, test_params.m1, test_params.μ1, 
                          test_params.m2, test_params.μ2, test_params.T; 
                          Φ=test_params.Φ, Φbar=test_params.Φbar)
            
            if ξ > EPS_SEGMENT
                B0_corr = B0_correction(λ, test_params.k, test_params.m1, test_params.m2,
                                       test_params.μ1, test_params.μ2, test_params.T,
                                       test_params.Φ, test_params.Φbar, ξ)
            else
                B0_corr = (0.0, 0.0)
            end
            
            B0_total_real = B0_result[1] + B0_corr[1]
            B0_contribution_P = factor * prefactor_P * B0_total_real
            total_contribution_P = A_contribution + B0_contribution_P
            A_weight_percent = abs(A_contribution / total_contribution_P) * 100
            
            println(@sprintf("%-8.2f %-20.10e %-20.10e %-20.10e %-20.2f", 
                            ξ, A_contribution, B0_contribution_P, total_contribution_P, A_weight_percent))
        end
        
        println("\n5. 标量通道的分析：")
        println(@sprintf("%-8s %-20s %-20s %-20s %-20s", "ξ", "A项贡献", "B0项贡献(S)", "总贡献(S)", "A项占比(%)"))
        println("-"^90)
        
        for ξ in ξ_values
            B0_result = B0(λ, test_params.k, test_params.m1, test_params.μ1, 
                          test_params.m2, test_params.μ2, test_params.T; 
                          Φ=test_params.Φ, Φbar=test_params.Φbar)
            
            if ξ > EPS_SEGMENT
                B0_corr = B0_correction(λ, test_params.k, test_params.m1, test_params.m2,
                                       test_params.μ1, test_params.μ2, test_params.T,
                                       test_params.Φ, test_params.Φbar, ξ)
            else
                B0_corr = (0.0, 0.0)
            end
            
            B0_total_real = B0_result[1] + B0_corr[1]
            B0_contribution_S = factor * prefactor_S * B0_total_real
            total_contribution_S = A_contribution + B0_contribution_S
            A_weight_percent = abs(A_contribution / total_contribution_S) * 100
            
            println(@sprintf("%-8.2f %-20.10e %-20.10e %-20.10e %-20.2f", 
                            ξ, A_contribution, B0_contribution_S, total_contribution_S, A_weight_percent))
        end
        
        # 结论
        println("\n" * "="^90)
        println("关键发现：")
        B0_result_ref = B0(λ, test_params.k, test_params.m1, test_params.μ1, 
                          test_params.m2, test_params.μ2, test_params.T; 
                          Φ=test_params.Φ, Φbar=test_params.Φbar)
        B0_contribution_P_ref = factor * prefactor_P * B0_result_ref[1]
        B0_contribution_S_ref = factor * prefactor_S * B0_result_ref[1]
        
        println(@sprintf("1. A项的绝对贡献：%.6e", A_contribution))
        println(@sprintf("2. B0项的零阶贡献(P通道)：%.6e (A项的 %.2f%%)", 
                        B0_contribution_P_ref, abs(B0_contribution_P_ref/A_contribution)*100))
        println(@sprintf("3. B0项的零阶贡献(S通道)：%.6e (A项的 %.2f%%)", 
                        B0_contribution_S_ref, abs(B0_contribution_S_ref/A_contribution)*100))
        println("\n结论：A项占据了极化函数的主导地位，因此 ξ 对 B0 的修正")
        println("      只能产生相对较小的整体变化（约0.01%量级）。")
        println("="^90 * "\n")
        
        # 验证测试
        @test abs(A_contribution) > 0
        @test abs(B0_contribution_P_ref) > 0
        @test abs(B0_contribution_S_ref) > 0
    end
end

println("\n" * "="^70)
println("PolarizationAniso 模块测试完成！")
println("="^70)


