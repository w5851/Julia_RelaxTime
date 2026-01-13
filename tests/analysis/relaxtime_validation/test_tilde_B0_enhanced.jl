"""
测试增强版 tilde_B0_correction_k_positive 函数

测试内容：
1. 不同积分策略的正确性
2. 聚簇 GL 参数敏感性
3. 边界情况处理
4. 诊断输出功能
5. 与参考值的精度对比

运行方式：
```julia
using Pkg
Pkg.activate(".")
include("tests/unit/test_tilde_B0_enhanced.jl")
```
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
using Test
using Printf

include("../../../src/relaxtime/OneLoopIntegralsAniso.jl")
include("../../../src/relaxtime/OneLoopIntegrals.jl")
using QuadGK: quadgk

using .OneLoopIntegrals: energy_cutoff

using .OneLoopIntegralsCorrection: tilde_B0_correction_k_positive
    IntegrationStrategy, STRATEGY_QUADGK, STRATEGY_INTERVAL_GL, STRATEGY_CLUSTER_GL,
    IntegrationDiagnostics, real_integrand_k_positive

# 标准测试参数
const TEST_PARAMS = (
    λ = -1.0,
    k = 0.01,
    m = 0.3,
    m_prime = 0.3,
    ξ = -0.2,
    T = 0.15,
    μ = 0.0,
    Φ = 0.0,
    Φbar = 0.0,
)

# 计算参考值
function compute_reference(; params=TEST_PARAMS)
    integrand(E) = real_integrand_k_positive(:quark, params.λ, params.k, params.m, params.m_prime, E, 
        params.ξ, params.T, params.μ, params.Φ, params.Φbar)
    Emin = params.m
    Emax = energy_cutoff(params.m)
    val, _ = quadgk(integrand, Emin, Emax; rtol=1e-12)
    return val
end

@testset "tilde_B0_correction_k_positive 增强版测试" begin
    
    @testset "基本功能测试" begin
        p = TEST_PARAMS
        result = tilde_B0_correction_k_positive(:quark, p.λ, p.k, p.m, p.m_prime, p.μ, p.T, p.Φ, p.Φbar, p.ξ)
        
        @test result isa Tuple{Float64, Float64}
        @test isfinite(result[1])
        @test isfinite(result[2])
    end
    
    @testset "策略选择测试" begin
        p = TEST_PARAMS
        ref = compute_reference()
        
        # STRATEGY_QUADGK
        r_quadgk = tilde_B0_correction_k_positive(:quark, p.λ, p.k, p.m, p.m_prime, p.μ, p.T, p.Φ, p.Φbar, p.ξ;
            strategy=STRATEGY_QUADGK)
        @test isfinite(r_quadgk[1])
        
        # STRATEGY_INTERVAL_GL
        r_interval = tilde_B0_correction_k_positive(:quark, p.λ, p.k, p.m, p.m_prime, p.μ, p.T, p.Φ, p.Φbar, p.ξ;
            strategy=STRATEGY_INTERVAL_GL, cluster_n=64)
        @test isfinite(r_interval[1])
        
        # STRATEGY_CLUSTER_GL (默认，应该最精确)
        r_cluster = tilde_B0_correction_k_positive(:quark, p.λ, p.k, p.m, p.m_prime, p.μ, p.T, p.Φ, p.Φbar, p.ξ;
            strategy=STRATEGY_CLUSTER_GL, cluster_beta=8.0, cluster_n=64)
        @test isfinite(r_cluster[1])
        
        # 聚簇 GL 应该比其他方法更精确
        err_cluster = abs((r_cluster[1] - ref) / ref)
        @test err_cluster < 1e-4  # 相对误差 < 0.01%
    end
    
    @testset "聚簇参数敏感性" begin
        p = TEST_PARAMS
        ref = compute_reference()
        
        # 测试不同 beta 值
        for beta in [2.0, 4.0, 8.0]
            result = tilde_B0_correction_k_positive(:quark, p.λ, p.k, p.m, p.m_prime, p.μ, p.T, p.Φ, p.Φbar, p.ξ;
                strategy=STRATEGY_CLUSTER_GL, cluster_beta=beta, cluster_n=64)
            @test isfinite(result[1])
        end
        
        # beta=8.0 应该给出最佳精度
        r_beta8 = tilde_B0_correction_k_positive(:quark, p.λ, p.k, p.m, p.m_prime, p.μ, p.T, p.Φ, p.Φbar, p.ξ;
            strategy=STRATEGY_CLUSTER_GL, cluster_beta=8.0, cluster_n=64)
        err_beta8 = abs((r_beta8[1] - ref) / ref)
        @test err_beta8 < 1e-4
    end
    
    @testset "节点数收敛性" begin
        p = TEST_PARAMS
        ref = compute_reference()
        
        errors = Float64[]
        for n in [16, 32, 64, 128]
            result = tilde_B0_correction_k_positive(:quark, p.λ, p.k, p.m, p.m_prime, p.μ, p.T, p.Φ, p.Φbar, p.ξ;
                strategy=STRATEGY_CLUSTER_GL, cluster_beta=8.0, cluster_n=n)
            push!(errors, abs((result[1] - ref) / ref))
        end
        
        # 误差应该随节点数增加而减小（大致趋势）
        @test errors[end] < errors[1]  # n=128 比 n=16 更精确
    end
    
    @testset "诊断输出测试" begin
        p = TEST_PARAMS
        
        real_part, imag_part, diag = tilde_B0_correction_k_positive(:quark, p.λ, p.k, p.m, p.m_prime, p.μ, p.T, p.Φ, p.Φbar, p.ξ;
            strategy=STRATEGY_CLUSTER_GL, diagnostics=true)
        
        @test diag isa IntegrationDiagnostics
        @test diag.strategy == STRATEGY_CLUSTER_GL
        @test diag.n_roots >= 0
        @test diag.n_intervals >= 1
        @test length(diag.roots) == diag.n_roots
        @test length(diag.intervals) == diag.n_intervals
        @test diag.real_part == real_part
        @test diag.imag_part == imag_part
        @test diag.elapsed_ms >= 0.0
    end
    
    @testset "不同参数组合" begin
        test_cases = [
            (λ=-1.0, k=0.01, m=0.3, m_prime=0.3),
            (λ=-1.0, k=0.1, m=0.3, m_prime=0.3),
            (λ=-0.5, k=0.05, m=0.3, m_prime=0.3),
            (λ=0.5, k=0.1, m=0.3, m_prime=0.35),
        ]
        
        for tc in test_cases
            result = tilde_B0_correction_k_positive(:quark, tc.λ, tc.k, tc.m, tc.m_prime, 
                0.0, 0.15, 0.0, 0.0, -0.2;
                strategy=STRATEGY_CLUSTER_GL)
            @test isfinite(result[1])
            @test isfinite(result[2])
        end
    end
    
    @testset "符号类型测试" begin
        p = TEST_PARAMS
        
        # :quark
        r_quark = tilde_B0_correction_k_positive(:quark, p.λ, p.k, p.m, p.m_prime, p.μ, p.T, p.Φ, p.Φbar, p.ξ)
        @test isfinite(r_quark[1])
        
        # :antiquark
        r_antiquark = tilde_B0_correction_k_positive(:antiquark, p.λ, p.k, p.m, p.m_prime, p.μ, p.T, p.Φ, p.Φbar, p.ξ)
        @test isfinite(r_antiquark[1])
    end
end

println("All tests passed!")

