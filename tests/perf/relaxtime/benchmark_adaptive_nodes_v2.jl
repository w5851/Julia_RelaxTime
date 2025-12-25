"""
自适应节点数性能基准测试 v2

使用 BenchmarkTools 进行更精确的计时
"""

using Printf
using BenchmarkTools

# 加载模块
include("../../../src/relaxtime/OneLoopIntegralsAniso.jl")
using .OneLoopIntegralsCorrection: tilde_B0_correction_k_positive, 
    STRATEGY_HYBRID, STRATEGY_QUADGK

# 测试参数集
const TEST_CASES = [
    (0.5, 0.3, 0.3, 0.3, 0.3, 0.15, 0.3, 0.3, 0.1, "标准参数"),
    (0.8, 0.5, 0.3, 0.3, 0.3, 0.15, 0.3, 0.3, 0.1, "较大λ和k"),
    (0.3, 0.1, 0.3, 0.3, 0.3, 0.15, 0.3, 0.3, 0.1, "较小k"),
    (0.5, 0.3, 0.5, 0.3, 0.3, 0.15, 0.3, 0.3, 0.1, "不等质量"),
    (1.0, 0.8, 0.3, 0.3, 0.3, 0.15, 0.3, 0.3, 0.1, "大λ大k"),
    (0.2, 0.15, 0.3, 0.3, 0.3, 0.15, 0.3, 0.3, 0.1, "小参数"),
    (0.6, 0.4, 0.3, 0.5, 0.3, 0.15, 0.3, 0.3, 0.1, "m' > m"),
    (0.4, 0.2, 0.3, 0.3, 0.3, 0.20, 0.3, 0.3, 0.1, "较高温度"),
]

# 有根的测试用例（需要特殊参数）
const TEST_CASES_WITH_ROOTS = [
    # 这些参数会产生 A±B=0 的根
    (0.5, 0.5, 0.3, 0.3, 0.3, 0.15, 0.3, 0.3, 0.1, "k=λ (可能有根)"),
    (0.3, 0.3, 0.3, 0.3, 0.3, 0.15, 0.3, 0.3, 0.1, "小k=λ"),
]

function main()
    println("=" ^ 70)
    println("自适应节点数性能基准测试 (BenchmarkTools)")
    println("=" ^ 70)
    println()
    
    # 获取诊断信息
    println("节点分配情况:")
    println("-" ^ 70)
    
    for (i, (λ, k, m, m_prime, μ, T, Φ, Φbar, ξ, desc)) in enumerate(TEST_CASES)
        _, _, diag = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
            strategy=STRATEGY_HYBRID, diagnostics=true)
        
        n_roots = diag.n_roots
        n_intervals = diag.n_intervals
        
        if n_intervals == 1
            total_nodes = n_roots == 0 ? 16 : 32
        else
            middle_intervals = max(0, n_intervals - 2)
            total_nodes = 2 * 16 + middle_intervals * 32
        end
        
        @printf("Case %d: %-15s | 根数=%d, 区间数=%d, 总节点≈%d\n",
            i, desc, n_roots, n_intervals, total_nodes)
    end
    
    println()
    println("性能测试 (使用 @btime):")
    println("-" ^ 70)
    
    # 选取代表性用例进行详细测试
    λ, k, m, m_prime, μ, T, Φ, Φbar, ξ = TEST_CASES[1][1:9]
    
    println("\n标准参数 - HYBRID 策略:")
    @btime tilde_B0_correction_k_positive(:quark, $λ, $k, $m, $m_prime, $μ, $T, $Φ, $Φbar, $ξ;
        strategy=STRATEGY_HYBRID)
    
    println("\n标准参数 - QuadGK 策略:")
    @btime tilde_B0_correction_k_positive(:quark, $λ, $k, $m, $m_prime, $μ, $T, $Φ, $Φbar, $ξ;
        strategy=STRATEGY_QUADGK)
    
    # 测试有根情况
    println("\n" * "=" ^ 70)
    println("有根情况测试:")
    println("-" ^ 70)
    
    for (λ, k, m, m_prime, μ, T, Φ, Φbar, ξ, desc) in TEST_CASES_WITH_ROOTS
        _, _, diag = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
            strategy=STRATEGY_HYBRID, diagnostics=true)
        
        @printf("\n%s: 根数=%d, 区间数=%d\n", desc, diag.n_roots, diag.n_intervals)
        
        if diag.n_roots > 0
            println("HYBRID:")
            @btime tilde_B0_correction_k_positive(:quark, $λ, $k, $m, $m_prime, $μ, $T, $Φ, $Φbar, $ξ;
                strategy=STRATEGY_HYBRID)
        end
    end
    
    println()
    println("精度验证:")
    println("-" ^ 70)
    
    for (i, (λ, k, m, m_prime, μ, T, Φ, Φbar, ξ, desc)) in enumerate(TEST_CASES)
        real_hybrid, _ = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
            strategy=STRATEGY_HYBRID)
        real_quadgk, _ = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
            strategy=STRATEGY_QUADGK, rtol=1e-10, atol=1e-12)
        
        rel_err = abs(real_hybrid - real_quadgk) / max(abs(real_quadgk), 1e-15)
        @printf("Case %d: %-15s | 相对误差: %.2e %s\n",
            i, desc, rel_err, rel_err < 1e-4 ? "✓" : "✗")
    end
    
    println()
    println("=" ^ 70)
end

main()
