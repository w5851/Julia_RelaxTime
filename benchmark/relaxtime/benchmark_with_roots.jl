"""
有根情况的性能基准测试
"""

using Printf
using BenchmarkTools

include("../../../src/relaxtime/OneLoopIntegralsAniso.jl")
using .OneLoopIntegralsCorrection: tilde_B0_correction_k_positive, 
    STRATEGY_HYBRID, STRATEGY_QUADGK

# 有根的测试用例
const ROOT_CASES = [
    (0.1, 0.2, 0.3, 0.3, "1根-小参数"),
    (0.1, 0.3, 0.3, 0.5, "1根-不等质量"),
    (0.2, 0.3, 0.3, 0.3, "1根-中等参数"),
    (0.3, 0.5, 0.3, 0.3, "1根-较大k"),
]

function main()
    μ, T, Φ, Φbar, ξ = 0.3, 0.15, 0.3, 0.3, 0.1
    
    println("=" ^ 70)
    println("有根情况性能基准测试")
    println("=" ^ 70)
    println()
    
    println("诊断信息:")
    println("-" ^ 70)
    
    for (λ, k, m, m_prime, desc) in ROOT_CASES
        _, _, diag = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
            strategy=STRATEGY_HYBRID, diagnostics=true)
        
        # 计算总节点数
        n_intervals = diag.n_intervals
        if n_intervals == 1
            total_nodes = diag.n_roots == 0 ? 16 : 32
        else
            middle_intervals = max(0, n_intervals - 2)
            total_nodes = 2 * 16 + middle_intervals * 32
        end
        
        @printf("%-20s | 根数=%d, 区间数=%d, 总节点≈%d\n",
            desc, diag.n_roots, n_intervals, total_nodes)
        @printf("                     | 根位置: %s\n", diag.roots)
    end
    
    println()
    println("性能测试:")
    println("-" ^ 70)
    
    for (λ, k, m, m_prime, desc) in ROOT_CASES
        println("\n$desc:")
        
        print("  HYBRID: ")
        @btime tilde_B0_correction_k_positive(:quark, $λ, $k, $m, $m_prime, $μ, $T, $Φ, $Φbar, $ξ;
            strategy=STRATEGY_HYBRID)
        
        print("  QuadGK: ")
        @btime tilde_B0_correction_k_positive(:quark, $λ, $k, $m, $m_prime, $μ, $T, $Φ, $Φbar, $ξ;
            strategy=STRATEGY_QUADGK)
    end
    
    println()
    println("精度验证:")
    println("-" ^ 70)
    
    for (λ, k, m, m_prime, desc) in ROOT_CASES
        real_hybrid, _ = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
            strategy=STRATEGY_HYBRID)
        real_quadgk, _ = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
            strategy=STRATEGY_QUADGK, rtol=1e-10, atol=1e-12)
        
        rel_err = abs(real_hybrid - real_quadgk) / max(abs(real_quadgk), 1e-15)
        @printf("%-20s | 相对误差: %.2e %s\n",
            desc, rel_err, rel_err < 1e-4 ? "✓" : "✗")
    end
    
    println()
    println("=" ^ 70)
end

main()
