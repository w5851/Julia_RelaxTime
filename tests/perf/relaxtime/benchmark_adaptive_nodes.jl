"""
自适应节点数性能基准测试

比较固定 32 节点 vs 自适应节点数 (16/32) 的性能差异
"""

using Printf
using QuadGK

# 加载模块
include("../../../src/relaxtime/OneLoopIntegralsAniso.jl")
using .OneLoopIntegralsCorrection: tilde_B0_correction_k_positive, 
    STRATEGY_HYBRID, STRATEGY_QUADGK, IntegrationDiagnostics

# 测试参数集
const TEST_CASES = [
    # (λ, k, m, m_prime, μ, T, Φ, Φbar, ξ, description)
    (0.5, 0.3, 0.3, 0.3, 0.3, 0.15, 0.3, 0.3, 0.1, "标准参数"),
    (0.8, 0.5, 0.3, 0.3, 0.3, 0.15, 0.3, 0.3, 0.1, "较大λ和k"),
    (0.3, 0.1, 0.3, 0.3, 0.3, 0.15, 0.3, 0.3, 0.1, "较小k"),
    (0.5, 0.3, 0.5, 0.3, 0.3, 0.15, 0.3, 0.3, 0.1, "不等质量"),
    (1.0, 0.8, 0.3, 0.3, 0.3, 0.15, 0.3, 0.3, 0.1, "大λ大k"),
    (0.2, 0.15, 0.3, 0.3, 0.3, 0.15, 0.3, 0.3, 0.1, "小参数"),
    (0.6, 0.4, 0.3, 0.5, 0.3, 0.15, 0.3, 0.3, 0.1, "m' > m"),
    (0.4, 0.2, 0.3, 0.3, 0.3, 0.20, 0.3, 0.3, 0.1, "较高温度"),
]

function benchmark_single(λ, k, m, m_prime, μ, T, Φ, Φbar, ξ; n_warmup=10, n_iter=100)
    sign_ = :quark
    
    # 预热
    for _ in 1:n_warmup
        tilde_B0_correction_k_positive(sign_, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
            strategy=STRATEGY_HYBRID)
    end
    
    # 计时
    times = Float64[]
    for _ in 1:n_iter
        t0 = time()
        tilde_B0_correction_k_positive(sign_, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
            strategy=STRATEGY_HYBRID)
        push!(times, (time() - t0) * 1000)  # ms
    end
    
    return mean(times), std(times)
end

function mean(x)
    return sum(x) / length(x)
end

function std(x)
    m = mean(x)
    return sqrt(sum((xi - m)^2 for xi in x) / (length(x) - 1))
end

function main()
    println("=" ^ 70)
    println("自适应节点数性能基准测试")
    println("=" ^ 70)
    println()
    
    # 获取诊断信息来查看节点分配
    println("节点分配情况:")
    println("-" ^ 70)
    
    for (i, (λ, k, m, m_prime, μ, T, Φ, Φbar, ξ, desc)) in enumerate(TEST_CASES)
        _, _, diag = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
            strategy=STRATEGY_HYBRID, diagnostics=true)
        
        n_roots = diag.n_roots
        n_intervals = diag.n_intervals
        
        # 计算总节点数（自适应）
        # 规则：中间区间（SING_BOTH）用 32 节点，其他用 16 节点
        if n_intervals == 1
            total_nodes = n_roots == 0 ? 16 : 32  # 无根用16，有根用32
        else
            # 第一个和最后一个区间用 16，中间区间用 32
            middle_intervals = max(0, n_intervals - 2)
            total_nodes = 2 * 16 + middle_intervals * 32
        end
        
        @printf("Case %d: %-15s | 根数=%d, 区间数=%d, 总节点≈%d\n",
            i, desc, n_roots, n_intervals, total_nodes)
    end
    
    println()
    println("性能测试:")
    println("-" ^ 70)
    @printf("%-20s | %12s | %12s\n", "测试用例", "平均时间(ms)", "标准差(ms)")
    println("-" ^ 70)
    
    total_time = 0.0
    for (i, (λ, k, m, m_prime, μ, T, Φ, Φbar, ξ, desc)) in enumerate(TEST_CASES)
        avg, sd = benchmark_single(λ, k, m, m_prime, μ, T, Φ, Φbar, ξ)
        total_time += avg
        @printf("%-20s | %12.4f | %12.4f\n", desc, avg, sd)
    end
    
    println("-" ^ 70)
    @printf("平均每次调用: %.4f ms\n", total_time / length(TEST_CASES))
    println()
    
    # 与 QuadGK 对比
    println("与 QuadGK 对比:")
    println("-" ^ 70)
    
    λ, k, m, m_prime, μ, T, Φ, Φbar, ξ = TEST_CASES[1][1:9]
    
    # HYBRID 时间
    _, _, diag_hybrid = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
        strategy=STRATEGY_HYBRID, diagnostics=true)
    
    # QuadGK 时间
    _, _, diag_quadgk = tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
        strategy=STRATEGY_QUADGK, diagnostics=true)
    
    # 多次测量取平均
    hybrid_times = Float64[]
    quadgk_times = Float64[]
    
    for _ in 1:100
        t0 = time()
        tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
            strategy=STRATEGY_HYBRID)
        push!(hybrid_times, (time() - t0) * 1000)
        
        t0 = time()
        tilde_B0_correction_k_positive(:quark, λ, k, m, m_prime, μ, T, Φ, Φbar, ξ;
            strategy=STRATEGY_QUADGK)
        push!(quadgk_times, (time() - t0) * 1000)
    end
    
    @printf("HYBRID:  %.4f ms (std: %.4f)\n", mean(hybrid_times), std(hybrid_times))
    @printf("QuadGK:  %.4f ms (std: %.4f)\n", mean(quadgk_times), std(quadgk_times))
    @printf("加速比:  %.2fx\n", mean(quadgk_times) / mean(hybrid_times))
    
    println()
    println("精度验证:")
    println("-" ^ 70)
    
    for (i, (λ, k, m, m_prime, μ, T, Φ, Φbar, ξ, desc)) in enumerate(TEST_CASES[1:4])
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
