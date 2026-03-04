# 分析被积函数在 E→m (左端点) 的奇点行为
# 以及这对有根情况下积分策略的影响

using Printf
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegralsAniso.jl"))
include(joinpath(@__DIR__, "../../../src/relaxtime/OneLoopIntegrals.jl"))
using QuadGK: quadgk

println("=" ^ 70)
println("分析被积函数在 E→m 的奇点行为")
println("=" ^ 70)

# 测试参数
ξ = -0.2; T = 0.15; μ = 0.0; Φ = 0.0; Φbar = 0.0

# 有根情况
tc = (name="有根1", λ=-1.0, k=0.01, m=0.3, m_prime=0.3)

Emin = tc.m
Emax = OneLoopIntegrals.energy_cutoff(tc.m)

println("\n1. 被积函数各部分在 E→m 的行为:")
println("-" ^ 70)

# 采样接近左端点
Es = [Emin + 1e-10, Emin + 1e-8, Emin + 1e-6, Emin + 1e-4, Emin + 1e-2, Emin + 0.1]

println(@sprintf("%-12s %12s %12s %12s %12s %12s", "E-m", "p", "coeff_x", "1/coeff_x²", "p*1/coeff_x²", "integrand"))
println("-" ^ 80)

for E in Es
    p = sqrt(E^2 - tc.m^2)
    coeff_x = 2 * p * tc.k
    denominator_const = tc.λ^2 + 2*tc.λ*E + tc.m^2 - tc.m_prime^2 - tc.k^2
    
    # real_integral_tool 的主要项
    term1 = -2*denominator_const/coeff_x^2
    
    # 完整被积函数
    integrand_val = OneLoopIntegralsCorrection.real_integrand_k_positive(
        :quark, tc.λ, tc.k, tc.m, tc.m_prime, E, ξ, T, μ, Φ, Φbar)
    
    println(@sprintf("%.2e %12.4e %12.4e %12.4e %12.4e %12.4e", 
        E - tc.m, p, coeff_x, 1/coeff_x^2, p/coeff_x^2, integrand_val))
end

println("\n2. 分析 p * (1/coeff_x²) 的极限行为:")
println("-" ^ 70)
println("coeff_x = 2*p*k, 所以 1/coeff_x² = 1/(4*p²*k²)")
println("p * 1/coeff_x² = p / (4*p²*k²) = 1/(4*p*k²)")
println("当 p→0 时，这一项 → ∞")
println("\n但是 p = sqrt(E²-m²) ≈ sqrt(2m*(E-m)) 当 E→m")
println("所以 1/p ≈ 1/sqrt(2m*(E-m)) ~ (E-m)^(-1/2)")
println("这是一个可积的奇点！(积分 ∫(E-m)^(-1/2) dE 收敛)")

println("\n3. 检查有根情况下的根位置:")
println("-" ^ 70)

roots = OneLoopIntegralsCorrection.find_roots_AB(tc.λ, tc.k, tc.m, tc.m_prime, Emin, Emax)
println("找到的根: ", roots)
println("积分区间: [", Emin, ", ", Emax, "]")

if !isempty(roots)
    println("\n根与左端点的距离:")
    for (i, r) in enumerate(roots)
        println(@sprintf("  根 %d: E=%.6f, 距离左端点=%.6f", i, r, r - Emin))
    end
end

println("\n4. 有根情况下不同策略的精度对比:")
println("-" ^ 70)

# 高精度参考值
integrand(E) = OneLoopIntegralsCorrection.real_integrand_k_positive(
    :quark, tc.λ, tc.k, tc.m, tc.m_prime, E, ξ, T, μ, Φ, Φbar)
ref, _ = quadgk(integrand, Emin, Emax; rtol=1e-12)
println("参考值 (QuadGK rtol=1e-12): ", ref)

# 当前 HYBRID 策略
result_hybrid = OneLoopIntegralsCorrection.tilde_B0_correction_k_positive(
    :quark, tc.λ, tc.k, tc.m, tc.m_prime, μ, T, Φ, Φbar, ξ;
    strategy=OneLoopIntegralsCorrection.STRATEGY_HYBRID, cluster_n=32, diagnostics=true)
println("\nHYBRID (n=32):")
println("  结果: ", result_hybrid[1])
println("  相对误差: ", abs((result_hybrid[1] - ref) / ref))
println("  区间数: ", result_hybrid[3].n_intervals)
println("  区间: ", result_hybrid[3].intervals)

println("\n5. 测试改进策略：考虑左端点奇点")
println("-" ^ 70)

# 手动测试：对第一个区间使用 SING_BOTH (考虑左端点的 p→0 奇点)
intervals = OneLoopIntegralsCorrection.build_intervals_from_roots(roots, Emin, Emax)
println("区间分割: ", intervals)

# 对比不同节点选择
n = 32
total_current = 0.0
total_improved = 0.0

for (idx, (a, b)) in enumerate(intervals)
    global total_current, total_improved
    # 当前策略
    if idx == 1
        current_sing = OneLoopIntegralsCorrection.SING_RIGHT
        # 改进：第一个区间左端是 E=m，有 p→0 奇点，应该用 SING_BOTH
        improved_sing = OneLoopIntegralsCorrection.SING_BOTH
    elseif idx == length(intervals)
        current_sing = OneLoopIntegralsCorrection.SING_LEFT
        improved_sing = OneLoopIntegralsCorrection.SING_LEFT
    else
        current_sing = OneLoopIntegralsCorrection.SING_BOTH
        improved_sing = OneLoopIntegralsCorrection.SING_BOTH
    end
    
    xs_cur, wx_cur = OneLoopIntegralsCorrection.hybrid_nodes(a, b, n, current_sing)
    xs_imp, wx_imp = OneLoopIntegralsCorrection.hybrid_nodes(a, b, n, improved_sing)
    
    int_cur = sum(wx_cur .* integrand.(xs_cur))
    int_imp = sum(wx_imp .* integrand.(xs_imp))
    
    total_current += int_cur
    total_improved += int_imp
    
    println(@sprintf("区间 %d [%.4f, %.4f]:", idx, a, b))
    println(@sprintf("  当前策略 (%s): %.10f", current_sing, int_cur))
    println(@sprintf("  改进策略 (%s): %.10f", improved_sing, int_imp))
end

println("\n总结:")
println(@sprintf("  当前策略总和: %.10f, 相对误差: %.2e", total_current, abs((total_current - ref) / ref)))
println(@sprintf("  改进策略总和: %.10f, 相对误差: %.2e", total_improved, abs((total_improved - ref) / ref)))

println("\n6. 更多有根测试用例:")
println("-" ^ 70)

test_cases = [
    (name="有根1", λ=-1.0, k=0.01, m=0.3, m_prime=0.3),
    (name="有根2", λ=-1.0, k=0.1, m=0.3, m_prime=0.3),
    (name="有根3", λ=-0.8, k=0.02, m=0.25, m_prime=0.25),
]

println(@sprintf("%-8s %12s %12s %12s", "Case", "当前误差", "改进误差", "改进?"))
println("-" ^ 50)

for tc in test_cases
    Emin = tc.m
    Emax = OneLoopIntegrals.energy_cutoff(tc.m)
    integrand(E) = OneLoopIntegralsCorrection.real_integrand_k_positive(
        :quark, tc.λ, tc.k, tc.m, tc.m_prime, E, ξ, T, μ, Φ, Φbar)
    ref, _ = quadgk(integrand, Emin, Emax; rtol=1e-12)
    
    roots = OneLoopIntegralsCorrection.find_roots_AB(tc.λ, tc.k, tc.m, tc.m_prime, Emin, Emax)
    intervals = OneLoopIntegralsCorrection.build_intervals_from_roots(roots, Emin, Emax)
    
    total_current = 0.0
    total_improved = 0.0
    
    for (idx, (a, b)) in enumerate(intervals)
        local current_sing, improved_sing
        if idx == 1
            current_sing = isempty(roots) ? OneLoopIntegralsCorrection.SING_NONE : 
                          (length(intervals) == 1 ? OneLoopIntegralsCorrection.SING_BOTH : OneLoopIntegralsCorrection.SING_RIGHT)
            # 改进：第一个区间始终考虑左端点的 p→0 奇点
            improved_sing = OneLoopIntegralsCorrection.SING_BOTH
        elseif idx == length(intervals)
            current_sing = OneLoopIntegralsCorrection.SING_LEFT
            improved_sing = OneLoopIntegralsCorrection.SING_BOTH  # 也考虑左端点
        else
            current_sing = OneLoopIntegralsCorrection.SING_BOTH
            improved_sing = OneLoopIntegralsCorrection.SING_BOTH
        end
        
        xs_cur, wx_cur = OneLoopIntegralsCorrection.hybrid_nodes(a, b, n, current_sing)
        xs_imp, wx_imp = OneLoopIntegralsCorrection.hybrid_nodes(a, b, n, improved_sing)
        
        total_current += sum(wx_cur .* integrand.(xs_cur))
        total_improved += sum(wx_imp .* integrand.(xs_imp))
    end
    
    err_cur = abs((total_current - ref) / ref)
    err_imp = abs((total_improved - ref) / ref)
    better = err_imp < err_cur ? "是" : "否"
    
    println(@sprintf("%-8s %12.2e %12.2e %12s", tc.name, err_cur, err_imp, better))
end
