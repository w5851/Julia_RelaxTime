"""
分析被积函数中 1/p 项的来源

问题：
1. 为什么只在各向异性中出现 1/p 项而各向同性下没有？
2. 各向异性中出现的 1/p 项导致的积分下限处的发散是否会被 p*coeff(...) 抵消？
"""

using Printf

# 加载模块
include("../../../src/relaxtime/OneLoopIntegralsAniso.jl")
include("../../../src/relaxtime/OneLoopIntegrals.jl")

using .OneLoopIntegralsCorrection
using .OneLoopIntegrals

println("=" ^ 70)
println("分析被积函数中 1/p 项的来源")
println("=" ^ 70)
println()

# 物理参数
m = 0.3
m_prime = 0.3
λ = 0.5
k = 0.3
μ = 0.3
T = 0.15
Φ = 0.3
Φbar = 0.3
ξ = 0.1

println("物理参数: m=$m, m'=$m_prime, λ=$λ, k=$k")
println()

println("=" ^ 70)
println("1. 各向同性 B0 的被积函数分析")
println("=" ^ 70)
println()

println("各向同性 B0 的被积函数 (k>0 情况):")
println("  f(E) = dist(E) * log|...| / k")
println()
println("其中 dist(E) 是分布函数，log|...| 是角度积分的结果。")
println("注意：这里没有显式的 1/p 项！")
println()

println("=" ^ 70)
println("2. 各向异性修正项的被积函数分析")
println("=" ^ 70)
println()

println("根据公式文档，各向异性修正项的被积函数为:")
println()
println("  k>0 情况:")
println("  f_aniso(E) = (1/4k³) * C(E)/p * (B² L - 4Bpk)")
println()
println("其中:")
println("  C(E) = ξβ * p²/(2E) * G(u)  -- 各向异性修正系数")
println("  B = λ² + 2λE - k² + m² - m'²")
println("  L = ln|(2pk + B)/(2pk - B)|")
println("  p = sqrt(E² - m²)")
println()

println("展开 C(E)/p:")
println("  C(E)/p = ξβ * p²/(2E) * G(u) / p = ξβ * p/(2E) * G(u)")
println()
println("所以被积函数中实际上是 p/E * G(u) * (B² L - 4Bpk)/(4k³)")
println("这里 p → 0 时，整个表达式 → 0（因为 p 在分子）")
println()

println("=" ^ 70)
println("3. 代码实现分析")
println("=" ^ 70)
println()

println("查看代码中的 real_integrand_k_positive 函数:")
println()
println("  real_ = 2.0*p * real_integral_tool(coeff_x, denominator_const) *")
println("          correction_cos_theta_coefficient(sign_, p, m, μ, T, Φ, Φbar, ξ)")
println()
println("其中:")
println("  coeff_x = 2pk  (A 系数)")
println("  denominator_const = λ² + 2λE + m² - m'² - k²  (B 系数)")
println()
println("  real_integral_tool(A, B) = -2B/A² + B²/A³ * ln|(A+B)/(A-B)|")
println("                          = -2B/(2pk)² + B²/(2pk)³ * ln|...|")
println("                          = -B/(2p²k²) + B²/(8p³k³) * ln|...|")
println()
println("所以:")
println("  2p * real_integral_tool = 2p * [-B/(2p²k²) + B²/(8p³k³) * ln|...|]")
println("                         = -B/(pk²) + B²/(4p²k³) * ln|...|")
println()
println("这里确实出现了 1/p 和 1/p² 项！")
println()

println("=" ^ 70)
println("4. 发散性分析")
println("=" ^ 70)
println()

println("当 E → m 时，p = sqrt(E² - m²) → 0")
println()
println("被积函数的行为:")
println("  2p * real_integral_tool * C(E)")
println("  = [-B/(pk²) + B²/(4p²k³) * ln|...|] * [ξβ * p²/(2E) * G(u)]")
println("  = -B * ξβ * p/(2Ek²) * G(u) + B² * ξβ/(8Ek³) * G(u) * ln|...|")
println()
println("第一项: ~ p → 0 (可积)")
println("第二项: ~ ln|...| 当 p → 0 时")
println()
println("关键问题：ln|(2pk + B)/(2pk - B)| 当 p → 0 时的行为")
println("  当 p → 0: ln|(B)/(−B)| = ln|−1| = ln(1) = 0 (如果 B ≠ 0)")
println("  但如果 B = 0（即存在根），则 ln|...| 发散")
println()

println("=" ^ 70)
println("5. 数值验证")
println("=" ^ 70)
println()

# 计算被积函数在 E → m 附近的值
Emin = m
Emax = OneLoopIntegrals.energy_cutoff(m)

println("E 接近 m 时被积函数的值:")
println("-" ^ 70)
@printf("%-15s | %-15s | %-15s | %-15s\n", "E - m", "p", "被积函数", "p * 被积函数/p")

for δE in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
    E = m + δE
    p = sqrt(E^2 - m^2)
    
    # 计算被积函数
    coeff_x = 2 * p * k
    denominator_const = λ^2 + 2*λ*E + m^2 - m_prime^2 - k^2
    
    # real_integral_tool 的值
    if abs(coeff_x) > 1e-15
        term1 = -2*denominator_const/coeff_x^2
        term2 = denominator_const^2/coeff_x^3
        log_arg = (coeff_x + denominator_const) / (coeff_x - denominator_const)
        term3 = log(abs(log_arg))
        real_int = term1 + term2 * term3
    else
        real_int = 0.0
    end
    
    # correction_cos_theta_coefficient 的近似值（简化）
    # 实际上 C(E) ~ ξ * p² / E，所以 C(E) → 0 当 p → 0
    C_approx = ξ * p^2 / E  # 简化近似
    
    integrand = 2 * p * real_int * C_approx
    
    @printf("%.2e       | %.6e | %.6e | %.6e\n", δE, p, integrand, integrand / p * p)
end

println()
println("=" ^ 70)
println("6. 结论")
println("=" ^ 70)
println()
println("1. 各向同性 B0 的被积函数不含显式 1/p 项，因为角度积分后")
println("   结果是 log|...| 形式，与 p 无关（除了 log 的参数）。")
println()
println("2. 各向异性修正项的被积函数确实包含 1/p 和 1/p² 项，")
println("   来源于 x² 积分 ∫x²/(Ax+B)dx 的结果。")
println()
println("3. 但是，各向异性修正系数 C(E) = ξβ * p²/(2E) * G(u)")
println("   包含 p² 因子，这会抵消 1/p² 项的发散。")
println()
println("4. 最终被积函数在 E → m (p → 0) 时的行为:")
println("   - 第一项 ~ p → 0 (可积)")
println("   - 第二项 ~ ln|...| → 有限值 (如果 B ≠ 0)")
println()
println("5. 因此，积分在下限 E = m 处是收敛的，不存在真正的发散。")
println("   但被积函数在 p → 0 附近变化剧烈，需要聚簇节点来提高精度。")
