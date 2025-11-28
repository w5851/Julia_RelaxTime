# 分析ss_to_ss散射过程中k0²-u<0的问题

include("../../src/Constants_PNJL.jl")
include("../../src/relaxtime/ScatteringAmplitude.jl")

using .Constants_PNJL
using .ScatteringAmplitude
using .ScatteringAmplitude.TotalPropagator

println("="^70)
println("ss→ss散射过程中k0²-u<0问题分析")
println("="^70)

# 设置参数（与测试相同）
T = 0.15  # fm⁻¹
Φ = 0.5
Φbar = 0.5
ξ = 1.0

m_u = 0.3  # fm⁻¹
m_s = 0.5  # fm⁻¹
μ_u = 0.25  # fm⁻¹
μ_s = 0.25  # fm⁻¹

# 简化的A值（测试用）
A_u = 0.05  # fm⁻²
A_s = 0.08  # fm⁻²

quark_params = (
    m = (u=m_u, d=m_u, s=m_s),
    μ = (u=μ_u, d=μ_u, s=μ_s),
    A = (u=A_u, d=A_u, s=A_s)
)

thermo_params = (T=T, Φ=Φ, Φbar=Φbar, ξ=ξ)

# 测试中使用的参数
s = 20.0  # fm⁻²
t = -0.2  # fm⁻²
process = :ss_to_ss

println("\n## 1. 输入参数")
println("-"^70)
println("散射过程: $process (s + s → s + s)")
println("s夸克质量: m_s = $m_s fm⁻¹")
println("Mandelstam s: $s fm⁻²")
println("Mandelstam t: $t fm⁻²")

# 提取质量
m1, m2, m3, m4 = get_quark_masses_for_process(process, quark_params)
println("\n粒子质量:")
println("  m₁ = $m1 fm⁻¹ (入射s夸克1)")
println("  m₂ = $m2 fm⁻¹ (入射s夸克2)")
println("  m₃ = $m3 fm⁻¹ (出射s夸克3)")
println("  m₄ = $m4 fm⁻¹ (出射s夸克4)")

# 计算u
u = m1^2 + m2^2 + m3^2 + m4^2 - s - t
println("\n## 2. Mandelstam变量约束")
println("-"^70)
println("约束: s + t + u = m₁² + m₂² + m₃² + m₄²")
println("      s + t + u = $(m1^2) + $(m2^2) + $(m3^2) + $(m4^2) = $(m1^2 + m2^2 + m3^2 + m4^2) fm⁻²")
println("\n计算u:")
println("  u = Σm² - s - t")
println("  u = $(m1^2 + m2^2 + m3^2 + m4^2) - $s - ($t)")
println("  u = $(u) fm⁻²")

if u > 0
    println("\n⚠️  注意: u = $(u) > 0")
    println("    对于qq散射，通常期望 u < 0 (类空动量转移)")
else
    println("\n✓ u = $(u) < 0 (正常)")
end

# 分析质心系能量
println("\n## 3. 质心系能量分析")
println("-"^70)

# 计算质心系能量
E_cm = sqrt(s)
println("质心系总能量: E_cm = √s = $(E_cm) fm⁻¹")

# 计算阈值能量
E_threshold = m1 + m2
E_threshold_squared = (m1 + m2)^2
println("阈值能量: E_th = m₁ + m₂ = $(E_threshold) fm⁻¹")
println("阈值能量平方: s_th = E_th² = $(E_threshold_squared) fm⁻²")

if s > E_threshold_squared
    println("✓ s = $s > s_th = $(E_threshold_squared) (能量充足)")
    excess_energy = E_cm - E_threshold
    println("  过阈能量: ΔE = $(excess_energy) fm⁻¹")
else
    println("✗ s = $s < s_th = $(E_threshold_squared) (低于阈值！)")
end

# u道质心系动量计算
println("\n## 4. u道质心系动量分析")
println("-"^70)
println("u道交换对应: p₁ - p₄ (粒子1和4交换介子)")
println("u道能量: k0 = |E₁ - E₄|")
println("u道动量: k = √(k0² - u)")

# 计算单粒子能量
E1 = sqrt(m1^2 + (E_cm/2)^2)  # 粗略估计
E4 = sqrt(m4^2 + (E_cm/2)^2)
println("\n粗略估计的质心系单粒子能量:")
println("  E₁ ≈ √(m₁² + p²) ≈ $(E1) fm⁻¹")
println("  E₄ ≈ √(m₄² + p²) ≈ $(E4) fm⁻¹")

k0_estimate = abs(E1 - E4)
delta_estimate = k0_estimate^2 - u
println("\n粗略估计:")
println("  k0 ≈ |E₁ - E₄| ≈ $(k0_estimate) fm⁻¹")
println("  k0² ≈ $(k0_estimate^2) fm⁻²")
println("  k0² - u ≈ $(k0_estimate^2) - $(u) = $(delta_estimate) fm⁻²")

if delta_estimate < 0
    println("\n⚠️  k0² - u < 0: 无法计算实数动量!")
    println("    这表示u道交换在这个运动学配置下不可达")
end

# 实际调用看看会发生什么
println("\n## 5. 实际调用calculate_cms_momentum")
println("-"^70)

try
    cms_u = calculate_cms_momentum(process, s, t, :u, quark_params; u=u)
    println("✓ 成功计算u道CMS动量:")
    println("  k0 = $(cms_u.k0) fm⁻¹")
    println("  k = $(cms_u.k) fm⁻¹")
catch e
    println("✗ 调用失败:")
    println("  错误: ", sprint(showerror, e))
end

# 分析根本原因
println("\n## 6. 根本原因分析")
println("-"^70)
println("问题根源:")
println("1. **u值为正**: u = $(u) > 0")
println("   - 对于ss→ss，Σm² = 4×$(m_s^2) = $(4*m_s^2) fm⁻²")
println("   - s = $s fm⁻² 很大（远超阈值）")
println("   - t = $t fm⁻² 很小（小角散射）")
println("   - 导致 u = $(4*m_s^2) - $s - ($t) = $(u) > 0")
println()
println("2. **物理意义**:")
println("   - qq散射中，通常 t < 0, u < 0 (两者都是类空动量转移)")
println("   - u > 0 意味着u道是类时的，对应能量转移而非动量转移")
println("   - 这在弹性散射中通常不会发生")
println()
println("3. **k0² - u < 0 的含义**:")
println("   - k0 来自于质心系能量差 |E₁ - E₄|")
println("   - 对于同种粒子弹性散射，|E₁ - E₄| 通常很小")
println("   - 当 u > k0² 时，√(k0² - u) 为虚数")
println("   - 这表示u道交换在这个运动学条件下是**禁戒的**")

# 参数建议
println("\n## 7. 参数选择建议")
println("-"^70)
println("要避免这个问题，需要确保 u < 0:")
println()

# 计算需要的t值
u_target = -0.5  # 目标u值（负数）
t_needed = 4*m_s^2 - s - u_target
println("方案1: 固定s=$s fm⁻², 增大|t|")
println("  要使 u = $u_target fm⁻²")
println("  需要 t = $(t_needed) fm⁻²")
println("  (即 |t| = $(abs(t_needed)) fm⁻²)")
println()

# 计算需要的s值
s_needed = 4*m_s^2 - t - u_target
println("方案2: 固定t=$t fm⁻², 减小s")
println("  要使 u = $u_target fm⁻²")
println("  需要 s = $(s_needed) fm⁻²")
println("  阈值 s_th = $(E_threshold_squared) fm⁻²")
if s_needed > E_threshold_squared
    println("  ✓ $s_needed > $(E_threshold_squared) (满足阈值条件)")
else
    println("  ✗ $s_needed < $(E_threshold_squared) (低于阈值)")
end
println()

println("方案3: 同时调整s和t")
println("  建议参数组合:")
for s_test in [5.0, 8.0, 10.0, 12.0]
    for t_test in [-0.5, -1.0, -1.5, -2.0]
        u_test = 4*m_s^2 - s_test - t_test
        if u_test < 0 && s_test > E_threshold_squared
            println("    s = $s_test fm⁻², t = $t_test fm⁻² → u = $(round(u_test, digits=2)) fm⁻² ✓")
            break
        end
    end
end

println("\n" * "="^70)
println("分析完成")
println("="^70)

