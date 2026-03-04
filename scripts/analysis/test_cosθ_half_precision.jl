"""
测试半区间角度积分节点的精度优势

对于对称被积函数 f(cosθ) = f(-cosθ)，利用对称性可以提高精度：
∫₋₁¹ f(cosθ) d(cosθ) = 2 ∫₀¹ f(cosθ) d(cosθ)
"""

include("../../src/integration/GaussLegendre.jl")
using .GaussLegendre

println("=" ^ 70)
println("半区间角度积分节点精度测试")
println("=" ^ 70)

println("\n节点配置:")
println("  全区间 [-1, 1]: ", length(DEFAULT_COSΘ_NODES), " 个节点")
println("  半区间 [0, 1]:  ", length(DEFAULT_COSΘ_HALF_NODES), " 个节点")

# ============================================================================
# 测试1: 对称函数 - 半区间应该更精确
# ============================================================================
println("\n" * repeat("=", 70))
println("测试1: 对称函数 f(x) = x²")
println(repeat("=", 70))

# 定义对称函数
f_symmetric(x) = x^2

# 精确值: ∫₋₁¹ x² dx = 2/3
exact_value = 2.0 / 3.0

# 方法1: 全区间积分
result_full = sum(DEFAULT_COSΘ_WEIGHTS .* f_symmetric.(DEFAULT_COSΘ_NODES))
error_full = abs(result_full - exact_value) / exact_value * 100

# 方法2: 半区间积分（利用对称性）
result_half = 2.0 * sum(DEFAULT_COSΘ_HALF_WEIGHTS .* f_symmetric.(DEFAULT_COSΘ_HALF_NODES))
error_half = abs(result_half - exact_value) / exact_value * 100

println("\n精确值: ", exact_value)
println("\n全区间积分 (32节点):")
println("  数值结果: ", result_full)
println("  相对误差: ", error_full, " %")
println("\n半区间对称积分 (32节点):")
println("  数值结果: ", result_half)
println("  相对误差: ", error_half, " %")
println("\n精度提升: ", error_full / error_half, "x")

# ============================================================================
# 测试2: 更复杂的对称函数
# ============================================================================
println("\n" * repeat("=", 70))
println("测试2: 对称函数 f(x) = exp(-x²)")
println(repeat("=", 70))

f_gaussian(x) = exp(-x^2)

# 精确值需要数值计算（用高精度）
using QuadGK
exact_value_2, _ = quadgk(f_gaussian, -1.0, 1.0, rtol=1e-12)

result_full_2 = sum(DEFAULT_COSΘ_WEIGHTS .* f_gaussian.(DEFAULT_COSΘ_NODES))
error_full_2 = abs(result_full_2 - exact_value_2) / exact_value_2 * 100

result_half_2 = 2.0 * sum(DEFAULT_COSΘ_HALF_WEIGHTS .* f_gaussian.(DEFAULT_COSΘ_HALF_NODES))
error_half_2 = abs(result_half_2 - exact_value_2) / exact_value_2 * 100

println("\n精确值: ", exact_value_2)
println("\n全区间积分 (32节点):")
println("  数值结果: ", result_full_2)
println("  相对误差: ", error_full_2, " %")
println("\n半区间对称积分 (32节点):")
println("  数值结果: ", result_half_2)
println("  相对误差: ", error_half_2, " %")
println("\n精度提升: ", error_full_2 / error_half_2, "x")

# ============================================================================
# 测试3: 非对称函数 - 不能使用半区间方法
# ============================================================================
println("\n" * repeat("=", 70))
println("测试3: 非对称函数 f(x) = x³ (验证半区间方法不适用)")
println(repeat("=", 70))

f_asymmetric(x) = x^3

# 精确值: ∫₋₁¹ x³ dx = 0
exact_value_3 = 0.0

result_full_3 = sum(DEFAULT_COSΘ_WEIGHTS .* f_asymmetric.(DEFAULT_COSΘ_NODES))

println("\n精确值: ", exact_value_3)
println("\n全区间积分 (32节点):")
println("  数值结果: ", result_full_3)
println("  绝对误差: ", abs(result_full_3 - exact_value_3))
println("\n⚠️  对于非对称函数，不能使用半区间方法")

# ============================================================================
# 测试4: 物理应用示例 - 费米分布的角度部分
# ============================================================================
println("\n" * repeat("=", 70))
println("测试4: 物理应用 - 类似费米分布的对称函数")
println(repeat("=", 70))

# 模拟一个物理上常见的对称函数
function f_physics(cosθ, p=5.0, T=0.15)
    E = p * sqrt(1 + cosθ^2)  # 模拟能量依赖
    return 1.0 / (1.0 + exp(E / T))  # 费米分布形式
end

# 使用高精度参考
exact_value_4, _ = quadgk(x -> f_physics(x), -1.0, 1.0, rtol=1e-12)

result_full_4 = sum(DEFAULT_COSΘ_WEIGHTS .* f_physics.(DEFAULT_COSΘ_NODES))
error_full_4 = abs(result_full_4 - exact_value_4) / exact_value_4 * 100

result_half_4 = 2.0 * sum(DEFAULT_COSΘ_HALF_WEIGHTS .* f_physics.(DEFAULT_COSΘ_HALF_NODES))
error_half_4 = abs(result_half_4 - exact_value_4) / exact_value_4 * 100

println("\n精确值: ", exact_value_4)
println("\n全区间积分 (32节点):")
println("  数值结果: ", result_full_4)
println("  相对误差: ", error_full_4, " %")
println("\n半区间对称积分 (32节点):")
println("  数值结果: ", result_half_4)
println("  相对误差: ", error_half_4, " %")
println("\n精度提升: ", error_full_4 / error_half_4, "x")

# ============================================================================
# 使用建议
# ============================================================================
println("\n" * repeat("=", 70))
println("使用建议")
println(repeat("=", 70))
println("""
✓ 使用 DEFAULT_COSΘ_HALF_NODES 的场景:
  - 被积函数关于 cosθ = 0 对称: f(cosθ) = f(-cosθ)
  - 需要更高精度但不想增加节点数
  - 典型例子: 各向同性分布、偶数阶勒让德多项式等

✗ 不适用的场景:
  - 非对称函数
  - 积分区间本身不对称
  - 需要完整的角度信息

使用方法:
  using .GaussLegendre: DEFAULT_COSΘ_HALF_NODES, DEFAULT_COSΘ_HALF_WEIGHTS
  
  # 对于对称函数 f
  result = 2.0 * sum(DEFAULT_COSΘ_HALF_WEIGHTS .* f.(DEFAULT_COSΘ_HALF_NODES))
  
  # 注意乘以因子 2
""")
println("=" ^ 70)

