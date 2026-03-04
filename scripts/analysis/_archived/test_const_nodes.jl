"""
测试新的 const 常量节点使用方式
"""

include("../../src/integration/GaussLegendre.jl")
using .GaussLegendre

println("=" ^ 70)
println("GaussLegendre 模块默认常量测试")
println("=" ^ 70)

println("\n默认节点数量:")
println("  动量节点: ", length(DEFAULT_MOMENTUM_NODES))
println("  角度节点: ", length(DEFAULT_COSΘ_NODES))

println("\n积分区间:")
println("  动量: [", minimum(DEFAULT_MOMENTUM_NODES), ", ", maximum(DEFAULT_MOMENTUM_NODES), "] fm⁻¹")
println("  cosθ: [", minimum(DEFAULT_COSΘ_NODES), ", ", maximum(DEFAULT_COSΘ_NODES), "]")

println("\n前3个动量节点:")
for i in 1:3
    println("  p[$i] = ", DEFAULT_MOMENTUM_NODES[i], " fm⁻¹")
end

println("\n验证简单积分:")
# ∫₀¹⁰ x² dx = 10³/3 ≈ 333.333
f(x) = x^2
result = sum(DEFAULT_MOMENTUM_WEIGHTS .* f.(DEFAULT_MOMENTUM_NODES))
exact = 10.0^3 / 3.0
error = abs(result - exact) / exact * 100
println("  ∫₀¹⁰ x² dx")
println("    数值结果: ", result)
println("    精确值:   ", exact)
println("    相对误差: ", error, " %")

println("\n✓ 所有常量可以直接使用，无需调用函数")
println("=" ^ 70)

