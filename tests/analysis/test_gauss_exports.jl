"""
验证 GaussLegendre 模块的所有导出符号
"""

include("../../src/integration/GaussLegendre.jl")
using .GaussLegendre

println("=" ^ 70)
println("GaussLegendre 模块导出符号列表")
println("=" ^ 70)

println("\n导出的符号:")
for name in names(GaussLegendre)
    println("  ✓ ", name)
end

println("\n" * "=" ^ 70)
println("常量检查")
println("=" ^ 70)

println("\n角度积分节点:")
println("  DEFAULT_COSΘ_NODES:       ", length(DEFAULT_COSΘ_NODES), " 个节点, 区间 [-1, 1]")
println("  DEFAULT_COSΘ_HALF_NODES:  ", length(DEFAULT_COSΘ_HALF_NODES), " 个节点, 区间 [0, 1]")

println("\n动量积分节点:")
println("  DEFAULT_MOMENTUM_NODES:   ", length(DEFAULT_MOMENTUM_NODES), " 个节点, 区间 [0, 10] fm⁻¹")

println("\n" * "=" ^ 70)
println("✓ 所有常量已成功导出并可用")
println("=" ^ 70)

