"""
简单测试：验证 polarization_aniso_cached 的缓存统计是否工作
"""

push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
push!(LOAD_PATH, joinpath(@__DIR__, "../src/relaxtime"))

include("../src/Constants_PNJL.jl")
include("../src/integration/GaussLegendre.jl")
include("../src/relaxtime/OneLoopIntegrals.jl")
include("../src/relaxtime/PolarizationCache.jl")

using .PolarizationCache: polarization_aniso_cached, reset_cache!, get_cache_stats
using .GaussLegendre: gauleg
using .OneLoopIntegrals: A

# 设置参数
T = 150.0 / 197.327
m_u = 300.0 / 197.327
μ_u = 0.0
Φ = 0.5
Φbar = 0.5
ξ = 0.0
k0 = 100.0 / 197.327
k_norm = 50.0 / 197.327

nodes, weights = gauleg(0.0, 20.0, 64)
A_u = A(m_u, μ_u, T, Φ, Φbar, nodes, weights)

println("="^70)
println("测试 polarization_aniso_cached 缓存统计")
println("="^70)

reset_cache!()
stats1 = get_cache_stats()
println("\n调用前:")
println("  总调用: $(stats1.total_calls)")
println("  命中: $(stats1.cache_hits)")
println("  未命中: $(stats1.cache_misses)")
println("  缓存大小: $(stats1.cache_size)")

# 第一次调用
println("\n第一次调用 polarization_aniso_cached...")
result1 = polarization_aniso_cached(:P, k0, k_norm, m_u, m_u, μ_u, μ_u, T, Φ, Φbar, ξ, A_u, A_u, 0)
stats2 = get_cache_stats()
println("第1次调用后:")
println("  总调用: $(stats2.total_calls)")
println("  命中: $(stats2.cache_hits)")
println("  未命中: $(stats2.cache_misses)")
println("  缓存大小: $(stats2.cache_size)")
println("  结果: $(result1)")

# 第二次调用（应该命中缓存）
println("\n第二次调用相同参数（应命中缓存）...")
result2 = polarization_aniso_cached(:P, k0, k_norm, m_u, m_u, μ_u, μ_u, T, Φ, Φbar, ξ, A_u, A_u, 0)
stats3 = get_cache_stats()
println("第2次调用后:")
println("  总调用: $(stats3.total_calls)")
println("  命中: $(stats3.cache_hits)")
println("  未命中: $(stats3.cache_misses)")
println("  缓存大小: $(stats3.cache_size)")
println("  结果: $(result2)")

# 验证结果
println("\n="^70)
if stats3.total_calls == 2 && stats3.cache_hits == 1
    println("✓ 缓存统计正常工作！")
else
    println("✗ 缓存统计异常！")
    println("  预期：total_calls=2, cache_hits=1")
    println("  实际：total_calls=$(stats3.total_calls), cache_hits=$(stats3.cache_hits)")
end
println("="^70)
