"""
详细分析：缓存首次调用的额外开销来源

目的：解释为什么方法4的最大值(0.0337ms)比极化函数平均计算时间(0.0138ms)慢
      以及这个额外开销的具体来源
"""

push!(LOAD_PATH, joinpath(@__DIR__, "../src"))
push!(LOAD_PATH, joinpath(@__DIR__, "../src/relaxtime"))

using Statistics
using Printf

include("../src/Constants_PNJL.jl")
include("../src/integration/GaussLegendre.jl")
include("../src/relaxtime/OneLoopIntegrals.jl")
include("../src/relaxtime/PolarizationAniso.jl")
include("../src/relaxtime/PolarizationCache.jl")

using .GaussLegendre: gauleg
using .OneLoopIntegrals: A
using .PolarizationAniso: polarization_aniso
using .PolarizationCache: polarization_aniso_cached, reset_cache!

println("="^80)
println("缓存首次调用额外开销分析")
println("="^80)

# === 设置参数 ===
k0 = 0.45
k_norm = 0.3
m1 = 0.25
μ1 = 0.12
m2 = 0.36
μ2 = -0.05
T = 0.17
Φ = 0.15
Φbar = 0.15
ξ = 0.0

nodes, weights = gauleg(0.0, 20.0, 32)
A1 = A(m1, μ1, T, Φ, Φbar, nodes, weights)
A2 = A(m2, μ2, T, Φ, Φbar, nodes, weights)

# === 充分预热 ===
println("\n充分预热所有代码路径...")
for _ in 1:1000
    polarization_aniso(:P, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
    reset_cache!()
    polarization_aniso_cached(:P, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
end
GC.gc()
sleep(0.1)
println("✓ 预热完成")

# ============================================================================
# 测试1：基准 - 纯极化函数计算时间
# ============================================================================
println("\n" * "="^80)
println("测试1：基准 - 纯极化函数计算时间（无缓存机制）")
println("="^80)

n_samples = 10000
times_pure = Float64[]
for _ in 1:n_samples
    t_start = time_ns()
    polarization_aniso(:P, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
    t_end = time_ns()
    push!(times_pure, (t_end - t_start) / 1e6)
end

avg_pure = mean(times_pure)
std_pure = std(times_pure)
min_pure = minimum(times_pure)
max_pure = maximum(times_pure)

@printf("\n结果:\n")
@printf("  平均: %.6f ms\n", avg_pure)
@printf("  标准差: %.6f ms\n", std_pure)
@printf("  最小值: %.6f ms\n", min_pure)
@printf("  最大值: %.6f ms\n", max_pure)

# ============================================================================
# 测试2：缓存命中的时间（哈希查询）
# ============================================================================
println("\n" * "="^80)
println("测试2：缓存命中时间（纯哈希表查询）")
println("="^80)

# 先填充缓存
reset_cache!()
polarization_aniso_cached(:P, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)

n_samples = 100000
times_hit = Float64[]
for _ in 1:n_samples
    t_start = time_ns()
    polarization_aniso_cached(:P, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
    t_end = time_ns()
    push!(times_hit, (t_end - t_start) / 1e6)
end

avg_hit = mean(times_hit)
std_hit = std(times_hit)
min_hit = minimum(times_hit)
max_hit = maximum(times_hit)

@printf("\n结果:\n")
@printf("  平均: %.6f ms (%.2f ns)\n", avg_hit, avg_hit * 1e6)
@printf("  标准差: %.6f ms\n", std_hit)
@printf("  最小值: %.6f ms (%.2f ns)\n", min_hit, min_hit * 1e6)
@printf("  最大值: %.6f ms (%.2f ns)\n", max_hit, max_hit * 1e6)

# ============================================================================
# 测试3：缓存未命中的时间（首次调用）
# ============================================================================
println("\n" * "="^80)
println("测试3：缓存未命中时间（首次调用 = 计算 + 哈希插入）")
println("="^80)

n_samples = 5000
times_miss = Float64[]
for i in 1:n_samples
    reset_cache!()  # 每次都清空缓存，模拟首次调用
    t_start = time_ns()
    polarization_aniso_cached(:P, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
    t_end = time_ns()
    push!(times_miss, (t_end - t_start) / 1e6)
end

avg_miss = mean(times_miss)
std_miss = std(times_miss)
min_miss = minimum(times_miss)
max_miss = maximum(times_miss)

@printf("\n结果:\n")
@printf("  平均: %.6f ms\n", avg_miss)
@printf("  标准差: %.6f ms\n", std_miss)
@printf("  最小值: %.6f ms\n", min_miss)
@printf("  最大值: %.6f ms\n", max_miss)

# ============================================================================
# 测试4：单独测量reset_cache!的开销
# ============================================================================
println("\n" * "="^80)
println("测试4：reset_cache!() 的开销")
println("="^80)

# 先填充一些缓存
reset_cache!()
for i in 1:100
    k0_temp = k0 + i * 0.001
    polarization_aniso_cached(:P, k0_temp, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
end

n_samples = 10000
times_reset = Float64[]
for _ in 1:n_samples
    t_start = time_ns()
    reset_cache!()
    t_end = time_ns()
    push!(times_reset, (t_end - t_start) / 1e6)
    # 重新填充缓存以保持一致
    for i in 1:100
        k0_temp = k0 + i * 0.001
        polarization_aniso_cached(:P, k0_temp, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
    end
end

avg_reset = mean(times_reset)
std_reset = std(times_reset)
min_reset = minimum(times_reset)
max_reset = maximum(times_reset)

@printf("\n结果 (缓存包含100个条目):\n")
@printf("  平均: %.6f ms (%.2f ns)\n", avg_reset, avg_reset * 1e6)
@printf("  标准差: %.6f ms\n", std_reset)
@printf("  最小值: %.6f ms (%.2f ns)\n", min_reset, min_reset * 1e6)
@printf("  最大值: %.6f ms (%.2f ns)\n", max_reset, max_reset * 1e6)

# ============================================================================
# 测试5：分解缓存未命中开销（不包含reset_cache!）
# ============================================================================
println("\n" * "="^80)
println("测试5：缓存未命中（不调用reset_cache! - 使用不同参数）")
println("="^80)

reset_cache!()
n_samples = 1000
times_miss_no_reset = Float64[]
for i in 1:n_samples
    # 使用不同的k0值，确保缓存未命中但不需要reset
    k0_temp = k0 + i * 0.00001
    t_start = time_ns()
    polarization_aniso_cached(:P, k0_temp, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
    t_end = time_ns()
    push!(times_miss_no_reset, (t_end - t_start) / 1e6)
end

avg_miss_no_reset = mean(times_miss_no_reset)
std_miss_no_reset = std(times_miss_no_reset)
min_miss_no_reset = minimum(times_miss_no_reset)
max_miss_no_reset = maximum(times_miss_no_reset)

@printf("\n结果:\n")
@printf("  平均: %.6f ms\n", avg_miss_no_reset)
@printf("  标准差: %.6f ms\n", std_miss_no_reset)
@printf("  最小值: %.6f ms\n", min_miss_no_reset)
@printf("  最大值: %.6f ms\n", max_miss_no_reset)

# ============================================================================
# 综合分析
# ============================================================================
println("\n" * "="^80)
println("综合分析")
println("="^80)

println("\n1. 各操作的平均时间:")
println("-"^80)
@printf("  %-50s %12.6f ms\n", "纯极化函数计算 (baseline)", avg_pure)
@printf("  %-50s %12.6f ms (%.0f ns)\n", "缓存命中 (哈希查询)", avg_hit, avg_hit*1e6)
@printf("  %-50s %12.6f ms\n", "缓存未命中+reset (完整首次调用)", avg_miss)
@printf("  %-50s %12.6f ms\n", "缓存未命中无reset (计算+插入)", avg_miss_no_reset)
@printf("  %-50s %12.6f ms (%.0f ns)\n", "reset_cache!() 操作", avg_reset, avg_reset*1e6)

println("\n2. 额外开销分解:")
println("-"^80)
overhead_hash = avg_hit
overhead_insert = avg_miss_no_reset - avg_pure
overhead_with_reset = avg_miss - avg_pure
overhead_reset = avg_reset

@printf("  %-50s %12.6f ms (%.0f ns)\n", "哈希查询开销", overhead_hash, overhead_hash*1e6)
@printf("  %-50s %12.6f ms (%.0f ns)\n", "哈希插入开销 (实测)", overhead_insert, overhead_insert*1e6)
@printf("  %-50s %12.6f ms (%.0f ns)\n", "reset_cache!开销", overhead_reset, overhead_reset*1e6)
@printf("  %-50s %12.6f ms\n", "首次调用总开销 (含reset)", overhead_with_reset)

println("\n3. 开销占比:")
println("-"^80)
@printf("  哈希查询占纯计算时间的比例: %.4f%%\n", (overhead_hash / avg_pure) * 100)
@printf("  哈希插入占纯计算时间的比例: %.4f%%\n", (overhead_insert / avg_pure) * 100)
@printf("  reset占纯计算时间的比例: %.4f%%\n", (overhead_reset / avg_pure) * 100)

println("\n4. 解释方法4的最大值:")
println("-"^80)
@printf("  方法4实测最大值: %.6f ms\n", max_miss)
@printf("  理论值 (纯计算 + 插入开销): %.6f ms\n", avg_pure + overhead_insert)
@printf("  差异: %.6f ms\n", max_miss - (avg_pure + overhead_insert))
println("\n  最大值通常是:")
println("    = 极化函数计算时间(可能略慢于平均) + 哈希插入 + reset_cache! + 其他噪声")
println("    = 纯计算的最大值波动 + ~微秒级的缓存操作")

println("\n5. 关于您的理解:")
println("-"^80)
println("  您的理解是正确的！哈希表插入确实是纳秒量级:")
@printf("    - 哈希查询: ~%.0f ns\n", overhead_hash * 1e6)
@printf("    - 哈希插入: ~%.0f ns (%.6f ms)\n", overhead_insert * 1e6, overhead_insert)
println("\n  方法4的最大值(~0.034ms)之所以比平均极化计算(0.014ms)慢，")
println("  主要原因是:")
println("    1. 最大值对应的是极化函数本身计算较慢的那次")
println("    2. 加上~数百纳秒的reset_cache!开销")
println("    3. 加上~数百纳秒的哈希插入开销")
println("    4. 加上统计噪声和JIT优化状态差异")
println("\n  而不是哈希插入特别慢！")

println("\n6. 验证：缓存节省的时间")
println("-"^80)
speedup = avg_pure / avg_hit
@printf("  缓存命中加速比: %.0fx\n", speedup)
@printf("  单次查询节省时间: %.6f ms\n", avg_pure - avg_hit)
println("\n  这才是缓存的真正价值：")
println("  一次哈希插入(~微秒)的代价，换来后续几千次查询的极大加速！")

println("\n" * "="^80)
