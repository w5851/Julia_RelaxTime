"""
分析测试5标准差异常大的原因

关键问题：为什么测试5的标准差(0.015ms)比测试1(0.002ms)大这么多？

假设：随着哈希表不断增长，哈希冲突和内存分配会导致性能劣化
"""

push!(LOAD_PATH, joinpath(@__DIR__, "../../src"))
push!(LOAD_PATH, joinpath(@__DIR__, "../../src/relaxtime"))

using Statistics
using Printf

include("../../src/Constants_PNJL.jl")
include("../../src/integration/GaussLegendre.jl")
include("../../src/relaxtime/OneLoopIntegrals.jl")
include("../../src/relaxtime/PolarizationAniso.jl")
include("../../src/relaxtime/PolarizationCache.jl")

using .GaussLegendre: gauleg
using .OneLoopIntegrals: A
using .PolarizationAniso: polarization_aniso
using .PolarizationCache: polarization_aniso_cached, reset_cache!, POLARIZATION_CACHE

println("="^80)
println("缓存大小对性能波动的影响分析")
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

# === 预热 ===
println("\n预热中...")
for _ in 1:500
    polarization_aniso(:P, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
    reset_cache!()
    polarization_aniso_cached(:P, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
end
GC.gc()
println("✓ 预热完成")

# ============================================================================
# 测试1：基准 - 纯极化函数（无缓存，多次重复）
# ============================================================================
println("\n" * "="^80)
println("基准测试：纯极化函数（相同参数，无缓存干扰）")
println("="^80)

n_samples = 10000
times_pure = Float64[]
for _ in 1:n_samples
    t_start = time_ns()
    polarization_aniso(:P, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
    t_end = time_ns()
    push!(times_pure, (t_end - t_start) / 1e6)
end

@printf("\n结果:\n")
@printf("  平均: %.6f ms\n", mean(times_pure))
@printf("  标准差: %.6f ms\n", std(times_pure))
@printf("  最小值: %.6f ms\n", minimum(times_pure))
@printf("  最大值: %.6f ms\n", maximum(times_pure))
@printf("  CV (变异系数): %.4f\n", std(times_pure) / mean(times_pure))

# ============================================================================
# 测试2：缓存未命中 - 缓存从小到大增长（模拟测试5）
# ============================================================================
println("\n" * "="^80)
println("测试：缓存未命中 - 缓存逐渐增长到1000条（模拟测试5）")
println("="^80)

reset_cache!()
n_samples = 1000
times_growing = Float64[]
cache_sizes = Int[]

for i in 1:n_samples
    k0_temp = k0 + i * 0.00001
    t_start = time_ns()
    polarization_aniso_cached(:P, k0_temp, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
    t_end = time_ns()
    push!(times_growing, (t_end - t_start) / 1e6)
    push!(cache_sizes, length(POLARIZATION_CACHE))
end

@printf("\n结果:\n")
@printf("  平均: %.6f ms\n", mean(times_growing))
@printf("  标准差: %.6f ms\n", std(times_growing))
@printf("  最小值: %.6f ms\n", minimum(times_growing))
@printf("  最大值: %.6f ms\n", maximum(times_growing))
@printf("  CV (变异系数): %.4f\n", std(times_growing) / mean(times_growing))
@printf("  最终缓存大小: %d\n", cache_sizes[end])

# ============================================================================
# 测试3：按缓存大小分段分析
# ============================================================================
println("\n" * "="^80)
println("分段分析：不同缓存大小阶段的性能")
println("="^80)

segments = [
    (1, 100, "前100次 (缓存: 1-100)"),
    (101, 300, "101-300次 (缓存: 101-300)"),
    (301, 600, "301-600次 (缓存: 301-600)"),
    (601, 1000, "601-1000次 (缓存: 601-1000)")
]

println("\n按缓存大小分段:")
println("-"^80)
@printf("%-30s %10s %10s %10s %10s\n", "阶段", "平均(ms)", "标准差(ms)", "最大值(ms)", "CV")
println("-"^80)

for (start_idx, end_idx, label) in segments
    segment_times = times_growing[start_idx:end_idx]
    @printf("%-30s %10.6f %10.6f %10.6f %10.4f\n", 
            label,
            mean(segment_times),
            std(segment_times),
            maximum(segment_times),
            std(segment_times) / mean(segment_times))
end

# ============================================================================
# 测试4：固定缓存大小（预填充后测试）
# ============================================================================
println("\n" * "="^80)
println("测试：固定缓存大小 - 预先填充1000条后继续测试")
println("="^80)

# 已经有1000条缓存了，继续添加1000条
times_fixed = Float64[]
for i in 1001:2000
    k0_temp = k0 + i * 0.00001
    t_start = time_ns()
    polarization_aniso_cached(:P, k0_temp, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
    t_end = time_ns()
    push!(times_fixed, (t_end - t_start) / 1e6)
end

@printf("\n结果 (缓存大小: 1000-2000):\n")
@printf("  平均: %.6f ms\n", mean(times_fixed))
@printf("  标准差: %.6f ms\n", std(times_fixed))
@printf("  最小值: %.6f ms\n", minimum(times_fixed))
@printf("  最大值: %.6f ms\n", maximum(times_fixed))
@printf("  CV (变异系数): %.4f\n", std(times_fixed) / mean(times_fixed))

# ============================================================================
# 测试5：极小缓存（每次只有1条）
# ============================================================================
println("\n" * "="^80)
println("测试：每次reset后只添加1条（最小缓存开销）")
println("="^80)

times_minimal = Float64[]
for i in 1:1000
    reset_cache!()
    k0_temp = k0 + i * 0.00001
    t_start = time_ns()
    polarization_aniso_cached(:P, k0_temp, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
    t_end = time_ns()
    push!(times_minimal, (t_end - t_start) / 1e6)
end

@printf("\n结果:\n")
@printf("  平均: %.6f ms\n", mean(times_minimal))
@printf("  标准差: %.6f ms\n", std(times_minimal))
@printf("  最小值: %.6f ms\n", minimum(times_minimal))
@printf("  最大值: %.6f ms\n", maximum(times_minimal))
@printf("  CV (变异系数): %.4f\n", std(times_minimal) / mean(times_minimal))

# ============================================================================
# 综合分析
# ============================================================================
println("\n" * "="^80)
println("综合分析：标准差差异的原因")
println("="^80)

println("\n标准差对比:")
println("-"^80)
@printf("%-40s %10.6f ms (CV=%.4f)\n", "基准(无缓存，相同参数)", std(times_pure), std(times_pure)/mean(times_pure))
@printf("%-40s %10.6f ms (CV=%.4f)\n", "缓存增长0→1000", std(times_growing), std(times_growing)/mean(times_growing))
@printf("%-40s %10.6f ms (CV=%.4f)\n", "固定缓存1000→2000", std(times_fixed), std(times_fixed)/mean(times_fixed))
@printf("%-40s %10.6f ms (CV=%.4f)\n", "每次reset(缓存=1)", std(times_minimal), std(times_minimal)/mean(times_minimal))

println("\n" * "="^80)
println("结论")
println("="^80)
println("\n测试5标准差大的主要原因:")
println("\n1. **哈希表动态扩容**")
println("   - Julia的Dict在增长时会重新分配内存并rehash")
println("   - 典型扩容点: 32, 64, 128, 256, 512, 1024...")
println("   - 扩容时的那次操作会显著变慢")
println()
println("2. **内存分配模式变化**")
println("   - 缓存从0增长到1000，内存使用持续增加")
println("   - 可能触发内存页分配、缓存行失效等")
println()
println("3. **哈希冲突增加**")
println("   - 缓存越大，哈希冲突概率越高")
println("   - 查找和插入时间略有增加")
println()
println("4. **CPU缓存效应**")
println("   - 小哈希表可能在CPU缓存中")
println("   - 大哈希表需要访问主内存")
println()
println("验证:")
println("  - 如果分段分析显示后期标准差更大 → 缓存大小影响")
println("  - 如果固定大小缓存标准差稳定 → 增长过程影响")
println("  - 如果每次reset标准差小 → 动态扩容是主因")

println("\n" * "="^80)

# ============================================================================
# 额外：可视化时间序列
# ============================================================================
println("\n" * "="^80)
println("时间序列可视化（检测扩容事件）")
println("="^80)

println("\n查找明显的性能尖峰（可能的扩容点）:")
println("-"^80)

# 找出超过平均值3倍标准差的点
threshold = mean(times_growing) + 3 * std(times_growing)
spikes = [(i, times_growing[i], cache_sizes[i]) for i in 1:length(times_growing) if times_growing[i] > threshold]

if !isempty(spikes)
    @printf("%-10s %-15s %-15s\n", "迭代", "耗时(ms)", "缓存大小")
    println("-"^80)
    for (idx, time, size) in spikes
        @printf("%-10d %-15.6f %-15d", idx, time, size)
        # 检查是否接近2的幂（典型扩容点）
        if size > 0
            log2_size = log2(size)
            if abs(log2_size - round(log2_size)) < 0.1
                print("  ← 可能的扩容点 (2^$(round(Int, log2_size)))")
            end
        end
        println()
    end
else
    println("未检测到显著的性能尖峰")
end

println("\n" * "="^80)

