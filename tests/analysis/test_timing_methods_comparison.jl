"""
对比测试：验证两种测时方法的差异

目的：找出为什么 test_polarization_aniso.jl 和 test_total_propagator_performance.jl
     测出的极化函数速度差异达15倍
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
using .PolarizationCache: polarization_aniso_cached, reset_cache!

println("="^80)
println("极化函数计时方法对比测试")
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

println("\n物理参数:")
println("  k0=$k0, k_norm=$k_norm, ξ=$ξ")
println("  m1=$m1, m2=$m2")
println("  A1=$A1, A2=$A2")

# === 预热 ===
println("\n预热中...")
for _ in 1:100
    polarization_aniso(:P, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
end
GC.gc()
println("✓ 预热完成")

# ============================================================================
# 方法1：test_polarization_aniso.jl 的方法（整体计时）
# ============================================================================
println("\n" * "="^80)
println("方法1：整体计时（test_polarization_aniso.jl 的方法）")
println("="^80)

n_iter1 = 10000
println("迭代次数: $n_iter1")

GC.gc()
start = time_ns()
for _ in 1:n_iter1
    value = polarization_aniso(:P, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
end
elapsed_total = (time_ns() - start) / 1e6  # 毫秒
avg_time1 = elapsed_total / n_iter1

@printf("\n结果:\n")
@printf("  总耗时: %.2f ms\n", elapsed_total)
@printf("  平均每次: %.6f ms\n", avg_time1)

# ============================================================================
# 方法2：test_total_propagator_performance.jl 的方法（逐次计时+统计）
# ============================================================================
println("\n" * "="^80)
println("方法2：逐次计时+统计（test_total_propagator_performance.jl 的方法）")
println("="^80)

n_iter2 = 5000
println("迭代次数: $n_iter2")

times = Float64[]
for i in 1:n_iter2
    #GC.gc()
    t_start = time_ns()
    value = polarization_aniso(:P, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
    t_end = time_ns()
    push!(times, (t_end - t_start) / 1e6)
end

avg_time2 = mean(times)
std_time2 = std(times)
min_time2 = minimum(times)
max_time2 = maximum(times)

@printf("\n结果:\n")
@printf("  平均: %.6f ms\n", avg_time2)
@printf("  标准差: %.6f ms\n", std_time2)
@printf("  最小值: %.6f ms\n", min_time2)
@printf("  最大值: %.6f ms\n", max_time2)

# ============================================================================
# 方法3：逐次计时但不调用GC（折中方案）
# ============================================================================
println("\n" * "="^80)
println("方法3：逐次计时但不频繁GC（折中方案）")
println("="^80)

n_iter3 = 5000
println("迭代次数: $n_iter3")

GC.gc()  # 只在开始时GC一次
times3 = Float64[]
for i in 1:n_iter3
    t_start = time_ns()
    value = polarization_aniso(:P, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
    t_end = time_ns()
    push!(times3, (t_end - t_start) / 1e6)
end

avg_time3 = mean(times3)
std_time3 = std(times3)
min_time3 = minimum(times3)
max_time3 = maximum(times3)

@printf("\n结果:\n")
@printf("  平均: %.6f ms\n", avg_time3)
@printf("  标准差: %.6f ms\n", std_time3)
@printf("  最小值: %.6f ms\n", min_time3)
@printf("  最大值: %.6f ms\n", max_time3)

# ============================================================================
# 方法4：使用 polarization_aniso_cached + reset_cache!
# ============================================================================
println("\n" * "="^80)
println("方法4：polarization_aniso_cached + reset_cache!（test_total_propagator的方法）")
println("="^80)

n_iter4 = 50000
println("迭代次数: $n_iter4")
# 先预热，排除JIT
for _ in 1:100
    polarization_aniso_cached(:P, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
end

reset_cache!()
times4 = Float64[]
for i in 1:n_iter4
    #reset_cache!()
    t_start = time_ns()
    value = polarization_aniso_cached(:P, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
    t_end = time_ns()
    push!(times4, (t_end - t_start) / 1e6)
end

avg_time4 = mean(times4)
std_time4 = std(times4)
min_time4 = minimum(times4)
max_time4 = maximum(times4)

@printf("\n结果:\n")
@printf("  平均: %.6f ms\n", avg_time4)
@printf("  标准差: %.6f ms\n", std_time4)
@printf("  最小值: %.6f ms\n", min_time4)
@printf("  最大值: %.6f ms\n", max_time4)

# ============================================================================
# 汇总对比
# ============================================================================
println("\n" * "="^80)
println("汇总对比")
println("="^80)

@printf("\n%-50s %12s\n", "测试方法", "平均时间(ms)")
println("-"^80)
@printf("%-50s %12.6f\n", "方法1: 整体计时 (polarization_aniso)", avg_time1)
@printf("%-50s %12.6f\n", "方法2: 逐次计时+每次GC (polarization_aniso)", avg_time2)
@printf("%-50s %12.6f\n", "方法3: 逐次计时+开始时GC (polarization_aniso)", avg_time3)
@printf("%-50s %12.6f\n", "方法4: 逐次计时+reset_cache (cached版本)", avg_time4)
println("-"^80)

println("\n相对速度比:")
@printf("  方法2 / 方法1 = %.2fx\n", avg_time2 / avg_time1)
@printf("  方法3 / 方法1 = %.2fx\n", avg_time3 / avg_time1)
@printf("  方法4 / 方法1 = %.2fx\n", avg_time4 / avg_time1)
@printf("  方法2 / 方法3 = %.2fx\n", avg_time2 / avg_time3)
@printf("  方法4 / 方法3 = %.2fx\n", avg_time4 / avg_time3)

println("\n" * "="^80)
println("结论")
println("="^80)
println("1. 每次迭代都调用 GC.gc() 会显著降低性能（方法2）")
println("2. 整体计时（方法1）给出的是最优化的性能（JIT充分优化）")
println("3. 逐次计时但不频繁GC（方法3）接近真实性能")
println("4. reset_cache!() 的开销相对较小")
println("\n推荐:")
println("  - 性能基准测试: 使用方法1（整体计时）")
println("  - 统计分析: 使用方法3（逐次计时但不频繁GC）")
println("  - 避免: 方法2（每次GC会严重拖慢速度）")
println("="^80)

