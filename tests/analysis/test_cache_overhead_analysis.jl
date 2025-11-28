"""
缓存开销详细分析：区分JIT编译、哈希表构造和哈希查询的开销
"""

push!(LOAD_PATH, joinpath(@__DIR__, "../../src"))
push!(LOAD_PATH, joinpath(@__DIR__, "../../src/relaxtime"))

using Statistics
using Printf

include("../../src/Constants_PNJL.jl")
include("../../src/integration/GaussLegendre.jl")
include("../../src/relaxtime/OneLoopIntegrals.jl")
include("../../src/relaxtime/PolarizationAniso.jl")
include("../../src/relaxtime/TotalPropagator.jl")

using .GaussLegendre: gauleg
using .OneLoopIntegrals: A
using .PolarizationAniso: polarization_aniso
using .TotalPropagator: reset_cache!, get_cache_stats

# 从TotalPropagator导入缓存版本
using .TotalPropagator: PolarizationCache
const polarization_aniso_cached = PolarizationCache.polarization_aniso_cached

println("="^80)
println("缓存开销详细分析")
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

# ============================================================================
# 测试1：JIT编译开销（冷启动）
# ============================================================================
println("\n" * "="^80)
println("测试1：JIT编译开销（冷启动 - 第一次调用）")
println("="^80)

reset_cache!()
GC.gc()

t_start = time_ns()
result1 = polarization_aniso_cached(:P, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
t_jit = (time_ns() - t_start) / 1e6

stats1 = get_cache_stats()
@printf("第一次调用时间（含JIT编译）: %.4f ms\n", t_jit)
@printf("缓存统计: total_calls=%d, cache_hits=%d, cache_size=%d\n", 
        stats1.total_calls, stats1.cache_hits, stats1.cache_size)

# ============================================================================
# 测试2：缓存未命中开销（哈希查询+计算+插入）
# ============================================================================
println("\n" * "="^80)
println("测试2：缓存未命中开销（哈希查询失败 + 实际计算 + 插入缓存）")
println("="^80)

# 第二次调用，JIT已编译，但参数不同会导致缓存未命中
reset_cache!()
GC.gc()

t_start = time_ns()
result2 = polarization_aniso_cached(:P, k0+0.01, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
t_miss = (time_ns() - t_start) / 1e6

stats2 = get_cache_stats()
@printf("缓存未命中时间（JIT已完成）: %.4f ms\n", t_miss)
@printf("缓存统计: total_calls=%d, cache_hits=%d, cache_size=%d\n", 
        stats2.total_calls, stats2.cache_hits, stats2.cache_size)

# ============================================================================
# 测试3：缓存命中开销（仅哈希查询）
# ============================================================================
println("\n" * "="^80)
println("测试3：缓存命中开销（仅哈希查询）")
println("="^80)

# 第三次调用，使用相同参数，会命中缓存
GC.gc()

t_start = time_ns()
result3 = polarization_aniso_cached(:P, k0+0.01, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
t_hit = (time_ns() - t_start) / 1e6

stats3 = get_cache_stats()
@printf("缓存命中时间: %.6f ms\n", t_hit)
@printf("缓存统计: total_calls=%d, cache_hits=%d, cache_size=%d\n", 
        stats3.total_calls, stats3.cache_hits, stats3.cache_size)

# ============================================================================
# 测试4：不带缓存的直接计算（基准）
# ============================================================================
println("\n" * "="^80)
println("测试4：不带缓存的直接计算（基准参考）")
println("="^80)

# 预热
for _ in 1:10
    polarization_aniso(:P, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
end

GC.gc()
times_direct = Float64[]
for _ in 1:1000
    t_start = time_ns()
    polarization_aniso(:P, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
    t_end = time_ns()
    push!(times_direct, (t_end - t_start) / 1e6)
end

avg_direct = mean(times_direct)
std_direct = std(times_direct)

@printf("直接计算平均时间: %.4f ms (标准差: %.4f ms)\n", avg_direct, std_direct)

# ============================================================================
# 测试5：哈希表操作的纯开销
# ============================================================================
println("\n" * "="^80)
println("测试5：Julia Dict 哈希表操作基准测试")
println("="^80)

# 创建一个与PolarizationKey类似结构的测试
struct TestKey
    channel::Symbol
    k0::Float64
    k_norm::Float64
    m1::Float64
    m2::Float64
    μ1::Float64
    μ2::Float64
    T::Float64
    Φ::Float64
    Φbar::Float64
    ξ::Float64
    A1::Float64
    A2::Float64
    num_s_quark::Int
end

# 自定义hash和==（与PolarizationCache类似）
function Base.hash(key::TestKey, h::UInt)
    round_val(x) = round(x, sigdigits=12)
    h = hash(key.channel, h)
    h = hash(round_val(key.k0), h)
    h = hash(round_val(key.k_norm), h)
    h = hash(round_val(key.m1), h)
    h = hash(round_val(key.m2), h)
    h = hash(round_val(key.μ1), h)
    h = hash(round_val(key.μ2), h)
    h = hash(round_val(key.T), h)
    h = hash(round_val(key.Φ), h)
    h = hash(round_val(key.Φbar), h)
    h = hash(round_val(key.ξ), h)
    h = hash(round_val(key.A1), h)
    h = hash(round_val(key.A2), h)
    h = hash(key.num_s_quark, h)
    return h
end

const EPS = 1e-12
function Base.:(==)(k1::TestKey, k2::TestKey)
    k1.channel != k2.channel && return false
    k1.num_s_quark != k2.num_s_quark && return false
    
    function approx_equal(a::Float64, b::Float64)
        abs_max = max(abs(a), abs(b))
        abs_max < 1e-15 && return true
        return abs(a - b) / abs_max < EPS
    end
    
    approx_equal(k1.k0, k2.k0) || return false
    approx_equal(k1.k_norm, k2.k_norm) || return false
    approx_equal(k1.m1, k2.m1) || return false
    approx_equal(k1.m2, k2.m2) || return false
    approx_equal(k1.μ1, k2.μ1) || return false
    approx_equal(k1.μ2, k2.μ2) || return false
    approx_equal(k1.T, k2.T) || return false
    approx_equal(k1.Φ, k2.Φ) || return false
    approx_equal(k1.Φbar, k2.Φbar) || return false
    approx_equal(k1.ξ, k2.ξ) || return false
    approx_equal(k1.A1, k2.A1) || return false
    approx_equal(k1.A2, k2.A2) || return false
    
    return true
end

# 测试哈希表插入
test_dict = Dict{TestKey, Tuple{Float64, Float64}}()
test_key = TestKey(:P, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)

times_insert = Float64[]
for i in 1:10000
    key_i = TestKey(:P, k0+i*1e-6, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
    t_start = time_ns()
    test_dict[key_i] = (1.0, 2.0)
    t_end = time_ns()
    push!(times_insert, (t_end - t_start) / 1e6)
end

@printf("哈希表插入平均时间: %.6f ms (标准差: %.6f ms)\n", 
        mean(times_insert), std(times_insert))

# 测试哈希表查询（命中）
times_lookup = Float64[]
for i in 1:10000
    key_i = TestKey(:P, k0+i*1e-6, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, 0)
    t_start = time_ns()
    result = get(test_dict, key_i, (0.0, 0.0))
    t_end = time_ns()
    push!(times_lookup, (t_end - t_start) / 1e6)
end

@printf("哈希表查询平均时间: %.6f ms (标准差: %.6f ms)\n", 
        mean(times_lookup), std(times_lookup))

# ============================================================================
# 汇总分析
# ============================================================================
println("\n" * "="^80)
println("开销分解汇总")
println("="^80)

@printf("\n%-50s %12s\n", "操作类型", "时间 (ms)")
println("-"^80)
@printf("%-50s %12.4f  ← 包含JIT编译\n", "首次调用（冷启动）", t_jit)
@printf("%-50s %12.4f\n", "缓存未命中（JIT后，含计算+插入）", t_miss)
@printf("%-50s %12.6f  ← 仅哈希查询\n", "缓存命中（仅查询）", t_hit)
@printf("%-50s %12.4f\n", "直接计算（无缓存，基准）", avg_direct)
println("-"^80)
@printf("%-50s %12.6f\n", "纯哈希表插入操作", mean(times_insert))
@printf("%-50s %12.6f\n", "纯哈希表查询操作", mean(times_lookup))
println("-"^80)

println("\n开销分析:")
jit_overhead = t_jit - t_miss
cache_overhead = t_miss - avg_direct
hash_query_overhead = t_hit

@printf("  JIT编译开销: %.4f ms (%.1f%%)\n", jit_overhead, (jit_overhead/t_jit)*100)
@printf("  缓存管理开销: %.6f ms (包含哈希计算+比较+插入)\n", cache_overhead)
@printf("  哈希查询开销: %.6f ms (命中时)\n", hash_query_overhead)
@printf("  实际计算时间: %.4f ms (基准)\n", avg_direct)

println("\n" * "="^80)
println("结论")
println("="^80)
println("1. **71ms的最大值主要是JIT编译开销，不是哈希表建立开销**")
@printf("   - JIT编译: ~%.1f ms\n", jit_overhead)
@printf("   - 实际首次计算: ~%.2f ms\n", t_miss)
println()
println("2. **哈希表操作开销非常小**")
@printf("   - 插入操作: ~%.6f ms (微秒级)\n", mean(times_insert))
@printf("   - 查询操作: ~%.6f ms (微秒级)\n", mean(times_lookup))
println()
println("3. **缓存效率分析**")
@printf("   - 缓存未命中: %.4f ms (计算+插入)\n", t_miss)
@printf("   - 缓存命中: %.6f ms (仅查询)\n", t_hit)
@printf("   - 加速比: %.0fx\n", t_miss/t_hit)
println()
println("4. **优化建议**")
println("   - 哈希表本身不需要优化（已经足够快）")
println("   - JIT编译是一次性开销，后续调用不受影响")
println("   - 缓存命中率才是关键性能指标")
println("="^80)

