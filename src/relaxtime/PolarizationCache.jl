"""
# PolarizationCache.jl

极化函数计算缓存模块，通过哈希表避免重复计算相同参数的极化函数。

## 使用场景
在单次输运系数计算中，多个散射过程/散射道可能需要相同参数的极化函数：
- 同一散射过程的t/u道可能共享Π
- 不同散射过程（如uu→uu和dd→dd）共享Π_uu
- 相同(k0, k_norm)但不同通道(:P vs :S)需分别缓存

## 性能优势
- 极化函数计算成本：~0.1-1ms（包含高斯积分）
- 哈希表查询成本：~1ns
- 单次输运系数计算可能需要10⁴-10⁶次极化函数调用
- 缓存命中率可达30%-70%（取决于计算网格密度）

## 设计原则
1. 使用浮点数容差比较（避免舍入误差导致缓存失效）
2. 线程安全（如需并行计算）
3. 单次计算会话内有效（不跨计算任务持久化）
4. 提供缓存统计功能（命中率、节省时间）
"""
module PolarizationCache

include("PolarizationAniso.jl")
using .PolarizationAniso: polarization_aniso

# 导出函数
export polarization_aniso_cached, reset_cache!, get_cache_stats

@inline function _channel_code(channel::Symbol)::UInt8
    channel === :P && return UInt8(0x01)
    channel === :S && return UInt8(0x02)
    throw(ArgumentError("Unknown polarization channel: $channel (expected :P or :S)"))
end

# 缓存量化精度：通过截断 Float64 的 mantissa（有效位）把“浮点容差”落到离散桶。
# 好处：
# 1) Dict 的 hash/== 一致性由“离散桶”保证
# 2) 避免 round(sigdigits=...) 带来的 log10/^ 等昂贵开销（profile 热点）
# 3) 量化后的相对误差上界约为 2^(-CACHE_MANTISSA_BITS)
const CACHE_MANTISSA_BITS = 40  # 2^-40 ≈ 9.09e-13，接近原先 1e-12 的容差目标

const _SIGN_EXP_MASK = UInt64(0xFFF0000000000000)
const _MANTISSA_MASK = UInt64(0x000FFFFFFFFFFFFF)
const _MANTISSA_KEEP_MASK = _MANTISSA_MASK & ~((UInt64(1) << (52 - CACHE_MANTISSA_BITS)) - 1)

@inline function _round_cache(x::Float64)
    # 这些特殊值保持不变
    (x == 0.0 || isnan(x) || isinf(x)) && return x
    bits = reinterpret(UInt64, x)
    sign_exp = bits & _SIGN_EXP_MASK
    mantissa = bits & _MANTISSA_MASK
    return reinterpret(Float64, sign_exp | (mantissa & _MANTISSA_KEEP_MASK))
end

# ----------------------------------------------------------------------------
# 缓存键结构
# ----------------------------------------------------------------------------

"""
    PolarizationKey

极化函数缓存的键结构，包含所有影响计算结果的参数。

# 字段
- `channel`: 通道类型（:P或:S）
- `k0`: 能量分量（MeV）
- `k_norm`: 三动量大小（MeV）
- `m1, m2`: 两个夸克质量（MeV）
- `μ1, μ2`: 两个夸克化学势（MeV）
- `T`: 温度（MeV）
- `Φ, Φbar`: Polyakov环期望值
- `ξ`: 各向异性参数
- `A1, A2`: 预计算的A函数值
- `num_s_quark`: 奇异夸克数量（0或1）

# 注意
使用hash和==运算符时会考虑浮点数容差（EPS_CACHE）
"""
struct PolarizationKey
    channel_code::UInt8
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
    num_s_quark::UInt8

    function PolarizationKey(channel::Symbol, k0::Float64, k_norm::Float64,
                             m1::Float64, m2::Float64, μ1::Float64, μ2::Float64,
                             T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64,
                             A1::Float64, A2::Float64, num_s_quark::Int)
        return new(
            _channel_code(channel),
            _round_cache(k0),
            _round_cache(k_norm),
            _round_cache(m1),
            _round_cache(m2),
            _round_cache(μ1),
            _round_cache(μ2),
            _round_cache(T),
            _round_cache(Φ),
            _round_cache(Φbar),
            _round_cache(ξ),
            _round_cache(A1),
            _round_cache(A2),
            UInt8(num_s_quark),
        )
    end
end

@inline function Base.isequal(a::PolarizationKey, b::PolarizationKey)
    return a.channel_code == b.channel_code &&
           isequal(a.k0, b.k0) &&
           isequal(a.k_norm, b.k_norm) &&
           isequal(a.m1, b.m1) &&
           isequal(a.m2, b.m2) &&
           isequal(a.μ1, b.μ1) &&
           isequal(a.μ2, b.μ2) &&
           isequal(a.T, b.T) &&
           isequal(a.Φ, b.Φ) &&
           isequal(a.Φbar, b.Φbar) &&
           isequal(a.ξ, b.ξ) &&
           isequal(a.A1, b.A1) &&
           isequal(a.A2, b.A2) &&
           a.num_s_quark == b.num_s_quark
end

@inline function Base.hash(k::PolarizationKey, h::UInt)
    h = hash(k.channel_code, h)
    h = hash(k.k0, h)
    h = hash(k.k_norm, h)
    h = hash(k.m1, h)
    h = hash(k.m2, h)
    h = hash(k.μ1, h)
    h = hash(k.μ2, h)
    h = hash(k.T, h)
    h = hash(k.Φ, h)
    h = hash(k.Φbar, h)
    h = hash(k.ξ, h)
    h = hash(k.A1, h)
    h = hash(k.A2, h)
    h = hash(k.num_s_quark, h)
    return h
end

# ----------------------------------------------------------------------------
# 全局缓存存储
# ----------------------------------------------------------------------------

"""全局缓存字典：PolarizationKey => (Π_real, Π_imag)"""
const POLARIZATION_CACHE = Dict{PolarizationKey, Tuple{Float64, Float64}}()

"""缓存统计：总调用次数"""
const CACHE_TOTAL_CALLS = Ref(0)

"""缓存统计：缓存命中次数"""
const CACHE_HIT_CALLS = Ref(0)

# ----------------------------------------------------------------------------
# 带缓存的极化函数接口
# ----------------------------------------------------------------------------

"""
    polarization_aniso_cached(channel, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1_value, A2_value, num_s_quark)

带缓存的极化函数计算，自动查询缓存或调用底层计算函数。

# 参数
与`PolarizationAniso.polarization_aniso`完全相同

# 返回值
返回元组 `(Π_real, Π_imag)`

# 性能
- 缓存命中：~1ns（哈希表查询）
- 缓存未命中：~0.1-1ms（调用polarization_aniso计算并缓存结果）

# 示例
```julia
using PolarizationCache

# 第一次调用：计算并缓存
Π_real, Π_imag = polarization_aniso_cached(:P, 100.0, 50.0, 5.0, 5.0, 300.0, 300.0, 150.0, 0.5, 0.5, 0.0, -50.0, -50.0, 0)

# 第二次调用相同参数：从缓存读取（快~10⁵倍）
Π_real, Π_imag = polarization_aniso_cached(:P, 100.0, 50.0, 5.0, 5.0, 300.0, 300.0, 150.0, 0.5, 0.5, 0.0, -50.0, -50.0, 0)

# 查看缓存统计
stats = get_cache_stats()
println("缓存命中率: \$(stats.hit_rate * 100)%")
```
"""
function polarization_aniso_cached(channel::Symbol, k0::Float64, k_norm::Float64, 
                                  m1::Float64, m2::Float64, μ1::Float64, μ2::Float64, 
                                  T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, 
                                  A1_value::Float64, A2_value::Float64, num_s_quark::Int)
    # 更新统计
    CACHE_TOTAL_CALLS[] += 1
    
    # 构造缓存键
    key = PolarizationKey(channel, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1_value, A2_value, num_s_quark)
    
    # 单次查找：命中则直接返回，避免 haskey + getindex 的双查找
    cached = get(POLARIZATION_CACHE, key, nothing)
    if cached !== nothing
        CACHE_HIT_CALLS[] += 1
        return cached
    end

    # 缓存未命中，计算并缓存
    result = polarization_aniso(channel, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1_value, A2_value, num_s_quark)
    POLARIZATION_CACHE[key] = result
    return result
end

# ----------------------------------------------------------------------------
# 缓存管理函数
# ----------------------------------------------------------------------------

"""
    reset_cache!()

清空极化函数缓存和统计信息。

# 使用场景
- 开始新的输运系数计算任务
- 改变物理参数（温度、化学势等）后
- 内存不足需要释放缓存

# 注意
通常在单次输运系数计算的开始和结束时调用

# 示例
```julia
# 开始新的计算任务
reset_cache!()

# ... 进行输运系数计算 ...

# 结束后查看统计并清空
stats = get_cache_stats()
println("本次计算缓存命中率: \$(stats.hit_rate * 100)%")
reset_cache!()
```
"""
function reset_cache!()
    empty!(POLARIZATION_CACHE)
    CACHE_TOTAL_CALLS[] = 0
    CACHE_HIT_CALLS[] = 0
    return nothing
end

"""
    get_cache_stats()

获取缓存统计信息。

# 返回值
返回NamedTuple，包含以下字段：
- `total_calls`: 总调用次数
- `cache_hits`: 缓存命中次数
- `cache_misses`: 缓存未命中次数
- `hit_rate`: 缓存命中率（0.0-1.0）
- `cache_size`: 当前缓存条目数

# 示例
```julia
stats = get_cache_stats()
println("总调用: \$(stats.total_calls)")
println("命中次数: \$(stats.cache_hits)")
println("未命中次数: \$(stats.cache_misses)")
println("命中率: \$(stats.hit_rate * 100)%")
println("缓存大小: \$(stats.cache_size) 条目")
```
"""
function get_cache_stats()
    total = CACHE_TOTAL_CALLS[]
    hits = CACHE_HIT_CALLS[]
    misses = total - hits
    hit_rate = total > 0 ? hits / total : 0.0
    cache_size = length(POLARIZATION_CACHE)
    
    return (
        total_calls = total,
        cache_hits = hits,
        cache_misses = misses,
        hit_rate = hit_rate,
        cache_size = cache_size
    )
end

end # module PolarizationCache
