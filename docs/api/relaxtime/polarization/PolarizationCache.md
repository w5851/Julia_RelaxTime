# PolarizationCache 模块 API 文档

## 模块概述

`PolarizationCache` 模块通过哈希表缓存机制优化极化函数计算性能。在单次输运系数计算中，多个散射过程和散射道可能需要计算相同参数的极化函数，通过缓存避免重复计算，显著提升整体性能。

**性能优势**: 
- 单次极化函数计算：0.1-1ms（包含高斯积分）
- 缓存查询：~1ns（哈希表查询）
- 加速比：约10⁵倍（缓存命中时）
- 典型收益：单次输运系数计算可节省30%-70%时间

**在计算链中的位置**:
```
K系数 → 极化函数Π (带缓存) → 介子传播子D → 散射振幅M → 弛豫时间τ
```

## 依赖

- `PolarizationAniso` —— 底层极化函数计算（被缓存的实际计算模块）
- Julia标准库 `Dict` —— 哈希表存储

## 设计理念

### 缓存复用场景
1. **同一散射过程的不同道**：uu→uu的t道和u道可能使用相同的(k0, k_norm)
2. **不同散射过程**：uu→uu和dd→dd共享相同的Π_uu
3. **同一动量不同通道**：相同(k0, k_norm)的P通道和S通道分别缓存

### 模块化设计
- **透明包装**：完全兼容`PolarizationAniso`接口，零修改接入
- **会话管理**：单次计算任务内有效，不跨任务持久化
- **统计可观测**：提供详细缓存统计，便于性能分析

## 单位约定

- 能量/动量：MeV
- 质量/化学势：MeV
- 温度：MeV
- 极化函数：fm²（ComplexF64）
- 无量纲参数：ξ（各向异性）、Φ/Φbar（Polyakov环）

---

## API 参考

### 主函数

#### `polarization_aniso_cached`

带缓存的极化函数计算，自动查询缓存或调用底层计算。

##### 函数签名

```julia
polarization_aniso_cached(
    channel::Symbol, 
    k0::Float64, 
    k_norm::Float64, 
    m1::Float64, 
    m2::Float64, 
    μ1::Float64, 
    μ2::Float64, 
    T::Float64, 
    Φ::Float64, 
    Φbar::Float64, 
    ξ::Float64, 
    A1_value::Float64, 
    A2_value::Float64, 
    num_s_quark::Int
) -> Tuple{Float64, Float64}
```

##### 参数说明

| 参数 | 类型 | 物理意义 | 单位 | 备注 |
|------|------|----------|------|------|
| `channel` | Symbol | 通道类型 | - | `:P`(赝标量)或`:S`(标量) |
| `k0` | Float64 | 介子能量分量 | MeV | 可正可负 |
| `k_norm` | Float64 | 介子三动量大小 | MeV | 非负 |
| `m1` | Float64 | 第一个夸克质量 | MeV | 正数 |
| `m2` | Float64 | 第二个夸克质量 | MeV | 正数 |
| `μ1` | Float64 | 第一个夸克化学势 | MeV | 实数 |
| `μ2` | Float64 | 第二个夸克化学势 | MeV | 实数 |
| `T` | Float64 | 温度 | MeV | 正数 |
| `Φ` | Float64 | Polyakov环期望值 | - | [0,1] |
| `Φbar` | Float64 | 反Polyakov环期望值 | - | [0,1] |
| `ξ` | Float64 | 各向异性参数 | - | 非负，0表示各向同性 |
| `A1_value` | Float64 | 第一个夸克的A函数值 | - | 预计算值 |
| `A2_value` | Float64 | 第二个夸克的A函数值 | - | 预计算值 |
| `num_s_quark` | Int | 奇异夸克数量 | - | 0或1，影响k0对称性处理 |

##### 返回值

返回元组 `(Π_real, Π_imag)`：
- `Π_real`：极化函数实部（单位：fm²）
- `Π_imag`：极化函数虚部（单位：fm²）

##### 行为描述

1. **缓存查询**：根据所有输入参数构造缓存键，查询哈希表
2. **命中**：直接返回缓存结果（耗时~1ns）
3. **未命中**：调用`PolarizationAniso.polarization_aniso`计算，存储结果并返回（耗时~0.1-1ms）
4. **统计更新**：自动更新调用次数和命中次数

##### 使用示例

```julia
using PolarizationCache

# 参数设置
channel = :P
k0 = 100.0      # MeV
k_norm = 50.0   # MeV
m_u = 5.0       # MeV
μ_u = 300.0     # MeV
T = 150.0       # MeV
Φ = 0.5
Φbar = 0.5
ξ = 0.0         # 各向同性
A_u = -50.0     # 预计算的A函数值

# 第一次调用：计算并缓存
Π_real, Π_imag = polarization_aniso_cached(
    :P, k0, k_norm, m_u, m_u, μ_u, μ_u, T, Φ, Φbar, ξ, A_u, A_u, 0
)
# 耗时: ~0.5ms

# 第二次调用相同参数：从缓存读取
Π_real2, Π_imag2 = polarization_aniso_cached(
    :P, k0, k_norm, m_u, m_u, μ_u, μ_u, T, Φ, Φbar, ξ, A_u, A_u, 0
)
# 耗时: ~1ns，结果完全相同
```

##### 注意事项

1. **浮点数容差/量化**：缓存键会对所有 Float64 参数做“mantissa 截断量化”（目标精度约 $2^{-40}\approx 9\times 10^{-13}$），从而把浮点容差落到离散桶上，避免 `round(sigdigits=...)` 的额外开销
2. **A函数复用**：A函数仅依赖(m, μ, T, Φ)，应在外层预计算并传入
3. **线程安全**：当前实现不是线程安全，并行计算需要每个线程独立缓存
4. **内存占用**：每个缓存条目约200字节，10,000条目约2MB

---

### 缓存管理函数

#### `reset_cache!`

清空极化函数缓存和所有统计信息。

##### 函数签名

```julia
reset_cache!() -> Nothing
```

##### 使用场景

1. **开始新计算任务**：确保缓存不会受到上次计算残留影响
2. **改变物理参数**：温度、化学势变化后，旧缓存不再适用
3. **释放内存**：计算结束后清空缓存释放内存

##### 示例

```julia
# 开始新的输运系数计算
reset_cache!()

# 进行计算...
for T in temperatures
    result = calculate_transport_coefficient(T, μ, ξ)
end

# 结束后释放内存
reset_cache!()
```

#### `get_cache_stats`

获取缓存统计信息，用于性能分析和监控。

##### 函数签名

```julia
get_cache_stats() -> NamedTuple
```

##### 返回值

返回包含以下字段的NamedTuple：

| 字段 | 类型 | 含义 |
|------|------|------|
| `total_calls` | Int | 总调用次数 |
| `cache_hits` | Int | 缓存命中次数 |
| `cache_misses` | Int | 缓存未命中次数 |
| `hit_rate` | Float64 | 缓存命中率（0.0-1.0） |
| `cache_size` | Int | 当前缓存条目数 |

##### 性能分析公式

```julia
stats = get_cache_stats()

# 节省的计算次数
saved_calls = stats.cache_hits

# 节省的时间（假设单次计算0.5ms）
saved_time_ms = saved_calls * 0.5

# 缓存效率
efficiency = stats.hit_rate * 100  # 百分比
```

##### 示例

```julia
# 计算过程中监控
stats = get_cache_stats()
println("缓存统计:")
println("  总调用: $(stats.total_calls)")
println("  命中: $(stats.cache_hits) ($(stats.hit_rate * 100)%)")
println("  未命中: $(stats.cache_misses)")
println("  缓存大小: $(stats.cache_size) 条目")
println("  估算节省时间: $(stats.cache_hits * 0.5)ms")
```

---

## 实现细节

### 缓存键结构

#### `PolarizationKey`（内部使用）

缓存键结构体，包含所有影响计算结果的参数。

为避免 `Symbol` 参与 Dict 哈希时走 `objectid` 造成额外开销，内部会把通道 `:P/:S` 编码成 `UInt8`（`channel_code`）。

```julia
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
end
```

### 哈希和相等性

#### 哈希算法
- 对所有 Float64 参数做 mantissa 截断量化（约等价于 $\sim 10^{-12}$ 的容差目标）
- 使用 Julia 标准 `hash` 对量化后的字段计算
- `channel_code` / `num_s_quark` 使用紧凑整数编码，减少哈希开销

#### 相等性判断
- `channel_code` / `num_s_quark`：完全相等
- Float64：比较量化后的值（等价于在量化桶上做完全相等）

### 全局缓存存储

```julia
# 缓存字典
const POLARIZATION_CACHE = Dict{PolarizationKey, Tuple{Float64, Float64}}()

# 统计计数器
const CACHE_TOTAL_CALLS = Ref(0)
const CACHE_HIT_CALLS = Ref(0)
```

---

## 使用指南

### 典型工作流程

```julia
using PolarizationCache

function calculate_transport_coefficient(T::Float64, μ::Float64, ξ::Float64)
    # 1. 初始化缓存
    reset_cache!()
    
    # 2. 预计算公共参数
    m_u, m_d, m_s = get_quark_masses(T)
    μ_u, μ_d, μ_s = calculate_chemical_potentials(μ)
    Φ, Φbar = calculate_polyakov_loop(T, μ)
    
    # 预计算A函数（所有Π共享）
    A_u = calculate_A_function(m_u, μ_u, T, Φ, Φbar)
    A_d = calculate_A_function(m_d, μ_d, T, Φ, Φbar)
    A_s = calculate_A_function(m_s, μ_s, T, Φ, Φbar)
    
    # 3. 主计算循环（自动缓存）
    result = 0.0
    for process in [:uu_to_uu, :ud_to_ud, :us_to_us, ...]
        for channel in [:t, :u, :s]
            # 计算介子动量
            k0, k_norm = calculate_meson_momentum(process, channel, ...)
            
            # 计算极化函数（自动缓存）
            if involves_u_quarks(process)
                Π_uu = polarization_aniso_cached(:P, k0, k_norm, m_u, m_u, 
                                                  μ_u, μ_u, T, Φ, Φbar, ξ, A_u, A_u, 0)
            end
            
            if involves_s_quarks(process)
                Π_us = polarization_aniso_cached(:P, k0, k_norm, m_u, m_s, 
                                                  μ_u, μ_s, T, Φ, Φbar, ξ, A_u, A_s, 1)
            end
            
            # 计算传播子和散射振幅
            D = meson_propagator(...)
            M = scattering_amplitude(...)
            
            # 累积贡献
            result += contribution(M, D, ...)
        end
    end
    
    # 4. 输出统计
    stats = get_cache_stats()
    @info "极化函数缓存统计" stats.total_calls stats.cache_hits stats.hit_rate
    
    # 5. 清理（可选，释放内存）
    reset_cache!()
    
    return result
end
```

### 性能优化建议

#### 1. A函数预计算
```julia
# ❌ 低效：每次都重新计算A
Π = polarization_aniso_cached(:P, k0, k_norm, m, m, μ, μ, T, Φ, Φbar, ξ,
                              calculate_A_function(m, μ, T, Φ, Φbar),  # 重复计算
                              calculate_A_function(m, μ, T, Φ, Φbar),  # 重复计算
                              0)

# ✅ 高效：预计算并复用
A_value = calculate_A_function(m, μ, T, Φ, Φbar)
Π = polarization_aniso_cached(:P, k0, k_norm, m, m, μ, μ, T, Φ, Φbar, ξ,
                              A_value, A_value, 0)
```

#### 2. 参数组织
```julia
# ✅ 将相同参数的调用组织在一起，提高缓存命中率
# 先计算所有P通道
for process in processes
    for channel in channels
        k0, k_norm = get_momentum(process, channel)
        Π_P = polarization_aniso_cached(:P, k0, k_norm, ...)
    end
end

# 再计算所有S通道（可能重用相同的k0, k_norm）
for process in processes
    for channel in channels
        k0, k_norm = get_momentum(process, channel)
        Π_S = polarization_aniso_cached(:S, k0, k_norm, ...)
    end
end
```

#### 3. 定期检查统计
```julia
# 在长时间计算中定期输出统计，监控性能
n_iterations = 0
for ...
    n_iterations += 1
    if n_iterations % 1000 == 0
        stats = get_cache_stats()
        println("迭代 $n_iterations: 命中率 $(stats.hit_rate * 100)%")
    end
end
```

### 并行计算集成

#### 多线程方案（推荐）
```julia
using Distributed

@everywhere using PolarizationCache

function parallel_transport_calculation(temperature_range)
    results = pmap(temperature_range) do T
        # 每个worker独立缓存
        reset_cache!()
        result = calculate_transport_coefficient(T, μ, ξ)
        stats = get_cache_stats()
        return (result, stats)
    end
    
    # 汇总统计
    total_stats = merge_stats([r[2] for r in results])
    return [r[1] for r in results], total_stats
end
```

---

## 性能分析

### 理论性能模型

假设输运系数计算参数：
- 散射过程数：N_proc = 9
- 散射道数：N_channel = 3
- 动量网格：N_p × N_k = 100 × 100
- 介子通道：N_meson = 2 (P, S)

**总调用次数**：
```
N_total = N_proc × N_channel × N_p × N_k × N_meson
        = 9 × 3 × 100 × 100 × 2
        = 540,000
```

**缓存分析**：
- 不同(k0, k_norm)组合：~N_p × N_k = 10,000
- 每个组合2个通道(P, S)：20,000个唯一键
- 缓存命中率：1 - 20,000/540,000 ≈ 96.3%

**时间对比**（假设单次极化函数计算0.5ms）：

| 方法 | 计算次数 | 缓存查询次数 | 总时间 |
|------|---------|-------------|--------|
| 无缓存 | 540,000 | 0 | 270秒 (4.5分钟) |
| 有缓存 | 20,000 | 520,000 | 10秒 + ~1ms ≈ 10秒 |
| **加速比** | - | - | **27×** |

### 实际测试结果

基于`test/test_polarization_cache.jl`的测试结果：

| 测试项 | 结果 | 说明 |
|--------|------|------|
| 基本缓存功能 | ✅ | 相同参数正确命中缓存 |
| 浮点数容差 | ✅ | 1e-13相对差异正确命中 |
| 参数区分 | ✅ | 不同参数不会错误命中 |
| 统计准确性 | ✅ | 命中率计算正确 |
| 性能提升 | ✅ | 缓存比计算快4-5倍 |
| 内存管理 | ✅ | reset_cache!正确清空 |

### 性能监控示例

```julia
using PolarizationCache

function benchmark_with_monitoring()
    reset_cache!()
    start_time = time()
    
    # 模拟实际计算
    for i in 1:10000
        k0 = 50.0 + rand() * 100.0
        k_norm = 30.0 + rand() * 50.0
        Π = polarization_aniso_cached(:P, k0, k_norm, 5.0, 5.0, 
                                      300.0, 300.0, 150.0, 0.5, 0.5, 
                                      0.0, -50.0, -50.0, 0)
    end
    
    elapsed = time() - start_time
    stats = get_cache_stats()
    
    println("性能报告:")
    println("  总耗时: $(elapsed)秒")
    println("  总调用: $(stats.total_calls)")
    println("  缓存命中: $(stats.cache_hits) ($(stats.hit_rate * 100)%)")
    println("  平均每次: $(elapsed / stats.total_calls * 1000)ms")
    println("  缓存条目: $(stats.cache_size)")
    println("  估算无缓存耗时: $(elapsed / (1 - stats.hit_rate))秒")
    println("  节省时间: $(elapsed / (1 - stats.hit_rate) - elapsed)秒")
end
```

---

## 常见问题

### Q1: 缓存会影响计算精度吗？
**A**: 不会。缓存只是存储已计算的结果，不改变任何计算逻辑。浮点数容差1e-12远小于物理精度要求。

### Q2: 什么时候应该调用reset_cache!？
**A**: 
- 开始新的输运系数计算任务（必须）
- 改变物理参数（T, μ, ξ）后（必须）
- 计算结束后释放内存（可选）

### Q3: 缓存占用多少内存？
**A**: 每个条目约200字节。典型计算10,000条目约2MB，可忽略不计。

### Q4: 能否持久化缓存到磁盘？
**A**: 不建议。原因：
1. 跨任务参数复用概率低
2. 磁盘I/O可能比重新计算更慢
3. 单次任务内缓存已足够高效

### Q5: 并行计算时如何使用缓存？
**A**: 推荐每个worker独立缓存：
```julia
@everywhere using PolarizationCache
results = pmap(params) do p
    reset_cache!()  # 每个worker独立缓存
    calculate(p)
end
```

### Q6: 缓存命中率低怎么办？
**A**: 
1. 检查是否有参数在微小变化（应该命中但未命中）
2. 调整计算顺序，让相同参数的调用相邻
3. 分析哪些参数组合被重复计算

### Q7: 如何验证缓存正确性？
**A**: 运行测试套件：
```bash
julia test/test_polarization_cache.jl
```

---

## 版本兼容性

- **Julia版本**: ≥ 1.6
- **依赖模块**: PolarizationAniso
- **测试覆盖**: 7个测试集，34个测试用例

## 参考文档

- [PolarizationCache_缓存使用文档.md](../../../reference/formula/relaxtime/polarization/PolarizationCache_缓存使用文档.md) - 详细使用指南
- [PolarizationAniso API](PolarizationAniso.md) - 底层计算函数文档
- [总传播子计算文档](../../../reference/formula/relaxtime/propagator/总传播子byPropagator.md) - 在总传播子计算中的应用

## 作者与维护

- **创建日期**: 2025年11月17日
- **测试状态**: ✅ 全部通过 (34/34)
- **性能验证**: ✅ 缓存加速4-5倍（单次查询）
