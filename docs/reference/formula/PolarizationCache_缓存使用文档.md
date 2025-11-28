# 极化函数缓存 - 使用文档

## 1. 模块标识
- **模块名称**: `PolarizationCache`
- **对应文件**: `src/relaxtime/PolarizationCache.jl`

## 2. 功能说明

### 2.1 核心功能
通过哈希表缓存极化函数计算结果，避免在单次输运系数计算中重复计算相同参数的极化函数。

### 2.2 性能优势
- **计算成本**：单次极化函数计算需要高斯积分，耗时约0.1-1ms
- **缓存成本**：哈希表查询耗时约1ns
- **加速比**：缓存命中可获得10⁵-10⁶倍加速
- **典型收益**：在输运系数计算中，缓存命中率可达30%-70%

### 2.3 使用场景
在以下情况下，极化函数会被重复计算相同参数：
- 同一散射过程的t/u道共享相同的介子动量
- 不同散射过程（如uu→uu和dd→dd）共享Π_uu
- 相同(k0, k_norm)但不同通道(:P vs :S)

## 3. API 接口

### 3.1 带缓存的极化函数计算

```julia
polarization_aniso_cached(channel, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1_value, A2_value, num_s_quark)
```

**参数**：与`PolarizationAniso.polarization_aniso`完全相同

**返回值**：`(Π_real, Π_imag)` - 复数极化函数的实部和虚部

**行为**：
1. 查询缓存：如果参数匹配（考虑浮点数容差），直接返回缓存结果
2. 计算并缓存：如果缓存未命中，调用底层函数计算并存储结果

### 3.2 缓存管理

```julia
# 清空缓存和统计信息
reset_cache!()

# 获取缓存统计
stats = get_cache_stats()
# 返回 NamedTuple:
#   - total_calls: 总调用次数
#   - cache_hits: 缓存命中次数
#   - cache_misses: 缓存未命中次数
#   - hit_rate: 缓存命中率（0.0-1.0）
#   - cache_size: 当前缓存条目数
```

## 4. 使用示例

### 4.1 基本使用

```julia
using PolarizationCache

# 在输运系数计算开始前重置缓存
reset_cache!()

# 正常调用（会自动缓存）
Π_real, Π_imag = polarization_aniso_cached(
    :P,          # 赝标量通道
    100.0,       # k0 (MeV)
    50.0,        # k_norm (MeV)
    5.0, 5.0,    # m1, m2 (MeV)
    300.0, 300.0,# μ1, μ2 (MeV)
    150.0,       # T (MeV)
    0.5, 0.5,    # Φ, Φbar
    0.0,         # ξ
    -50.0, -50.0,# A1, A2
    0            # num_s_quark
)

# 重复调用相同参数会从缓存读取
Π_real2, Π_imag2 = polarization_aniso_cached(:P, 100.0, 50.0, ...) # 快~10⁵倍

# 查看缓存统计
stats = get_cache_stats()
println("缓存命中率: $(stats.hit_rate * 100)%")
```

### 4.2 输运系数计算中的使用

```julia
using PolarizationCache

function calculate_transport_coefficient(T, μ, ξ, ...)
    # 开始新计算任务
    reset_cache!()
    
    # 预计算A函数（所有极化函数共享）
    A_u = calculate_A_function(m_u, μ_u, T, Φ, Φbar)
    A_s = calculate_A_function(m_s, μ_s, T, Φ, Φbar)
    
    # 循环计算多个散射过程和散射道
    for process in [:uu_to_uu, :ud_to_ud, :us_to_us, ...]
        for channel in [:t, :u, :s]
            # 计算介子动量
            k0, k_norm = calculate_meson_momentum(process, channel, ...)
            
            # 计算极化函数（自动缓存）
            Π_uu = polarization_aniso_cached(:P, k0, k_norm, m_u, m_u, μ_u, μ_u, T, Φ, Φbar, ξ, A_u, A_u, 0)
            Π_us = polarization_aniso_cached(:P, k0, k_norm, m_u, m_s, μ_u, μ_s, T, Φ, Φbar, ξ, A_u, A_s, 1)
            
            # 计算传播子和散射振幅
            D = meson_propagator(...)
            M = scattering_amplitude(...)
            
            # ... 累积贡献 ...
        end
    end
    
    # 输出缓存统计
    stats = get_cache_stats()
    println("极化函数调用: $(stats.total_calls)次")
    println("缓存命中: $(stats.cache_hits)次 ($(stats.hit_rate * 100)%)")
    println("实际计算: $(stats.cache_misses)次")
    println("节省时间: 约$(stats.cache_hits * 0.5)ms")  # 假设单次0.5ms
    
    # 计算结束后清空缓存（释放内存）
    reset_cache!()
    
    return result
end
```

### 4.3 批量计算时的缓存效率

```julia
using PolarizationCache

# 模拟批量计算不同温度下的输运系数
temperatures = [100.0, 120.0, 140.0, 160.0, 180.0]

for T in temperatures
    println("\n计算 T = $T MeV")
    reset_cache!()  # 每个温度独立缓存
    
    # 模拟计算100个不同动量点
    for i in 1:100
        k0 = 50.0 + i * 2.0
        k_norm = 30.0 + i * 1.0
        
        # P通道和S通道都需要计算
        Π_P = polarization_aniso_cached(:P, k0, k_norm, ...)
        Π_S = polarization_aniso_cached(:S, k0, k_norm, ...)
    end
    
    stats = get_cache_stats()
    println("  总调用: $(stats.total_calls)")
    println("  命中率: $(stats.hit_rate * 100)%")
    println("  缓存大小: $(stats.cache_size) 条目")
end
```

## 5. 技术细节

### 5.1 缓存键设计
使用`PolarizationKey`结构体作为哈希表键，包含所有影响计算结果的14个参数：
- `channel`：通道类型
- `k0, k_norm`：四动量
- `m1, m2, μ1, μ2`：夸克质量和化学势
- `T, Φ, Φbar, ξ`：热力学参数
- `A1, A2`：预计算的A函数值
- `num_s_quark`：奇异夸克数量

### 5.2 浮点数容差
- 使用相对容差 `EPS_CACHE = 1e-12` 比较浮点数参数
- 哈希时保留12位有效数字避免舍入误差导致缓存失效
- 判等时：`|a - b| / max(|a|, |b|) < 1e-12`

### 5.3 线程安全
**注意**：当前实现**不是线程安全**的。如果需要并行计算：
- 方案1：每个线程维护独立缓存（推荐）
- 方案2：使用线程安全的Dict（需修改代码）

### 5.4 内存管理
- 缓存条目数取决于计算的参数组合数
- 典型单次计算：10³-10⁴个条目，内存占用~1-10MB
- 建议：每次计算任务结束后调用`reset_cache!()`释放内存

## 6. 性能预期

### 6.1 理论分析
假设输运系数计算需要：
- 散射过程数：9种（uu, ud, us, dd, ds, ss, ...)
- 散射道：3种（t, u, s）
- 动量积分点：100×100网格
- 介子通道：2种（P, S）

**不使用缓存**：
- 极化函数调用：9 × 3 × 100 × 100 × 2 = 540,000次
- 计算时间：540,000 × 0.5ms = 270秒 = 4.5分钟

**使用缓存（假设50%命中率）**：
- 实际计算：270,000次
- 缓存读取：270,000次（几乎不耗时）
- 计算时间：约135秒 = 2.25分钟
- **节省时间：50%**

### 6.2 实际测试建议
在实际代码中添加性能监控：

```julia
# 开始计算
reset_cache!()
start_time = time()

# ... 进行输运系数计算 ...

# 结束统计
elapsed = time() - start_time
stats = get_cache_stats()

println("总耗时: $(elapsed)秒")
println("极化函数调用: $(stats.total_calls)次")
println("缓存命中率: $(stats.hit_rate * 100)%")
println("平均每次: $(elapsed / stats.total_calls * 1000)ms")
println("预估无缓存耗时: $(elapsed / (1 - stats.hit_rate))秒")
println("节省时间: $(elapsed / (1 - stats.hit_rate) - elapsed)秒")
```

## 7. 常见问题

### Q1: 缓存会影响计算精度吗？
不会。缓存只是存储已计算的结果，不改变任何计算逻辑。浮点数容差`1e-12`远小于物理精度要求。

### Q2: 何时需要调用reset_cache!？
- 开始新的输运系数计算任务
- 改变物理参数（温度、化学势等）后
- 内存不足时

### Q3: 缓存占用多少内存？
每个缓存条目约200字节（14个Float64 + 1个Symbol + 1个Int + 2个结果Float64）。
10,000个条目约2MB，通常可忽略。

### Q4: 能否持久化缓存到磁盘？
不建议。原因：
1. 极化函数依赖14个参数，跨任务复用概率低
2. 磁盘I/O可能比重新计算更慢
3. 单次计算任务内的缓存已足够高效

### Q5: 并行计算时如何使用缓存？
推荐方案：每个线程/进程独立缓存

```julia
using Distributed

@everywhere using PolarizationCache

# 每个worker独立缓存
results = pmap(temperatures) do T
    reset_cache!()  # 线程本地缓存
    result = calculate_at_temperature(T)
    stats = get_cache_stats()
    return (result, stats)
end
```

## 8. 与其他模块的集成

### 8.1 替换现有代码
只需将：
```julia
using .PolarizationAniso: polarization_aniso
Π = polarization_aniso(...)
```

改为：
```julia
using .PolarizationCache: polarization_aniso_cached, reset_cache!, get_cache_stats
Π = polarization_aniso_cached(...)
```

### 8.2 向后兼容
`PolarizationCache.jl`内部调用`PolarizationAniso.jl`，完全兼容原有接口和功能。

### 8.3 测试验证
运行测试套件验证缓存功能：
```bash
julia test/test_polarization_cache.jl
```

预期输出：
- ✅ 基本缓存功能
- ✅ 浮点数容差
- ✅ 不同参数区分
- ✅ 缓存统计
- ✅ 性能对比（>100x加速）
- ✅ 重置功能
- ✅ 实际使用场景模拟
