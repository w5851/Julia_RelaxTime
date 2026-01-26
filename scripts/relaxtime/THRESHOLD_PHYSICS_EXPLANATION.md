# 阈值物理的正确理解

**日期**: 2026-01-26
**问题**: 为什么用 max 而不是 min？上界在哪里？

---

## 问题1：为什么是 max 而不是 min？

### 代码

```julia
s_threshold_initial = (mi + mj)^2
s_threshold_final = (mc + md)^2
s_threshold = max(s_threshold_initial, s_threshold_final)
if s < s_threshold - EPS_KINEMATIC
    return 0.0
end
```

### 物理解释

对于 2→2 散射过程 **i + j → c + d**：

**质心系能量守恒**：
```
√s = E_i + E_j = E_c + E_d
```

**两个独立的约束**：

1. **初态粒子必须能够存在**：
   ```
   √s ≥ m_i + m_j  →  s ≥ (m_i + m_j)²
   ```
   - 这是初态的静止质量
   - 如果 s < (m_i + m_j)²，初态粒子根本无法存在

2. **末态粒子必须能够产生**：
   ```
   √s ≥ m_c + m_d  →  s ≥ (m_c + m_d)²
   ```
   - 这是末态的静止质量
   - 如果 s < (m_c + m_d)²，末态粒子无法产生

**关键**：这两个条件必须**同时**满足！

因此：
```
s ≥ (m_i + m_j)²  AND  s ≥ (m_c + m_d)²
```

等价于：
```
s ≥ max((m_i + m_j)², (m_c + m_d)²)
```

### 为什么不是 min？

**反例：ssbar → uubar**

如果用 `min`：
```
s_threshold = min((2m_s)², (2m_u)²) = (2m_u)² = 265 MeV²
```

在 s = 1000 MeV² 时（大于265但小于171,113）：
- ✓ 能量足够产生 uubar（需要265 MeV²）
- ✗ **但初态 ssbar 根本不存在！**（需要171,113 MeV²）

**物理上不可能**：你不能从不存在的粒子产生其他粒子！

### 正确的理解

**阈值是下界**，但是**两个约束的下界**：

```
s_threshold = max(s_threshold_initial, s_threshold_final)
```

这确保了：
1. ✅ 初态粒子能够存在
2. ✅ 末态粒子能够产生
3. ✅ 能量守恒

---

## 问题2：质心能量的上界在哪里？

### 物理上的上界

在实际计算中，质心能量 s 确实有上界，但这个上界**不是来自运动学约束**，而是来自：

1. **动量截断 Λ**（如果使用）
2. **分布函数的衰减**
3. **实验/计算的能量范围**

### 代码中的处理

#### 1. 在 `calculate_t_bounds` 中

```julia
# 检查运动学可行性
if sqrt_arg1 < -EPS_KINEMATIC
    s_threshold_ij = (mi + mj)^2
    error("Initial state kinematic violation: s = $s < threshold = $s_threshold_ij")
end
if sqrt_arg2 < -EPS_KINEMATIC
    s_threshold_cd = (mc + md)^2
    error("Final state kinematic violation: s = $s < threshold = $s_threshold_cd")
end
```

**这里只检查下界**，因为：
- 运动学上，s 可以任意大（只要有足够能量）
- 没有运动学上的上界

#### 2. 在 `AverageScatteringRate` 中

```julia
# 在 design_w0cdf_s_grid 函数中
s_bo = max((mi + mj)^2, (mc + md)^2)  # 下界
s_up = if p_cutoff !== nothing
    min((sqrt(mi^2 + p_cutoff^2) + sqrt(mj^2 + p_cutoff^2))^2,
        (sqrt(mc^2 + p_cutoff^2) + sqrt(md^2 + p_cutoff^2))^2)
else
    Inf  # 无上界
end
```

**这里有上界处理**：
- 如果指定了动量截断 `p_cutoff`（如 Λ），则计算对应的 s 上界
- 如果没有截断，s_up = Inf（无上界）

#### 3. 在 σ(s) 缓存中

```julia
# 在 get_sigma 函数中
if s < cache.s_vals[1] || s > cache.s_vals[end]
    return 0.0  # 超出缓存范围返回0
end
```

**隐式上界**：
- σ(s) 缓存有一个有限的 s 范围
- 超出范围时返回 0
- 这是一个**实用的上界**，不是物理上的

### 物理意义

**运动学上**：
- ✅ 有下界：s ≥ max((m_i+m_j)², (m_c+m_d)²)
- ❌ 无上界：s 可以任意大

**实际计算中**：
- ✅ 有下界：运动学阈值
- ✅ 有上界：
  - 动量截断 Λ → s_max ≈ (2√(m²+Λ²))²
  - 分布函数衰减 → 高 s 贡献可忽略
  - 缓存范围 → 超出返回 0

---

## 总结

### 关于 max vs min

**正确**：`s_threshold = max(s_threshold_initial, s_threshold_final)`

**原因**：
- 初态和末态都必须满足各自的阈值条件
- 取 max 确保两个条件同时满足
- 取 min 会导致物理上不可能的情况

### 关于上界

**运动学上**：
- 没有上界（s 可以任意大）

**实际计算中**：
- 在 `design_w0cdf_s_grid` 中处理
- 如果有动量截断 Λ：s_up = min((2√(m²+Λ²))²)
- 如果无截断：s_up = Inf
- σ(s) 缓存超出范围返回 0

**代码位置**：
```julia
# AverageScatteringRate.jl, design_w0cdf_s_grid 函数
s_up = if p_cutoff !== nothing
    min((sqrt(mi^2 + p_cutoff^2) + sqrt(mj^2 + p_cutoff^2))^2,
        (sqrt(mc^2 + p_cutoff^2) + sqrt(md^2 + p_cutoff^2))^2)
else
    Inf
end

# 跳过超出上界的点
s >= s_up && continue
```

---

## 图示

```
运动学约束：

s_threshold_initial = (m_i + m_j)²  ←─ 初态阈值
s_threshold_final = (m_c + m_d)²    ←─ 末态阈值

s_threshold = max(初态, 末态)       ←─ 总阈值（下界）

         ↓
    [s_threshold, ∞)                ←─ 运动学允许的范围
         ↓
    [s_threshold, s_up]             ←─ 实际计算范围（如果有Λ）
         ↓
    [s_min_cache, s_max_cache]      ←─ σ(s)缓存范围
```

---

*文档时间: 2026-01-26*
*结论: max 是正确的，上界在 design_w0cdf_s_grid 中处理*

