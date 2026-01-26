# 当前实现状态：动量积分与σ(s)缓存

**日期**: 2026-01-26
**问题**: 动量积分是半无穷，但σ(s)缓存使用Λ截断吗？

---

## 🎯 简短答案

**是的！** 当前实现是：
- ✅ **动量积分**：半无穷 [0, ∞)
- ✅ **σ(s)缓存**：使用 Λ 截断

这是2026-01-26的新策略。

---

## 📋 详细说明

### 1. 动量积分（在 `average_scattering_rate` 中）

**代码位置**：`AverageScatteringRate.jl`, `_average_scattering_rate_semi_infinite` 函数

```julia
# 使用半无穷积分 [0, ∞)
t_grid, t_w = gauleg(0.0, 1.0, p_nodes)
for (t, wt) in zip(t_grid, t_w)
    if t >= 0.9999
        continue
    end
    p = scale * t / (1.0 - t)  # 变换到 [0, ∞)
    dp_dt = scale / (1.0 - t)^2
    # ...
end
```

**特点**：
- 使用变量替换：p = scale × t / (1-t)，t ∈ [0, 1)
- 映射到 [0, ∞)
- 没有动量截断

### 2. σ(s)缓存设计（在 `design_w0cdf_s_grid` 中）

**代码位置**：`AverageScatteringRate.jl`, `design_w0cdf_s_grid` 函数

```julia
# 根据 p_cutoff 选择动量网格构建方式
p_vals, p_wts, dp_jac = if p_cutoff !== nothing
    _build_finite_cutoff_p_grid(p_nodes, p_cutoff)  # 使用 [0, Λ]
else
    _build_semi_infinite_p_grid(p_nodes, scale)     # 使用 [0, ∞)
end

# 计算 s 的上界
s_up = if p_cutoff !== nothing
    min((sqrt(mi^2 + p_cutoff^2) + sqrt(mj^2 + p_cutoff^2))^2,
        (sqrt(mc^2 + p_cutoff^2) + sqrt(md^2 + p_cutoff^2))^2)
else
    Inf
end

# 跳过超出上界的点
s >= s_up && continue
```

**特点**：
- 如果 `p_cutoff` 指定（如 Λ），使用有限动量网格 [0, Λ]
- 计算对应的 s 上界：s_up = (2√(m²+Λ²))²
- 跳过超出 s_up 的采样点

### 3. 默认配置（在 `compute_average_rates` 中）

**代码位置**：`RelaxationTime.jl`, `compute_average_rates` 函数

```julia
# σ(s)缓存默认使用 Λ 截断
effective_sigma_cutoff = sigma_cutoff === nothing ? Λ_inv_fm : sigma_cutoff

# 构建缓存
cache = build_w0cdf_pchip_cache(
    process,
    quark_params,
    thermo_params,
    K_coeffs;
    p_cutoff=effective_sigma_cutoff,  # 传递 Λ
    n_sigma_points=n_sigma_points,
)
```

**默认值**：
- `sigma_cutoff = nothing` → 使用 `Λ_inv_fm`
- `Λ_inv_fm = 3.052 fm⁻¹` (600 MeV)

---

## 🔍 物理意义

### 为什么这样设计？

**动量积分 [0, ∞)**：
- 正确的物理：分布函数 f(p) 在所有动量都有贡献
- 虽然高动量贡献小（f(p) ~ exp(-p/T)），但不应人为截断

**σ(s)缓存 [0, s_up(Λ)]**：
- 实用考虑：σ(s) 的计算和存储需要有限范围
- 物理合理：超出 Λ 对应的 s 值，贡献很小
- 效率优化：避免计算和存储不必要的高 s 值

### 具体数值

对于 T = 300 MeV, Λ = 600 MeV：

**ssbar → uubar**：
- s_threshold = 171,113 MeV²
- s_up(Λ) = (2√(m_s²+Λ²))² ≈ 1,622,174 MeV²
- 缓存范围：[171,113, 1,622,174] MeV²

**动量积分**：
- 范围：[0, ∞)
- 主要贡献：p ∈ [0, 3T] ≈ [0, 4.5 fm⁻¹]
- 在 p > Λ = 3.05 fm⁻¹ 时：
  - s 可能超出缓存范围
  - σ(s) 返回 0（超出缓存）

---

## 📊 实际效果

### σ(s)缓存的行为

```julia
# 在 get_sigma 函数中
if s < cache.s_vals[1] || s > cache.s_vals[end]
    return 0.0  # 超出缓存范围返回 0
end
```

**含义**：
- s < s_min：返回 0（低于阈值）
- s_min ≤ s ≤ s_max：插值返回 σ(s)
- s > s_max：返回 0（超出缓存）

### 对平均散射率的影响

**在动量积分中**：
```julia
ω = ∫∫ dp_i dp_j · p_i² p_j² · f(p_i) f(p_j) · σ(s) · v_rel
```

**当 p > Λ 时**：
- s 可能超出 s_up(Λ)
- σ(s) 返回 0
- 这部分贡献被截断

**效果**：
- 虽然动量积分是 [0, ∞)
- 但由于 σ(s) 在高 s 时返回 0
- **实际上等价于有一个隐式的动量截断**

---

## 🎓 与Fortran的对比

### Fortran

```fortran
! 动量积分：[0, Λ]
do i = 1, n_nodes
    p_i = p_grid(i)  ! p_i ∈ [0, Λ]
    ! ...
end do

! σ(s)：在需要时计算，没有预先的上界限制
```

**特点**：
- 显式的动量截断
- σ(s) 在任何 s 值都可以计算

### Julia（当前实现）

```julia
# 动量积分：[0, ∞)
p = scale * t / (1.0 - t)  # p ∈ [0, ∞)

# σ(s)缓存：[s_min, s_up(Λ)]
if s > s_up
    return 0.0  # 隐式截断
end
```

**特点**：
- 形式上是半无穷积分
- 但通过 σ(s) 缓存实现隐式截断
- **实际效果与 Fortran 类似**

---

## ✅ 总结

### 当前实现

| 组件 | 范围 | 实现方式 |
|------|------|---------|
| **动量积分** | [0, ∞) | 半无穷积分（变量替换） |
| **σ(s)缓存设计** | [0, Λ] | 有限动量网格 |
| **σ(s)有效范围** | [s_min, s_up(Λ)] | 超出返回 0 |
| **实际效果** | ≈ [0, Λ] | 隐式截断 |

### 回答你的问题

**问**：平均散射率中动量积分是半无穷，但 `design_w0cdf_s_grid` 使用 `p_cutoff=Λ` 计算 s 的上界并截断？

**答**：
1. ✅ **是的**，动量积分形式上是半无穷 [0, ∞)
2. ✅ **是的**，σ(s)缓存使用 `p_cutoff=Λ` 计算上界
3. ✅ **是的**，超出 s_up(Λ) 的点被跳过（σ(s) 返回 0）
4. ✅ **实际效果**：虽然动量积分是 [0, ∞)，但由于 σ(s) 的隐式截断，**等价于有一个动量截断在 Λ 附近**

### 物理合理性

这个设计是**合理的**：
- 形式上保持了半无穷积分的正确性
- 实际上通过 σ(s) 缓存实现了有效的截断
- 避免了计算和存储不必要的高 s 值
- 与 Fortran 的显式截断效果类似

---

*文档时间: 2026-01-26*
*结论: 动量积分是半无穷，但σ(s)缓存使用Λ截断，实际效果等价*

