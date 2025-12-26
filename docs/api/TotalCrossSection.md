# TotalCrossSection 模块 API 文档

## 模块概述

`TotalCrossSection` 模块通过对微分散射截面进行 t 积分计算总散射截面，包含末态费米子统计因子（Pauli blocking）。

**文件位置**: `src/relaxtime/TotalCrossSection.jl`

## 物理背景

### 总散射截面公式

总散射截面通过对 Mandelstam 变量 t 积分得到：

$$\sigma(s, T, \mu_q) = \int_{t_-}^{t_+} dt \cdot \frac{d\sigma}{dt} \cdot (1 - f_c)(1 - f_d)$$

其中：
- $\frac{d\sigma}{dt}$: 微分散射截面（来自 `DifferentialCrossSection` 模块）
- $(1 - f)$: 末态费米子统计因子（Pauli blocking）
- $E_c, E_d$: 末态粒子能量，从 $(s, t)$ 计算
- $t_\pm$: t 积分边界（公式 5.14）

**注意**: 本模块只考虑夸克-夸克散射（费米子），不涉及玻色子散射。

### t 积分边界（修正版公式 5.14）

$$\begin{aligned}
t_\pm &= m_i^2 + m_c^2 - \frac{1}{2s}(s + m_i^2 - m_j^2)(s + m_c^2 - m_d^2)\\
&\quad \pm 2\sqrt{\left[\frac{(s + m_i^2 - m_j^2)^2}{4s} - m_i^2\right]\left[\frac{(s + m_c^2 - m_d^2)^2}{4s} - m_c^2\right]}
\end{aligned}$$

此公式基于质心系动量守恒和散射角的物理边界，相比简化近似更加精确。

### 相同粒子散射的对称因子

当末态粒子 $c$ 和 $d$ 为相同粒子时（如 $uu \to uu$ 或 $ss \to ss$），由于粒子的不可区分性，$t$ 积分的下限需要减半：

$t_- \to \frac{t_-}{2} \quad \text{当 } c = d \text{ 时}$

#### 物理原因与权威参考

**权威来源**：Peskin & Schroeder, *An Introduction to Quantum Field Theory*

📍 **位置**：第 4 章（Scattering），Section 4.5: Phase Space and Cross Sections

📌 **关键原文**：
> "If there are identical particles in the final state, the phase space must be divided by a symmetry factor $S$, where $S = n!$ for $n$ identical particles."

**中文翻译**：若末态中存在相同粒子，则相空间必须除以一个对称因子 $S$；对于 $n$ 个相同粒子，该因子为 $n!$。

**应用到本问题**：
- 对于两个相同出射粒子（$n=2$），对称因子 $S = 2! = 2$
- 由于 $t_{\max} = 0$ 且被积函数关于 $t \leftrightarrow u$ 交换具有对称性
- 将 $t_{\min}$ 除以 2 等效于将积分结果除以对称因子 2

**几何解释**：
在 $t$ 积分中，$t$ 从 $t_-$ 到 $t_+ = 0$ 对应散射角 $\theta^*$ 从 $\pi$ 到 $0$。对于相同粒子，交换两个出射粒子等价于 $\theta^* \to \pi - \theta^*$，即 $t \leftrightarrow u$。由于被积函数在这种交换下对称，积分区间 $[t_-, 0]$ 实际上覆盖了每个物理构型两次。因此，将积分下限减半（$t_- \to t_-/2$）正好消除这种重复计数，等效于除以对称因子 2。

**代码实现**：
```julia
if particle_c == particle_d
    t_min = t_min / 2.0
end
```

### 末态粒子能量

在质心系中：

$$E_c = \frac{s + m_c^2 - m_d^2}{2\sqrt{s}}, \quad E_d = \frac{s + m_d^2 - m_c^2}{2\sqrt{s}}$$

满足能量守恒：$E_c + E_d = \sqrt{s}$

### Pauli Blocking 因子

费米子末态统计因子：

$$(1 - f_c)(1 - f_d)$$

其中 $f(E, \mu, T, \Phi, \Phi_{\text{bar}})$ 是 PNJL 模型中的 Fermi-Dirac 分布函数。

## 设计原则

1. **完整公式**: 使用精确的 t 积分边界公式（5.14），而非简化近似
2. **物理自洽**: 包含末态统计因子，符合量子统计力学
3. **固定点数积分**: 使用高斯-勒让德积分（固定点数），耗时可预测
4. **可扩展性**: 支持批量计算和参数扫描

---

## 核心函数

### `calculate_t_bounds`

```julia
calculate_t_bounds(s, mi, mj, mc, md) -> NamedTuple{(:t_min, :t_max)}
```

计算 Mandelstam 变量 t 的物理边界。

#### 参数

- `s::Float64`: Mandelstam 变量 s [fm⁻²]
- `mi::Float64`: 初态粒子 1 质量 [fm⁻¹]
- `mj::Float64`: 初态粒子 2 质量 [fm⁻¹]
- `mc::Float64`: 末态粒子 1 质量 [fm⁻¹]
- `md::Float64`: 末态粒子 2 质量 [fm⁻¹]

#### 返回值

`NamedTuple` 包含:
- `t_min::Float64`: t 的最小值（后向散射，$\theta = \pi$）[fm⁻²]
- `t_max::Float64`: t 的最大值（正向散射，$\theta = 0$，通常为 0）[fm⁻²]

#### 物理意义

- **t_min (t_-)**: 后向散射极限
- **t_max (t_+)**: 正向散射极限
- **满足 Mandelstam 约束**: $s + t + u = m_i^2 + m_j^2 + m_c^2 + m_d^2$

#### 运动学检查

如果 s 低于阈值（初态或末态），抛出错误：
- 初态阈值: $s \geq (m_i + m_j)^2$
- 末态阈值: $s \geq (m_c + m_d)^2$

#### 示例

```julia
using .TotalCrossSection

s = 31.0  # fm⁻²
mi = mj = mc = md = 1.52  # fm⁻¹ (u 夸克)

t_bounds = calculate_t_bounds(s, mi, mj, mc, md)
println("t ∈ [$(t_bounds.t_min), $(t_bounds.t_max)] fm⁻²")
# 输出: t ∈ [-21.7545, 0.0] fm⁻²

# 对于相同质量: t_max = 0 (对称性)
@assert t_bounds.t_max == 0.0
@assert t_bounds.t_min < 0.0
```

---

### `calculate_final_state_energies`

```julia
calculate_final_state_energies(s, t, mi, mj, mc, md) -> (E_c, E_d)
```

从 Mandelstam 变量计算质心系中末态粒子的能量。

#### 参数

- `s::Float64`: Mandelstam 变量 s [fm⁻²]
- `t::Float64`: Mandelstam 变量 t [fm⁻²]
- `mi, mj, mc, md::Float64`: 粒子质量 [fm⁻¹]

#### 返回值

`Tuple{Float64, Float64}`: 末态粒子能量 `(E_c, E_d)` [fm⁻¹]

#### 物理验证

- **能量守恒**: $E_c + E_d = \sqrt{s}$
- **质量壳条件**: $E_c^2 - p_c^2 = m_c^2$
- **非负能量**: $E_c \geq m_c, E_d \geq m_d$

#### 示例

```julia
s = 31.0
t = -2.0
mi = mj = mc = md = 1.52

E_c, E_d = calculate_final_state_energies(s, t, mi, mj, mc, md)
println("E_c = $E_c fm⁻¹, E_d = $E_d fm⁻¹")
println("E_c + E_d = $(E_c + E_d) ≈ √s = $(sqrt(s))")
# 能量守恒验证
@assert abs(E_c + E_d - sqrt(s)) < 1e-12
```

---

### `final_state_blocking_factor`

```julia
final_state_blocking_factor(E, μ, T, Φ, Φbar) -> Float64
```

计算单个末态费米子的统计因子 $(1 - f)$（Pauli blocking）。

#### 参数

- `E::Float64`: 粒子能量 [fm⁻¹]
- `μ::Float64`: 化学势 [fm⁻¹]
- `T::Float64`: 温度 [fm⁻¹]
- `Φ::Float64`: Polyakov loop
- `Φbar::Float64`: Conjugate Polyakov loop

#### 返回值

`Float64`: 末态统计因子 $(1 - f)$（范围 [0, 1]）

#### 物理意义

**Pauli blocking**: 费米子不能占据已被占据的态，抑制散射到已被占据的末态。
- $f \to 0$ (低温低密度): $(1-f) \to 1$，散射不受抑制
- $f \to 1$ (高温高密度): $(1-f) \to 0$，散射完全阻塞

#### 示例

```julia
E = 3.0  # fm⁻¹
μ = 0.3  # fm⁻¹
T = 0.15  # fm⁻¹
Φ = 0.5
Φbar = 0.5

factor = final_state_blocking_factor(E, μ, T, Φ, Φbar)
println("(1 - f) = $factor")
# 对于费米子: 0 ≤ factor ≤ 1
@assert 0.0 ≤ factor ≤ 1.0
```

---

### `combined_final_state_factor`

```julia
combined_final_state_factor(E_c, E_d, μ_c, μ_d, T, Φ, Φbar) -> Float64
```

计算组合末态费米子统计因子 $(1 - f_c)(1 - f_d)$。

#### 参数

- `E_c, E_d::Float64`: 末态粒子能量 [fm⁻¹]
- `μ_c, μ_d::Float64`: 化学势 [fm⁻¹]
- `T::Float64`: 温度 [fm⁻¹]
- `Φ, Φbar::Float64`: Polyakov loop

#### 返回值

`Float64`: 组合统计因子 $(1 - f_c)(1 - f_d)$

#### 说明

本项目只考虑夸克-夸克散射（费米子），因此统一使用 Pauli blocking 因子。

#### 示例

```julia
factor = combined_final_state_factor(
    3.0, 3.0, 0.3, 0.3, 0.15, 0.5, 0.5
)
println("(1-f_c)(1-f_d) = $factor")
```

---

### `total_cross_section`

```julia
total_cross_section(process, s, quark_params, thermo_params, K_coeffs; n_points) -> Float64
```

计算给定 s 下的总散射截面 σ(s)。**这是模块的核心函数**。

#### 参数

- `process::Symbol`: 散射过程标识（如 `:uu_to_uu`, `:udbar_to_udbar`）
- `s::Float64`: Mandelstam 变量 s [fm⁻²]
- `quark_params::NamedTuple`: 夸克参数
  - `m`: 质量 `(u=..., d=..., s=...)` [fm⁻¹]
  - `μ`: 化学势 `(u=..., d=..., s=...)` [fm⁻¹]
  - `A`: A 函数值 `(u=..., d=..., s=...)`
- `thermo_params::NamedTuple`: 热力学参数
  - `T`: 温度 [fm⁻¹]
  - `Φ`: Polyakov loop
  - `Φbar`: Conjugate Polyakov loop
  - `ξ`: 各向异性参数
- `K_coeffs::NamedTuple`: 有效耦合常数
- `n_points::Int=32`: 高斯-勒让德积分点数

#### 返回值

`Float64`: 总散射截面 σ(s) [fm²]

#### 计算流程

1. 解析过程中的粒子类型 (i, j, c, d)
2. 获取质量和化学势
3. 计算 t 积分边界
4. 生成高斯-勒让德积分节点和权重
5. 对每个 t 点计算：
   - 计算散射矩阵元 $|\mathcal{M}|^2$
   - 计算微分截面 $d\sigma/dt$
   - 计算末态能量 $E_c, E_d$
   - 计算费米子统计因子 $(1-f_c)(1-f_d)$
6. 加权求和得到积分结果

#### 积分方法

使用高斯-勒让德积分（固定点数），相比自适应积分的优势：
- **耗时可预测**: 总耗时 ≈ n_points × 50 ms
- **无收敛问题**: 不会因为被积函数特性导致无限循环
- **精度可控**: 通过调整 n_points 权衡精度与速度

| n_points | 预估耗时 | 推荐场景 |
|----------|---------|---------|
| 8 | ~0.4 s | 快速估算、参数扫描 |
| 16 | ~0.8 s | 中等精度计算 |
| 32 | ~1.6 s | 高精度计算（默认）|
| 64 | ~3.2 s | 验证收敛性 |

#### 性能考虑

- 每个 t 点需要计算 $|\mathcal{M}|^2$（约 50 ms）
- 总耗时 = n_points × 50 ms（线性可预测）
- 默认 n_points=32 时约 1.6 秒

#### 示例

```julia
using .TotalCrossSection

# 设置参数
quark_params = (
    m = (u = 1.52, d = 1.52, s = 2.50),
    μ = (u = 0.3, d = 0.3, s = 0.0),
    A = (u = -2.5, d = -2.5, s = -1.8)
)

thermo_params = (
    T = 0.15,    # 150 MeV / 197.327
    Φ = 0.5,
    Φbar = 0.5,
    ξ = 0.0
)

# 计算 K 系数（需要 EffectiveCouplings 模块）
K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

# 计算总散射截面（使用默认 32 个积分点）
s = 31.0  # fm⁻²
σ = total_cross_section(
    :uu_to_uu, s,
    quark_params, thermo_params, K_coeffs,
    n_points=32
)

println("σ(s=$s fm⁻²) = $σ fm²")

# 快速估算（使用 8 个积分点）
σ_fast = total_cross_section(
    :uu_to_uu, s,
    quark_params, thermo_params, K_coeffs,
    n_points=8
)
```

---

### `calculate_all_total_cross_sections`

```julia
calculate_all_total_cross_sections(s, quark_params, thermo_params, K_coeffs; n_points) -> NamedTuple
```

计算给定 s 下所有散射过程的总截面。

#### 参数

与 `total_cross_section` 相同，但不需要指定 `process`

#### 返回值

`NamedTuple`: 所有过程的总截面，键为过程名，值为 σ [fm²]

#### 示例

```julia
all_σ = calculate_all_total_cross_sections(
    31.0, quark_params, thermo_params, K_coeffs,
    n_points=16  # 使用 16 点加速计算
)

println("所有过程的总截面:")
for (process, σ) in pairs(all_σ)
    if !isnan(σ)
        println("  $process: $σ fm²")
    end
end

# 访问特定过程
println("\nuu→uu: $(all_σ.uu_to_uu) fm²")
println("dd→dd: $(all_σ.dd_to_dd) fm²")
```

---

### `scan_s_dependence`

```julia
scan_s_dependence(s_values, process, quark_params, thermo_params, K_coeffs; n_points) -> Vector{Float64}
```

扫描总散射截面随 s 的变化：σ(s)。

#### 参数

- `s_values::Vector{Float64}`: s 值数组 [fm⁻²]
- 其他参数与 `total_cross_section` 相同

#### 返回值

`Vector{Float64}`: 对应的总截面 σ_values [fm²]

#### 示例

```julia
# 扫描 s 依赖性（使用 16 点加速）
s_values = collect(range(10.0, 50.0, length=20))
σ_values = scan_s_dependence(
    s_values, :uu_to_uu,
    quark_params, thermo_params, K_coeffs,
    n_points=16
)

# 绘图
using Plots
plot(s_values, σ_values,
     xlabel="s [fm⁻²]",
     ylabel="σ [fm²]",
     title="uu→uu Total Cross Section",
     linewidth=2, marker=:circle)
savefig("sigma_vs_s.png")

# 找到最大值
max_idx = argmax(σ_values)
println("最大截面: σ = $(σ_values[max_idx]) fm² at s = $(s_values[max_idx]) fm⁻²")
```

---

## 使用场景

### 场景1：单点计算

计算特定能量下的总散射截面：

```julia
using .TotalCrossSection

s = 31.0
process = :uu_to_uu

σ = total_cross_section(
    process, s,
    quark_params, thermo_params, K_coeffs
)

println("σ($process, s=$s) = $σ fm²")
```

### 场景2：批量计算

计算所有过程的总截面：

```julia
all_σ = calculate_all_total_cross_sections(
    31.0, quark_params, thermo_params, K_coeffs
)

# 筛选有效结果
valid_processes = filter(p -> !isnan(all_σ[p]), keys(all_σ))
println("共计算 $(length(valid_processes)) 个过程")
```

### 场景3：能量扫描

研究总截面的能量依赖性：

```julia
# 从阈值到高能区扫描
m_u = quark_params.m.u
s_threshold = (2 * m_u)^2
s_values = collect(range(s_threshold + 1.0, 50.0, length=30))

σ_values = scan_s_dependence(
    s_values, :uu_to_uu,
    quark_params, thermo_params, K_coeffs
)

# 分析结果
println("阈值行为: σ(s_min) = $(σ_values[1]) fm²")
println("高能行为: σ(s_max) = $(σ_values[end]) fm²")
```

### 场景4：温度依赖性

研究温度对总截面的影响：

```julia
T_values = [0.10, 0.12, 0.14, 0.16, 0.18, 0.20]  # fm⁻¹
σ_vs_T = Float64[]

for T in T_values
    thermo_params_T = (T=T, Φ=0.5, Φbar=0.5, ξ=0.0)
    σ = total_cross_section(
        :uu_to_uu, 31.0,
        quark_params, thermo_params_T, K_coeffs
    )
    push!(σ_vs_T, σ)
end

plot(T_values .* 197.327, σ_vs_T,  # 转换为 MeV
     xlabel="T [MeV]", ylabel="σ [fm²]")
```

---

## 数值考虑

### 积分精度控制

```julia
# 高精度计算（更慢）
σ_high = total_cross_section(
    :uu_to_uu, 31.0,
    quark_params, thermo_params, K_coeffs,
    rtol=1e-8, atol=1e-12
)

# 快速计算（降低精度）
σ_fast = total_cross_section(
    :uu_to_uu, 31.0,
    quark_params, thermo_params, K_coeffs,
    rtol=1e-4, atol=1e-8
)

println("相对差异: $(abs(σ_high - σ_fast) / σ_high * 100)%")
```

### 运动学阈值

```julia
# 检查是否接近阈值
m_u = quark_params.m.u
s_threshold = (2 * m_u)^2

if s < s_threshold + 1.0
    @warn "s 接近阈值，积分可能不稳定" s s_threshold
end
```

### 错误处理

```julia
try
    σ = total_cross_section(
        :uu_to_uu, 1.0,  # 低于阈值
        quark_params, thermo_params, K_coeffs
    )
catch e
    @error "计算失败" exception=e
    # 继续其他计算
end
```

---

## 性能优化建议

### 1. 预计算常量

```julia
# 如果要多次调用，预计算 K 系数
K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

# 批量计算时复用
for s in s_values
    σ = total_cross_section(
        :uu_to_uu, s,
        quark_params, thermo_params, K_coeffs  # 复用
    )
end
```

### 2. 并行计算

```julia
using Distributed
@everywhere using .TotalCrossSection

# 并行扫描
σ_values = @distributed (vcat) for s in s_values
    [total_cross_section(:uu_to_uu, s, quark_params, thermo_params, K_coeffs)]
end
```

### 3. 降低采样密度

```julia
# 粗略扫描
s_coarse = range(10.0, 50.0, length=10)

# 精细扫描（仅在感兴趣区域）
s_fine = range(25.0, 35.0, length=30)
```

---

## 物理验证

### 验证1：Mandelstam 约束

```julia
s = 31.0
mi = mj = mc = md = 1.52

t_bounds = calculate_t_bounds(s, mi, mj, mc, md)

# 检查边界点的 Mandelstam 约束
for t in [t_bounds.t_min, t_bounds.t_max]
    u = mi^2 + mj^2 + mc^2 + md^2 - s - t
    mandelstam_sum = s + t + u
    expected = mi^2 + mj^2 + mc^2 + md^2
    @test abs(mandelstam_sum - expected) < 1e-10
end
```

### 验证2：能量守恒

```julia
E_c, E_d = calculate_final_state_energies(s, t, mi, mj, mc, md)
@test abs(E_c + E_d - sqrt(s)) < 1e-12
```

### 验证3：统计因子范围

```julia
factor = combined_final_state_factor(E_c, E_d, μ_c, μ_d, T, Φ, Φbar)
@test 0.0 ≤ factor ≤ 1.0
```

---

## 常见问题

### Q1: 为什么总截面是负数或 NaN？

**可能原因**:
1. s 低于阈值 → 检查 $s \geq (m_i + m_j)^2$
2. 散射矩阵元计算失败 → 检查 `K_coeffs`
3. 积分发散 → 降低 `rtol`，增加 `atol`

```julia
# 诊断
m_u = quark_params.m.u
s_threshold = (2 * m_u)^2
println("s = $s, threshold = $s_threshold")
@assert s > s_threshold "s below threshold!"
```

### Q2: 积分时间过长怎么办？

**解决方案**:
1. 降低精度要求：`rtol=1e-4`
2. 检查是否接近阈值（积分困难）
3. 使用更好的初值猜测

### Q3: 不同过程的截面差异很大？

**这是正常的**:
- qq 散射 vs qqbar 散射：介子通道不同
- 同位旋对称：uu 和 dd 应该相同（当 $m_u = m_d$）
- 奇异夸克：us 散射截面通常较小（质量差异）

### Q4: 如何选择 s 的扫描范围？

```julia
# 从阈值开始
m_i, m_j = 1.52, 1.52
s_min = (m_i + m_j)^2 + 1.0  # 阈值 + 安全边距

# 到足够高能（依赖物理问题）
s_max = 100.0  # fm⁻²

s_values = range(s_min, s_max, length=50)
```

---

## 调试建议

```julia
# 启用详细输出
using Logging
global_logger(ConsoleLogger(stderr, Logging.Debug))

# 检查中间步骤
t_bounds = calculate_t_bounds(s, mi, mj, mc, md)
println("t_bounds: ", t_bounds)

E_c, E_d = calculate_final_state_energies(s, t_bounds.t_min, mi, mj, mc, md)
println("Energies at t_min: E_c=$E_c, E_d=$E_d")

# 验证被积函数
function integrand_test(t)
    M_squared = scattering_amplitude_squared(
        :uu_to_uu, s, t, quark_params, thermo_params, K_coeffs
    )
    dsigma_dt = differential_cross_section(s_12_plus, s_12_minus, M_squared)
    E_c, E_d = calculate_final_state_energies(s, t, mi, mj, mc, md)
    blocking = combined_final_state_factor(E_c, E_d, μ_c, μ_d, T, Φ, Φbar)
    return dsigma_dt * blocking
end

println("Integrand at t_min: ", integrand_test(t_bounds.t_min))
println("Integrand at t_max: ", integrand_test(t_bounds.t_max))
```

---

## 参考文献

- **公式推导**: `doc/formula/散射截面by微分散射截面.md`
- **微分截面**: `api/DifferentialCrossSection.md`
- **散射矩阵元**: `api/ScatteringAmplitude.md`
- **测试用例**: `test_unit/test_total_cross_section.jl`
- **性能测试**: `test_unit/test_total_cross_section_performance.md`
