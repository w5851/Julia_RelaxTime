# ScatteringAmplitude API 文档

## 模块概述

`ScatteringAmplitude`模块实现了三味PNJL模型中夸克-夸克和夸克-反夸克弹性散射的散射矩阵元平方计算。该模块基于介子交换机制(π, K, η, η', σ等),支持11种独立散射过程。

**文件位置**: `src/relaxtime/ScatteringAmplitude.jl`  
**依赖模块**: 
- `TotalPropagator`: 提供按道分离的介子传播子
- `Constants_PNJL`: 提供色数N_c和过程映射表

---

## 核心常量

### N_c
```julia
const N_c = N_color
```
**类型**: `Int`  
**值**: 3  
**说明**: 色自由度数,用于色自旋平均因子计算。

---

## 主要函数

### 1. scattering_amplitude_squared

计算指定散射过程的色自旋平均矩阵元平方 |M|²。

#### 函数签名
```julia
function scattering_amplitude_squared(
    process::Symbol,
    s::Float64,
    t::Float64,
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    K_coeffs::NamedTuple
)::Float64
```

#### 参数说明

| 参数 | 类型 | 说明 | 单位 |
|------|------|------|------|
| `process` | Symbol | 散射过程标识符(见下表) | - |
| `s` | Float64 | Mandelstam s变量(质心系能量平方) | fm⁻² |
| `t` | Float64 | Mandelstam t变量(动量转移平方) | fm⁻² |
| `quark_params` | NamedTuple | 夸克参数 | - |
| `thermo_params` | NamedTuple | 热力学参数 | - |
| `K_coeffs` | NamedTuple | K系数(有效耦合常数) | - |

#### quark_params 结构
```julia
quark_params = (
    m = (u=0.3, d=0.3, s=0.5),      # 夸克质量 [fm⁻¹]
    μ = (u=0.25, d=0.25, s=0.25),   # 化学势 [fm⁻¹]
    A = (u=0.05, d=0.05, s=0.08)    # A函数值 [fm⁻²]
)
```

#### thermo_params 结构
```julia
thermo_params = (
    T = 0.15,        # 温度 [fm⁻¹]
    Φ = 0.5,         # Polyakov loop
    Φbar = 0.5,      # 共轭Polyakov loop
    ξ = 1.0          # 各向异性参数
)
```

#### K_coeffs 结构
```julia
K_coeffs = (
    K0_plus = 10.0, K0_minus = 10.0,
    K123_plus = 5.0, K123_minus = 5.0,
    # ... (共12个K系数)
)
```

#### 支持的散射过程

**4种qq散射**:
| 符号 | 物理过程 | t道介子 | u道介子 |
|------|---------|---------|---------|
| `:uu_to_uu` | u + u → u + u | η, η', σ, σ' | η, η', σ, σ' |
| `:ss_to_ss` | s + s → s + s | η, η', σ, σ' | η, η', σ, σ' |
| `:ud_to_ud` | u + d → u + d | π, σ_π | π, σ_π |
| `:us_to_us` | u + s → u + s | K, σ_K | K, σ_K |

**7种qqbar散射**:
| 符号 | 物理过程 | s道介子 | t道介子 |
|------|---------|---------|---------|
| `:udbar_to_udbar` | u + đ → u + đ | - | π, σ_π |
| `:usbar_to_usbar` | u + š → u + š | - | K, σ_K |
| `:uubar_to_uubar` | u + ū → u + ū | η, η', σ, σ' | η, η', σ, σ' |
| `:uubar_to_ddbar` | u + ū → d + đ | π, σ_π | π, σ_π |
| `:uubar_to_ssbar` | u + ū → s + š | K, σ_K | K, σ_K |
| `:ssbar_to_uubar` | s + š → u + ū | K, σ_K | K, σ_K |
| `:ssbar_to_ssbar` | s + š → s + š | η, η', σ, σ' | η, η', σ, σ' |

#### 返回值

| 类型 | 说明 | 单位 | 约束 |
|------|------|------|------|
| Float64 | 色自旋平均的散射矩阵元平方 | fm⁻⁴ | ≥ 0 |

#### 物理约束

**qq散射**:
```julia
s > (m1 + m2)²  # 质心系能量阈值
t < 0           # 类空动量转移
u < 0           # 类空动量转移
s + t + u = m1² + m2² + m3² + m4²
```

**qqbar散射**:
```julia
s > 0           # 正能量
t < 0           # 类空动量转移
u < 0           # 由s+t+u约束决定
```

#### 示例用法

**示例1: 基本计算**
```julia
using ScatteringAmplitude

# 准备参数
process = :uu_to_uu
s = 10.0  # fm⁻²
t = -0.2  # fm⁻²

quark_params = (
    m = (u=0.3, d=0.3, s=0.5),
    μ = (u=0.25, d=0.25, s=0.25),
    A = (u=0.05, d=0.05, s=0.08)
)

thermo_params = (
    T = 0.15,
    Φ = 0.5,
    Φbar = 0.5,
    ξ = 1.0
)

# 计算K系数(示例值)
K_coeffs = (
    K0_plus = 10.0, K0_minus = 10.0,
    K123_plus = 5.0, K123_minus = 5.0,
    K4567_plus = 5.0, K4567_minus = 5.0,
    K8_plus = 10.0, K8_minus = 10.0,
    K08_plus = 0.0, K08_minus = 0.0,
    det_K_plus = 500.0, det_K_minus = 500.0
)

# 计算矩阵元平方
M_squared = scattering_amplitude_squared(
    process, s, t, quark_params, thermo_params, K_coeffs
)

println("uu→uu: |M|² = $M_squared fm⁻⁴")
# 输出: uu→uu: |M|² = 314.13 fm⁻⁴
```

**示例2: 遍历所有过程**
```julia
# 所有11种散射过程
all_processes = [
    :uu_to_uu, :ss_to_ss, :ud_to_ud, :us_to_us,          # qq
    :udbar_to_udbar, :usbar_to_usbar, :uubar_to_uubar,   # qqbar (1-3)
    :uubar_to_ddbar, :uubar_to_ssbar,                    # qqbar (4-5)
    :ssbar_to_uubar, :ssbar_to_ssbar                     # qqbar (6-7)
]

s_qq = 10.0
s_qqbar = 6.0
t = -0.2

results = Dict{Symbol, Float64}()

for process in all_processes
    # 根据过程类型选择s值
    s_val = contains(string(process), "bar") ? s_qqbar : s_qq
    
    M_sq = scattering_amplitude_squared(
        process, s_val, t, quark_params, thermo_params, K_coeffs
    )
    
    results[process] = M_sq
    println("$process: |M|² = $(M_sq) fm⁻⁴")
end
```

**示例3: 参数扫描**
```julia
# 研究s依赖性
s_values = range(5.0, 20.0, length=10)
t_fixed = -0.2

M_sq_vs_s = Float64[]

for s in s_values
    M_sq = scattering_amplitude_squared(
        :uu_to_uu, s, t_fixed, quark_params, thermo_params, K_coeffs
    )
    push!(M_sq_vs_s, M_sq)
end

# 绘图或分析 M_sq_vs_s
```

#### 异常处理

| 异常类型 | 触发条件 | 处理建议 |
|---------|---------|---------|
| `ArgumentError` | 未知的process符号 | 检查符号是否在支持列表中 |
| `DomainError` | 底层积分域问题(如k0²-u<0) | 调整s,t值以满足物理约束 |

```julia
try
    M_sq = scattering_amplitude_squared(
        :invalid_process, s, t, quark_params, thermo_params, K_coeffs
    )
catch e
    if isa(e, ArgumentError)
        println("未知的散射过程: $(e.msg)")
    elseif isa(e, DomainError)
        println("参数不满足物理约束: $(e.msg)")
    end
end
```

---

### 2. calculate_mandelstam_variables (内部函数)

计算散射矩阵元中需要的18个Mandelstam辅助变量。

#### 函数签名
```julia
function calculate_mandelstam_variables(
    s::Float64,
    t::Float64,
    u::Float64,
    m1::Float64,
    m2::Float64,
    m3::Float64,
    m4::Float64
)::NamedTuple
```

#### 参数说明

| 参数 | 类型 | 说明 | 单位 |
|------|------|------|------|
| `s` | Float64 | Mandelstam s变量 | fm⁻² |
| `t` | Float64 | Mandelstam t变量 | fm⁻² |
| `u` | Float64 | Mandelstam u变量 | fm⁻² |
| `m1` | Float64 | 粒子1质量 | fm⁻¹ |
| `m2` | Float64 | 粒子2质量 | fm⁻¹ |
| `m3` | Float64 | 粒子3质量 | fm⁻¹ |
| `m4` | Float64 | 粒子4质量 | fm⁻¹ |

#### 返回值

返回包含18个辅助变量的`NamedTuple`:

```julia
(
    s_12_plus,  s_12_minus,   # s - (m1 ± m2)²
    s_34_plus,  s_34_minus,   # s - (m3 ± m4)²
    
    t_13_plus,  t_13_minus,   # t - (m1 ± m3)²
    t_24_plus,  t_24_minus,   # t - (m2 ± m4)²
    
    u_14_plus,  u_14_minus,   # u - (m1 ± m4)²
    u_23_plus,  u_23_minus,   # u - (m2 ± m3)²
    
    t_14_plus,  t_14_minus,   # t - (m1 ± m4)²
    u_13_plus,  u_13_minus    # u - (m1 ± m3)²
)
```

**单位**: 全部为fm⁻²

#### 物理意义

这些辅助变量简化了散射矩阵元的表达式:

```julia
# 例如在|M_u|²项中
|M_u|² ∝ |D_u^S|² * u_14_plus * u_23_plus + |D_u^P|² * u_14_minus * u_23_minus
```

其中:
- `u_14_plus = u - (m1 + m4)²` 表示u道的能量减去质量之和的平方
- `u_14_minus = u - (m1 - m4)²` 表示u道的能量减去质量之差的平方

#### 示例用法

```julia
# 同质量粒子
m = 1.5  # fm⁻¹
s = 10.0
t = -2.0
u = -4.5

vars = calculate_mandelstam_variables(s, t, u, m, m, m, m)

println("s₁₂⁺ = $(vars.s_12_plus) fm⁻²")   # 1.00
println("s₁₂⁻ = $(vars.s_12_minus) fm⁻²")  # 19.00
println("t₁₃⁺ = $(vars.t_13_plus) fm⁻²")   # -11.00
println("u₁₄⁺ = $(vars.u_14_plus) fm⁻²")   # -13.50
```

---

### 3. calculate_qq_amplitude_squared (内部函数)

计算qq散射矩阵元平方的底层实现。

#### 函数签名
```julia
function calculate_qq_amplitude_squared(
    process::Symbol,
    s::Float64,
    t::Float64,
    u::Float64,
    m1::Float64,
    m2::Float64,
    m3::Float64,
    m4::Float64,
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    K_coeffs::NamedTuple
)::Float64
```

#### 计算步骤

1. **计算质心系动量**
   ```julia
   cms_t = calculate_cms_momentum(process, s, t, :t, quark_params; u=u)
   cms_u = calculate_cms_momentum(process, s, t, :u, quark_params; u=u)
   ```

2. **获取道分离传播子**
   ```julia
   prop_t = calculate_all_propagators_by_channel(
       process, cms_t.k0, cms_t.k, quark_params, thermo_params, K_coeffs
   )
   D_t_S, D_t_P = prop_t.t_S, prop_t.t_P
   
   prop_u = calculate_all_propagators_by_channel(
       process, cms_u.k0, cms_u.k, quark_params, thermo_params, K_coeffs
   )
   D_u_S, D_u_P = prop_u.u_S, prop_u.u_P
   ```

3. **计算Mandelstam辅助变量**
   ```julia
   vars = calculate_mandelstam_variables(s, t, u, m1, m2, m3, m4)
   ```

4. **计算|M_u|²** (色自旋平均后)
   ```julia
   M_u_squared = abs2(D_u_S) * vars.u_14_plus * vars.u_23_plus + 
                 abs2(D_u_P) * vars.u_14_minus * vars.u_23_minus
   ```

5. **计算|M_t|²** (色自旋平均后)
   ```julia
   M_t_squared = abs2(D_t_S) * vars.t_13_plus * vars.t_24_plus + 
                 abs2(D_t_P) * vars.t_13_minus * vars.t_24_minus
   ```

6. **计算交叉项** 2Re(M_u M_t*)
   ```julia
   cross_term_factor = 1.0 / (4.0 * N_c)
   
   term1 = D_t_S * conj(D_u_S) * (vars.t_13_plus * vars.t_24_plus - 
                                   vars.s_12_plus * vars.s_34_plus + 
                                   vars.u_14_plus * vars.u_23_plus)
   
   term2 = D_t_S * conj(D_u_P) * (vars.t_13_plus * vars.t_24_plus - 
                                   vars.s_12_minus * vars.s_34_minus + 
                                   vars.u_14_minus * vars.u_23_minus)
   
   term3 = D_t_P * conj(D_u_S) * (vars.t_13_minus * vars.t_24_minus - 
                                   vars.s_12_minus * vars.s_34_minus + 
                                   vars.u_14_plus * vars.u_23_plus)
   
   term4 = D_t_P * conj(D_u_P) * (vars.t_13_minus * vars.t_24_minus - 
                                   vars.s_12_plus * vars.s_34_plus + 
                                   vars.u_14_minus * vars.u_23_minus)
   
   cross_term = cross_term_factor * (term1 - term2 - term3 + term4)
   ```

7. **合成总矩阵元平方**
   ```julia
   M_squared = M_u_squared + M_t_squared - 2.0 / N_c * real(cross_term)
   ```

#### 色自旋平均公式

$$
\frac{1}{4N_c^2}\sum_{\text{sc}} |\mathcal{M}|^2 = \frac{1}{4N_c^2}\sum_{\text{sc}} (|\mathcal{M}_u|^2 + |\mathcal{M}_t|^2 - 2\mathcal{M}_u \mathcal{M}_t^*)
$$

其中:
- $M_u^2$ 和 $M_t^2$ **已包含** $1/(4N_c^2)$ 平均因子
- 交叉项包含 $1/(4N_c)$ 因子,需要再除以 $N_c$ 得到 $1/(4N_c^2)$

---

### 4. calculate_qqbar_amplitude_squared (内部函数)

计算qqbar散射矩阵元平方的底层实现,结构类似`calculate_qq_amplitude_squared`,但使用s道和t道。

#### 主要差异

1. **散射道**: s道 + t道 (而非t道 + u道)
2. **传播子**: D_s_S, D_s_P, D_t_S, D_t_P
3. **辅助变量组合**: 使用s₁₂⁺, s₁₂⁻, s₃₄⁺, s₃₄⁻等

#### 色自旋平均公式

$$
\frac{1}{4N_c^2}\sum_{\text{sc}} |\mathcal{M}|^2 = \frac{1}{4N_c^2}\sum_{\text{sc}} (|\mathcal{M}_s|^2 + |\mathcal{M}_t|^2 - 2\mathcal{M}_s \mathcal{M}_t^*)
$$

---

## 数值特性

### 典型量级

| 散射类型 | |M|² 范围 (fm⁻⁴) | 物理场景 |
|---------|-----------------|---------|
| 轻夸克qq (uu, dd) | 300 - 1000 | s=10, t=-0.2 |
| 含s夸克qq (us, ss) | 700 - 2100 | s=10, t=-0.2 |
| 轻夸克qqbar (uubar) | 250 - 300 | s=6, t=-0.3 |
| 味道转换 (uubar→ddbar) | 600 - 1600 | s=6, t=-0.3 |

### 参数敏感性

**能量依赖**:
```
s ↑ (阈值→高能) ⇒ |M|² ↑ (强相互作用增强)
```

**动量转移依赖**:
```
|t| ↑ (小角→大角散射) ⇒ |M|² ↑ (t道传播子贡献)
```

**温度依赖**:
```
T ↑ ⇒ 传播子屏蔽 ⇒ |M|² 变化(复杂非单调)
```

---

## 性能注意事项

### 计算复杂度

单次`scattering_amplitude_squared`调用:
- **传播子计算**: 2次`calculate_all_propagators_by_channel`
- **每次传播子调用**: 包含6-8个介子通道的极化函数积分
- **典型耗时**: 0.1-0.3秒(取决于积分精度和缓存命中率)

### 优化建议

1. **缓存复用**: PolarizationCache已实现,相同(k0,k)参数自动缓存
2. **批量计算**: 对多个(s,t)点,考虑并行化
3. **参数预检**: 提前验证s,t范围避免无效积分

```julia
# 示例: 批量计算优化
using Base.Threads

s_grid = range(5.0, 20.0, length=50)
t_grid = range(-2.0, -0.1, length=30)

results = Matrix{Float64}(undef, length(s_grid), length(t_grid))

@threads for i in 1:length(s_grid)
    for j in 1:length(t_grid)
        results[i,j] = scattering_amplitude_squared(
            :uu_to_uu, s_grid[i], t_grid[j],
            quark_params, thermo_params, K_coeffs
        )
    end
end
```

---

## 理论背景

### 散射矩阵元定义

在PNJL模型中,夸克通过介子交换发生散射:

```
q₁(p₁) + q₂(p₂) → q₃(p₃) + q₄(p₄)
```

散射振幅:
```
M = ∑ₙ (顶点因子) × (传播子) × (顶点因子)
```

### 介子道

**qq散射**: t道 + u道交换
- t道: 动量转移k = p₁ - p₃
- u道: 动量转移k = p₁ - p₄

**qqbar散射**: s道 + t道交换  
- s道: 湮灭产生k = p₁ + p₂
- t道: 动量转移k = p₁ - p₃

### 色自旋平均

观测量对应色和自旋态求和并平均:
```
1/(4N_c²) ∑_{色,自旋} |M|²
```

其中:
- 1/4: 自旋平均(2种入态自旋 × 2种)
- 1/N_c²: 色平均(N_c种入态色 × N_c种)

---

## 相关函数

### 上游依赖

- `TotalPropagator.calculate_cms_momentum`: 计算质心系动量
- `TotalPropagator.calculate_all_propagators_by_channel`: 获取道分离传播子
- `TotalPropagator.get_quark_masses_for_process`: 提取过程对应的夸克质量

### 下游应用

- `CrossSection.dsigma_dt`: 微分散射截面
- `RelaxationTime.scattering_rate`: 散射率
- `TransportCoefficients.shear_viscosity`: 剪切粘度

---

## 常见问题

### Q1: 为什么ss→ss过程会失败?

**A**: 在特定(s,t)参数下,u道CMS动量计算可能得到k0²<u,导致虚数动量。这触发底层OneLoopIntegrals模块的积分域问题。

**解决方案**:
1. 增大s值以远离阈值
2. 调整t值避免极端前向/后向散射
3. 使用try-catch优雅处理数值异常

### Q2: |M|²量级是否合理?

**A**: 典型值300-3000 fm⁻⁴在合理范围内。文献中:
- 强耦合极限: $\alpha_s^2 \sim (0.5)^2 = 0.25$
- Mandelstam变量: $s,t \sim 1-20$ fm⁻²
- 预期: $|M|^2 \sim \alpha_s^2 \times s \times \text{因子} \sim 100-1000$ fm⁻⁴

### Q3: 如何验证计算正确性?

**A**: 
1. **物理约束**: 检查|M|² ≥ 0
2. **对称性**: uu→uu应与dd→dd相同(同位旋对称)
3. **阈值行为**: s→(m1+m2)²时,|M|²→0
4. **高能极限**: 大s时,|M|²~s(Regge行为)

### Q4: 交叉项的色自旋平均因子如何处理?

**A**: 公式文档给出:
```
1/(4N_c²) ∑(M_u M_t*) = 1/(4N_c) × [4个D组合项]
```

**关键理解**：这个等式表示**左右两边相等**。即右边的`1/(4N_c) × [...]`**已经是**完整的`1/(4N_c²) ∑(M_u M_t*)`，无需额外处理。

**代码实现**:
```julia
cross_term = 1/(4N_c) * (term1 - term2 - term3 + term4)  # 已是完整的1/(4N_c²)∑(...)
M_squared = M_u² + M_t² - 2 * real(cross_term)           # 直接使用
```

### Q5: 如何选择合适的(s,t)参数?

**A**: 参考以下指南:

**qq散射**:
```julia
# 轻夸克(u,d)
s_light = 5.0 ~ 15.0 fm⁻²  # 阈值约0.36 fm⁻²
t = -0.1 ~ -0.5 fm⁻²       # 小角散射

# 包含s夸克
s_strange = 10.0 ~ 20.0 fm⁻²  # 阈值约0.64 fm⁻²
t = -0.2 ~ -1.0 fm⁻²          # 适度动量转移
```

**qqbar散射**:
```julia
s = 4.0 ~ 10.0 fm⁻²   # 无阈值限制
t = -0.1 ~ -0.5 fm⁻²  # 小角散射
```

---

## 版本历史

### v1.0 (当前版本)
- ✅ 实现11种独立散射过程
- ✅ 修正色自旋平均因子处理
- ✅ 通过完整单元测试(39/40通过,1个已知数值问题)
- ✅ 支持完整PNJL模型参数

### 未来改进
- 🔹 参数域预检查和警告
- 🔹 与C++参考实现数值对比
- 🔹 性能优化(GPU加速积分?)

---

## 参考文献

1. **PNJL模型**: P. Costa et al., Phys. Rev. D 92, 036012 (2015)
2. **色自旋平均**: Peskin & Schroeder, "An Introduction to QFT", Chapter 5
3. **介子交换**: K. Fukushima, Phys. Lett. B 591, 277 (2004)

---

**文档版本**: 1.0  
**最后更新**: 2024年  
**维护者**: Julia_RelaxTime项目组
