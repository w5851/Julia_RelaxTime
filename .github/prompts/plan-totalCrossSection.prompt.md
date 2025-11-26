# 计划：Part3 - 总散射截面计算模块 TotalCrossSection.jl

## 问题诊断

### 当前 t 积分边界实现的错误

**当前实现** (`DifferentialCrossSection.jl` 第274行)：
```julia
t_min = -(s_12_plus * s_34_plus) / (4.0 * s)
t_max = 0.0
```

**正确公式** (来自 `doc/formula/散射截面by微分散射截面.md`)：
```math
t_± = m_i^2 + m_c^2 - (1/2s)·(s + m_i^2 - m_j^2)(s + m_c^2 - m_d^2)
      ± 2√{[(s + m_i^2 - m_j^2)²/(4s) - m_i²]·[(s + m_c^2 - m_d^2)²/(4s) - m_c²]}
```

**问题分析**：
1. ❌ 当前实现缺少常数项 $m_i^2 + m_c^2$
2. ❌ 当前实现缺少对称项 $-\frac{1}{2s}(s + m_i^2 - m_j^2)(s + m_c^2 - m_d^2)$
3. ❌ 当前实现只计算了 $t_-$（使用负号），且公式也不完整
4. ✅ `t_max = 0` 的假设在某些情况下不正确

### 模块职责重新分配

**当前设计问题**：
- `calculate_t_bounds` 放在 `DifferentialCrossSection.jl` 中
- 但它实际上是 **总散射截面积分** 的准备工作
- 微分散射截面本身不需要知道积分边界

**建议调整**：
1. ✅ 将 `calculate_t_bounds` **移动**到 `TotalCrossSection.jl`
2. ✅ 修正公式为完整的 $t_\pm$ 计算
3. ✅ `DifferentialCrossSection.jl` 保持纯粹：只计算单点 $d\sigma/dt$

---

## Part3 实现计划

### 概述

实现总散射截面计算模块，对微分散射截面进行 t 积分，并考虑末态统计因子（Pauli blocking / Bose enhancement）。

**物理公式**：
```math
σ(s, T, μ_q) = ∫_{t_-}^{t_+} dt · (dσ/dt) · (1 ± f(E_c, T, μ)) · (1 ± f(E_d, T, μ))
```

其中：
- $(1 \pm f)$: 末态统计因子
  - $+$: 玻色子（Bose enhancement）
  - $-$: 费米子（Pauli blocking）
- $E_c, E_d$: 末态粒子能量（需从 s, t 计算）

---

## 步骤 1: 修正 t 积分边界计算

### 1.1 从 DifferentialCrossSection.jl 移除 `calculate_t_bounds`

**操作**：
1. 删除 `DifferentialCrossSection.jl` 中的 `calculate_t_bounds` 函数（第254-278行）
2. 移除导出：删除 `export calculate_t_bounds`（第57行）
3. 更新模块文档，移除相关示例

**理由**：
- 职责单一原则：微分截面模块不应负责积分边界
- 避免循环依赖：`TotalCrossSection` 会依赖 `DifferentialCrossSection`

### 1.2 在 TotalCrossSection.jl 中实现正确的 `calculate_t_bounds`

**新实现**：
```julia
function calculate_t_bounds(
    s::Float64,
    mi::Float64,  # 初态粒子1质量
    mj::Float64,  # 初态粒子2质量
    mc::Float64,  # 末态粒子1质量
    md::Float64   # 末态粒子2质量
)::NamedTuple
    """
    计算 Mandelstam 变量 t 的物理边界 t_±
    
    公式 (5.14):
    t_± = m_i² + m_c² - (1/2s)·(s + m_i² - m_j²)(s + m_c² - m_d²)
          ± 2√{[(s + m_i² - m_j²)²/(4s) - m_i²]·[(s + m_c² - m_d²)²/(4s) - m_c²]}
    """
    
    # 常数项
    constant_term = mi^2 + mc^2
    
    # 对称项
    term_ij = s + mi^2 - mj^2
    term_cd = s + mc^2 - md^2
    symmetric_term = -(term_ij * term_cd) / (2.0 * s)
    
    # 平方根项（质心系动量的函数）
    # sqrt_arg1 = (term_ij²/(4s) - mi²) = [(s - (mi+mj)²)(s - (mi-mj)²)] / (4s)
    sqrt_arg1 = (term_ij^2 / (4.0 * s)) - mi^2
    sqrt_arg2 = (term_cd^2 / (4.0 * s)) - mc^2
    
    # 检查运动学可行性
    if sqrt_arg1 < 0.0
        error("Initial state kinematic violation: s = $s, threshold = $((mi+mj)^2)")
    end
    if sqrt_arg2 < 0.0
        error("Final state kinematic violation: s = $s, threshold = $((mc+md)^2)")
    end
    
    sqrt_term = 2.0 * sqrt(sqrt_arg1 * sqrt_arg2)
    
    # 计算 t_± = constant + symmetric ± sqrt
    t_minus = constant_term + symmetric_term - sqrt_term  # 后向散射
    t_plus  = constant_term + symmetric_term + sqrt_term  # 正向散射
    
    return (t_min = t_minus, t_max = t_plus)
end
```

**验证关系**：
- 使用 Mandelstam 约束：$s + t + u = m_i^2 + m_j^2 + m_c^2 + m_d^2$
- 质心系动量：$p_{\text{cm}}^2 = \frac{(s - (m_i+m_j)^2)(s - (m_i-m_j)^2)}{4s}$

---

## 步骤 2: 计算末态粒子能量

### 2.1 从 Mandelstam 变量计算末态能量

**物理背景**：
在质心系中，已知 $(s, t)$ 后可以计算末态粒子能量。

**实现函数**：
```julia
function calculate_final_state_energies(
    s::Float64,
    t::Float64,
    mi::Float64,
    mj::Float64,
    mc::Float64,
    md::Float64
)::Tuple{Float64, Float64}
    """
    从 Mandelstam 变量计算质心系中末态粒子能量
    
    返回: (E_c, E_d) [fm⁻¹]
    """
    # 质心系总能量
    sqrt_s = sqrt(s)
    
    # 末态粒子能量（质心系）
    # E_c = (s + m_c² - m_d²) / (2√s)
    # E_d = (s + m_d² - m_c²) / (2√s)
    E_c = (s + mc^2 - md^2) / (2.0 * sqrt_s)
    E_d = (s + md^2 - mc^2) / (2.0 * sqrt_s)
    
    return (E_c, E_d)
end
```

**物理验证**：
- 能量守恒：$E_c + E_d = \sqrt{s}$
- 质量壳：$E_c^2 - p_c^2 = m_c^2$

### 2.2 确定粒子统计类型

**实现函数**：
```julia
function particle_statistics(flavor::Symbol)::Symbol
    """
    返回粒子的统计类型
    
    返回:
    - :fermion: 费米子（夸克）
    - :boson: 玻色子（反夸克、胶子）
    """
    # 对于 PNJL 模型中的夸克散射
    # 夸克是费米子，反夸克也按费米子处理（交换对称性）
    if flavor in [:u, :d, :s, :ubar, :dbar, :sbar]
        return :fermion
    else
        error("Unknown particle flavor: $flavor")
    end
end
```

---

## 步骤 3: 实现末态统计因子

### 3.1 分布函数计算

**复用现有模块**：
```julia
include("../src/relaxtime/OneLoopIntegrals.jl")
using .OneLoopIntegrals: distribution_value
```

**包装函数**：
```julia
function final_state_blocking_factor(
    E::Float64,
    μ::Float64,
    T::Float64,
    Φ::Float64,
    Φbar::Float64,
    statistics::Symbol
)::Float64
    """
    计算末态统计因子 (1 ± f)
    
    参数:
    - statistics: :fermion 或 :boson
    
    返回:
    - fermion: 1 - f (Pauli blocking)
    - boson:   1 + f (Bose enhancement)
    """
    # 计算分布函数 f(E, μ, T, Φ, Φbar)
    # 使用 :pnjl 模式，:plus 符号（正能粒子）
    f = distribution_value(:pnjl, :plus, E, μ, T, Φ, Φbar)
    
    if statistics == :fermion
        return 1.0 - f  # Pauli blocking
    elseif statistics == :boson
        return 1.0 + f  # Bose enhancement
    else
        error("Unknown statistics: $statistics")
    end
end
```

### 3.2 双粒子末态因子

**组合函数**：
```julia
function combined_final_state_factor(
    E_c::Float64,
    E_d::Float64,
    μ_c::Float64,
    μ_d::Float64,
    T::Float64,
    Φ::Float64,
    Φbar::Float64,
    stat_c::Symbol,
    stat_d::Symbol
)::Float64
    """
    计算组合末态因子: (1 ± f_c)(1 ± f_d)
    """
    factor_c = final_state_blocking_factor(E_c, μ_c, T, Φ, Φbar, stat_c)
    factor_d = final_state_blocking_factor(E_d, μ_d, T, Φ, Φbar, stat_d)
    
    return factor_c * factor_d
end
```

---

## 步骤 4: 实现总散射截面核心函数

### 4.1 核心积分函数

```julia
function total_cross_section(
    process::Symbol,
    s::Float64,
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    K_coeffs::NamedTuple;
    rtol::Float64=1e-6,
    atol::Float64=1e-10
)::Float64
    """
    计算给定 s 下的总散射截面 σ(s)
    
    公式:
    σ(s) = ∫_{t_-}^{t_+} dt · (dσ/dt) · (1 ± f_c)(1 ± f_d)
    
    参数:
    - process: 散射过程标识
    - s: Mandelstam 变量 s [fm⁻²]
    - rtol, atol: 积分容差
    
    返回: σ(s) [fm²]
    """
    
    # 步骤1: 获取过程信息
    process_info = SCATTERING_MESON_MAP[process]
    particles = process_info[:particles]  # (i, j, c, d)
    
    # 步骤2: 获取质量和化学势
    mi = get_mass_for_flavor(particles[1], quark_params)
    mj = get_mass_for_flavor(particles[2], quark_params)
    mc = get_mass_for_flavor(particles[3], quark_params)
    md = get_mass_for_flavor(particles[4], quark_params)
    
    μ_c = get_chemical_potential(particles[3], quark_params)
    μ_d = get_chemical_potential(particles[4], quark_params)
    
    # 步骤3: 获取统计类型
    stat_c = particle_statistics(particles[3])
    stat_d = particle_statistics(particles[4])
    
    # 步骤4: 计算 t 积分边界
    t_bounds = calculate_t_bounds(s, mi, mj, mc, md)
    
    # 步骤5: 预计算 s 依赖的 Mandelstam 变量
    s_12_plus = s - (mi + mj)^2
    s_12_minus = s - (mi - mj)^2
    
    # 步骤6: 对 t 积分
    T = thermo_params.T
    Φ = thermo_params.Φ
    Φbar = thermo_params.Φbar
    
    σ, err = quadgk(t_bounds.t_min, t_bounds.t_max, 
                    rtol=rtol, atol=atol) do t
        # 6.1 计算微分截面
        M_squared = scattering_amplitude_squared(
            process, s, t, quark_params, thermo_params, K_coeffs
        )
        
        dsigma_dt = differential_cross_section(
            s_12_plus, s_12_minus, M_squared
        )
        
        # 6.2 计算末态能量
        E_c, E_d = calculate_final_state_energies(s, t, mi, mj, mc, md)
        
        # 6.3 计算末态统计因子
        blocking_factor = combined_final_state_factor(
            E_c, E_d, μ_c, μ_d, T, Φ, Φbar, stat_c, stat_d
        )
        
        # 6.4 返回被积函数
        return dsigma_dt * blocking_factor
    end
    
    return σ
end
```

### 4.2 辅助函数

```julia
function get_mass_for_flavor(flavor::Symbol, quark_params::NamedTuple)::Float64
    if flavor in [:u, :ubar]
        return quark_params.m.u
    elseif flavor in [:d, :dbar]
        return quark_params.m.d
    elseif flavor in [:s, :sbar]
        return quark_params.m.s
    else
        error("Unknown flavor: $flavor")
    end
end

function get_chemical_potential(flavor::Symbol, quark_params::NamedTuple)::Float64
    # 反粒子化学势为负
    if flavor == :u
        return quark_params.μ.u
    elseif flavor == :ubar
        return -quark_params.μ.u
    elseif flavor == :d
        return quark_params.μ.d
    elseif flavor == :dbar
        return -quark_params.μ.d
    elseif flavor == :s
        return quark_params.μ.s
    elseif flavor == :sbar
        return -quark_params.μ.s
    else
        error("Unknown flavor: $flavor")
    end
end
```

---

## 步骤 5: 批量计算与 s 扫描

### 5.1 批量计算所有过程

```julia
function calculate_all_total_cross_sections(
    s::Float64,
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    K_coeffs::NamedTuple;
    rtol::Float64=1e-6,
    atol::Float64=1e-10
)::NamedTuple
    """
    计算给定 s 下所有散射过程的总截面
    
    返回: NamedTuple{Symbol, Float64}
    """
    results = Dict{Symbol, Float64}()
    
    for process in keys(SCATTERING_MESON_MAP)
        try
            σ = total_cross_section(
                process, s, quark_params, thermo_params, K_coeffs,
                rtol=rtol, atol=atol
            )
            results[process] = σ
        catch e
            @warn "Failed to calculate $process" exception=e
            results[process] = NaN
        end
    end
    
    return NamedTuple(results)
end
```

### 5.2 s 依赖性扫描

```julia
function scan_s_dependence(
    s_values::Vector{Float64},
    process::Symbol,
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    K_coeffs::NamedTuple;
    rtol::Float64=1e-6,
    atol::Float64=1e-10
)::Vector{Float64}
    """
    扫描总散射截面随 s 的变化: σ(s)
    
    返回: σ_values [fm²]
    """
    σ_values = Float64[]
    
    for s in s_values
        σ = total_cross_section(
            process, s, quark_params, thermo_params, K_coeffs,
            rtol=rtol, atol=atol
        )
        push!(σ_values, σ)
    end
    
    return σ_values
end
```

---

## 步骤 6: 测试用例设计

### 6.1 t 积分边界修正验证

**测试文件**: `test_unit/test_t_bounds_correction.jl`

```julia
@testset "t 积分边界修正验证" begin
    s = 31.0
    mi = mj = mc = md = 1.52  # uu→uu, 相同质量
    
    # 使用新公式计算
    t_bounds = calculate_t_bounds(s, mi, mj, mc, md)
    
    # 验证对称性（相同质量）
    @test abs(t_bounds.t_min + t_bounds.t_max) < 1e-10
    
    # 验证 Mandelstam 约束
    # 在 t=t_± 时，对应前向/后向散射
    # 使用 s + t + u = Σm² 验证
    u_min = (mi^2 + mj^2 + mc^2 + md^2) - s - t_bounds.t_max
    u_max = (mi^2 + mj^2 + mc^2 + md^2) - s - t_bounds.t_min
    
    @test u_min >= (mi - md)^2  # u 的运动学下限
    @test u_max >= (mi - md)^2
    
    println("新公式 t_bounds:")
    println("  t_min = $(t_bounds.t_min) fm⁻²")
    println("  t_max = $(t_bounds.t_max) fm⁻²")
    println("  对应 u_min = $u_min, u_max = $u_max fm⁻²")
end
```

### 6.2 末态能量计算验证

```julia
@testset "末态能量计算" begin
    s = 31.0
    t = -2.0
    mi = mj = mc = md = 1.52
    
    E_c, E_d = calculate_final_state_energies(s, t, mi, mj, mc, md)
    
    # 验证能量守恒
    @test E_c + E_d ≈ sqrt(s) rtol=1e-12
    
    # 验证质量壳条件（需要计算动量）
    sqrt_s = sqrt(s)
    p_cm_sq = ((s - (mi+mj)^2) * (s - (mi-mj)^2)) / (4.0 * s)
    
    @test E_c >= mc  # 能量至少等于静质量
    @test E_d >= md
    
    println("末态能量:")
    println("  E_c = $E_c fm⁻¹")
    println("  E_d = $E_d fm⁻¹")
    println("  E_c + E_d = $(E_c + E_d) = √s = $(sqrt_s)")
end
```

### 6.3 总散射截面基本计算

```julia
@testset "总散射截面计算" begin
    s = 31.0
    process = :uu_to_uu
    
    σ = total_cross_section(
        process, s, quark_params, thermo_params, K_coeffs,
        rtol=1e-4, atol=1e-8
    )
    
    @test σ > 0.0
    @test !isnan(σ) && !isinf(σ)
    
    println("uu→uu 总散射截面:")
    println("  σ(s=$s) = $σ fm²")
    
    # 估算量级：如果 dσ/dt ~ 0.1 fm², Δt ~ 4 fm⁻²
    # 则 σ ~ 0.4 fm²（粗略估计）
end
```

### 6.4 s 依赖性测试

```julia
@testset "s 依赖性扫描" begin
    # 选择 s 范围：从阈值附近到高能
    mi = mj = mc = md = quark_params.m.u
    s_threshold = (mi + mj)^2
    
    s_values = range(s_threshold * 1.1, s_threshold * 5.0, length=10)
    
    σ_values = scan_s_dependence(
        collect(s_values), :uu_to_uu,
        quark_params, thermo_params, K_coeffs,
        rtol=1e-4
    )
    
    # 验证物理行为
    @test all(σ -> σ > 0, σ_values)
    
    # 验证阈值附近行为（通常截面增大）
    # 注意：具体行为依赖于矩阵元的能量依赖性
    
    println("\ns 依赖性:")
    for (s, σ) in zip(s_values, σ_values)
        println("  s = $(round(s, digits=2)) fm⁻²  →  σ = $(round(σ, digits=4)) fm²")
    end
end
```

---

## 步骤 7: 模块结构与导出

```julia
module TotalCrossSection

"""
# TotalCrossSection.jl

总散射截面计算模块，对微分散射截面进行 t 积分。

公式 (5.13):
σ(s) = ∫_{t_-}^{t_+} dt · (dσ/dt) · (1 ± f_c)(1 ± f_d)

## 核心功能
- 修正的 t 积分边界计算（公式 5.14）
- 末态粒子能量计算
- 末态统计因子（Pauli blocking / Bose enhancement）
- 总散射截面 σ(s) 计算
- s 依赖性扫描
"""

using QuadGK

include("../Constants_PNJL.jl")
include("ScatteringAmplitude.jl")
include("DifferentialCrossSection.jl")
include("OneLoopIntegrals.jl")

using .Constants_PNJL: SCATTERING_MESON_MAP
using .ScatteringAmplitude: scattering_amplitude_squared
using .DifferentialCrossSection: differential_cross_section
using .OneLoopIntegrals: distribution_value

export calculate_t_bounds
export calculate_final_state_energies
export particle_statistics
export total_cross_section
export calculate_all_total_cross_sections
export scan_s_dependence

# ... (实现函数)

end  # module TotalCrossSection
```

---

## 步骤 8: 文档编写

### 8.1 API 文档

**文件**: `api/TotalCrossSection.md`

**内容**：
- 物理公式完整推导
- t_± 公式的详细说明
- 末态统计因子的物理意义
- 与 Part4（平均散射率）的关系
- 使用示例和典型结果

### 8.2 公式文档补充

**文件**: `doc/formula/散射截面by微分散射截面.md`

**补充内容**：
- t_± 公式的推导过程
- 质心系能量-动量关系
- 散射角与 Mandelstam 变量的关系

---

## 实现优先级

### 高优先级（必须）

1. **修正 t 边界计算**（步骤1）
   - 移除 `DifferentialCrossSection.jl` 中的旧实现
   - 在 `TotalCrossSection.jl` 中实现正确公式

2. **核心积分函数**（步骤4）
   - `total_cross_section` 主函数
   - 末态能量和统计因子

3. **基本测试**（步骤6.1-6.3）
   - t 边界验证
   - 末态能量验证
   - 总截面计算验证

### 中优先级（推荐）

4. **批量计算**（步骤5）
   - 所有过程的总截面
   - s 依赖性扫描

5. **完整测试**（步骤6.4）
   - s 依赖性测试
   - 物理行为验证

### 低优先级（可选）

6. **性能优化**
   - 自适应积分策略
   - 缓存 t 边界
   - 并行化 s 扫描

7. **详细文档**（步骤8）
   - API 文档
   - 公式推导补充

---

## 预期成果

完成后将具备：

1. ✅ **正确的 t 积分边界**: 使用公式 (5.14)，不是简化近似
2. ✅ **完整的物理模型**: 包含末态统计因子（Pauli blocking）
3. ✅ **总散射截面 σ(s)**: 作为 s 的函数，为 Part4 准备
4. ✅ **职责清晰**: 模块分工明确，无冗余功能
5. ✅ **充分验证**: 通过 Mandelstam 约束和能量守恒检验

这将为 Part4（平均散射率）的四重积分奠定坚实基础。

---

## 注意事项

### 数值积分考虑

1. **奇点处理**: 
   - 接近 $t_\pm$ 边界时微分截面可能发散
   - 使用自适应积分 `quadgk` 处理

2. **积分容差**:
   - 默认 `rtol=1e-6`, `atol=1e-10`
   - 可根据精度需求调整

3. **计算成本**:
   - 每个 t 点需要计算散射矩阵元（~50 ms）
   - 积分通常需要 10-100 个点
   - 单个 σ(s) 计算约 1-10 秒

### 物理合理性检查

1. **阈值行为**: σ(s) 在 $s \to (m_i+m_j)^2$ 时应符合阈值定律
2. **高能行为**: σ(s) 在高能时的渐近行为
3. **正定性**: σ(s) > 0 对所有物理参数
4. **统计因子**: $(1-f) \in [0, 1]$ 对费米子

### 与 Part4 的接口

Part4（平均散射率）将使用 σ(s) 作为输入：
```julia
Γ = ∫∫∫∫ dpi dθi dpj dθj · fi(pi) fj(pj) · σ(s(pi, pj, θij)) · v_rel
```

因此需要确保：
- σ(s) 在广泛的 s 范围内数值稳定
- 提供高效的 s 扫描功能
- 考虑向量化/并行化以加速四重积分
