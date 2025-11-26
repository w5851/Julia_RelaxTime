# 计划：微分散射截面计算模块 DifferentialCrossSection.jl

## 概述
实现相对论性Boltzmann动力学理论中的微分散射截面计算，连接已实现的散射矩阵元模块与后续的总散射截面/弛豫时间计算。

---

## 步骤 1: 理解公式与依赖关系

**公式**: 
$$\frac{d\sigma}{dt} = \frac{1}{16\pi s_{12}^+ s_{12}^-} \cdot \frac{1}{4N_C^2} \sum_{\text{spin, color}} |\mathcal{M}|^2$$

**关键点**:
1. **已实现部分**: $\frac{1}{4N_C^2} \sum |\mathcal{M}|^2$ 已由 `ScatteringAmplitude.jl` 中的 `scattering_amplitude_squared` 函数计算
2. **需实现部分**: 
   - 计算运动学因子 $\frac{1}{16\pi s_{12}^+ s_{12}^-}$
   - 其中 $s_{12}^+ = s - (m_1 + m_2)^2$，$s_{12}^- = s - (m_1 - m_2)^2$
3. **复用**: `calculate_mandelstam_variables` 已计算所需的 $s_{12}^\pm$ 变量

---

## 步骤 2: 创建核心计算函数

**文件**: `src/relaxtime/DifferentialCrossSection.jl`

**设计思路**: 利用已存在的 `calculate_mandelstam_variables` 函数预计算 s±, t±, u± 变量，避免内部重复计算

**核心函数（方案A - 推荐）**:
```julia
function differential_cross_section(
    process::Symbol,
    s_12_plus::Float64,
    s_12_minus::Float64,
    M_squared::Float64
) -> Float64
    """
    从预计算的运动学变量和散射矩阵元计算微分截面
    
    参数:
    - process: 散射过程标识（用于日志/调试）
    - s_12_plus: s - (m₁ + m₂)² [fm⁻²]
    - s_12_minus: s - (m₁ - m₂)² [fm⁻²]
    - M_squared: 已计算的 |M|² [fm⁻⁴]
    
    返回: dσ/dt [fm²]
    """
    kinematic_factor = 1.0 / (16π * s_12_plus * s_12_minus)
    return kinematic_factor * M_squared
end
```

**便捷包装函数（方案B - 一步到位）**:
```julia
function differential_cross_section(
    process::Symbol,
    s::Float64,
    t::Float64,
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    K_coeffs::NamedTuple
) -> Float64
    """
    一步计算微分截面（内部调用预计算）
    
    实现逻辑:
    1. 获取四个粒子质量
    2. 计算 u = m₁² + m₂² + m₃² + m₄² - s - t
    3. 调用 calculate_mandelstam_variables(s, t, u, m₁, m₂, m₃, m₄)
    4. 提取 s_12_plus, s_12_minus
    5. 调用 scattering_amplitude_squared 获取 |M|²
    6. 使用方案A的核心函数计算截面
    """
    # 步骤1: 获取质量
    m1, m2, m3, m4 = get_quark_masses_for_process(process, quark_params)
    
    # 步骤2-3: 计算所有Mandelstam辅助变量
    u = m1^2 + m2^2 + m3^2 + m4^2 - s - t
    mandelstam_vars = calculate_mandelstam_variables(s, t, u, m1, m2, m3, m4)
    
    # 步骤4-5: 获取矩阵元
    M_squared = scattering_amplitude_squared(
        process, s, t, quark_params, thermo_params, K_coeffs
    )
    
    # 步骤6: 计算截面
    return differential_cross_section(
        process, mandelstam_vars.s_12_plus, mandelstam_vars.s_12_minus, M_squared
    )
end
```

**推荐用法**:
```julia
# 高性能场景（批量计算时复用 mandelstam_vars）
mandelstam_vars = calculate_mandelstam_variables(s, t, u, m1, m2, m3, m4)
M_squared = scattering_amplitude_squared(...)
dsigma_dt = differential_cross_section(
    process, mandelstam_vars.s_12_plus, mandelstam_vars.s_12_minus, M_squared
)

# 便捷场景（单次计算）
dsigma_dt = differential_cross_section(
    process, s, t, quark_params, thermo_params, K_coeffs
)
```

**单位**: fm² (散射截面的自然单位)

---

## 步骤 3: 添加批量计算功能

**函数**: `calculate_all_differential_cross_sections`

**目的**: 对所有13种散射过程计算微分截面

**返回**: NamedTuple包含所有过程的 $d\sigma/dt$ 值

**优化实现（复用已计算的矩阵元）**:
```julia
function calculate_all_differential_cross_sections(
    s::Float64, t::Float64,
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    K_coeffs::NamedTuple
) -> NamedTuple
    """
    批量计算所有散射过程的微分截面
    
    优化策略: 先批量计算所有 |M|²，然后逐个添加运动学因子
    """
    # 步骤1: 批量计算所有散射矩阵元平方
    M_squared_all = calculate_all_scattering_amplitudes_squared(
        s, t, quark_params, thermo_params, K_coeffs
    )
    
    # 步骤2: 为每个过程添加运动学因子
    results = Dict{Symbol, Float64}()
    
    for (process, M_squared) in pairs(M_squared_all)
        # 获取质量
        m1, m2, m3, m4 = get_quark_masses_for_process(process, quark_params)
        
        # 计算运动学变量（仅需 s±）
        s_12_plus = s - (m1 + m2)^2
        s_12_minus = s - (m1 - m2)^2
        
        # 计算微分截面
        results[process] = differential_cross_section(
            process, s_12_plus, s_12_minus, M_squared
        )
    end
    
    return NamedTuple(results)
end
```

**性能优势**:
1. 复用 `calculate_all_scattering_amplitudes_squared` 的批量计算能力
2. 传播子计算只执行一次（最耗时部分）
3. 运动学因子计算开销很小（简单算术）

---

## 步骤 4: 验证物理约束与边界条件

**运动学检查**:
1. **阈值条件**: $s \geq (m_1 + m_2)^2$ 确保 $s_{12}^+ \geq 0$
2. **避免奇点**: 
   - 当 $s_{12}^+ \to 0$ 时，截面发散（阈值行为）
   - 当 $s_{12}^- \to 0$ 时（$m_1 = m_2$ 的特殊情况需谨慎处理）
3. **t 道约束**: $t_{\min} \leq t \leq t_{\max}$，其中:
   - $t_{\min} = -\frac{(s - (m_1+m_2)^2)(s - (m_3+m_4)^2)}{4s}$
   - $t_{\max} = 0$ (对于正向散射)

**建议实现**:
```julia
function check_kinematic_threshold(s::Float64, m1::Float64, m2::Float64)
    s_threshold = (m1 + m2)^2
    if s < s_threshold
        @warn "s below threshold" s=s threshold=s_threshold
        return false
    end
    
    s_plus = s - s_threshold
    if s_plus < 1e-12  # 接近阈值
        @warn "s very close to threshold, cross section may diverge"
    end
    
    return true
end
```

---

## 步骤 5: 模块结构与导出

**模块声明**:
```julia
module DifferentialCrossSection

"""
# DifferentialCrossSection.jl

微分散射截面计算模块，连接散射矩阵元与总散射截面/弛豫时间。

公式: dσ/dt = [1/(16π s₁₂⁺ s₁₂⁻)] × [1/(4Nc²) Σ|M|²]
其中 Σ|M|² 已由 ScatteringAmplitude 模块计算。

## 设计特点
- 双层接口: 核心函数（接受预计算变量）+ 便捷包装（一步到位）
- 批量优化: 复用 calculate_all_scattering_amplitudes_squared
- 运动学检查: 自动验证阈值条件和数值稳定性
"""

include("../Constants_PNJL.jl")
include("ScatteringAmplitude.jl")

using .Constants_PNJL: SCATTERING_MESON_MAP
using .ScatteringAmplitude: scattering_amplitude_squared,
                            calculate_all_scattering_amplitudes_squared,
                            get_quark_masses_for_process,
                            calculate_mandelstam_variables

# 预计算常数因子
const KINEMATIC_PREFACTOR = 1.0 / (16π)

# 导出函数
export differential_cross_section  # 核心函数（重载）
export calculate_all_differential_cross_sections  # 批量计算
export check_kinematic_threshold  # 运动学检查
export calculate_t_bounds  # t积分边界
```

**函数重载策略**:
```julia
# 方法1: 核心计算（高性能）
differential_cross_section(process::Symbol, s_12_plus::Float64, 
                          s_12_minus::Float64, M_squared::Float64)

# 方法2: 便捷包装（易用性）
differential_cross_section(process::Symbol, s::Float64, t::Float64,
                          quark_params::NamedTuple, thermo_params::NamedTuple,
                          K_coeffs::NamedTuple)
```

---

## 步骤 6: 测试用例设计

**测试文件**: `test_unit/test_differential_cross_section.jl`

**测试内容**:
1. **基本计算测试**:
   - 使用已知参数计算单个过程的 $d\sigma/dt$
   - 验证单位正确性（fm²）
   - 与散射矩阵元测试使用相同参数验证一致性

2. **运动学约束测试**:
   - 测试阈值附近行为 ($s \approx (m_1+m_2)^2$)
   - 验证截面为正 ($d\sigma/dt > 0$)
   - 检查 $s_{12}^-$ 为零的退化情况（$m_1 = m_2$）

3. **批量计算测试**:
   - 验证所有13种过程都能正确计算
   - 检查不同过程的截面量级合理性
   - 对比 uu 与 ss 散射的差异（质量效应）

4. **物理一致性测试**:
   - 验证同位旋对称: $d\sigma/dt(uu) = d\sigma/dt(dd)$（当 $m_u=m_d$）
   - 验证电荷共轭: $d\sigma/dt(\text{dubar}) = d\sigma/dt(\text{udbar})$
   - 检查 t→0 极限行为（正向散射增强）

**示例测试代码**:
```julia
@testset "微分散射截面基本计算" begin
    s = 31.0  # fm⁻²
    t = -2.0
    
    # 测试1: 使用便捷包装函数
    dsigma_dt = differential_cross_section(
        :uu_to_uu, s, t, quark_params, thermo_params, K_coeffs
    )
    
    @test dsigma_dt > 0.0
    @test !isnan(dsigma_dt) && !isinf(dsigma_dt)
    println("uu→uu: dσ/dt = ", dsigma_dt, " fm²")
    
    # 测试2: 验证与手动计算的一致性
    m1, m2, m3, m4 = get_quark_masses_for_process(:uu_to_uu, quark_params)
    u = m1^2 + m2^2 + m3^2 + m4^2 - s - t
    mandelstam_vars = calculate_mandelstam_variables(s, t, u, m1, m2, m3, m4)
    
    M_squared = scattering_amplitude_squared(
        :uu_to_uu, s, t, quark_params, thermo_params, K_coeffs
    )
    
    # 使用预计算变量的核心函数
    dsigma_dt_manual = differential_cross_section(
        :uu_to_uu, mandelstam_vars.s_12_plus, mandelstam_vars.s_12_minus, M_squared
    )
    
    @test dsigma_dt ≈ dsigma_dt_manual rtol=1e-12
    
    # 测试3: 验证运动学因子公式
    expected = M_squared / (16π * mandelstam_vars.s_12_plus * mandelstam_vars.s_12_minus)
    @test dsigma_dt_manual ≈ expected rtol=1e-10
end

@testset "预计算变量的性能优势" begin
    s = 31.0
    t = -2.0
    process = :uu_to_uu
    
    # 预计算方案（高性能）
    m1, m2, m3, m4 = get_quark_masses_for_process(process, quark_params)
    u = m1^2 + m2^2 + m3^2 + m4^2 - s - t
    mandelstam_vars = calculate_mandelstam_variables(s, t, u, m1, m2, m3, m4)
    M_squared = scattering_amplitude_squared(
        process, s, t, quark_params, thermo_params, K_coeffs
    )
    
    # 多次调用核心函数（模拟批量计算场景）
    results = Float64[]
    for i in 1:100
        dsigma = differential_cross_section(
            process, mandelstam_vars.s_12_plus, mandelstam_vars.s_12_minus, M_squared
        )
        push!(results, dsigma)
    end
    
    @test all(r -> r == results[1], results)  # 所有结果相同
    println("预计算方案: 100次调用，无重复计算")
end
```

---

## 步骤 7: 文档编写

**API文档**: `api/DifferentialCrossSection.md`

**内容**:
- 函数签名与参数说明
- 返回值物理意义（单位 fm²）
- 运动学约束与适用范围
- 使用示例
- 与总散射截面的积分关系: $\sigma_{\text{total}} = \int_{t_{\min}}^{t_{\max}} \frac{d\sigma}{dt} dt$

**公式文档**: 已存在 `doc/formula/微分散射截面by散射矩阵元.md`，需补充:
- 完整的运动学因子推导
- $s_{12}^\pm$ 的物理意义（相对动量）
- t 积分边界的推导

---

## 进一步考虑

**性能优化**:
1. ✅ **已采用**: 设计双层接口（核心函数 + 包装函数），支持预计算变量复用
2. ✅ **已采用**: 批量计算时先统一调用 `calculate_all_scattering_amplitudes_squared`
3. 预计算常数因子 `const KINEMATIC_PREFACTOR = 1.0 / (16π)` 作为模块常量
4. 可选：对于固定质量组合，缓存 s± 计算结果

**扩展功能**:
1. 添加 t 积分边界计算函数 `calculate_t_bounds(s, m1, m2, m3, m4)`
   ```julia
   function calculate_t_bounds(s::Float64, m1::Float64, m2::Float64, m3::Float64, m4::Float64)
       s_12_plus = s - (m1 + m2)^2
       s_34_plus = s - (m3 + m4)^2
       t_min = -(s_12_plus * s_34_plus) / (4*s)
       t_max = 0.0  # 正向散射
       return (t_min=t_min, t_max=t_max)
   end
   ```
2. 为后续总截面计算准备积分器接口
3. 支持打印详细的运动学信息（调试用）

**错误处理**:
1. 在核心函数中检查 `s_12_plus` 和 `s_12_minus` 的符号
   ```julia
   if s_12_plus <= 0.0
       error("Kinematic threshold violation: s_12_plus = $s_12_plus ≤ 0")
   end
   if abs(s_12_minus) < 1e-14  # 处理 m1 ≈ m2 的退化情况
       @warn "s_12_minus ≈ 0, using regularization"
       s_12_minus = max(abs(s_12_minus), 1e-14)
   end
   ```
2. 捕获运动学约束违反
3. 提供清晰的警告信息

**接口设计优势总结**:
- **灵活性**: 用户可选择便捷函数或高性能函数
- **可组合性**: 核心函数可与其他模块（如积分器）无缝集成
- **可测试性**: 分离运动学计算和矩阵元计算，便于单元测试
- **性能**: 批量计算时避免重复的传播子计算（主要瓶颈）

---

## 实现优先级

**高优先级（必须）**:
- 步骤2：核心计算函数 `differential_cross_section`
- 步骤4：运动学阈值检查
- 步骤6：基本测试用例

**中优先级（推荐）**:
- 步骤3：批量计算功能
- 步骤6：物理一致性测试

**低优先级（可选）**:
- 步骤7：详细文档
- 进一步考虑中的扩展功能

---

## 预期成果

完成后将具备：
1. ✅ 单个散射过程的微分截面计算
2. ✅ 所有13种过程的批量计算
3. ✅ 完整的运动学约束验证
4. ✅ 与散射矩阵元模块的无缝集成
5. ✅ 为总散射截面积分做好准备

这将为后续的弛豫时间计算奠定坚实基础。
