module TotalCrossSection

"""
# TotalCrossSection.jl

总散射截面计算模块，对微分散射截面进行 t 积分，包含末态统计因子。

## 物理公式

总散射截面通过对 Mandelstam 变量 t 积分得到：

σ(s, T, μ_q) = ∫_{t_-}^{t_+} dt · (dσ/dt) · (1 - f_c)(1 - f_d)

其中：
- dσ/dt: 微分散射截面（来自 DifferentialCrossSection 模块）
- (1 - f): 末态费米子统计因子（Pauli blocking）
- E_c, E_d: 末态粒子能量，从 (s, t) 计算
- t_±: t 积分边界（公式 5.14）

注：本模块只考虑夸克-夸克散射（费米子），不涉及玻色子散射。

## t 积分边界（修正版）

t_± = m_i² + m_c² - (1/2s)·(s + m_i² - m_j²)(s + m_c² - m_d²)
      ± 2√{[(s + m_i² - m_j²)²/(4s) - m_i²]·[(s + m_c² - m_d²)²/(4s) - m_c²]}

此公式基于质心系动量守恒和散射角的物理边界。

## 核心功能

1. **修正的 t 积分边界**: 使用完整公式 (5.14)，而非简化近似
2. **末态粒子能量**: 从 Mandelstam 变量计算质心系能量
3. **末态统计因子**: Pauli blocking (费米子专用)
4. **总散射截面 σ(s)**: 高斯-勒让德数值积分（固定点数，可预测耗时）
5. **批量计算与扫描**: 支持多过程、s 依赖性分析

## 使用示例

```julia
using .TotalCrossSection

# 单个过程的总散射截面
σ = total_cross_section(
    :uu_to_uu, 31.0,
    quark_params, thermo_params, K_coeffs,
    n_points=32  # 高斯积分点数
)

# 扫描 s 依赖性
s_values = collect(range(10.0, 50.0, length=20))
σ_values = scan_s_dependence(
    s_values, :uu_to_uu,
    quark_params, thermo_params, K_coeffs
)

# 计算所有过程
all_σ = calculate_all_total_cross_sections(
    31.0, quark_params, thermo_params, K_coeffs
)
```

## 参考文档

- doc/formula/散射截面by微分散射截面.md: 公式推导
- api/TotalCrossSection.md: API 详细文档
"""

# 导入依赖模块
include(joinpath(@__DIR__, "..", "Constants_PNJL.jl"))
include(joinpath(@__DIR__, "..", "integration", "GaussLegendre.jl"))
include(joinpath(@__DIR__, "ScatteringAmplitude.jl"))
include(joinpath(@__DIR__, "DifferentialCrossSection.jl"))
include(joinpath(@__DIR__, "OneLoopIntegrals.jl"))

using .Constants_PNJL: SCATTERING_MESON_MAP
using .GaussLegendre: gauleg
using .ScatteringAmplitude: scattering_amplitude_squared
using .DifferentialCrossSection: differential_cross_section
using .OneLoopIntegrals: distribution_value

# 默认积分点数
const DEFAULT_T_INTEGRAL_POINTS = 32

# 导出函数
export calculate_t_bounds
export calculate_final_state_energies
export final_state_blocking_factor
export combined_final_state_factor
export total_cross_section
export calculate_all_total_cross_sections
export scan_s_dependence

# 数值容差
const EPS_KINEMATIC = 1e-14  # 运动学阈值容差

# ----------------------------------------------------------------------------
# 1. t 积分边界计算（修正版，公式 5.14）
# ----------------------------------------------------------------------------

"""
    calculate_t_bounds(s, mi, mj, mc, md) -> NamedTuple

计算 Mandelstam 变量 t 的物理边界 t_±。

# 公式 (5.14)

t_± = m_i² + m_c² - (1/2s)·(s + m_i² - m_j²)(s + m_c² - m_d²)
      ± 2√{[(s + m_i² - m_j²)²/(4s) - m_i²]·[(s + m_c² - m_d²)²/(4s) - m_c²]}

其中：
- 常数项: m_i² + m_c²
- 对称项: -(1/2s)·(s + m_i² - m_j²)(s + m_c² - m_d²)
- 动量项: 2√(p_cm_i² · p_cm_c²)

# 参数

- `s::Float64`: Mandelstam 变量 s [fm⁻²]
- `mi::Float64`: 初态粒子 1 质量 [fm⁻¹]
- `mj::Float64`: 初态粒子 2 质量 [fm⁻¹]
- `mc::Float64`: 末态粒子 1 质量 [fm⁻¹]
- `md::Float64`: 末态粒子 2 质量 [fm⁻¹]

# 返回

`(t_min, t_max)`: t 的物理边界 [fm⁻²]

# 物理意义

- t_min (t_-): 后向散射极限
- t_max (t_+): 正向散射极限
- 满足 Mandelstam 约束: s + t + u = m_i² + m_j² + m_c² + m_d²

# 运动学检查

如果 s 低于阈值（初态或末态），抛出错误：
- 初态阈值: s ≥ (m_i + m_j)²
- 末态阈值: s ≥ (m_c + m_d)²

# 示例

```julia
s = 31.0  # fm⁻²
mi = mj = mc = md = 1.52  # fm⁻¹ (u 夸克)

t_bounds = calculate_t_bounds(s, mi, mj, mc, md)
# 对于相同质量: t_min ≈ -t_max (对称性)
```
"""
function calculate_t_bounds(
    s::Float64,
    mi::Float64,
    mj::Float64,
    mc::Float64,
    md::Float64
)::NamedTuple
    # 常数项
    constant_term = mi^2 + mc^2
    
    # 对称项的中间变量
    term_ij = s + mi^2 - mj^2
    term_cd = s + mc^2 - md^2
    symmetric_term = -(term_ij * term_cd) / (2.0 * s)
    
    # 平方根项（质心系动量的平方）
    # sqrt_arg1 = (term_ij²/(4s) - mi²) = p_cm_i²
    # sqrt_arg2 = (term_cd²/(4s) - mc²) = p_cm_c²
    sqrt_arg1 = (term_ij^2 / (4.0 * s)) - mi^2
    sqrt_arg2 = (term_cd^2 / (4.0 * s)) - mc^2
    
    # 检查运动学可行性
    if sqrt_arg1 < -EPS_KINEMATIC
        s_threshold_ij = (mi + mj)^2
        error("Initial state kinematic violation: s = $s < threshold = $s_threshold_ij")
    end
    if sqrt_arg2 < -EPS_KINEMATIC
        s_threshold_cd = (mc + md)^2
        error("Final state kinematic violation: s = $s < threshold = $s_threshold_cd")
    end
    
    # 处理接近阈值的情况（避免负数开方）
    sqrt_arg1 = max(0.0, sqrt_arg1)
    sqrt_arg2 = max(0.0, sqrt_arg2)
    
    # 计算动量项
    sqrt_term = 2.0 * sqrt(sqrt_arg1 * sqrt_arg2)
    
    # 最终结果: t_± = constant + symmetric ± sqrt
    t_minus = constant_term + symmetric_term - sqrt_term  # 后向散射
    t_plus  = constant_term + symmetric_term + sqrt_term  # 正向散射
    
    return (t_min = t_minus, t_max = t_plus)
end

# ----------------------------------------------------------------------------
# 2. 末态粒子能量计算
# ----------------------------------------------------------------------------

"""
    calculate_final_state_energies(s, t, mi, mj, mc, md) -> (E_c, E_d)

从 Mandelstam 变量计算质心系中末态粒子的能量。

# 公式

在质心系中：
- E_c = (s + m_c² - m_d²) / (2√s)
- E_d = (s + m_d² - m_c²) / (2√s)

# 参数

- `s::Float64`: Mandelstam 变量 s [fm⁻²]
- `t::Float64`: Mandelstam 变量 t [fm⁻²]
- `mi, mj, mc, md::Float64`: 粒子质量 [fm⁻¹]

# 返回

`(E_c, E_d)`: 末态粒子能量 [fm⁻¹]

# 物理验证

- 能量守恒: E_c + E_d = √s
- 质量壳条件: E_c² - p_c² = m_c²
- E_c ≥ m_c, E_d ≥ m_d

# 示例

```julia
s = 31.0
t = -2.0
mi = mj = mc = md = 1.52

E_c, E_d = calculate_final_state_energies(s, t, mi, mj, mc, md)
@assert E_c + E_d ≈ sqrt(s)  # 能量守恒
```
"""
function calculate_final_state_energies(
    s::Float64,
    t::Float64,
    mi::Float64,
    mj::Float64,
    mc::Float64,
    md::Float64
)::Tuple{Float64, Float64}
    # 质心系总能量
    sqrt_s = sqrt(s)
    
    # 末态粒子能量（质心系）
    E_c = (s + mc^2 - md^2) / (2.0 * sqrt_s)
    E_d = (s + md^2 - mc^2) / (2.0 * sqrt_s)
    
    return (E_c, E_d)
end

# ----------------------------------------------------------------------------
# 3. 末态统计因子（费米子专用）
# ----------------------------------------------------------------------------

"""
    final_state_blocking_factor(E, μ, T, Φ, Φbar) -> Float64

计算末态费米子统计因子 (1 - f)。

# 公式

费米子: 1 - f(E, μ, T, Φ, Φbar) (Pauli blocking)

其中 f 是 Fermi-Dirac 分布函数。

# 参数

- `E::Float64`: 粒子能量 [fm⁻¹]
- `μ::Float64`: 化学势 [fm⁻¹]
- `T::Float64`: 温度 [fm⁻¹]
- `Φ::Float64`: Polyakov loop
- `Φbar::Float64`: Conjugate Polyakov loop

# 返回

`(1 - f)`: 末态统计因子（无量纲，范围 [0, 1]）

# 物理意义

Pauli blocking: 费米子不能占据已被占据的态，抑制散射到已被占据的末态。

# 说明

本项目只考虑夸克-夸克散射（费米子），不涉及玻色子（介子）散射。

# 示例

```julia
E = 3.0  # fm⁻¹
μ = 0.3  # fm⁻¹
T = 0.15  # fm⁻¹
Φ = 0.5
Φbar = 0.5

factor = final_state_blocking_factor(E, μ, T, Φ, Φbar)
# 对于费米子: 0 ≤ factor ≤ 1
```
"""
function final_state_blocking_factor(
    E::Float64,
    μ::Float64,
    T::Float64,
    Φ::Float64,
    Φbar::Float64
)::Float64
    # 计算 Fermi-Dirac 分布函数 f(E, μ, T, Φ, Φbar)
    # 使用 :pnjl 模式，:plus 符号（正能粒子）
    f = distribution_value(:pnjl, :plus, E, μ, T, Φ, Φbar)
    
    # Pauli blocking (费米子)
    return 1.0 - f
end

"""
    combined_final_state_factor(E_c, E_d, μ_c, μ_d, T, Φ, Φbar) -> Float64

计算组合末态费米子统计因子 (1 - f_c)(1 - f_d)。

# 参数

- `E_c, E_d::Float64`: 末态粒子能量 [fm⁻¹]
- `μ_c, μ_d::Float64`: 化学势 [fm⁻¹]
- `T::Float64`: 温度 [fm⁻¹]
- `Φ, Φbar::Float64`: Polyakov loop

# 返回

`(1 - f_c)(1 - f_d)`: 组合统计因子（费米子）

# 说明

本项目只考虑夸克-夸克散射（费米子），因此统一使用 Pauli blocking 因子。

# 示例

```julia
factor = combined_final_state_factor(
    3.0, 3.0, 0.3, 0.3, 0.15, 0.5, 0.5
)
```
"""
function combined_final_state_factor(
    E_c::Float64,
    E_d::Float64,
    μ_c::Float64,
    μ_d::Float64,
    T::Float64,
    Φ::Float64,
    Φbar::Float64
)::Float64
    factor_c = final_state_blocking_factor(E_c, μ_c, T, Φ, Φbar)
    factor_d = final_state_blocking_factor(E_d, μ_d, T, Φ, Φbar)
    
    return factor_c * factor_d
end

# ----------------------------------------------------------------------------
# 5. 辅助函数：粒子、质量和化学势
# ----------------------------------------------------------------------------

"""
    parse_particles_from_process(process) -> Tuple{Symbol, Symbol, Symbol, Symbol}

从过程名称解析出粒子类型 (i, j, c, d)。

# 参数

- `process::Symbol`: 过程标识（如 :uu_to_uu, :ud_to_dū, :us_to_us）

# 返回

`(particle_i, particle_j, particle_c, particle_d)`: 粒子味标识

# 示例

```julia
i, j, c, d = parse_particles_from_process(:uu_to_uu)
# 返回: (:u, :u, :u, :u)

i, j, c, d = parse_particles_from_process(:ud_to_dū)
# 返回: (:u, :d, :d, :ubar)
```
"""
function parse_particles_from_process(process::Symbol)::Tuple{Symbol, Symbol, Symbol, Symbol}
    # 将过程名称转换为字符串
    process_str = string(process)
    
    # 分割 "_to_"
    parts = split(process_str, "_to_")
    if length(parts) != 2
        error("Invalid process format: $process (expected format: 'ab_to_cd')")
    end
    
    initial_str = String(parts[1])
    final_str = String(parts[2])
    
    # 解析初态粒子
    particle_i, particle_j = parse_particle_pair(initial_str)
    
    # 解析末态粒子
    particle_c, particle_d = parse_particle_pair(final_str)
    
    return (particle_i, particle_j, particle_c, particle_d)
end

"""
    parse_particle_pair(pair_str) -> Tuple{Symbol, Symbol}

从字符串解析粒子对（如 "uu" → (:u, :u), "dū" → (:d, :ubar)）。
"""
function parse_particle_pair(pair_str::String)::Tuple{Symbol, Symbol}
    # 处理反粒子标记 ū, đ, s̄
    pair_str = replace(pair_str, "ū" => "ubar")
    pair_str = replace(pair_str, "đ" => "dbar")  
    pair_str = replace(pair_str, "s̄" => "sbar")
    
    # 匹配模式：两个字母或 letterbar
    particles = Symbol[]
    
    i = 1
    while i <= length(pair_str)
        if i + 3 <= length(pair_str) && pair_str[i:i+3] == "ubar"
            push!(particles, :ubar)
            i += 4
        elseif i + 3 <= length(pair_str) && pair_str[i:i+3] == "dbar"
            push!(particles, :dbar)
            i += 4
        elseif i + 3 <= length(pair_str) && pair_str[i:i+3] == "sbar"
            push!(particles, :sbar)
            i += 4
        elseif pair_str[i] in ['u', 'd', 's']
            push!(particles, Symbol(pair_str[i]))
            i += 1
        else
            error("Unknown particle symbol at position $i in '$pair_str'")
        end
    end
    
    if length(particles) != 2
        error("Expected 2 particles in '$pair_str', found $(length(particles))")
    end
    
    return (particles[1], particles[2])
end

"""
    get_mass_for_flavor(flavor, quark_params) -> Float64

从 quark_params 获取指定味的夸克质量。

# 参数

- `flavor::Symbol`: 夸克味标识
- `quark_params::NamedTuple`: 夸克参数

# 返回

质量 [fm⁻¹]
"""
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

"""
    get_chemical_potential(flavor, quark_params) -> Float64

从 quark_params 获取指定味的化学势。

# 参数

- `flavor::Symbol`: 夸克味标识
- `quark_params::NamedTuple`: 夸克参数

# 返回

化学势 [fm⁻¹]

# 说明

反粒子化学势为负：μ(q̄) = -μ(q)
"""
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

# ----------------------------------------------------------------------------
# 6. 核心函数：总散射截面计算
# ----------------------------------------------------------------------------

"""
    total_cross_section(process, s, quark_params, thermo_params, K_coeffs; n_points) -> Float64

计算给定 s 下的总散射截面 σ(s)。

# 公式

σ(s) = ∫_{t_-}^{t_+} dt · (dσ/dt) · (1 - f_c)(1 - f_d)

其中：
- dσ/dt: 微分散射截面
- (1 - f_c)(1 - f_d): 末态费米子统计因子（Pauli blocking）
- t_±: 由 calculate_t_bounds 计算

# 参数

- `process::Symbol`: 散射过程标识（如 :uu_to_uu）
- `s::Float64`: Mandelstam 变量 s [fm⁻²]
- `quark_params::NamedTuple`: 夸克参数（质量、化学势）
- `thermo_params::NamedTuple`: 热力学参数（T, Φ, Φbar）
- `K_coeffs::NamedTuple`: 有效耦合常数
- `n_points::Int=32`: 高斯-勒让德积分点数（默认 32 点）

# 返回

`σ(s)`: 总散射截面 [fm²]

# 计算流程

1. 获取过程信息（粒子类型）
2. 获取质量和化学势
3. 计算 t 积分边界
4. 生成高斯-勒让德积分节点和权重
5. 对每个 t 点计算被积函数：
   - 计算散射矩阵元 |M|²
   - 计算微分截面 dσ/dt
   - 计算末态能量 E_c, E_d
   - 计算费米子统计因子 (1-f_c)(1-f_d)
6. 加权求和得到积分结果

# 说明

- 本模块只考虑夸克-夸克散射（费米子），不涉及玻色子散射
- 使用高斯-勒让德积分（固定点数），计算时间可预测
- 积分点数 n_points 越大精度越高，但耗时线性增加

# 性能考虑

- 每个 t 点需要计算 |M|²（~50 ms）
- 总耗时 ≈ n_points × 50 ms
- n_points=32 时，约 1.6 秒；n_points=16 时，约 0.8 秒

# 示例

```julia
σ = total_cross_section(
    :uu_to_uu, 31.0,
    quark_params, thermo_params, K_coeffs,
    n_points=32
)
println("σ(s=31) = \$σ fm²")
```
"""
function total_cross_section(
    process::Symbol,
    s::Float64,
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    K_coeffs::NamedTuple;
    n_points::Int=DEFAULT_T_INTEGRAL_POINTS
)::Float64
    # 步骤1: 解析过程中的粒子
    particle_i, particle_j, particle_c, particle_d = parse_particles_from_process(process)
    
    # 步骤2: 获取质量和化学势
    mi = get_mass_for_flavor(particle_i, quark_params)
    mj = get_mass_for_flavor(particle_j, quark_params)
    mc = get_mass_for_flavor(particle_c, quark_params)
    md = get_mass_for_flavor(particle_d, quark_params)
    
    μ_c = get_chemical_potential(particle_c, quark_params)
    μ_d = get_chemical_potential(particle_d, quark_params)
    
    # 步骤3: 计算 t 积分边界
    t_bounds = calculate_t_bounds(s, mi, mj, mc, md)
    t_min = t_bounds.t_min
    t_max = t_bounds.t_max
    
    # 步骤4: 生成高斯-勒让德积分节点和权重
    t_nodes, t_weights = gauleg(t_min, t_max, n_points)
    
    # 步骤5: 预计算 s 依赖的 Mandelstam 变量
    s_12_plus = s - (mi + mj)^2
    s_12_minus = s - (mi - mj)^2
    
    # 步骤6: 热力学参数
    T = thermo_params.T
    Φ = thermo_params.Φ
    Φbar = thermo_params.Φbar
    
    # 步骤7: 对 t 积分（高斯求积）
    σ = 0.0
    for i in 1:n_points
        t = t_nodes[i]
        w = t_weights[i]
        
        # 7.1 计算散射矩阵元
        M_squared = scattering_amplitude_squared(
            process, s, t, quark_params, thermo_params, K_coeffs
        )
        
        # 7.2 计算微分截面
        dsigma_dt = differential_cross_section(
            s_12_plus, s_12_minus, M_squared
        )
        
        # 7.3 计算末态能量
        E_c, E_d = calculate_final_state_energies(s, t, mi, mj, mc, md)
        
        # 7.4 计算末态统计因子（费米子）
        blocking_factor = combined_final_state_factor(
            E_c, E_d, μ_c, μ_d, T, Φ, Φbar
        )
        
        # 7.5 累加
        σ += w * dsigma_dt * blocking_factor
    end
    
    return σ
end

# ----------------------------------------------------------------------------
# 7. 批量计算与 s 扫描
# ----------------------------------------------------------------------------

"""
    calculate_all_total_cross_sections(s, quark_params, thermo_params, K_coeffs; n_points) -> NamedTuple

计算给定 s 下所有散射过程的总截面。

# 参数

与 `total_cross_section` 相同，但不需要指定 `process`

# 返回

`NamedTuple{Symbol, Float64}`: 所有过程的总截面

# 示例

```julia
all_σ = calculate_all_total_cross_sections(
    31.0, quark_params, thermo_params, K_coeffs,
    n_points=32
)
println("uu→uu: \$(all_σ.uu_to_uu) fm²")
println("dd→dd: \$(all_σ.dd_to_dd) fm²")
```
"""
function calculate_all_total_cross_sections(
    s::Float64,
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    K_coeffs::NamedTuple;
    n_points::Int=DEFAULT_T_INTEGRAL_POINTS
)::NamedTuple
    results = Dict{Symbol, Float64}()
    
    for process in keys(SCATTERING_MESON_MAP)
        try
            σ = total_cross_section(
                process, s, quark_params, thermo_params, K_coeffs,
                n_points=n_points
            )
            results[process] = σ
        catch e
            @warn "Failed to calculate $process" exception=e
            results[process] = NaN
        end
    end
    
    return NamedTuple(results)
end

"""
    scan_s_dependence(s_values, process, quark_params, thermo_params, K_coeffs; n_points) -> Vector{Float64}

扫描总散射截面随 s 的变化：σ(s)。

# 参数

- `s_values::Vector{Float64}`: s 值数组 [fm⁻²]
- 其他参数与 `total_cross_section` 相同

# 返回

`σ_values::Vector{Float64}`: 对应的总截面 [fm²]

# 示例

```julia
s_values = collect(range(10.0, 50.0, length=20))
σ_values = scan_s_dependence(
    s_values, :uu_to_uu,
    quark_params, thermo_params, K_coeffs,
    n_points=32
)

using Plots
plot(s_values, σ_values, xlabel="s [fm⁻²]", ylabel="σ [fm²]")
```
"""
function scan_s_dependence(
    s_values::Vector{Float64},
    process::Symbol,
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    K_coeffs::NamedTuple;
    n_points::Int=DEFAULT_T_INTEGRAL_POINTS
)::Vector{Float64}
    σ_values = Float64[]
    
    for s in s_values
        try
            σ = total_cross_section(
                process, s, quark_params, thermo_params, K_coeffs,
                n_points=n_points
            )
            push!(σ_values, σ)
        catch e
            @warn "Failed to calculate σ at s=$s" exception=e
            push!(σ_values, NaN)
        end
    end
    
    return σ_values
end

end  # module TotalCrossSection
