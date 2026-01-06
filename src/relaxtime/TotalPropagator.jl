"""
# TotalPropagator.jl

总传播子计算模块，支持任意散射道组合及介子组合的总传播子计算。

在PNJL模型中，夸克弹性散射通过多种介子交换实现。本模块将单个介子传播子
按味因子加权求和，得到特定散射过程的总传播子，用于计算散射矩阵元。

公式参考 doc/formula/总传播子byPropagator.md

## 核心设计原则
- 支持任意散射过程（如 :uu_to_uu, :us_to_us, :ud_to_ud 等）
- 自动识别散射道（t/u/s）对应的味因子
- 灵活处理一般介子（π, K, σ_π, σ_K）和混合介子（η/η', σ/σ'）
- 性能优化：内部调用带缓存的极化函数，自动复用相同参数的计算结果
"""
module TotalPropagator

include("../Constants_PNJL.jl")
include("MesonPropagator.jl")
include("PolarizationCache.jl")

using .Constants_PNJL: SCATTERING_MESON_MAP
using .MesonPropagator: meson_propagator_simple, meson_propagator_mixed, calculate_coupling_matrix
using .PolarizationCache: polarization_aniso_cached, reset_cache!, get_cache_stats

export total_propagator_simple, total_propagator_mixed, get_flavor_factor
export calculate_all_propagators, calculate_all_propagators_by_channel, total_propagator_auto
export calculate_cms_momentum, get_quark_masses_for_process
export reset_cache!, get_cache_stats  # 导出缓存管理函数

# ----------------------------------------------------------------------------
# 味因子查询表（表5.3）
# ----------------------------------------------------------------------------

const FLAVOR_FACTOR_TABLE = Dict{Tuple{Symbol, Symbol}, Float64}(
    (:u, :u) => 1.0,
    (:u, :d) => √2,
    (:u, :s) => √2,
    (:d, :u) => √2,
    (:d, :d) => -1.0,
    (:d, :s) => √2,
    (:s, :u) => √2,
    (:s, :d) => √2,
    (:s, :s) => 2.0
)

# ----------------------------------------------------------------------------
# 辅助函数：质量提取
# ----------------------------------------------------------------------------

"""
    get_quark_masses_for_process(process::Symbol, quark_params::NamedTuple) -> NTuple{4, Float64}

根据散射过程提取四个粒子的质量。

# 参数
- `process`: 散射过程符号（如 :uu_to_uu, :uubar_to_ssbar）
- `quark_params`: 夸克参数NamedTuple，结构：
  ```julia
  (m = (u=m_u, d=m_d, s=m_s), ...)
  ```

# 返回值
返回元组 (m1, m2, m3, m4)，对应散射过程 q1+q2→q3+q4 的四个粒子质量（单位：fm⁻¹）

# 质量映射规则
- 夸克和反夸克使用相同的质量（如u和ū都用m_u）
- 从process符号解析粒子类型：如 :uubar_to_ssbar → (u, ū, s, š)
- 映射到quark_params.m中的对应质量

# 示例
```julia
quark_params = (m = (u=0.3, d=0.3, s=0.5), ...)
m1, m2, m3, m4 = get_quark_masses_for_process(:uubar_to_ssbar, quark_params)
# 返回: (0.3, 0.3, 0.5, 0.5)  # (m_u, m_u, m_s, m_s)
```
"""
function get_quark_masses_for_process(process::Symbol, quark_params::NamedTuple)
    # 解析散射过程
    q1, q2, q3, q4 = parse_scattering_process(process)
    
    # 提取味类型（忽略反粒子标记）
    flavor1 = extract_quark_flavor(q1)
    flavor2 = extract_quark_flavor(q2)
    flavor3 = extract_quark_flavor(q3)
    flavor4 = extract_quark_flavor(q4)
    
    # 获取对应质量
    m1 = quark_params.m[flavor1]
    m2 = quark_params.m[flavor2]
    m3 = quark_params.m[flavor3]
    m4 = quark_params.m[flavor4]
    
    return (m1, m2, m3, m4)
end

# ----------------------------------------------------------------------------
# 辅助函数：质心系动量计算
# ----------------------------------------------------------------------------

"""
    calculate_cms_momentum(process::Symbol, s::Float64, t::Float64, channel::Symbol,
                          quark_params::NamedTuple; u::Union{Float64, Nothing}=nothing) 
                          -> NamedTuple

计算质心系下介子的四动量，用于传播子计算。

# 物理背景
在散射过程 q₁+q₂→q₃+q₄ 中，介子交换发生在不同的散射道（s道、t道、u道）。
介子的四动量(k0, k)由质心系运动学和Mandelstam变量唯一确定。

# 质心系计算流程
1. 计算四个粒子在质心系的能量：
   - E1 = (s + m1² - m2²) / (2√s)
   - E2 = (s - m1² + m2²) / (2√s)
   - E3 = (s + m3² - m4²) / (2√s)
   - E4 = (s - m3² + m4²) / (2√s)

2. 根据散射道类型计算介子动量：
   - **s道**（q₁+q₂→介子→q₃+q₄）: k0 = |E1+E2|, k = 0
   - **t道**（q₁→q₃+介子）: k0 = |E1-E3|, k = √(k0²-t)
   - **u道**（q₁→q₄+介子）: k0 = |E1-E4|, k = √(k0²-u)

# 参数
- `process`: 散射过程符号（如 :uu_to_uu, :uubar_to_uubar）
- `s`: Mandelstam变量s（单位：fm⁻²）
- `t`: Mandelstam变量t（单位：fm⁻²）
- `channel`: 散射道类型（:s, :t, :u）
- `quark_params`: 夸克参数NamedTuple，用于提取质量
- `u`: 可选的Mandelstam变量u（单位：fm⁻²）。如果为`nothing`，则自动计算 u = Σm² - s - t

# 返回值
返回NamedTuple: `(k0 = k0, k = k)`
- `k0`: 介子能量分量（单位：fm⁻¹）
- `k`: 介子三动量大小（单位：fm⁻¹）

# 边界条件处理
- 当 k0²-t ≥ 0 或 k0²-u ≥ 0：正常计算 k = √(k0²-t/u)
- 当 -1e-12 < k0²-t/u < 0：数值误差，设 k = 0（发出警告）
- 当 k0²-t/u < -1e-12：数值误差，设 k = 0（发出警告）

# Mandelstam约束
s + t + u = m1² + m2² + m3² + m4²

如果不提供u，函数将自动计算以保证该约束严格满足。

# 示例
```julia
quark_params = (m = (u=0.3, d=0.3, s=0.5), ...)
s = 4.0  # fm⁻²
t = -0.5  # fm⁻²

# 自动计算u
result = calculate_cms_momentum(:uu_to_uu, s, t, :t, quark_params)
println("k0 = ", result.k0, " fm⁻¹")
println("k = ", result.k, " fm⁻¹")

# 使用预计算的u
u = sum_masses_squared - s - t
result = calculate_cms_momentum(:uu_to_uu, s, t, :t, quark_params; u=u)
```
"""
function calculate_cms_momentum(process::Symbol, s::Float64, t::Float64, channel::Symbol,
                               quark_params::NamedTuple; u::Union{Float64, Nothing}=nothing)
    # 1. 提取四个粒子质量
    m1, m2, m3, m4 = get_quark_masses_for_process(process, quark_params)
    
    # 2. 计算u（如果未提供）
    if u === nothing
        u = m1^2 + m2^2 + m3^2 + m4^2 - s - t
    end
    
    # 3. 计算质心系能量
    sqrt_s = sqrt(s)
    E1 = (s + m1^2 - m2^2) / (2.0 * sqrt_s)
    E2 = (s - m1^2 + m2^2) / (2.0 * sqrt_s)
    E3 = (s + m3^2 - m4^2) / (2.0 * sqrt_s)
    E4 = (s - m3^2 + m4^2) / (2.0 * sqrt_s)
    
    # 4. 根据散射道计算介子四动量
    local k0::Float64, k::Float64
    
    if channel == :s
        # s道：k0 = |E1 + E2|, k = 0
        k0 = abs(E1 + E2)
        k = 0.0
        
    elseif channel == :t
        # t道：k0 = |E1 - E3|, k = √(k0² - t)
        k0 = abs(E1 - E3)
        delta = k0^2 - t
        
        if delta >= 0.0
            k = sqrt(delta)
        elseif delta > -1e-12
            # 数值误差范围内，设为0
            k = 0.0
            @warn "calculate_cms_momentum: k0²-t = $delta ∈ (-1e-12, 0), setting k=0" maxlog=10
        else
            # 明显的数值问题，设为0并警告
            k = 0.0
            @warn "calculate_cms_momentum: k0²-t = $delta < -1e-12, setting k=0. Check Mandelstam variables." maxlog=10
        end
        
    elseif channel == :u
        # u道：k0 = |E1 - E4|, k = √(k0² - u)
        k0 = abs(E1 - E4)
        delta = k0^2 - u
        
        if delta >= 0.0
            k = sqrt(delta)
        elseif delta > -1e-12
            # 数值误差范围内，设为0
            k = 0.0
            @warn "calculate_cms_momentum: k0²-u = $delta ∈ (-1e-12, 0), setting k=0" maxlog=10
        else
            # 明显的数值问题，设为0并警告
            k = 0.0
            @warn "calculate_cms_momentum: k0²-u = $delta < -1e-12, setting k=0. Check Mandelstam variables." maxlog=10
        end
        
    else
        error("Unknown channel: $channel. Use :s, :t, or :u")
    end
    
    return (k0 = k0, k = k)
end

# ----------------------------------------------------------------------------
# 辅助函数
# ----------------------------------------------------------------------------

function extract_quark_flavor(quark::Symbol)
    s = string(quark)
    if endswith(s, "bar")
        return Symbol(s[1:end-3])
    else
        return quark
    end
end

function get_flavor_factor(quark1::Symbol, quark2::Symbol)
    flavor1 = extract_quark_flavor(quark1)
    flavor2 = extract_quark_flavor(quark2)
    
    key = (flavor1, flavor2)
    if haskey(FLAVOR_FACTOR_TABLE, key)
        return FLAVOR_FACTOR_TABLE[key]
    else
        error("Unknown flavor combination: ($flavor1, $flavor2)")
    end
end

"""
    parse_particle_pair_str(pair_str::AbstractString) -> Tuple{Symbol, Symbol}

将过程符号里形如 `"ud"`, `"uubar"`, `"ubarsbar"` 的二粒子字符串解析为两个粒子 `Symbol`。

支持的 token：
- 夸克：`u`, `d`, `s`
- 反夸克：`ubar`, `dbar`, `sbar`

约束：
- 解析结果必须恰好是 2 个粒子，否则抛错。

示例：
```julia
parse_particle_pair_str("ud")       # (:u, :d)
parse_particle_pair_str("uubar")    # (:u, :ubar)
parse_particle_pair_str("ubarubar") # (:ubar, :ubar)
parse_particle_pair_str("ubardbar") # (:ubar, :dbar)
```
"""
function parse_particle_pair_str(pair_str::AbstractString)::Tuple{Symbol, Symbol}
    # Token-based parser that supports both quarks and antiquarks, e.g.
    # "ud" -> (:u,:d)
    # "uubar" -> (:u,:ubar)
    # "ubarubar" -> (:ubar,:ubar)
    # "ubardbar" -> (:ubar,:dbar)
    s = String(pair_str)
    particles = Symbol[]
    i = firstindex(s)
    while i <= lastindex(s)
        rest = SubString(s, i)
        if startswith(rest, "ubar")
            push!(particles, :ubar)
            i = nextind(s, i, 4)
        elseif startswith(rest, "dbar")
            push!(particles, :dbar)
            i = nextind(s, i, 4)
        elseif startswith(rest, "sbar")
            push!(particles, :sbar)
            i = nextind(s, i, 4)
        else
            c = s[i]
            if c == 'u'
                push!(particles, :u)
            elseif c == 'd'
                push!(particles, :d)
            elseif c == 's'
                push!(particles, :s)
            else
                error("Invalid particle token at position $i in '$s'")
            end
            i = nextind(s, i)
        end
    end
    length(particles) == 2 || error("Invalid particle pair '$s' (expected 2 particles)")
    return (particles[1], particles[2])
end

"""
    parse_scattering_process(process::Symbol) -> Tuple{Symbol, Symbol, Symbol, Symbol}

将散射过程 key（例如 `:uu_to_uu`, `:uubar_to_ssbar`, `:ubardbar_to_ubardbar`）解析为
四个粒子符号 `(q1, q2, q3, q4)`，对应反应 `q1 + q2 → q3 + q4`。

说明：
- 输入格式要求为 `q1q2_to_q3q4`（以 `"_to_"` 分隔）。
- 粒子对字符串的解析由 `parse_particle_pair_str` 完成，因此同时支持夸克与反夸克 token。

该解析结果会被用于：
- `get_quark_masses_for_process`：提取四个粒子的质量（忽略 bar 标记取味）
- `get_flavor_factors_for_channel`：为 s/t/u 道计算味因子组合
- 混合介子传播子的味结构选择
"""
function parse_scattering_process(process::Symbol)
    s = string(process)
    parts = split(s, "_to_")
    length(parts) == 2 || error("Invalid process format: $process. Expected: q1q2_to_q3q4")

    q1, q2 = parse_particle_pair_str(parts[1])
    q3, q4 = parse_particle_pair_str(parts[2])
    return (q1, q2, q3, q4)
end

function get_flavor_factors_for_channel(process::Symbol, channel::Symbol)
    q1, q2, q3, q4 = parse_scattering_process(process)
    
    if channel == :t
        T1 = get_flavor_factor(q1, q3)
        T2 = get_flavor_factor(q2, q4)
    elseif channel == :u
        T1 = get_flavor_factor(q1, q4)
        T2 = get_flavor_factor(q2, q3)
    elseif channel == :s
        T1 = get_flavor_factor(q1, q2)
        T2 = get_flavor_factor(q3, q4)
    else
        error("Unknown channel: $channel. Use :t, :u, or :s")
    end
    
    return (T1, T2)
end

# ----------------------------------------------------------------------------
# 辅助函数：介子类型与极化函数参数映射
# ----------------------------------------------------------------------------

"""
    get_polarization_params(meson::Symbol, quark_params::NamedTuple) -> NamedTuple

根据介子类型返回计算极化函数所需的参数。

# 参数映射表
| 介子 | 通道 | m1 | m2 | μ1 | μ2 | A1 | A2 | num_s_quark |
|------|------|----|----|----|----|----|----|-------------|
| :pi | :P | m_u | m_u | μ_u | μ_u | A_u | A_u | 0 |
| :K | :P | m_u | m_s | μ_u | μ_s | A_u | A_s | 1 |
| :sigma_pi | :S | m_u | m_u | μ_u | μ_u | A_u | A_u | 0 |
| :sigma_K | :S | m_u | m_s | μ_u | μ_s | A_u | A_s | 1 |

# 参数
- `meson`: 介子类型（:pi, :K, :sigma_pi, :sigma_K）
- `quark_params`: 夸克参数NamedTuple，结构为：
  ```julia
  (m = (u=..., d=..., s=...), 
   μ = (u=..., d=..., s=...), 
   A = (u=..., d=..., s=...))
  ```

# 返回值
返回NamedTuple：(channel, m1, m2, μ1, μ2, A1, A2, num_s_quark)
"""
function get_polarization_params(meson::Symbol, quark_params::NamedTuple)
    if meson == :pi
        return (channel = :P,
                m1 = quark_params.m.u, m2 = quark_params.m.u,
                μ1 = quark_params.μ.u, μ2 = quark_params.μ.u,
                A1 = quark_params.A.u, A2 = quark_params.A.u,
                num_s_quark = 0)
    elseif meson == :K
        return (channel = :P,
                m1 = quark_params.m.u, m2 = quark_params.m.s,
                μ1 = quark_params.μ.u, μ2 = quark_params.μ.s,
                A1 = quark_params.A.u, A2 = quark_params.A.s,
                num_s_quark = 1)
    elseif meson == :sigma_pi
        return (channel = :S,
                m1 = quark_params.m.u, m2 = quark_params.m.u,
                μ1 = quark_params.μ.u, μ2 = quark_params.μ.u,
                A1 = quark_params.A.u, A2 = quark_params.A.u,
                num_s_quark = 0)
    elseif meson == :sigma_K
        return (channel = :S,
                m1 = quark_params.m.u, m2 = quark_params.m.s,
                μ1 = quark_params.μ.u, μ2 = quark_params.μ.s,
                A1 = quark_params.A.u, A2 = quark_params.A.s,
                num_s_quark = 1)
    else
        error("Unknown meson type: $meson. Use :pi, :K, :sigma_pi, or :sigma_K")
    end
end

# ----------------------------------------------------------------------------
# 主函数
# ----------------------------------------------------------------------------

"""
    total_propagator_simple(process, channel, meson_list, k0, k_norm, 
                           quark_params, thermo_params, K_coeffs) -> ComplexF64

计算一般介子的总传播子，包含味因子。内部自动调用带缓存的极化函数。

# 公式
D_total = T1 * sum(D_meson) * T2

其中每个介子传播子通过内部计算的极化函数获得。

# 参数
- `process`: 散射过程符号（如 :uu_to_uu）
- `channel`: 散射道类型（:t, :u, :s）
- `meson_list`: 介子列表（如 [:pi, :K]）
- `k0`: 四动量能量分量（单位：fm⁻¹）
- `k_norm`: 三动量大小（单位：fm⁻¹）
- `quark_params`: 夸克参数NamedTuple，结构：
  ```julia
  (m = (u=m_u, d=m_d, s=m_s),   # 夸克质量（fm⁻¹）
   μ = (u=μ_u, d=μ_d, s=μ_s),   # 夸克化学势（fm⁻¹）
   A = (u=A_u, d=A_d, s=A_s))   # A函数值（fm）
  ```
- `thermo_params`: 热力学参数NamedTuple，结构：
  ```julia
  (T=T, Φ=Φ, Φbar=Φbar, ξ=ξ)    # 温度（fm⁻¹）、Polyakov环、各向异性参数
  ```
- `K_coeffs`: K系数NamedTuple（通过`EffectiveCouplings.calculate_effective_couplings`获取）

# 返回值
总传播子（ComplexF64，单位 fm²）

# 示例
```julia
using .EffectiveCouplings: calculate_effective_couplings, calculate_G_from_A
using .OneLoopIntegrals: A

# 设置参数
T = 150.0 / 197.327  # 150 MeV → fm⁻¹
k0 = 100.0 / 197.327
k_norm = 50.0 / 197.327

# 计算夸克参数
m_u = 300.0 / 197.327
m_s = 500.0 / 197.327
μ_u = 0.0
μ_s = 0.0
Φ = 0.5
Φbar = 0.5

A_u = A(T, μ_u, m_u, Φ, Φbar)
A_s = A(T, μ_s, m_s, Φ, Φbar)

quark_params = (
    m = (u=m_u, d=m_u, s=m_s),
    μ = (u=μ_u, d=μ_u, s=μ_s),
    A = (u=A_u, d=A_u, s=A_s)
)

thermo_params = (T=T, Φ=Φ, Φbar=Φbar, ξ=0.0)

# 计算K系数
G_u = calculate_G_from_A(A_u)
G_s = calculate_G_from_A(A_s)
K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

# 计算总传播子
D = total_propagator_simple(:uu_to_uu, :t, [:pi, :K], 
                           k0, k_norm, quark_params, thermo_params, K_coeffs)
```
"""
function total_propagator_simple(process::Symbol, channel::Symbol, meson_list::Vector{Symbol},
                                k0::Float64, k_norm::Float64,
                                quark_params::NamedTuple, thermo_params::NamedTuple,
                                K_coeffs::NamedTuple)
    T1, T2 = get_flavor_factors_for_channel(process, channel)
    
    D_sum = ComplexF64(0.0, 0.0)
    for meson in meson_list
        # 获取该介子的极化函数参数
        pol_params = get_polarization_params(meson, quark_params)
        
        # 调用带缓存的极化函数
        Π_real, Π_imag = polarization_aniso_cached(
            pol_params.channel, k0, k_norm,
            pol_params.m1, pol_params.m2,
            pol_params.μ1, pol_params.μ2,
            thermo_params.T, thermo_params.Φ, thermo_params.Φbar, thermo_params.ξ,
            pol_params.A1, pol_params.A2,
            pol_params.num_s_quark
        )
        
        Π = ComplexF64(Π_real, Π_imag)
        
        # 计算介子传播子
        D_meson = meson_propagator_simple(meson, K_coeffs, Π)
        D_sum += D_meson
    end
    
    D_total = T1 * D_sum * T2
    return D_total
end

"""
    total_propagator_mixed(process, channel, meson_channel, k0, k_norm,
                          quark_params, thermo_params, K_coeffs) -> ComplexF64

计算混合介子（η/η'或σ/σ'）的总传播子。内部自动调用带缓存的极化函数。

注意：混合介子传播子已包含味结构，无需额外乘味因子。

# 参数
- `process`: 散射过程符号（如 :uu_to_uu）
- `channel`: 散射道类型（:t, :u, :s）
- `meson_channel`: 介子通道类型（:P为赝标量η/η'，:S为标量σ/σ'）
- `k0`: 四动量能量分量（单位：fm⁻¹）
- `k_norm`: 三动量大小（单位：fm⁻¹）
- `quark_params`: 夸克参数NamedTuple，结构：
  ```julia
  (m = (u=m_u, d=m_d, s=m_s),   # 夸克质量（fm⁻¹）
   μ = (u=μ_u, d=μ_d, s=μ_s),   # 夸克化学势（fm⁻¹）
   A = (u=A_u, d=A_d, s=A_s))   # A函数值（fm）
  ```
- `thermo_params`: 热力学参数NamedTuple，结构：
  ```julia
  (T=T, Φ=Φ, Φbar=Φbar, ξ=ξ)    # 温度（fm⁻¹）、Polyakov环、各向异性参数
  ```
- `K_coeffs`: K系数NamedTuple（通过`EffectiveCouplings.calculate_effective_couplings`获取）

# 返回值
混合介子传播子（ComplexF64，单位 fm²）

# 极化函数计算规则
- Π_uu：使用 (m1=m_u, m2=m_u, μ1=μ_u, μ2=μ_u, A1=A_u, A2=A_u, num_s_quark=0)
- Π_ss：使用 (m1=m_s, m2=m_s, μ1=μ_s, μ2=μ_s, A1=A_s, A2=A_s, num_s_quark=2)

# 示例
```julia
# 计算η/η'传播子（赝标量通道）
D_eta = total_propagator_mixed(:uu_to_uu, :t, :P, 
                               k0, k_norm, quark_params, thermo_params, K_coeffs)

# 计算σ/σ'传播子（标量通道）
D_sigma = total_propagator_mixed(:uu_to_uu, :t, :S,
                                 k0, k_norm, quark_params, thermo_params, K_coeffs)
```
"""
function total_propagator_mixed(process::Symbol, channel::Symbol, meson_channel::Symbol,
                               k0::Float64, k_norm::Float64,
                               quark_params::NamedTuple, thermo_params::NamedTuple,
                               K_coeffs::NamedTuple)
    # 计算Π_uu（u夸克极化函数）
    Π_uu_real, Π_uu_imag = polarization_aniso_cached(
        meson_channel, k0, k_norm,
        quark_params.m.u, quark_params.m.u,
        quark_params.μ.u, quark_params.μ.u,
        thermo_params.T, thermo_params.Φ, thermo_params.Φbar, thermo_params.ξ,
        quark_params.A.u, quark_params.A.u,
        0  # num_s_quark = 0 for uu
    )
    Π_uu = ComplexF64(Π_uu_real, Π_uu_imag)
    
    # 计算Π_ss（s夸克极化函数）
    Π_ss_real, Π_ss_imag = polarization_aniso_cached(
        meson_channel, k0, k_norm,
        quark_params.m.s, quark_params.m.s,
        quark_params.μ.s, quark_params.μ.s,
        thermo_params.T, thermo_params.Φ, thermo_params.Φbar, thermo_params.ξ,
        quark_params.A.s, quark_params.A.s,
        2  # num_s_quark = 2 for ss
    )
    Π_ss = ComplexF64(Π_ss_real, Π_ss_imag)
    
    # 计算M矩阵
    M_matrix = calculate_coupling_matrix(Π_uu, Π_ss, K_coeffs, meson_channel)
    
    # 从K_coeffs中获取det_K
    det_K = (meson_channel == :P) ? K_coeffs.det_K_plus : K_coeffs.det_K_minus
    
    # 计算混合介子传播子
    q1, q2, q3, q4 = parse_scattering_process(process)
    D_mixed = meson_propagator_mixed(det_K, M_matrix, q1, q2, q3, q4, channel)
    
    return D_mixed
end

# ----------------------------------------------------------------------------
# 自动化接口：根据散射过程自动计算传播子
# ----------------------------------------------------------------------------

"""
    total_propagator_auto(process, channel, k0, k_norm, 
                         quark_params, thermo_params, K_coeffs) -> ComplexF64

根据散射过程和散射道自动查询介子列表并计算总传播子。

# 功能说明
1. 从常量映射表 `SCATTERING_MESON_MAP` 自动查询介子列表
2. 分别计算一般介子和混合介子的贡献
3. 将所有贡献相加得到该散射道的总传播子

# 介子计算规则
- **一般介子**（π, K, σ_π, σ_K）：调用 `total_propagator_simple`
- **混合介子η/η'**（`:mixed_P => true`）：调用 `total_propagator_mixed(..., :P)` 一次
- **混合介子σ/σ'**（`:mixed_S => true`）：调用 `total_propagator_mixed(..., :S)` 一次

# 参数
- `process`: 散射过程符号（如 :uu_to_uu, :uubar_to_uubar）
- `channel`: 散射道类型（:t, :u, :s）
- `k0`: 四动量能量分量（单位：fm⁻¹）
- `k_norm`: 三动量大小（单位：fm⁻¹）
- `quark_params`: 夸克参数NamedTuple
- `thermo_params`: 热力学参数NamedTuple
- `K_coeffs`: K系数NamedTuple

# 返回值
总传播子（ComplexF64，单位 fm²）

# 示例
```julia
# 计算 u+u→u+u 散射的t道总传播子
D_t = total_propagator_auto(:uu_to_uu, :t, 
                            k0, k_norm, quark_params, thermo_params, K_coeffs)
```

# 错误处理
- 如果散射过程不存在：抛出错误并列出所有支持的散射过程
- 如果散射道不适用：抛出错误（如对qq散射请求s道）
"""
function total_propagator_auto(process::Symbol, channel::Symbol,
                               k0::Float64, k_norm::Float64,
                               quark_params::NamedTuple, thermo_params::NamedTuple,
                               K_coeffs::NamedTuple)
    # 验证散射过程是否存在
    if !haskey(SCATTERING_MESON_MAP, process)
        supported = join(sort(collect(keys(SCATTERING_MESON_MAP))), ", ")
        error("Unknown scattering process: $process. Supported processes: $supported")
    end
    
    process_info = SCATTERING_MESON_MAP[process]
    scattering_type = process_info[:type]
    
    # 验证散射道是否适用于该散射类型
    if !haskey(process_info[:channels], channel)
        if scattering_type == :qq
            error("Scattering process $process (qq type) only has t and u channels, not $channel")
        else
            error("Scattering process $process (qqbar type) only has t and s channels, not $channel")
        end
    end
    
    channel_info = process_info[:channels][channel]
    
    # 初始化总传播子
    D_total = ComplexF64(0.0, 0.0)
    
    # 1. 计算一般介子的贡献
    simple_mesons = channel_info[:simple]
    if !isempty(simple_mesons)
        D_simple = total_propagator_simple(
            process, channel, simple_mesons,
            k0, k_norm, quark_params, thermo_params, K_coeffs
        )
        D_total += D_simple
    end
    
    # 2. 计算赝标量混合介子（η/η'）的贡献
    if channel_info[:mixed_P]
        D_mixed_P = total_propagator_mixed(
            process, channel, :P,
            k0, k_norm, quark_params, thermo_params, K_coeffs
        )
        D_total += D_mixed_P
    end
    
    # 3. 计算标量混合介子（σ/σ'）的贡献
    if channel_info[:mixed_S]
        D_mixed_S = total_propagator_mixed(
            process, channel, :S,
            k0, k_norm, quark_params, thermo_params, K_coeffs
        )
        D_total += D_mixed_S
    end
    
    return D_total
end

"""
    calculate_all_propagators(process, k0, k_norm, 
                             quark_params, thermo_params, K_coeffs) -> NamedTuple

自动计算指定散射过程所有相关散射道的总传播子。

# 功能说明
根据散射类型自动判断需要计算的散射道：
- **qq散射**（夸克-夸克）：计算t道和u道
- **qqbar散射**（夸克-反夸克）：计算t道和s道

# 参数
- `process`: 散射过程符号（如 :uu_to_uu, :uubar_to_uubar）
- `k0`: 四动量能量分量（单位：fm⁻¹）
- `k_norm`: 三动量大小（单位：fm⁻¹）
- `quark_params`: 夸克参数NamedTuple
- `thermo_params`: 热力学参数NamedTuple
- `K_coeffs`: K系数NamedTuple

# 返回值
返回NamedTuple：
- **qq散射**：`(t = D_t::ComplexF64, u = D_u::ComplexF64)`
- **qqbar散射**：`(t = D_t::ComplexF64, s = D_s::ComplexF64)`

单位：fm²

# 示例
```julia
# 夸克-夸克散射（返回t和u道）
result = calculate_all_propagators(:uu_to_uu, 
                                   k0, k_norm, quark_params, thermo_params, K_coeffs)
println("D_t = ", result.t)
println("D_u = ", result.u)

# 夸克-反夸克散射（返回t和s道）
result = calculate_all_propagators(:uubar_to_uubar,
                                   k0, k_norm, quark_params, thermo_params, K_coeffs)
println("D_t = ", result.t)
println("D_s = ", result.s)
```

# 性能优化
该函数会自动利用极化函数缓存。对于对称的散射道（如uu→uu的t道和u道），
相同的极化函数只计算一次。

# 错误处理
如果散射过程不存在，抛出错误并列出所有支持的散射过程。
"""
function calculate_all_propagators(process::Symbol,
                                  k0::Float64, k_norm::Float64,
                                  quark_params::NamedTuple, thermo_params::NamedTuple,
                                  K_coeffs::NamedTuple)
    # 验证散射过程是否存在
    if !haskey(SCATTERING_MESON_MAP, process)
        supported = join(sort(collect(keys(SCATTERING_MESON_MAP))), ", ")
        error("Unknown scattering process: $process. Supported processes: $supported")
    end
    
    process_info = SCATTERING_MESON_MAP[process]
    scattering_type = process_info[:type]
    
    if scattering_type == :qq
        # 夸克-夸克散射：计算t道和u道
        D_t = total_propagator_auto(process, :t, k0, k_norm, quark_params, thermo_params, K_coeffs)
        D_u = total_propagator_auto(process, :u, k0, k_norm, quark_params, thermo_params, K_coeffs)
        return (t = D_t, u = D_u)
    else  # :qqbar
        # 夸克-反夸克散射：计算t道和s道
        D_t = total_propagator_auto(process, :t, k0, k_norm, quark_params, thermo_params, K_coeffs)
        D_s = total_propagator_auto(process, :s, k0, k_norm, quark_params, thermo_params, K_coeffs)
        return (t = D_t, s = D_s)
    end
end

"""
    calculate_all_propagators_by_channel(process, k0, k_norm, 
                                        quark_params, thermo_params, K_coeffs) -> NamedTuple

自动计算指定散射过程所有相关散射道的总传播子，并按标量(S)和赝标量(P)通道分离。

# 功能说明
根据散射类型自动判断需要计算的散射道，并将每个散射道的传播子按通道分离：
- **qq散射**（夸克-夸克）：计算t道和u道，各分离为S和P通道
- **qqbar散射**（夸克-反夸克）：计算t道和s道，各分离为S和P通道

# 通道分离逻辑
- **P通道（赝标量）**：累加π + K + η/η'混合介子的贡献
- **S通道（标量）**：累加σ_π + σ_K + σ/σ'混合介子的贡献

每个介子仅在其物理对应的通道中计算，无跨通道计算。

# 参数
- `process`: 散射过程符号（如 :uu_to_uu, :uubar_to_uubar）
- `k0`: 四动量能量分量（单位：fm⁻¹）
- `k_norm`: 三动量大小（单位：fm⁻¹）
- `quark_params`: 夸克参数NamedTuple
- `thermo_params`: 热力学参数NamedTuple
- `K_coeffs`: K系数NamedTuple

# 返回值
返回NamedTuple，按散射类型不同：
- **qq散射**：`(t_S = D_t^S, t_P = D_t^P, u_S = D_u^S, u_P = D_u^P)`
- **qqbar散射**：`(t_S = D_t^S, t_P = D_t^P, s_S = D_s^S, s_P = D_s^P)`

所有值为ComplexF64类型，单位：fm²

# 示例
```julia
# 夸克-夸克散射（返回t和u道的S/P分量）
result = calculate_all_propagators_by_channel(:uu_to_uu, 
                                              k0, k_norm, quark_params, thermo_params, K_coeffs)
println("D_t^S = ", result.t_S)  # t道标量通道
println("D_t^P = ", result.t_P)  # t道赝标量通道
println("D_u^S = ", result.u_S)  # u道标量通道
println("D_u^P = ", result.u_P)  # u道赝标量通道

# 夸克-反夸克散射（返回t和s道的S/P分量）
result = calculate_all_propagators_by_channel(:uubar_to_uubar,
                                              k0, k_norm, quark_params, thermo_params, K_coeffs)
println("D_t^S = ", result.t_S)  # t道标量通道
println("D_t^P = ", result.t_P)  # t道赝标量通道
println("D_s^S = ", result.s_S)  # s道标量通道
println("D_s^P = ", result.s_P)  # s道赝标量通道
```

# 性能优化
该函数会自动利用极化函数缓存。相同参数的极化函数只计算一次。

# 错误处理
如果散射过程不存在，抛出错误并列出所有支持的散射过程。
"""
function calculate_all_propagators_by_channel(process::Symbol,
                                              k0::Float64, k_norm::Float64,
                                              quark_params::NamedTuple, thermo_params::NamedTuple,
                                              K_coeffs::NamedTuple)
    # 验证散射过程是否存在
    if !haskey(SCATTERING_MESON_MAP, process)
        supported = join(sort(collect(keys(SCATTERING_MESON_MAP))), ", ")
        error("Unknown scattering process: $process. Supported processes: $supported")
    end
    
    process_info = SCATTERING_MESON_MAP[process]
    scattering_type = process_info[:type]
    
    if scattering_type == :qq
        # 夸克-夸克散射：计算t道和u道的S/P分量
        D_t_S, D_t_P = calculate_propagator_by_channel(process, :t, k0, k_norm, quark_params, thermo_params, K_coeffs)
        D_u_S, D_u_P = calculate_propagator_by_channel(process, :u, k0, k_norm, quark_params, thermo_params, K_coeffs)
        return (t_S = D_t_S, t_P = D_t_P, u_S = D_u_S, u_P = D_u_P)
    else  # :qqbar
        # 夸克-反夸克散射：计算t道和s道的S/P分量
        D_t_S, D_t_P = calculate_propagator_by_channel(process, :t, k0, k_norm, quark_params, thermo_params, K_coeffs)
        D_s_S, D_s_P = calculate_propagator_by_channel(process, :s, k0, k_norm, quark_params, thermo_params, K_coeffs)
        return (t_S = D_t_S, t_P = D_t_P, s_S = D_s_S, s_P = D_s_P)
    end
end

"""
    calculate_propagator_by_channel(process, channel, k0, k_norm,
                                    quark_params, thermo_params, K_coeffs) -> (ComplexF64, ComplexF64)

计算指定散射道的标量(S)和赝标量(P)通道传播子。

# 功能说明
根据散射过程和散射道，分别累加标量和赝标量通道的介子贡献：
- **P通道（赝标量）**：累加π、K介子 + η/η'混合介子
- **S通道（标量）**：累加σ_π、σ_K介子 + σ/σ'混合介子

# 参数
- `process`: 散射过程符号（如 :uu_to_uu）
- `channel`: 散射道类型（:t, :u, :s）
- `k0`: 四动量能量分量（单位：fm⁻¹）
- `k_norm`: 三动量大小（单位：fm⁻¹）
- `quark_params`: 夸克参数NamedTuple
- `thermo_params`: 热力学参数NamedTuple
- `K_coeffs`: K系数NamedTuple

# 返回值
返回元组 (D_S, D_P)，分别为标量和赝标量通道的传播子（ComplexF64，单位：fm²）

# 实现逻辑
1. 从SCATTERING_MESON_MAP获取该散射道的介子列表和混合介子标记
2. 根据介子类型分类：
   - 赝标量介子（:pi, :K）→ 累加到P通道
   - 标量介子（:sigma_pi, :sigma_K）→ 累加到S通道
3. 如果存在混合介子：
   - mixed_P=true → 计算η/η'贡献并累加到P通道
   - mixed_S=true → 计算σ/σ'贡献并累加到S通道
4. 返回两个通道的累加结果
"""
function calculate_propagator_by_channel(process::Symbol, channel::Symbol,
                                        k0::Float64, k_norm::Float64,
                                        quark_params::NamedTuple, thermo_params::NamedTuple,
                                        K_coeffs::NamedTuple)
    # 获取散射道信息
    process_info = SCATTERING_MESON_MAP[process]
    channel_info = process_info[:channels][channel]
    
    # 初始化S和P通道的传播子
    D_S = ComplexF64(0.0, 0.0)
    D_P = ComplexF64(0.0, 0.0)
    
    # 获取味因子
    T1, T2 = get_flavor_factors_for_channel(process, channel)
    
    # 1. 处理一般介子
    simple_mesons = channel_info[:simple]
    for meson in simple_mesons
        # 获取该介子的极化函数参数
        pol_params = get_polarization_params(meson, quark_params)
        
        # 调用带缓存的极化函数
        Π_real, Π_imag = polarization_aniso_cached(
            pol_params.channel, k0, k_norm,
            pol_params.m1, pol_params.m2,
            pol_params.μ1, pol_params.μ2,
            thermo_params.T, thermo_params.Φ, thermo_params.Φbar, thermo_params.ξ,
            pol_params.A1, pol_params.A2,
            pol_params.num_s_quark
        )
        
        Π = ComplexF64(Π_real, Π_imag)
        
        # 计算介子传播子
        D_meson = meson_propagator_simple(meson, K_coeffs, Π)
        
        # 根据介子类型累加到对应通道
        if meson in [:pi, :K]
            # 赝标量介子 → P通道
            D_P += T1 * D_meson * T2
        elseif meson in [:sigma_pi, :sigma_K]
            # 标量介子 → S通道
            D_S += T1 * D_meson * T2
        else
            error("Unknown meson type in channel separation: $meson")
        end
    end
    
    # 2. 处理赝标量混合介子（η/η'）
    if channel_info[:mixed_P]
        D_mixed_P = total_propagator_mixed(
            process, channel, :P,
            k0, k_norm, quark_params, thermo_params, K_coeffs
        )
        D_P += D_mixed_P
    end
    
    # 3. 处理标量混合介子（σ/σ'）
    if channel_info[:mixed_S]
        D_mixed_S = total_propagator_mixed(
            process, channel, :S,
            k0, k_norm, quark_params, thermo_params, K_coeffs
        )
        D_S += D_mixed_S
    end
    
    return (D_S, D_P)
end

end # module TotalPropagator
