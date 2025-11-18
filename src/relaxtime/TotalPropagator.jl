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
export calculate_all_propagators, total_propagator_auto
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

function parse_scattering_process(process::Symbol)
    s = string(process)
    parts = split(s, "_to_")
    if length(parts) != 2
        error("Invalid process format: $process. Expected: q1q2_to_q3q4")
    end
    
    initial = parts[1]
    final = parts[2]
    
    # 解析初态
    if occursin("bar", initial)
        initial_clean = replace(initial, "bar" => "")
        if length(initial_clean) == 2 && endswith(initial, "bar")
            q1 = Symbol(initial_clean[1:1])
            q2 = Symbol(initial_clean[2:2] * "bar")
        else
            error("Invalid initial state: $initial")
        end
    else
        if length(initial) != 2
            error("Invalid initial state: $initial")
        end
        q1 = Symbol(initial[1:1])
        q2 = Symbol(initial[2:2])
    end
    
    # 解析末态
    if occursin("bar", final)
        final_clean = replace(final, "bar" => "")
        if length(final_clean) == 2 && endswith(final, "bar")
            q3 = Symbol(final_clean[1:1])
            q4 = Symbol(final_clean[2:2] * "bar")
        else
            error("Invalid final state: $final")
        end
    else
        if length(final) != 2
            error("Invalid final state: $final")
        end
        q3 = Symbol(final[1:1])
        q4 = Symbol(final[2:2])
    end
    
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

end # module TotalPropagator
