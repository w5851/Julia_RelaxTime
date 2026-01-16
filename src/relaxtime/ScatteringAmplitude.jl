module ScatteringAmplitude

"""
# ScatteringAmplitude.jl

散射矩阵元计算模块，支持PNJL模型中的夸克弹性散射过程。

在三味PNJL模型中，夸克-夸克、夸克-反夸克弹性散射通过标量和赝标量介子交换
传递相互作用。本模块实现从总传播子到散射矩阵元平方|M|²的计算，包含色和
自旋平均因子1/(4N_c²)。

公式参考 doc/formula/ScatteringAmplitude_散射矩阵元by总传播子.md

## 核心设计原则
- 支持所有11个散射过程（4个qq散射 + 7个qqbar散射）
- 直接计算|M|²而无需单独构建振幅M
- 按标量(S)和赝标量(P)通道分离传播子
- 自动计算质心系介子四动量
- 正确处理不同过程的Mandelstam变量组合
"""

include("../Constants_PNJL.jl")
include("../utils/ParticleSymbols.jl")
include("TotalPropagator.jl")

using .Constants_PNJL: N_color, SCATTERING_MESON_MAP
using .TotalPropagator: calculate_cms_momentum
using .ParticleSymbols: get_quark_masses_for_process

export scattering_amplitude_squared
export calculate_mandelstam_variables
export calculate_all_scattering_amplitudes_squared

# 使用N_color作为N_c
const N_c = N_color

# ----------------------------------------------------------------------------
# 辅助函数：Mandelstam辅助变量计算
# ----------------------------------------------------------------------------

"""
    calculate_mandelstam_variables(s, t, u, m1, m2, m3, m4) -> NamedTuple

计算散射矩阵元所需的18个Mandelstam辅助变量。

# 辅助变量定义
对于散射过程 q₁ + q₂ → q₃ + q₄，定义：
- s_ij^± = s - (m_i ± m_j)²
- t_ij^± = t - (m_i ± m_j)²
- u_ij^± = u - (m_i ± m_j)²

其中 +/- 符号对应质量的和或差。

# 参数
- `s`: Mandelstam变量 s（单位：fm⁻²）
- `t`: Mandelstam变量 t（单位：fm⁻²）
- `u`: Mandelstam变量 u（单位：fm⁻²）
- `m1`, `m2`, `m3`, `m4`: 四个粒子的质量（单位：fm⁻¹）

# 返回值
返回NamedTuple包含18个辅助变量（单位：fm⁻²）：
```julia
(s_12_plus, s_12_minus, s_34_plus, s_34_minus,
 t_13_plus, t_13_minus, t_24_plus, t_24_minus,
 u_14_plus, u_14_minus, u_23_plus, u_23_minus,
 t_14_plus, t_14_minus, t_23_plus, t_23_minus,
 u_13_plus, u_13_minus)
```

# Mandelstam约束
输入应满足：s + t + u = m₁² + m₂² + m₃² + m₄²

# 示例
```julia
s = 4.0  # fm⁻²
t = -0.5
u = 0.36 - s - t  # 假设所有质量为0.3 fm⁻¹
m = 0.3  # fm⁻¹

vars = calculate_mandelstam_variables(s, t, u, m, m, m, m)
println("s_12^+ = ", vars.s_12_plus)
println("t_13^- = ", vars.t_13_minus)
```
"""
function calculate_mandelstam_variables(s::Float64, t::Float64, u::Float64,
                                       m1::Float64, m2::Float64, m3::Float64, m4::Float64)
    # s 相关变量
    s_12_plus = s - (m1 + m2)^2
    s_12_minus = s - (m1 - m2)^2
    s_34_plus = s - (m3 + m4)^2
    s_34_minus = s - (m3 - m4)^2
    
    # t 相关变量
    t_13_plus = t - (m1 + m3)^2
    t_13_minus = t - (m1 - m3)^2
    t_24_plus = t - (m2 + m4)^2
    t_24_minus = t - (m2 - m4)^2
    t_14_plus = t - (m1 + m4)^2
    t_14_minus = t - (m1 - m4)^2
    t_23_plus = t - (m2 + m3)^2
    t_23_minus = t - (m2 - m3)^2
    
    # u 相关变量
    u_14_plus = u - (m1 + m4)^2
    u_14_minus = u - (m1 - m4)^2
    u_23_plus = u - (m2 + m3)^2
    u_23_minus = u - (m2 - m3)^2
    u_13_plus = u - (m1 + m3)^2
    u_13_minus = u - (m1 - m3)^2
    
    return (
        s_12_plus = s_12_plus, s_12_minus = s_12_minus,
        s_34_plus = s_34_plus, s_34_minus = s_34_minus,
        t_13_plus = t_13_plus, t_13_minus = t_13_minus,
        t_24_plus = t_24_plus, t_24_minus = t_24_minus,
        u_14_plus = u_14_plus, u_14_minus = u_14_minus,
        u_23_plus = u_23_plus, u_23_minus = u_23_minus,
        t_14_plus = t_14_plus, t_14_minus = t_14_minus,
        t_23_plus = t_23_plus, t_23_minus = t_23_minus,
        u_13_plus = u_13_plus, u_13_minus = u_13_minus
    )
end

# ----------------------------------------------------------------------------
# 主函数：散射矩阵元平方
# ----------------------------------------------------------------------------

"""
    scattering_amplitude_squared(process, s, t, quark_params, thermo_params, K_coeffs) -> Float64

计算散射矩阵元平方 |M|²（包含色和自旋平均因子 1/(4N_c²)）。

# 物理背景
在三味PNJL模型中，夸克弹性散射通过介子交换实现。散射矩阵元平方是计算
散射截面和弛豫时间的核心输入。

# 计算流程
1. 根据散射过程提取四个粒子质量
2. 计算Mandelstam变量 u = Σm² - s - t
3. 计算质心系介子四动量 (k0, k) 用于传播子计算
4. 获取通道分离的传播子（S和P通道）
5. 计算18个Mandelstam辅助变量
6. 根据散射类型（qq或qqbar）应用相应公式计算|M|²

# 参数
- `process`: 散射过程符号（如 :uu_to_uu, :uubar_to_ssbar）
- `s`: Mandelstam变量 s（单位：fm⁻²）
- `t`: Mandelstam变量 t（单位：fm⁻²）
- `quark_params`: 夸克参数NamedTuple，结构：
  ```julia
  (m = (u=m_u, d=m_d, s=m_s),   # 夸克质量（fm⁻¹）
   μ = (u=μ_u, d=μ_d, s=μ_s),   # 夸克化学势（fm⁻¹）
   A = (u=A_u, d=A_u, s=A_s))   # A函数值（fm）
  ```
- `thermo_params`: 热力学参数NamedTuple，结构：
  ```julia
  (T=T, Φ=Φ, Φbar=Φbar, ξ=ξ)    # 温度（fm⁻¹）、Polyakov环、各向异性参数
  ```
- `K_coeffs`: K系数NamedTuple（通过`EffectiveCouplings.calculate_effective_couplings`获取）

# 返回值
散射矩阵元平方 |M|²（单位：fm⁻⁴，包含色自旋平均因子1/(4N_c²)）

# 散射类型与公式
- **qq散射**（4种）：|M|² = |M_u|² + |M_t|² - 2Re(M_u M_t*)
  - 涉及t道和u道传播子
  - 过程：:uu_to_uu, :ss_to_ss, :ud_to_ud, :us_to_us

- **qqbar散射**（7种）：|M|² = |M_s|² + |M_t|² - 2Re(M_s M_t*)
  - 涉及s道和t道传播子
  - 过程：:udbar_to_udbar, :usbar_to_usbar, :uubar_to_uubar, 
          :uubar_to_ddbar, :uubar_to_ssbar, :ssbar_to_uubar, :ssbar_to_ssbar

# Mandelstam约束
自动计算 u = m₁² + m₂² + m₃² + m₄² - s - t 以满足约束条件。

# 示例
```julia
using .EffectiveCouplings: calculate_effective_couplings, calculate_G_from_A
using .OneLoopIntegrals: A
using .Constants_PNJL

# 设置参数
T = 150.0 / 197.327  # 150 MeV → fm⁻¹
m_u = 300.0 / 197.327
m_s = 500.0 / 197.327
μ_u = 0.0
μ_s = 0.0
Φ = 0.5
Φbar = 0.5
ξ = 0.0

# 计算夸克参数
A_u = A(T, μ_u, m_u, Φ, Φbar)
A_s = A(T, μ_s, m_s, Φ, Φbar)
G_u = calculate_G_from_A(A_u)
G_s = calculate_G_from_A(A_s)

quark_params = (
    m = (u=m_u, d=m_u, s=m_s),
    μ = (u=μ_u, d=μ_u, s=μ_s),
    A = (u=A_u, d=A_u, s=A_s)
)
thermo_params = (T=T, Φ=Φ, Φbar=Φbar, ξ=ξ)

# 计算K系数
G_fm2 = Constants_PNJL.G_GeV_inv2 / (197.327^2)
K_fm5 = Constants_PNJL.K_GeV_inv5 / (197.327^5)
K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

# 计算散射矩阵元平方
s = 4.0  # fm⁻²
t = -0.5  # fm⁻²
M_squared = scattering_amplitude_squared(
    :uu_to_uu, s, t, quark_params, thermo_params, K_coeffs
)
println("|M|² = ", M_squared, " fm⁻⁴")
```

# 注意事项
1. 确保Mandelstam变量满足物理约束（s, t > 0 for qq; s > 0, t < 0 for qqbar）
2. 传播子的虚部会导致|M|²为复数，但物理上应取实部
3. 温度T和化学势μ必须与quark_params中的A函数值一致
4. 各向异性参数ξ影响极化函数，进而影响传播子和散射矩阵元
"""
function scattering_amplitude_squared(process::Symbol, s::Float64, t::Float64,
                                     quark_params::NamedTuple, thermo_params::NamedTuple,
                                     K_coeffs::NamedTuple)
    # 1. 提取四个粒子质量
    m1, m2, m3, m4 = get_quark_masses_for_process(process, quark_params)
    
    # 2. 计算 u（Mandelstam约束）
    u = m1^2 + m2^2 + m3^2 + m4^2 - s - t
    
    # 3. 获取散射过程信息
    if !haskey(SCATTERING_MESON_MAP, process)
        supported = join(sort(collect(keys(SCATTERING_MESON_MAP))), ", ")
        error("Unknown scattering process: $process. Supported processes: $supported")
    end
    
    process_info = SCATTERING_MESON_MAP[process]
    scattering_type = process_info[:type]
    
    # 4. 根据散射类型计算|M|²
    if scattering_type == :qq
        return calculate_qq_amplitude_squared(
            process, s, t, u, m1, m2, m3, m4,
            quark_params, thermo_params, K_coeffs
        )
    else  # :qqbar
        return calculate_qqbar_amplitude_squared(
            process, s, t, u, m1, m2, m3, m4,
            quark_params, thermo_params, K_coeffs
        )
    end
end

# ----------------------------------------------------------------------------
# 内部函数：qq散射矩阵元平方
# ----------------------------------------------------------------------------

"""
    calculate_qq_amplitude_squared(process, s, t, u, m1, m2, m3, m4,
                                   quark_params, thermo_params, K_coeffs) -> Float64

计算夸克-夸克散射的矩阵元平方。

# 公式（包含色自旋平均）
|M|² = |M_u|² + |M_t|² - 2Re(M_u M_t*)

其中：
- |M_u|² = |D_u^S|² u_14^+ u_23^+ + |D_u^P|² u_14^- u_23^-
- |M_t|² = |D_t^S|² t_13^+ t_24^+ + |D_t^P|² t_13^- t_24^-
- 交叉项包含4个D_t^(S/P) D_u^(S/P)* 的组合

# 计算流程
1. 计算t道质心系动量 (k0_t, k_t)
2. 计算u道质心系动量 (k0_u, k_u)
3. 获取t道和u道的S/P通道传播子
4. 计算Mandelstam辅助变量
5. 应用公式计算|M|²
"""
function calculate_qq_amplitude_squared(process::Symbol, s::Float64, t::Float64, u::Float64,
                                       m1::Float64, m2::Float64, m3::Float64, m4::Float64,
                                       quark_params::NamedTuple, thermo_params::NamedTuple,
                                       K_coeffs::NamedTuple)
    # 1. 计算质心系动量
    cms_t = calculate_cms_momentum(process, s, t, :t, quark_params; u=u)
    cms_u = calculate_cms_momentum(process, s, t, :u, quark_params; u=u)
    
    # 2. 获取t道传播子（按需只算 t 道，避免同时计算 u 道的浪费）
    D_t_S, D_t_P = TotalPropagator.calculate_propagator_by_channel(
        process, :t, cms_t.k0, cms_t.k, quark_params, thermo_params, K_coeffs
    )
    
    # 3. 获取u道传播子（按需只算 u 道）
    D_u_S, D_u_P = TotalPropagator.calculate_propagator_by_channel(
        process, :u, cms_u.k0, cms_u.k, quark_params, thermo_params, K_coeffs
    )
    
    # 4. 计算Mandelstam辅助变量
    vars = calculate_mandelstam_variables(s, t, u, m1, m2, m3, m4)
    
    # 5. 计算色自旋平均后的|M_u|² (根据公式已包含1/(4N_c²)平均因子)
    # 1/(4N_c²) ∑|M_u|² = |D_u^S|² u₁₄⁺ u₂₃⁺ + |D_u^P|² u₁₄⁻ u₂₃⁻
    M_u_squared = abs2(D_u_S) * vars.u_14_plus * vars.u_23_plus + 
                  abs2(D_u_P) * vars.u_14_minus * vars.u_23_minus
    
    # 6. 计算色自旋平均后的|M_t|² (根据公式已包含1/(4N_c²)平均因子)
    # 1/(4N_c²) ∑|M_t|² = |D_t^S|² t₁₃⁺ t₂₄⁺ + |D_t^P|² t₁₃⁻ t₂₄⁻
    M_t_squared = abs2(D_t_S) * vars.t_13_plus * vars.t_24_plus + 
                  abs2(D_t_P) * vars.t_13_minus * vars.t_24_minus
    
    # 7. 计算色自旋平均后的交叉项
    # 公式：1/(4N_c²) ∑(M_u M_t*) = 1/(4N_c) × [4个D组合项]
    # 左右两边相等，即右边已经是完整的 1/(4N_c²) ∑(M_u M_t*)
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
    
    # 8. 总矩阵元平方
    # 公式：|M|²_total = [1/(4N_c²)∑|M_u|²] + [1/(4N_c²)∑|M_t|²] - 2Re[1/(4N_c²)∑(M_u M_t*)]
    # 注意：cross_term已经通过cross_term_factor = 1/(4N_c)完整计算了1/(4N_c²)∑(M_u M_t*)
    M_squared = M_u_squared + M_t_squared - 2.0 * real(cross_term)
    
    return M_squared
end

# ----------------------------------------------------------------------------
# 内部函数：qqbar散射矩阵元平方
# ----------------------------------------------------------------------------

"""
    calculate_qqbar_amplitude_squared(process, s, t, u, m1, m2, m3, m4,
                                      quark_params, thermo_params, K_coeffs) -> Float64

计算夸克-反夸克散射的矩阵元平方。

# 公式（包含色自旋平均）
|M|² = |M_s|² + |M_t|² - 2Re(M_s M_t*)

其中：
- |M_s|² = |D_s^S|² s_12^+ s_34^+ + |D_s^P|² s_12^- s_34^-
- |M_t|² = |D_t^S|² t_13^+ t_24^+ + |D_t^P|² t_13^- t_24^-
- 交叉项包含4个D_s^(S/P) D_t^(S/P)* 的组合

# 计算流程
1. 计算s道质心系动量 (k0_s, k_s)
2. 计算t道质心系动量 (k0_t, k_t)
3. 获取s道和t道的S/P通道传播子
4. 计算Mandelstam辅助变量
5. 应用公式计算|M|²
"""
function calculate_qqbar_amplitude_squared(process::Symbol, s::Float64, t::Float64, u::Float64,
                                          m1::Float64, m2::Float64, m3::Float64, m4::Float64,
                                          quark_params::NamedTuple, thermo_params::NamedTuple,
                                          K_coeffs::NamedTuple)
    # 1. 计算质心系动量
    cms_s = calculate_cms_momentum(process, s, t, :s, quark_params; u=u)
    cms_t = calculate_cms_momentum(process, s, t, :t, quark_params; u=u)
    
    # 2. 获取s道传播子（按需只算 s 道）
    D_s_S, D_s_P = TotalPropagator.calculate_propagator_by_channel(
        process, :s, cms_s.k0, cms_s.k, quark_params, thermo_params, K_coeffs
    )
    
    # 3. 获取t道传播子（按需只算 t 道）
    D_t_S, D_t_P = TotalPropagator.calculate_propagator_by_channel(
        process, :t, cms_t.k0, cms_t.k, quark_params, thermo_params, K_coeffs
    )
    
    # 4. 计算Mandelstam辅助变量
    vars = calculate_mandelstam_variables(s, t, u, m1, m2, m3, m4)
    
    # 5. 计算色自旋平均后的|M_s|² (根据公式已包含1/(4N_c²)平均因子)
    # 1/(4N_c²) ∑|M_s|² = |D_s^S|² s₁₂⁺ s₃₄⁺ + |D_s^P|² s₁₂⁻ s₃₄⁻
    M_s_squared = abs2(D_s_S) * vars.s_12_plus * vars.s_34_plus + 
                  abs2(D_s_P) * vars.s_12_minus * vars.s_34_minus
    
    # 6. 计算色自旋平均后的|M_t|² (根据公式已包含1/(4N_c²)平均因子)
    # 1/(4N_c²) ∑|M_t|² = |D_t^S|² t₁₃⁺ t₂₄⁺ + |D_t^P|² t₁₃⁻ t₂₄⁻
    M_t_squared = abs2(D_t_S) * vars.t_13_plus * vars.t_24_plus + 
                  abs2(D_t_P) * vars.t_13_minus * vars.t_24_minus
    
    # 7. 计算色自旋平均后的交叉项
    # 公式：1/(4N_c²) ∑(M_s M_t*) = 1/(4N_c) × [4个D组合项]
    # 左右两边相等，即右边已经是完整的 1/(4N_c²) ∑(M_s M_t*)
    cross_term_factor = 1.0 / (4.0 * N_c)
    
    term1 = D_s_S * conj(D_t_S) * (vars.s_12_plus * vars.s_34_plus - 
                                    vars.u_14_plus * vars.u_23_plus + 
                                    vars.t_13_plus * vars.t_24_plus)
    
    term2 = D_s_S * conj(D_t_P) * (vars.s_12_plus * vars.s_34_plus - 
                                    vars.u_14_minus * vars.u_23_minus + 
                                    vars.t_13_minus * vars.t_24_minus)
    
    term3 = D_s_P * conj(D_t_S) * (vars.s_12_minus * vars.s_34_minus - 
                                    vars.u_14_minus * vars.u_23_minus + 
                                    vars.t_13_plus * vars.t_24_plus)
    
    term4 = D_s_P * conj(D_t_P) * (vars.s_12_minus * vars.s_34_minus - 
                                    vars.u_14_plus * vars.u_23_plus + 
                                    vars.t_13_minus * vars.t_24_minus)
    
    cross_term = cross_term_factor * (term1 - term2 - term3 + term4)
    
    # 8. 总矩阵元平方
    # 公式：|M|²_total = [1/(4N_c²)∑|M_s|²] + [1/(4N_c²)∑|M_t|²] - 2Re[1/(4N_c²)∑(M_s M_t*)]
    # 注意：cross_term已经通过cross_term_factor = 1/(4N_c)完整计算了1/(4N_c²)∑(M_s M_t*)
    M_squared = M_s_squared + M_t_squared - 2.0 * real(cross_term)
    
    return M_squared
end

# ----------------------------------------------------------------------------
# 批量计算函数
# ----------------------------------------------------------------------------

"""
    calculate_all_scattering_amplitudes_squared(s, t, quark_params, thermo_params, K_coeffs) -> NamedTuple

批量计算所有13种散射过程的散射矩阵元平方 |M|²。

这是一个便利函数，一次性计算所有支持的散射过程，返回包含每个过程结果的NamedTuple。
适用于需要同时分析多个散射过程的场景，如计算总散射率或弛豫时间。

# 参数
- `s::Float64`: Mandelstam变量s（质心系能量平方，单位：fm⁻²）
- `t::Float64`: Mandelstam变量t（动量转移平方，单位：fm⁻²）
- `quark_params::NamedTuple`: 夸克参数，包含m(质量), μ(化学势), A(单圈积分值)
- `thermo_params::NamedTuple`: 热力学参数，包含T(温度), Φ(Polyakov环), Φbar(共轭Polyakov环), ξ(各向异性参数)
- `K_coeffs::NamedTuple`: 有效耦合常数，由`calculate_effective_couplings`计算

# 返回值
NamedTuple，包含所有13种散射过程的|M|²值（单位：fm⁻⁴）：

**夸克-夸克散射（4种）：**
- `uu_to_uu`: u+u→u+u
- `ss_to_ss`: s+s→s+s
- `ud_to_ud`: u+d→u+d
- `us_to_us`: u+s→u+s

**夸克-反夸克散射（9种）：**
- `udbar_to_udbar`: u+đ→u+đ
- `usbar_to_usbar`: u+s̄→u+s̄
- `dubar_to_dubar`: d+ū→d+ū （与udbar等价）
- `subar_to_subar`: s+ū→s+ū （与usbar等价）
- `uubar_to_uubar`: u+ū→u+ū
- `uubar_to_ddbar`: u+ū→d+đ
- `uubar_to_ssbar`: u+ū→s+s̄
- `ssbar_to_uubar`: s+s̄→u+ū
- `ssbar_to_ssbar`: s+s̄→s+s̄

# 物理约束
- 所有返回值均满足 |M|² ≥ 0（物理必要条件）
- dubar_to_dubar 与 udbar_to_udbar 在同位旋对称下数值相同
- subar_to_subar 与 usbar_to_usbar 在同位旋对称下数值相同

# 性能说明
- 每个过程独立计算，无相互依赖
- 利用极化函数缓存机制加速重复计算
- 典型耗时：~0.5-1ms（13个过程总计，取决于缓存命中率）

# 示例
```julia
# 准备参数
s = 8.0  # fm⁻²
t = -0.3  # fm⁻²
quark_params = (
    m = (u=0.3, d=0.3, s=0.5),
    μ = (u=0.0, d=0.0, s=0.0),
    A = (u=0.05, d=0.05, s=0.08)
)
thermo_params = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

# 批量计算
results = calculate_all_scattering_amplitudes_squared(s, t, quark_params, thermo_params, K_coeffs)

# 访问结果
println("uu→uu: |M|² = ", results.uu_to_uu, " fm⁻⁴")
println("ss→ss: |M|² = ", results.ss_to_ss, " fm⁻⁴")

# 计算总散射率（示例）
total_rate = sum(values(results))
```

# 参见
- `scattering_amplitude_squared`: 单个散射过程的计算
- `calculate_effective_couplings`: 有效耦合常数计算
"""
function calculate_all_scattering_amplitudes_squared(
    s::Float64, 
    t::Float64,
    quark_params::NamedTuple, 
    thermo_params::NamedTuple,
    K_coeffs::NamedTuple
)::NamedTuple
    
    # 获取所有散射过程（按固定顺序）
    all_processes = [
        # qq散射（4种）
        :uu_to_uu,
        :ss_to_ss,
        :ud_to_ud,
        :us_to_us,
        # qqbar散射（9种）
        :udbar_to_udbar,
        :usbar_to_usbar,
        :dubar_to_dubar,
        :subar_to_subar,
        :uubar_to_uubar,
        :uubar_to_ddbar,
        :uubar_to_ssbar,
        :ssbar_to_uubar,
        :ssbar_to_ssbar
    ]
    
    # 计算每个过程的矩阵元平方
    results = Dict{Symbol, Float64}()
    for process in all_processes
        results[process] = scattering_amplitude_squared(
            process, s, t, quark_params, thermo_params, K_coeffs
        )
    end
    
    # 转换为NamedTuple（保持固定顺序）
    return NamedTuple(results)
end

end # module ScatteringAmplitude
