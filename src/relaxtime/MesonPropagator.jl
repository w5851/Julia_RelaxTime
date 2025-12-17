"""
# MesonPropagator.jl

介子传播子计算模块，支持一般介子（π、K、σ_π、σ_K）和混合介子（η/η'、σ/σ'）的完整计算。

公式参考 `doc/formula/Propagator_传播子byPolarization.md`（第8节）。

## 核心设计原则
- 所有函数接收预计算的极化函数Π和K系数作为参数
- 各向异性修正通过预计算的Π值间接影响传播子
- 批量计算时可复用同一组K系数以提升性能
"""
module MesonPropagator

include("../Constants_PNJL.jl")
include("EffectiveCouplings.jl")

using .Constants_PNJL: λ₀, λ₈, ψ_u, ψ_d, ψ_s, ψbar_u, ψbar_d, ψbar_s
using .EffectiveCouplings: coupling_matrix_determinant
using LinearAlgebra: det

# 传播子分母正则：避免在共振/极点附近出现 NaN/Inf。
# 注意：这不是物理“宽度”建模，只是数值稳健性保护；值应尽量小。
const PROPAGATOR_DENOM_EPS = 1e-12

@inline function _isfinite_complex(z::ComplexF64)::Bool
    return isfinite(real(z)) && isfinite(imag(z))
end

export meson_propagator_simple, meson_propagator_mixed
export calculate_coupling_matrix, extract_flavor, get_quark_wavefunction, calculate_current_vector

# ----------------------------------------------------------------------------
# 辅助函数：味提取和波函数选择
# ----------------------------------------------------------------------------

"""
    extract_flavor(q::Symbol) -> (Symbol, Bool)

从夸克/反夸克符号中提取味类型和反粒子标志。

# 参数
- `q`: 夸克/反夸克类型（`:u`, `:d`, `:s`, `:ubar`, `:dbar`, `:sbar`）

# 返回值
返回元组 `(flavor, is_bar)`：
- `flavor`: 味类型（`:u`, `:d`, `:s`）
- `is_bar`: 是否为反夸克（`true`表示反夸克）

# 示例
```julia
extract_flavor(:u) # 返回 (:u, false)
extract_flavor(:ubar) # 返回 (:u, true)
```
"""
@inline function extract_flavor(q::Symbol)
    s = string(q)
    if endswith(s, "bar")
        flavor = Symbol(s[1:end-3])  # 去除"bar"后缀
        return (flavor, true)
    else
        return (q, false)
    end
end

"""
    get_quark_wavefunction(flavor::Symbol, is_bar::Bool) -> Union{Vector{Float64}, Matrix{Float64}}

根据味类型和反粒子标志返回对应的夸克波函数。

# 参数
- `flavor`: 味类型（`:u`, `:d`, `:s`）
- `is_bar`: 是否为反夸克（`false`为夸克，`true`为反夸克）

# 返回值
- `is_bar=false`: 返回列向量ψ（3×1向量）
- `is_bar=true`: 返回行向量ψ̄（1×3矩阵）

# 物理意义
- ψ_u = [1, 0, 0]：u夸克在味空间的表示
- ψbar_u = [1 0 0]：u反夸克在味空间的表示（1×3矩阵）

# 示例
```julia
ψu = get_quark_wavefunction(:u, false)  # 返回列向量 [1.0, 0.0, 0.0]
ψbar_u = get_quark_wavefunction(:u, true)  # 返回行向量 [1.0 0.0 0.0]
```
"""
@inline function get_quark_wavefunction(flavor::Symbol, is_bar::Bool)
    if is_bar
        if flavor == :u
            return ψbar_u
        elseif flavor == :d
            return ψbar_d
        elseif flavor == :s
            return ψbar_s
        else
            error("Unknown flavor: $flavor")
        end
    else
        if flavor == :u
            return ψ_u
        elseif flavor == :d
            return ψ_d
        elseif flavor == :s
            return ψ_s
        else
            error("Unknown flavor: $flavor")
        end
    end
end

# ----------------------------------------------------------------------------
# 辅助函数：流算符向量计算
# ----------------------------------------------------------------------------

"""
    calculate_current_vector(q1::Symbol, q2::Symbol, channel::Symbol) -> Vector{Float64}

根据散射道映射表自动选择场算符，计算流算符向量。

# 散射道映射表
| 散射道 | ψ | ψ̄ | 公式 |
|-------|---|---|------|
| t道   | q1 | q2 | J = [ψ̄₂·λ₀·ψ₁, ψ̄₂·λ₈·ψ₁] |
| s道   | q1 | q2 | J = [ψ̄₂·λ₀·ψ₁, ψ̄₂·λ₈·ψ₁] |
| u道   | q1 | q2 | J = [ψ̄₂·λ₀·ψ₁, ψ̄₂·λ₈·ψ₁] |

注：散射道映射在调用前已通过channel确定q1和q2的对应关系。

# 参数
- `q1`: 第一个场算符对应的夸克/反夸克（`:u`, `:d`, `:s`, `:ubar`, `:dbar`, `:sbar`）
- `q2`: 第二个场算符对应的夸克/反夸克
- `channel`: 散射道类型（`:t`, `:s`, `:u`）

# 返回值
返回2×1列向量 J = [J₀, J₈]

# 物理意义
流算符向量编码了介子的味量子数：
- J₀：单态分量（与λ₀对应）
- J₈：八重态分量（与λ₈对应）

# 示例
```julia
# u+d→u+d 散射（t道）
J = calculate_current_vector(:u, :d, :t)
```
"""
function calculate_current_vector(q1::Symbol, q2::Symbol, channel::Symbol)
    # 提取味和反粒子标志
    flavor1, is_bar1 = extract_flavor(q1)
    flavor2, is_bar2 = extract_flavor(q2)
    
    # 计算流算符向量 J = [ψ̄·λ₀·ψ, ψ̄·λ₈·ψ]
    # 根据散射道映射表：
    # - 第一个粒子(q1)对应ψ（列向量）
    # - 第二个粒子(q2)对应ψ̄（行向量）
    # 无论是夸克还是反夸克，都使用相同味道的波函数
    
    # 获取波函数（总是用相同的味）
    ψ1 = get_quark_wavefunction(flavor1, false)  # 列向量
    ψbar2 = get_quark_wavefunction(flavor2, true)  # 行向量
    
    # 矩阵乘法顺序：ψbar2 (1×3) * λ (3×3) * ψ1 (3×1) = 标量
    J0 = (ψbar2 * λ₀ * ψ1)[1]  # 提取标量值
    J8 = (ψbar2 * λ₈ * ψ1)[1]
    
    return [J0, J8]
end

# ----------------------------------------------------------------------------
# 辅助函数：耦合矩阵计算
# ----------------------------------------------------------------------------

"""
    calculate_coupling_matrix(Π_uu::ComplexF64, Π_ss::ComplexF64, K_coeffs::NamedTuple, 
                              channel::Symbol) -> Matrix{ComplexF64}

计算混合介子传播子所需的2×2复数耦合矩阵M。

# 物理背景
在η-η'（赝标量）或σ-σ'（标量）混合介子系统中，耦合矩阵M描述了SU(3)单态和八重态的耦合。
矩阵元素由有效耦合系数K和极化函数Π共同决定。

# 公式
先从 `K_coeffs` 取出 (k00,k88,k08) 并计算
detK = k00*k88 - k08^2

再构造 2×2 复数矩阵 M（实现里用 m00/m08/m88 记号）：
- m00 = k00 - (4/3) * detK * (Π_uu + 2Π_ss)
- m08 = m80 = k08 + (4/3)*sqrt(2) * detK * (Π_uu - Π_ss)
- m88 = k88 - (4/3) * detK * (2Π_uu + Π_ss)

其中±表示通道类型：
- +:赝标量通道（P通道），用于η/η'
- -:标量通道（S通道），用于σ/σ'

# 参数
- `Π_uu`: u/d夸克的极化函数（单位：fm²）
- `Π_ss`: s夸克的极化函数（单位：fm²）
- `K_coeffs`: 预计算的K系数NamedTuple（通过`EffectiveCouplings.calculate_effective_couplings`获取）
- `channel`: 通道类型（`:P`为赝标量，`:S`为标量）

# 返回值
返回2×2复数矩阵M（单位：无量纲）

# 注意事项
1. K系数应预先通过`EffectiveCouplings.calculate_effective_couplings`计算
2. 函数内部自动根据channel选择对应的K系数（:P用K⁺，:S用K⁻）
3. M₀₈系数使用`4*sqrt(2)/3` ≈ 1.8856，不是4/(3√2) ≈ 0.9428
4. 利用M₀₈ = M₈₀的对称性减少计算量

# 使用示例
```julia
# 预计算K系数
K_coeffs = EffectiveCouplings.calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

# 预计算极化函数
Π_uu = 1e-5 + 0.0im  # fm²（ComplexF64）
Π_ss = 8e-6 + 0.0im  # fm²

# 计算赝标量通道（η/η'）的M矩阵
M_P = calculate_coupling_matrix(Π_uu, Π_ss, K_coeffs, :P)

# 计算标量通道（σ/σ'）的M矩阵
M_S = calculate_coupling_matrix(Π_uu, Π_ss, K_coeffs, :S)
```
"""
@fastmath function calculate_coupling_matrix(Π_uu::ComplexF64, Π_ss::ComplexF64, 
                                             K_coeffs::NamedTuple, channel::Symbol)
    # 根据通道选择对应的K系数
    if channel == :P
        # 赝标量通道（η/η'）使用K⁺系数
        K0 = K_coeffs.K0_plus
        K8 = K_coeffs.K8_plus
        K08 = K_coeffs.K08_plus
        det_K = K_coeffs.det_K_plus
    elseif channel == :S
        # 标量通道（σ/σ'）使用K⁻系数
        K0 = K_coeffs.K0_minus
        K8 = K_coeffs.K8_minus
        K08 = K_coeffs.K08_minus
        det_K = K_coeffs.det_K_minus
    else
        error("Unknown channel: $channel. Use :P (pseudoscalar) or :S (scalar)")
    end

    M00 = K0 - (4.0 / 3.0) * (Π_uu + 2.0 * Π_ss) * det_K
    M08 = K08 + (4.0 / 3.0) * sqrt(2.0) * (Π_uu - Π_ss) * det_K
    M88 = K8 - (4.0 / 3.0) * (2.0 * Π_uu + Π_ss) * det_K

    return ComplexF64[M00 M08; M08 M88]
end

# ----------------------------------------------------------------------------
# 主函数：一般介子传播子
# ----------------------------------------------------------------------------

"""
    meson_propagator_simple(meson_type::Symbol, K_coeffs::NamedTuple, Π::ComplexF64) -> ComplexF64

计算一般介子（π、K、σ_π、σ_K）的传播子。

# 公式
传播子采用：D(p0,p) = 2*K_α / (1 - 4*K_α*Π^±_process(p0,p))

其中：
- α：介子通道（π通道用123，K通道用4567）
- ±：赝标量P用+，标量S用-
- process：散射过程
# 介子类型映射
| 介子 | 标量/赝标量 | 散射过程 | K系数 | 极化函数 |
|------|------|----------|-------|---------|
| π | P | ūu/d̄d | K123_plus | Π_{uu}^P |
| K | P | ūs/d̄s | K4567_plus | Π_{us}^P |
| σ_π | S | ūu/d̄d | K123_minus | Π_{uu}^S |
| σ_K | S | ūs/d̄s | K4567_minus | Π_{us}^S |

# 参数
- `meson_type`: 介子类型（`:pi`, `:K`, `:sigma_pi`, `:sigma_K`）
- `K_coeffs`: 预计算的K系数NamedTuple（通过`EffectiveCouplings.calculate_effective_couplings`获取）
- `Π`: 预计算的极化函数（ComplexF64，单位：fm²）

# 返回值
返回传播子D（ComplexF64，单位：fm²）

# 注意事项
1. K介子应使用Π_{us}而非Π_{uu}
2. 函数自动根据meson_type选择正确的K系数
3. 各向异性修正通过预计算的Π值间接影响传播子
4. K系数应预先通过`EffectiveCouplings.calculate_effective_couplings`计算，批量计算时可复用

# 使用示例
```julia
using .EffectiveCouplings: calculate_effective_couplings

# 预计算K系数（所有介子共用）
K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

# 预计算极化函数
Π_uu_P = 1e-5 + 0.0im  # ūu极化函数（赝标量通道）
Π_us_P = 9e-6 + 0.0im  # ūs极化函数（赝标量通道）
Π_uu_S = 1.2e-5 + 0.0im  # ūu极化函数（标量通道）

# 计算不同介子的传播子
D_pi = meson_propagator_simple(:pi, K_coeffs, Π_uu_P)
D_K = meson_propagator_simple(:K, K_coeffs, Π_us_P)
D_sigma_pi = meson_propagator_simple(:sigma_pi, K_coeffs, Π_uu_S)
```
"""
@inline @fastmath function meson_propagator_simple(meson_type::Symbol, K_coeffs::NamedTuple, 
                                                   Π::ComplexF64)
    # 根据介子类型自动选择K系数
    if meson_type == :pi
        K = K_coeffs.K123_plus  # π是赝标量P通道，用K⁺
    elseif meson_type == :K
        K = K_coeffs.K4567_plus  # K是赝标量P通道，用K⁺
    elseif meson_type == :sigma_pi
        K = K_coeffs.K123_minus  # σ_π是标量S通道，用K⁻
    elseif meson_type == :sigma_K
        K = K_coeffs.K4567_minus  # σ_K是标量S通道，用K⁻
    else
        error("Unknown meson type: $meson_type. Use :pi, :K, :sigma_pi, or :sigma_K")
    end
    
    # D = 2K / (1 - 4KΠ)
    denom = 1.0 - 4.0 * K * Π
    if !_isfinite_complex(denom)
        return 0.0 + 0.0im
    end
    if abs(denom) < PROPAGATOR_DENOM_EPS
        denom += complex(0.0, PROPAGATOR_DENOM_EPS)
    end
    return (2.0 * K) / denom
end

# ----------------------------------------------------------------------------
# 主函数：混合介子传播子
# ----------------------------------------------------------------------------

"""
    meson_propagator_mixed(det_K::Float64, M_matrix::Matrix{ComplexF64}, 
                          q1::Symbol, q2::Symbol, q3::Symbol, q4::Symbol, 
                          channel::Symbol) -> ComplexF64

计算混合介子（η/η'、σ/σ'）的传播子。

# 物理背景
混合介子的传播子涉及SU(3)单态和八重态的混合。
传播子形式为：D = 2det(K)/det(M) × J^T M J'，其中J和J'为流算符向量。

# 散射道映射表
| 散射道 | 物理过程 | ψ | ψ̄ | ψ' | ψ̄' |
|-------|---------|---|---|----|----|
| t道   | q₁+q₂→q₃+q₄ | q1 | q3 | q2 | q4 |
| s道   | q₁+q̄₂→q₃+q̄₄ | q1 | q2 | q3 | q4 |
| u道   | q₁+q₂→q₃+q₄ | q1 | q4 | q2 | q3 |

# 公式
D(q²) = [2det(K) / det(M)] × J^T M J'

其中：
- det(K): 耦合矩阵行列式
- M: 2×2耦合矩阵（包含极化函数修正）
- J, J': 流算符向量（由散射过程的味量子数决定）

# 参数
- `det_K`: 预计算的耦合矩阵行列式（单位：fm⁴，通过`EffectiveCouplings.coupling_matrix_determinant`获取）
- `M_matrix`: 预计算的2×2复数耦合矩阵（通过`calculate_coupling_matrix`获取）
- `q1, q2, q3, q4`: 散射过程q₁+q₂→q₃+q₄的夸克/反夸克类型（`:u`, `:d`, `:s`, `:ubar`, `:dbar`, `:sbar`）
- `channel`: 散射道类型（`:t`, `:s`, `:u`）

# 返回值
返回传播子D（ComplexF64，单位：fm²）

# 散射过程约束
- q1和q3必定是夸克（`:u`, `:d`, `:s`）
- q2和q4要么全为夸克，要么全为反夸克
- t道和u道：纯夸克散射
- s道：夸克-反夸克湮灭

# 注意事项
1. det_K应预先通过`coupling_matrix_determinant(K0, K8, K08)`计算
2. M矩阵应预先通过`calculate_coupling_matrix`计算
3. 函数内部自动根据channel映射场算符
4. 批量计算时可复用det_K和M矩阵

# 使用示例
```julia
using .EffectiveCouplings: calculate_effective_couplings, coupling_matrix_determinant

# 预计算K系数
K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

# 预计算det_K（赝标量通道，用于η/η'）
det_K_P = coupling_matrix_determinant(K_coeffs.K0_plus, K_coeffs.K8_plus, K_coeffs.K08_plus)

# 预计算极化函数和M矩阵
Π_uu_P = 1e-5 + 0.0im
Π_ss_P = 8e-6 + 0.0im
M_P = calculate_coupling_matrix(Π_uu_P, Π_ss_P, K_coeffs, :P)

# 计算η/η'传播子（u+d→u+d散射，t道）
D_eta = meson_propagator_mixed(det_K_P, M_P, :u, :d, :u, :d, :t)

# 计算σ/σ'传播子（u+d̅→u+d̅散射，s道）
det_K_S = coupling_matrix_determinant(K_coeffs.K0_minus, K_coeffs.K8_minus, K_coeffs.K08_minus)
M_S = calculate_coupling_matrix(Π_uu_P, Π_ss_P, K_coeffs, :S)
D_sigma = meson_propagator_mixed(det_K_S, M_S, :u, :dbar, :u, :dbar, :s)
```
"""
@fastmath function meson_propagator_mixed(det_K::Float64, M_matrix::Matrix{ComplexF64},
                                          q1::Symbol, q2::Symbol, q3::Symbol, q4::Symbol,
                                          channel::Symbol)
    # 根据散射道映射表选择场算符
    if channel == :t
        # t道：q₁+q₂→q₃+q₄（直接传递）
        # ψ=q1, ψ̄=q3, ψ'=q2, ψ̄'=q4
        J = calculate_current_vector(q1, q3, channel)
        J_prime = calculate_current_vector(q2, q4, channel)
    elseif channel == :s
        # s道：q₁+q̄₂→q₃+q̄₄（湮灭-产生）
        # ψ=q1, ψ̄=q2, ψ'=q3, ψ̄'=q4
        J = calculate_current_vector(q1, q2, channel)
        J_prime = calculate_current_vector(q3, q4, channel)
    elseif channel == :u
        # u道：q₁+q₂→q₃+q₄（交叉传递）
        # ψ=q1, ψ̄=q4, ψ'=q2, ψ̄'=q3
        J = calculate_current_vector(q1, q4, channel)
        J_prime = calculate_current_vector(q2, q3, channel)
    else
        error("Unknown channel: $channel. Use :t, :s, or :u")
    end
    
    # 计算M矩阵的行列式
    det_M = det(M_matrix)
    if !_isfinite_complex(det_M)
        return 0.0 + 0.0im
    end
    if abs(det_M) < PROPAGATOR_DENOM_EPS
        det_M += complex(0.0, PROPAGATOR_DENOM_EPS)
    end
    
    # 计算传播子 D = 2det(K)/det(M) × J^T M J'
    # 矩阵乘法：J^T (2×1 → 1×2) × M (2×2) × J' (2×1) = 标量
    result = (transpose(J) * M_matrix * J_prime)[1]  # 提取标量值
    if !_isfinite_complex(result)
        return 0.0 + 0.0im
    end
    
    # 乘以前置因子
    D = 2.0 * det_K / det_M * result
    
    return D
end

end # module MesonPropagator
