"""
# EffectiveCouplings.jl

有效耦合常数计算模块，基于单圈积分A函数计算介子传播子所需的有效耦合系数。

公式参考 `doc/formula/K_有效耦合常数byA.md`。
"""
module EffectiveCouplings

include("../Constants_PNJL.jl")
include("OneLoopIntegrals.jl")
include("OneLoopIntegralsAniso.jl")
using .Constants_PNJL: N_color, ħc_MeV_fm
using .OneLoopIntegrals: A
using .OneLoopIntegralsCorrection: A_aniso

export calculate_effective_couplings, coupling_matrix_determinant
export calculate_G_from_A

# ----------------------------------------------------------------------------
# 转换函数：从A函数计算G^f
# ----------------------------------------------------------------------------

"""
    calculate_G_from_A(A_f::Float64, m_f::Float64; Nc::Int=3) -> Float64

从单圈积分A_f计算夸克凝聚相关函数G^f。

# 公式
```math
G^f = -\\frac{N_c}{4\\pi^2} m_f A_f(T, \\mu)
```

# 参数
- `A_f`: 单圈积分A函数的值（单位：fm）
- `m_f`: 该味夸克的质量（单位：fm⁻¹）
- `Nc`: 色数（默认3）

# 返回值
- `G_f`: 夸克凝聚相关函数（单位：无量纲）

# 示例
```julia
using .OneLoopIntegrals: A
T_inv_fm = 150.0 / 197.327  # 150 MeV
μ_inv_fm = 0.0
m_inv_fm = 300.0 / 197.327  # 300 MeV
Φ = 0.5
Φbar = 0.5

A_u = A(T_inv_fm, μ_inv_fm, m_inv_fm, Φ, Φbar)
G_u = calculate_G_from_A(A_u)
```
"""
@inline @fastmath function calculate_G_from_A(A_f::Float64, m_f::Float64; Nc::Int=3)
    # 与 C++ set_pse() 中 G[i] 的构造一致（去掉了外层的 K_f 因子，
    # 该因子会在 calculate_effective_couplings 的 K * G_f 中引入）
    return -Nc / (4.0 * π^2) * (m_f * A_f)
end

@inline @fastmath function calculate_G_from_A(A_f::Float64; Nc::Int=3)
    throw(ArgumentError("calculate_G_from_A(A_f) is deprecated; call calculate_G_from_A(A_f, m_f)"))
end

# ----------------------------------------------------------------------------
# 辅助函数：计算耦合矩阵行列式
# ----------------------------------------------------------------------------

"""
    coupling_matrix_determinant(K0::Float64, K8::Float64, K08::Float64) -> Float64

计算混合介子传播子所需的耦合矩阵行列式（也称为det_K）。

# 物理背景
在η-η'介子系统（赝标量通道）或σ-σ'介子系统（标量通道）中，存在SU(3)单态和八重态的混合。
混合介子的传播子需要对2×2矩阵求逆，其行列式det K直接影响传播子的极点位置和介子质量。

# 公式
det(K^±) = K_0^± K_8^± - (K_{08}^±)^2

其中±表示通道类型：
- +：赝标量通道（P通道），用于η/η'混合介子
- -：标量通道（S通道），用于σ/σ'混合介子

# 参数
参数名使用通用命名，适用于两种通道：
- `K0`: 单态通道耦合系数（K0_plus用于P通道，K0_minus用于S通道）（单位：fm²）
- `K8`: 八重态通道耦合系数（K8_plus用于P通道，K8_minus用于S通道）（单位：fm²）
- `K08`: 混合通道耦合系数（K08_plus用于P通道，K08_minus用于S通道）（单位：fm²）

# 返回值
- `det_K`: 耦合矩阵行列式（单位：fm⁴）

# 物理约束
对于物理上有意义的传播子，要求det K > 0，否则传播子会出现非物理的极点。

# 使用示例
```julia
# 计算赝标量通道（η/η'）的行列式
K_coeffs = calculate_effective_couplings(G, K, G_u, G_s)
det_K_P = coupling_matrix_determinant(K_coeffs.K0_plus, K_coeffs.K8_plus, K_coeffs.K08_plus)

# 计算标量通道（σ/σ'）的行列式
det_K_S = coupling_matrix_determinant(K_coeffs.K0_minus, K_coeffs.K8_minus, K_coeffs.K08_minus)
```

# 注意事项
1. 此函数适用于赝标量和标量两种通道，计算公式相同，只需传入相应的K系数即可
2. det K的符号直接关系到介子传播子的因果性
3. 在某些极端参数下（如非常高的温度或密度），det K可能变负，表示模型失效
"""
@inline @fastmath function coupling_matrix_determinant(K0::Float64, K8::Float64, 
                                                       K08::Float64)
    return K0 * K8 - K08^2
end

# ----------------------------------------------------------------------------
# 主函数：计算有效耦合系数
# ----------------------------------------------------------------------------

"""
    calculate_effective_couplings(G::Float64, K::Float64, G_u::Float64, G_s::Float64) -> NamedTuple

计算3味PNJL模型中的10个有效耦合系数。

# 物理背景
在随机相位近似(RPA)下，介子传播子由有效耦合系数K_α^±描述，这些系数依赖于原始耦合常数G、K
和夸克凝聚相关函数G^μ（u/d味）、G^s（s味）。

# 公式
K_0^± = G ∓ (1/3)K(2G^μ + G^s)
K_{123}^± = G ± (1/2)KG^s
K_{4567}^± = G ± (1/2)KG^μ
K_8^± = G ± (1/6)K(4G^μ - G^s)
K_{08}^± = ±(1/6)√2 K(G^μ - G^s)

其中：
- K_0: 单态通道（σ/f_0）
- K_{123}: π介子通道（π⁰, π±）
- K_{4567}: K介子通道（K⁰, K±）
- K_8: 八重态通道（η_8）
- K_{08}: 混合通道（η-η'混合）

# 参数
- `G`: 四夸克相互作用耦合常数（单位：fm²）
- `K`: 't Hooft六夸克相互作用耦合常数（单位：fm⁵）
- `G_u`: u/d夸克凝聚相关函数（无量纲），通过`calculate_G_from_A(A_u)`计算
- `G_s`: s夸克凝聚相关函数（无量纲），通过`calculate_G_from_A(A_s)`计算

# 返回值
返回一个NamedTuple，包含10个有效耦合系数和2个行列式（单位：fm²和fm⁴）：
(K0_plus, K0_minus, K123_plus, K123_minus, K4567_plus, K4567_minus,
 K8_plus, K8_minus, K08_plus, K08_minus, det_K_plus, det_K_minus)

其中：
- `*_plus`: 赝标量通道（P通道）
- `*_minus`: 标量通道（S通道）

# 物理意义
- 赝标量通道（+）：对应自旋0、负宇称的介子（π、K、η、η'等）
- 标量通道（-）：对应自旋0、正宇称的介子（σ_π、σ_K、σ, σ'等）

# 使用示例
```julia
using .Constants_PNJL: G_fm2, K_fm5, ħc_MeV_fm
using .OneLoopIntegrals: A

# 设置物理参数
T_MeV = 150.0
μ_MeV = 0.0
m_u_MeV = 300.0
m_s_MeV = 500.0

T_inv_fm = T_MeV / ħc_MeV_fm
μ_inv_fm = μ_MeV / ħc_MeV_fm
m_u_inv_fm = m_u_MeV / ħc_MeV_fm
m_s_inv_fm = m_s_MeV / ħc_MeV_fm

Φ = 0.5
Φbar = 0.5

# 计算A函数
A_u = A(T_inv_fm, μ_inv_fm, m_u_inv_fm, Φ, Φbar)
A_s = A(T_inv_fm, μ_inv_fm, m_s_inv_fm, Φ, Φbar)

# 计算G^f
G_u = calculate_G_from_A(A_u)
G_s = calculate_G_from_A(A_s)

# 计算有效耦合系数
K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

println("K₀⁺ = ", K_coeffs.K0_plus, " fm²")
println("K₁₂₃⁻ = ", K_coeffs.K123_minus, " fm² (π介子)")
println("det(K⁺) = ", K_coeffs.det_K_plus, " fm⁴ (η/η'混合)")
println("det(K⁻) = ", K_coeffs.det_K_minus, " fm⁴ (σ/σ'混合)")
```

# 注意事项
1. 输入的G和K应使用fm为单位的自然单位制
2. G_u和G_s是无量纲量，通过A函数计算得到
3. 在手征极限下（G^μ = G^s = 0），所有K_α^±均退化为G
4. K_{08}^±的符号差异来自于η-η'混合的非对角项
"""
@inline @fastmath function calculate_effective_couplings(G::Float64, K::Float64, 
                                                         G_u::Float64, G_s::Float64)
    # K₀± = G ∓ (1/3)K(2G^μ + G^s)
    term_0 = (1.0 / 3.0) * K * (2.0 * G_u + G_s)
    K0_plus = G - term_0
    K0_minus = G + term_0
    
    # K₁₂₃± = G ± (1/2)KG^s（π介子通道）
    term_123 = 0.5 * K * G_s
    K123_plus = G + term_123
    K123_minus = G - term_123
    
    # K₄₅₆₇± = G ± (1/2)KG^μ（K介子通道）
    term_4567 = 0.5 * K * G_u
    K4567_plus = G + term_4567
    K4567_minus = G - term_4567
    
    # K₈± = G ± (1/6)K(4G^μ - G^s)
    term_8 = (1.0 / 6.0) * K * (4.0 * G_u - G_s)
    K8_plus = G + term_8
    K8_minus = G - term_8
    
    # K₀₈± = ±(1/6)√2 K(G^μ - G^s)（混合通道）
    term_08 = (1.0 / 6.0) * sqrt(2.0) * K * (G_u - G_s)
    K08_plus = term_08
    K08_minus = -term_08
    
    # 计算耦合矩阵行列式（用于混合介子传播子）
    det_K_plus = coupling_matrix_determinant(K0_plus, K8_plus, K08_plus)
    det_K_minus = coupling_matrix_determinant(K0_minus, K8_minus, K08_minus)
    
    return (K0_plus=K0_plus, K0_minus=K0_minus,
            K123_plus=K123_plus, K123_minus=K123_minus,
            K4567_plus=K4567_plus, K4567_minus=K4567_minus,
            K8_plus=K8_plus, K8_minus=K8_minus,
            K08_plus=K08_plus, K08_minus=K08_minus,
            det_K_plus=det_K_plus, det_K_minus=det_K_minus)
end

end # module EffectiveCouplings
