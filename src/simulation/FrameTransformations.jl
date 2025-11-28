"""
    FrameTransformations

洛伦兹变换工具模块，用于质心系(CMS)和实验室系(Lab)之间的坐标变换。

提供功能：
- Källen函数计算
- 质心系参数计算（速度β、洛伦兹因子γ、Mandelstam变量s）
- 仿射变换矩阵构造
- 四动量boost变换
"""
module FrameTransformations

using LinearAlgebra

export kallen_lambda, calculate_cms_parameters, build_affine_transform, boost_to_lab, boost_energy

"""
    kallen_lambda(a, b, c)

计算Källen函数: λ(a,b,c) = a² + b² + c² - 2ab - 2ac - 2bc

# 参数
- `a, b, c`: 输入参数（通常为Mandelstam变量或质量平方）

# 返回
- Källen函数的值

# 示例
```julia
λ = kallen_lambda(s, m3^2, m4^2)
```
"""
function kallen_lambda(a::Real, b::Real, c::Real)
    return a^2 + b^2 + c^2 - 2*a*b - 2*a*c - 2*b*c
end

"""
    calculate_cms_parameters(p1, p2, m1, m2, m3, m4)

计算质心系参数，包括速度、洛伦兹因子、Mandelstam变量s和质心系动量。

# 参数
- `p1, p2`: 入射粒子三动量向量 [fm⁻¹]
- `m1, m2, m3, m4`: 四个粒子的质量 [fm⁻¹]

# 返回
返回NamedTuple包含：
- `s`: Mandelstam变量s [fm⁻²]
- `p_star`: 质心系动量模 [fm⁻¹]
- `E3_star, E4_star`: 质心系出射粒子能量 [fm⁻¹]
- `beta`: 质心系速度向量（3-vector，无量纲）
- `gamma`: 洛伦兹因子（无量纲）
- `beta2`: β²（无量纲）

# 异常
- 如果 s ≤ (m3 + m4)²，抛出AssertionError（低于阈值）

# 示例
```julia
cms = calculate_cms_parameters([0.5, 0, 1.8], [-0.5, 0, -1.8], 1.52, 1.52, 1.52, 1.52)
```
"""
function calculate_cms_parameters(p1::Vector{<:Real}, p2::Vector{<:Real}, 
                                  m1::Real, m2::Real, 
                                  m3::Real, m4::Real)
    # 计算入射粒子能量
    E1 = sqrt(dot(p1, p1) + m1^2)
    E2 = sqrt(dot(p2, p2) + m2^2)
    
    # 总三动量和总能量
    Ptot = p1 .+ p2
    Etot = E1 + E2
    
    # Mandelstam s
    s = Etot^2 - dot(Ptot, Ptot)
    
    # 阈值检查
    threshold = (m3 + m4)^2
    if s <= threshold + 1e-12  # 添加小量避免数值误差
        error("Below threshold: s = $s fm⁻², threshold = $threshold fm⁻²")
    end
    
    # 质心系速度（向量）
    β = Ptot / Etot
    β2 = dot(β, β)
    
    # 洛伦兹因子
    γ = 1 / sqrt(1 - β2)
    
    # 质心系动量模和能量
    sqrt_s = sqrt(s)
    lambda_val = kallen_lambda(s, m3^2, m4^2)
    
    if lambda_val < 0
        error("Negative Källen function: λ = $lambda_val (s=$s, m3²=$(m3^2), m4²=$(m4^2))")
    end
    
    p_star = sqrt(lambda_val) / (2 * sqrt_s)
    E3_star = sqrt(p_star^2 + m3^2)
    E4_star = sqrt_s - E3_star
    
    return (s=s, p_star=p_star, E3_star=E3_star, E4_star=E4_star,
            beta=β, gamma=γ, beta2=β2)
end

"""
    build_affine_transform(β, γ, E_star)

构造从质心系到实验室系的仿射变换矩阵和偏移向量。

变换公式：p_lab = A * p_cms + b

其中：
- A = I + (γ-1)/β² * β⊗β  （当β² < 1e-14时，A = I）
- b = γ * E_star * β

# 参数
- `β`: 质心系速度向量（3-vector）
- `γ`: 洛伦兹因子
- `E_star`: 质心系粒子能量 [fm⁻¹]

# 返回
- `A`: 3×3仿射变换矩阵
- `b`: 3维偏移向量 [fm⁻¹]

# 示例
```julia
A, b = build_affine_transform(cms.beta, cms.gamma, cms.E3_star)
```
"""
function build_affine_transform(β::Vector{<:Real}, γ::Real, E_star::Real)
    β2 = dot(β, β)
    
    # 小速度情况：接近静止系
    if β2 < 1e-14
        A = Matrix{Float64}(I, 3, 3)
    else
        # A = I + (γ-1)/β² * β⊗β
        A = Matrix{Float64}(I, 3, 3) + ((γ - 1) / β2) * (β * β')
    end
    
    # 偏移向量
    b = γ * E_star * β
    
    return A, b
end

"""
    boost_to_lab(p_star, A, b)

将质心系三动量boost到实验室系。

# 参数
- `p_star`: 质心系三动量向量 [fm⁻¹]
- `A`: 仿射变换矩阵（3×3）
- `b`: 偏移向量（3-vector）[fm⁻¹]

# 返回
- 实验室系三动量向量 [fm⁻¹]

# 示例
```julia
p3_cms = cms.p_star * [sin(θ)*cos(φ), sin(θ)*sin(φ), cos(θ)]
p3_lab = boost_to_lab(p3_cms, A, b)
```
"""
function boost_to_lab(p_star::Vector{<:Real}, A::Matrix{<:Real}, b::Vector{<:Real})
    return A * p_star + b
end

"""
    boost_energy(E_star, p_star, β, γ)

计算boost后的能量（从质心系到实验室系）。

# 参数
- `E_star`: 质心系能量 [fm⁻¹]
- `p_star`: 质心系三动量向量 [fm⁻¹]
- `β`: 质心系速度向量
- `γ`: 洛伦兹因子

# 返回
- 实验室系能量 [fm⁻¹]

# 示例
```julia
E3_lab = boost_energy(cms.E3_star, p3_cms, cms.beta, cms.gamma)
```
"""
function boost_energy(E_star::Real, p_star::Vector{<:Real}, 
                      β::Vector{<:Real}, γ::Real)
    return γ * (E_star + dot(p_star, β))
end

end # module
