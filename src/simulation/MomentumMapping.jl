"""
    MomentumMapping

主控模块，实现2→2散射的动量映射计算。

功能：
- 根据入射动量和质心系方向计算出射动量
- 完整的运动学信息存储
- 物理约束验证（四动量守恒、质壳条件）
"""
module MomentumMapping

using LinearAlgebra

# 导入依赖模块
include(joinpath(@__DIR__, "FrameTransformations.jl"))
include(joinpath(@__DIR__, "EllipsoidCalculation.jl"))

using .FrameTransformations
using .EllipsoidCalculation

export ScatteringKinematics, calculate_outgoing_momenta, validate_kinematics

"""
    ScatteringKinematics

散射运动学完整信息结构体。

# 字段
## 入射态 (实验室系)
- `p1_lab, p2_lab::Vector{Float64}`: 入射粒子三动量 [fm⁻¹]
- `E1, E2::Float64`: 入射粒子能量 [fm⁻¹]

## 质心系参数
- `s::Float64`: Mandelstam变量s [fm⁻²]
- `p_star::Float64`: 质心系动量模 [fm⁻¹]
- `beta::Vector{Float64}`: 质心系速度向量
- `gamma::Float64`: 洛伦兹因子

## 出射态 (实验室系)
- `p3_lab, p4_lab::Vector{Float64}`: 出射粒子三动量 [fm⁻¹]
- `E3, E4::Float64`: 出射粒子能量 [fm⁻¹]

## 质心系角度
- `theta_star, phi_star::Float64`: 出射粒子3在质心系的方向角 [rad]

## 椭球参数
- `ellipsoid::EllipsoidParams`: 可行域椭球参数
- `affine_A::Matrix{Float64}`: 仿射变换矩阵
- `affine_b::Vector{Float64}`: 偏移向量
"""
struct ScatteringKinematics
    # 入射态
    p1_lab::Vector{Float64}
    p2_lab::Vector{Float64}
    E1::Float64
    E2::Float64
    
    # 质心系参数
    s::Float64
    p_star::Float64
    beta::Vector{Float64}
    gamma::Float64
    
    # 出射态 (实验室系)
    p3_lab::Vector{Float64}
    p4_lab::Vector{Float64}
    E3::Float64
    E4::Float64
    
    # 质心系角度
    theta_star::Float64
    phi_star::Float64
    
    # 椭球参数
    ellipsoid::EllipsoidParams
    affine_A::Matrix{Float64}
    affine_b::Vector{Float64}
end

"""
    calculate_outgoing_momenta(p1, p2, m1, m2, m3, m4, theta_star, phi_star)

根据入射动量和质心系出射方向，计算完整的散射运动学。

# 参数
- `p1, p2`: 入射粒子三动量向量 [fm⁻¹]
- `m1, m2, m3, m4`: 四个粒子的质量 [fm⁻¹]
- `theta_star`: 质心系极角 [rad], 范围 [0, π]
- `phi_star`: 质心系方位角 [rad], 范围 [0, 2π]

# 返回
- `ScatteringKinematics`: 完整运动学信息

# 计算流程
1. 计算质心系参数（速度、洛伦兹因子、Mandelstam s）
2. 构造仿射变换矩阵A和偏移向量b
3. 根据质心系方向(θ*, φ*)计算质心系出射动量p3*
4. Boost到实验室系得到p3_lab
5. 通过四动量守恒得到p4_lab
6. 计算椭球参数

# 异常
- 如果s低于阈值，抛出error
- 如果Källen函数为负，抛出error

# 示例
```julia
result = calculate_outgoing_momenta(
    [0.5, 0.0, 1.8], [-0.5, 0.0, -1.8],  # p1, p2
    1.52, 1.52, 1.52, 1.52,              # m1, m2, m3, m4
    π/4, π/6                             # theta_star, phi_star
)
```
"""
function calculate_outgoing_momenta(
    p1::Vector{<:Real}, p2::Vector{<:Real},
    m1::Real, m2::Real, m3::Real, m4::Real,
    theta_star::Real, phi_star::Real
)
    # 确保输入是Float64类型
    p1 = Vector{Float64}(p1)
    p2 = Vector{Float64}(p2)
    
    # 1. 计算质心系参数
    cms = calculate_cms_parameters(p1, p2, m1, m2, m3, m4)
    
    # 2. 构造仿射变换
    A, b = build_affine_transform(cms.beta, cms.gamma, cms.E3_star)
    
    # 3. 计算质心系出射动量（粒子3）
    p3_cms = cms.p_star * [
        sin(theta_star) * cos(phi_star), 
        sin(theta_star) * sin(phi_star), 
        cos(theta_star)
    ]
    
    # 4. Boost到实验室系
    p3_lab = boost_to_lab(p3_cms, A, b)
    E3 = boost_energy(cms.E3_star, p3_cms, cms.beta, cms.gamma)
    
    # 5. 四动量守恒得到粒子4
    E1 = sqrt(dot(p1, p1) + m1^2)
    E2 = sqrt(dot(p2, p2) + m2^2)
    p4_lab = p1 + p2 - p3_lab
    E4 = E1 + E2 - E3
    
    # 6. 计算椭球参数
    ellipsoid = calculate_ellipsoid_parameters(A, cms.p_star, center_offset=b)
    
    return ScatteringKinematics(
        p1, p2, E1, E2,
        cms.s, cms.p_star, cms.beta, cms.gamma,
        p3_lab, p4_lab, E3, E4,
        theta_star, phi_star,
        ellipsoid, A, b
    )
end

"""
    validate_kinematics(result, m1, m2, m3, m4; tol=1e-10)

验证散射运动学的物理约束条件。

# 参数
- `result::ScatteringKinematics`: 运动学计算结果
- `m1, m2, m3, m4`: 粒子质量 [fm⁻¹]
- `tol`: 数值容差（默认1e-10）

# 检查项
1. **能量守恒**: |E1 + E2 - E3 - E4| < tol
2. **动量守恒**: |p1 + p2 - p3 - p4| < tol
3. **质壳条件** (粒子3): |E3² - p3² - m3²| < tol
4. **质壳条件** (粒子4): |E4² - p4² - m4²| < tol

# 返回
- `Bool`: 所有检查是否通过
- `Dict`: 各项检查的详细结果

# 异常
- 如果任何检查失败且tol < 1e-8，抛出AssertionError

# 示例
```julia
is_valid, details = validate_kinematics(result, 1.52, 1.52, 1.52, 1.52)
```
"""
function validate_kinematics(result::ScatteringKinematics, 
                             m1::Real, m2::Real, 
                             m3::Real, m4::Real; 
                             tol::Real=1e-10)
    # 能量守恒
    ΔE = abs(result.E1 + result.E2 - result.E3 - result.E4)
    
    # 动量守恒
    Δp = norm(result.p1_lab + result.p2_lab - result.p3_lab - result.p4_lab)
    
    # 质壳条件 (粒子3)
    shell3 = abs(result.E3^2 - dot(result.p3_lab, result.p3_lab) - m3^2)
    
    # 质壳条件 (粒子4)
    shell4 = abs(result.E4^2 - dot(result.p4_lab, result.p4_lab) - m4^2)
    
    # 检查结果
    checks = Dict(
        "energy_conservation" => (ΔE, ΔE < tol),
        "momentum_conservation" => (Δp, Δp < tol),
        "on_shell_p3" => (shell3, shell3 < tol),
        "on_shell_p4" => (shell4, shell4 < tol)
    )
    
    all_passed = all(check[2] for check in values(checks))
    
    # 如果容差很小且检查失败，抛出错误
    if !all_passed && tol < 1e-8
        error_msg = "Physics constraints violated:\n"
        for (name, (value, passed)) in checks
            if !passed
                error_msg *= "  $name: value = $value (tol = $tol)\n"
            end
        end
        error(error_msg)
    end
    
    return all_passed, checks
end

"""
    print_kinematics_info(result; show_details=true)

打印散射运动学信息到控制台。

# 参数
- `result::ScatteringKinematics`: 运动学结果
- `show_details::Bool`: 是否显示详细信息（默认true）

# 输出格式
包括质心系参数、椭球参数、出射动量等信息。

# 示例
```julia
print_kinematics_info(result)
```
"""
function print_kinematics_info(result::ScatteringKinematics; show_details::Bool=true)
    ħc = 197.327  # MeV·fm
    
    println("=" ^ 60)
    println("散射运动学信息")
    println("=" ^ 60)
    
    println("\n质心系参数:")
    println("  s = $(round(result.s, digits=6)) fm⁻² ($(round(result.s * ħc^2, digits=2)) MeV²)")
    println("  √s = $(round(sqrt(result.s), digits=6)) fm⁻¹ ($(round(sqrt(result.s) * ħc, digits=2)) MeV)")
    println("  p* = $(round(result.p_star, digits=6)) fm⁻¹")
    println("  β = $(round(norm(result.beta), digits=6))")
    println("  γ = $(round(result.gamma, digits=6))")
    
    println("\n椭球参数:")
    println("  中心: $(round.(result.ellipsoid.center, digits=4)) fm⁻¹")
    println("  半轴长: $(round.(result.ellipsoid.half_lengths, digits=4)) fm⁻¹")
    
    if show_details
        println("\n出射动量 (实验室系):")
        println("  θ* = $(round(rad2deg(result.theta_star), digits=2))°")
        println("  φ* = $(round(rad2deg(result.phi_star), digits=2))°")
        println("  p₃ = $(round.(result.p3_lab, digits=4)) fm⁻¹")
        println("  E₃ = $(round(result.E3, digits=4)) fm⁻¹ ($(round(result.E3 * ħc, digits=2)) MeV)")
        println("  p₄ = $(round.(result.p4_lab, digits=4)) fm⁻¹")
        println("  E₄ = $(round(result.E4, digits=4)) fm⁻¹ ($(round(result.E4 * ħc, digits=2)) MeV)")
    end
    
    println("=" ^ 60)
end

"""
    calculate_mandelstam_t(result)

从运动学结果计算Mandelstam变量t。

# 参数
- `result::ScatteringKinematics`: 运动学结果

# 返回
- `Float64`: Mandelstam t = (p1 - p3)² [fm⁻²]

# 示例
```julia
t = calculate_mandelstam_t(result)
```
"""
function calculate_mandelstam_t(result::ScatteringKinematics)
    Δp = result.p1_lab - result.p3_lab
    ΔE = result.E1 - result.E3
    return ΔE^2 - dot(Δp, Δp)
end

"""
    calculate_mandelstam_u(result, m1, m2, m3, m4)

从运动学结果计算Mandelstam变量u。

# 参数
- `result::ScatteringKinematics`: 运动学结果
- `m1, m2, m3, m4`: 粒子质量 [fm⁻¹]

# 返回
- `Float64`: Mandelstam u = (p1 - p4)² [fm⁻²]

或通过约束关系: u = m1² + m2² + m3² + m4² - s - t

# 示例
```julia
u = calculate_mandelstam_u(result, 1.52, 1.52, 1.52, 1.52)
```
"""
function calculate_mandelstam_u(result::ScatteringKinematics,
                               m1::Real, m2::Real, m3::Real, m4::Real)
    # 方法1：直接计算
    Δp = result.p1_lab - result.p4_lab
    ΔE = result.E1 - result.E4
    u_direct = ΔE^2 - dot(Δp, Δp)
    
    # 方法2：通过约束关系（用于验证）
    t = calculate_mandelstam_t(result)
    u_constraint = m1^2 + m2^2 + m3^2 + m4^2 - result.s - t
    
    # 检查一致性
    if abs(u_direct - u_constraint) > 1e-8
        @warn "Mandelstam u calculation mismatch" u_direct u_constraint
    end
    
    return u_direct
end

end # module
