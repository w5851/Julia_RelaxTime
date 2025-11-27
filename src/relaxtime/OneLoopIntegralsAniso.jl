"""
单圈积分在各向异性一阶修正下的修正项B0_correction和A_correction的计算函数
"""
module OneLoopIntegralsCorrection

export B0_correction, A_correction,
    A_aniso

include("../QuarkDistribution_Aniso.jl")
include("OneLoopIntegrals.jl")
using QuadGK: quadgk
using .PNJLQuarkDistributions_Aniso: correction_cos_theta_coefficient, distribution_aniso
using .OneLoopIntegrals: internal_momentum, EPS_K, DEFAULT_RTOL, DEFAULT_ATOL,
    energy_cutoff, singularity_k_positive, const_integral_term_A

""" k>0时
对x积分所得结果的Sokhotski–Plemelj formula实部部分(柯西主值积分)
被积函数展开到一阶,与x相关的项只多了x^2
对x的积分的被积函数可以抽象为x^2/(Ax+B)
积分变量x=cosθ从-1积到1
A对应的形参名:coeff_x
B对应的形参名:denominator_const
"""
function real_integral_tool(coeff_x::Float64, denominator_const::Float64)
    term1 = -2*denominator_const/coeff_x^2
    term2 = denominator_const^2/coeff_x^3
    log_arg = (coeff_x + denominator_const) / (coeff_x - denominator_const)
    term3 = log(abs(log_arg))
    result = term1 + term2 * term3
    return result
end

"""对x积分的虚部部分"""
function imag_integral_tool(coeff_x::Float64, denominator_const::Float64)
    return π *denominator_const ^ 2 / coeff_x^3
end

"""阶跃函数,确保奇点出现在积分范围内时才贡献虚部"""
@inline function heaviside_step(coeff_x::Float64, denominator_const::Float64, x_min::Float64, x_max::Float64)
    # 计算奇点位置
    x = -denominator_const / coeff_x 
    return (x >= x_min && x <= x_max) ? 1.0 : 0.0
end

"""coeff_x和denominator_const的计算函数"""
@inline function compute_coefficients(λ::Float64, k::Float64, m::Float64, m_prime::Float64, E::Float64)
    if k>EPS_K # k>0
        denominator_const = λ^2+2*λ*E + m^2 - m_prime^2-k^2
        p = internal_momentum(E,m)
        coeff_x = 2*p*k
        return coeff_x, denominator_const
    else # k=0
        denominator_const = λ^2+m^2-m_prime^2
        coeff_E = 2*λ
        return coeff_E, denominator_const
    end
end

"""k=0时的积分实部被积函数"""
function real_integrand_k_zero(sign_::Symbol, λ::Float64, m::Float64, m_prime::Float64, E::Float64,
    ξ::Float64, T::Float64, μ::Float64, Φ::Float64, Φbar::Float64)
    coeff_E, denominator_const = compute_coefficients(λ, 0.0, m, m_prime, E)
    p = internal_momentum(E,m)
    # 含各向异性修正项
    real_ = (4.0/3.0)*p/(coeff_E*E+denominator_const)*correction_cos_theta_coefficient(sign_, p, m, μ, T, Φ, Φbar, ξ)
    return real_
end

"""k>0时的积分实部被积函数"""
function real_integrand_k_positive(sign_::Symbol, λ::Float64, k::Float64, m::Float64, m_prime::Float64, E::Float64,
    ξ::Float64, T::Float64, μ::Float64, Φ::Float64, Φbar::Float64)
    coeff_x, denominator_const = compute_coefficients(λ, k, m, m_prime, E)
    p = internal_momentum(E,m)
    # 含各向异性修正项
    real_ = 2.0*p*real_integral_tool(coeff_x, denominator_const)*
        correction_cos_theta_coefficient(sign_, p, m, μ, T, Φ, Φbar, ξ)
    return real_
end

"""k=0时的积分虚部函数"""
function imag_integrand_k_zero(sign_::Symbol, λ::Float64, m::Float64, m_prime::Float64,
    ξ::Float64, T::Float64, μ::Float64, Φ::Float64, Φbar::Float64)
    coeff_E, denominator_const = compute_coefficients(λ, 0.0, m, m_prime, 0.0)
    E_pole = -denominator_const / coeff_E
    Θ = heaviside_step(coeff_E, denominator_const, m, energy_cutoff(m))
    if Θ == 0.0
        return 0.0
    else
        p_pole = internal_momentum(E_pole,m)
        # 含各向异性修正项
        imag_ = π*(4.0/3.0)*p_pole/coeff_E*
            correction_cos_theta_coefficient(sign_, p_pole, m, μ, T, Φ, Φbar, ξ)
        return imag_
    end
end

"""k>0时的积分虚部被积函数"""
function imag_integrand_k_positive(sign_::Symbol, λ::Float64, k::Float64, m::Float64, m_prime::Float64, E::Float64,
    ξ::Float64, T::Float64, μ::Float64, Φ::Float64, Φbar::Float64)

    coeff_x, denominator_const = compute_coefficients(λ, k, m, m_prime, E)
    
    p = internal_momentum(E,m)
    # 含各向异性修正项
    imag_ = 2.0*p*imag_integral_tool(coeff_x, denominator_const)*
            correction_cos_theta_coefficient(sign_, p, m, μ, T, Φ, Φbar, ξ)
    return imag_
end

"""k=0时的 B0分量 含各向异性修正项的积分计算"""
function tilde_B0_correction_k_zero(sign_::Symbol, λ::Float64, m::Float64, m_prime::Float64, μ::Float64, T::Float64,
    Φ::Float64, Φbar::Float64, ξ::Float64; rtol::Float64=DEFAULT_RTOL, atol::Float64=DEFAULT_ATOL)
    Emin = m
    Emax = energy_cutoff(m)

    integrand_real(E) = real_integrand_k_zero(sign_, λ, m, m_prime, E,
        ξ, T, μ, Φ, Φbar) # 闭包被积函数-实部
    real_part, _ = quadgk(integrand_real, Emin, Emax; rtol=rtol, atol=atol)
    # 计算虚部
    imag_part = imag_integrand_k_zero(sign_, λ, m, m_prime,
        ξ, T, μ, Φ, Φbar)

    return real_part , imag_part
end

"""k>0时的 B0分量 含各向异性修正项的积分计算"""
function tilde_B0_correction_k_positive(sign_::Symbol, λ::Float64, k::Float64, m::Float64, m_prime::Float64, μ::Float64, T::Float64,
    Φ::Float64, Φbar::Float64, ξ::Float64; rtol::Float64=DEFAULT_RTOL, atol::Float64=DEFAULT_ATOL)
    Emin = m
    Emax = energy_cutoff(m)

    intervals, sign_type = singularity_k_positive(λ, k, m, m_prime, Emin, Emax)

    integrand_real(E) = real_integrand_k_positive(sign_, λ, k, m, m_prime, E,
        ξ, T, μ, Φ, Φbar) # 闭包被积函数-实部
    integrand_imag(E) = imag_integrand_k_positive(sign_, λ, k, m, m_prime, E,
        ξ, T, μ, Φ, Φbar) # 闭包被积函数-虚部

    real_part, _ = quadgk(integrand_real, Emin, Emax; rtol=rtol, atol=atol)
    
    # 根据区间类型计算虚部
    imag_part = 0.0
    if !isempty(intervals)
        for (E1, E2) in intervals
            imag_part_segment, _ = quadgk(integrand_imag, E1, E2; rtol=rtol, atol=atol)
            imag_part += imag_part_segment
        end
        imag_part *= sign(λ)
    end

    return real_part , imag_part
end

"""含各向异性修正项的 B0分量 积分计算"""
function tilde_B0_correction(sign_::Symbol, λ::Float64, k::Float64, m::Float64, m_prime::Float64, μ::Float64, T::Float64,
    Φ::Float64, Φbar::Float64, ξ::Float64; rtol::Float64=DEFAULT_RTOL, atol::Float64=DEFAULT_ATOL)
    if k > EPS_K
        return tilde_B0_correction_k_positive(sign_, λ, k, m, m_prime, μ, T,
            Φ, Φbar, ξ; rtol=rtol, atol=atol)
    else
        return tilde_B0_correction_k_zero(sign_, λ, m, m_prime, μ, T,
            Φ, Φbar, ξ; rtol=rtol, atol=atol)
    end
end

"""
    B0_correction(λ, k, m1, m2, μ1, μ2, T, Φ, Φbar, ξ; rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL)
动量各向异性下B0的一阶修正
"""
function B0_correction(λ::Float64, k::Float64, m1::Float64, m2::Float64, μ1::Float64, μ2::Float64, T::Float64,
    Φ::Float64, Φbar::Float64, ξ::Float64; rtol::Float64=DEFAULT_RTOL, atol::Float64=DEFAULT_ATOL)

    tol_kwargs = (; rtol=rtol, atol=atol)

    term1 = tilde_B0_correction(:quark, -λ, k, m1, m2, μ1, T, Φ, Φbar, ξ; tol_kwargs...)
    term2 = tilde_B0_correction(:antiquark, λ, k, m1, m2, μ1, T, Φ, Φbar, ξ; tol_kwargs...)
    term3 = tilde_B0_correction(:quark, λ, k, m2, m1, μ2, T, Φ, Φbar, ξ; tol_kwargs...)
    term4 = tilde_B0_correction(:antiquark, -λ, k, m2, m1, μ2, T, Φ, Φbar, ξ; tol_kwargs...)

    real_part = term1[1] - term2[1] + term3[1] - term4[1]
    imag_part = term1[2] - term2[2] + term3[2] - term4[2]
    return real_part, imag_part
end

# -------------------------------------
"""
    A_correction(m, μ, T, Φ, Φbar, ξ, nodes_p, weights_p)
计算单线积分函数A对ξ的一阶修正项(需要加上各向同性项才能进行后续计算，即零阶项),需要传入预生成的动量的积分节点与权重
可以使用 GaussLegendre 模块中的 DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS (积分上限 10.0 fm⁻¹)
"""
function A_correction(m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64,
    ξ::Float64, nodes_p::Vector{Float64}, weights_p::Vector{Float64})
    integral = 0.0
    @inbounds for i in eachindex(nodes_p)
        node_p = nodes_p[i]
        weight_p = weights_p[i]
        E = sqrt(node_p^2 + m^2)
        quark_correction = correction_cos_theta_coefficient(:quark, node_p, m, μ, T, Φ, Φbar, ξ)
        antiquark_correction = correction_cos_theta_coefficient(:antiquark, node_p, m, μ, T, Φ, Φbar, ξ)
        integral += weight_p * node_p^2 / E * (quark_correction + antiquark_correction)
    end
    return 4.0 * integral / 3.0
end

"""
    A_aniso(m, μ, T, Φ, Φbar, ξ, nodes_p, weights_p, nodes_cosθ, weights_cosθ)
计算单线积分函数A在动量各向异性下的完整形式,需要传入预生成的动量和角度的积分节点与权重
可以使用 GaussLegendre 模块中的常量：
- DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS (p ∈ [0, 10] fm⁻¹)
- DEFAULT_COSΘ_NODES, DEFAULT_COSΘ_WEIGHTS (cosθ ∈ [-1, 1])
"""
function A_aniso(m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64,
    ξ::Float64, nodes_p::Vector{Float64}, weights_p::Vector{Float64},
    nodes_cosθ::Vector{Float64}, weights_cosθ::Vector{Float64})
    # 只计算分布函数的积分部分，常数项单独处理
    integral = 0.0
    @inbounds @simd for i in eachindex(nodes_p)
        node_p = nodes_p[i]
        weight_p = weights_p[i]
        # 对角度进行积分
        angle_integral_quark = 0.0
        angle_integral_antiquark = 0.0
        E = sqrt(node_p^2 + m^2)  # 使用各向同性能量作为分母
        @inbounds @simd for j in eachindex(nodes_cosθ)
            cosθ = nodes_cosθ[j]
            weight_cosθ = weights_cosθ[j]
            dist_quark = distribution_aniso(:quark, node_p, m, μ, T, Φ, Φbar, ξ, cosθ)
            dist_antiquark = distribution_aniso(:antiquark, node_p, m, μ, T, Φ, Φbar, ξ, cosθ)
            angle_integral_quark += weight_cosθ * node_p^2 / E * dist_quark
            angle_integral_antiquark += weight_cosθ * node_p^2 / E * dist_antiquark
        end
        integral += weight_p * (angle_integral_quark + angle_integral_antiquark)
    end
    # 分布函数项：angle integration给出因子2（cosθ权重和），φ积分贡献2π→归一化为2，总共4
    # 常数项：无角度依赖，需要完整4π因子
    return 2.0 * integral - 4.0 * const_integral_term_A(m)
end
end # module OneLoopIntegralsCorrection