"""
各向异性下对夸克有效分布函数对ξ的一阶展开修正项:
在ξ很小的情况下, 使用线性近似.

具体公式见文档doc/formula/PNJL_夸克有效分布函数_动量各向异性.md
"""
module PNJLQuarkDistributions_Aniso

export distribution_aniso_correction, distribution_aniso, correction_cos_theta_coefficient

include("QuarkDistribution.jl")
using .PNJLQuarkDistributions: quark_distribution, antiquark_distribution
"""计算PNJL模型中夸克有效分布函数对能量的导数"""
function quark_df_dE(E_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64)
    # 计算温度的倒数
    β_fm = 1 / T_inv_fm

    # 计算指数项
    exp_term = clamp(exp(-(E_inv_fm - μ_inv_fm) * β_fm), 1e-200, 1e200)
    exp_term2 = exp_term * exp_term
    exp_term3 = exp_term2 * exp_term

    # 计算分子和分母
    numerator = Φ * exp_term + 2 * Φbar * exp_term2 + exp_term3
    denominator = 1 + 3 * Φ * exp_term + 3 * Φbar * exp_term2 + exp_term3
    d_numerator = -(Φ * exp_term + 4 * Φbar * exp_term2 + 3 * exp_term3)
    d_denominator = -(3 * Φ * exp_term + 6 * Φbar * exp_term2 + 3 * exp_term3)

    # x=β(E - μ)
    df_dx = d_numerator * denominator - numerator * d_denominator
    df_dx /= denominator ^ 2

    # 计算导数
    df_dE = -β_fm * df_dx

    return df_dE
end

"""计算PNJL模型中反夸克有效分布函数对能量的导数"""
function antiquark_df_dE(E_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64)
    # 计算温度的倒数
    β_fm = 1 / T_inv_fm

    # 计算指数项
    exp_term = clamp(exp(-(E_inv_fm + μ_inv_fm) * β_fm), 1e-200, 1e200)
    exp_term2 = exp_term * exp_term
    exp_term3 = exp_term2 * exp_term

    # 计算分子和分母
    numerator = Φbar * exp_term + 2 * Φ * exp_term2 + exp_term3
    denominator = 1 + 3 * Φbar * exp_term + 3 * Φ * exp_term2 + exp_term3
    d_numerator = -(Φbar * exp_term + 4 * Φ * exp_term2 + 3 * exp_term3)
    d_denominator = -(3 * Φbar * exp_term + 6 * Φ * exp_term2 + 3 * exp_term3)

    # x=β(E + μ)
    df_dx = d_numerator * denominator - numerator * d_denominator
    df_dx /= denominator ^ 2

    # 计算导数
    df_dE = -β_fm * df_dx

    return df_dE    
end

# -------------------------------------
"""
    correction_cos_theta_coefficient(sign_, p_inv_fm, m_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, ξ)
计算PNJL模型中分布函数一阶修正项中cosθ的系数
"""
function correction_cos_theta_coefficient(sign_::Symbol, p_inv_fm::Float64, m_inv_fm::Float64, 
    μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64, ξ::Float64)
    E_inv_fm = sqrt(p_inv_fm^2 + m_inv_fm^2)
    coeff = 0.5 * ξ * (p_inv_fm^2) / E_inv_fm
    if sign_ === :quark
        df_dE = quark_df_dE(E_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar) 
    elseif sign_ === :antiquark
        df_dE = antiquark_df_dE(E_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar)
    else
        throw(ArgumentError("sign_ must be :quark or :antiquark"))
    end

    return coeff*df_dE
end

# ------------------------------------
"""计算PNJL模型中夸克有效分布函数的各向异性一阶修正项"""
function quark_distribution_aniso_correction(p_inv_fm::Float64, m_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, cosθ::Float64)
    E_inv_fm = sqrt(p_inv_fm^2 + m_inv_fm^2)
    df_dE = quark_df_dE(E_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar)

    return 0.5*ξ*(p_inv_fm*cosθ)^2 / E_inv_fm * df_dE
end

"""计算PNJL模型中夸克有效分布函数的动量各向异性下完整Romatschke-Strickland形式"""
function quark_distribution_aniso(p_inv_fm::Float64, m_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, cosθ::Float64)
    # 计算各向异性修正后的能量
    E_aniso_inv_fm = sqrt(p_inv_fm^2 + m_inv_fm^2 + ξ * (p_inv_fm * cosθ)^2)
    
    # 计算分布函数
    return quark_distribution(E_aniso_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar)
end

# ------------------------------------
"""计算PNJL模型中反夸克有效分布函数的各向异性一阶修正项"""
function antiquark_distribution_aniso_correction(p_inv_fm::Float64, m_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, cosθ::Float64)
    E_inv_fm = sqrt(p_inv_fm^2 + m_inv_fm^2)
    df_dE = antiquark_df_dE(E_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar)

    return 0.5*ξ*(p_inv_fm*cosθ)^2 / E_inv_fm * df_dE
end

"""计算PNJL模型中反夸克有效分布函数的动量各向异性下完整Romatschke-Strickland形式"""
function antiquark_distribution_aniso(p_inv_fm::Float64, m_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, cosθ::Float64)
    # 计算各向异性修正后的能量
    E_aniso_inv_fm = sqrt(p_inv_fm^2 + m_inv_fm^2 + ξ * (p_inv_fm * cosθ)^2)
    
    # 计算分布函数
    return antiquark_distribution(E_aniso_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar)
end

"""根据符号选择计算夸克或反夸克有效分布函数的各向异性一阶修正项"""
function distribution_aniso_correction(sign_::Symbol, p_inv_fm::Float64, m_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, cosθ::Float64)
    if sign_ === :quark
        return quark_distribution_aniso_correction(p_inv_fm, m_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, ξ, cosθ)
    elseif sign_ === :antiquark
        return antiquark_distribution_aniso_correction(p_inv_fm, m_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, ξ, cosθ)
    else
        throw(ArgumentError("sign_ must be :quark or :antiquark"))
    end
end

"""根据符号选择计算夸克或反夸克有效分布函数的动量各向异性下完整Romatschke-Strickland形式"""
function distribution_aniso(sign_::Symbol, p_inv_fm::Float64, m_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, cosθ::Float64)
    if sign_ === :quark
        return quark_distribution_aniso(p_inv_fm, m_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, ξ, cosθ)
    elseif sign_ === :antiquark
        return antiquark_distribution_aniso(p_inv_fm, m_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar, ξ, cosθ)
    else
        throw(ArgumentError("sign_ must be :quark or :antiquark"))
    end
end

end # module PNJLQuarkDistributions_Aniso