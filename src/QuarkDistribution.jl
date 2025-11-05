module PNJLQuarkDistributions
export quark_distribution, antiquark_distribution, 
        quark_distribution_antiderivative, quark_distribution_integral,
        antiquark_distribution_antiderivative, antiquark_distribution_integral
export distribution
      
"""PNJL模型中夸克有效分布函数"""
@fastmath function quark_distribution(E_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64)
    # 计算温度的倒数
    β_fm = 1 / T_inv_fm

    # 计算指数项
    exp_term = clamp(exp(-(E_inv_fm - μ_inv_fm) * β_fm), 1e-200, 1e200)
    exp_term2 = exp_term*exp_term
    exp_term3 = exp_term2*exp_term

    # 计算分子
    numerator = Φ * exp_term + 2 * Φbar * exp_term2 + exp_term3

    # 计算分母
    denominator = 1 + 3 * Φ * exp_term + 3 * Φbar * exp_term2 + exp_term3

    # 返回分布函数
    return numerator / denominator
end

"""PNJL模型中反夸克有效分布函数"""
@fastmath function antiquark_distribution(E_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64)
    # 计算温度的倒数
    β_fm = 1 / T_inv_fm

    # 计算指数项
    exp_term = clamp(exp(-(E_inv_fm + μ_inv_fm) * β_fm), 1e-200, 1e200)
    exp_term2 = exp_term * exp_term
    exp_term3 = exp_term2 * exp_term

    # 计算分子
    numerator = Φbar * exp_term + 2 * Φ * exp_term2 + exp_term3

    # 计算分母
    denominator = 1 + 3 * Φbar * exp_term + 3 * Φ * exp_term2 + exp_term3

    # 返回分布函数
    return numerator / denominator
end

# -----------------------
"""PNJL模型中夸克有效分布函数的原函数"""
function quark_distribution_antiderivative(E_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64)
    # 计算温度的倒数（β）
    β_fm = 1 / T_inv_fm

    # 计算指数项
    exp_term = clamp(exp(-(E_inv_fm - μ_inv_fm) * β_fm), 1e-200, 1e200)
    exp_term2 = exp_term * exp_term
    exp_term3 = exp_term2 * exp_term

    # 计算对数项
    log_term = 1 + 3 * Φ * exp_term + 3 * Φbar * exp_term2 + exp_term3

    # 返回原函数（针对能量 E 的不定积分）, 常数因子选择使得 d/dE antiderivative = quark_distribution
    return -T_inv_fm * log(log_term) / 3
end

"""计算夸克有效分布函数在给定能量区间的积分"""
@inline function quark_distribution_integral(E_min_inv_fm::Float64, E_max_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64)
    return quark_distribution_antiderivative(E_max_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar) - 
           quark_distribution_antiderivative(E_min_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar)
end

# ----------------------
"""PNJL模型中反夸克有效分布函数的原函数"""
function antiquark_distribution_antiderivative(E_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64)
    # 计算温度的倒数（β）
    β_fm = 1 / T_inv_fm

    # 计算指数项
    exp_term = clamp(exp(-(E_inv_fm + μ_inv_fm) * β_fm), 1e-200, 1e200)
    exp_term2 = exp_term * exp_term
    exp_term3 = exp_term2 * exp_term

    # 计算对数项
    log_term = 1 + 3 * Φbar * exp_term + 3 * Φ * exp_term2 + exp_term3

    # 返回原函数（针对能量 E 的不定积分）, 常数因子选择使得 d/dE antiderivative = antiquark_distribution
    return -T_inv_fm * log(log_term) / 3
end

"""计算反夸克有效分布函数在给定能量区间的积分"""
@inline function antiquark_distribution_integral(E_min_inv_fm::Float64, E_max_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64)
    return antiquark_distribution_antiderivative(E_max_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar) - 
           antiquark_distribution_antiderivative(E_min_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar)
end

"""根据符号选择计算夸克或反夸克有效分布函数"""
function distribution(sign_::Symbol, E_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64)
    if sign_ === :quark
        return quark_distribution(E_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar)
    elseif sign_ === :antiquark
        return antiquark_distribution(E_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar)
    else
        throw(ArgumentError("sign_ must be :quark or :antiquark"))
    end
end
end # module PNJLQuarkDistributions


module NJLQuarkDistributions
export quark_distribution, antiquark_distribution

function quark_distribution(E_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64)
    β_fm = 1 / T_inv_fm
    exp_term = clamp(exp((E_inv_fm - μ_inv_fm) * β_fm), 1e-200, 1e200)
    return 1.0 / (1.0 + exp_term)
end

function antiquark_distribution(E_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64)
    β_fm = 1 / T_inv_fm
    exp_term = clamp(exp((-E_inv_fm - μ_inv_fm) * β_fm), 1e-200, 1e200)
    return 1.0 / (1.0 + exp_term)
end

end # module NJLQuarkDistributions