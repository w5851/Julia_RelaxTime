module PNJLQuarkDistributions
export quark_distribution, antiquark_distribution, 
        quark_distribution_antiderivative, quark_distribution_integral,
        antiquark_distribution_antiderivative, antiquark_distribution_integral
export distribution

@inline function _safe_log_term(x::Real; min_val::Real=1e-300)
    return log(ifelse(x <= min_val, oftype(x, min_val), x))
end

@inline function _pnjl_quark_distribution_core(x::Real, Φ::Real, Φbar::Real)
    if x >= 0.0
        z = exp(-x)
        z2 = z * z
        z3 = z2 * z
        numerator = Φ * z + 2.0 * Φbar * z2 + z3
        denominator = 1.0 + 3.0 * Φ * z + 3.0 * Φbar * z2 + z3
        return numerator / denominator
    end

    y = exp(x)
    y2 = y * y
    y3 = y2 * y
    numerator = Φ * y2 + 2.0 * Φbar * y + 1.0
    denominator = y3 + 3.0 * Φ * y2 + 3.0 * Φbar * y + 1.0
    return numerator / denominator
end

@inline function _pnjl_antiquark_distribution_core(x::Real, Φ::Real, Φbar::Real)
    if x >= 0.0
        z = exp(-x)
        z2 = z * z
        z3 = z2 * z
        numerator = Φbar * z + 2.0 * Φ * z2 + z3
        denominator = 1.0 + 3.0 * Φbar * z + 3.0 * Φ * z2 + z3
        return numerator / denominator
    end

    y = exp(x)
    y2 = y * y
    y3 = y2 * y
    numerator = Φbar * y2 + 2.0 * Φ * y + 1.0
    denominator = y3 + 3.0 * Φbar * y2 + 3.0 * Φ * y + 1.0
    return numerator / denominator
end
      
"""PNJL模型中夸克有效分布函数
注：也许将硬截断clamp改为分子分母同时除以可能出现的最大项exp_term3会更好(确保计算过程中不出现很大的数，但是可能会出现很小的数)，
但目前这样足够稳定且高效"""
@fastmath function quark_distribution(E_inv_fm::Real, μ_inv_fm::Real, T_inv_fm::Real, Φ::Real, Φbar::Real)
    x = (E_inv_fm - μ_inv_fm) / T_inv_fm
    return _pnjl_quark_distribution_core(x, Φ, Φbar)
end

"""PNJL模型中反夸克有效分布函数"""
@fastmath function antiquark_distribution(E_inv_fm::Real, μ_inv_fm::Real, T_inv_fm::Real, Φ::Real, Φbar::Real)
    x = (E_inv_fm + μ_inv_fm) / T_inv_fm
    return _pnjl_antiquark_distribution_core(x, Φ, Φbar)
end

# -----------------------
"""PNJL模型中夸克有效分布函数的原函数"""
function quark_distribution_antiderivative(E_inv_fm::Real, μ_inv_fm::Real, T_inv_fm::Real, Φ::Real, Φbar::Real)
    # 计算温度的倒数（β）
    β_fm = 1 / T_inv_fm

    # 计算指数项
    exp_term = clamp(exp(-(E_inv_fm - μ_inv_fm) * β_fm), 1e-200, 1e200)
    exp_term2 = exp_term * exp_term
    exp_term3 = exp_term2 * exp_term

    # 计算对数项
    log_term = 1 + 3 * Φ * exp_term + 3 * Φbar * exp_term2 + exp_term3

    # 返回原函数（针对能量 E 的不定积分）, 常数因子选择使得 d/dE antiderivative = quark_distribution
    return -T_inv_fm * _safe_log_term(log_term) / 3
end

"""计算夸克有效分布函数在给定能量区间的积分"""
@inline function quark_distribution_integral(E_min_inv_fm::Real, E_max_inv_fm::Real, μ_inv_fm::Real, T_inv_fm::Real, Φ::Real, Φbar::Real)
    return quark_distribution_antiderivative(E_max_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar) - 
           quark_distribution_antiderivative(E_min_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar)
end

# ----------------------
"""PNJL模型中反夸克有效分布函数的原函数"""
function antiquark_distribution_antiderivative(E_inv_fm::Real, μ_inv_fm::Real, T_inv_fm::Real, Φ::Real, Φbar::Real)
    # 计算温度的倒数（β）
    β_fm = 1 / T_inv_fm

    # 计算指数项
    exp_term = clamp(exp(-(E_inv_fm + μ_inv_fm) * β_fm), 1e-200, 1e200)
    exp_term2 = exp_term * exp_term
    exp_term3 = exp_term2 * exp_term

    # 计算对数项
    log_term = 1 + 3 * Φbar * exp_term + 3 * Φ * exp_term2 + exp_term3

    # 返回原函数（针对能量 E 的不定积分）, 常数因子选择使得 d/dE antiderivative = antiquark_distribution
    return -T_inv_fm * _safe_log_term(log_term) / 3
end

"""计算反夸克有效分布函数在给定能量区间的积分"""
@inline function antiquark_distribution_integral(E_min_inv_fm::Real, E_max_inv_fm::Real, μ_inv_fm::Real, T_inv_fm::Real, Φ::Real, Φbar::Real)
    return antiquark_distribution_antiderivative(E_max_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar) - 
           antiquark_distribution_antiderivative(E_min_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar)
end

"""根据符号选择计算夸克或反夸克有效分布函数"""
function distribution(sign_::Symbol, E_inv_fm::Real, μ_inv_fm::Real, T_inv_fm::Real, Φ::Real, Φbar::Real)
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
