"""PNJLDistributions (models-side)

为 models 子系统提供最小的 PNJL 分布函数 / 各向异性 RS 形式接口，供输运 provider 复用。

约定：单位均为 fm⁻¹。

注意：这是“最小 API 抽取”切片；实现目前与仓库既有公式保持一致，
用于让 `Models.transport_provider(:models)` 不再依赖根目录的 legacy 分布模块命名。
"""

@inline function _clamp_exp(x)
    return clamp(exp(x), 1e-200, 1e200)
end

"""Polyakov-aware effective quark distribution f_q(E; μ,T,Φ,Φbar)."""
@fastmath function pnjl_quark_distribution(E_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64)
    β_fm = 1.0 / T_inv_fm
    exp_term = _clamp_exp(-(E_inv_fm - μ_inv_fm) * β_fm)
    exp_term2 = exp_term * exp_term
    exp_term3 = exp_term2 * exp_term

    numerator = Φ * exp_term + 2.0 * Φbar * exp_term2 + exp_term3
    denominator = 1.0 + 3.0 * Φ * exp_term + 3.0 * Φbar * exp_term2 + exp_term3
    return numerator / denominator
end

"""Polyakov-aware effective antiquark distribution f_\u0304q(E; μ,T,Φ,Φbar)."""
@fastmath function pnjl_antiquark_distribution(E_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64)
    β_fm = 1.0 / T_inv_fm
    exp_term = _clamp_exp(-(E_inv_fm + μ_inv_fm) * β_fm)
    exp_term2 = exp_term * exp_term
    exp_term3 = exp_term2 * exp_term

    numerator = Φbar * exp_term + 2.0 * Φ * exp_term2 + exp_term3
    denominator = 1.0 + 3.0 * Φbar * exp_term + 3.0 * Φ * exp_term2 + exp_term3
    return numerator / denominator
end

"""Romatschke–Strickland form: f_q(p,cosθ; ξ) via E_ξ = sqrt(p^2+m^2+ξ(p cosθ)^2)."""
@inline function pnjl_quark_distribution_aniso(
    p_inv_fm::Float64,
    m_inv_fm::Float64,
    μ_inv_fm::Float64,
    T_inv_fm::Float64,
    Φ::Float64,
    Φbar::Float64,
    ξ::Float64,
    cosθ::Float64,
)
    E_aniso_inv_fm = sqrt(p_inv_fm * p_inv_fm + m_inv_fm * m_inv_fm + ξ * (p_inv_fm * cosθ)^2)
    return pnjl_quark_distribution(E_aniso_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar)
end

"""Romatschke–Strickland form: f_\u0304q(p,cosθ; ξ)."""
@inline function pnjl_antiquark_distribution_aniso(
    p_inv_fm::Float64,
    m_inv_fm::Float64,
    μ_inv_fm::Float64,
    T_inv_fm::Float64,
    Φ::Float64,
    Φbar::Float64,
    ξ::Float64,
    cosθ::Float64,
)
    E_aniso_inv_fm = sqrt(p_inv_fm * p_inv_fm + m_inv_fm * m_inv_fm + ξ * (p_inv_fm * cosθ)^2)
    return pnjl_antiquark_distribution(E_aniso_inv_fm, μ_inv_fm, T_inv_fm, Φ, Φbar)
end
