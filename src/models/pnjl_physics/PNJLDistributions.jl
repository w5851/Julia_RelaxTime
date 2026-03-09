"""PNJLDistributions (models-side)

为 models 子系统提供最小的 PNJL 分布函数 / 各向异性 RS 形式接口，供输运 provider 复用。

约定：单位均为 fm⁻¹。

注意：这是“最小 API 抽取”切片；实现目前与仓库既有公式保持一致，
用于让 `Models.transport_provider(model)` 不再依赖根目录的 legacy 分布模块命名。
"""

@inline function _pnjl_quark_distribution_core(x::Float64, Φ::Float64, Φbar::Float64)
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

@inline function _pnjl_antiquark_distribution_core(x::Float64, Φ::Float64, Φbar::Float64)
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

"""Polyakov-aware effective quark distribution f_q(E; μ,T,Φ,Φbar)."""
@fastmath function pnjl_quark_distribution(E_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64)
    x = (E_inv_fm - μ_inv_fm) / T_inv_fm
    return _pnjl_quark_distribution_core(x, Φ, Φbar)
end

"""Polyakov-aware effective antiquark distribution f_\u0304q(E; μ,T,Φ,Φbar)."""
@fastmath function pnjl_antiquark_distribution(E_inv_fm::Float64, μ_inv_fm::Float64, T_inv_fm::Float64, Φ::Float64, Φbar::Float64)
    x = (E_inv_fm + μ_inv_fm) / T_inv_fm
    return _pnjl_antiquark_distribution_core(x, Φ, Φbar)
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
