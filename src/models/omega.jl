"""omega.jl

Ω (grand potential) 组装入口（新架构）。

约定：
- `x_state` 的前 3 个分量是凝聚参数 φ = (φu, φd, φs)
- 若 `x_state` 长度 ≥ 5：第 4、5 个分量分别为 Φ、Φbar
- 若 `x_state` 仅有 3 个分量：默认 Φ=Φbar=1（适合 NJL）
"""

using StaticArrays

export omega, omega_components, grand_potential

@inline function _extract_state(x_state)
    st = meanfield_state(x_state)
    return st.phi, st.Phi, st.PhiBar
end

"""返回 (chi, poly, vac, therm, masses, omega) 的 NamedTuple。"""
function omega_components(
    model::AbstractQCDModel,
    x_state,
    T,
    mu_vec;
    p_num::Int=64,
    t_num::Int=8,
    xi=0.0,
    kwargs...
)
    φ, Φ, Φbar = _extract_state(x_state)
    mu_vec = normalize_mu_vec(mu_vec)

    masses = calculate_mass_vec(model, φ)
    vac = vacuum_contribution(model, masses; x_state=x_state, T=T, mu_vec=mu_vec, Phi=Φ, PhiBar=Φbar, kwargs...)
    poly = polyakov_potential(model, Φ, Φbar, T; kwargs...)
    therm = thermal_contribution(model, masses, Φ, Φbar, mu_vec, T; p_num=p_num, t_num=t_num, xi=xi, kwargs...)
    chi = calculate_chiral(model, φ)

    ω = chi + poly + vac + therm
    return (chi=chi, poly=poly, vac=vac, therm=therm, masses=masses, omega=ω)
end

"""返回 Ω 的标量值。"""
function omega(
    model::AbstractQCDModel,
    x_state,
    T,
    mu_vec;
    p_num::Int=64,
    t_num::Int=8,
    xi=0.0,
    kwargs...
)
    return omega_components(model, x_state, T, mu_vec; p_num=p_num, t_num=t_num, xi=xi, kwargs...).omega
end

"""别名：grand_potential ≡ omega。"""
grand_potential(args...; kwargs...) = omega(args...; kwargs...)
