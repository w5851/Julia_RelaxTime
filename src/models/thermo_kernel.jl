"""thermo_kernel.jl

主域公共热力学计算骨架（模型无关，基于多重派发的 `omega` 入口）。

约定：
- 压强 `P = -Ω`
- `rho` 通过 `∂P/∂μ`（ForwardDiff）
- `thermo` 返回 `(pressure, rho_norm, entropy, energy)`
"""

using StaticArrays
using ForwardDiff

export model_pressure, model_rho, model_thermo

@inline function _apply_pressure_shift(
    pressure::T,
    pressure_shift_fn,
    model,
    x_state,
    mu_vec,
    T_fm,
) where {T}
    pressure_shift_fn === nothing && return pressure
    return pressure + pressure_shift_fn(model, x_state, mu_vec, T_fm)
end

@inline function _rho0_ref()
    if isdefined(Main, :Constants_PNJL) && isdefined(Main.Constants_PNJL, :ρ0_inv_fm3)
        return getproperty(Main.Constants_PNJL, :ρ0_inv_fm3)
    end
    return 0.16
end

"""基于模型 `omega` 的统一压强入口：`P = -Ω`。"""
@inline function model_pressure(
    model::AbstractQCDModel,
    x_state,
    mu_vec,
    T_fm;
    p_num::Int=24,
    t_num::Int=8,
    xi=0.0,
    pressure_shift_fn=nothing,
    kwargs...
)
    pressure = -omega(model, x_state, T_fm, mu_vec; p_num=p_num, t_num=t_num, xi=xi, kwargs...)
    return _apply_pressure_shift(pressure, pressure_shift_fn, model, x_state, mu_vec, T_fm)
end

"""统一数密度向量入口：`ρ_i = ∂P/∂μ_i`。"""
function model_rho(
    model::AbstractQCDModel,
    x_state,
    mu_vec,
    T_fm;
    p_num::Int=24,
    t_num::Int=8,
    xi=0.0,
    pressure_shift_fn=nothing,
    kwargs...
)
    μ0 = normalize_mu_vec(mu_vec)
    pressure_mu = μ -> model_pressure(model, x_state, μ, T_fm; p_num=p_num, t_num=t_num, xi=xi, pressure_shift_fn=pressure_shift_fn, kwargs...)
    grad = ForwardDiff.gradient(pressure_mu, μ0)
    grad_type = typeof(grad[1])
    return SVector{3, grad_type}(Tuple(grad))
end

"""统一热力学量入口：返回 `(pressure, rho_norm, entropy, energy)`。"""
function model_thermo(
    model::AbstractQCDModel,
    x_state,
    mu_vec,
    T_fm;
    p_num::Int=24,
    t_num::Int=8,
    xi=0.0,
    pressure_shift_fn=nothing,
    kwargs...
)
    require_capability(model, :model_thermo)
    ρ0 = _rho0_ref()
    μ0 = normalize_mu_vec(mu_vec)
    rho_vec = model_rho(model, x_state, μ0, T_fm; p_num=p_num, t_num=t_num, xi=xi, pressure_shift_fn=pressure_shift_fn, kwargs...)
    rho_norm = sum(rho_vec) / (3.0 * ρ0)

    pressure_T = τ -> model_pressure(model, x_state, μ0, τ; p_num=p_num, t_num=t_num, xi=xi, pressure_shift_fn=pressure_shift_fn, kwargs...)
    entropy = ForwardDiff.derivative(pressure_T, T_fm)

    pressure = model_pressure(model, x_state, μ0, T_fm; p_num=p_num, t_num=t_num, xi=xi, pressure_shift_fn=pressure_shift_fn, kwargs...)
    energy = -pressure + sum(μ0 .* rho_vec) + T_fm * entropy

    return pressure, rho_norm, entropy, energy
end
