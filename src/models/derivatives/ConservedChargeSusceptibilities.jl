"""
    ConservedChargeSusceptibilities

守恒荷静态涨落对象的首版实现。

当前阶段覆盖：
- baryon 方向 `chi_B(...; order=n)`
- baryon 累积量与比值量
- 二阶 `B/Q/S` 与 mixed susceptibilities
- flavor 基底下总压强对化学势的一阶/二阶导数

实现约定：
- 主路线使用 AD。
- 内部主对象优先取 `P(T, μ_u, μ_d, μ_s)` 对 flavor 化学势的导数。
- `χ_n` 的 `T` 缩放因子在导数完成后解析施加，避免直接对 `P/T^4` 做底层高阶 AD。
"""
module ConservedChargeSusceptibilities

using ForwardDiff
using StaticArrays

using ..PNJLCore: DEFAULT_MOMENTUM_COUNT, DEFAULT_THETA_COUNT
using ..HigherOrderDerivatives: nth_derivative, susceptibility_scale

export flavor_pressure_derivatives, conserved_charge_susceptibility
export chi_BQS, cumulant_BQS
export chi_B, chi1_B, chi2_B, chi3_B, chi4_B
export chi_Q, chi1_Q, chi2_Q, chi3_Q, chi4_Q
export chi_S, chi1_S, chi2_S, chi3_S, chi4_S
export chi11_BQ, chi11_BS, chi11_QS
export cumulant_B, baryon_Ssigma, baryon_kappa_sigma2

const _MODEL_CACHE = Dict{Symbol, Any}()
const _FLAVOR_IMPLICIT_GAP_SOLVER_CACHE = Dict{Tuple{Float64, Int, Int, UInt}, Any}()
const _BQS_TO_FLAVOR = @SMatrix [
    1.0 / 3.0   2.0 / 3.0   0.0
    1.0 / 3.0  -1.0 / 3.0   0.0
    1.0 / 3.0  -1.0 / 3.0  -1.0
]

@inline function _validate_baryon_inputs(T_fm::Real, order::Int)
    T_fm > 0 || throw(ArgumentError("T_fm must be positive, got $T_fm"))
    order >= 1 || throw(ArgumentError("order must be >= 1, got $order"))
    return nothing
end

@inline function _get_model(model_kind::Symbol)
    return get!(_MODEL_CACHE, model_kind) do
        Main.Models.create_model(model_kind)
    end
end

@inline function _flavor_solver(; model, xi::Real, p_num::Int, t_num::Int)
    key = (Float64(xi), Int(p_num), Int(t_num), Base.objectid(model))
    if haskey(_FLAVOR_IMPLICIT_GAP_SOLVER_CACHE, key)
        return _FLAVOR_IMPLICIT_GAP_SOLVER_CACHE[key]
    end

    f = Main.Models.create_flavor_mu_implicit_gap_solver(model; xi=xi, p_num=p_num, t_num=t_num)
    _FLAVOR_IMPLICIT_GAP_SOLVER_CACHE[key] = f
    return f
end

@inline _normalize_flavor_mu(mu_vec) = Main.Models.normalize_mu_vec(mu_vec)

@inline function _flavor_mu_from_bqs(muB_fm::Real, muQ_fm::Real, muS_fm::Real)
    return SVector(
        muB_fm / 3 + 2 * muQ_fm / 3,
        muB_fm / 3 - muQ_fm / 3,
        muB_fm / 3 - muQ_fm / 3 - muS_fm,
    )
end

@inline function _single_axis_order(orders::NTuple{3, Int})
    axis = 0
    order = 0
    for i in 1:3
        if orders[i] != 0
            axis == 0 || return nothing
            axis = i
            order = orders[i]
        end
    end
    axis == 0 && return nothing
    return (axis, order)
end

@inline function _flavor_direction_from_axis(axis::Int)
    1 <= axis <= 3 || throw(ArgumentError("axis must be 1(B), 2(Q), or 3(S), got $axis"))
    return SVector{3}(_BQS_TO_FLAVOR[1, axis], _BQS_TO_FLAVOR[2, axis], _BQS_TO_FLAVOR[3, axis])
end

@inline function _baryon_muT_scale(T_fm::Real, order::Int)
    return T_fm^(order - 4)
end

function _pressure_from_flavor_mu(
    T_fm::Real,
    mu_vec;
    xi::Real=0.0,
    p_num::Int=DEFAULT_MOMENTUM_COUNT,
    t_num::Int=DEFAULT_THETA_COUNT,
    model=nothing,
)
    T_fm > 0 || throw(ArgumentError("T_fm must be positive, got $T_fm"))
    μ = _normalize_flavor_mu(mu_vec)
    m = model === nothing ? _get_model(:PNJL) : model
    solver = _flavor_solver(model=m, xi=xi, p_num=p_num, t_num=t_num)
    θ = [T_fm, μ[1], μ[2], μ[3]]
    x, _ = solver(θ)
    x_sv = SVector{5}(Tuple(x))
    return Main.Models.model_pressure(m, x_sv, μ, T_fm; p_num=p_num, t_num=t_num, xi=xi)
end

function _single_axis_susceptibility(
    T_fm::Real,
    muB_fm::Real,
    muQ_fm::Real,
    muS_fm::Real,
    axis::Int,
    order::Int;
    xi::Real=0.0,
    p_num::Int=DEFAULT_MOMENTUM_COUNT,
    t_num::Int=DEFAULT_THETA_COUNT,
    model=nothing,
)
    T_fm > 0 || throw(ArgumentError("T_fm must be positive, got $T_fm"))
    order >= 1 || throw(ArgumentError("order must be >= 1, got $order"))

    μ0 = _flavor_mu_from_bqs(muB_fm, muQ_fm, muS_fm)
    direction = _flavor_direction_from_axis(axis)
    m = model === nothing ? _get_model(:PNJL) : model
    pressure_axis = δ -> _pressure_from_flavor_mu(T_fm, μ0 + δ * direction; xi=xi, p_num=p_num, t_num=t_num, model=m)
    δ0 = zero(promote_type(typeof(T_fm), typeof(muB_fm), typeof(muQ_fm), typeof(muS_fm)))
    return nth_derivative(pressure_axis, δ0, order) * susceptibility_scale(T_fm, order)
end

"""
    flavor_pressure_derivatives(T_fm, mu_vec; order=1|2, kwargs...) -> NamedTuple

计算总压强 `P(T, μ_u, μ_d, μ_s)` 对 flavor 化学势的导数。
"""
function flavor_pressure_derivatives(
    T_fm::Real,
    mu_vec;
    order::Int=1,
    xi::Real=0.0,
    p_num::Int=DEFAULT_MOMENTUM_COUNT,
    t_num::Int=DEFAULT_THETA_COUNT,
    model=nothing,
)
    T_fm > 0 || throw(ArgumentError("T_fm must be positive, got $T_fm"))
    μ0 = _normalize_flavor_mu(mu_vec)
    m = model === nothing ? _get_model(:PNJL) : model
    pressure_mu = μ -> _pressure_from_flavor_mu(T_fm, μ; xi=xi, p_num=p_num, t_num=t_num, model=m)

    pressure = pressure_mu(μ0)
    grad_mu = ForwardDiff.gradient(pressure_mu, μ0)

    if order == 1
        return (pressure=pressure, grad_mu=grad_mu)
    elseif order == 2
        hessian_mu = ForwardDiff.hessian(pressure_mu, μ0)
        return (pressure=pressure, grad_mu=grad_mu, hessian_mu=hessian_mu)
    else
        throw(ArgumentError("flavor_pressure_derivatives currently supports order=1 or 2, got $order"))
    end
end

@inline function _second_order_bqs_indices(orders::NTuple{3, Int})
    idxs = Int[]
    for (i, n) in enumerate(orders)
        for _ in 1:n
            push!(idxs, i)
        end
    end
    length(idxs) == 2 || throw(ArgumentError("orders=$orders is not a total second-order tuple"))
    return (idxs[1], idxs[2])
end

"""
    conserved_charge_susceptibility(T_fm, muB_fm, muQ_fm, muS_fm; orders=(i,j,k), kwargs...) -> Float64

返回 `χ_{ijk}^{BQS}` 的当前首版实现：
- 支持总二阶 `(2,0,0)`, `(0,2,0)`, `(0,0,2)`, `(1,1,0)`, `(1,0,1)`, `(0,1,1)`
- 支持纯 `B/Q/S` 方向 `χ_n, n=1..4`
"""
function conserved_charge_susceptibility(
    T_fm::Real,
    muB_fm::Real,
    muQ_fm::Real,
    muS_fm::Real;
    orders::NTuple{3, Int},
    xi::Real=0.0,
    p_num::Int=DEFAULT_MOMENTUM_COUNT,
    t_num::Int=DEFAULT_THETA_COUNT,
    model=nothing,
)
    total_order = orders[1] + orders[2] + orders[3]
    pure_axis = _single_axis_order(orders)

    if pure_axis !== nothing && 1 <= pure_axis[2] <= 4
        axis, order = pure_axis
        return _single_axis_susceptibility(T_fm, muB_fm, muQ_fm, muS_fm, axis, order; xi=xi, p_num=p_num, t_num=t_num, model=model)
    elseif total_order == 2
        μ = _flavor_mu_from_bqs(muB_fm, muQ_fm, muS_fm)
        derivs = flavor_pressure_derivatives(T_fm, μ; order=2, xi=xi, p_num=p_num, t_num=t_num, model=model)
        h_bqs = transpose(_BQS_TO_FLAVOR) * derivs.hessian_mu * _BQS_TO_FLAVOR
        i, j = _second_order_bqs_indices(orders)
        return h_bqs[i, j] * susceptibility_scale(T_fm, total_order)
    else
        throw(ArgumentError("conserved_charge_susceptibility currently supports total second order and pure single-axis B/Q/S order 1..4, got orders=$orders"))
    end
end

"""
    chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(i,j,k), kwargs...) -> Float64

统一的守恒荷广义磁化率接口。

示例：
```julia
chi2Q = Models.chi_BQS(T_fm, 0.0, 0.0, 0.0; orders=(0, 2, 0))
chi11BQ = Models.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(1, 1, 0))
```
"""
@inline function chi_BQS(
    T_fm::Real,
    muB_fm::Real,
    muQ_fm::Real,
    muS_fm::Real;
    orders::NTuple{3, Int},
    kwargs...
)
    return conserved_charge_susceptibility(T_fm, muB_fm, muQ_fm, muS_fm; orders=orders, kwargs...)
end

"""
    cumulant_BQS(T_fm, muB_fm, muQ_fm, muS_fm, V; orders=(i,j,k), kwargs...) -> Float64

统一的守恒荷累积量接口：
`C_{ijk}^{BQS} = V * T^3 * χ_{ijk}^{BQS}`。

示例：
```julia
c2Q = Models.cumulant_BQS(T_fm, 0.0, 0.0, 0.0, V; orders=(0, 2, 0))
c11BQ = Models.cumulant_BQS(T_fm, muB_fm, muQ_fm, muS_fm, V; orders=(1, 1, 0))
```
"""
@inline function cumulant_BQS(
    T_fm::Real,
    muB_fm::Real,
    muQ_fm::Real,
    muS_fm::Real,
    V::Real;
    orders::NTuple{3, Int},
    kwargs...
)
    T_fm > 0 || throw(ArgumentError("T_fm must be positive, got $T_fm"))
    return V * T_fm^3 * chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=orders, kwargs...)
end

"""
    chi_B(T_fm, muB_fm; order=1, kwargs...)

返回 baryon 方向广义磁化率 chi_n^B。

定义对象是 (P / T^4) 对 (mu_B / T) 的 n 阶导数。
当前实现建立在对称化学势路径 mu_u = mu_d = mu_s = mu_B / 3 之上，
因此适用于 M1 baryon 首版静态涨落 API。
"""
function chi_B(
    T_fm::Real,
    muB_fm::Real;
    order::Int=1,
    xi::Real=0.0,
    p_num::Int=DEFAULT_MOMENTUM_COUNT,
    t_num::Int=DEFAULT_THETA_COUNT,
    model=nothing,
)
    _validate_baryon_inputs(T_fm, order)
    return _single_axis_susceptibility(T_fm, muB_fm, zero(muB_fm), zero(muB_fm), 1, order; xi=xi, p_num=p_num, t_num=t_num, model=model)
end

@inline chi1_B(T_fm::Real, muB_fm::Real; kwargs...) = chi_B(T_fm, muB_fm; order=1, kwargs...)
@inline chi2_B(T_fm::Real, muB_fm::Real; kwargs...) = chi_B(T_fm, muB_fm; order=2, kwargs...)
@inline chi3_B(T_fm::Real, muB_fm::Real; kwargs...) = chi_B(T_fm, muB_fm; order=3, kwargs...)
@inline chi4_B(T_fm::Real, muB_fm::Real; kwargs...) = chi_B(T_fm, muB_fm; order=4, kwargs...)
@inline chi_Q(T_fm::Real, muB_fm::Real=0.0, muQ_fm::Real=0.0, muS_fm::Real=0.0; order::Int=1, kwargs...) = conserved_charge_susceptibility(T_fm, muB_fm, muQ_fm, muS_fm; orders=(0, order, 0), kwargs...)
@inline chi1_Q(T_fm::Real, muB_fm::Real=0.0, muQ_fm::Real=0.0, muS_fm::Real=0.0; kwargs...) = chi_Q(T_fm, muB_fm, muQ_fm, muS_fm; order=1, kwargs...)
@inline chi2_Q(T_fm::Real, muB_fm::Real=0.0, muQ_fm::Real=0.0, muS_fm::Real=0.0; kwargs...) = chi_Q(T_fm, muB_fm, muQ_fm, muS_fm; order=2, kwargs...)
@inline chi3_Q(T_fm::Real, muB_fm::Real=0.0, muQ_fm::Real=0.0, muS_fm::Real=0.0; kwargs...) = chi_Q(T_fm, muB_fm, muQ_fm, muS_fm; order=3, kwargs...)
@inline chi4_Q(T_fm::Real, muB_fm::Real=0.0, muQ_fm::Real=0.0, muS_fm::Real=0.0; kwargs...) = chi_Q(T_fm, muB_fm, muQ_fm, muS_fm; order=4, kwargs...)
@inline chi_S(T_fm::Real, muB_fm::Real=0.0, muQ_fm::Real=0.0, muS_fm::Real=0.0; order::Int=1, kwargs...) = conserved_charge_susceptibility(T_fm, muB_fm, muQ_fm, muS_fm; orders=(0, 0, order), kwargs...)
@inline chi1_S(T_fm::Real, muB_fm::Real=0.0, muQ_fm::Real=0.0, muS_fm::Real=0.0; kwargs...) = chi_S(T_fm, muB_fm, muQ_fm, muS_fm; order=1, kwargs...)
@inline chi2_S(T_fm::Real, muB_fm::Real=0.0, muQ_fm::Real=0.0, muS_fm::Real=0.0; kwargs...) = chi_S(T_fm, muB_fm, muQ_fm, muS_fm; order=2, kwargs...)
@inline chi3_S(T_fm::Real, muB_fm::Real=0.0, muQ_fm::Real=0.0, muS_fm::Real=0.0; kwargs...) = chi_S(T_fm, muB_fm, muQ_fm, muS_fm; order=3, kwargs...)
@inline chi4_S(T_fm::Real, muB_fm::Real=0.0, muQ_fm::Real=0.0, muS_fm::Real=0.0; kwargs...) = chi_S(T_fm, muB_fm, muQ_fm, muS_fm; order=4, kwargs...)
@inline chi11_BQ(T_fm::Real, muB_fm::Real=0.0, muQ_fm::Real=0.0, muS_fm::Real=0.0; kwargs...) = conserved_charge_susceptibility(T_fm, muB_fm, muQ_fm, muS_fm; orders=(1, 1, 0), kwargs...)
@inline chi11_BS(T_fm::Real, muB_fm::Real=0.0, muQ_fm::Real=0.0, muS_fm::Real=0.0; kwargs...) = conserved_charge_susceptibility(T_fm, muB_fm, muQ_fm, muS_fm; orders=(1, 0, 1), kwargs...)
@inline chi11_QS(T_fm::Real, muB_fm::Real=0.0, muQ_fm::Real=0.0, muS_fm::Real=0.0; kwargs...) = conserved_charge_susceptibility(T_fm, muB_fm, muQ_fm, muS_fm; orders=(0, 1, 1), kwargs...)

"""
    cumulant_B(T_fm, muB_fm, V; order=1, kwargs...)

返回 baryon 方向累积量 C_n^B = V * T^3 * chi_n^B。
"""
function cumulant_B(
    T_fm::Real,
    muB_fm::Real,
    V::Real;
    order::Int=1,
    kwargs...
)
    _validate_baryon_inputs(T_fm, order)
    return cumulant_BQS(T_fm, muB_fm, zero(muB_fm), zero(muB_fm), V; orders=(order, 0, 0), kwargs...)
end

"""
    baryon_Ssigma(T_fm, muB_fm; kwargs...)

返回 Ssigma = chi_3^B / chi_2^B。
"""
function baryon_Ssigma(T_fm::Real, muB_fm::Real; kwargs...)
    chi2 = chi2_B(T_fm, muB_fm; kwargs...)
    chi3 = chi3_B(T_fm, muB_fm; kwargs...)
    return chi3 / chi2
end

"""
    baryon_kappa_sigma2(T_fm, muB_fm; kwargs...)

返回 kappa_sigma2 = chi_4^B / chi_2^B。
"""
function baryon_kappa_sigma2(T_fm::Real, muB_fm::Real; kwargs...)
    chi2 = chi2_B(T_fm, muB_fm; kwargs...)
    chi4 = chi4_B(T_fm, muB_fm; kwargs...)
    return chi4 / chi2
end

end # module ConservedChargeSusceptibilities
