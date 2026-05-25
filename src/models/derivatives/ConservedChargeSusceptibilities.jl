"""
    ConservedChargeSusceptibilities

守恒荷静态涨落对象的首版实现。

当前阶段覆盖：
- `B/Q/S` 单轴 `chi_*(...; order=n)`
- baryon 累积量与比值量
- 二阶 `B/Q/S` 与任意总阶数的 mixed BQS susceptibilities
- flavor 基底下总压强对化学势的一阶/二阶导数

实现约定：
- 主路线使用 AD/TaylorDiff。
- 内部主对象优先取 `P(T, μ_u, μ_d, μ_s)` 对 flavor 化学势的导数。
- `χ_n` 的 `T` 缩放因子在导数完成后解析施加，避免直接对 `P/T^4` 做底层高阶 AD。
"""
module ConservedChargeSusceptibilities

using StaticArrays

using ..PNJLCore: DEFAULT_MOMENTUM_COUNT, DEFAULT_THETA_COUNT
using ..PNJLChiBTaylorDiff: chi_BQS_mixed_taylorjet, chi_direction_taylordiff
using ..PNJLChiBTaylorDiff: nth_derivative_from_series, pressure_series_direction
using ..Models: create_model, normalize_mu_vec

export flavor_pressure_derivatives, conserved_charge_susceptibility
export chi_BQS, cumulant_BQS
export chi_B, chi1_B, chi2_B, chi3_B, chi4_B
export chi_Q, chi1_Q, chi2_Q, chi3_Q, chi4_Q
export chi_S, chi1_S, chi2_S, chi3_S, chi4_S
export chi11_BQ, chi11_BS, chi11_QS
export cumulant_B, baryon_Ssigma, baryon_kappa_sigma2

const _MODEL_CACHE = Dict{Symbol, Any}()
const _BQS_TO_FLAVOR = @SMatrix [
    1.0 / 3.0   2.0 / 3.0   0.0
    1.0 / 3.0  -1.0 / 3.0   0.0
    1.0 / 3.0  -1.0 / 3.0  -1.0
]
const _DERIVATIVE_BACKENDS = (:auto, :taylordiff, :mixedjet)

@inline function _validate_baryon_inputs(T_fm::Real, order::Int)
    T_fm > 0 || throw(ArgumentError("T_fm must be positive, got $T_fm"))
    order >= 1 || throw(ArgumentError("order must be >= 1, got $order"))
    return nothing
end

@inline function _validate_derivative_backend(derivative_backend::Symbol)
    derivative_backend === :forwarddiff && throw(ArgumentError("derivative_backend=:forwarddiff has been retired from conserved charge susceptibilities; use derivative_backend=:auto, :taylordiff, or :mixedjet. Single-axis B/Q/S routes use TaylorDiff fast path and mixed BQS routes use MixedTaylorJet."))
    derivative_backend in _DERIVATIVE_BACKENDS && return derivative_backend
    throw(ArgumentError("derivative_backend must be one of $(_DERIVATIVE_BACKENDS), got $(derivative_backend)"))
end

@inline function _resolve_single_axis_backend(axis::Int, order::Int, derivative_backend::Symbol)
    _ = axis
    _ = order
    backend = _validate_derivative_backend(derivative_backend)
    backend === :mixedjet && return :taylordiff
    return backend === :auto ? :taylordiff : backend
end

@inline function _get_model(model_kind::Symbol)
    return get!(_MODEL_CACHE, model_kind) do
        create_model(model_kind)
    end
end

@inline _normalize_flavor_mu(mu_vec) = normalize_mu_vec(mu_vec)

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
    derivative_backend::Symbol=:auto,
    series_iterations::Union{Nothing, Int}=nothing,
    linear_solve::Symbol=:auto,
    series_residual_tol::Real=1e-7,
)
    T_fm > 0 || throw(ArgumentError("T_fm must be positive, got $T_fm"))
    order >= 1 || throw(ArgumentError("order must be >= 1, got $order"))

    backend = _resolve_single_axis_backend(axis, order, derivative_backend)
    m = model === nothing ? _get_model(:PNJL) : model
    backend === :taylordiff || error("unreachable derivative backend: $backend")
    base_mu_vec = _flavor_mu_from_bqs(muB_fm, muQ_fm, muS_fm)
    direction = _flavor_direction_from_axis(axis)
    return chi_direction_taylordiff(
        T_fm,
        base_mu_vec,
        direction;
        order=order,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        model=m,
        series_iterations=series_iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )
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
    derivative_backend::Symbol=:auto,
    series_iterations::Union{Nothing, Int}=nothing,
    linear_solve::Symbol=:auto,
    series_residual_tol::Real=1e-7,
)
    T_fm > 0 || throw(ArgumentError("T_fm must be positive, got $T_fm"))
    μ0 = _normalize_flavor_mu(mu_vec)
    m = model === nothing ? _get_model(:PNJL) : model
    _validate_derivative_backend(derivative_backend)
    order == 1 || order == 2 || throw(ArgumentError("flavor_pressure_derivatives currently supports order=1 or 2, got $order"))

    e1 = SVector{3}(1.0, 0.0, 0.0)
    e2 = SVector{3}(0.0, 1.0, 0.0)
    e3 = SVector{3}(0.0, 0.0, 1.0)
    series_kwargs = (; xi=xi, p_num=p_num, t_num=t_num, model=m, series_iterations=series_iterations, linear_solve=linear_solve, series_residual_tol=series_residual_tol)

    s1 = pressure_series_direction(T_fm, μ0, e1; order=order, series_kwargs...)
    s2 = pressure_series_direction(T_fm, μ0, e2; order=order, series_kwargs...)
    s3 = pressure_series_direction(T_fm, μ0, e3; order=order, series_kwargs...)

    pressure = nth_derivative_from_series(s1, 0)
    grad_mu = Vector{Float64}(undef, 3)
    grad_mu[1] = nth_derivative_from_series(s1, 1)
    grad_mu[2] = nth_derivative_from_series(s2, 1)
    grad_mu[3] = nth_derivative_from_series(s3, 1)

    if order == 1
        return (pressure=pressure, grad_mu=grad_mu)
    end

    hessian_mu = Matrix{Float64}(undef, 3, 3)
    diag11 = nth_derivative_from_series(s1, 2)
    diag22 = nth_derivative_from_series(s2, 2)
    diag33 = nth_derivative_from_series(s3, 2)
    hessian_mu[1, 1] = diag11
    hessian_mu[2, 2] = diag22
    hessian_mu[3, 3] = diag33

    s12 = pressure_series_direction(T_fm, μ0, e1 + e2; order=order, series_kwargs...)
    s13 = pressure_series_direction(T_fm, μ0, e1 + e3; order=order, series_kwargs...)
    s23 = pressure_series_direction(T_fm, μ0, e2 + e3; order=order, series_kwargs...)
    hessian_mu[1, 2] = hessian_mu[2, 1] = (nth_derivative_from_series(s12, 2) - diag11 - diag22) / 2
    hessian_mu[1, 3] = hessian_mu[3, 1] = (nth_derivative_from_series(s13, 2) - diag11 - diag33) / 2
    hessian_mu[2, 3] = hessian_mu[3, 2] = (nth_derivative_from_series(s23, 2) - diag22 - diag33) / 2
    return (pressure=pressure, grad_mu=grad_mu, hessian_mu=hessian_mu)
end

"""
    conserved_charge_susceptibility(T_fm, muB_fm, muQ_fm, muS_fm; orders=(i,j,k), kwargs...) -> Float64

返回 `χ_{ijk}^{BQS}`：
- 纯 `B/Q/S` 单方向由 `:auto`/`:taylordiff` 路由到单变量 TaylorDiff fast path
- mixed BQS 组合由 `:auto`/`:taylordiff`/`:mixedjet` 路由到内部 multivariate Taylor jet
- `:forwarddiff` legacy reference/fallback 已下线，会抛出迁移错误
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
    derivative_backend::Symbol=:auto,
    series_iterations::Union{Nothing, Int}=nothing,
    linear_solve::Symbol=:auto,
    series_residual_tol::Real=1e-7,
)
    pure_axis = _single_axis_order(orders)

    if pure_axis !== nothing
        axis, order = pure_axis
        backend = _resolve_single_axis_backend(axis, order, derivative_backend)
        return _single_axis_susceptibility(
            T_fm,
            muB_fm,
            muQ_fm,
            muS_fm,
            axis,
            order;
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            model=model,
            derivative_backend=backend,
            series_iterations=series_iterations,
            linear_solve=linear_solve,
            series_residual_tol=series_residual_tol,
        )
    else
        _validate_derivative_backend(derivative_backend)
        return chi_BQS_mixed_taylorjet(
            T_fm,
            muB_fm,
            muQ_fm,
            muS_fm;
            orders=orders,
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            model=model,
            series_iterations=series_iterations,
            linear_solve=linear_solve,
            series_residual_tol=series_residual_tol,
        )
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
    derivative_backend::Symbol=:auto,
    series_iterations::Union{Nothing, Int}=nothing,
    linear_solve::Symbol=:auto,
    series_residual_tol::Real=1e-7,
)
    _validate_baryon_inputs(T_fm, order)
    return _single_axis_susceptibility(
        T_fm,
        muB_fm,
        zero(muB_fm),
        zero(muB_fm),
        1,
        order;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        model=model,
        derivative_backend=derivative_backend,
        series_iterations=series_iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )
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
