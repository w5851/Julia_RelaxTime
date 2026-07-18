"""
    PNJLChiBTaylorDiff

TaylorDiff-backed single-direction PNJL susceptibility helpers.

The route is:
1. solve the primal PNJL gap equation at the target point;
2. build the primal Jacobian `J0 = dF/dx`;
3. represent one flavor-space chemical-potential direction and `x(δ)` as univariate TaylorDiff polynomials;
4. run explicit Taylor-series Newton steps in the gap residual;
5. extract pressure derivatives from the resulting pressure series.
"""
module PNJLChiBTaylorDiff

using ForwardDiff
using LinearAlgebra
using StaticArrays
using TaylorDiff

using ..TaylorDiffForwardDiffCompat
using ..MixedTaylorJets
using ..Conditions: GapParams, gap_core_residual!
using ..Models: AbstractPNJLModel, calculate_mass_vec, create_model, solve_gap, state_vector
using ..Models: cached_nodes, thermal_p_max_inv_fm, model_pressure

export pressure_series_B, chi_B_taylordiff, chi_B_taylordiff_all
export gap_series_parameter_direction, pressure_series_parameter_direction
export pressure_jet_bqs, chi_BQS_mixed_taylorjet
export taylor_constant, taylor_variable, nth_derivative_from_series
export PNJLTaylorBaseContext, build_taylor_base_context

const _MODEL_CACHE = Dict{Symbol, Any}()
const _ALLOWED_LINEAR_SOLVES = (:auto, :refactor_each_order, :factorized_each_order, :factorized_batched)
const _BQS_TO_FLAVOR = @SMatrix [
    1.0 / 3.0   2.0 / 3.0   0.0
    1.0 / 3.0  -1.0 / 3.0   0.0
    1.0 / 3.0  -1.0 / 3.0  -1.0
]

"""
    PNJLTaylorBaseContext

Branch-locked zero-order state and reusable linearization for PNJL Taylor
series derivatives. The context is tied to one model instance and one exact
`(T, mu_vec, xi, p_num, t_num)` point.
"""
struct PNJLTaylorBaseContext{M, F}
    model::M
    T_fm::Float64
    mu_vec::SVector{3, Float64}
    xi::Float64
    p_num::Int
    t_num::Int
    x0::SVector{5, Float64}
    J0::Matrix{Float64}
    factorization::F
    base_state_source::Symbol
    base_state_polished::Bool
    branch_locked::Bool
    base_residual_norm::Float64
    primal_solve_count::Int
    jacobian_factorization_count::Int
    polish_iteration_count::Int
end

@inline function taylor_constant(x::Real, order::Int)
    order >= 1 || throw(ArgumentError("order must be >= 1, got $(order)"))
    return _taylor_constant(Float64(x), Val(order))
end

@inline _taylor_constant(x::Float64, ::Val{P}) where {P} = TaylorDiff.TaylorScalar{P}(x)

@inline function taylor_variable(x::Real, order::Int)
    order >= 1 || throw(ArgumentError("order must be >= 1, got $(order)"))
    return _taylor_variable(Float64(x), Val(order))
end

@inline _taylor_variable(x::Float64, ::Val{P}) where {P} = TaylorDiff.TaylorScalar{P}(x, 1.0)

@inline function nth_derivative_from_series(y::TaylorDiff.TaylorScalar, order::Int)
    order >= 0 || throw(ArgumentError("order must be >= 0, got $(order)"))
    order == 0 && return Float64(TaylorDiff.value(y))
    return Float64(TaylorDiff.extract_derivative(y, Val(order)))
end

@inline function nth_derivative_from_series(y::MixedTaylorJet{D, N, L}, orders::NTuple{D, Int}) where {D, N, L}
    all(n -> n >= 0, orders) || throw(ArgumentError("orders must be non-negative, got $(orders)"))
    sum(orders) <= N || throw(ArgumentError("requested derivative order $(sum(orders)) exceeds jet order $N"))
    return Float64(jet_extract_derivative(y, orders))
end

@inline function _as_svector3(x::SVector{3, <:Real}, name::AbstractString)
    return SVector{3, Float64}(Float64(x[1]), Float64(x[2]), Float64(x[3]))
end

@inline function _as_svector3(x::NTuple{3, <:Real}, name::AbstractString)
    return SVector{3, Float64}(Float64(x[1]), Float64(x[2]), Float64(x[3]))
end

@inline function _as_svector3(x::AbstractVector{<:Real}, name::AbstractString)
    length(x) == 3 || throw(ArgumentError("$name must have length 3, got $(length(x))"))
    return SVector{3, Float64}(Float64(x[1]), Float64(x[2]), Float64(x[3]))
end

@inline function _validate_flavor_direction(direction)
    d = _as_svector3(direction, "direction")
    any(x -> !iszero(x), d) || throw(ArgumentError("direction must not be the zero vector"))
    all(isfinite, d) || throw(ArgumentError("direction must contain only finite values"))
    return d
end

@inline function _validate_parameter_direction(T_direction::Real, mu_direction)
    T_dir = Float64(T_direction)
    isfinite(T_dir) || throw(ArgumentError("T_direction must be finite, got $(T_direction)"))
    mu_dir = _as_svector3(mu_direction, "mu_direction")
    all(isfinite, mu_dir) || throw(ArgumentError("mu_direction must contain only finite values"))
    (!iszero(T_dir) || any(x -> !iszero(x), mu_dir)) ||
        throw(ArgumentError("combined TaylorDiff parameter direction must not be zero"))
    return T_dir, mu_dir
end

@inline function _get_model()
    return get!(_MODEL_CACHE, :PNJL) do
        create_model(:PNJL)
    end
end

@inline function _validate_common_inputs(T_fm::Real, order::Int, p_num::Int, t_num::Int)
    T_fm > 0 || throw(ArgumentError("T_fm must be positive, got $(T_fm)"))
    order >= 1 || throw(ArgumentError("order must be >= 1, got $(order)"))
    p_num > 0 || throw(ArgumentError("p_num must be positive, got $(p_num)"))
    t_num > 0 || throw(ArgumentError("t_num must be positive, got $(t_num)"))
    return nothing
end

@inline function _validate_linear_solve(linear_solve::Symbol)
    linear_solve in _ALLOWED_LINEAR_SOLVES && return linear_solve
    throw(ArgumentError("linear_solve must be one of $(_ALLOWED_LINEAR_SOLVES), got $(linear_solve)"))
end

@inline function _resolve_linear_solve(linear_solve::Symbol, order::Int)
    validated = _validate_linear_solve(linear_solve)
    validated !== :auto && return validated
    return order >= 16 ? :factorized_each_order : :refactor_each_order
end

@inline function _thermal_nodes(model::AbstractPNJLModel, p_num::Int, t_num::Int)
    return cached_nodes(p_num, t_num; p_max_inv_fm=thermal_p_max_inv_fm(model))
end

@inline function _gap_params(
    model::AbstractPNJLModel,
    T_fm,
    xi,
    p_num::Int,
    t_num::Int;
    thermo_quadrature_kwargs::NamedTuple=(;),
)
    return GapParams(
        T_fm,
        _thermal_nodes(model, p_num, t_num),
        xi;
        p_num=p_num,
        t_num=t_num,
        model_kind=:PNJL,
        thermo_quadrature_kwargs...,
    )
end

@inline function _state_svector(x)
    return SVector{5}(x[1], x[2], x[3], x[4], x[5])
end

@inline function _flavor_mu_from_bqs(muB_fm, muQ_fm, muS_fm)
    return SVector{3}(
        muB_fm / 3 + 2 * muQ_fm / 3,
        muB_fm / 3 - muQ_fm / 3,
        muB_fm / 3 - muQ_fm / 3 - muS_fm,
    )
end

@inline function _flavor_direction_from_bqs(direction_bqs)
    d = _as_svector3(direction_bqs, "direction_bqs")
    return _BQS_TO_FLAVOR * d
end

function _gap_residual(
    model::AbstractPNJLModel,
    x_state::SVector{5},
    mu_vec::SVector{3},
    T_fm,
    xi,
    p_num::Int,
    t_num::Int,
    ;
    thermo_quadrature_kwargs::NamedTuple=(;),
)
    params = _gap_params(
        model, T_fm, xi, p_num, t_num; thermo_quadrature_kwargs=thermo_quadrature_kwargs,
    )
    Tout = promote_type(eltype(x_state), eltype(mu_vec), typeof(T_fm), typeof(xi))
    out = Vector{Tout}(undef, 5)
    gap_core_residual!(out, model, x_state, mu_vec, params)
    return SVector{5, Tout}(Tuple(out))
end

function _primal_state(
    model::AbstractPNJLModel,
    T_fm::Real,
    mu_vec0::SVector{3, Float64},
    xi::Real,
    p_num::Int,
    t_num::Int,
    ;
    thermo_quadrature_kwargs::NamedTuple=(;),
)
    st = solve_gap(
        model,
        Float64(T_fm),
        mu_vec0;
        xi=Float64(xi),
        p_num=p_num,
        t_num=t_num,
        thermo_quadrature_kwargs...,
    )
    return SVector{5, Float64}(Tuple(state_vector(st)))
end

@inline function _as_state5(x, name::AbstractString)
    length(x) == 5 || throw(ArgumentError("$name must have length 5, got $(length(x))"))
    state = SVector{5, Float64}(ntuple(i -> Float64(x[i]), Val(5)))
    all(isfinite, state) || throw(ArgumentError("$name must contain only finite values"))
    return state
end

@inline function _masses_from_state(model::AbstractPNJLModel, x_state::SVector{5, Float64})
    phi = SVector{3, Float64}(x_state[1], x_state[2], x_state[3])
    return SVector{3, Float64}(Tuple(calculate_mass_vec(model, phi)))
end

@inline function _primal_residual(
    model::AbstractPNJLModel,
    x_state::SVector{5, Float64},
    T_fm::Real,
    mu_vec::SVector{3, Float64},
    xi::Real,
    p_num::Int,
    t_num::Int,
    ;
    thermo_quadrature_kwargs::NamedTuple=(;),
)
    residual = _gap_residual(
        model,
        x_state,
        mu_vec,
        Float64(T_fm),
        Float64(xi),
        p_num,
        t_num;
        thermo_quadrature_kwargs=thermo_quadrature_kwargs,
    )
    residual_norm = maximum(abs, residual)
    isfinite(residual_norm) || throw(ArgumentError("PNJL TaylorDiff base-state residual is not finite"))
    return residual, Float64(residual_norm)
end

function _assert_same_branch(
    model::AbstractPNJLModel,
    initial_state::SVector{5, Float64},
    candidate_state::SVector{5, Float64};
    expected_masses=nothing,
    state_rtol::Real=0.20,
    mass_rtol::Real=0.10,
)
    initial_masses = expected_masses === nothing ?
        _masses_from_state(model, initial_state) :
        _as_svector3(expected_masses, "base_masses")
    candidate_masses = _masses_from_state(model, candidate_state)

    state_drift = maximum(abs(candidate_state[i] - initial_state[i]) / max(abs(initial_state[i]), 1e-3) for i in 1:5)
    mass_drift = maximum(abs(candidate_masses[i] - initial_masses[i]) / max(abs(initial_masses[i]), 1e-6) for i in 1:3)
    if state_drift > Float64(state_rtol) || mass_drift > Float64(mass_rtol)
        throw(ArgumentError("local PNJL base-state polish would leave the supplied branch: state_relative_drift=$(state_drift), mass_relative_drift=$(mass_drift)"))
    end
    return candidate_masses
end

function _polish_base_state(
    model::AbstractPNJLModel,
    initial_state::SVector{5, Float64},
    T_fm::Float64,
    mu_vec::SVector{3, Float64},
    xi::Float64,
    p_num::Int,
    t_num::Int;
    residual_tol::Float64,
    max_iterations::Int,
    expected_masses=nothing,
)
    x_state = initial_state
    residual, residual_norm = _primal_residual(model, x_state, T_fm, mu_vec, xi, p_num, t_num)
    residual_norm <= residual_tol && return x_state, residual_norm, 0

    for iteration in 1:max_iterations
        J = _jacobian_primal(model, x_state, T_fm, mu_vec, xi, p_num, t_num)
        step = J \ Vector(residual)
        all(isfinite, step) || throw(ArgumentError("local PNJL base-state polish produced a non-finite Newton step"))

        accepted = false
        alpha = 1.0
        for _ in 1:12
            candidate = x_state - alpha * SVector{5, Float64}(Tuple(step))
            _assert_same_branch(model, initial_state, candidate; expected_masses=expected_masses)
            candidate_residual, candidate_norm = _primal_residual(model, candidate, T_fm, mu_vec, xi, p_num, t_num)
            if candidate_norm < residual_norm
                x_state = candidate
                residual = candidate_residual
                residual_norm = candidate_norm
                accepted = true
                break
            end
            alpha *= 0.5
        end
        accepted || throw(ArgumentError("local PNJL base-state polish could not reduce the residual without leaving the supplied branch"))
        residual_norm <= residual_tol && return x_state, residual_norm, iteration
    end

    throw(ArgumentError("local PNJL base-state polish did not reach tolerance: residual_norm=$(residual_norm), tolerance=$(residual_tol), iterations=$(max_iterations)"))
end

function _jacobian_primal(
    model::AbstractPNJLModel,
    x0::SVector{5, Float64},
    T_fm::Real,
    mu_vec0::SVector{3, Float64},
    xi::Real,
    p_num::Int,
    t_num::Int,
    ;
    thermo_quadrature_kwargs::NamedTuple=(;),
)
    residual = x -> Vector(_gap_residual(
        model,
        _state_svector(x),
        mu_vec0,
        Float64(T_fm),
        Float64(xi),
        p_num,
        t_num;
        thermo_quadrature_kwargs=thermo_quadrature_kwargs,
    ))
    J0 = ForwardDiff.jacobian(residual, x0)
    all(isfinite, J0) || throw(ArgumentError("PNJL TaylorDiff gap Jacobian contains non-finite entries"))
    return J0
end

"""
    build_taylor_base_context(T_fm, mu_vec0; kwargs...) -> PNJLTaylorBaseContext

Build one reusable PNJL zero-order state, Jacobian, and factorization. When
`base_state` is supplied, no global gap solve is performed. A residual-only
local Newton polish may refine that same branch; a large state or mass drift is
rejected instead of silently selecting another phase branch.
"""
function build_taylor_base_context(
    T_fm::Real,
    mu_vec0;
    xi::Real=0.0,
    p_num::Int=8,
    t_num::Int=4,
    model=nothing,
    base_state=nothing,
    base_masses=nothing,
    base_state_source::Symbol=base_state === nothing ? :internal_gap_solve : :explicit_base_state,
    base_residual_tol::Real=1e-7,
    allow_local_polish::Bool=true,
    polish_max_iterations::Int=4,
)
    _validate_common_inputs(T_fm, 1, p_num, t_num)
    base_residual_tol > 0 || throw(ArgumentError("base_residual_tol must be positive, got $(base_residual_tol)"))
    polish_max_iterations >= 0 || throw(ArgumentError("polish_max_iterations must be non-negative, got $(polish_max_iterations)"))

    m = model === nothing ? _get_model() : model
    m isa AbstractPNJLModel || throw(ArgumentError("TaylorDiff base context requires an AbstractPNJLModel, got $(typeof(m))"))
    T0 = Float64(T_fm)
    mu0 = _as_svector3(mu_vec0, "mu_vec0")
    xi0 = Float64(xi)
    all(isfinite, mu0) || throw(ArgumentError("mu_vec0 must contain only finite values"))
    isfinite(xi0) || throw(ArgumentError("xi must be finite, got $(xi)"))

    primal_solve_count = base_state === nothing ? 1 : 0
    initial_state = base_state === nothing ?
        _primal_state(m, T0, mu0, xi0, p_num, t_num) :
        _as_state5(base_state, "base_state")

    if base_masses !== nothing
        expected = _as_svector3(base_masses, "base_masses")
        actual = _masses_from_state(m, initial_state)
        all(isapprox(actual[i], expected[i]; rtol=1e-8, atol=1e-10) for i in 1:3) ||
            throw(ArgumentError("base_masses are inconsistent with base_state for the resolved PNJL model"))
    end

    _, initial_residual_norm = _primal_residual(m, initial_state, T0, mu0, xi0, p_num, t_num)
    x0, base_residual_norm, polish_iteration_count = if initial_residual_norm <= Float64(base_residual_tol)
        initial_state, initial_residual_norm, 0
    elseif allow_local_polish && polish_max_iterations > 0
        _polish_base_state(
            m,
            initial_state,
            T0,
            mu0,
            xi0,
            p_num,
            t_num;
            residual_tol=Float64(base_residual_tol),
            max_iterations=polish_max_iterations,
            expected_masses=base_masses,
        )
    else
        throw(ArgumentError("PNJL TaylorDiff base_state residual exceeds tolerance and local polish is disabled: residual_norm=$(initial_residual_norm), tolerance=$(base_residual_tol)"))
    end

    _assert_same_branch(m, initial_state, x0; expected_masses=base_masses)
    J0 = _jacobian_primal(m, x0, T0, mu0, xi0, p_num, t_num)
    Jfac = factorize(J0)
    return PNJLTaylorBaseContext(
        m,
        T0,
        mu0,
        xi0,
        p_num,
        t_num,
        x0,
        J0,
        Jfac,
        base_state_source,
        polish_iteration_count > 0,
        true,
        base_residual_norm,
        primal_solve_count,
        1,
        polish_iteration_count,
    )
end

@inline function _validate_base_context(
    ctx::PNJLTaylorBaseContext,
    model::AbstractPNJLModel,
    T_fm::Real,
    mu_vec::SVector{3, Float64},
    xi::Real,
    p_num::Int,
    t_num::Int,
)
    ctx.model === model || throw(ArgumentError("PNJL TaylorDiff base context model does not match the requested model instance"))
    ctx.T_fm == Float64(T_fm) || throw(ArgumentError("PNJL TaylorDiff base context T_fm does not match the requested point"))
    ctx.mu_vec == mu_vec || throw(ArgumentError("PNJL TaylorDiff base context mu_vec does not match the requested point"))
    ctx.xi == Float64(xi) || throw(ArgumentError("PNJL TaylorDiff base context xi does not match the requested point"))
    ctx.p_num == p_num || throw(ArgumentError("PNJL TaylorDiff base context p_num does not match the requested numerics"))
    ctx.t_num == t_num || throw(ArgumentError("PNJL TaylorDiff base context t_num does not match the requested numerics"))
    return ctx
end

@inline function _coefficient(x::TaylorDiff.TaylorScalar, n::Int)
    n == 0 && return Float64(TaylorDiff.value(x))
    return Float64(TaylorDiff.partials(x)[n])
end

@inline _coefficient(x::MixedTaylorJet, pos::Int) = Float64(jet_coefficient(x, pos))
@inline _coefficient(x::Real, n::Int) = n == 0 ? Float64(x) : 0.0

@inline function _from_coefficients(coeffs, ::Val{P}) where {P}
    return TaylorDiff.TaylorScalar(Float64(coeffs[1]), ntuple(i -> Float64(coeffs[i + 1]), Val(P)))
end

@inline function _series_from_coeff_matrix(coeffs::AbstractMatrix{<:Real}, ::Val{P}) where {P}
    return SVector{5}(ntuple(i -> _from_coefficients(@view(coeffs[i, :]), Val(P)), 5))
end

@inline function _series_from_base_direction(
    base_mu_vec::SVector{3, Float64},
    direction::SVector{3, Float64},
    delta::TaylorDiff.TaylorScalar,
    ::Val{P},
) where {P}
    return SVector{3}(
        base_mu_vec[1] + direction[1] * delta,
        base_mu_vec[2] + direction[2] * delta,
        base_mu_vec[3] + direction[3] * delta,
    )
end

@inline function _series_from_parameter_direction(
    T_fm::Real,
    base_mu_vec::SVector{3, Float64},
    T_direction::Float64,
    mu_direction::SVector{3, Float64},
    delta::TaylorDiff.TaylorScalar,
    ::Val{P},
) where {P}
    T_series = taylor_constant(T_fm, P) + T_direction * delta
    mu_vec = _series_from_base_direction(base_mu_vec, mu_direction, delta, Val(P))
    return T_series, mu_vec
end

function _rhs_coefficients(rhs::SVector{5}, ::Val{P}) where {P}
    rhs_coeffs = Matrix{Float64}(undef, 5, P + 1)
    @inbounds for n in 0:P
        col = n + 1
        for i in 1:5
            rhs_coeffs[i, col] = _coefficient(rhs[i], n)
        end
    end
    return rhs_coeffs
end

function _solve_series_linear_refactor_each_order(J0, rhs::SVector{5}, ::Val{P}) where {P}
    coeffs = Matrix{Float64}(undef, 5, P + 1)
    rhs_n = Vector{Float64}(undef, 5)
    @inbounds for n in 0:P
        for i in 1:5
            rhs_n[i] = _coefficient(rhs[i], n)
        end
        sol_n = J0 \ rhs_n
        for i in 1:5
            coeffs[i, n + 1] = sol_n[i]
        end
    end
    return _series_from_coeff_matrix(coeffs, Val(P))
end

function _solve_series_linear_factorized_each_order(Jfac, rhs::SVector{5}, ::Val{P}) where {P}
    coeffs = Matrix{Float64}(undef, 5, P + 1)
    rhs_n = Vector{Float64}(undef, 5)
    @inbounds for n in 0:P
        for i in 1:5
            rhs_n[i] = _coefficient(rhs[i], n)
        end
        sol_n = Jfac \ rhs_n
        for i in 1:5
            coeffs[i, n + 1] = sol_n[i]
        end
    end
    return _series_from_coeff_matrix(coeffs, Val(P))
end

function _solve_series_linear_factorized_batched(Jfac, rhs::SVector{5}, ::Val{P}) where {P}
    rhs_coeffs = _rhs_coefficients(rhs, Val(P))
    sol_coeffs = Jfac \ rhs_coeffs
    return _series_from_coeff_matrix(sol_coeffs, Val(P))
end

@inline function _solve_series_linear(linear_operator, rhs::SVector{5}, ::Val{P}, linear_solve::Symbol) where {P}
    if linear_solve === :factorized_batched
        return _solve_series_linear_factorized_batched(linear_operator, rhs, Val(P))
    elseif linear_solve === :factorized_each_order
        return _solve_series_linear_factorized_each_order(linear_operator, rhs, Val(P))
    end
    return _solve_series_linear_refactor_each_order(linear_operator, rhs, Val(P))
end

function _series_residual_norm(residual::SVector{5}, order::Int)
    max_norm = 0.0
    @inbounds for i in 1:5
        for n in 0:order
            max_norm = max(max_norm, abs(_coefficient(residual[i], n)))
        end
    end
    return max_norm
end

function _series_residual_scaled_norm(residual::SVector{5}, x_series::SVector{5}, order::Int)
    max_norm = 0.0
    @inbounds for n in 0:order
        x_scale = maximum(abs(_coefficient(x_series[i], n)) for i in 1:5)
        scale = max(1.0, x_scale)
        for i in 1:5
            max_norm = max(max_norm, abs(_coefficient(residual[i], n)) / scale)
        end
    end
    return max_norm
end

function _assert_series_residual!(
    residual::SVector{5},
    x_series::SVector{5},
    order::Int,
    residual_norm::Float64,
    series_residual_tol::Real,
)
    isfinite(residual_norm) || throw(ArgumentError("PNJL TaylorDiff gap series residual is not finite"))
    primal_residual_norm = maximum(abs(_coefficient(residual[i], 0)) for i in 1:5)
    scaled_residual_norm = _series_residual_scaled_norm(residual, x_series, order)
    tol = Float64(series_residual_tol)
    (primal_residual_norm <= tol && scaled_residual_norm <= tol) && return nothing
    throw(ArgumentError("PNJL TaylorDiff gap series did not converge: residual_norm=$(residual_norm), scaled_residual_norm=$(scaled_residual_norm), primal_residual_norm=$(primal_residual_norm), tolerance=$(series_residual_tol), order=$(order)"))
end

@inline function _solve_numeric_linear(linear_operator, rhs::AbstractVector{<:Real})
    return linear_operator \ rhs
end

function _solve_gap_series_parameter_direction(
    model::AbstractPNJLModel,
    T_fm::Real,
    mu_vec0,
    T_direction::Real,
    mu_direction;
    order::Int,
    xi::Real,
    p_num::Int,
    t_num::Int,
    iterations::Int=0,
    linear_solve::Symbol=:auto,
    series_residual_tol::Real=1e-7,
    base_context::Union{Nothing, PNJLTaylorBaseContext}=nothing,
    thermo_quadrature_kwargs::NamedTuple=(;),
)
    _validate_common_inputs(T_fm, order, p_num, t_num)
    iterations >= 0 || throw(ArgumentError("series_iterations must be non-negative, got $(iterations)"))
    resolved_linear_solve = if base_context !== nothing && linear_solve === :auto
        :factorized_each_order
    else
        _resolve_linear_solve(linear_solve, order)
    end

    base_mu_vec = _as_svector3(mu_vec0, "mu_vec0")
    T_dir, direction_vec = _validate_parameter_direction(T_direction, mu_direction)

    context = if base_context === nothing
        nothing
    else
        isempty(thermo_quadrature_kwargs) || throw(ArgumentError(
            "custom thermo quadrature is not supported with a prebuilt PNJLTaylorBaseContext",
        ))
        _validate_base_context(base_context, model, T_fm, base_mu_vec, xi, p_num, t_num)
    end
    x0 = context === nothing ?
        _primal_state(model, T_fm, base_mu_vec, xi, p_num, t_num; thermo_quadrature_kwargs=thermo_quadrature_kwargs) :
        context.x0
    J0 = context === nothing ?
        _jacobian_primal(model, x0, T_fm, base_mu_vec, xi, p_num, t_num; thermo_quadrature_kwargs=thermo_quadrature_kwargs) :
        context.J0
    linear_operator = if resolved_linear_solve === :refactor_each_order
        J0
    elseif context === nothing
        factorize(J0)
    else
        context.factorization
    end

    delta = taylor_variable(0.0, order)
    T_series, mu_vec = _series_from_parameter_direction(
        T_fm,
        base_mu_vec,
        T_dir,
        direction_vec,
        delta,
        Val(order),
    )

    coeffs = zeros(Float64, 5, order + 1)
    @inbounds for i in 1:5
        coeffs[i, 1] = x0[i]
    end

    rhs_n = Vector{Float64}(undef, 5)
    for n in 1:order
        x_series_n = _series_from_coeff_matrix(coeffs, Val(order))
        residual_n = _gap_residual(model, x_series_n, mu_vec, T_series, Float64(xi), p_num, t_num;
            thermo_quadrature_kwargs=thermo_quadrature_kwargs)
        @inbounds for i in 1:5
            rhs_n[i] = _coefficient(residual_n[i], n)
        end
        sol_n = _solve_numeric_linear(linear_operator, rhs_n)
        @inbounds for i in 1:5
            coeffs[i, n + 1] = -sol_n[i]
        end
    end

    x_series = _series_from_coeff_matrix(coeffs, Val(order))

    for _ in 1:iterations
        residual = _gap_residual(model, x_series, mu_vec, T_series, Float64(xi), p_num, t_num;
            thermo_quadrature_kwargs=thermo_quadrature_kwargs)
        step = _solve_series_linear(linear_operator, residual, Val(order), resolved_linear_solve)
        x_series = x_series - step
    end

    residual = _gap_residual(model, x_series, mu_vec, T_series, Float64(xi), p_num, t_num;
        thermo_quadrature_kwargs=thermo_quadrature_kwargs)
    residual_norm = _series_residual_norm(residual, order)
    _assert_series_residual!(residual, x_series, order, residual_norm, series_residual_tol)

    return (
        x_state=x_series,
        mu_vec=mu_vec,
        T_series=T_series,
        residual=residual,
        residual_norm=residual_norm,
        iterations=iterations,
        linear_solve=resolved_linear_solve,
        base_context=context,
    )
end

function _solve_gap_series_direction(
    model::AbstractPNJLModel,
    T_fm::Real,
    mu_vec0,
    direction;
    kwargs...,
)
    direction_vec = _validate_flavor_direction(direction)
    return _solve_gap_series_parameter_direction(
        model,
        T_fm,
        mu_vec0,
        0.0,
        direction_vec;
        kwargs...,
    )
end

function _solve_gap_series(
    model::AbstractPNJLModel,
    T_fm::Real,
    muB_fm::Real;
    muQ_fm::Real=0.0,
    muS_fm::Real=0.0,
    kwargs...,
)
    mu_vec0 = _flavor_mu_from_bqs(Float64(muB_fm), Float64(muQ_fm), Float64(muS_fm))
    direction = _flavor_direction_from_bqs((1.0, 0.0, 0.0))
    return _solve_gap_series_direction(model, T_fm, mu_vec0, direction; kwargs...)
end

"""
    pressure_series_direction(T_fm, mu_vec0, direction; order, kwargs...) -> TaylorScalar

Evaluate PNJL pressure along an arbitrary flavor-space direction `direction`
starting from the base flavor chemical potential vector `mu_vec0`.
"""
function pressure_series_direction(
    T_fm::Real,
    mu_vec0,
    direction;
    order::Int,
    xi::Real=0.0,
    p_num::Int=8,
    t_num::Int=4,
    model=nothing,
    series_iterations::Union{Nothing, Int}=nothing,
    linear_solve::Symbol=:auto,
    series_residual_tol::Real=1e-7,
)
    _validate_common_inputs(T_fm, order, p_num, t_num)
    m = model === nothing ? _get_model() : model
    m isa AbstractPNJLModel || throw(ArgumentError("TaylorDiff chi_B backend requires an AbstractPNJLModel, got $(typeof(m))"))
    iterations = series_iterations === nothing ? 0 : Int(series_iterations)
    result = _solve_gap_series_direction(
        m,
        T_fm,
        mu_vec0,
        direction;
        order=order,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        iterations=iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )
    pressure = model_pressure(m, result.x_state, result.mu_vec, result.T_series; p_num=p_num, t_num=t_num, xi=xi)
    isfinite(TaylorDiff.value(pressure)) || throw(ArgumentError("PNJL TaylorDiff pressure series value is not finite"))
    all(isfinite, TaylorDiff.partials(pressure)) || throw(ArgumentError("PNJL TaylorDiff pressure series contains non-finite coefficients"))
    return pressure
end

"""
    gap_series_parameter_direction(T_fm, mu_vec0, T_direction, mu_direction; order, kwargs...)

Solve the PNJL gap equation as a univariate Taylor series along the combined
parameter ray

`T(δ) = T_fm + T_direction * δ`,
`mu_vec(δ) = mu_vec0 + mu_direction * δ`.

The returned named tuple contains Taylor-series `x_state`, `T_series`, and
`mu_vec` fields. This is the shared low-level route for pressure, mass, thermo,
and conserved-charge derivatives.
"""
function gap_series_parameter_direction(
    T_fm::Real,
    mu_vec0,
    T_direction::Real,
    mu_direction;
    order::Int,
    xi::Real=0.0,
    p_num::Int=8,
    t_num::Int=4,
    model=nothing,
    series_iterations::Union{Nothing, Int}=nothing,
    linear_solve::Symbol=:auto,
    series_residual_tol::Real=1e-7,
    base_context::Union{Nothing, PNJLTaylorBaseContext}=nothing,
    thermo_quadrature_kwargs::NamedTuple=(;),
)
    _validate_common_inputs(T_fm, order, p_num, t_num)
    m = model === nothing ? _get_model() : model
    m isa AbstractPNJLModel || throw(ArgumentError("TaylorDiff gap-series backend requires an AbstractPNJLModel, got $(typeof(m))"))
    iterations = series_iterations === nothing ? 0 : Int(series_iterations)
    return _solve_gap_series_parameter_direction(
        m,
        T_fm,
        mu_vec0,
        T_direction,
        mu_direction;
        order=order,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        iterations=iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
        base_context=base_context,
        thermo_quadrature_kwargs=thermo_quadrature_kwargs,
    )
end

"""
    pressure_series_parameter_direction(T_fm, mu_vec0, T_direction, mu_direction; order, kwargs...)

Evaluate PNJL pressure on the same combined parameter ray used by
`gap_series_parameter_direction`.
"""
function pressure_series_parameter_direction(
    T_fm::Real,
    mu_vec0,
    T_direction::Real,
    mu_direction;
    order::Int,
    xi::Real=0.0,
    p_num::Int=8,
    t_num::Int=4,
    model=nothing,
    series_iterations::Union{Nothing, Int}=nothing,
    linear_solve::Symbol=:auto,
    series_residual_tol::Real=1e-7,
    base_context::Union{Nothing, PNJLTaylorBaseContext}=nothing,
)
    result = gap_series_parameter_direction(
        T_fm,
        mu_vec0,
        T_direction,
        mu_direction;
        order=order,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        model=model,
        series_iterations=series_iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
        base_context=base_context,
    )
    m = model === nothing ? _get_model() : model
    pressure = model_pressure(m, result.x_state, result.mu_vec, result.T_series; p_num=p_num, t_num=t_num, xi=xi)
    isfinite(TaylorDiff.value(pressure)) || throw(ArgumentError("PNJL TaylorDiff pressure series value is not finite"))
    all(isfinite, TaylorDiff.partials(pressure)) || throw(ArgumentError("PNJL TaylorDiff pressure series contains non-finite coefficients"))
    return pressure
end

@inline function _active_bqs_axes(orders::NTuple{3, Int})
    all(n -> n >= 0, orders) || throw(ArgumentError("orders must be non-negative, got $(orders)"))
    axes = Int[]
    for i in 1:3
        orders[i] > 0 && push!(axes, i)
    end
    return Tuple(axes)
end

@inline function _jet_from_coeffs(coeffs::AbstractMatrix{<:Real}, row::Int, ::Type{MixedTaylorJet{D, N, L}}) where {D, N, L}
    return MixedTaylorJet{D, N, L}(ntuple(j -> Float64(coeffs[row, j]), Val(L)))
end

@inline function _jet_series_from_coeff_matrix(coeffs::AbstractMatrix{<:Real}, ::Type{MixedTaylorJet{D, N, L}}) where {D, N, L}
    return SVector{5}(ntuple(i -> _jet_from_coeffs(coeffs, i, MixedTaylorJet{D, N, L}), 5))
end

function _solve_jet_series_linear(linear_operator, rhs::SVector{5}, ::Type{MixedTaylorJet{D, N, L}}) where {D, N, L}
    rhs_coeffs = Matrix{Float64}(undef, 5, L)
    @inbounds for pos in 1:L
        for i in 1:5
            rhs_coeffs[i, pos] = _coefficient(rhs[i], pos)
        end
    end
    sol_coeffs = linear_operator \ rhs_coeffs
    return _jet_series_from_coeff_matrix(sol_coeffs, MixedTaylorJet{D, N, L})
end

function _jet_residual_norm(residual::SVector{5}, ::Type{MixedTaylorJet{D, N, L}}) where {D, N, L}
    max_norm = 0.0
    @inbounds for i in 1:5
        for pos in 1:L
            max_norm = max(max_norm, abs(_coefficient(residual[i], pos)))
        end
    end
    return max_norm
end

function _jet_residual_scaled_norm(residual::SVector{5}, x_series::SVector{5}, ::Type{MixedTaylorJet{D, N, L}}) where {D, N, L}
    max_norm = 0.0
    @inbounds for pos in 1:L
        x_scale = maximum(abs(_coefficient(x_series[i], pos)) for i in 1:5)
        scale = max(1.0, x_scale)
        for i in 1:5
            max_norm = max(max_norm, abs(_coefficient(residual[i], pos)) / scale)
        end
    end
    return max_norm
end

function _assert_jet_residual!(
    residual::SVector{5},
    x_series::SVector{5},
    ::Type{MixedTaylorJet{D, N, L}},
    residual_norm::Float64,
    series_residual_tol::Real,
) where {D, N, L}
    isfinite(residual_norm) || throw(ArgumentError("PNJL mixed Taylor jet gap residual is not finite"))
    primal_residual_norm = maximum(abs(_coefficient(residual[i], 1)) for i in 1:5)
    scaled_residual_norm = _jet_residual_scaled_norm(residual, x_series, MixedTaylorJet{D, N, L})
    tol = Float64(series_residual_tol)
    (primal_residual_norm <= tol && scaled_residual_norm <= tol) && return nothing
    throw(ArgumentError("PNJL mixed Taylor jet gap series did not converge: residual_norm=$(residual_norm), scaled_residual_norm=$(scaled_residual_norm), primal_residual_norm=$(primal_residual_norm), tolerance=$(series_residual_tol), order=$(N), dim=$(D)"))
end

function _mixed_mu_vec(
    base_mu_vec::SVector{3, Float64},
    active_axes::NTuple{D, Int},
    deltas::NTuple{D, <:MixedTaylorJet},
    ::Type{MixedTaylorJet{D, N, L}},
) where {D, N, L}
    return SVector{3}(ntuple(flavor -> begin
        value = convert(MixedTaylorJet{D, N, L}, base_mu_vec[flavor])
        for j in 1:D
            value = value + _BQS_TO_FLAVOR[flavor, active_axes[j]] * deltas[j]
        end
        value
    end, Val(3)))
end

function _solve_gap_mixed_jet_bqs_impl(
    ::Val{D},
    ::Val{N},
    model::AbstractPNJLModel,
    T_fm::Real,
    muB_fm::Real,
    muQ_fm::Real,
    muS_fm::Real,
    active_axes::NTuple{D, Int};
    xi::Real,
    p_num::Int,
    t_num::Int,
    iterations::Int=0,
    linear_solve::Symbol=:auto,
    series_residual_tol::Real=1e-7,
) where {D, N}
    _validate_common_inputs(T_fm, N, p_num, t_num)
    D >= 2 || throw(ArgumentError("mixed Taylor jet backend expects at least two active BQS axes, got D=$D"))
    iterations >= 0 || throw(ArgumentError("series_iterations must be non-negative, got $(iterations)"))
    resolved_linear_solve = _resolve_linear_solve(linear_solve, N)

    base_mu_vec = _flavor_mu_from_bqs(Float64(muB_fm), Float64(muQ_fm), Float64(muS_fm))
    x0 = _primal_state(model, T_fm, base_mu_vec, xi, p_num, t_num)
    J0 = _jacobian_primal(model, x0, T_fm, base_mu_vec, xi, p_num, t_num)
    linear_operator = resolved_linear_solve === :refactor_each_order ? J0 : factorize(J0)

    deltas = ntuple(i -> jet_variable(0.0, i, D, N), Val(D))
    JetT = typeof(deltas[1])
    T_series = convert(JetT, Float64(T_fm))
    mu_vec = _mixed_mu_vec(base_mu_vec, active_axes, deltas, JetT)

    levels = jet_degree_positions(JetT)
    L = length(jet_basis(JetT))
    coeffs = zeros(Float64, 5, L)
    @inbounds for i in 1:5
        coeffs[i, 1] = x0[i]
    end

    for degree in 1:N
        x_series_n = _jet_series_from_coeff_matrix(coeffs, JetT)
        residual_n = _gap_residual(model, x_series_n, mu_vec, T_series, Float64(xi), p_num, t_num)
        positions = levels[degree + 1]
        rhs_mat = Matrix{Float64}(undef, 5, length(positions))
        @inbounds for (col, pos) in enumerate(positions)
            for i in 1:5
                rhs_mat[i, col] = _coefficient(residual_n[i], pos)
            end
        end
        sol_mat = linear_operator \ rhs_mat
        @inbounds for (col, pos) in enumerate(positions)
            for i in 1:5
                coeffs[i, pos] = -sol_mat[i, col]
            end
        end
    end

    x_series = _jet_series_from_coeff_matrix(coeffs, JetT)

    for _ in 1:iterations
        residual = _gap_residual(model, x_series, mu_vec, T_series, Float64(xi), p_num, t_num)
        step = _solve_jet_series_linear(linear_operator, residual, JetT)
        x_series = x_series - step
    end

    residual = _gap_residual(model, x_series, mu_vec, T_series, Float64(xi), p_num, t_num)
    residual_norm = _jet_residual_norm(residual, JetT)
    _assert_jet_residual!(residual, x_series, JetT, residual_norm, series_residual_tol)

    return (
        x_state=x_series,
        mu_vec=mu_vec,
        T_series=T_series,
        residual=residual,
        residual_norm=residual_norm,
        iterations=iterations,
        linear_solve=resolved_linear_solve,
        active_axes=active_axes,
    )
end

@inline function _solve_gap_mixed_jet_bqs(
    model::AbstractPNJLModel,
    T_fm::Real,
    muB_fm::Real,
    muQ_fm::Real,
    muS_fm::Real,
    active_axes::NTuple{2, Int};
    order::Int,
    kwargs...,
)
    return _solve_gap_mixed_jet_bqs_impl(Val(2), Val(order), model, T_fm, muB_fm, muQ_fm, muS_fm, active_axes; kwargs...)
end

@inline function _solve_gap_mixed_jet_bqs(
    model::AbstractPNJLModel,
    T_fm::Real,
    muB_fm::Real,
    muQ_fm::Real,
    muS_fm::Real,
    active_axes::NTuple{3, Int};
    order::Int,
    kwargs...,
)
    return _solve_gap_mixed_jet_bqs_impl(Val(3), Val(order), model, T_fm, muB_fm, muQ_fm, muS_fm, active_axes; kwargs...)
end

"""
    pressure_jet_bqs(T_fm, muB_fm, muQ_fm, muS_fm; orders, kwargs...) -> MixedTaylorJet

Evaluate PNJL pressure as a multivariate Taylor jet in the active B/Q/S axes
specified by `orders`. This is intended for mixed derivatives; single-axis
public routes keep using the univariate TaylorDiff fast path.
"""
function pressure_jet_bqs(
    T_fm::Real,
    muB_fm::Real,
    muQ_fm::Real,
    muS_fm::Real;
    orders::NTuple{3, Int},
    xi::Real=0.0,
    p_num::Int=8,
    t_num::Int=4,
    model=nothing,
    series_iterations::Union{Nothing, Int}=nothing,
    linear_solve::Symbol=:auto,
    series_residual_tol::Real=1e-7,
)
    total_order = sum(orders)
    total_order >= 1 || throw(ArgumentError("orders must have positive total order, got $(orders)"))
    active_axes = _active_bqs_axes(orders)
    length(active_axes) >= 2 || throw(ArgumentError("pressure_jet_bqs is for mixed derivatives; got single-axis orders=$(orders)"))
    m = model === nothing ? _get_model() : model
    m isa AbstractPNJLModel || throw(ArgumentError("mixed Taylor jet backend requires an AbstractPNJLModel, got $(typeof(m))"))
    iterations = series_iterations === nothing ? 0 : Int(series_iterations)
    result = _solve_gap_mixed_jet_bqs(
        m,
        T_fm,
        muB_fm,
        muQ_fm,
        muS_fm,
        active_axes;
        order=total_order,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        iterations=iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )
    pressure = model_pressure(m, result.x_state, result.mu_vec, result.T_series; p_num=p_num, t_num=t_num, xi=xi)
    isfinite(pressure) || throw(ArgumentError("PNJL mixed Taylor jet pressure contains non-finite coefficients"))
    return pressure
end

function chi_BQS_mixed_taylorjet(
    T_fm::Real,
    muB_fm::Real,
    muQ_fm::Real,
    muS_fm::Real;
    orders::NTuple{3, Int},
    xi::Real=0.0,
    p_num::Int=8,
    t_num::Int=4,
    model=nothing,
    series_iterations::Union{Nothing, Int}=nothing,
    linear_solve::Symbol=:auto,
    series_residual_tol::Real=1e-7,
)
    total_order = sum(orders)
    total_order >= 1 || throw(ArgumentError("orders must have positive total order, got $(orders)"))
    active_axes = _active_bqs_axes(orders)
    if length(active_axes) == 1
        axis = active_axes[1]
        base_mu_vec = _flavor_mu_from_bqs(Float64(muB_fm), Float64(muQ_fm), Float64(muS_fm))
        direction = _flavor_direction_from_bqs(ntuple(i -> i == axis ? 1.0 : 0.0, 3))
        return chi_direction_taylordiff(
            T_fm,
            base_mu_vec,
            direction;
            order=total_order,
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            model=model,
            series_iterations=series_iterations,
            linear_solve=linear_solve,
            series_residual_tol=series_residual_tol,
        )
    end

    pressure = pressure_jet_bqs(
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
    active_orders = ntuple(i -> orders[active_axes[i]], length(active_axes))
    return nth_derivative_from_series(pressure, active_orders) * Float64(T_fm)^(total_order - 4)
end

"""
    pressure_series_B(T_fm, muB_fm; order, kwargs...) -> TaylorScalar

Evaluate PNJL pressure along the baryon direction with `mu_B` represented as a
TaylorDiff univariate series. `muQ_fm` and `muS_fm` are held constant when
provided.
"""
function pressure_series_B(
    T_fm::Real,
    muB_fm::Real;
    muQ_fm::Real=0.0,
    muS_fm::Real=0.0,
    order::Int,
    xi::Real=0.0,
    p_num::Int=8,
    t_num::Int=4,
    model=nothing,
    series_iterations::Union{Nothing, Int}=nothing,
    linear_solve::Symbol=:auto,
    series_residual_tol::Real=1e-7,
)
    base_mu_vec = _flavor_mu_from_bqs(Float64(muB_fm), Float64(muQ_fm), Float64(muS_fm))
    direction = _flavor_direction_from_bqs((1.0, 0.0, 0.0))
    return pressure_series_direction(
        T_fm,
        base_mu_vec,
        direction;
        order=order,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        model=model,
        series_iterations=series_iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )
end

"""
    chi_direction_taylordiff(T_fm, mu_vec0, direction; order, kwargs...) -> Float64

TaylorDiff-based susceptibility along an arbitrary flavor-space direction.
"""
function chi_direction_taylordiff(
    T_fm::Real,
    mu_vec0,
    direction;
    order::Int=1,
    xi::Real=0.0,
    p_num::Int=8,
    t_num::Int=4,
    model=nothing,
    series_iterations::Union{Nothing, Int}=nothing,
    linear_solve::Symbol=:auto,
    series_residual_tol::Real=1e-7,
)
    pressure = pressure_series_direction(
        T_fm,
        mu_vec0,
        direction;
        order=order,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        model=model,
        series_iterations=series_iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )
    return nth_derivative_from_series(pressure, order) * Float64(T_fm)^(order - 4)
end

function chi_direction_taylordiff_all(
    T_fm::Real,
    mu_vec0,
    direction;
    max_order::Int=4,
    xi::Real=0.0,
    p_num::Int=8,
    t_num::Int=4,
    model=nothing,
    series_iterations::Union{Nothing, Int}=nothing,
    linear_solve::Symbol=:auto,
    series_residual_tol::Real=1e-7,
)
    pressure = pressure_series_direction(
        T_fm,
        mu_vec0,
        direction;
        order=max_order,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        model=model,
        series_iterations=series_iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )
    return [
        nth_derivative_from_series(pressure, n) * Float64(T_fm)^(n - 4)
        for n in 1:max_order
    ]
end

"""
    chi_B_taylordiff(T_fm, muB_fm; order, kwargs...) -> Float64

TaylorDiff-based baryon susceptibility:

`chi_n^B = T^(n - 4) d^n P / d mu_B^n`.
"""
function chi_B_taylordiff(
    T_fm::Real,
    muB_fm::Real;
    muQ_fm::Real=0.0,
    muS_fm::Real=0.0,
    order::Int=1,
    xi::Real=0.0,
    p_num::Int=8,
    t_num::Int=4,
    model=nothing,
    series_iterations::Union{Nothing, Int}=nothing,
    linear_solve::Symbol=:auto,
    series_residual_tol::Real=1e-7,
)
    base_mu_vec = _flavor_mu_from_bqs(Float64(muB_fm), Float64(muQ_fm), Float64(muS_fm))
    direction = _flavor_direction_from_bqs((1.0, 0.0, 0.0))
    pressure = pressure_series_direction(
        T_fm,
        base_mu_vec,
        direction;
        order=order,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        model=model,
        series_iterations=series_iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )
    return nth_derivative_from_series(pressure, order) * Float64(T_fm)^(order - 4)
end

function chi_B_taylordiff_all(
    T_fm::Real,
    muB_fm::Real;
    muQ_fm::Real=0.0,
    muS_fm::Real=0.0,
    max_order::Int=4,
    xi::Real=0.0,
    p_num::Int=8,
    t_num::Int=4,
    model=nothing,
    series_iterations::Union{Nothing, Int}=nothing,
    linear_solve::Symbol=:auto,
    series_residual_tol::Real=1e-7,
)
    base_mu_vec = _flavor_mu_from_bqs(Float64(muB_fm), Float64(muQ_fm), Float64(muS_fm))
    direction = _flavor_direction_from_bqs((1.0, 0.0, 0.0))
    pressure = pressure_series_direction(
        T_fm,
        base_mu_vec,
        direction;
        order=max_order,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        model=model,
        series_iterations=series_iterations,
        linear_solve=linear_solve,
        series_residual_tol=series_residual_tol,
    )
    return [
        nth_derivative_from_series(pressure, n) * Float64(T_fm)^(n - 4)
        for n in 1:max_order
    ]
end

end # module PNJLChiBTaylorDiff
