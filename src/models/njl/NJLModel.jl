"""NJLModel

将 `NJLCore` 封装为新架构下的具体模型类型，并实现最小接口。

当前阶段只实现：
- 参数持有与构造
- 质量方程 `calculate_mass_vec`
- 手征项 `calculate_chiral`
- Polyakov 势 `polyakov_potential`（NJL=0）
"""

include(joinpath(@__DIR__, "NJLCore.jl"))

using Base.MathConstants: π
using StaticArrays

const _GAUSSLEGENDRE_PATH = normpath(joinpath(@__DIR__, "..", "..", "integration", "GaussLegendre.jl"))
if !isdefined(@__MODULE__, :GaussLegendre)
    include(_GAUSSLEGENDRE_PATH)
end
using .GaussLegendre:
    gauleg,
    DEFAULT_COSΘ_HALF_NODES,
    DEFAULT_COSΘ_HALF_WEIGHTS,
    DEFAULT_MOMENTUM_NODES,
    DEFAULT_MOMENTUM_WEIGHTS

const DEFAULT_THETA_COUNT = length(DEFAULT_COSΘ_HALF_NODES)
const DEFAULT_MOMENTUM_COUNT = length(DEFAULT_MOMENTUM_NODES)

const _NODE_CACHE = Dict{Tuple{Int, Int}, NTuple{3, Matrix{Float64}}}()

@inline function _theta_nodes(t_num::Int)
    if t_num == DEFAULT_THETA_COUNT
        return DEFAULT_COSΘ_HALF_NODES, DEFAULT_COSΘ_HALF_WEIGHTS .* 2.0
    end
    nodes, weights = gauleg(0.0, 1.0, t_num)
    return nodes, weights .* 2.0
end

@inline function _momentum_nodes(p_num::Int)
    if p_num == DEFAULT_MOMENTUM_COUNT
        return DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
    end
    return gauleg(0.0, 10.0, p_num)
end

@inline function _cached_nodes(p_num::Int, t_num::Int)
    key = (p_num, t_num)
    return get!(_NODE_CACHE, key) do
        momentum_nodes, momentum_weights = _momentum_nodes(p_num)
        cosθ_nodes, cosθ_weights = _theta_nodes(t_num)

        p_mesh = repeat(momentum_nodes, 1, t_num)
        cosθ_mesh = repeat(cosθ_nodes', p_num, 1)
        weight_mesh = momentum_weights * cosθ_weights'
        coefficients = weight_mesh .* p_mesh.^2 ./ (2 * π)^2
        return (p_mesh, cosθ_mesh, coefficients)
    end
end

@inline function _energy_anisotropic(mass, p, xi, cosθ)
    return sqrt(p^2 + mass^2 + xi * (p * cosθ)^2)
end

@inline function _log1pexp(x)
    # For ForwardDiff.Dual, comparisons like `x > 0` are not defined.
    # In NJL usage here, the argument is typically negative (-(E±μ)/T),
    # so the simple form is stable enough and AD-friendly.
    return log1p(exp(x))
end

@inline function _vacuum_integral(Λ, mass)
    mass_abs = abs(mass)
    epsilon = one(mass_abs) * 1e-12
    mass_safe = mass_abs + epsilon
    sqrt_term = sqrt(Λ^2 + mass_safe^2)
    poly_part = Λ * sqrt_term * (2 * Λ^2 + mass_safe^2)
    log_term = mass_safe^4 * log((Λ + sqrt_term) / mass_safe)
    return (poly_part - log_term) / (16 * π^2)
end

"""三味 NJL 模型（仅包含 NJLCore 参数）。"""
struct NJLModel <: AbstractNJLModel
    params::NJLCore.NJLParams
end

"""从 config profile 构造 NJLModel。"""
function NJLModel(; profile::String=get(ENV, "NJL_PARAM_PROFILE", "default"))
    return NJLModel(NJLCore.njl_params(profile=profile))
end

# -------------------------
# Interface implementations
# -------------------------

function calculate_mass_vec(model::NJLModel, φ; kwargs...)
    return NJLCore.calculate_mass_vec(model.params, φ)
end

function calculate_chiral(model::NJLModel, φ; kwargs...)
    return NJLCore.chiral_potential(model.params, φ)
end

function polyakov_potential(::NJLModel, Φ, Φbar, T; kwargs...)
    return zero(promote_type(typeof(Φ), typeof(Φbar), typeof(T)))
end

function vacuum_contribution(model::NJLModel, masses; kwargs...)
    p = model.params
    if p.N_flavor != 3
        throw(ArgumentError("vacuum_contribution currently supports N_flavor=3, got $(p.N_flavor)"))
    end
    Tm = eltype(masses)
    Λ = p.Lambda_inv_fm * one(Tm)
    total = zero(Tm)
    @inbounds for i in 1:3
        total += _vacuum_integral(Λ, masses[i])
    end
    return -2 * convert(Tm, p.N_color) * total
end

function thermal_contribution(
    model::NJLModel,
    masses,
    Φ,
    Φbar,
    mu_vec,
    T;
    p_num::Int=DEFAULT_MOMENTUM_COUNT,
    t_num::Int=DEFAULT_THETA_COUNT,
    xi=0.0,
    kwargs...
)
    p = model.params
    if p.N_color != 3
        throw(ArgumentError("thermal_contribution currently supports N_color=3, got $(p.N_color)"))
    end
    if p.N_flavor != 3
        throw(ArgumentError("thermal_contribution currently supports N_flavor=3, got $(p.N_flavor)"))
    end

    # NJL: ignore Φ / Φbar (they are only meaningful for PNJL family)
    _ = (Φ, Φbar, kwargs)

    p_mesh, cosθ_mesh, coeff = _cached_nodes(p_num, t_num)
    # IMPORTANT: during implicit differentiation we may see mixed ForwardDiff.Dual tags
    # (e.g. Duals seeded for φ inside ForwardDiff.gradient, and Duals seeded for (T, μ)
    # outside). Avoid explicit `convert(Tm, ...)` which can force invalid conversions;
    # instead, promote via normal arithmetic and accumulate in a type that can hold
    # nested duals.
    Tp = typeof(masses[1] + T + mu_vec[1] + xi)
    total = zero(Tp)
    T_prom = T + zero(Tp)

    @inbounds for flavor in 1:3
        m = masses[flavor]
        μ = mu_vec[flavor] + zero(Tp)
        for idx in eachindex(p_mesh)
            E = _energy_anisotropic(m, p_mesh[idx], xi, cosθ_mesh[idx])
            a = -(E - μ) / T_prom
            b = -(E + μ) / T_prom
            total += (_log1pexp(a) + _log1pexp(b)) * coeff[idx]
        end
    end

    # Ω_thermal = -2 * Nc * T * ∑_i ∫ d^3p/(2π)^3 [ln(1+e^{-(E-μ)/T}) + ln(1+e^{-(E+μ)/T})]
    return -2 * p.N_color * T_prom * total
end

@inline function _fermi_dirac(x)
    return inv(one(x) + exp(x))
end

"""返回 (quark, antiquark) 数密度（单位 fm⁻³）。

约定：返回 `(quark=SVector{3}(...), antiquark=SVector{3}(...))`。
"""
function number_densities(
    model::NJLModel,
    x_state,
    T,
    mu_vec;
    p_num::Int=DEFAULT_MOMENTUM_COUNT,
    t_num::Int=DEFAULT_THETA_COUNT,
    xi=0.0,
    kwargs...
)
    _ = kwargs

    st = x_state isa MeanFieldState ? x_state : MeanFieldState(x_state)
    φ = st.phi
    masses = calculate_mass_vec(model, φ)
    mu_vec = normalize_mu_vec(mu_vec)

    p_mesh, cosθ_mesh, coeff = _cached_nodes(p_num, t_num)

    p = model.params
    Tp = typeof(masses[1] + T + mu_vec[1] + xi)
    T_prom = T + zero(Tp)
    pref = 2 * p.N_color

    acc_q = MVector{3, Tp}(zero(Tp), zero(Tp), zero(Tp))
    acc_aq = MVector{3, Tp}(zero(Tp), zero(Tp), zero(Tp))

    @inbounds for flavor in 1:3
        m = masses[flavor]
        μ = mu_vec[flavor] + zero(Tp)
        total_q = zero(Tp)
        total_aq = zero(Tp)
        for idx in eachindex(p_mesh)
            E = _energy_anisotropic(m, p_mesh[idx], xi, cosθ_mesh[idx])
            f_q = _fermi_dirac((E - μ) / T_prom)
            f_aq = _fermi_dirac((E + μ) / T_prom)
            total_q += f_q * coeff[idx]
            total_aq += f_aq * coeff[idx]
        end
        acc_q[flavor] = pref * total_q
        acc_aq[flavor] = pref * total_aq
    end

    return (quark=SVector{3}(acc_q), antiquark=SVector{3}(acc_aq))
end
