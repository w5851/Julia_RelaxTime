"""NJL2Model

两味 NJL 模型（u,d），基于多重派发接入 models 子系统。

兼容策略：保持 `MeanFieldState.phi` 为 3 分量，约定 φ = (φu, φd, 0)。
"""

include(joinpath(@__DIR__, "NJL2Core.jl"))

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

const _NJL2_NODE_CACHE = Dict{Tuple{Int, Int}, NTuple{3, Matrix{Float64}}}()

@inline function _njl2_theta_nodes(t_num::Int)
    if t_num == DEFAULT_THETA_COUNT
        return DEFAULT_COSΘ_HALF_NODES, DEFAULT_COSΘ_HALF_WEIGHTS .* 2.0
    end
    nodes, weights = gauleg(0.0, 1.0, t_num)
    return nodes, weights .* 2.0
end

@inline function _njl2_momentum_nodes(p_num::Int)
    if p_num == DEFAULT_MOMENTUM_COUNT
        return DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS
    end
    return gauleg(0.0, 10.0, p_num)
end

@inline function _njl2_cached_nodes(p_num::Int, t_num::Int)
    key = (p_num, t_num)
    return get!(_NJL2_NODE_CACHE, key) do
        momentum_nodes, momentum_weights = _njl2_momentum_nodes(p_num)
        cosθ_nodes, cosθ_weights = _njl2_theta_nodes(t_num)

        p_mesh = repeat(momentum_nodes, 1, t_num)
        cosθ_mesh = repeat(cosθ_nodes', p_num, 1)
        weight_mesh = momentum_weights * cosθ_weights'
        coefficients = weight_mesh .* p_mesh.^2 ./ (2 * π)^2
        return (p_mesh, cosθ_mesh, coefficients)
    end
end

@inline function _njl2_energy_anisotropic(mass, p, xi, cosθ)
    return sqrt(p^2 + mass^2 + xi * (p * cosθ)^2)
end

@inline _njl2_log1pexp(x) = log1p(exp(x))

@inline function _njl2_vacuum_integral(Λ, mass)
    mass_abs = abs(mass)
    epsilon = one(mass_abs) * 1e-12
    mass_safe = mass_abs + epsilon
    sqrt_term = sqrt(Λ^2 + mass_safe^2)
    poly_part = Λ * sqrt_term * (2 * Λ^2 + mass_safe^2)
    log_term = mass_safe^4 * log((Λ + sqrt_term) / mass_safe)
    return (poly_part - log_term) / (16 * π^2)
end

struct NJL2Model <: AbstractNJLModel
    params::NJL2Core.NJL2Params
end

function NJL2Model(; profile::String=get(ENV, "NJL2_PARAM_PROFILE", "default"))
    return NJL2Model(NJL2Core.njl2_params(profile=profile))
end

function calculate_mass_vec(model::NJL2Model, φ; kwargs...)
    _ = kwargs
    return NJL2Core.calculate_mass_vec(model.params, φ)
end

function calculate_chiral(model::NJL2Model, φ; kwargs...)
    _ = kwargs
    return NJL2Core.chiral_potential(model.params, φ)
end

function polyakov_potential(::NJL2Model, Φ, Φbar, T; kwargs...)
    _ = kwargs
    return zero(promote_type(typeof(Φ), typeof(Φbar), typeof(T)))
end

function vacuum_contribution(model::NJL2Model, masses; kwargs...)
    _ = kwargs
    p = model.params
    p.N_flavor == 2 || throw(ArgumentError("NJL2Model expects N_flavor=2, got $(p.N_flavor)"))

    Tm = eltype(masses)
    Λ = p.Lambda_inv_fm * one(Tm)
    total = zero(Tm)
    @inbounds for i in 1:2
        total += _njl2_vacuum_integral(Λ, masses[i])
    end
    return -2 * convert(Tm, p.N_color) * total
end

function thermal_contribution(
    model::NJL2Model,
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
    _ = (Φ, Φbar, kwargs)

    p = model.params
    p.N_flavor == 2 || throw(ArgumentError("NJL2Model expects N_flavor=2, got $(p.N_flavor)"))

    p_mesh, cosθ_mesh, coeff = _njl2_cached_nodes(p_num, t_num)
    Tp = typeof(masses[1] + T + mu_vec[1] + xi)
    total = zero(Tp)
    T_prom = T + zero(Tp)

    @inbounds for flavor in 1:2
        m = masses[flavor]
        μ = mu_vec[flavor] + zero(Tp)
        for idx in eachindex(p_mesh)
            E = _njl2_energy_anisotropic(m, p_mesh[idx], xi, cosθ_mesh[idx])
            a = -(E - μ) / T_prom
            b = -(E + μ) / T_prom
            total += (_njl2_log1pexp(a) + _njl2_log1pexp(b)) * coeff[idx]
        end
    end

    return -2 * p.N_color * T_prom * total
end

@inline _njl2_fermi_dirac(x) = inv(one(x) + exp(x))

function number_densities(
    model::NJL2Model,
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

    p_mesh, cosθ_mesh, coeff = _njl2_cached_nodes(p_num, t_num)

    p = model.params
    Tp = typeof(masses[1] + T + mu_vec[1] + xi)
    T_prom = T + zero(Tp)
    pref = 2 * p.N_color

    acc_q = MVector{2, Tp}(zero(Tp), zero(Tp))
    acc_aq = MVector{2, Tp}(zero(Tp), zero(Tp))

    @inbounds for flavor in 1:2
        m = masses[flavor]
        μ = mu_vec[flavor] + zero(Tp)
        total_q = zero(Tp)
        total_aq = zero(Tp)
        for idx in eachindex(p_mesh)
            E = _njl2_energy_anisotropic(m, p_mesh[idx], xi, cosθ_mesh[idx])
            f_q = _njl2_fermi_dirac((E - μ) / T_prom)
            f_aq = _njl2_fermi_dirac((E + μ) / T_prom)
            total_q += f_q * coeff[idx]
            total_aq += f_aq * coeff[idx]
        end
        acc_q[flavor] = pref * total_q
        acc_aq[flavor] = pref * total_aq
    end

    return (
        quark=SVector{3, Tp}(acc_q[1], acc_q[2], zero(Tp)),
        antiquark=SVector{3, Tp}(acc_aq[1], acc_aq[2], zero(Tp)),
    )
end
