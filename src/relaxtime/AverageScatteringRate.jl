module AverageScatteringRate

"""
Thermal-weighted, positive cross-section interpolation with an explicit error
budget. Momentum integration is semi-infinite; the physical cross-section
support retains the model cutoff. Use Models.TransportSamplingProfiles for the
qualified background, A integrals, thermal average and terminal transport.
"""

const _VALIDATION_UTILS_PATH = normpath(joinpath(@__DIR__, "..", "utils", "ValidationUtils.jl"))
if !isdefined(Main, :ValidationUtils)
    Base.include(Main, _VALIDATION_UTILS_PATH)
end

using Base.Threads: ReentrantLock

using Main.Constants_PNJL: Λ_inv_fm
using ..GaussLegendre: gauleg
using ..GaussLegendre: standard_nodes_weights
using ..AFieldBuilder: ensure_quark_params_has_A
import ..TotalCrossSection
using ..TotalCrossSection: total_cross_section
using ..TotalCrossSection: parse_particles_from_process
using ..ParticleSymbols: is_antiquark
using Main.PNJLQuarkDistributions: quark_distribution, antiquark_distribution
using Main.PNJLQuarkDistributions_Aniso: quark_distribution_aniso, antiquark_distribution_aniso

using Main.ParameterTypes: QuarkParams, ThermoParams
using Main.ParameterAdapters: normalize_quark_input, normalize_thermo_input
using Main.ValidationUtils: validate_grid_weight_pair

include("ThermalAngularQuadrature.jl")
using .ThermalAngularQuadrature: angular_groups
include("ThermalWeightedSampling.jl")

export average_scattering_rate

const DEFAULT_P_NODES = 48
const DEFAULT_ANGLE_NODES = 36
const DEFAULT_PHI_NODES = 36
const DEFAULT_N_SIGMA_POINTS = 24
# Automatic A construction retains its own quadrature convention.
const DEFAULT_A_COS_NODES = 4
const DQ = 6.0
const TWO_PI = 2.0 * π
const PROPAGATOR_XI_POLICIES = (:match_thermo, :isotropic)
const SIGMA_CACHE_POLICIES = (:thermal_weighted,)
const RELATIVE_VELOCITY_POLICIES = (:cm, :medium_moller)

@inline function distribution_with_anisotropy(flavor::Symbol, p::Float64, m::Float64, μ::Float64,
    T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, cosθ::Float64)
    if ξ == 0.0
        E = energy_from_p(p, m)
        return is_antiquark(flavor) ? antiquark_distribution(E, μ, T, Φ, Φbar) : quark_distribution(E, μ, T, Φ, Φbar)
    else
        return is_antiquark(flavor) ? antiquark_distribution_aniso(p, m, μ, T, Φ, Φbar, ξ, cosθ) : quark_distribution_aniso(p, m, μ, T, Φ, Φbar, ξ, cosθ)
    end
end

@inline function energy_from_p(p::Float64, m::Float64)
    return sqrt(p * p + m * m)
end

@inline function _normalize_propagator_xi_policy(policy::Symbol)::Symbol
    policy in PROPAGATOR_XI_POLICIES ||
        throw(ArgumentError("unknown propagator_xi_policy=$(policy); expected one of $(PROPAGATOR_XI_POLICIES)"))
    return policy
end

@inline function _normalize_sigma_cache_policy(policy::Symbol)::Symbol
    policy in SIGMA_CACHE_POLICIES ||
        throw(ArgumentError("unknown sigma_cache_policy=$(policy); expected one of $(SIGMA_CACHE_POLICIES)"))
    return policy
end

@inline function _normalize_relative_velocity_policy(policy)::Symbol
    normalized = policy isa Symbol ? policy : Symbol(policy)
    normalized in RELATIVE_VELOCITY_POLICIES ||
        throw(ArgumentError("unknown relative_velocity=$(policy); expected one of $(RELATIVE_VELOCITY_POLICIES)"))
    return normalized
end

@inline function _with_xi(thermo_params::NamedTuple, ξ::Float64)::NamedTuple
    return merge(thermo_params, (ξ=ξ,))
end

@inline function _base_quark_params_without_A(quark_params::NamedTuple)::NamedTuple
    return (m=quark_params.m, μ=quark_params.μ)
end

@inline function _resolve_propagator_thermo_params(thermo_params::NamedTuple, policy::Symbol)::NamedTuple
    policy === :isotropic && return _with_xi(thermo_params, 0.0)
    return thermo_params
end

function _resolve_propagator_quark_params(
    quark_params::NamedTuple,
    propagator_thermo_params::NamedTuple,
    policy::Symbol,
    propagator_quark_params::Union{Nothing,NamedTuple,QuarkParams},
)::NamedTuple
    if propagator_quark_params !== nothing
        return normalize_quark_input(propagator_quark_params)
    elseif policy === :isotropic
        return ensure_quark_params_has_A(
            _base_quark_params_without_A(quark_params),
            propagator_thermo_params;
            use_aniso=false,
            warn_on_auto=false,
        )
    end
    return quark_params
end


@inline function get_mass(flavor::Symbol, quark_params::NamedTuple)
    if flavor === :u || flavor === :ubar; return quark_params.m.u
    elseif flavor === :d || flavor === :dbar; return quark_params.m.d
    elseif flavor === :s || flavor === :sbar; return quark_params.m.s
    else; error("Unknown flavor $flavor") end
end

@inline get_mass(flavor::Symbol, quark_params::QuarkParams) = get_mass(flavor, normalize_quark_input(quark_params))

@inline function get_mu(flavor::Symbol, quark_params::NamedTuple)
    # Convention: always return the quark chemical potential μ_q (positive sign).
    # The particle/antiparticle distinction is handled by using
    # `quark_distribution*` vs `antiquark_distribution*`.
    if flavor === :u || flavor === :ubar; return quark_params.μ.u
    elseif flavor === :d || flavor === :dbar; return quark_params.μ.d
    elseif flavor === :s || flavor === :sbar; return quark_params.μ.s
    else; error("Unknown flavor $flavor") end
end

@inline get_mu(flavor::Symbol, quark_params::QuarkParams) = get_mu(flavor, normalize_quark_input(quark_params))

const _FNV64_OFFSET = UInt64(0xcbf29ce484222325)
const _FNV64_PRIME = UInt64(0x00000100000001b3)

@inline _fingerprint_value(x::Nothing) = nothing
@inline _fingerprint_value(x::Bool) = x
@inline _fingerprint_value(x::Symbol) = x
@inline _fingerprint_value(x::Unsigned) = x
@inline _fingerprint_value(x::Integer) = Int(x)
@inline _fingerprint_value(x::AbstractFloat) = Float64(x)
@inline _fingerprint_value(x::Real) = Float64(x)
_fingerprint_value(x::Tuple) = map(_fingerprint_value, x)

function _fingerprint_value(x::NamedTuple)
    names = sort(collect(keys(x)); by=String)
    return Tuple(name => _fingerprint_value(getproperty(x, name)) for name in names)
end

@inline function _fnv64_mix(h::UInt64, x::UInt64)::UInt64
    return (h ⊻ x) * _FNV64_PRIME
end

function _float_vector_fingerprint_hash(values::Vector{Float64})::UInt64
    h = _FNV64_OFFSET
    for value in values
        h = _fnv64_mix(h, reinterpret(UInt64, value))
    end
    return h
end

function _s_grid_fingerprint(values::Vector{Float64})
    return (
        length=length(values),
        first=isempty(values) ? nothing : values[1],
        last=isempty(values) ? nothing : values[end],
        hash=_float_vector_fingerprint_hash(values),
    )
end

function _quark_fingerprint(quark_params::NamedTuple)
    return (
        m=_fingerprint_value(quark_params.m),
        μ=_fingerprint_value(quark_params.μ),
        A=hasproperty(quark_params, :A) ? _fingerprint_value(quark_params.A) : nothing,
    )
end

function _thermo_fingerprint(thermo_params::NamedTuple)
    return (
        T=_fingerprint_value(thermo_params.T),
        Φ=_fingerprint_value(thermo_params.Φ),
        Φbar=_fingerprint_value(thermo_params.Φbar),
        ξ=hasproperty(thermo_params, :ξ) ? _fingerprint_value(thermo_params.ξ) : 0.0,
    )
end

const DEFAULT_SEMI_INF_SCALE = 10.0  # 半无穷积分的尺度参数
const DEFAULT_SEMI_INF_NODES = 32    # 半无穷积分的节点数

const _INTERVAL_GRID_CACHE = Dict{Tuple{Float64,Float64,Int},NTuple{2,Vector{Float64}}}()
const _INTERVAL_GRID_LOCK = ReentrantLock()
const _SEMI_INF_GRID_CACHE = Dict{Tuple{Int,Float64},NTuple{2,Vector{Float64}}}()
const _SEMI_INF_GRID_LOCK = ReentrantLock()

@inline function _build_interval_grid(a::Float64, b::Float64, n::Int)
    nodes_std, weights_std = standard_nodes_weights(n)
    nodes = similar(nodes_std)
    weights = similar(weights_std)
    scale = 0.5 * (b - a)
    shift = 0.5 * (b + a)
    @inbounds @simd for i in eachindex(nodes_std)
        nodes[i] = scale * nodes_std[i] + shift
        weights[i] = scale * weights_std[i]
    end
    return nodes, weights
end

function _cached_interval_grid(a::Float64, b::Float64, n::Int)
    if a == 0.0 && b == 1.0 && n == DEFAULT_ANGLE_NODES
        return _DEFAULT_HALF_COS_GRID
    elseif a == -1.0 && b == 1.0 && n == DEFAULT_ANGLE_NODES
        return _DEFAULT_FULL_COS_GRID
    elseif a == 0.0 && b == TWO_PI && n == DEFAULT_PHI_NODES
        return _DEFAULT_PHI_GRID
    end

    key = (a, b, n)
    lock(_INTERVAL_GRID_LOCK)
    try
        return get!(_INTERVAL_GRID_CACHE, key) do
            _build_interval_grid(a, b, n)
        end
    finally
        unlock(_INTERVAL_GRID_LOCK)
    end
end

function _build_semi_infinite_momentum_grid(p_nodes::Int, scale::Float64)
    t_grid, t_w = gauleg(0.0, 1.0, p_nodes)
    p_vals = Float64[]
    quadrature_wts = Float64[]
    sizehint!(p_vals, length(t_grid))
    sizehint!(quadrature_wts, length(t_w))

    for (t, wt) in zip(t_grid, t_w)
        if t >= 0.9999
            continue
        end
        inv_gap = 1.0 / (1.0 - t)
        p = scale * t * inv_gap
        push!(p_vals, p)
        push!(quadrature_wts, wt * scale * inv_gap^2)
    end
    return p_vals, quadrature_wts
end

function _cached_semi_infinite_momentum_grid(p_nodes::Int, scale::Float64)
    if p_nodes == DEFAULT_SEMI_INF_NODES && scale == DEFAULT_SEMI_INF_SCALE
        return _DEFAULT_DENSITY_SEMI_INF_GRID
    elseif p_nodes == DEFAULT_P_NODES && scale == DEFAULT_SEMI_INF_SCALE
        return _DEFAULT_ASR_SEMI_INF_GRID
    end

    key = (p_nodes, scale)
    lock(_SEMI_INF_GRID_LOCK)
    try
        return get!(_SEMI_INF_GRID_CACHE, key) do
            _build_semi_infinite_momentum_grid(p_nodes, scale)
        end
    finally
        unlock(_SEMI_INF_GRID_LOCK)
    end
end

const _DEFAULT_HALF_COS_GRID = _build_interval_grid(0.0, 1.0, DEFAULT_ANGLE_NODES)
const _DEFAULT_FULL_COS_GRID = _build_interval_grid(-1.0, 1.0, DEFAULT_ANGLE_NODES)
const _DEFAULT_PHI_GRID = _build_interval_grid(0.0, TWO_PI, DEFAULT_PHI_NODES)
const _DEFAULT_DENSITY_SEMI_INF_GRID = _build_semi_infinite_momentum_grid(DEFAULT_SEMI_INF_NODES, DEFAULT_SEMI_INF_SCALE)
const _DEFAULT_ASR_SEMI_INF_GRID = _build_semi_infinite_momentum_grid(DEFAULT_P_NODES, DEFAULT_SEMI_INF_SCALE)

"""
    number_density(flavor, m, μ, T, Φ, Φbar, ξ; kwargs...)

计算夸克/反夸克数密度，使用半无穷积分 [0, ∞)。

积分变换: p = scale * t / (1-t), dp = scale / (1-t)^2 dt
其中 t ∈ [0, 1)

# Arguments
- `flavor::Symbol`: 夸克味道 (:u, :d, :s, :ubar, :dbar, :sbar)
- `m::Float64`: 夸克质量 (fm⁻¹)
- `μ::Float64`: 化学势 (fm⁻¹)
- `T::Float64`: 温度 (fm⁻¹)
- `Φ::Float64`: Polyakov loop
- `Φbar::Float64`: Polyakov loop conjugate
- `ξ::Float64`: 各向异性参数

# Keyword Arguments
- `p_nodes::Int`: 动量积分节点数 (默认32)
- `angle_nodes::Int`: 角度积分节点数 (默认36)
- `scale::Float64`: 半无穷积分尺度参数 (默认10.0)
"""
function number_density(flavor::Symbol, m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64;
    p_nodes::Int=DEFAULT_SEMI_INF_NODES, angle_nodes::Int=DEFAULT_ANGLE_NODES,
    p_grid::Union{Nothing,Vector{Float64}}=nothing, p_w::Union{Nothing,Vector{Float64}}=nothing,
    cos_grid::Union{Nothing,Vector{Float64}}=nothing, cos_w::Union{Nothing,Vector{Float64}}=nothing,
    scale::Float64=DEFAULT_SEMI_INF_SCALE)
    validate_grid_weight_pair("number_density", "p_grid", p_grid, "p_w", p_w)
    validate_grid_weight_pair("number_density", "cos_grid", cos_grid, "cos_w", cos_w)

    cos_grid === nothing && ((cos_grid, cos_w) = _cached_interval_grid(0.0, 1.0, angle_nodes))

    # Momentum integration:
    # - Default: semi-infinite [0, ∞) via t∈[0,1) mapping.
    # - If p_grid/p_w are provided: integrate directly on that finite grid.
    if p_grid !== nothing && p_w !== nothing
        p_vals = p_grid
        integral = 0.0
        for (p, wp) in zip(p_vals, p_w)
            for (cθ, wθ) in zip(cos_grid, cos_w)
                f = distribution_with_anisotropy(flavor, p, m, μ, T, Φ, Φbar, ξ, cθ)
                integral += wp * wθ * p^2 * f
            end
        end
        return DQ * integral / (2.0 * π^2)
    else
        p_vals, quadrature_wts = _cached_semi_infinite_momentum_grid(p_nodes, scale)
    end

    integral = 0.0
    for (p, wp) in zip(p_vals, quadrature_wts)
        for (cθ, wθ) in zip(cos_grid, cos_w)
            f = distribution_with_anisotropy(flavor, p, m, μ, T, Φ, Φbar, ξ, cθ)
            integral += wp * wθ * p^2 * f
        end
    end
    # ρ = d_q / (2π^2) ∫ p^2 dp ∫_0^1 dcosθ f
    return DQ * integral / (2.0 * π^2)
end

@inline function _resolve_sigma_support_bounds(
    mi::Float64,
    mj::Float64,
    mc::Float64,
    md::Float64,
    sigma_cutoff::Union{Nothing,Float64},
    p_grid::Union{Nothing,Vector{Float64}},
    p_w::Union{Nothing,Vector{Float64}},
)::Tuple{Float64,Float64,Float64}
    has_finite_p_cut = (p_grid !== nothing && p_w !== nothing)
    Λ = if sigma_cutoff !== nothing
        sigma_cutoff
    elseif has_finite_p_cut
        maximum(p_grid)
    else
        NaN
    end

    s_bo = max((mi + mj)^2, (mc + md)^2)
    s_up = if !isnan(Λ)
        min((sqrt(mi^2 + Λ^2) + sqrt(mj^2 + Λ^2))^2,
            (sqrt(mc^2 + Λ^2) + sqrt(md^2 + Λ^2))^2)
    else
        Inf
    end
    return Λ, s_bo, s_up
end

@inline function _relative_velocity_cm(s::Float64, mi::Float64, mj::Float64)
    # Exactly equal incoming masses admit a factored threshold expression.
    # Do not round nearly equal masses into this branch.
    if mi == mj
        return 2.0 * sqrt(max(0.0, (s - 4.0 * mi^2) / s))
    end
    s_rt = sqrt(s)
    Ei_cm = (s + mi^2 - mj^2) / (2.0 * s_rt)
    Ej_cm = (s - mi^2 + mj^2) / (2.0 * s_rt)
    p_cm = sqrt(max(0.0, (s - (mi + mj)^2) * (s - (mi - mj)^2))) / (2.0 * s_rt)
    numerator = (Ei_cm * Ej_cm + p_cm * p_cm)^2 - (mi * mj)^2
    return numerator > 0.0 ? sqrt(numerator) / (Ei_cm * Ej_cm) : 0.0
end

"""Møller relative speed evaluated with the incoming medium-frame energies.

The invariant numerator is ``sqrt(lambda(s,mi^2,mj^2))/2``.  The current
historical cache uses the same numerator divided by the two centre-of-mass
energies; the paper's medium-frame definition divides by the actual incoming
energies at each thermal event.  Keeping this as a separate helper prevents
the two conventions from being mixed after the event has been compressed to
its invariant ``s`` value.
"""
@inline function _relative_velocity_medium_moller(
    s::Float64, mi::Float64, mj::Float64, ei::Float64, ej::Float64,
)
    denominator = 2.0 * ei * ej
    denominator > 0.0 || return 0.0
    lambda = (s - (mi + mj)^2) * (s - (mi - mj)^2)
    return sqrt(max(0.0, lambda)) / denominator
end

@inline function _relative_velocity_for_mode(
    policy::Symbol, s::Float64, mi::Float64, mj::Float64, ei::Float64, ej::Float64,
)
    policy === :cm && return _relative_velocity_cm(s, mi, mj)
    policy === :medium_moller && return _relative_velocity_medium_moller(s, mi, mj, ei, ej)
    throw(ArgumentError("unknown relative_velocity=$(policy); expected one of $(RELATIVE_VELOCITY_POLICIES)"))
end

@inline function _relative_velocity_for_mode(
    policy, s::Float64, mi::Float64, mj::Float64, ei::Float64, ej::Float64,
)
    return _relative_velocity_for_mode(_normalize_relative_velocity_policy(policy), s, mi, mj, ei, ej)
end

include("ThermalWeightedCache.jl")

"""
    average_scattering_rate(process, quark_params, thermo_params, K_coeffs; kwargs...)

Build or consume a ThermalWeightedCache. The cache identity includes the
physical context, all effective nodes and weights, cutoff and error budget.
Unconverged construction throws; no unchecked table or legacy subtraction
fallback is used. Custom grids require independent outer-quadrature validation.
"""
function average_scattering_rate(
    process::Symbol,
    quark_params::Union{NamedTuple, QuarkParams},
    thermo_params::Union{NamedTuple, ThermoParams},
    K_coeffs::NamedTuple;
    p_nodes::Int=DEFAULT_P_NODES,
    angle_nodes::Int=DEFAULT_ANGLE_NODES,
    phi_nodes::Int=DEFAULT_PHI_NODES,
    p_grid::Union{Nothing,Vector{Float64}}=nothing, p_w::Union{Nothing,Vector{Float64}}=nothing,
    cos_grid::Union{Nothing,Vector{Float64}}=nothing, cos_w::Union{Nothing,Vector{Float64}}=nothing,
    phi_grid::Union{Nothing,Vector{Float64}}=nothing, phi_w::Union{Nothing,Vector{Float64}}=nothing,
    cs_cache::Union{Nothing,ThermalWeightedCache}=nothing,
    n_sigma_points::Int=DEFAULT_N_SIGMA_POINTS,
    scale::Float64=DEFAULT_SEMI_INF_SCALE,
    density_p_grid::Union{Nothing,Vector{Float64}}=nothing,
    density_p_w::Union{Nothing,Vector{Float64}}=nothing,
    density_p_nodes::Int=DEFAULT_SEMI_INF_NODES,
    density_scale::Float64=DEFAULT_SEMI_INF_SCALE,
    apply_s_domain_cut::Bool=true,
    sigma_cutoff::Union{Nothing,Float64}=nothing,  # 新增：σ(s)有效范围的动量截断
    interpolation_mode::Symbol=:auto,
    propagator_xi_policy::Symbol=:match_thermo,
    propagator_quark_params::Union{Nothing,NamedTuple,QuarkParams}=nothing,
    sigma_cache_policy::Symbol=:thermal_weighted,
    angular_reduction::Symbol=:grouped,
    sampling_config::ThermalSamplingConfig=DEFAULT_THERMAL_SAMPLING_CONFIG,
)::Float64
    _normalize_sigma_cache_policy(sigma_cache_policy)
    interpolation_mode in (:auto, :thermal_weighted) ||
        throw(ArgumentError("only thermal-weighted interpolation is supported"))
    angular_reduction === :grouped ||
        throw(ArgumentError("thermal-weighted sampling requires exact grouped angular quadrature"))
    return _thermal_weighted_average(process,quark_params,thermo_params,K_coeffs,cs_cache;
        p_nodes=p_nodes,angle_nodes=angle_nodes,phi_nodes=phi_nodes,
        p_grid=p_grid,p_w=p_w,cos_grid=cos_grid,cos_w=cos_w,phi_grid=phi_grid,phi_w=phi_w,
        scale=scale,density_p_grid=density_p_grid,density_p_w=density_p_w,
        density_p_nodes=density_p_nodes,density_scale=density_scale,
        sigma_cutoff=sigma_cutoff,apply_s_domain_cut=apply_s_domain_cut,n_sigma_points=n_sigma_points,
        propagator_xi_policy=propagator_xi_policy,propagator_quark_params=propagator_quark_params,
        sampling_config=sampling_config)
end

end # module AverageScatteringRate
