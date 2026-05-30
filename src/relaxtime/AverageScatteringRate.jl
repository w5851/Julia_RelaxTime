module AverageScatteringRate

"""
# AverageScatteringRate Module

平均散射率计算模块（各向异性可选）。

## Features

- 积分采用 Gauss-Legendre：动量 32 节点，角度 4 节点（可覆盖 7 阶多项式）。
- 支持各向异性分布 `quark_distribution_aniso(..., ξ, cosθ)`；`ξ=0` 时退化为各向同性。
- 散射截面使用 `TotalCrossSection.total_cross_section`，可按需预计算并插值。

## Cross-section Cache

- 本仓库已将生产默认策略固定为 **w0cdf 取点 + PCHIP 插值**。
- `CrossSectionCache` 仅用于承载预计算的 σ(s) 表，并用 PCHIP 做插值；当质心能量 s 超出缓存覆盖区间时，直接返回 0（而不是钳制到边界），以避免不可达的 s 区域产生伪贡献。
- `CrossSectionCache(process)` 创建空缓存；通过 `precompute_cross_section!` 填充后即可用于 `average_scattering_rate`。

## Dual Interface Pattern

This module supports **both struct and NamedTuple parameters**:

```julia
# Using structs (recommended)
using Main.ParameterTypes: QuarkParams, ThermoParams

q = QuarkParams(m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
t = ThermoParams(0.15, 0.5, 0.5, 0.0)
rate = average_scattering_rate(:uu_to_uu, q, t, K_coeffs)

# Using NamedTuples (backward compatible)
q_nt = (m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
t_nt = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
rate = average_scattering_rate(:uu_to_uu, q_nt, t_nt, K_coeffs)
```

Both produce identical results. Internal normalization ensures type stability and zero overhead.
"""

# Include-once helper
# Dependencies loaded by RelaxTime.jl entry point

const _VALIDATION_UTILS_PATH = normpath(joinpath(@__DIR__, "..", "utils", "ValidationUtils.jl"))
if !isdefined(Main, :ValidationUtils)
    Base.include(Main, _VALIDATION_UTILS_PATH)
end

using LinearAlgebra
using Statistics
using Base.Threads: ReentrantLock

using Main.Constants_PNJL: Λ_inv_fm
using ..GaussLegendre: gauleg
using ..GaussLegendre: standard_nodes_weights
import ..TotalCrossSection
using ..TotalCrossSection: total_cross_section
using ..TotalCrossSection: parse_particles_from_process
using ..ParticleSymbols: is_antiquark
using Main.PNJLQuarkDistributions: quark_distribution, antiquark_distribution
using Main.PNJLQuarkDistributions_Aniso: quark_distribution_aniso, antiquark_distribution_aniso

using Main.ParameterTypes: QuarkParams, ThermoParams
using Main.ParameterAdapters: normalize_quark_input, normalize_thermo_input
using Main.ValidationUtils: validate_grid_weight_pair

export average_scattering_rate, CrossSectionCache, precompute_cross_section!, build_w0cdf_pchip_cache

const DEFAULT_P_NODES = 20
const DEFAULT_ANGLE_NODES = 4  # cosθ节点数
const DEFAULT_PHI_NODES = 8    # φ节点数
const DEFAULT_W0CDF_P_NODES = 14
const DEFAULT_W0CDF_ANGLE_NODES = DEFAULT_ANGLE_NODES
const DEFAULT_W0CDF_PHI_NODES = DEFAULT_PHI_NODES
const DEFAULT_SIGMA_GRID_N = 60
const DQ = 6.0 # 简并度d_q=2*N_c=6
const TWO_PI = 2.0 * π

# --------------------------- 工具函数 ---------------------------

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

# -------------------- 截面缓存与插值 --------------------
mutable struct CrossSectionData
    s_vals::Vector{Float64}
    sigma_vals::Vector{Float64}
end

mutable struct CachedInterpolation
    pchip_slopes::Vector{Float64}
    pchip_dirty::Bool
    peak_ratio::Float64
    peak_s::Float64
    peak_dirty::Bool
    local_s_vals::Vector{Float64}
    local_sigma_vals::Vector{Float64}
    local_dirty::Bool
    local_n_points::Int
    local_upper_s::Float64
end

mutable struct AsymptoticConfig
    asym_enabled::Bool
    asym_s0::Float64
    asym_A::Float64
    asym_requested::Bool
end

mutable struct CrossSectionCache
    process::Symbol
    data::CrossSectionData
    interpolation::CachedInterpolation
    asymptotic::AsymptoticConfig
    fingerprint::Union{Nothing,NamedTuple}
end

CrossSectionData() = CrossSectionData(Float64[], Float64[])
CachedInterpolation() = CachedInterpolation(Float64[], true, 0.0, 0.0, true, Float64[], Float64[], true, 0, 0.0)
AsymptoticConfig() = AsymptoticConfig(false, 0.0, 0.0, false)
CrossSectionCache(process::Symbol, data::CrossSectionData, interpolation::CachedInterpolation, asymptotic::AsymptoticConfig) =
    CrossSectionCache(process, data, interpolation, asymptotic, nothing)

function Base.propertynames(::CrossSectionCache, private::Bool=false)
    names = (:process, :data, :interpolation, :asymptotic, :fingerprint,
        :s_vals, :sigma_vals, :pchip_slopes, :pchip_dirty,
        :peak_ratio, :peak_s, :peak_dirty,
        :local_s_vals, :local_sigma_vals, :local_dirty, :local_n_points, :local_upper_s,
        :asym_enabled, :asym_s0, :asym_A, :asym_requested)
    return private ? names : names
end

function Base.getproperty(cache::CrossSectionCache, name::Symbol)
    if name === :process || name === :data || name === :interpolation || name === :asymptotic || name === :fingerprint
        return getfield(cache, name)
    elseif name === :s_vals || name === :sigma_vals
        return getfield(getfield(cache, :data), name)
    elseif name === :pchip_slopes || name === :pchip_dirty
        return getfield(getfield(cache, :interpolation), name)
    elseif name === :peak_ratio || name === :peak_s || name === :peak_dirty
        return getfield(getfield(cache, :interpolation), name)
    elseif name === :local_s_vals || name === :local_sigma_vals || name === :local_dirty || name === :local_n_points || name === :local_upper_s
        return getfield(getfield(cache, :interpolation), name)
    elseif name === :asym_enabled || name === :asym_s0 || name === :asym_A || name === :asym_requested
        return getfield(getfield(cache, :asymptotic), name)
    end
    return getfield(cache, name)
end

function Base.setproperty!(cache::CrossSectionCache, name::Symbol, value)
    if name === :s_vals || name === :sigma_vals
        return setfield!(getfield(cache, :data), name, value)
    elseif name === :pchip_slopes || name === :pchip_dirty
        return setfield!(getfield(cache, :interpolation), name, value)
    elseif name === :peak_ratio || name === :peak_s || name === :peak_dirty
        return setfield!(getfield(cache, :interpolation), name, value)
    elseif name === :local_s_vals || name === :local_sigma_vals || name === :local_dirty || name === :local_n_points || name === :local_upper_s
        return setfield!(getfield(cache, :interpolation), name, value)
    elseif name === :asym_enabled || name === :asym_s0 || name === :asym_A || name === :asym_requested
        return setfield!(getfield(cache, :asymptotic), name, value)
    end
    return setfield!(cache, name, value)
end

CrossSectionCache(process::Symbol) = CrossSectionCache(process, CrossSectionData(), CachedInterpolation(), AsymptoticConfig())

const _CACHE_FINGERPRINT_VERSION = 1
const _CACHE_FINGERPRINT_CONTEXT_FIELDS = (
    :version,
    :process,
    :quark,
    :thermo,
    :K_coeffs,
    :n_points,
    :threshold_subtraction,
    :asym_window,
    :asym_fit_min_points,
    :asym_extra_points,
)
const _CACHE_FINGERPRINT_WARNED = IdDict{CrossSectionCache,Bool}()
const _CACHE_FINGERPRINT_WARN_LOCK = ReentrantLock()
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

function _cross_section_cache_fingerprint(
    process::Symbol,
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    K_coeffs::NamedTuple;
    n_points::Int,
    threshold_subtraction::Bool,
    asym_window::Float64,
    asym_fit_min_points::Int,
    asym_extra_points::Int,
    s_grid::Vector{Float64},
    grid_metadata::NamedTuple=(kind=:explicit,),
)
    return (
        version=_CACHE_FINGERPRINT_VERSION,
        process=process,
        quark=_quark_fingerprint(quark_params),
        thermo=_thermo_fingerprint(thermo_params),
        K_coeffs=_fingerprint_value(K_coeffs),
        n_points=Int(n_points),
        threshold_subtraction=Bool(threshold_subtraction),
        asym_window=Float64(asym_window),
        asym_fit_min_points=Int(asym_fit_min_points),
        asym_extra_points=Int(asym_extra_points),
        s_grid=_s_grid_fingerprint(s_grid),
        grid=_fingerprint_value(grid_metadata),
    )
end

function _explicit_grid_metadata(input_s_grid::Vector{Float64})
    return (
        kind=:explicit,
        input_s_grid=_s_grid_fingerprint(input_s_grid),
    )
end

function _w0cdf_grid_metadata(
    input_s_grid::Vector{Float64};
    N::Int,
    design_p_nodes::Int,
    design_angle_nodes::Int,
    design_phi_nodes::Int,
    p_cutoff::Union{Nothing,Float64},
    scale::Float64,
)
    return (
        kind=:w0cdf,
        N=Int(N),
        design_p_nodes=Int(design_p_nodes),
        design_angle_nodes=Int(design_angle_nodes),
        design_phi_nodes=Int(design_phi_nodes),
        p_cutoff=_fingerprint_value(p_cutoff),
        scale=Float64(scale),
        input_s_grid=_s_grid_fingerprint(input_s_grid),
    )
end

@inline function _fingerprint_grid_kind(fingerprint::NamedTuple)
    grid = fingerprint.grid
    return grid isa Tuple ? _fingerprint_grid_value(grid, :kind) : :unknown
end

function _fingerprint_grid_value(grid::Tuple, key::Symbol)
    for pair in grid
        pair.first === key && return pair.second
    end
    return nothing
end

function _fingerprint_context_mismatch(actual::NamedTuple, expected::NamedTuple)
    for field in _CACHE_FINGERPRINT_CONTEXT_FIELDS
        hasproperty(actual, field) || return field
        getproperty(actual, field) == getproperty(expected, field) || return field
    end
    return nothing
end

function _w0cdf_grid_mismatch(actual::NamedTuple, sigma_cutoff::Union{Nothing,Float64}, scale::Float64)
    _fingerprint_grid_kind(actual) === :w0cdf || return nothing

    expected = (
        kind=:w0cdf,
        N=DEFAULT_SIGMA_GRID_N,
        design_p_nodes=DEFAULT_W0CDF_P_NODES,
        design_angle_nodes=DEFAULT_W0CDF_ANGLE_NODES,
        design_phi_nodes=DEFAULT_W0CDF_PHI_NODES,
        p_cutoff=_fingerprint_value(sigma_cutoff),
        scale=Float64(scale),
    )
    for field in keys(expected)
        actual_value = _fingerprint_grid_value(actual.grid, field)
        actual_value == getproperty(expected, field) || return field
    end
    return nothing
end

function _warn_missing_cache_fingerprint_once(cache::CrossSectionCache)
    should_warn = false
    lock(_CACHE_FINGERPRINT_WARN_LOCK)
    try
        if !haskey(_CACHE_FINGERPRINT_WARNED, cache)
            _CACHE_FINGERPRINT_WARNED[cache] = true
            should_warn = true
        end
    finally
        unlock(_CACHE_FINGERPRINT_WARN_LOCK)
    end

    if should_warn
        @warn "CrossSectionCache has no fingerprint; allowing legacy/manual cache reuse. Rebuild the cache with precompute_cross_section! or build_w0cdf_pchip_cache for strict validation, or pass require_cache_fingerprint=true to reject fingerprint-less caches." process=cache.process
    end
    return nothing
end

function _validate_cross_section_cache!(
    cache::CrossSectionCache,
    process::Symbol,
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    K_coeffs::NamedTuple;
    n_points::Int,
    threshold_subtraction::Bool,
    asym_window::Float64,
    asym_fit_min_points::Int,
    asym_extra_points::Int,
    sigma_cutoff::Union{Nothing,Float64},
    scale::Float64,
    require_cache_fingerprint::Bool,
)
    cache.process === process || throw(ArgumentError("CrossSectionCache process mismatch: cache is for $(cache.process), but current process is $(process)"))
    isempty(cache.s_vals) && return cache

    if cache.fingerprint === nothing
        if require_cache_fingerprint
            throw(ArgumentError("CrossSectionCache for $(process) has no fingerprint; rebuild it with precompute_cross_section! or build_w0cdf_pchip_cache before strict reuse"))
        end
        _warn_missing_cache_fingerprint_once(cache)
        return cache
    end

    cache.fingerprint.s_grid == _s_grid_fingerprint(cache.s_vals) ||
        throw(ArgumentError("CrossSectionCache fingerprint for $(process) no longer matches the cached s_grid; rebuild the cache before reuse"))

    expected = _cross_section_cache_fingerprint(
        process,
        quark_params,
        thermo_params,
        K_coeffs;
        n_points=n_points,
        threshold_subtraction=threshold_subtraction,
        asym_window=asym_window,
        asym_fit_min_points=asym_fit_min_points,
        asym_extra_points=asym_extra_points,
        s_grid=cache.s_vals,
        grid_metadata=(kind=:validation_context,),
    )

    mismatch = _fingerprint_context_mismatch(cache.fingerprint, expected)
    if mismatch !== nothing
        throw(ArgumentError("CrossSectionCache fingerprint mismatch for $(process): $(mismatch) differs from current parameters; rebuild the cache for the current quark_params, thermo_params, K_coeffs, n_sigma_points, and threshold_subtraction settings"))
    end

    grid_mismatch = _w0cdf_grid_mismatch(cache.fingerprint, sigma_cutoff, scale)
    if grid_mismatch !== nothing
        throw(ArgumentError("CrossSectionCache fingerprint mismatch for $(process): w0cdf grid field $(grid_mismatch) differs from the current average_scattering_rate cache strategy; rebuild with build_w0cdf_pchip_cache using the current sigma_cutoff/default grid settings"))
    end
    return cache
end

function CrossSectionCache(process::Symbol, s_vals::Vector{Float64}, sigma_vals::Vector{Float64})
    length(s_vals) == length(sigma_vals) || error("CrossSectionCache: s_vals and sigma_vals length mismatch")
    cache = CrossSectionCache(process, CrossSectionData(s_vals, sigma_vals), CachedInterpolation(), AsymptoticConfig())
    _ensure_pchip_slopes!(cache)
    return cache
end

function insert_sigma!(cache::CrossSectionCache, s::Float64, σ::Float64)
    idx = searchsortedfirst(cache.s_vals, s)
    insert!(cache.s_vals, idx, s)
    insert!(cache.sigma_vals, idx, σ)
    cache.pchip_dirty = true
    cache.peak_dirty = true
    cache.local_dirty = true
end

function _ensure_peak_profile!(cache::CrossSectionCache)
    if !cache.peak_dirty
        return
    end
    n = length(cache.sigma_vals)
    if n == 0
        cache.peak_ratio = 0.0
        cache.peak_s = 0.0
        cache.peak_dirty = false
        return
    end
    idx = argmax(cache.sigma_vals)
    μ = mean(cache.sigma_vals)
    cache.peak_ratio = (μ > 0.0 && isfinite(μ)) ? (cache.sigma_vals[idx] / μ) : 0.0
    cache.peak_s = cache.s_vals[idx]
    cache.peak_dirty = false
end

@inline function _signmatch(a::Float64, b::Float64)
    return (a == 0.0 && b == 0.0) || (a > 0.0 && b > 0.0) || (a < 0.0 && b < 0.0)
end

function _ensure_pchip_slopes!(cache::CrossSectionCache)
    if !cache.pchip_dirty && length(cache.pchip_slopes) == length(cache.s_vals)
        return
    end
    n = length(cache.s_vals)
    cache.pchip_slopes = zeros(Float64, n)
    cache.pchip_dirty = false
    n <= 1 && return

    x = cache.s_vals
    y = cache.sigma_vals
    h = [x[i+1] - x[i] for i in 1:(n-1)]
    d = [(y[i+1] - y[i]) / h[i] for i in 1:(n-1)]

    if n == 2
        cache.pchip_slopes[1] = d[1]
        cache.pchip_slopes[2] = d[1]
        return
    end

    # Endpoints (Fritsch–Carlson)
    h1, h2 = h[1], h[2]
    d1, d2 = d[1], d[2]
    m1 = ((2h1 + h2) * d1 - h1 * d2) / (h1 + h2)
    if !_signmatch(m1, d1)
        m1 = 0.0
    elseif (!_signmatch(d1, d2)) && (abs(m1) > abs(3d1))
        m1 = 3d1
    end
    cache.pchip_slopes[1] = m1

    hn1, hn2 = h[n-2], h[n-1]
    dn1, dn2 = d[n-2], d[n-1]
    mn = ((2hn2 + hn1) * dn2 - hn2 * dn1) / (hn1 + hn2)
    if !_signmatch(mn, dn2)
        mn = 0.0
    elseif (!_signmatch(dn2, dn1)) && (abs(mn) > abs(3dn2))
        mn = 3dn2
    end
    cache.pchip_slopes[n] = mn

    # Interior slopes
    for i in 2:(n-1)
        if d[i-1] == 0.0 || d[i] == 0.0 || !_signmatch(d[i-1], d[i])
            cache.pchip_slopes[i] = 0.0
            continue
        end
        w1 = 2h[i] + h[i-1]
        w2 = h[i] + 2h[i-1]
        cache.pchip_slopes[i] = (w1 + w2) / (w1 / d[i-1] + w2 / d[i])
    end
end

@inline function _pchip_eval(x1::Float64, x2::Float64, y1::Float64, y2::Float64, m1::Float64, m2::Float64, x::Float64)
    h = x2 - x1
    t = (x - x1) / h
    t2 = t * t
    t3 = t2 * t
    h00 = 2t3 - 3t2 + 1
    h10 = t3 - 2t2 + t
    h01 = -2t3 + 3t2
    h11 = t3 - t2
    return h00 * y1 + h10 * h * m1 + h01 * y2 + h11 * h * m2
end

function precompute_cross_section!(cache::CrossSectionCache, s_grid::Vector{Float64},
    quark_params::Union{NamedTuple, QuarkParams}, thermo_params::Union{NamedTuple, ThermoParams}, K_coeffs::NamedTuple;
    n_points::Int=TotalCrossSection.DEFAULT_T_INTEGRAL_POINTS,
    threshold_subtraction::Bool=false,
    asym_window::Float64=0.6,
    asym_fit_min_points::Int=8,
    asym_extra_points::Int=10)
    quark_nt = normalize_quark_input(quark_params)
    thermo_nt = normalize_thermo_input(thermo_params)
    return _precompute_cross_section_core!(cache, s_grid, quark_nt, thermo_nt, K_coeffs;
        n_points=n_points,
        threshold_subtraction=threshold_subtraction,
        asym_window=asym_window,
        asym_fit_min_points=asym_fit_min_points,
        asym_extra_points=asym_extra_points)
end

function _precompute_cross_section_core!(cache::CrossSectionCache, s_grid::Vector{Float64},
    quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple;
    n_points::Int=TotalCrossSection.DEFAULT_T_INTEGRAL_POINTS,
    threshold_subtraction::Bool=false,
    asym_window::Float64=0.6,
    asym_fit_min_points::Int=8,
    asym_extra_points::Int=10,
    fingerprint_grid_metadata::Union{Nothing,NamedTuple}=nothing)
    # compute raw σ(s) for grid
    raw = Float64[]
    # prepare containers for potential extra samples (defined regardless of threshold_subtraction)
    extra_s = Float64[]
    extra_raw = Float64[]
    for s in s_grid
        σ = total_cross_section(cache.process, s, quark_params, thermo_params, K_coeffs; n_points=n_points)
        push!(raw, σ)
    end

    # optional threshold asymptotic subtraction
    cache.asym_enabled = false
    cache.asym_requested = threshold_subtraction
    if threshold_subtraction
        # compute s_th from process masses
        pi_sym, pj_sym, pc_sym, pd_sym = parse_particles_from_process(cache.process)
        mi = get_mass(pi_sym, quark_params)
        mj = get_mass(pj_sym, quark_params)
        mc = get_mass(pc_sym, quark_params)
        md = get_mass(pd_sym, quark_params)
        s_th = max((mi + mj)^2, (mc + md)^2)

        # select points near threshold for fitting
        idxs = findall(i -> isfinite(raw[i]) && raw[i] > 0.0 && (s_grid[i] - s_th) > 1e-12 && (s_grid[i] - s_th) <= asym_window, 1:length(s_grid))

        # Prepare containers for possible extra samples
        extra_s = Float64[]
        extra_raw = Float64[]

        # If insufficient points, compute extra high-resolution samples near threshold
        if length(idxs) < asym_fit_min_points
            extra_s = collect(range(s_th + 1e-12, stop = s_th + asym_window, length = asym_extra_points))
            for s_ex in extra_s
                σ_ex = total_cross_section(cache.process, s_ex, quark_params, thermo_params, K_coeffs; n_points=n_points)
                push!(extra_raw, σ_ex)
            end
        end

        # combine original grid and extra samples for fitting and (later) for caching
        s_all = vcat(s_grid, extra_s)
        raw_all = vcat(raw, extra_raw)
        idxs_all = findall(i -> isfinite(raw_all[i]) && raw_all[i] > 0.0 && (s_all[i] - s_th) > 1e-12 && (s_all[i] - s_th) <= asym_window, 1:length(s_all))
        if length(idxs_all) >= asym_fit_min_points
            vals = [raw_all[i] * sqrt(s_all[i] - s_th) for i in idxs_all]
            A_est = median(vals)
            if isfinite(A_est) && A_est > 0.0
                cache.asym_enabled = true
                cache.asym_s0 = s_th
                cache.asym_A = A_est
            end
        end
    end

    # store sigma_vals as regularized (raw - asym) if asym enabled
    # If extra samples were taken, include them in the cached grid so interpolation benefits
    if !isempty(extra_s)
        # Build mapping from s -> raw (prefer original grid values if duplicated)
        raw_map = Dict{Float64, Float64}()
        for (s, σ) in zip(s_grid, raw)
            raw_map[s] = σ
        end
        for (s, σ) in zip(extra_s, extra_raw)
            # only set if not already present (avoid overwriting original grid values)
            if !haskey(raw_map, s)
                raw_map[s] = σ
            end
        end
        s_used = sort(collect(keys(raw_map)))
        cache.s_vals = copy(s_used)
        sigma_list = Float64[]
        for s in s_used
            σ = raw_map[s]
            if cache.asym_enabled
                σ_asym = (s > cache.asym_s0) ? cache.asym_A / sqrt(max(1e-16, s - cache.asym_s0)) : 0.0
                push!(sigma_list, max(0.0, σ - σ_asym))
            else
                push!(sigma_list, σ)
            end
        end
        cache.sigma_vals = sigma_list
    else
        cache.s_vals = copy(s_grid)
        if cache.asym_enabled
            reg = Float64[]
            for (s, σ) in zip(s_grid, raw)
                σ_asym = (s > cache.asym_s0) ? cache.asym_A / sqrt(max(1e-16, s - cache.asym_s0)) : 0.0
                push!(reg, max(0.0, σ - σ_asym))
            end
            cache.sigma_vals = reg
        else
            cache.sigma_vals = copy(raw)
        end
    end
    cache.pchip_dirty = true
    cache.peak_dirty = true
    cache.local_dirty = true
    _ensure_pchip_slopes!(cache)
    grid_metadata = fingerprint_grid_metadata === nothing ? _explicit_grid_metadata(s_grid) : fingerprint_grid_metadata
    cache.fingerprint = _cross_section_cache_fingerprint(
        cache.process,
        quark_params,
        thermo_params,
        K_coeffs;
        n_points=n_points,
        threshold_subtraction=threshold_subtraction,
        asym_window=asym_window,
        asym_fit_min_points=asym_fit_min_points,
        asym_extra_points=asym_extra_points,
        s_grid=cache.s_vals,
        grid_metadata=grid_metadata,
    )
    return cache
end

function interpolate_sigma(cache::CrossSectionCache, s::Float64)
    n = length(cache.s_vals)
    if n == 0
        return nothing
    elseif n == 1
        return cache.sigma_vals[1]
    elseif s == cache.s_vals[1]
        return cache.sigma_vals[1]
    elseif s == cache.s_vals[end]
        return cache.sigma_vals[end]
    elseif s < cache.s_vals[1] || s > cache.s_vals[end]
        return nothing
    end

    idx = searchsortedfirst(cache.s_vals, s)
    s1, s2 = cache.s_vals[idx-1], cache.s_vals[idx]
    σ1, σ2 = cache.sigma_vals[idx-1], cache.sigma_vals[idx]

    # Standardized interpolation: :pchip
    _ensure_pchip_slopes!(cache)
    m1 = cache.pchip_slopes[idx-1]
    m2 = cache.pchip_slopes[idx]
    y = _pchip_eval(s1, s2, σ1, σ2, m1, m2, s)
    return isfinite(y) ? max(0.0, y) : 0.0
end

@inline function interpolate_sigma_linear(cache::CrossSectionCache, s::Float64)
    n = length(cache.s_vals)
    if n == 0
        return nothing
    elseif n == 1
        return cache.sigma_vals[1]
    elseif s == cache.s_vals[1]
        return cache.sigma_vals[1]
    elseif s == cache.s_vals[end]
        return cache.sigma_vals[end]
    elseif s < cache.s_vals[1] || s > cache.s_vals[end]
        return nothing
    end

    idx = searchsortedfirst(cache.s_vals, s)
    s1, s2 = cache.s_vals[idx-1], cache.s_vals[idx]
    y1, y2 = cache.sigma_vals[idx-1], cache.sigma_vals[idx]
    t = (s - s1) / (s2 - s1)
    y = y1 + t * (y2 - y1)
    return isfinite(y) ? max(0.0, y) : 0.0
end

function get_sigma(cache::CrossSectionCache, s::Float64,
    quark_params::Union{NamedTuple, QuarkParams}, thermo_params::Union{NamedTuple, ThermoParams}, K_coeffs::NamedTuple;
    n_points::Int=TotalCrossSection.DEFAULT_T_INTEGRAL_POINTS)
    quark_nt = normalize_quark_input(quark_params)
    thermo_nt = normalize_thermo_input(thermo_params)
    return _get_sigma_core(cache, s, quark_nt, thermo_nt, K_coeffs; n_points=n_points)
end

function _get_sigma_core(cache::CrossSectionCache, s::Float64,
    quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple;
    n_points::Int=TotalCrossSection.DEFAULT_T_INTEGRAL_POINTS,
    interpolation_mode::Symbol=:pchip)
    if interpolation_mode == :direct
        return total_cross_section(cache.process, s, quark_params, thermo_params, K_coeffs; n_points=n_points)
    end

    # Only cached PCHIP interpolation is supported.
    n = length(cache.s_vals)
    n == 0 && error("CrossSectionCache has no points; precompute σ(s) first")
    if n == 1
        return cache.sigma_vals[1]
    elseif s < cache.s_vals[1] || s > cache.s_vals[end]
        # outside cached window: return only asymptotic part if available
        if cache.asym_enabled && s > cache.asym_s0
            return cache.asym_A / sqrt(max(1e-16, s - cache.asym_s0))
        end
        return 0.0
    elseif s == cache.s_vals[1]
        return cache.sigma_vals[1]
    elseif s == cache.s_vals[end]
        return cache.sigma_vals[end]
    end

    if interpolation_mode == :hybrid_threshold
        _ensure_peak_profile!(cache)
        if cache.process in (:uubar_to_uubar, :uubar_to_ddbar)
            peak_offset = cache.peak_s - cache.asym_s0
            if cache.asym_enabled && cache.peak_ratio >= 20.0 && peak_offset >= 0.0 && peak_offset <= 0.8 && s <= (cache.asym_s0 + 0.8)
                # Lightweight error fuse: evaluate two nearby quadratures.
                n_lo = max(n_points, 8)
                n_hi = max(n_lo + 4, 12)
                σ_lo = total_cross_section(cache.process, s, quark_params, thermo_params, K_coeffs; n_points=n_lo)
                σ_hi = total_cross_section(cache.process, s, quark_params, thermo_params, K_coeffs; n_points=n_hi)
                rel = abs(σ_hi - σ_lo) / max(abs(σ_hi), 1e-12)
                if rel > 0.05
                    n_ref = max(n_hi + 4, 16)
                    return total_cross_section(cache.process, s, quark_params, thermo_params, K_coeffs; n_points=n_ref)
                end
                return σ_hi
            end
        end
    end

    val = if interpolation_mode == :linear
        interpolate_sigma_linear(cache, s)
    elseif interpolation_mode == :hybrid_threshold
        interpolate_sigma_linear(cache, s)
    else
        interpolate_sigma(cache, s)
    end
    val === nothing && error("interpolation failed inside cache window")
    # add back analytic asymptotic if enabled
    if cache.asym_enabled && s > cache.asym_s0
        return val + cache.asym_A / sqrt(max(1e-16, s - cache.asym_s0))
    end
    return val
end

# -------------------- w0cdf σ-grid design (internal) --------------------
function _build_semi_infinite_p_grid(p_nodes::Int, scale::Float64)
    t_grid, t_w = gauleg(0.0, 1.0, p_nodes)
    p_vals = Float64[]
    p_wts = Float64[]
    dp_jac = Float64[]
    for (t, wt) in zip(t_grid, t_w)
        t >= 0.9999 && continue
        p = scale * t / (1.0 - t)
        dp_dt = scale / (1.0 - t)^2
        push!(p_vals, p)
        push!(p_wts, wt)
        push!(dp_jac, dp_dt)
    end
    return p_vals, p_wts, dp_jac
end

function _build_finite_cutoff_p_grid(p_nodes::Int, p_max::Float64)
    p_vals, p_wts = gauleg(0.0, p_max, p_nodes)
    dp_jac = ones(Float64, length(p_vals))
    return p_vals, p_wts, dp_jac
end

"""
    design_w0cdf_s_grid(process, quark_params, thermo_params; kwargs...)

基于 ω 积分权重的 CDF 设计 σ(s) 采样网格。

# Arguments
- `process::Symbol`: 散射过程
- `quark_params`: 夸克参数，可以是 `QuarkParams` 结构体或 NamedTuple (m, μ)
- `thermo_params`: 热力学参数，可以是 `ThermoParams` 结构体或 NamedTuple (T, Φ, Φbar, ξ)

# Keyword Arguments
- `N::Int`: 采样点数（默认 60）
- `p_nodes::Int`: 动量积分节点数（默认 14）
- `angle_nodes::Int`: 角度积分节点数（默认 4）
- `phi_nodes::Int`: 方位角积分节点数（默认 8）
- `p_cutoff::Union{Nothing,Float64}`: 动量截断（fm⁻¹）。
  - `nothing`（默认）：使用半无穷积分 [0, ∞)
  - 指定值（如 `Λ_inv_fm`）：使用有限截断 [0, p_cutoff]
- `scale::Float64`: 半无穷积分的尺度参数（默认 10.0，仅当 p_cutoff=nothing 时使用）

# Returns
- `s_grid::Vector{Float64}`: σ(s) 采样点的 s 值数组
"""
function design_w0cdf_s_grid(
    process::Symbol,
    quark_params::Union{NamedTuple, QuarkParams},
    thermo_params::Union{NamedTuple, ThermoParams};
    N::Int=DEFAULT_SIGMA_GRID_N,
    p_nodes::Int=DEFAULT_W0CDF_P_NODES,
    angle_nodes::Int=DEFAULT_W0CDF_ANGLE_NODES,
    phi_nodes::Int=DEFAULT_W0CDF_PHI_NODES,
    p_cutoff::Union{Nothing,Float64}=nothing,
    scale::Float64=DEFAULT_SEMI_INF_SCALE,
)
    quark_nt = normalize_quark_input(quark_params)
    thermo_nt = normalize_thermo_input(thermo_params)
    return _design_w0cdf_s_grid_core(process, quark_nt, thermo_nt;
        N=N,
        p_nodes=p_nodes,
        angle_nodes=angle_nodes,
        phi_nodes=phi_nodes,
        p_cutoff=p_cutoff,
        scale=scale)
end

function _design_w0cdf_s_grid_core(
    process::Symbol,
    quark_params::NamedTuple,
    thermo_params::NamedTuple;
    N::Int=DEFAULT_SIGMA_GRID_N,
    p_nodes::Int=DEFAULT_W0CDF_P_NODES,
    angle_nodes::Int=DEFAULT_W0CDF_ANGLE_NODES,
    phi_nodes::Int=DEFAULT_W0CDF_PHI_NODES,
    p_cutoff::Union{Nothing,Float64}=nothing,
    scale::Float64=DEFAULT_SEMI_INF_SCALE,
)
    pi_sym, pj_sym, pc_sym, pd_sym = parse_particles_from_process(process)
    mi = get_mass(pi_sym, quark_params)
    mj = get_mass(pj_sym, quark_params)
    μi = get_mu(pi_sym, quark_params)
    μj = get_mu(pj_sym, quark_params)
    mc = get_mass(pc_sym, quark_params)
    md = get_mass(pd_sym, quark_params)

    T = thermo_params.T
    Φ = thermo_params.Φ
    Φbar = thermo_params.Φbar
    ξ = hasproperty(thermo_params, :ξ) ? thermo_params.ξ : 0.0

    # 根据 p_cutoff 选择动量网格构建方式
    p_vals, p_wts, dp_jac = if p_cutoff !== nothing
        _build_finite_cutoff_p_grid(p_nodes, p_cutoff)
    else
        _build_semi_infinite_p_grid(p_nodes, scale)
    end
    
    cos_grid, cos_w = gauleg(-1.0, 1.0, angle_nodes)
    phi_grid, phi_w = gauleg(0.0, TWO_PI, phi_nodes)

    sqrt_s_samples = Float64[]
    weights = Float64[]

    s_bo = max((mi + mj)^2, (mc + md)^2)
    # 如果指定了 p_cutoff，则限制 s 的上限
    s_up = if p_cutoff !== nothing
        min((sqrt(mi^2 + p_cutoff^2) + sqrt(mj^2 + p_cutoff^2))^2,
            (sqrt(mc^2 + p_cutoff^2) + sqrt(md^2 + p_cutoff^2))^2)
    else
        Inf
    end

    for (p_i, w_pi, dp_i) in zip(p_vals, p_wts, dp_jac)
        Ei = energy_from_p(p_i, mi)
        for (p_j, w_pj, dp_j) in zip(p_vals, p_wts, dp_jac)
            Ej = energy_from_p(p_j, mj)
            for (cθi, w_cθi) in zip(cos_grid, cos_w)
                sθi = sqrt(max(1.0 - cθi * cθi, 0.0))
                f_i = distribution_with_anisotropy(pi_sym, p_i, mi, μi, T, Φ, Φbar, ξ, cθi)
                f_i == 0.0 && continue
                for (cθj, w_cθj) in zip(cos_grid, cos_w)
                    sθj = sqrt(max(1.0 - cθj * cθj, 0.0))
                    f_j = distribution_with_anisotropy(pj_sym, p_j, mj, μj, T, Φ, Φbar, ξ, cθj)
                    f_j == 0.0 && continue
                    for (φ, wφ) in zip(phi_grid, phi_w)
                        cosΘ = cθi * cθj + sθi * sθj * cos(φ)
                        s = mi^2 + mj^2 + 2.0 * (Ei * Ej - p_i * p_j * cosΘ)
                        s <= s_bo && continue
                        # 如果指定了 p_cutoff，跳过超出 s_up 的点
                        s >= s_up && continue

                        s_rt = sqrt(s)
                        Ei_cm = (s + mi^2 - mj^2) / (2.0 * s_rt)
                        Ej_cm = (s - mi^2 + mj^2) / (2.0 * s_rt)
                        pi_cm = sqrt(max(0.0, (s - (mi + mj)^2) * (s - (mi - mj)^2))) / (2.0 * s_rt)
                        v_rel_num = (Ei_cm * Ej_cm + pi_cm * pi_cm)^2 - (mi * mj)^2
                        v_rel = v_rel_num > 0.0 ? sqrt(v_rel_num) / (Ei_cm * Ej_cm) : 0.0
                        (v_rel == 0.0 || v_rel > 2.0) && continue

                        w0 = w_pi * w_pj * w_cθi * w_cθj * wφ * (p_i^2) * (p_j^2) * f_i * f_j * v_rel * dp_i * dp_j
                        (isfinite(w0) && w0 > 0.0) || continue
                        push!(sqrt_s_samples, s_rt)
                        push!(weights, w0)
                    end
                end
            end
        end
    end

    isempty(weights) && error("w0cdf design produced empty weights for process $process")

    order = sortperm(sqrt_s_samples)
    xs = sqrt_s_samples[order]
    ws = weights[order]
    cdf = cumsum(ws)
    tot = cdf[end]
    tot <= 0.0 && error("w0cdf design produced non-positive total weight")

    # Keep exactly N points, but guarantee endpoints for safe clamping.
    grid = Vector{Float64}(undef, N)
    grid[1] = xs[1]
    for i in 2:(N-1)
        q = (i - 0.5) / N * tot
        idx = searchsortedfirst(cdf, q)
        idx = clamp(idx, 1, length(xs))
        grid[i] = xs[idx]
    end
    grid[end] = xs[end]

    s_grid = grid .^ 2
    return s_grid
end

"""
    build_w0cdf_pchip_cache(process, quark_params, thermo_params, K_coeffs; kwargs...)

构建基于 w0cdf 设计的 σ(s) 缓存。

# Arguments
- `process::Symbol`: 散射过程
- `quark_params`: 夸克参数，可以是 `QuarkParams` 结构体或 NamedTuple (m, μ)
- `thermo_params`: 热力学参数，可以是 `ThermoParams` 结构体或 NamedTuple (T, Φ, Φbar, ξ)
- `K_coeffs::NamedTuple`: 有效耦合系数

# Keyword Arguments
- `N::Int`: 采样点数（默认 60）
- `design_p_nodes::Int`: w0cdf 设计时的动量节点数（默认 14）
- `design_angle_nodes::Int`: w0cdf 设计时的角度节点数（默认 4）
- `design_phi_nodes::Int`: w0cdf 设计时的方位角节点数（默认 8）
- `p_cutoff::Union{Nothing,Float64}`: 动量截断（fm⁻¹）。
  - `nothing`：使用半无穷积分 [0, ∞)
  - 指定值（如 `Λ_inv_fm`）：使用有限截断 [0, p_cutoff]，**推荐用于生产**
- `scale::Float64`: 半无穷积分的尺度参数（默认 10.0，仅当 p_cutoff=nothing 时使用）
- `n_sigma_points::Int`: σ(s) 计算时的 t 积分点数
"""
function build_w0cdf_pchip_cache(
    process::Symbol,
    quark_params::Union{NamedTuple, QuarkParams},
    thermo_params::Union{NamedTuple, ThermoParams},
    K_coeffs::NamedTuple;
    N::Int=DEFAULT_SIGMA_GRID_N,
    design_p_nodes::Int=DEFAULT_W0CDF_P_NODES,
    design_angle_nodes::Int=DEFAULT_W0CDF_ANGLE_NODES,
    design_phi_nodes::Int=DEFAULT_W0CDF_PHI_NODES,
    p_cutoff::Union{Nothing,Float64}=nothing,
    scale::Float64=DEFAULT_SEMI_INF_SCALE,
    n_sigma_points::Int=TotalCrossSection.DEFAULT_T_INTEGRAL_POINTS,
    threshold_subtraction::Bool=false,
    asym_window::Float64=0.6,
    asym_fit_min_points::Int=8,
    asym_extra_points::Int=10,
)
    quark_nt = normalize_quark_input(quark_params)
    thermo_nt = normalize_thermo_input(thermo_params)
    s_grid = _design_w0cdf_s_grid_core(
        process,
        quark_nt,
        thermo_nt;
        N=N,
        p_nodes=design_p_nodes,
        angle_nodes=design_angle_nodes,
        phi_nodes=design_phi_nodes,
        p_cutoff=p_cutoff,
        scale=scale,
    )
    cache = CrossSectionCache(process)
    _precompute_cross_section_core!(cache, s_grid, quark_nt, thermo_nt, K_coeffs;
        n_points=n_sigma_points,
        threshold_subtraction=threshold_subtraction,
        asym_window=asym_window,
        asym_fit_min_points=asym_fit_min_points,
        asym_extra_points=asym_extra_points,
        fingerprint_grid_metadata=_w0cdf_grid_metadata(
            s_grid;
            N=N,
            design_p_nodes=design_p_nodes,
            design_angle_nodes=design_angle_nodes,
            design_phi_nodes=design_phi_nodes,
            p_cutoff=p_cutoff,
            scale=scale,
        ))
    _ensure_pchip_slopes!(cache)
    return cache
end

# -------------------- ρ 计算（各向异性） --------------------
# 半无穷积分的默认参数
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
- `angle_nodes::Int`: 角度积分节点数 (默认2)
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

# -------------------- 平均散射率主函数 --------------------

@inline function _resolve_auto_threshold_subtraction(
    threshold_subtraction::Bool,
    mi::Float64,
    mj::Float64,
    mc::Float64,
    md::Float64,
)::Bool
    lhs = mi^2 + mj^2
    rhs = mc^2 + md^2
    scale = max(abs(lhs), abs(rhs), 1.0)
    tol = 64 * eps(scale)
    if lhs > rhs + tol
        return true
    end
    return threshold_subtraction
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

@inline function _ensure_default_angle_grids(
    cos_grid::Union{Nothing,Vector{Float64}},
    cos_w::Union{Nothing,Vector{Float64}},
    phi_grid::Union{Nothing,Vector{Float64}},
    phi_w::Union{Nothing,Vector{Float64}},
    angle_nodes::Int,
    phi_nodes::Int,
)
    cos_grid === nothing && ((cos_grid, cos_w) = _cached_interval_grid(-1.0, 1.0, angle_nodes))
    phi_grid === nothing && ((phi_grid, phi_w) = _cached_interval_grid(0.0, TWO_PI, phi_nodes))
    return cos_grid, cos_w, phi_grid, phi_w
end

"""
    average_scattering_rate(process, quark_params, thermo_params, K_coeffs; kwargs...)

计算平均散射率（采用半无穷积分 [0, ∞) 的实现）。

# Arguments
- `process::Symbol`: 散射过程 (如 :ud_to_ud)
- `quark_params`: 夸克参数，可以是 `QuarkParams` 结构体或 NamedTuple (m, μ, A)
- `thermo_params`: 热力学参数，可以是 `ThermoParams` 结构体或 NamedTuple (T, Φ, Φbar, ξ)
- `K_coeffs::NamedTuple`: 有效耦合系数

# Keyword Arguments
- `p_nodes::Int`: 动量积分节点数 (默认6)
- `angle_nodes::Int`: 角度积分节点数 (默认2)
- `phi_nodes::Int`: 方位角积分节点数 (默认4)
- `n_sigma_points::Int`: σ(s) 计算时的 t 积分点数
- `require_cache_fingerprint::Bool`: 为 `true` 时拒绝无指纹的外部 σ(s) 缓存
- `scale::Float64`: 半无穷积分尺度参数
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
    cs_cache::Union{Nothing,CrossSectionCache}=nothing,
    n_sigma_points::Int=TotalCrossSection.DEFAULT_T_INTEGRAL_POINTS,
    scale::Float64=DEFAULT_SEMI_INF_SCALE,
    density_p_grid::Union{Nothing,Vector{Float64}}=nothing,
    density_p_w::Union{Nothing,Vector{Float64}}=nothing,
    density_p_nodes::Int=DEFAULT_SEMI_INF_NODES,
    density_scale::Float64=DEFAULT_SEMI_INF_SCALE,
    apply_s_domain_cut::Bool=true,
    sigma_cutoff::Union{Nothing,Float64}=nothing,  # 新增：σ(s)有效范围的动量截断
    # 允许在未传入 cs_cache 时自动构建并配置阈值减法
    threshold_subtraction::Bool=false,
    asym_window::Float64=0.6,
    asym_fit_min_points::Int=8,
    asym_extra_points::Int=10,
    interpolation_mode::Symbol=:pchip,
    require_cache_fingerprint::Bool=false,
    band_edges::Union{Nothing,Vector{Float64}}=nothing,
    band_omega_out::Union{Nothing,Base.RefValue{Vector{Float64}}}=nothing,
    band_omega_sigma_out::Union{Nothing,Base.RefValue{Vector{Float64}}}=nothing,
)::Float64
    quark_params = normalize_quark_input(quark_params)
    thermo_params = normalize_thermo_input(thermo_params)
    # 解析粒子、质量、化学势
    pi_sym, pj_sym, pc_sym, pd_sym = parse_particles_from_process(process)
    mi = get_mass(pi_sym, quark_params); mj = get_mass(pj_sym, quark_params)
    mc = get_mass(pc_sym, quark_params); md = get_mass(pd_sym, quark_params)
    μi = get_mu(pi_sym, quark_params); μj = get_mu(pj_sym, quark_params)
    μc = get_mu(pc_sym, quark_params); μd = get_mu(pd_sym, quark_params)

    T = thermo_params.T; Φ = thermo_params.Φ; Φbar = thermo_params.Φbar
    ξ = hasproperty(thermo_params, :ξ) ? thermo_params.ξ : 0.0

    # Build the finalized σ-cache strategy by default.
    # 如果指定了 sigma_cutoff，则使用有限截断的 w0cdf 设计
    thr_for_build = _resolve_auto_threshold_subtraction(threshold_subtraction, mi, mj, mc, md)
    if cs_cache === nothing
        cs_cache = build_w0cdf_pchip_cache(
            process,
            quark_params,
            thermo_params,
            K_coeffs;
            N=DEFAULT_SIGMA_GRID_N,
            design_p_nodes=DEFAULT_W0CDF_P_NODES,
            design_angle_nodes=DEFAULT_W0CDF_ANGLE_NODES,
            design_phi_nodes=DEFAULT_W0CDF_PHI_NODES,
            p_cutoff=sigma_cutoff,  # 传递 sigma_cutoff 作为 p_cutoff
            scale=scale,
            n_sigma_points=n_sigma_points,
            threshold_subtraction=thr_for_build,
            asym_window=asym_window,
            asym_fit_min_points=asym_fit_min_points,
            asym_extra_points=asym_extra_points,
        )
    elseif interpolation_mode !== :direct
        _validate_cross_section_cache!(
            cs_cache,
            process,
            quark_params,
            thermo_params,
            K_coeffs;
            n_points=n_sigma_points,
            threshold_subtraction=thr_for_build,
            asym_window=asym_window,
            asym_fit_min_points=asym_fit_min_points,
            asym_extra_points=asym_extra_points,
            sigma_cutoff=sigma_cutoff,
            scale=scale,
            require_cache_fingerprint=require_cache_fingerprint,
        )
    end

    # 使用半无穷积分 [0, ∞)
    return _average_scattering_rate_semi_infinite(
        process, pi_sym, pj_sym, mi, mj, μi, μj, T, Φ, Φbar, ξ,
        quark_params, thermo_params, K_coeffs,
        p_nodes, angle_nodes, phi_nodes, scale,
        p_grid, p_w, cos_grid, cos_w, phi_grid, phi_w,
        cs_cache, n_sigma_points,
        density_p_grid, density_p_w, density_p_nodes, density_scale,
        mc, md, apply_s_domain_cut, sigma_cutoff, interpolation_mode,
        band_edges, band_omega_out, band_omega_sigma_out,
    )
end

# 半无穷积分版本（保留原有实现）
function _average_scattering_rate_semi_infinite(
    process, pi_sym, pj_sym, mi, mj, μi, μj, T, Φ, Φbar, ξ,
    quark_params, thermo_params, K_coeffs,
    p_nodes, angle_nodes, phi_nodes, scale,
    p_grid, p_w, cos_grid, cos_w, phi_grid, phi_w,
    cs_cache, n_sigma_points,
    density_p_grid, density_p_w, density_p_nodes, density_scale,
    mc, md, apply_s_domain_cut, sigma_cutoff,
    interpolation_mode::Symbol,
    band_edges::Union{Nothing,Vector{Float64}},
    band_omega_out::Union{Nothing,Base.RefValue{Vector{Float64}}},
    band_omega_sigma_out::Union{Nothing,Base.RefValue{Vector{Float64}}},
)
    validate_grid_weight_pair("average_scattering_rate", "p_grid", p_grid, "p_w", p_w)
    validate_grid_weight_pair("average_scattering_rate", "cos_grid", cos_grid, "cos_w", cos_w)
    validate_grid_weight_pair("average_scattering_rate", "phi_grid", phi_grid, "phi_w", phi_w)

    # Momentum integration:
    # - Default: semi-infinite [0, ∞) via t∈[0,1) mapping.
    # - If p_grid/p_w are provided: integrate directly on that finite grid.
    if p_grid !== nothing && p_w !== nothing
        p_vals = p_grid
        quadrature_wts = p_w
    else
        p_vals, quadrature_wts = _cached_semi_infinite_momentum_grid(p_nodes, scale)
    end

    # 使用 sigma_cutoff 参数（如果提供）来确定 s 范围，否则从动量网格推断。
    Λ, s_bo, s_up = _resolve_sigma_support_bounds(mi, mj, mc, md, sigma_cutoff, p_grid, p_w)
    cos_grid, cos_w, phi_grid, phi_w = _ensure_default_angle_grids(cos_grid, cos_w, phi_grid, phi_w, angle_nodes, phi_nodes)

    # 数密度（用于归一化）
    # - 默认：半无穷积分 [0,∞)
    # - 若提供 density_p_grid/density_p_w：则直接使用该有限动量网格
    ρ_i = number_density(pi_sym, mi, μi, T, Φ, Φbar, ξ;
        p_nodes=density_p_nodes, angle_nodes=angle_nodes,
        p_grid=density_p_grid, p_w=density_p_w,
        scale=density_scale)
    ρ_j = number_density(pj_sym, mj, μj, T, Φ, Φbar, ξ;
        p_nodes=density_p_nodes, angle_nodes=angle_nodes,
        p_grid=density_p_grid, p_w=density_p_w,
        scale=density_scale)
    if ρ_i == 0.0 || ρ_j == 0.0
        return 0.0
    end

    # 原始公式的前因子: DQ²/(2π)⁵ = DQ²/(32π⁵)
    prefactor = (DQ^2) / (32.0 * π^5 * ρ_i * ρ_j)
    # 当提供了 sigma_cutoff 或有限动量网格时，应用 s 范围截断
    apply_s_cut = apply_s_domain_cut && (!isnan(Λ))
    s_th = max((mi + mj)^2, (mc + md)^2)

    ω = _omega_integral_5d(
        process, pi_sym, pj_sym,
        mi, mj, μi, μj, T, Φ, Φbar, ξ,
        quark_params, thermo_params, K_coeffs,
        p_vals, quadrature_wts,
        cos_grid, cos_w,
        phi_grid, phi_w,
        cs_cache, n_sigma_points,
        apply_s_cut, s_bo, s_up,
        interpolation_mode,
        s_th,
        band_edges,
        band_omega_out,
        band_omega_sigma_out,
    )

    return prefactor * ω
end

@inline function _band_index(ds::Float64, edges::Vector{Float64})
    n = length(edges)
    if n < 2
        return 0
    end
    if ds < edges[1]
        return 0
    end
    bi = searchsortedlast(edges, ds)
    bi <= 0 && return 0
    if bi >= n
        return n - 1
    end
    return bi
end

function _omega_integral_5d(
    process::Symbol,
    pi_sym::Symbol,
    pj_sym::Symbol,
    mi::Float64,
    mj::Float64,
    μi::Float64,
    μj::Float64,
    T::Float64,
    Φ::Float64,
    Φbar::Float64,
    ξ::Float64,
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    K_coeffs::NamedTuple,
    p_vals::Vector{Float64},
    quadrature_wts::Vector{Float64},
    cos_grid::Vector{Float64},
    cos_w::Vector{Float64},
    phi_grid::Vector{Float64},
    phi_w::Vector{Float64},
    cs_cache::CrossSectionCache,
    n_sigma_points::Int,
    apply_s_cut::Bool,
    s_bo::Float64,
    s_up::Float64,
    interpolation_mode::Symbol,
    s_th::Float64,
    band_edges::Union{Nothing,Vector{Float64}},
    band_omega_out::Union{Nothing,Base.RefValue{Vector{Float64}}},
    band_omega_sigma_out::Union{Nothing,Base.RefValue{Vector{Float64}}},
)::Float64
    ω = 0.0
    band_omega = nothing
    band_omega_sigma = nothing
    if band_edges !== nothing && length(band_edges) >= 2
        nbin = length(band_edges) - 1
        band_omega = zeros(Float64, nbin)
        band_omega_sigma = zeros(Float64, nbin)
    end

    for (p_i, w_pi) in zip(p_vals, quadrature_wts)
        Ei = energy_from_p(p_i, mi)

        for (p_j, w_pj) in zip(p_vals, quadrature_wts)
            Ej = energy_from_p(p_j, mj)

            for (cθi, w_cθi) in zip(cos_grid, cos_w)
                sθi = sqrt(max(1.0 - cθi * cθi, 0.0))
                for (cθj, w_cθj) in zip(cos_grid, cos_w)
                    sθj = sqrt(max(1.0 - cθj * cθj, 0.0))
                    f_i = distribution_with_anisotropy(pi_sym, p_i, mi, μi, T, Φ, Φbar, ξ, cθi)
                    f_j = distribution_with_anisotropy(pj_sym, p_j, mj, μj, T, Φ, Φbar, ξ, cθj)
                    if f_i == 0.0 || f_j == 0.0
                        continue
                    end
                    for (φ, wφ) in zip(phi_grid, phi_w)
                        cosΘ = cθi * cθj + sθi * sθj * cos(φ)
                        s = mi^2 + mj^2 + 2.0 * (Ei * Ej - p_i * p_j * cosΘ)
                        if apply_s_cut
                            if (s <= s_bo) || (s >= s_up)
                                continue
                            end
                        else
                            if s <= (mi + mj)^2
                                continue
                            end
                        end

                        s_rt = sqrt(s)
                        Ei_cm = (s + mi^2 - mj^2) / (2.0 * s_rt)
                        Ej_cm = (s - mi^2 + mj^2) / (2.0 * s_rt)
                        pi_cm = sqrt(max(0.0, (s - (mi + mj)^2) * (s - (mi - mj)^2))) / (2.0 * s_rt)
                        pj_cm = pi_cm
                        v_rel_num = (Ei_cm * Ej_cm + pi_cm * pj_cm)^2 - (mi * mj)^2
                        v_rel = v_rel_num > 0.0 ? sqrt(v_rel_num) / (Ei_cm * Ej_cm) : 0.0
                        if v_rel == 0.0 || v_rel > 2.0
                            continue
                        end

                        σ = _get_sigma_core(cs_cache, s, quark_params, thermo_params, K_coeffs;
                            n_points=n_sigma_points, interpolation_mode=interpolation_mode)
                        base = w_pi * w_pj * w_cθi * w_cθj * wφ * (p_i^2) * (p_j^2) * f_i * f_j * v_rel
                        ω += base * σ
                        if band_omega !== nothing
                            ds = s - s_th
                            bi = _band_index(ds, band_edges)
                            if bi > 0
                                @inbounds band_omega[bi] += base
                                @inbounds band_omega_sigma[bi] += base * σ
                            end
                        end
                    end
                end
            end
        end
    end

    if band_omega !== nothing
        if band_omega_out !== nothing
            band_omega_out[] = band_omega
        end
        if band_omega_sigma_out !== nothing
            band_omega_sigma_out[] = band_omega_sigma
        end
    end

    return ω
end

end # module
