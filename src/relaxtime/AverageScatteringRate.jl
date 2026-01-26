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

include("../Constants_PNJL.jl")
include("../integration/GaussLegendre.jl")
include("TotalCrossSection.jl")
include("../QuarkDistribution.jl")
include("../QuarkDistribution_Aniso.jl")

using LinearAlgebra

using .Constants_PNJL: Λ_inv_fm
using .GaussLegendre: gauleg
using .TotalCrossSection: total_cross_section
using .TotalCrossSection: parse_particles_from_process
using .TotalCrossSection.ScatteringAmplitude.ParticleSymbols: is_antiquark
using .PNJLQuarkDistributions: quark_distribution, antiquark_distribution
using .PNJLQuarkDistributions_Aniso: quark_distribution_aniso, antiquark_distribution_aniso

# Import parameter types from Main
if !isdefined(Main, :ParameterTypes)
    Base.include(Main, joinpath(@__DIR__, "..", "ParameterTypes.jl"))
end
using Main.ParameterTypes: QuarkParams, ThermoParams, as_namedtuple

export average_scattering_rate, CrossSectionCache, precompute_cross_section!, build_w0cdf_pchip_cache

# Normalization helpers for dual interface support
@inline _nt_quark(q) = q isa QuarkParams ? as_namedtuple(q) : q
@inline _nt_thermo(t) = t isa ThermoParams ? as_namedtuple(t) : t

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



@inline function get_mass(flavor::Symbol, quark_params::Union{NamedTuple, QuarkParams})
    quark_params = _nt_quark(quark_params)
    if flavor === :u || flavor === :ubar; return quark_params.m.u
    elseif flavor === :d || flavor === :dbar; return quark_params.m.d
    elseif flavor === :s || flavor === :sbar; return quark_params.m.s
    else; error("Unknown flavor $flavor") end
end

@inline function get_mu(flavor::Symbol, quark_params::Union{NamedTuple, QuarkParams})
    quark_params = _nt_quark(quark_params)
    # Convention: always return the quark chemical potential μ_q (positive sign).
    # The particle/antiparticle distinction is handled by using
    # `quark_distribution*` vs `antiquark_distribution*`.
    if flavor === :u || flavor === :ubar; return quark_params.μ.u
    elseif flavor === :d || flavor === :dbar; return quark_params.μ.d
    elseif flavor === :s || flavor === :sbar; return quark_params.μ.s
    else; error("Unknown flavor $flavor") end
end

# -------------------- 截面缓存与插值 --------------------
mutable struct CrossSectionCache
    process::Symbol
    s_vals::Vector{Float64}
    sigma_vals::Vector{Float64}

    # Cached slopes for :pchip (same length as s_vals).
    pchip_slopes::Vector{Float64}
    pchip_dirty::Bool
end

CrossSectionCache(process::Symbol) = CrossSectionCache(process, Float64[], Float64[], Float64[], true)

function CrossSectionCache(process::Symbol, s_vals::Vector{Float64}, sigma_vals::Vector{Float64})
    length(s_vals) == length(sigma_vals) || error("CrossSectionCache: s_vals and sigma_vals length mismatch")
    cache = CrossSectionCache(process, s_vals, sigma_vals, Float64[], true)
    _ensure_pchip_slopes!(cache)
    return cache
end

function insert_sigma!(cache::CrossSectionCache, s::Float64, σ::Float64)
    idx = searchsortedfirst(cache.s_vals, s)
    insert!(cache.s_vals, idx, s)
    insert!(cache.sigma_vals, idx, σ)
    cache.pchip_dirty = true
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
    n_points::Int=TotalCrossSection.DEFAULT_T_INTEGRAL_POINTS)
    quark_params = _nt_quark(quark_params)
    thermo_params = _nt_thermo(thermo_params)
    for s in s_grid
        σ = total_cross_section(cache.process, s, quark_params, thermo_params, K_coeffs; n_points=n_points)
        insert_sigma!(cache, s, σ)
    end
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

function get_sigma(cache::CrossSectionCache, s::Float64,
    quark_params::Union{NamedTuple, QuarkParams}, thermo_params::Union{NamedTuple, ThermoParams}, K_coeffs::NamedTuple;
    n_points::Int=TotalCrossSection.DEFAULT_T_INTEGRAL_POINTS)
    quark_params = _nt_quark(quark_params)
    thermo_params = _nt_thermo(thermo_params)
    # Only cached PCHIP interpolation is supported.
    n = length(cache.s_vals)
    n == 0 && error("CrossSectionCache has no points; precompute σ(s) first")
    if n == 1
        return cache.sigma_vals[1]
    elseif s < cache.s_vals[1] || s > cache.s_vals[end]
        return 0.0
    elseif s == cache.s_vals[1]
        return cache.sigma_vals[1]
    elseif s == cache.s_vals[end]
        return cache.sigma_vals[end]
    end
    val = interpolate_sigma(cache, s)
    val === nothing && error("interpolation failed inside cache window")
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
    quark_params = _nt_quark(quark_params)
    thermo_params = _nt_thermo(thermo_params)
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
)
    quark_params = _nt_quark(quark_params)
    thermo_params = _nt_thermo(thermo_params)
    s_grid = design_w0cdf_s_grid(
        process,
        quark_params,
        thermo_params;
        N=N,
        p_nodes=design_p_nodes,
        angle_nodes=design_angle_nodes,
        phi_nodes=design_phi_nodes,
        p_cutoff=p_cutoff,
        scale=scale,
    )
    cache = CrossSectionCache(process)
    precompute_cross_section!(cache, s_grid, quark_params, thermo_params, K_coeffs; n_points=n_sigma_points)
    _ensure_pchip_slopes!(cache)
    return cache
end

# -------------------- ρ 计算（各向异性） --------------------
# 半无穷积分的默认参数
const DEFAULT_SEMI_INF_SCALE = 10.0  # 半无穷积分的尺度参数
const DEFAULT_SEMI_INF_NODES = 32    # 半无穷积分的节点数

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
    
    cos_grid === nothing && ((cos_grid, cos_w) = gauleg(0.0, 1.0, angle_nodes))

    # Momentum integration:
    # - Default: semi-infinite [0, ∞) via t∈[0,1) mapping.
    # - If p_grid/p_w are provided: integrate directly on that finite grid.
    p_vals = Float64[]
    p_wts = Float64[]
    dp_jac = Float64[]
    if p_grid !== nothing && p_w !== nothing
        length(p_grid) == length(p_w) || error("number_density: p_grid and p_w length mismatch")
        p_vals = p_grid
        p_wts = p_w
        dp_jac = ones(Float64, length(p_grid))
    else
        # 使用 [0, 1) 上的Gauss-Legendre节点，通过变换映射到 [0, ∞)
        # 变换: p = scale * t / (1-t), dp/dt = scale / (1-t)^2
        t_grid, t_w = gauleg(0.0, 1.0, p_nodes)
        for (t, wt) in zip(t_grid, t_w)
            # 避免 t=1 的奇点
            if t >= 0.9999
                continue
            end
            p = scale * t / (1.0 - t)
            dp_dt = scale / (1.0 - t)^2
            push!(p_vals, p)
            push!(p_wts, wt)
            push!(dp_jac, dp_dt)
        end
    end

    integral = 0.0
    for (p, wp, dp) in zip(p_vals, p_wts, dp_jac)
        for (cθ, wθ) in zip(cos_grid, cos_w)
            f = distribution_with_anisotropy(flavor, p, m, μ, T, Φ, Φbar, ξ, cθ)
            integral += wp * wθ * p^2 * f * dp
        end
    end
    # ρ = d_q / (2π^2) ∫ p^2 dp ∫_0^1 dcosθ f
    return DQ * integral / (2.0 * π^2)
end

# -------------------- 平均散射率主函数 --------------------
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
    sigma_cutoff::Union{Nothing,Float64}=nothing  # 新增：σ(s)有效范围的动量截断
)::Float64
    quark_params = _nt_quark(quark_params)
    thermo_params = _nt_thermo(thermo_params)
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
    cs_cache === nothing && (cs_cache = build_w0cdf_pchip_cache(
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
    ))

    # 使用半无穷积分 [0, ∞)
    return _average_scattering_rate_semi_infinite(
        process, pi_sym, pj_sym, mi, mj, μi, μj, T, Φ, Φbar, ξ,
        quark_params, thermo_params, K_coeffs,
        p_nodes, angle_nodes, phi_nodes, scale,
        p_grid, p_w, cos_grid, cos_w, phi_grid, phi_w,
        cs_cache, n_sigma_points,
        density_p_grid, density_p_w, density_p_nodes, density_scale,
        mc, md, apply_s_domain_cut, sigma_cutoff
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
    mc, md, apply_s_domain_cut, sigma_cutoff
)
    # Momentum integration:
    # - Default: semi-infinite [0, ∞) via t∈[0,1) mapping.
    # - If p_grid/p_w are provided: integrate directly on that finite grid.
    p_vals = Float64[]
    p_wts = Float64[]
    dp_jac = Float64[]
    if p_grid !== nothing && p_w !== nothing
        length(p_grid) == length(p_w) || error("average_scattering_rate: p_grid and p_w length mismatch")
        p_vals = p_grid
        p_wts = p_w
        dp_jac = ones(Float64, length(p_grid))
    else
        t_grid, t_w = gauleg(0.0, 1.0, p_nodes)
        for (t, wt) in zip(t_grid, t_w)
            if t >= 0.9999
                continue
            end
            p = scale * t / (1.0 - t)
            dp_dt = scale / (1.0 - t)^2
            push!(p_vals, p)
            push!(p_wts, wt)
            push!(dp_jac, dp_dt)
        end
    end

    # If numerator uses a finite cutoff grid, optionally apply the same s-domain cuts as Fortran:
    #   s_bo = max((mi+mj)^2, (mc+md)^2) * (1+1e-3)
    #   s_up = min((Ei_max+Ej_max)^2, (Ec_max+Ed_max)^2) with Ei_max = sqrt(mi^2+Λ^2)
    has_finite_p_cut = (p_grid !== nothing && p_w !== nothing)
    # 使用 sigma_cutoff 参数（如果提供）来确定 s 范围，否则从动量网格推断
    Λ = if sigma_cutoff !== nothing
        sigma_cutoff
    elseif has_finite_p_cut
        maximum(p_vals)
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
    # 使用完整积分区间 [-1,1] 和 [0,2π]
    cos_grid === nothing && ((cos_grid, cos_w) = gauleg(-1.0, 1.0, angle_nodes))
    phi_grid === nothing && ((phi_grid, phi_w) = gauleg(0.0, TWO_PI, phi_nodes))

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
    ω = _omega_integral_5d(
        process, pi_sym, pj_sym,
        mi, mj, μi, μj, T, Φ, Φbar, ξ,
        quark_params, thermo_params, K_coeffs,
        p_vals, p_wts, dp_jac,
        cos_grid, cos_w,
        phi_grid, phi_w,
        cs_cache, n_sigma_points,
        apply_s_cut, s_bo, s_up,
    )

    return prefactor * ω
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
    p_wts::Vector{Float64},
    dp_jac::Vector{Float64},
    cos_grid::Vector{Float64},
    cos_w::Vector{Float64},
    phi_grid::Vector{Float64},
    phi_w::Vector{Float64},
    cs_cache::CrossSectionCache,
    n_sigma_points::Int,
    apply_s_cut::Bool,
    s_bo::Float64,
    s_up::Float64,
)::Float64
    ω = 0.0

    for (p_i, w_pi, dp_i) in zip(p_vals, p_wts, dp_jac)
        Ei = energy_from_p(p_i, mi)

        for (p_j, w_pj, dp_j) in zip(p_vals, p_wts, dp_jac)
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

                        σ = get_sigma(cs_cache, s, quark_params, thermo_params, K_coeffs; n_points=n_sigma_points)
                        ω += w_pi * w_pj * w_cθi * w_cθj * wφ * (p_i^2) * (p_j^2) * f_i * f_j * v_rel * σ * dp_i * dp_j
                    end
                end
            end
        end
    end

    return ω
end

end # module
