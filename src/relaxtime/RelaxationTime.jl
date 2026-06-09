module RelaxationTime

"""
# RelaxationTime Module

基于平均散射率公式的弛豫时间计算。

## Features

- 仅需计算一次所需的平均散射率，即可复用至六种夸克种类（u、d、s、ubar、dbar、sbar）。
- 基于同位旋对称性，d和dbar共享u与ubar的计算结果，因此仅需显式评估四种夸克味。
- 粒子数密度 `densities` 由调用方提供（可预先计算或通过其他方式插值）。
- 若部分平均速率已存在，则通过 `existing_rates` 提供；缺失过程将自动补全。
- 同时返回 tau 和 tau_inv 以及可复用的平均速率。

## Dual Interface Pattern

This module supports **both struct and NamedTuple parameters** through a dual interface:

- **Structs (Recommended)**: Use `QuarkParams` and `ThermoParams` for type safety and better IDE support
- **NamedTuples (Backward Compatible)**: Existing NamedTuple-based code continues to work without modification

### Example with Structs

```julia
using Main.ParameterTypes: QuarkParams, ThermoParams

quark_params = QuarkParams(
    m=(u=1.52, d=1.52, s=3.04),
    μ=(u=0.3, d=0.3, s=0.3)
)
thermo_params = ThermoParams(0.15, 0.5, 0.5, 0.0)

result = relaxation_times(quark_params, thermo_params, K_coeffs; densities=densities)
```

### Example with NamedTuples

```julia
quark_params = (m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
thermo_params = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)

result = relaxation_times(quark_params, thermo_params, K_coeffs; densities=densities)
```

Both approaches produce identical results. See `docs/guides/PARAMETER_STRUCT_MIGRATION.md` for details.

## Internal Normalization

The module uses shared normalization helpers from ParameterAdapters to convert
struct inputs to NamedTuples at function boundaries. This ensures:
- Type stability in internal implementation
- Zero runtime overhead (helpers are inlined)
- Backward compatibility with existing code
"""

# Dependencies loaded by RelaxTime.jl entry point

const _VALIDATION_UTILS_PATH = normpath(joinpath(@__DIR__, "..", "utils", "ValidationUtils.jl"))
if !isdefined(Main, :ValidationUtils)
    Base.include(Main, _VALIDATION_UTILS_PATH)
end

using Main.ParameterTypes: QuarkParams, ThermoParams
using Main.ParameterAdapters: normalize_quark_input, normalize_thermo_input,
    normalize_symbol_mapping_input, lookup_symbol_value
using Main.ValidationUtils: validate_grid_weight_pair

using ..AverageScatteringRate: average_scattering_rate, CrossSectionCache,
    DEFAULT_P_NODES, DEFAULT_ANGLE_NODES, DEFAULT_PHI_NODES,
    DEFAULT_SIGMA_GRID_N, build_w0cdf_pchip_cache
using ..TotalCrossSection: DEFAULT_T_INTEGRAL_POINTS
using Main.Constants_PNJL: SCATTERING_PROCESS_KEYS, Λ_inv_fm
using ..AFieldBuilder: ensure_quark_params_has_A as _ensure_A

export relaxation_rates, relaxation_times, compute_average_rates, REQUIRED_PROCESSES

# Single source of truth for supported scattering processes.
# This list is derived from `Constants_PNJL.SCATTERING_MESON_MAP` keys.
const REQUIRED_PROCESSES = SCATTERING_PROCESS_KEYS

const RATE_ALIASES = (
    dubar_to_dubar=(:dubar_to_dubar, :udbar_to_udbar),
    subar_to_subar=(:subar_to_subar, :usbar_to_usbar),
    ubardbar_to_ubardbar=(:ubardbar_to_ubardbar, :ud_to_ud),
    ubarubar_to_ubarubar=(:ubarubar_to_ubarubar, :uu_to_uu),
    ubarsbar_to_ubarsbar=(:ubarsbar_to_ubarsbar, :us_to_us),
    sbarsbar_to_sbarsbar=(:sbarsbar_to_sbarsbar, :ss_to_ss),
)

@inline function ensure_quark_params_has_A(
    quark_params,
    thermo_params;
    p_nodes::Int=16,
    p_max::Float64=20.0,
    cos_nodes::Int=DEFAULT_ANGLE_NODES,
    use_aniso::Bool=true,
)::NamedTuple
    quark_nt = normalize_quark_input(quark_params)
    thermo_nt = normalize_thermo_input(thermo_params)
    return _ensure_A(
        quark_nt,
        thermo_nt;
        p_nodes=p_nodes,
        p_max=p_max,
        cos_nodes=cos_nodes,
        use_aniso=use_aniso,
        warn_on_auto=true,
    )
end

@inline density_value(densities::NamedTuple, key::Symbol) = lookup_symbol_value(densities, key, "densities")

@inline _has_rate(rates::NamedTuple, key::Symbol) = hasproperty(rates, key)
@inline _get_rate(rates::NamedTuple, key::Symbol) = getproperty(rates, key)

@inline function _resolve_rate_key(rates::NamedTuple, key::Symbol)::Symbol
    aliases = get(RATE_ALIASES, key, (key,))
    for candidate in aliases
        _has_rate(rates, candidate) && return candidate
    end
    return key
end

@inline function _rate_value_core(rates::NamedTuple, key::Symbol)
    resolved_key = _resolve_rate_key(rates, key)
    _has_rate(rates, resolved_key) || throw(ArgumentError("average rate for $(key) not found"))
    return _get_rate(rates, resolved_key)
end

@inline function _density_map(densities::NamedTuple)
    return (
        u=density_value(densities, :u),
        d=density_value(densities, :d),
        s=density_value(densities, :s),
        ubar=density_value(densities, :ubar),
        dbar=density_value(densities, :dbar),
        sbar=density_value(densities, :sbar),
    )
end

@inline function _rate_map(rates::NamedTuple)
    return (
        uu=rate_value(rates, :uu_to_uu),
        ud=rate_value(rates, :ud_to_ud),
        us=rate_value(rates, :us_to_us),
        udbar=rate_value(rates, :udbar_to_udbar),
        dubar=rate_value(rates, :dubar_to_dubar),
        uubar=rate_value(rates, :uubar_to_uubar),
        uubar_ddbar=rate_value(rates, :uubar_to_ddbar),
        usbar=rate_value(rates, :usbar_to_usbar),
        subar=rate_value(rates, :subar_to_subar),
        uubar_ssbar=rate_value(rates, :uubar_to_ssbar),
        ss=rate_value(rates, :ss_to_ss),
        ssbar_uubar=rate_value(rates, :ssbar_to_uubar),
        ssbar=rate_value(rates, :ssbar_to_ssbar),
        ubardbar=rate_value(rates, :ubardbar_to_ubardbar),
        ubarubar=rate_value(rates, :ubarubar_to_ubarubar),
        ubarsbar=rate_value(rates, :ubarsbar_to_ubarsbar),
        sbarsbar=rate_value(rates, :sbarsbar_to_sbarsbar),
    )
end

@inline function _species_omega_u(n, w)
    return n.u * (w.uu + w.ud) +
           n.ubar * (w.uubar + w.uubar_ddbar + w.uubar_ssbar + w.udbar) +
           n.s * w.us +
           n.sbar * w.usbar
end

@inline function _species_omega_s(n, w)
    return 2.0 * n.u * w.us +
           2.0 * n.ubar * w.usbar +
           n.s * w.ss +
           n.sbar * (w.ssbar + 2.0 * w.ssbar_uubar)
end

@inline function _species_omega_ubar(n, w)
    return n.u * (w.uubar + w.uubar_ddbar + w.uubar_ssbar + w.dubar) +
           n.ubar * (w.ubardbar + w.ubarubar) +
           n.s * w.subar +
           n.sbar * w.ubarsbar
end

@inline function _species_omega_sbar(n, w)
    return 2.0 * n.u * w.usbar +
           2.0 * n.ubar * w.ubarsbar +
           n.sbar * w.sbarsbar +
           n.s * (w.ssbar + 2.0 * w.ssbar_uubar)
end

@inline function _warn_and_clamp_nonnegative(omega_u::Real, omega_s::Real, omega_ubar::Real, omega_sbar::Real)
    if omega_u < -1e-12 || omega_s < -1e-12 || omega_ubar < -1e-12 || omega_sbar < -1e-12
        @warn "negative relaxation rate encountered; clamping to 0" omega_u=omega_u omega_s=omega_s omega_ubar=omega_ubar omega_sbar=omega_sbar
    end
    return max(omega_u, 0.0), max(omega_s, 0.0), max(omega_ubar, 0.0), max(omega_sbar, 0.0)
end

@inline function rate_value(rates, key::Symbol)
    return _rate_value_core(normalize_symbol_mapping_input(rates, "rates"), key)
end

function rate_lookup(args...)
    Base.depwarn("rate_lookup is deprecated; use rate_value instead.", :rate_lookup)
    return rate_value(args...)
end

"""
    compute_average_rates(quark_params, thermo_params, K_coeffs; kwargs...)

Compute missing averaged scattering rates while reusing any existing results or cross-section caches.

# Arguments
- `quark_params`: Quark parameters, either a `QuarkParams` struct or a NamedTuple with fields `m` and `μ`
- `thermo_params`: Thermodynamic parameters, either a `ThermoParams` struct or a NamedTuple with fields `T`, `Φ`, `Φbar`, `ξ`
- `K_coeffs`: Coupling coefficients as a NamedTuple
- `existing_rates`: Optional pre-computed rates to reuse
- `cs_caches`: Dictionary of cross-section caches for performance
- `p_nodes`, `angle_nodes`, `phi_nodes`: Integration node counts
- `p_grid`, `p_w`, `cos_grid`, `cos_w`, `phi_grid`, `phi_w`: Custom integration grids and weights
- `n_sigma_points`: Number of points for cross-section interpolation
- `sigma_grid_n`: Number of w0cdf sample points when building σ(s) caches automatically
- `sigma_cutoff`: Momentum cutoff for σ(s) effective range (defaults to Λ)
- `threshold_subtraction`, `asym_window`, `asym_fit_min_points`, `asym_extra_points`: Threshold-asymptotic cache controls forwarded to `average_scattering_rate`
- `interpolation_mode`: σ(s) evaluation mode forwarded to `average_scattering_rate`
- `require_cache_fingerprint`: Reject externally supplied σ(s) caches that do not carry fingerprint metadata
- `propagator_xi_policy`: `:match_thermo` keeps the current behavior; `:isotropic` evaluates σ(s)/propagators at `ξ=0` while retaining external distribution/density `ξ`
- `sigma_cache_policy`: σ(s) cache construction/evaluation policy; `:default` preserves the existing strategy, `:validated_anchored` enables a diagnostic anchored grid with local threshold addback

# Returns
A NamedTuple containing average scattering rates for all required processes.

# Examples
```julia
# Using structs (recommended)
q = QuarkParams(m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
t = ThermoParams(0.15, 0.5, 0.5, 0.0)
K = (K_pi=1.0, K_K=1.0, K_eta=1.0)
rates = compute_average_rates(q, t, K)

# Using NamedTuples (backward compatible)
q_nt = (m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
t_nt = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
rates = compute_average_rates(q_nt, t_nt, K)
```
"""
function compute_average_rates(
    quark_params::Union{NamedTuple, QuarkParams},
    thermo_params::Union{NamedTuple, ThermoParams},
    K_coeffs::NamedTuple;
    existing_rates::Union{Nothing,NamedTuple,AbstractDict}=nothing,
    cs_caches::Dict{Symbol,CrossSectionCache}=Dict{Symbol,CrossSectionCache}(),
    p_nodes::Int=DEFAULT_P_NODES,
    angle_nodes::Int=DEFAULT_ANGLE_NODES,
    phi_nodes::Int=DEFAULT_PHI_NODES,
    p_grid::Union{Nothing,Vector{Float64}}=nothing,
    p_w::Union{Nothing,Vector{Float64}}=nothing,
    cos_grid::Union{Nothing,Vector{Float64}}=nothing,
    cos_w::Union{Nothing,Vector{Float64}}=nothing,
    phi_grid::Union{Nothing,Vector{Float64}}=nothing,
    phi_w::Union{Nothing,Vector{Float64}}=nothing,
    n_sigma_points::Int=DEFAULT_T_INTEGRAL_POINTS,
    sigma_grid_n::Int=DEFAULT_SIGMA_GRID_N,
    sigma_cutoff::Union{Nothing,Float64}=nothing,  # σ(s)有效范围的动量截断，默认使用 Λ
    threshold_subtraction::Bool=false,
    asym_window::Float64=0.6,
    asym_fit_min_points::Int=8,
    asym_extra_points::Int=10,
    interpolation_mode::Symbol=:pchip,
    require_cache_fingerprint::Bool=false,
    propagator_xi_policy::Symbol=:match_thermo,
    propagator_quark_params::Union{Nothing,NamedTuple,QuarkParams}=nothing,
    sigma_cache_policy::Symbol=:default,
)::NamedTuple
    quark_nt = normalize_quark_input(quark_params)
    thermo_nt = normalize_thermo_input(thermo_params)
    existing_rates_nt = existing_rates === nothing ? nothing : normalize_symbol_mapping_input(existing_rates, "existing_rates")

    return _compute_average_rates_core(
        quark_nt,
        thermo_nt,
        K_coeffs;
        existing_rates=existing_rates_nt,
        cs_caches=cs_caches,
        p_nodes=p_nodes,
        angle_nodes=angle_nodes,
        phi_nodes=phi_nodes,
        p_grid=p_grid,
        p_w=p_w,
        cos_grid=cos_grid,
        cos_w=cos_w,
        phi_grid=phi_grid,
        phi_w=phi_w,
        n_sigma_points=n_sigma_points,
        sigma_grid_n=sigma_grid_n,
        sigma_cutoff=sigma_cutoff,
        threshold_subtraction=threshold_subtraction,
        asym_window=asym_window,
        asym_fit_min_points=asym_fit_min_points,
        asym_extra_points=asym_extra_points,
        interpolation_mode=interpolation_mode,
        require_cache_fingerprint=require_cache_fingerprint,
        propagator_xi_policy=propagator_xi_policy,
        propagator_quark_params=propagator_quark_params,
        sigma_cache_policy=sigma_cache_policy,
    )
end

function _compute_average_rates_core(
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    K_coeffs::NamedTuple;
    existing_rates::Union{Nothing,NamedTuple}=nothing,
    cs_caches::Dict{Symbol,CrossSectionCache}=Dict{Symbol,CrossSectionCache}(),
    p_nodes::Int=DEFAULT_P_NODES,
    angle_nodes::Int=DEFAULT_ANGLE_NODES,
    phi_nodes::Int=DEFAULT_PHI_NODES,
    p_grid::Union{Nothing,Vector{Float64}}=nothing,
    p_w::Union{Nothing,Vector{Float64}}=nothing,
    cos_grid::Union{Nothing,Vector{Float64}}=nothing,
    cos_w::Union{Nothing,Vector{Float64}}=nothing,
    phi_grid::Union{Nothing,Vector{Float64}}=nothing,
    phi_w::Union{Nothing,Vector{Float64}}=nothing,
    n_sigma_points::Int=DEFAULT_T_INTEGRAL_POINTS,
    sigma_grid_n::Int=DEFAULT_SIGMA_GRID_N,
    sigma_cutoff::Union{Nothing,Float64}=nothing,
    threshold_subtraction::Bool=false,
    asym_window::Float64=0.6,
    asym_fit_min_points::Int=8,
    asym_extra_points::Int=10,
    interpolation_mode::Symbol=:pchip,
    require_cache_fingerprint::Bool=false,
    propagator_xi_policy::Symbol=:match_thermo,
    propagator_quark_params::Union{Nothing,NamedTuple,QuarkParams}=nothing,
    sigma_cache_policy::Symbol=:default,
)::NamedTuple
    rates = Dict{Symbol,Float64}()
    if existing_rates !== nothing
        for (k, v) in pairs(existing_rates)
            rates[Symbol(k)] = v
        end
    end

    quark_params = ensure_quark_params_has_A(quark_params, thermo_params)

    # New strategy (2026-01-26):
    # - numerator momentum integrals p_i,p_j use semi-infinite [0,∞) integration
    # - σ(s) cache uses Λ cutoff: σ(s) = 0 when s exceeds Λ-based threshold
    # - number densities remain semi-infinite inside AverageScatteringRate
    validate_grid_weight_pair("compute_average_rates", "p_grid", p_grid, "p_w", p_w)
    validate_grid_weight_pair("compute_average_rates", "cos_grid", cos_grid, "cos_w", cos_w)
    validate_grid_weight_pair("compute_average_rates", "phi_grid", phi_grid, "phi_w", phi_w)
    
    # σ(s)缓存默认使用 Λ 截断
    effective_sigma_cutoff = sigma_cutoff === nothing ? Λ_inv_fm : sigma_cutoff
    
    # 动量积分不再使用截断，p_grid保持为nothing以触发半无穷积分
    # if p_grid === nothing: 将在 average_scattering_rate 中使用半无穷积分

    for process in REQUIRED_PROCESSES
        if haskey(rates, process)
            continue
        end
        # Reuse an existing populated cache if available; otherwise do NOT build it here.
        # Passing `nothing` as `cs_cache` allows `average_scattering_rate` to auto-build
        # the cache with the appropriate threshold-subtraction logic when needed.
        cache = haskey(cs_caches, process) ? cs_caches[process] : nothing
        cs_cache_arg = (cache === nothing || isempty(cache.s_vals)) ? nothing : cache

        rates[process] = average_scattering_rate(
            process,
            quark_params,
            thermo_params,
            K_coeffs;
            p_nodes=p_nodes,
            angle_nodes=angle_nodes,
            phi_nodes=phi_nodes,
            p_grid=p_grid,
            p_w=p_w,
            cos_grid=cos_grid,
            cos_w=cos_w,
            phi_grid=phi_grid,
            phi_w=phi_w,
            cs_cache=cs_cache_arg,
            n_sigma_points=n_sigma_points,
            sigma_grid_n=sigma_grid_n,
            sigma_cutoff=effective_sigma_cutoff,
            threshold_subtraction=threshold_subtraction,
            asym_window=asym_window,
            asym_fit_min_points=asym_fit_min_points,
            asym_extra_points=asym_extra_points,
            interpolation_mode=interpolation_mode,
            require_cache_fingerprint=require_cache_fingerprint,
            propagator_xi_policy=propagator_xi_policy,
            propagator_quark_params=propagator_quark_params,
            sigma_cache_policy=sigma_cache_policy,
        )
    end

    return NamedTuple(rates)
end

@inline _compute_average_rates_compat(args...; kwargs...) = _compute_average_rates_core(args...; kwargs...)
@inline _compute_average_rates_nt(args...; kwargs...) = _compute_average_rates_compat(args...; kwargs...)

# tau_i^-1 = sum_j rho_j * \bar{w}_{ij}
function relaxation_rates(
    densities::Union{NamedTuple,AbstractDict},
    rates::Union{NamedTuple,AbstractDict}
)::NamedTuple
    densities_nt = normalize_symbol_mapping_input(densities, "densities")
    rates_nt = normalize_symbol_mapping_input(rates, "rates")
    return _relaxation_rates_core(densities_nt, rates_nt)
end

function _relaxation_rates_core(
    densities::NamedTuple,
    rates::NamedTuple,
)::NamedTuple
    n = _density_map(densities)
    w = _rate_map(rates)

    omega_u = _species_omega_u(n, w)
    omega_s = _species_omega_s(n, w)
    # anti-u (shared with anti-d)
    # Matches Fortran: tau_lb = 1 / (
    #   n_u*(w6+w7+w9+wa5) + n_ub*(wa1+wa2) + n_s*wa6 + n_sb*wa3 )
    omega_ubar = _species_omega_ubar(n, w)
    # anti-s
    # Matches Fortran: tau_sb = 1 / (
    #   2*n_u*w8 + 2*n_ub*wa3 + n_sb*wa4 + n_s*(w11+2*w10) )
    omega_sbar = _species_omega_sbar(n, w)

    omega_u, omega_s, omega_ubar, omega_sbar = _warn_and_clamp_nonnegative(omega_u, omega_s, omega_ubar, omega_sbar)

    return (
        u = omega_u,
        d = omega_u,     # isospin symmetry
        s = omega_s,
        ubar = omega_ubar,
        dbar = omega_ubar, # isospin symmetry
        sbar = omega_sbar,
    )
end

@inline function safe_inv(x::Float64)
    return x == 0.0 ? Inf : 1.0 / x
end

const REQUIRED_RATE_KEYS_FOR_TAU = (
    :uu_to_uu,
    :ud_to_ud,
    :us_to_us,
    :udbar_to_udbar,
    :dubar_to_dubar,
    :uubar_to_uubar,
    :uubar_to_ddbar,
    :usbar_to_usbar,
    :subar_to_subar,
    :uubar_to_ssbar,
    :ss_to_ss,
    :ssbar_to_uubar,
    :ssbar_to_ssbar,
    :ubardbar_to_ubardbar,
    :ubarubar_to_ubarubar,
    :ubarsbar_to_ubarsbar,
    :sbarsbar_to_sbarsbar,
)

@inline function can_compute_tau_from_existing_rates(rates)::Bool
    try
        for k in REQUIRED_RATE_KEYS_FOR_TAU
            rate_value(rates, k)
        end
        return true
    catch err
        err isa InterruptException && rethrow()
        return false
    end
end

"""
    relaxation_times(quark_params, thermo_params, K_coeffs; kwargs...)

Calculate quark relaxation times based on average scattering rates.

This is the main entry point for computing relaxation times. It returns tau, tau_inv, 
and the averaged rates for reuse.

# Arguments
- `quark_params`: Quark parameters, either a `QuarkParams` struct or a NamedTuple with fields `m` and `μ`
- `thermo_params`: Thermodynamic parameters, either a `ThermoParams` struct or a NamedTuple with fields `T`, `Φ`, `Φbar`, `ξ`
- `K_coeffs`: Coupling coefficients as a NamedTuple
- `densities`: Particle number densities (required keyword argument)
- `existing_rates`: Optional pre-computed rates to reuse
- `cs_caches`: Dictionary of cross-section caches for performance
- `p_nodes`, `angle_nodes`, `phi_nodes`: Integration node counts
- `p_grid`, `p_w`, `cos_grid`, `cos_w`, `phi_grid`, `phi_w`: Custom integration grids and weights
- `n_sigma_points`: Number of points for cross-section interpolation
- `sigma_grid_n`: Number of w0cdf sample points when building σ(s) caches automatically
- `sigma_cutoff`: Momentum cutoff for σ(s) effective range
- `threshold_subtraction`, `asym_window`, `asym_fit_min_points`, `asym_extra_points`: Threshold-asymptotic cache controls forwarded to `average_scattering_rate`
- `interpolation_mode`: σ(s) evaluation mode forwarded to `average_scattering_rate`
- `require_cache_fingerprint`: Reject externally supplied σ(s) caches that do not carry fingerprint metadata
- `propagator_xi_policy`: `:match_thermo` keeps the current behavior; `:isotropic` evaluates σ(s)/propagators at `ξ=0` while retaining external distribution/density `ξ`
- `sigma_cache_policy`: σ(s) cache construction/evaluation policy; `:default` preserves the existing strategy, `:validated_anchored` enables a diagnostic anchored grid with local threshold addback

# Returns
A NamedTuple with fields:
- `tau`: Relaxation times for each quark flavor (u, d, s, ubar, dbar, sbar)
- `tau_inv`: Inverse relaxation times (scattering rates)
- `rates`: Average scattering rates for all processes

# Examples
```julia
# Using structs (recommended)
q = QuarkParams(m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
t = ThermoParams(0.15, 0.5, 0.5, 0.0)
K = (K_pi=1.0, K_K=1.0, K_eta=1.0)
densities = (u=0.1, d=0.1, s=0.05, ubar=0.1, dbar=0.1, sbar=0.05)
result = relaxation_times(q, t, K; densities=densities)

# Using NamedTuples (backward compatible)
q_nt = (m=(u=1.52, d=1.52, s=3.04), μ=(u=0.3, d=0.3, s=0.3))
t_nt = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
result = relaxation_times(q_nt, t_nt, K; densities=densities)
```
"""
function relaxation_times(
    quark_params::Union{NamedTuple, QuarkParams},
    thermo_params::Union{NamedTuple, ThermoParams},
    K_coeffs::NamedTuple;
    densities::Union{NamedTuple,AbstractDict},
    existing_rates::Union{Nothing,NamedTuple,AbstractDict}=nothing,
    cs_caches::Dict{Symbol,CrossSectionCache}=Dict{Symbol,CrossSectionCache}(),
    p_nodes::Int=DEFAULT_P_NODES,
    angle_nodes::Int=DEFAULT_ANGLE_NODES,
    phi_nodes::Int=DEFAULT_PHI_NODES,
    p_grid::Union{Nothing,Vector{Float64}}=nothing,
    p_w::Union{Nothing,Vector{Float64}}=nothing,
    cos_grid::Union{Nothing,Vector{Float64}}=nothing,
    cos_w::Union{Nothing,Vector{Float64}}=nothing,
    phi_grid::Union{Nothing,Vector{Float64}}=nothing,
    phi_w::Union{Nothing,Vector{Float64}}=nothing,
    n_sigma_points::Int=DEFAULT_T_INTEGRAL_POINTS,
    sigma_grid_n::Int=DEFAULT_SIGMA_GRID_N,
    sigma_cutoff::Union{Nothing,Float64}=nothing,
    threshold_subtraction::Bool=false,
    asym_window::Float64=0.6,
    asym_fit_min_points::Int=8,
    asym_extra_points::Int=10,
    interpolation_mode::Symbol=:pchip,
    require_cache_fingerprint::Bool=false,
    propagator_xi_policy::Symbol=:match_thermo,
    propagator_quark_params::Union{Nothing,NamedTuple,QuarkParams}=nothing,
    sigma_cache_policy::Symbol=:default,
)::NamedTuple
    quark_nt = normalize_quark_input(quark_params)
    thermo_nt = normalize_thermo_input(thermo_params)
    densities_nt = normalize_symbol_mapping_input(densities, "densities")
    existing_rates_nt = existing_rates === nothing ? nothing : normalize_symbol_mapping_input(existing_rates, "existing_rates")
    
    rates = if existing_rates_nt !== nothing && can_compute_tau_from_existing_rates(existing_rates_nt)
        existing_rates_nt
    else
        _compute_average_rates_core(
            quark_nt,
            thermo_nt,
            K_coeffs;
            existing_rates=existing_rates_nt,
            cs_caches=cs_caches,
            p_nodes=p_nodes,
            angle_nodes=angle_nodes,
            phi_nodes=phi_nodes,
            p_grid=p_grid,
            p_w=p_w,
            cos_grid=cos_grid,
            cos_w=cos_w,
            phi_grid=phi_grid,
            phi_w=phi_w,
            n_sigma_points=n_sigma_points,
            sigma_grid_n=sigma_grid_n,
            sigma_cutoff=sigma_cutoff,
            threshold_subtraction=threshold_subtraction,
            asym_window=asym_window,
            asym_fit_min_points=asym_fit_min_points,
            asym_extra_points=asym_extra_points,
            interpolation_mode=interpolation_mode,
            require_cache_fingerprint=require_cache_fingerprint,
            propagator_xi_policy=propagator_xi_policy,
            propagator_quark_params=propagator_quark_params,
            sigma_cache_policy=sigma_cache_policy,
        )
    end

    tau_inv = _relaxation_rates_core(densities_nt, rates)
    tau = (
        u = safe_inv(tau_inv.u),
        d = safe_inv(tau_inv.d),
        s = safe_inv(tau_inv.s),
        ubar = safe_inv(tau_inv.ubar),
        dbar = safe_inv(tau_inv.dbar),
        sbar = safe_inv(tau_inv.sbar),
    )

    return (tau=tau, tau_inv=tau_inv, rates=rates)
end

function _read_sigma_table_csv(path::AbstractString)
    s_vals = Float64[]
    sigma_vals = Float64[]
    for raw in eachline(path)
        line = strip(raw)
        isempty(line) && continue
        startswith(line, "#") && continue
        # Support either comma-separated or whitespace-separated formats.
        line = replace(line, ',' => ' ')
        parts = split(line)
        length(parts) < 2 && continue
        s_try = tryparse(Float64, parts[1])
        σ_try = tryparse(Float64, parts[2])
        (s_try === nothing || σ_try === nothing) && continue
        s = s_try
        σ = σ_try
        push!(s_vals, s)
        push!(sigma_vals, σ)
    end

    isempty(s_vals) && error("sigma table file has no data rows: $(path)")

    p = sortperm(s_vals)
    s_sorted = s_vals[p]
    σ_sorted = sigma_vals[p]

    # De-duplicate identical s values by keeping the last occurrence.
    s_out = Float64[]
    σ_out = Float64[]
    for (s, σ) in zip(s_sorted, σ_sorted)
        if !isempty(s_out) && s == s_out[end]
            σ_out[end] = σ
        else
            push!(s_out, s)
            push!(σ_out, σ)
        end
    end
    return (s_out, σ_out)
end

"""
    load_cross_section_caches_from_dir(dir) -> Dict{Symbol,CrossSectionCache}

从目录加载每个散射过程的 σ(s) 表（CSV），并构造 `cs_caches` 以注入到 `relaxation_times`。

本仓库的生产默认策略固定为 w0cdf+PCHIP，因此这里加载出的缓存会用 PCHIP 插值；当质心能量 s 超出缓存覆盖区间时，σ(s) 直接返回 0（而不是钳制到边界）。
运行时不会触发任何新的 σ(s) 计算。

目录内每个过程支持以下文件名之一：
- `sigma_<process>.csv`（推荐）
- `<process>.csv`

每个 CSV 的数据行格式为：
- `s,sigma` 或 `s sigma`（允许 # 开头注释行）
"""
function load_cross_section_caches_from_dir(dir::AbstractString)::Dict{Symbol,CrossSectionCache}
    isdir(dir) || error("sigma cache directory not found: $(dir)")

    cs_caches = Dict{Symbol,CrossSectionCache}()
    for process in REQUIRED_PROCESSES
        path1 = joinpath(dir, "sigma_$(process).csv")
        path2 = joinpath(dir, "$(process).csv")
        path = isfile(path1) ? path1 : (isfile(path2) ? path2 : "")
        isempty(path) && error("missing sigma table for $(process) under $(dir) (expected $(path1) or $(path2))")

        s_vals, σ_vals = _read_sigma_table_csv(path)
        cache = CrossSectionCache(process)
        cache.s_vals = s_vals
        cache.sigma_vals = σ_vals
        cache.pchip_dirty = true
        cs_caches[process] = cache
    end
    return cs_caches
end

end # module RelaxationTime
