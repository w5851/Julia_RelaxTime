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

The module uses internal normalization helpers `_nt_quark` and `_nt_thermo` to convert
struct inputs to NamedTuples at function boundaries. This ensures:
- Type stability in internal implementation
- Zero runtime overhead (helpers are inlined)
- Backward compatibility with existing code
"""

include("AverageScatteringRate.jl")
include("TotalCrossSection.jl")
include("OneLoopIntegrals.jl")
include("../Constants_PNJL.jl")

# Ensure shared parameter types are loaded for cross-module reuse
if !isdefined(Main, :ParameterTypes)
    Base.include(Main, joinpath(@__DIR__, "..", "ParameterTypes.jl"))
end

using Main.ParameterTypes: QuarkParams, ThermoParams, as_namedtuple

using .AverageScatteringRate: average_scattering_rate, CrossSectionCache,
    DEFAULT_P_NODES, DEFAULT_ANGLE_NODES, DEFAULT_PHI_NODES,
    build_w0cdf_pchip_cache
using .OneLoopIntegrals: A
using .TotalCrossSection: DEFAULT_T_INTEGRAL_POINTS
using .Constants_PNJL: SCATTERING_PROCESS_KEYS, Λ_inv_fm

export relaxation_rates, relaxation_times, compute_average_rates, REQUIRED_PROCESSES

@inline _nt_quark(q) = q isa QuarkParams ? as_namedtuple(q) : q
@inline _nt_thermo(t) = t isa ThermoParams ? as_namedtuple(t) : t

# Single source of truth for supported scattering processes.
# This list is derived from `Constants_PNJL.SCATTERING_MESON_MAP` keys.
const REQUIRED_PROCESSES = SCATTERING_PROCESS_KEYS

@inline function ensure_quark_params_has_A(quark_params, thermo_params)::NamedTuple
    quark_params = _nt_quark(quark_params)
    thermo_params = _nt_thermo(thermo_params)
    # Many low-level scattering routines require `quark_params.A` for polarization functions.
    # Older callers/tests may only provide (m, μ). In that case we compute A on-demand here.
    if hasproperty(quark_params, :A)
        return quark_params
    end
    hasproperty(quark_params, :m) || error("quark_params is missing :m")
    hasproperty(quark_params, :μ) || error("quark_params is missing :μ")
    hasproperty(thermo_params, :T) || error("thermo_params is missing :T")
    hasproperty(thermo_params, :Φ) || error("thermo_params is missing :Φ")
    hasproperty(thermo_params, :Φbar) || error("thermo_params is missing :Φbar")

    # A 的热部分对高温/轻质量参数更敏感：这里使用更稳健的上限/节点配置。
    nodes_p, weights_p = AverageScatteringRate.gauleg(0.0, 20.0, 16)
    A_u = A(quark_params.m.u, quark_params.μ.u, thermo_params.T, thermo_params.Φ, thermo_params.Φbar, nodes_p, weights_p)
    A_d = A(quark_params.m.d, quark_params.μ.d, thermo_params.T, thermo_params.Φ, thermo_params.Φbar, nodes_p, weights_p)
    A_s = A(quark_params.m.s, quark_params.μ.s, thermo_params.T, thermo_params.Φ, thermo_params.Φbar, nodes_p, weights_p)

    return merge(quark_params, (A=(u=A_u, d=A_d, s=A_s),))
end

@inline function density_lookup(densities, key::Symbol)
    if densities isa NamedTuple
        hasproperty(densities, key) || error("densities is missing $(key)")
        return getproperty(densities, key)
    elseif densities isa AbstractDict
        haskey(densities, key) || error("densities is missing $(key)")
        return densities[key]
    else
        error("densities must be a NamedTuple or Dict with keys :u, :d, :s, :ubar, :dbar, :sbar")
    end
end

@inline function rate_lookup(rates, key::Symbol)
    if rates isa NamedTuple
        # Backward/ergonomic aliases (isospin / charge conjugation).
        # These help older callers/tests that only provide a reduced rate set.
        if key === :dubar_to_dubar && !hasproperty(rates, :dubar_to_dubar) && hasproperty(rates, :udbar_to_udbar)
            key = :udbar_to_udbar
        elseif key === :subar_to_subar && !hasproperty(rates, :subar_to_subar) && hasproperty(rates, :usbar_to_usbar)
            key = :usbar_to_usbar
        elseif key === :ubardbar_to_ubardbar && !hasproperty(rates, :ubardbar_to_ubardbar) && hasproperty(rates, :ud_to_ud)
            key = :ud_to_ud
        elseif key === :ubarubar_to_ubarubar && !hasproperty(rates, :ubarubar_to_ubarubar) && hasproperty(rates, :uu_to_uu)
            key = :uu_to_uu
        elseif key === :ubarsbar_to_ubarsbar && !hasproperty(rates, :ubarsbar_to_ubarsbar) && hasproperty(rates, :us_to_us)
            key = :us_to_us
        elseif key === :sbarsbar_to_sbarsbar && !hasproperty(rates, :sbarsbar_to_sbarsbar) && hasproperty(rates, :ss_to_ss)
            key = :ss_to_ss
        end
        hasproperty(rates, key) || error("average rate for $(key) not found")
        return getproperty(rates, key)
    elseif rates isa AbstractDict
        if key === :dubar_to_dubar && !haskey(rates, :dubar_to_dubar) && haskey(rates, :udbar_to_udbar)
            key = :udbar_to_udbar
        elseif key === :subar_to_subar && !haskey(rates, :subar_to_subar) && haskey(rates, :usbar_to_usbar)
            key = :usbar_to_usbar
        elseif key === :ubardbar_to_ubardbar && !haskey(rates, :ubardbar_to_ubardbar) && haskey(rates, :ud_to_ud)
            key = :ud_to_ud
        elseif key === :ubarubar_to_ubarubar && !haskey(rates, :ubarubar_to_ubarubar) && haskey(rates, :uu_to_uu)
            key = :uu_to_uu
        elseif key === :ubarsbar_to_ubarsbar && !haskey(rates, :ubarsbar_to_ubarsbar) && haskey(rates, :us_to_us)
            key = :us_to_us
        elseif key === :sbarsbar_to_sbarsbar && !haskey(rates, :sbarsbar_to_sbarsbar) && haskey(rates, :ss_to_ss)
            key = :ss_to_ss
        end
        haskey(rates, key) || error("average rate for $(key) not found")
        return rates[key]
    else
        error("rates must be a NamedTuple or Dict")
    end
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
- `sigma_cutoff`: Momentum cutoff for σ(s) effective range (defaults to Λ)

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
    sigma_cutoff::Union{Nothing,Float64}=nothing  # σ(s)有效范围的动量截断，默认使用 Λ
)::NamedTuple
    # Normalize inputs at function entry
    quark_params = _nt_quark(quark_params)
    thermo_params = _nt_thermo(thermo_params)
    
    rates = Dict{Symbol,Float64}()
    if existing_rates !== nothing
        for (k, v) in pairs(existing_rates)
            rates[Symbol(k)] = v
        end
    end

    quark_params = ensure_quark_params_has_A(quark_params, thermo_params)

    # Unified standard (default):
    # - numerator momentum integrals p_i,p_j use the PNJL cutoff p∈[0,Λ]
    # - σ(s) cache also uses Λ cutoff for consistency
    # - number densities remain semi-infinite inside AverageScatteringRate
    if (p_grid === nothing) != (p_w === nothing)
        error("compute_average_rates: p_grid and p_w must be provided together")
    end
    
    # 默认使用 Λ 截断
    effective_sigma_cutoff = sigma_cutoff === nothing ? Λ_inv_fm : sigma_cutoff
    
    if p_grid === nothing
        # NOTE: `Λ_inv_fm` is the PNJL momentum cutoff Λ in units of fm⁻¹.
        # Do NOT invert it; otherwise the integration upper bound becomes ~0.3 fm⁻¹
        # instead of ~3 fm⁻¹, suppressing phase space and inflating τ.
        Λ = Λ_inv_fm
        # Use the requested p_nodes as the cutoff-grid resolution.
        p_grid, p_w = AverageScatteringRate.gauleg(0.0, Λ, p_nodes)
    end

    for process in REQUIRED_PROCESSES
        if haskey(rates, process)
            continue
        end
        cache = get!(cs_caches, process) do
            CrossSectionCache(process)
        end

        # If the cache is still empty, build the default designed σ-grid + PCHIP cache.
        # `average_scattering_rate` assumes any provided cache is already populated.
        # 使用 effective_sigma_cutoff 确保 σ(s) 缓存范围与动量积分范围一致
        if isempty(cache.s_vals)
            cache = build_w0cdf_pchip_cache(
                process,
                quark_params,
                thermo_params,
                K_coeffs;
                p_cutoff=effective_sigma_cutoff,
                n_sigma_points=n_sigma_points,
            )
            cs_caches[process] = cache
        end

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
            cs_cache=cache,
            n_sigma_points=n_sigma_points,
            sigma_cutoff=sigma_cutoff,
        )
    end

    return NamedTuple(rates)
end

# tau_i^-1 = sum_j rho_j * \bar{w}_{ij}
function relaxation_rates(
    densities::Union{NamedTuple,AbstractDict},
    rates::Union{NamedTuple,AbstractDict}
)::NamedTuple
    n_u = density_lookup(densities, :u)
    n_d = density_lookup(densities, :d)
    n_s = density_lookup(densities, :s)
    n_ubar = density_lookup(densities, :ubar)
    n_dbar = density_lookup(densities, :dbar)
    n_sbar = density_lookup(densities, :sbar)

    w_uu = rate_lookup(rates, :uu_to_uu)
    w_ud = rate_lookup(rates, :ud_to_ud)
    w_us = rate_lookup(rates, :us_to_us)
    w_udbar = rate_lookup(rates, :udbar_to_udbar)
    w_dubar = rate_lookup(rates, :dubar_to_dubar)
    w_uubar = rate_lookup(rates, :uubar_to_uubar)
    w_uubar_ddbar = rate_lookup(rates, :uubar_to_ddbar)
    w_usbar = rate_lookup(rates, :usbar_to_usbar)
    w_subar = rate_lookup(rates, :subar_to_subar)
    w_uubar_ssbar = rate_lookup(rates, :uubar_to_ssbar)
    w_ss = rate_lookup(rates, :ss_to_ss)
    w_ssbar_uubar = rate_lookup(rates, :ssbar_to_uubar)
    w_ssbar = rate_lookup(rates, :ssbar_to_ssbar)

    # Additional rates for antiquark relaxation times
    w_ubardbar = rate_lookup(rates, :ubardbar_to_ubardbar)
    w_ubarubar = rate_lookup(rates, :ubarubar_to_ubarubar)
    w_ubarsbar = rate_lookup(rates, :ubarsbar_to_ubarsbar)
    w_sbarsbar = rate_lookup(rates, :sbarsbar_to_sbarsbar)

    # u quark (shared with d by isospin symmetry)
    omega_u = n_u * (w_uu + w_ud) +
              n_ubar * (w_uubar + w_uubar_ddbar + w_uubar_ssbar + w_udbar) +
              n_s * w_us +
              n_sbar * w_usbar

    # s quark
    omega_s = 2.0 * n_u * w_us +
              2.0 * n_ubar * w_usbar +
              n_s * w_ss +
              n_sbar * (w_ssbar + 2.0 * w_ssbar_uubar)

    # anti-u (shared with anti-d)
    # Matches Fortran: tau_lb = 1 / (
    #   n_u*(w6+w7+w9+wa5) + n_ub*(wa1+wa2) + n_s*wa6 + n_sb*wa3 )
    omega_ubar = n_u * (w_uubar + w_uubar_ddbar + w_uubar_ssbar + w_dubar) +
                 n_ubar * (w_ubardbar + w_ubarubar) +
                 n_s * w_subar +
                 n_sbar * w_ubarsbar

    # anti-s
    # Matches Fortran: tau_sb = 1 / (
    #   2*n_u*w8 + 2*n_ub*wa3 + n_sb*wa4 + n_s*(w11+2*w10) )
    omega_sbar = 2.0 * n_u * w_usbar +
                 2.0 * n_ubar * w_ubarsbar +
                 n_sbar * w_sbarsbar +
                 n_s * (w_ssbar + 2.0 * w_ssbar_uubar)

    # 数值保护：平均速率应非负；若出现明显负值，提示但仍钳制到 0。
    if omega_u < -1e-12 || omega_s < -1e-12 || omega_ubar < -1e-12 || omega_sbar < -1e-12
        @warn "negative relaxation rate encountered; clamping to 0" omega_u=omega_u omega_s=omega_s omega_ubar=omega_ubar omega_sbar=omega_sbar
    end
    omega_u = max(omega_u, 0.0)
    omega_s = max(omega_s, 0.0)
    omega_ubar = max(omega_ubar, 0.0)
    omega_sbar = max(omega_sbar, 0.0)

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
            rate_lookup(rates, k)
        end
        return true
    catch
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
- `sigma_cutoff`: Momentum cutoff for σ(s) effective range

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
    sigma_cutoff::Union{Nothing,Float64}=nothing  # 新增：σ(s)有效范围的动量截断
)::NamedTuple
    # Normalize inputs at function entry
    quark_params = _nt_quark(quark_params)
    thermo_params = _nt_thermo(thermo_params)
    # Normalize inputs at function entry
    quark_params = _nt_quark(quark_params)
    thermo_params = _nt_thermo(thermo_params)
    
    rates = if existing_rates !== nothing && can_compute_tau_from_existing_rates(existing_rates)
        existing_rates isa NamedTuple ? existing_rates : (; (Symbol(k) => v for (k, v) in pairs(existing_rates))...)
    else
        compute_average_rates(
            quark_params,
            thermo_params,
            K_coeffs;
            existing_rates=existing_rates,
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
            sigma_cutoff=sigma_cutoff,
        )
    end

    tau_inv = relaxation_rates(densities, rates)
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

end # module
