module RelaxationTime

"""
基于平均散射率公式的弛豫时间计算。
- 仅需计算一次所需的平均散射率，即可复用至六种夸克种类（u、d、s、ubar、dbar、sbar）。
- 基于同位旋对称性，d和dbar共享u与ubar的计算结果，因此仅需显式评估四种夸克味。
- 粒子数密度 `densities` 由调用方提供（可预先计算或通过其他方式插值）。
- 若部分平均速率已存在，则通过 `existing_rates` 提供；缺失过程将自动补全。

同时返回 tau 和 tau_inv 以及可复用的平均速率。
"""

include("AverageScatteringRate.jl")
include("TotalCrossSection.jl")

using .AverageScatteringRate: average_scattering_rate, CrossSectionCache,
    DEFAULT_P_NODES, DEFAULT_ANGLE_NODES, DEFAULT_PHI_NODES
using .TotalCrossSection: DEFAULT_T_INTEGRAL_POINTS

export relaxation_rates, relaxation_times, compute_average_rates, REQUIRED_PROCESSES

# Processes needed by the documented tau expressions
const REQUIRED_PROCESSES = (
    :uu_to_uu,
    :ud_to_ud,
    :us_to_us,
    :usbar_to_usbar,
    :uubar_to_uubar,
    :uubar_to_ddbar,
    :uubar_to_ssbar,
    :udbar_to_udbar,
    :ss_to_ss,
    :ssbar_to_ssbar,
    :ssbar_to_uubar,
)

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
        hasproperty(rates, key) || error("average rate for $(key) not found")
        return getproperty(rates, key)
    elseif rates isa AbstractDict
        haskey(rates, key) || error("average rate for $(key) not found")
        return rates[key]
    else
        error("rates must be a NamedTuple or Dict")
    end
end

# Compute missing averaged scattering rates while reusing any existing results or cross-section caches.
function compute_average_rates(
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
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
    n_sigma_points::Int=DEFAULT_T_INTEGRAL_POINTS
)::NamedTuple
    rates = Dict{Symbol,Float64}()
    if existing_rates !== nothing
        for (k, v) in pairs(existing_rates)
            rates[Symbol(k)] = v
        end
    end

    for process in REQUIRED_PROCESSES
        if haskey(rates, process)
            continue
        end
        cache = get!(cs_caches, process) do
            CrossSectionCache(process)
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
    w_usbar = rate_lookup(rates, :usbar_to_usbar)
    w_uubar = rate_lookup(rates, :uubar_to_uubar)
    w_uubar_ddbar = rate_lookup(rates, :uubar_to_ddbar)
    w_uubar_ssbar = rate_lookup(rates, :uubar_to_ssbar)
    w_udbar = rate_lookup(rates, :udbar_to_udbar)
    w_ss = rate_lookup(rates, :ss_to_ss)
    w_ssbar = rate_lookup(rates, :ssbar_to_ssbar)
    w_ssbar_uubar = rate_lookup(rates, :ssbar_to_uubar)

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
    omega_ubar = n_u * (w_uubar + w_udbar + w_usbar + w_uubar_ddbar) +
                 n_ubar * (w_uubar + w_udbar) +
                 n_s * w_usbar +
                 n_sbar * w_usbar

    # anti-s
    omega_sbar = 2.0 * n_u * w_usbar +
                 2.0 * n_ubar * w_usbar +
                 n_s * w_ssbar +
                 n_sbar * (w_ssbar + 2.0 * w_ssbar_uubar)

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

# Main entry: returns tau, tau_inv, and the averaged rates for reuse.
function relaxation_times(
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
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
    n_sigma_points::Int=DEFAULT_T_INTEGRAL_POINTS
)::NamedTuple
    rates = compute_average_rates(
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
    )

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

end # module
