module AverageScatteringRate

"""
平均散射率计算模块（各向异性可选）。
- 积分采用 Gauss-Legendre：动量 32 节点，角度 4 节点（可覆盖 7 阶多项式）。
- 支持各向异性分布 `quark_distribution_aniso(..., ξ, cosθ)`；`ξ=0` 时退化为各向同性。
- 散射截面使用 `TotalCrossSection.total_cross_section`，可按需预计算并插值。
"""

include("../Constants_PNJL.jl")
include("../integration/GaussLegendre.jl")
include("TotalCrossSection.jl")
include("../QuarkDistribution.jl")
include("../QuarkDistribution_Aniso.jl")

using .Constants_PNJL: Λ_inv_fm
using .GaussLegendre: gauleg
using .TotalCrossSection: total_cross_section
using .PNJLQuarkDistributions: quark_distribution, antiquark_distribution
using .PNJLQuarkDistributions_Aniso: quark_distribution_aniso, antiquark_distribution_aniso

export average_scattering_rate, CrossSectionCache, precompute_cross_section!

const DEFAULT_P_NODES = 32
const DEFAULT_ANGLE_NODES = 4
const DQ = 6.0 # 简并度d_q=2*N_c=6
const TWO_PI = 2.0 * π

# --------------------------- 工具函数 ---------------------------

@inline function is_antiquark(flavor::Symbol)::Bool
    return flavor === :ubar || flavor === :dbar || flavor === :sbar
end

@inline function distribution_with_anisotropy(flavor::Symbol, p::Float64, m::Float64, μ::Float64,
    T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, cosθ::Float64)
    if ξ == 0.0
        return is_antiquark(flavor) ? antiquark_distribution(p, μ, T, Φ, Φbar) : quark_distribution(p, μ, T, Φ, Φbar)
    else
        return is_antiquark(flavor) ? antiquark_distribution_aniso(p, m, μ, T, Φ, Φbar, ξ, cosθ) : quark_distribution_aniso(p, m, μ, T, Φ, Φbar, ξ, cosθ)
    end
end

@inline function energy_from_p(p::Float64, m::Float64)
    return sqrt(p * p + m * m)
end

# 解析 process: 复用 TotalCrossSection 里约定的字符串方案
@inline function parse_particle_pair(pair_str::String)::Tuple{Symbol, Symbol}
    pair_str = replace(pair_str, "ū" => "ubar")
    pair_str = replace(pair_str, "đ" => "dbar")
    pair_str = replace(pair_str, "s̄" => "sbar")
    particles = Symbol[]
    i = 1
    while i <= length(pair_str)
        if i + 3 <= length(pair_str) && pair_str[i:i+3] == "ubar"
            push!(particles, :ubar); i += 4
        elseif i + 3 <= length(pair_str) && pair_str[i:i+3] == "dbar"
            push!(particles, :dbar); i += 4
        elseif i + 3 <= length(pair_str) && pair_str[i:i+3] == "sbar"
            push!(particles, :sbar); i += 4
        elseif pair_str[i] in ['u','d','s']
            push!(particles, Symbol(pair_str[i])); i += 1
        else
            error("Unknown particle symbol at position $i in '$pair_str'")
        end
    end
    length(particles) == 2 || error("Expected 2 particles in '$pair_str'")
    return (particles[1], particles[2])
end

@inline function parse_particles_from_process(process::Symbol)::Tuple{Symbol,Symbol,Symbol,Symbol}
    parts = split(string(process), "_to_")
    length(parts) == 2 || error("Invalid process format: $process")
    i,j = parse_particle_pair(String(parts[1]))
    c,d = parse_particle_pair(String(parts[2]))
    return (i,j,c,d)
end

@inline function get_mass(flavor::Symbol, quark_params::NamedTuple)
    if flavor in [:u,:ubar]; return quark_params.m.u
    elseif flavor in [:d,:dbar]; return quark_params.m.d
    elseif flavor in [:s,:sbar]; return quark_params.m.s
    else; error("Unknown flavor $flavor") end
end

@inline function get_mu(flavor::Symbol, quark_params::NamedTuple)
    if flavor == :u; return quark_params.μ.u
    elseif flavor == :ubar; return -quark_params.μ.u
    elseif flavor == :d; return quark_params.μ.d
    elseif flavor == :dbar; return -quark_params.μ.d
    elseif flavor == :s; return quark_params.μ.s
    elseif flavor == :sbar; return -quark_params.μ.s
    else; error("Unknown flavor $flavor") end
end

# -------------------- 截面缓存与插值 --------------------
mutable struct CrossSectionCache
    process::Symbol
    s_vals::Vector{Float64}
    sigma_vals::Vector{Float64}
end

function CrossSectionCache(process::Symbol)
    return CrossSectionCache(process, Float64[], Float64[])
end

function insert_sigma!(cache::CrossSectionCache, s::Float64, σ::Float64)
    idx = searchsortedfirst(cache.s_vals, s)
    insert!(cache.s_vals, idx, s)
    insert!(cache.sigma_vals, idx, σ)
end

function precompute_cross_section!(cache::CrossSectionCache, s_grid::Vector{Float64},
    quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple;
    n_points::Int=TotalCrossSection.DEFAULT_T_INTEGRAL_POINTS)
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
    elseif s <= cache.s_vals[1]
        return cache.sigma_vals[1]
    elseif s >= cache.s_vals[end]
        return cache.sigma_vals[end]
    else
        idx = searchsortedfirst(cache.s_vals, s)
        s1, s2 = cache.s_vals[idx-1], cache.s_vals[idx]
        σ1, σ2 = cache.sigma_vals[idx-1], cache.sigma_vals[idx]
        w = (s - s1)/(s2 - s1)
        return σ1 + w*(σ2 - σ1)
    end
end

function get_sigma(cache::CrossSectionCache, s::Float64,
    quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple;
    n_points::Int=TotalCrossSection.DEFAULT_T_INTEGRAL_POINTS)
    σ_interp = interpolate_sigma(cache, s)
    if σ_interp !== nothing
        return σ_interp
    else
        σ = total_cross_section(cache.process, s, quark_params, thermo_params, K_coeffs; n_points=n_points)
        insert_sigma!(cache, s, σ)
        return σ
    end
end

# -------------------- ρ 计算（各向异性） --------------------
function number_density(flavor::Symbol, m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64;
    p_nodes::Int=DEFAULT_P_NODES, angle_nodes::Int=DEFAULT_ANGLE_NODES)
    p_grid, p_w = gauleg(0.0, Λ_inv_fm, p_nodes)
    cos_grid, cos_w = gauleg(0.0, 1.0, angle_nodes)
    integral = 0.0
    for (p, wp) in zip(p_grid, p_w)
        for (cθ, wθ) in zip(cos_grid, cos_w)
            f = distribution_with_anisotropy(flavor, p, m, μ, T, Φ, Φbar, ξ, cθ)
            integral += wp * wθ * p^2 * f
        end
    end
    # ρ = d_q / (2π^2) ∫ p^2 dp ∫_0^1 dcosθ f
    return DQ * integral / (2.0 * π^2)
end

# -------------------- 平均散射率主函数 --------------------
function average_scattering_rate(
    process::Symbol,
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    K_coeffs::NamedTuple;
    p_nodes::Int=DEFAULT_P_NODES,
    angle_nodes::Int=DEFAULT_ANGLE_NODES,
    phi_nodes::Int=DEFAULT_ANGLE_NODES,
    cs_cache::CrossSectionCache=CrossSectionCache(process),
    n_sigma_points::Int=TotalCrossSection.DEFAULT_T_INTEGRAL_POINTS
)::Float64
    # 解析粒子、质量、化学势
    pi_sym, pj_sym, pc_sym, pd_sym = parse_particles_from_process(process)
    mi = get_mass(pi_sym, quark_params); mj = get_mass(pj_sym, quark_params)
    mc = get_mass(pc_sym, quark_params); md = get_mass(pd_sym, quark_params)
    μi = get_mu(pi_sym, quark_params); μj = get_mu(pj_sym, quark_params)
    μc = get_mu(pc_sym, quark_params); μd = get_mu(pd_sym, quark_params)

    T = thermo_params.T; Φ = thermo_params.Φ; Φbar = thermo_params.Φbar
    ξ = hasproperty(thermo_params, :ξ) ? thermo_params.ξ : 0.0

    # 预备节点
    p_grid, p_w = gauleg(0.0, Λ_inv_fm, p_nodes)
    cos_grid, cos_w = gauleg(0.0, 1.0, angle_nodes)
    phi_grid, phi_w = gauleg(0.0, Float64(π), phi_nodes)

    # 数密度（用于归一化）
    ρ_i = number_density(pi_sym, mi, μi, T, Φ, Φbar, ξ; p_nodes=p_nodes, angle_nodes=angle_nodes)
    ρ_j = number_density(pj_sym, mj, μj, T, Φ, Φbar, ξ; p_nodes=p_nodes, angle_nodes=angle_nodes)
    if ρ_i == 0.0 || ρ_j == 0.0
        return 0.0
    end

    prefactor = (DQ^2) / (4.0 * π^5 * ρ_i * ρ_j)
    ω = 0.0

    for (p_i, w_pi) in zip(p_grid, p_w)
        Ei = energy_from_p(p_i, mi)
        for (p_j, w_pj) in zip(p_grid, p_w)
            Ej = energy_from_p(p_j, mj)
            for (cθi, w_cθi) in zip(cos_grid, cos_w)
                sθi = sqrt(max(1.0 - cθi*cθi, 0.0))
                for (cθj, w_cθj) in zip(cos_grid, cos_w)
                    sθj = sqrt(max(1.0 - cθj*cθj, 0.0))
                    f_i = distribution_with_anisotropy(pi_sym, p_i, mi, μi, T, Φ, Φbar, ξ, cθi)
                    f_j = distribution_with_anisotropy(pj_sym, p_j, mj, μj, T, Φ, Φbar, ξ, cθj)
                    if f_i == 0.0 || f_j == 0.0
                        continue
                    end
                    for (φ, wφ) in zip(phi_grid, phi_w)
                        cosΘ = cθi * cθj + sθi * sθj * cos(φ)
                        p_dot = Ei * Ej - p_i * p_j * cosΘ
                        s = mi^2 + mj^2 + 2.0 * p_dot
                        if s <= (mi + mj)^2
                            continue
                        end
                        v_rel_num = p_dot^2 - mi^2 * mj^2
                        v_rel = v_rel_num > 0 ? sqrt(v_rel_num) / (Ei * Ej) : 0.0
                        if v_rel == 0.0
                            continue
                        end
                        σ = get_sigma(cs_cache, s, quark_params, thermo_params, K_coeffs; n_points=n_sigma_points)
                        ω += w_pi * w_pj * w_cθi * w_cθj * wφ * (p_i^2) * (p_j^2) * f_i * f_j * v_rel * σ
                    end
                end
            end
        end
    end

    return prefactor * ω
end

end # module
