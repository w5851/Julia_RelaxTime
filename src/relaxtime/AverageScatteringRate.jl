module AverageScatteringRate

"""
平均散射率计算模块（各向异性可选）。
- 积分采用 Gauss-Legendre：动量 32 节点，角度 4 节点（可覆盖 7 阶多项式）。
- 支持各向异性分布 `quark_distribution_aniso(..., ξ, cosθ)`；`ξ=0` 时退化为各向同性。
- 散射截面使用 `TotalCrossSection.total_cross_section`，可按需预计算并插值。

Cross-section cache notes:
- `CrossSectionCache(process; compute_missing=true, rtol=1e-2, max_refine=12)` 默认会在需要时调用
    `TotalCrossSection.total_cross_section` 来补点（并在局部做自适应细分以提升插值可靠性）。
- 若你只想用“手动插入/预计算的 σ(s) 点”做线性插值（例如单元测试、性能 smoke，或用常数 σ
    进行回归测试），可用 `compute_missing=false`：此时 `get_sigma` **不会**触发任何真实截面计算。
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
using .PNJLQuarkDistributions: quark_distribution, antiquark_distribution
using .PNJLQuarkDistributions_Aniso: quark_distribution_aniso, antiquark_distribution_aniso

export average_scattering_rate, CrossSectionCache, precompute_cross_section!

const DEFAULT_P_NODES = 6
const DEFAULT_ANGLE_NODES = 4  # cosθ节点数
const DEFAULT_PHI_NODES = 8    # φ节点数
const DQ = 6.0 # 简并度d_q=2*N_c=6
const TWO_PI = 2.0 * π

# --------------------------- 工具函数 ---------------------------

@inline function is_antiquark(flavor::Symbol)::Bool
    return flavor === :ubar || flavor === :dbar || flavor === :sbar
end

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
    # Convention: always return the quark chemical potential μ_q (positive sign).
    # The particle/antiparticle distinction is handled by using
    # `quark_distribution*` vs `antiquark_distribution*`.
    if flavor in [:u, :ubar]; return quark_params.μ.u
    elseif flavor in [:d, :dbar]; return quark_params.μ.d
    elseif flavor in [:s, :sbar]; return quark_params.μ.s
    else; error("Unknown flavor $flavor") end
end

# -------------------- 截面缓存与插值 --------------------
mutable struct CrossSectionCache
    process::Symbol
    s_vals::Vector{Float64}
    sigma_vals::Vector{Float64}

    # If false, `get_sigma` will only interpolate/clamp from existing points
    # and will not call the expensive exact cross-section calculation.
    compute_missing::Bool

    # Adaptive interpolation controls
    rtol::Float64
    max_refine::Int

    # If true, disable interpolation entirely and always compute exact σ(s)
    # at the queried s (while still memoizing exact points in the cache).
    direct::Bool

    # If false, never insert new σ(s) points (useful for direct mode where
    # queries are effectively unique and caching just wastes memory).
    memoize::Bool
end

function CrossSectionCache(
    process::Symbol;
    compute_missing::Bool=true,
    rtol::Float64=1e-2,
    max_refine::Int=12,
    direct::Bool=false,
    memoize::Bool=true,
)
    return CrossSectionCache(process, Float64[], Float64[], compute_missing, rtol, max_refine, direct, memoize)
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
    elseif s == cache.s_vals[1]
        return cache.sigma_vals[1]
    elseif s == cache.s_vals[end]
        return cache.sigma_vals[end]
    elseif s < cache.s_vals[1] || s > cache.s_vals[end]
        return nothing
    else
        idx = searchsortedfirst(cache.s_vals, s)
        s1, s2 = cache.s_vals[idx-1], cache.s_vals[idx]
        σ1, σ2 = cache.sigma_vals[idx-1], cache.sigma_vals[idx]
        w = (s - s1)/(s2 - s1)
        return σ1 + w*(σ2 - σ1)
    end
end

@inline function relerr(a::Float64, b::Float64)
    return abs(a - b) / max(1e-12, abs(b))
end

function sigma_at!(cache::CrossSectionCache, s::Float64,
    quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple;
    n_points::Int)

    idx = searchsortedfirst(cache.s_vals, s)
    if idx <= length(cache.s_vals) && cache.s_vals[idx] == s
        return cache.sigma_vals[idx]
    end

    try
        σ = total_cross_section(cache.process, s, quark_params, thermo_params, K_coeffs; n_points=n_points)
        if !isfinite(σ) || σ < 0.0
            @warn "sigma_at! produced non-finite/negative sigma; using 0" process=cache.process s=s sigma=σ
            return 0.0
        end
        if cache.memoize
            insert_sigma!(cache, s, σ)
        end
        return σ
    catch err
        @warn "sigma_at! failed; using 0" process=cache.process s=s error=err
        return 0.0
    end
end

function adaptive_interpolate_sigma!(cache::CrossSectionCache, s::Float64,
    quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple;
    n_points::Int)

    n = length(cache.s_vals)
    n == 0 && return sigma_at!(cache, s, quark_params, thermo_params, K_coeffs; n_points=n_points)

    # Always compute exact values outside the current cache window
    if s <= cache.s_vals[1] || s >= cache.s_vals[end]
        return sigma_at!(cache, s, quark_params, thermo_params, K_coeffs; n_points=n_points)
    end

    idx = searchsortedfirst(cache.s_vals, s)
    if idx <= length(cache.s_vals) && cache.s_vals[idx] == s
        return cache.sigma_vals[idx]
    end

    # Bracket s by cached neighbors
    sL = cache.s_vals[idx-1]
    sR = cache.s_vals[idx]
    σL = cache.sigma_vals[idx-1]
    σR = cache.sigma_vals[idx]

    # Refine only along the branch that contains the query s.
    # Criterion: midpoint deviation from linear interpolation within rtol.
    for _ in 1:cache.max_refine
        sM = 0.5 * (sL + sR)
        σM = sigma_at!(cache, sM, quark_params, thermo_params, K_coeffs; n_points=n_points)
        σM_lin = 0.5 * (σL + σR)
        if relerr(σM, σM_lin) <= cache.rtol
            break
        end
        if s < sM
            sR = sM
            σR = σM
        else
            sL = sM
            σL = σM
        end
    end

    w = (s - sL) / (sR - sL)
    return σL + w * (σR - σL)
end

function get_sigma(cache::CrossSectionCache, s::Float64,
    quark_params::NamedTuple, thermo_params::NamedTuple, K_coeffs::NamedTuple;
    n_points::Int=TotalCrossSection.DEFAULT_T_INTEGRAL_POINTS)
    if cache.direct
        if cache.memoize
            return sigma_at!(cache, s, quark_params, thermo_params, K_coeffs; n_points=n_points)
        end
        # Direct exact evaluation without storing points.
        try
            σ = total_cross_section(cache.process, s, quark_params, thermo_params, K_coeffs; n_points=n_points)
            if !isfinite(σ) || σ < 0.0
                @warn "direct sigma produced non-finite/negative sigma; using 0" process=cache.process s=s sigma=σ
                return 0.0
            end
            return σ
        catch err
            @warn "direct sigma failed; using 0" process=cache.process s=s error=err
            return 0.0
        end
    end
    if !cache.compute_missing
        n = length(cache.s_vals)
        n == 0 && error("CrossSectionCache has no points; cannot interpolate")
        if n == 1
            return cache.sigma_vals[1]
        elseif s <= cache.s_vals[1]
            return cache.sigma_vals[1]
        elseif s >= cache.s_vals[end]
            return cache.sigma_vals[end]
        else
            # Pure linear interpolation between cached neighbors
            idx = searchsortedfirst(cache.s_vals, s)
            s1, s2 = cache.s_vals[idx-1], cache.s_vals[idx]
            σ1, σ2 = cache.sigma_vals[idx-1], cache.sigma_vals[idx]
            w = (s - s1) / (s2 - s1)
            return σ1 + w * (σ2 - σ1)
        end
    end

    # Strategy:
    # - If cache already has points, use adaptive local refinement to improve
    #   linear interpolation accuracy near sharp features (thresholds/resonances).
    # - Otherwise compute exactly and insert.
    return adaptive_interpolate_sigma!(cache, s, quark_params, thermo_params, K_coeffs; n_points=n_points)
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
    
    # 使用 [0, 1) 上的Gauss-Legendre节点，通过变换映射到 [0, ∞)
    # 变换: p = scale * t / (1-t), dp/dt = scale / (1-t)^2
    t_grid, t_w = gauleg(0.0, 1.0, p_nodes)
    cos_grid === nothing && ((cos_grid, cos_w) = gauleg(0.0, 1.0, angle_nodes))
    
    integral = 0.0
    for (t, wt) in zip(t_grid, t_w)
        # 避免 t=1 的奇点
        if t >= 0.9999
            continue
        end
        # 变换到半无穷区间
        p = scale * t / (1.0 - t)
        dp_dt = scale / (1.0 - t)^2
        
        for (cθ, wθ) in zip(cos_grid, cos_w)
            f = distribution_with_anisotropy(flavor, p, m, μ, T, Φ, Φbar, ξ, cθ)
            integral += wt * wθ * p^2 * f * dp_dt
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
- `quark_params::NamedTuple`: 夸克参数 (m, μ, A)
- `thermo_params::NamedTuple`: 热力学参数 (T, Φ, Φbar, ξ)
- `K_coeffs::NamedTuple`: 有效耦合系数

# Keyword Arguments
- `p_nodes::Int`: 动量积分节点数 (默认6)
- `angle_nodes::Int`: 角度积分节点数 (默认2)
- `phi_nodes::Int`: 方位角积分节点数 (默认4)
- `scale::Float64`: 半无穷积分尺度参数
"""
function average_scattering_rate(
    process::Symbol,
    quark_params::NamedTuple,
    thermo_params::NamedTuple,
    K_coeffs::NamedTuple;
    p_nodes::Int=DEFAULT_P_NODES,
    angle_nodes::Int=DEFAULT_ANGLE_NODES,
    phi_nodes::Int=DEFAULT_PHI_NODES,
    p_grid::Union{Nothing,Vector{Float64}}=nothing, p_w::Union{Nothing,Vector{Float64}}=nothing,
    cos_grid::Union{Nothing,Vector{Float64}}=nothing, cos_w::Union{Nothing,Vector{Float64}}=nothing,
    phi_grid::Union{Nothing,Vector{Float64}}=nothing, phi_w::Union{Nothing,Vector{Float64}}=nothing,
    cs_cache::CrossSectionCache=CrossSectionCache(process),
    n_sigma_points::Int=TotalCrossSection.DEFAULT_T_INTEGRAL_POINTS,
    scale::Float64=DEFAULT_SEMI_INF_SCALE
)::Float64
    # 解析粒子、质量、化学势
    pi_sym, pj_sym, pc_sym, pd_sym = parse_particles_from_process(process)
    mi = get_mass(pi_sym, quark_params); mj = get_mass(pj_sym, quark_params)
    mc = get_mass(pc_sym, quark_params); md = get_mass(pd_sym, quark_params)
    μi = get_mu(pi_sym, quark_params); μj = get_mu(pj_sym, quark_params)
    μc = get_mu(pc_sym, quark_params); μd = get_mu(pd_sym, quark_params)

    T = thermo_params.T; Φ = thermo_params.Φ; Φbar = thermo_params.Φbar
    ξ = hasproperty(thermo_params, :ξ) ? thermo_params.ξ : 0.0

    # 使用半无穷积分 [0, ∞)
    return _average_scattering_rate_semi_infinite(
        process, pi_sym, pj_sym, mi, mj, μi, μj, T, Φ, Φbar, ξ,
        quark_params, thermo_params, K_coeffs,
        p_nodes, angle_nodes, phi_nodes, scale,
        p_grid, p_w, cos_grid, cos_w, phi_grid, phi_w,
        cs_cache, n_sigma_points
    )
end

# 半无穷积分版本（保留原有实现）
function _average_scattering_rate_semi_infinite(
    process, pi_sym, pj_sym, mi, mj, μi, μj, T, Φ, Φbar, ξ,
    quark_params, thermo_params, K_coeffs,
    p_nodes, angle_nodes, phi_nodes, scale,
    p_grid, p_w, cos_grid, cos_w, phi_grid, phi_w,
    cs_cache, n_sigma_points
)
    # 使用 [0, 1) 上的Gauss-Legendre节点，通过变换映射到 [0, ∞)
    t_grid, t_w = gauleg(0.0, 1.0, p_nodes)
    # 使用完整积分区间 [-1,1] 和 [0,2π]
    cos_grid === nothing && ((cos_grid, cos_w) = gauleg(-1.0, 1.0, angle_nodes))
    phi_grid === nothing && ((phi_grid, phi_w) = gauleg(0.0, TWO_PI, phi_nodes))

    # 数密度（用于归一化）- 使用半无穷积分
    ρ_i = number_density(pi_sym, mi, μi, T, Φ, Φbar, ξ; p_nodes=DEFAULT_SEMI_INF_NODES, angle_nodes=angle_nodes, scale=scale)
    ρ_j = number_density(pj_sym, mj, μj, T, Φ, Φbar, ξ; p_nodes=DEFAULT_SEMI_INF_NODES, angle_nodes=angle_nodes, scale=scale)
    if ρ_i == 0.0 || ρ_j == 0.0
        return 0.0
    end

    # 原始公式的前因子: DQ²/(2π)⁵ = DQ²/(32π⁵)
    prefactor = (DQ^2) / (32.0 * π^5 * ρ_i * ρ_j)
    ω = 0.0

    for (t_i, w_ti) in zip(t_grid, t_w)
        # 避免 t=1 的奇点
        if t_i >= 0.9999
            continue
        end
        # 变换到半无穷区间
        p_i = scale * t_i / (1.0 - t_i)
        dp_i_dt = scale / (1.0 - t_i)^2
        Ei = energy_from_p(p_i, mi)
        
        for (t_j, w_tj) in zip(t_grid, t_w)
            if t_j >= 0.9999
                continue
            end
            p_j = scale * t_j / (1.0 - t_j)
            dp_j_dt = scale / (1.0 - t_j)^2
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
                        s = mi^2 + mj^2 + 2.0 * (Ei * Ej - p_i * p_j * cosΘ)
                        if s <= (mi + mj)^2
                            continue
                        end
                        # 使用质心系能量计算v_rel (参考公式文档)
                        s_rt = sqrt(s)
                        Ei_cm = (s + mi^2 - mj^2) / (2.0 * s_rt)
                        Ej_cm = (s - mi^2 + mj^2) / (2.0 * s_rt)
                        pi_cm = sqrt(max(0.0, (s - (mi + mj)^2) * (s - (mi - mj)^2))) / (2.0 * s_rt)
                        pj_cm = pi_cm  # 质心系中动量大小相等
                        v_rel_num = (Ei_cm * Ej_cm + pi_cm * pj_cm)^2 - (mi * mj)^2
                        v_rel = v_rel_num > 0.0 ? sqrt(v_rel_num) / (Ei_cm * Ej_cm) : 0.0
                        if v_rel == 0.0 || v_rel > 2.0
                            continue
                        end
                        σ = get_sigma(cs_cache, s, quark_params, thermo_params, K_coeffs; n_points=n_sigma_points)
                        # 注意：需要乘以雅可比行列式 dp_i_dt * dp_j_dt
                        ω += w_ti * w_tj * w_cθi * w_cθj * wφ * (p_i^2) * (p_j^2) * f_i * f_j * v_rel * σ * dp_i_dt * dp_j_dt
                    end
                end
            end
        end
    end

    return prefactor * ω
end

end # module
