module TransportCoefficients

"""
在弛豫时间近似（RTA）下计算夸克物质输运系数：
- 剪切粘滞系数 η
- 体粘滞系数 ζ
- 电导率 σ

单位约定（自然单位 ħ=c=k_B=1）：
- 温度/化学势/质量/动量均使用 fm⁻¹
- τ 使用 fm

实现约定：
- τ 与动量 p 无关（直接使用 `RelaxationTime.relaxation_times` 的输出）
- ξ=0：各向同性，仅对动量 p 积分（乘以 4π）
- ξ≠0：分布函数采用 RS 形式 `*_distribution_aniso(p,...,ξ,cosθ)`，对 cosθ 再做一层积分（φ 角已积掉，乘以 2π）
"""

include("../Constants_PNJL.jl")
include("../integration/GaussLegendre.jl")
include("../QuarkDistribution.jl")
include("../QuarkDistribution_Aniso.jl")

using Base.MathConstants: π

using .Constants_PNJL: N_color
using .GaussLegendre: gauleg, DEFAULT_MOMENTUM_NODES, DEFAULT_MOMENTUM_WEIGHTS, DEFAULT_COSΘ_NODES, DEFAULT_COSΘ_WEIGHTS
using .PNJLQuarkDistributions: quark_distribution, antiquark_distribution
using .PNJLQuarkDistributions_Aniso: quark_distribution_aniso, antiquark_distribution_aniso

export shear_viscosity, bulk_viscosity, electric_conductivity, transport_coefficients

const TWO_PI = 2.0 * π

@inline energy_from_p(p::Float64, m::Float64) = sqrt(p * p + m * m)

@inline function degeneracy_default()::Float64
    return 2.0 * Float64(N_color)
end

@inline function default_charges()
    return (u = 2.0 / 3.0, d = -1.0 / 3.0, s = -1.0 / 3.0)
end

@inline function fermi_factor(f::Float64)
    return f * (1.0 - f)
end

@inline function distribution_for_species(
    species::Symbol,
    p::Float64,
    m::Float64,
    μ::Float64,
    T::Float64,
    Φ::Float64,
    Φbar::Float64,
    ξ::Float64,
    cosθ::Float64
)
    if ξ == 0.0
        E = energy_from_p(p, m)
        if species in (:u, :d, :s)
            return quark_distribution(E, μ, T, Φ, Φbar)
        else
            return antiquark_distribution(E, μ, T, Φ, Φbar)
        end
    else
        if species in (:u, :d, :s)
            return quark_distribution_aniso(p, m, μ, T, Φ, Φbar, ξ, cosθ)
        else
            return antiquark_distribution_aniso(p, m, μ, T, Φ, Φbar, ξ, cosθ)
        end
    end
end

@inline function mass_for_species(species::Symbol, quark_params::NamedTuple)::Float64
    if species in (:u, :d)
        return species === :u ? quark_params.m.u : quark_params.m.d
    elseif species === :s
        return quark_params.m.s
    elseif species in (:ubar, :dbar)
        return species === :ubar ? quark_params.m.u : quark_params.m.d
    elseif species === :sbar
        return quark_params.m.s
    else
        error("Unknown species: $species")
    end
end

@inline function mu_for_species(species::Symbol, quark_params::NamedTuple)::Float64
    # 注意：这里的 μ 取“味道的夸克化学势”本身（正值或按调用方约定），
    # 反粒子通过分布函数 `antiquark_distribution` 处理 E+μ 结构。
    if species in (:u, :ubar)
        return quark_params.μ.u
    elseif species in (:d, :dbar)
        return quark_params.μ.d
    elseif species in (:s, :sbar)
        return quark_params.μ.s
    else
        error("Unknown species: $species")
    end
end

@inline function tau_for_species(species::Symbol, tau::NamedTuple)::Float64
    hasproperty(tau, species) || error("tau is missing :$species")
    return getproperty(tau, species)
end

function _p_nodes_weights(p_nodes::Int, p_max::Float64, p_grid, p_w)
    if p_grid !== nothing
        p_w !== nothing || error("p_w must be provided when p_grid is provided")
        return p_grid, p_w
    end
    return gauleg(0.0, p_max, p_nodes)
end

function _cos_nodes_weights(cos_nodes::Int, cos_grid, cos_w)
    if cos_grid !== nothing
        cos_w !== nothing || error("cos_w must be provided when cos_grid is provided")
        return cos_grid, cos_w
    end
    return gauleg(-1.0, 1.0, cos_nodes)
end

"""剪切粘滞系数 η（RTA）。"""
function shear_viscosity(
    quark_params::NamedTuple,
    thermo_params::NamedTuple;
    tau::NamedTuple,
    degeneracy::Float64=degeneracy_default(),
    p_nodes::Int=length(DEFAULT_MOMENTUM_NODES),
    p_max::Float64=10.0,
    p_grid::Union{Nothing,Vector{Float64}}=nothing,
    p_w::Union{Nothing,Vector{Float64}}=nothing,
    cos_nodes::Int=length(DEFAULT_COSΘ_NODES),
    cos_grid::Union{Nothing,Vector{Float64}}=nothing,
    cos_w::Union{Nothing,Vector{Float64}}=nothing
)::Float64
    T = thermo_params.T
    Φ = thermo_params.Φ
    Φbar = thermo_params.Φbar
    ξ = get(thermo_params, :ξ, 0.0)

    nodes_p, weights_p = _p_nodes_weights(p_nodes, p_max, p_grid, p_w)

    pref_measure_iso = 4.0 * π / (2.0 * π)^3
    pref_measure_aniso = 2.0 * π / (2.0 * π)^3

    species_list = (:u, :d, :s, :ubar, :dbar, :sbar)

    if ξ == 0.0
        acc = 0.0
        @inbounds for (p, wp) in zip(nodes_p, weights_p)
            p2 = p * p
            for sp in species_list
                m = mass_for_species(sp, quark_params)
                μ = mu_for_species(sp, quark_params)
                E = energy_from_p(p, m)
                f = distribution_for_species(sp, p, m, μ, T, Φ, Φbar, 0.0, 0.0)
                ff = fermi_factor(f)
                τ = tau_for_species(sp, tau)
                acc += wp * (p2 * p2) / (E * E) * (degeneracy * τ * ff)
            end
        end
        integral = pref_measure_iso * acc
        return (1.0 / (15.0 * T)) * integral
    else
        nodes_cos, weights_cos = _cos_nodes_weights(cos_nodes, cos_grid, cos_w)
        acc = 0.0
        @inbounds for (p, wp) in zip(nodes_p, weights_p)
            p2 = p * p
            for (c, wc) in zip(nodes_cos, weights_cos)
                for sp in species_list
                    m = mass_for_species(sp, quark_params)
                    μ = mu_for_species(sp, quark_params)
                    E = energy_from_p(p, m)
                    f = distribution_for_species(sp, p, m, μ, T, Φ, Φbar, ξ, c)
                    ff = fermi_factor(f)
                    τ = tau_for_species(sp, tau)
                    acc += wp * wc * (p2 * p2) / (E * E) * (degeneracy * τ * ff)
                end
            end
        end
        integral = pref_measure_aniso * acc
        return (1.0 / (15.0 * T)) * integral
    end
end

"""电导率 σ（RTA）。"""
function electric_conductivity(
    quark_params::NamedTuple,
    thermo_params::NamedTuple;
    tau::NamedTuple,
    charges::NamedTuple=default_charges(),
    degeneracy::Float64=degeneracy_default(),
    p_nodes::Int=length(DEFAULT_MOMENTUM_NODES),
    p_max::Float64=10.0,
    p_grid::Union{Nothing,Vector{Float64}}=nothing,
    p_w::Union{Nothing,Vector{Float64}}=nothing,
    cos_nodes::Int=length(DEFAULT_COSΘ_NODES),
    cos_grid::Union{Nothing,Vector{Float64}}=nothing,
    cos_w::Union{Nothing,Vector{Float64}}=nothing
)::Float64
    T = thermo_params.T
    Φ = thermo_params.Φ
    Φbar = thermo_params.Φbar
    ξ = get(thermo_params, :ξ, 0.0)

    nodes_p, weights_p = _p_nodes_weights(p_nodes, p_max, p_grid, p_w)

    pref_measure_iso = 4.0 * π / (2.0 * π)^3
    pref_measure_aniso = 2.0 * π / (2.0 * π)^3

    function q2_for_species(sp::Symbol)
        if sp in (:u, :ubar)
            return charges.u^2
        elseif sp in (:d, :dbar)
            return charges.d^2
        elseif sp in (:s, :sbar)
            return charges.s^2
        else
            error("Unknown species: $sp")
        end
    end

    species_list = (:u, :d, :s, :ubar, :dbar, :sbar)

    if ξ == 0.0
        acc = 0.0
        @inbounds for (p, wp) in zip(nodes_p, weights_p)
            p2 = p * p
            for sp in species_list
                m = mass_for_species(sp, quark_params)
                μ = mu_for_species(sp, quark_params)
                E = energy_from_p(p, m)
                f = distribution_for_species(sp, p, m, μ, T, Φ, Φbar, 0.0, 0.0)
                ff = fermi_factor(f)
                τ = tau_for_species(sp, tau)
                acc += wp * p2 * (q2_for_species(sp) / (E * E)) * (degeneracy * τ * ff)
            end
        end
        integral = pref_measure_iso * acc
        return (2.0 / (3.0 * T)) * integral
    else
        nodes_cos, weights_cos = _cos_nodes_weights(cos_nodes, cos_grid, cos_w)
        acc = 0.0
        @inbounds for (p, wp) in zip(nodes_p, weights_p)
            p2 = p * p
            for (c, wc) in zip(nodes_cos, weights_cos)
                for sp in species_list
                    m = mass_for_species(sp, quark_params)
                    μ = mu_for_species(sp, quark_params)
                    E = energy_from_p(p, m)
                    f = distribution_for_species(sp, p, m, μ, T, Φ, Φbar, ξ, c)
                    ff = fermi_factor(f)
                    τ = tau_for_species(sp, tau)
                    acc += wp * wc * p2 * (q2_for_species(sp) / (E * E)) * (degeneracy * τ * ff)
                end
            end
        end
        integral = pref_measure_aniso * acc
        return (2.0 / (3.0 * T)) * integral
    end
end

"""体粘滞系数 ζ（RTA）。

注意：此实现按 docs/reference/formula/输运系数by弛豫时间.md 给出的表达式直译。
输入 `bulk_coeffs` 建议使用 `PNJL.ThermoDerivatives.bulk_derivative_coeffs` 的返回值。
"""
function bulk_viscosity(
    quark_params::NamedTuple,
    thermo_params::NamedTuple;
    tau::NamedTuple,
    bulk_coeffs::NamedTuple,
    degeneracy::Float64=degeneracy_default(),
    p_nodes::Int=length(DEFAULT_MOMENTUM_NODES),
    p_max::Float64=10.0,
    p_grid::Union{Nothing,Vector{Float64}}=nothing,
    p_w::Union{Nothing,Vector{Float64}}=nothing,
    cos_nodes::Int=length(DEFAULT_COSΘ_NODES),
    cos_grid::Union{Nothing,Vector{Float64}}=nothing,
    cos_w::Union{Nothing,Vector{Float64}}=nothing
)::Float64
    T = thermo_params.T
    Φ = thermo_params.Φ
    Φbar = thermo_params.Φbar
    ξ = get(thermo_params, :ξ, 0.0)

    dP_depsilon_n = bulk_coeffs.dP_depsilon_n
    dP_dn_epsilon = bulk_coeffs.dP_dn_epsilon
    dM_dT = bulk_coeffs.dM_dT
    dM_dmu = bulk_coeffs.dM_dmu

    nodes_p, weights_p = _p_nodes_weights(p_nodes, p_max, p_grid, p_w)

    pref_measure_iso = 4.0 * π / (2.0 * π)^3
    pref_measure_aniso = 2.0 * π / (2.0 * π)^3

    # 这里按三味道质量与导数来使用（与 ThermoDerivatives.mass_derivatives 一致：SVector{3}）
    function flavor_index(sp::Symbol)
        if sp in (:u, :d, :ubar, :dbar)
            return 1
        elseif sp in (:s, :sbar)
            return 3
        else
            error("Unknown species: $sp")
        end
    end

    # 只按文档里的求和 a=u,d,s；把反粒子项显式加入
    flavors = (:u, :d, :s)

    function one_flavor_pair_contrib(flavor::Symbol, p::Float64, cosθ::Float64)
        sp_q = flavor
        sp_aq = Symbol(string(flavor, "bar"))

        m = mass_for_species(sp_q, quark_params)
        μ = mu_for_species(sp_q, quark_params)
        E = energy_from_p(p, m)

        idx = flavor_index(sp_q)
        dE_dT_val = (m / E) * dM_dT[idx]
        dE_dmu_val = (m / E) * dM_dmu[idx]

        f_q = distribution_for_species(sp_q, p, m, μ, T, Φ, Φbar, ξ, cosθ)
        f_aq = distribution_for_species(sp_aq, p, m, μ, T, Φ, Φbar, ξ, cosθ)

        w_q = degeneracy * tau_for_species(sp_q, tau) * fermi_factor(f_q)
        w_aq = degeneracy * tau_for_species(sp_aq, tau) * fermi_factor(f_aq)

        bracket = (p * p) / (3.0 * E) - dP_depsilon_n * (E - T * dE_dT_val - μ * dE_dmu_val) + dP_dn_epsilon * dE_dmu_val

        term_sum = (m * m / E) * (w_q + w_aq) * bracket
        term_diff = (m * m / E) * (w_q - w_aq) * dP_dn_epsilon

        return term_sum - term_diff
    end

    if ξ == 0.0
        acc = 0.0
        @inbounds for (p, wp) in zip(nodes_p, weights_p)
            p2 = p * p
            for fl in flavors
                acc += wp * p2 * one_flavor_pair_contrib(fl, p, 0.0)
            end
        end
        integral = pref_measure_iso * acc
        return -(1.0 / (3.0 * T)) * integral
    else
        nodes_cos, weights_cos = _cos_nodes_weights(cos_nodes, cos_grid, cos_w)
        acc = 0.0
        @inbounds for (p, wp) in zip(nodes_p, weights_p)
            p2 = p * p
            for (c, wc) in zip(nodes_cos, weights_cos)
                for fl in flavors
                    acc += wp * wc * p2 * one_flavor_pair_contrib(fl, p, c)
                end
            end
        end
        integral = pref_measure_aniso * acc
        return -(1.0 / (3.0 * T)) * integral
    end
end

"""一次性计算 (η, ζ, σ)。"""
function transport_coefficients(
    quark_params::NamedTuple,
    thermo_params::NamedTuple;
    tau::NamedTuple,
    bulk_coeffs::Union{Nothing,NamedTuple}=nothing,
    charges::NamedTuple=default_charges(),
    degeneracy::Float64=degeneracy_default(),
    p_nodes::Int=length(DEFAULT_MOMENTUM_NODES),
    p_max::Float64=10.0,
    p_grid::Union{Nothing,Vector{Float64}}=nothing,
    p_w::Union{Nothing,Vector{Float64}}=nothing,
    cos_nodes::Int=length(DEFAULT_COSΘ_NODES),
    cos_grid::Union{Nothing,Vector{Float64}}=nothing,
    cos_w::Union{Nothing,Vector{Float64}}=nothing
)::NamedTuple
    eta = shear_viscosity(
        quark_params,
        thermo_params;
        tau=tau,
        degeneracy=degeneracy,
        p_nodes=p_nodes,
        p_max=p_max,
        p_grid=p_grid,
        p_w=p_w,
        cos_nodes=cos_nodes,
        cos_grid=cos_grid,
        cos_w=cos_w,
    )

    sigma = electric_conductivity(
        quark_params,
        thermo_params;
        tau=tau,
        charges=charges,
        degeneracy=degeneracy,
        p_nodes=p_nodes,
        p_max=p_max,
        p_grid=p_grid,
        p_w=p_w,
        cos_nodes=cos_nodes,
        cos_grid=cos_grid,
        cos_w=cos_w,
    )

    zeta = bulk_coeffs === nothing ? NaN : bulk_viscosity(
        quark_params,
        thermo_params;
        tau=tau,
        bulk_coeffs=bulk_coeffs,
        degeneracy=degeneracy,
        p_nodes=p_nodes,
        p_max=p_max,
        p_grid=p_grid,
        p_w=p_w,
        cos_nodes=cos_nodes,
        cos_grid=cos_grid,
        cos_w=cos_w,
    )

    return (eta=eta, zeta=zeta, sigma=sigma)
end

end # module
