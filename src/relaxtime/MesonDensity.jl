raw"""
MesonDensity

介子数密度最小实现入口。当前阶段仅覆盖稳定粒子极限：

- 玻色分布函数
- `π/K` 聚合通道默认简并因子
- 稳定粒子极限数密度
- `K/π` 比值与最小温度扫描

后续 BU / BW / 各向异性扩展在此基础上演进。
"""
module MesonDensity

using ..GaussLegendre: gauleg
using ..AFieldBuilder: ensure_quark_params_has_A
using ..EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using ..PolarizationAniso: polarization_with_width
using ..MesonPropagator: meson_propagator_simple
using Main.Constants_PNJL: G_fm2, K_fm5

export DEFAULT_MESON_DENSITY_Q_NODES
export DEFAULT_PHASE_SHIFT_Q_MAX, DEFAULT_PHASE_SHIFT_Q_NODES
export DEFAULT_PHASE_SHIFT_OMEGA_MAX, DEFAULT_PHASE_SHIFT_OMEGA_NODES
export bose_distribution, meson_degeneracy
export stable_meson_number_density, stable_kpi_ratio, stable_kpi_scan
export strict_bw_meson_number_density, strict_bw_meson_density_summary
export phase_shift_meson_number_density, phase_shift_meson_density_summary

const DEFAULT_MESON_DENSITY_Q_NODES = 256
const DEFAULT_PHASE_SHIFT_Q_MAX = 12.0
const DEFAULT_PHASE_SHIFT_Q_NODES = 48
const DEFAULT_PHASE_SHIFT_OMEGA_MAX = 10.0
const DEFAULT_PHASE_SHIFT_OMEGA_NODES = 48
const DEFAULT_PHASE_SHIFT_ETA = 1e-6

@inline function _require_nonnegative(name::AbstractString, value::Float64)
    value >= 0.0 && return
    throw(ArgumentError("$(name) must be nonnegative, got $(value)"))
end

@inline function _default_qmax(mass::Float64, T::Float64, μ::Float64)::Float64
    gap = max(mass - μ, 0.0)
    return max(8.0, 20.0 * T + 10.0 * gap, 8.0 * mass + 10.0 * T)
end

"""
    meson_degeneracy(meson; charge_resolved=false) -> Int

返回当前主线下的 `π/K` 简并因子。

- 聚合通道（默认）：`d_π = 3`、`d_K = 4`
- 电荷分辨通道：`d = 1`
"""
@inline function meson_degeneracy(meson::Symbol; charge_resolved::Bool=false)::Int
    if meson === :pi
        return charge_resolved ? 1 : 3
    elseif meson === :K
        return charge_resolved ? 1 : 4
    else
        throw(ArgumentError("Unsupported meson $(meson). Use :pi or :K."))
    end
end

raw"""
    bose_distribution(E, μ, T) -> Float64

玻色分布函数：

```math
g(E) = 1 / (\exp((E-\mu)/T) - 1)
```
"""
function bose_distribution(E::Float64, μ::Float64, T::Float64)::Float64
    _require_nonnegative("temperature T", T)
    T == 0.0 && return 0.0
    E > μ || throw(ArgumentError("Bose distribution requires E > μ to avoid pole, got E=$(E), μ=$(μ)"))

    exponent = (E - μ) / T
    exponent > 700.0 && return 0.0
    return 1.0 / expm1(exponent)
end

raw"""
    stable_meson_number_density(mass, T; μ=0.0, degeneracy=1,
                                qmax=nothing, num_q_nodes=DEFAULT_MESON_DENSITY_Q_NODES) -> Float64

稳定粒子极限介子数密度：

```math
n_M = d_M \\int_0^\\infty \\frac{dq\\,q^2}{2\\pi^2}
      \\frac{1}{\\exp((E_M-\\mu_M)/T)-1},
\\qquad
E_M = \\sqrt{q^2 + m_M^2}.
```
"""
function stable_meson_number_density(mass::Float64, T::Float64;
                                     μ::Float64=0.0,
                                     degeneracy::Integer=1,
                                     qmax::Union{Nothing,Float64}=nothing,
                                     num_q_nodes::Int=DEFAULT_MESON_DENSITY_Q_NODES)::Float64
    _require_nonnegative("mass", mass)
    _require_nonnegative("temperature T", T)
    degeneracy > 0 || throw(ArgumentError("degeneracy must be positive, got $(degeneracy)"))
    num_q_nodes > 1 || throw(ArgumentError("num_q_nodes must be > 1, got $(num_q_nodes)"))
    mass > μ || throw(ArgumentError("Stable boson density requires mass > μ, got mass=$(mass), μ=$(μ)"))
    T == 0.0 && return 0.0

    q_upper = qmax === nothing ? _default_qmax(mass, T, μ) : Float64(qmax)
    q_upper > 0.0 || throw(ArgumentError("qmax must be positive, got $(q_upper)"))

    nodes, weights = gauleg(0.0, q_upper, num_q_nodes)
    integral = 0.0
    @inbounds for i in eachindex(nodes, weights)
        q = nodes[i]
        E = hypot(q, mass)
        integral += weights[i] * q^2 * bose_distribution(E, μ, T)
    end
    return degeneracy * integral / (2.0 * π^2)
end

"""
    stable_kpi_ratio(m_pi, m_K, T; μ_pi=0.0, μ_K=0.0, d_pi=3, d_K=4, kwargs...) -> Float64

返回稳定粒子极限的 `K/π` 数密度比值。
"""
function stable_kpi_ratio(m_pi::Float64, m_K::Float64, T::Float64;
                          μ_pi::Float64=0.0, μ_K::Float64=0.0,
                          d_pi::Integer=3, d_K::Integer=4,
                          kwargs...)::Float64
    n_pi = stable_meson_number_density(m_pi, T; μ=μ_pi, degeneracy=d_pi, kwargs...)
    n_K = stable_meson_number_density(m_K, T; μ=μ_K, degeneracy=d_K, kwargs...)
    return iszero(n_pi) ? NaN : n_K / n_pi
end

"""
    stable_kpi_scan(temperatures; m_pi, m_K, μ_pi=0.0, μ_K=0.0, d_pi=3, d_K=4, kwargs...)
        -> NamedTuple

对一组温度执行稳定粒子极限 `π/K` 数密度扫描，返回：

- `temperatures`
- `n_pi`
- `n_K`
- `kpi_ratio`
"""
function stable_kpi_scan(temperatures::AbstractVector{<:Real};
                         m_pi::Float64, m_K::Float64,
                         μ_pi::Float64=0.0, μ_K::Float64=0.0,
                         d_pi::Integer=3, d_K::Integer=4,
                         kwargs...)
    Ts = Float64[Float64(T) for T in temperatures]
    n_pi = Vector{Float64}(undef, length(Ts))
    n_K = Vector{Float64}(undef, length(Ts))
    ratios = Vector{Float64}(undef, length(Ts))

    for i in eachindex(Ts)
        Ti = Ts[i]
        n_pi[i] = stable_meson_number_density(m_pi, Ti; μ=μ_pi, degeneracy=d_pi, kwargs...)
        n_K[i] = stable_meson_number_density(m_K, Ti; μ=μ_K, degeneracy=d_K, kwargs...)
        ratios[i] = iszero(n_pi[i]) ? NaN : n_K[i] / n_pi[i]
    end

    return (temperatures=Ts, n_pi=n_pi, n_K=n_K, kpi_ratio=ratios)
end

@inline function _strict_bw_kernel(Δω::Float64, gamma::Float64)::Float64
    half_gamma = gamma / 2.0
    return half_gamma / (Δω^2 + half_gamma^2)
end

raw"""
    strict_bw_meson_number_density(mass, gamma, T; μ=0.0, degeneracy=1,
                                   qmax=DEFAULT_PHASE_SHIFT_Q_MAX,
                                   q_nodes=DEFAULT_PHASE_SHIFT_Q_NODES,
                                   omega_max=DEFAULT_PHASE_SHIFT_OMEGA_MAX,
                                   omega_nodes=DEFAULT_PHASE_SHIFT_OMEGA_NODES,
                                   gamma_zero_tol=1e-12) -> NamedTuple

reduced strict BW 单通道介子数密度：

```math
n_M^{BW,red}(T)
= d_M \int_0^\infty \frac{dq\,q^2}{2\pi^2}
\int_0^{\omega_{\max}} \frac{d\Delta\omega}{2\pi}
g\!\left(\sqrt{q^2+m_M^2}+\Delta\omega\right)
\frac{\Gamma_M/2}{\Delta\omega^2+\Gamma_M^2/4}.
```

这里采用当前 Stage-1 reduced strict BW 口径：

- `E_M(q)=sqrt(q^2+m_M^2)`
- `Γ_M(q)=Γ_M`
- `ω = E_M(q) + Δω`
"""
function strict_bw_meson_number_density(
    mass::Float64,
    gamma::Float64,
    T::Float64;
    μ::Float64=0.0,
    degeneracy::Integer=1,
    qmax::Float64=DEFAULT_PHASE_SHIFT_Q_MAX,
    q_nodes::Int=DEFAULT_PHASE_SHIFT_Q_NODES,
    omega_max::Float64=DEFAULT_PHASE_SHIFT_OMEGA_MAX,
    omega_nodes::Int=DEFAULT_PHASE_SHIFT_OMEGA_NODES,
    gamma_zero_tol::Float64=1e-12,
)
    _require_nonnegative("mass", mass)
    _require_nonnegative("gamma", gamma)
    _require_nonnegative("temperature T", T)
    degeneracy > 0 || throw(ArgumentError("degeneracy must be positive, got $(degeneracy)"))
    _require_positive_node_count("q_nodes", q_nodes)
    _require_positive_node_count("omega_nodes", omega_nodes)
    qmax > 0.0 || throw(ArgumentError("qmax must be positive, got $(qmax)"))
    omega_max > 0.0 || throw(ArgumentError("omega_max must be positive, got $(omega_max)"))
    gamma_zero_tol >= 0.0 || throw(ArgumentError("gamma_zero_tol must be nonnegative, got $(gamma_zero_tol)"))
    mass > μ || throw(ArgumentError("strict BW helper requires mass > μ, got mass=$(mass), μ=$(μ)"))

    T == 0.0 && return (
        density=0.0,
        q_integral_estimate=0.0,
        omega_shell_at_qmax=0.0,
        qmax=qmax,
        q_nodes=q_nodes,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        degeneracy=Int(degeneracy),
        gamma=gamma,
        mass=mass,
        mode=:strict_bw_reduced,
    )

    if gamma <= gamma_zero_tol
        density = stable_meson_number_density(
            mass,
            T;
            μ=μ,
            degeneracy=degeneracy,
            qmax=qmax,
            num_q_nodes=q_nodes,
        )
        return (
            density=density,
            q_integral_estimate=density / Float64(degeneracy),
            omega_shell_at_qmax=0.0,
            qmax=qmax,
            q_nodes=q_nodes,
            omega_max=omega_max,
            omega_nodes=omega_nodes,
            degeneracy=Int(degeneracy),
            gamma=gamma,
            mass=mass,
            mode=:stable_limit,
        )
    end

    q_grid, q_w = gauleg(0.0, qmax, q_nodes)
    dω_grid, dω_w = gauleg(0.0, omega_max, omega_nodes)

    q_shell_weighted_sum = 0.0
    q_shell_at_qmax = NaN
    @inbounds for iq in eachindex(q_grid, q_w)
        q = q_grid[iq]
        E = hypot(q, mass)
        omega_val = 0.0
        for iω in eachindex(dω_grid, dω_w)
            Δω = dω_grid[iω]
            ω = E + Δω
            omega_val += dω_w[iω] * bose_distribution(ω, μ, T) * _strict_bw_kernel(Δω, gamma)
        end
        omega_val /= (2.0 * π)
        q_shell = (q^2 / (2.0 * π^2)) * omega_val
        q_shell_weighted_sum += q_w[iq] * q_shell
        if iq == length(q_grid)
            q_shell_at_qmax = q_shell
        end
    end

    return (
        density=Float64(degeneracy) * q_shell_weighted_sum,
        q_integral_estimate=q_shell_weighted_sum,
        omega_shell_at_qmax=q_shell_at_qmax,
        qmax=qmax,
        q_nodes=q_nodes,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        degeneracy=Int(degeneracy),
        gamma=gamma,
        mass=mass,
        mode=:strict_bw_reduced,
    )
end

"""
    strict_bw_meson_density_summary(pi_mass, pi_gamma, k_mass, k_gamma, T; kwargs...) -> NamedTuple

聚合 `π/K` reduced strict BW 数密度摘要。
"""
function strict_bw_meson_density_summary(
    pi_mass::Float64,
    pi_gamma::Float64,
    k_mass::Float64,
    k_gamma::Float64,
    T::Float64;
    μ_pi::Float64=0.0,
    μ_K::Float64=0.0,
    d_pi::Integer=meson_degeneracy(:pi),
    d_K::Integer=meson_degeneracy(:K),
    qmax::Float64=DEFAULT_PHASE_SHIFT_Q_MAX,
    q_nodes::Int=DEFAULT_PHASE_SHIFT_Q_NODES,
    omega_max::Float64=DEFAULT_PHASE_SHIFT_OMEGA_MAX,
    omega_nodes::Int=DEFAULT_PHASE_SHIFT_OMEGA_NODES,
    gamma_zero_tol::Float64=1e-12,
)
    pi_density = strict_bw_meson_number_density(
        pi_mass,
        pi_gamma,
        T;
        μ=μ_pi,
        degeneracy=Int(d_pi),
        qmax=qmax,
        q_nodes=q_nodes,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        gamma_zero_tol=gamma_zero_tol,
    )
    k_density = strict_bw_meson_number_density(
        k_mass,
        k_gamma,
        T;
        μ=μ_K,
        degeneracy=Int(d_K),
        qmax=qmax,
        q_nodes=q_nodes,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        gamma_zero_tol=gamma_zero_tol,
    )

    return (
        n_pi=pi_density.density,
        n_K=k_density.density,
        kpi_ratio=iszero(pi_density.density) ? NaN : k_density.density / pi_density.density,
        pi_density=pi_density,
        k_density=k_density,
        qmax=qmax,
        q_nodes=q_nodes,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        gamma_zero_tol=gamma_zero_tol,
    )
end

@inline function _require_phase_shift_isotropic_xi(ξ::Float64)
    abs(ξ) <= 1e-12 && return
    throw(ArgumentError(
        "current phase-shift meson-density helper only supports xi = 0 isotropic reduction, got xi=$(ξ)"
    ))
end

@inline function _require_positive_node_count(name::AbstractString, n::Int)
    n > 1 && return
    throw(ArgumentError("$(name) must be > 1, got $(n)"))
end

@inline function _complex_phase(z::ComplexF64)::Float64
    return atan(imag(z), real(z))
end

function _unwrap_phases(phases::Vector{Float64})
    out = similar(phases)
    isempty(phases) && return out
    out[1] = phases[1]
    shift = 0.0
    for i in 2:length(phases)
        Δ = phases[i] - phases[i - 1]
        if Δ > π
            shift -= 2π
        elseif Δ < -π
            shift += 2π
        end
        out[i] = phases[i] + shift
    end
    return out
end

function _simple_meson_pol_params(meson::Symbol, qp)
    if meson === :pi
        return (
            channel=:P,
            m1=Float64(qp.m.u), m2=Float64(qp.m.u),
            μ1=Float64(qp.μ.u), μ2=Float64(qp.μ.u),
            A1=Float64(qp.A.u), A2=Float64(qp.A.u),
            num_s_quark=0,
        )
    elseif meson === :K
        return (
            channel=:P,
            m1=Float64(qp.m.u), m2=Float64(qp.m.s),
            μ1=Float64(qp.μ.u), μ2=Float64(qp.μ.s),
            A1=Float64(qp.A.u), A2=Float64(qp.A.s),
            num_s_quark=1,
        )
    end
    throw(ArgumentError("Unsupported simple meson: $(meson)"))
end

function _build_k_coeffs(qp)
    G_u = calculate_G_from_A(Float64(qp.A.u), Float64(qp.m.u))
    G_s = calculate_G_from_A(Float64(qp.A.s), Float64(qp.m.s))
    return calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)
end

function _propagator_phase(meson::Symbol, ω::Float64, q::Float64, qp, tp, K_coeffs; eta::Float64)
    pol = _simple_meson_pol_params(meson, qp)
    Π_re, Π_im = polarization_with_width(
        pol.channel, ω, 2.0 * eta, q,
        pol.m1, pol.m2,
        pol.μ1, pol.μ2,
        Float64(tp.T), Float64(tp.Φ), Float64(tp.Φbar), Float64(tp.ξ),
        pol.A1, pol.A2, pol.num_s_quark,
    )
    Π = ComplexF64(Π_re, Π_im)
    D = meson_propagator_simple(meson, K_coeffs, Π)
    return _complex_phase(D)
end

"""
    phase_shift_meson_number_density(meson, quark_params, thermo_params; kwargs...) -> NamedTuple

当前最小 Phase-E3 口径下的单通道相移介子数密度辅助函数。

约束：
- 仅支持 `xi = 0`
- 当前聚合通道只支持 `:pi` 与 `:K`
- 数值积分固定采用 GL + 硬截断
"""
function phase_shift_meson_number_density(
    meson::Symbol,
    quark_params,
    thermo_params;
    degeneracy::Integer=meson_degeneracy(meson),
    qmax::Float64=DEFAULT_PHASE_SHIFT_Q_MAX,
    q_nodes::Int=DEFAULT_PHASE_SHIFT_Q_NODES,
    omega_min::Float64=0.05,
    omega_max::Float64=DEFAULT_PHASE_SHIFT_OMEGA_MAX,
    omega_nodes::Int=DEFAULT_PHASE_SHIFT_OMEGA_NODES,
    eta::Float64=DEFAULT_PHASE_SHIFT_ETA,
)
    degeneracy > 0 || throw(ArgumentError("degeneracy must be positive, got $(degeneracy)"))
    _require_positive_node_count("q_nodes", q_nodes)
    _require_positive_node_count("omega_nodes", omega_nodes)
    qmax > 0.0 || throw(ArgumentError("qmax must be positive, got $(qmax)"))
    omega_max > omega_min || throw(ArgumentError("omega_max must exceed omega_min"))

    tp = thermo_params
    _require_nonnegative("temperature T", Float64(tp.T))
    _require_phase_shift_isotropic_xi(Float64(tp.ξ))
    Float64(tp.T) > 0.0 || return (
        meson=meson,
        density=0.0,
        q_integral_estimate=0.0,
        omega_shell_at_qmax=0.0,
        qmax=qmax,
        q_nodes=q_nodes,
        omega_min=omega_min,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        degeneracy=Int(degeneracy),
    )

    qp = ensure_quark_params_has_A(quark_params, tp)
    K_coeffs = _build_k_coeffs(qp)
    q_grid, q_w = gauleg(0.0, qmax, q_nodes)
    omega_grid, omega_w = gauleg(omega_min, omega_max, omega_nodes)

    q_shell_weighted_sum = 0.0
    q_shell_at_qmax = NaN
    @inbounds for iq in eachindex(q_grid, q_w)
        q = q_grid[iq]
        phases = [_propagator_phase(meson, ω, q, qp, tp, K_coeffs; eta=eta) for ω in omega_grid]
        phase_unwrapped = _unwrap_phases(phases)
        omega_val = 0.0
        for iω in eachindex(omega_grid, omega_w, phase_unwrapped)
            gω = bose_distribution(Float64(omega_grid[iω]), 0.0, Float64(tp.T))
            omega_val += omega_w[iω] * gω * (1.0 + gω) * phase_unwrapped[iω]
        end
        omega_val /= (2.0 * π)
        q_shell = (q^2 / (2.0 * π^2)) * omega_val
        q_shell_weighted_sum += q_w[iq] * q_shell
        if iq == length(q_grid)
            q_shell_at_qmax = q_shell
        end
    end

    density = (Float64(degeneracy) / Float64(tp.T)) * q_shell_weighted_sum
    return (
        meson=meson,
        density=density,
        q_integral_estimate=q_shell_weighted_sum,
        omega_shell_at_qmax=q_shell_at_qmax,
        qmax=qmax,
        q_nodes=q_nodes,
        omega_min=omega_min,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        degeneracy=Int(degeneracy),
    )
end

"""
    phase_shift_meson_density_summary(quark_params, thermo_params; kwargs...) -> NamedTuple

当前最小 Phase-E3 口径下的 `π/K` 聚合通道相移介子数密度总结。
"""
function phase_shift_meson_density_summary(
    quark_params,
    thermo_params;
    pi_channel::Symbol=:pi,
    k_channel::Symbol=:K,
    d_pi::Integer=meson_degeneracy(:pi),
    d_K::Integer=meson_degeneracy(:K),
    qmax::Float64=DEFAULT_PHASE_SHIFT_Q_MAX,
    q_nodes::Int=DEFAULT_PHASE_SHIFT_Q_NODES,
    omega_min::Float64=0.05,
    omega_max::Float64=DEFAULT_PHASE_SHIFT_OMEGA_MAX,
    omega_nodes::Int=DEFAULT_PHASE_SHIFT_OMEGA_NODES,
    eta::Float64=DEFAULT_PHASE_SHIFT_ETA,
)
    pi_density = phase_shift_meson_number_density(
        pi_channel, quark_params, thermo_params;
        degeneracy=Int(d_pi),
        qmax=qmax,
        q_nodes=q_nodes,
        omega_min=omega_min,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        eta=eta,
    )
    k_density = phase_shift_meson_number_density(
        k_channel, quark_params, thermo_params;
        degeneracy=Int(d_K),
        qmax=qmax,
        q_nodes=q_nodes,
        omega_min=omega_min,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        eta=eta,
    )

    n_pi = Float64(pi_density.density)
    n_K = Float64(k_density.density)
    return (
        T_fm=Float64(thermo_params.T),
        xi=Float64(thermo_params.ξ),
        pi_channel=pi_channel,
        k_channel=k_channel,
        d_pi=Int(d_pi),
        d_K=Int(d_K),
        n_pi=n_pi,
        n_K=n_K,
        kpi_ratio=iszero(n_pi) ? NaN : n_K / n_pi,
        qmax=qmax,
        q_nodes=q_nodes,
        omega_min=omega_min,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        eta=eta,
        pi_density=pi_density,
        k_density=k_density,
    )
end

end # module MesonDensity
