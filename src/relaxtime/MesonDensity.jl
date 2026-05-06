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

using ForwardDiff
using ..GaussLegendre: gauleg
using ..AFieldBuilder: ensure_quark_params_has_A
using ..EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using ..MesonMass: solve_meson_mass
using ..PolarizationAniso: polarization_with_width
using ..MesonPropagator: meson_propagator_simple
using Main.Constants_PNJL: G_fm2, K_fm5

export DEFAULT_MESON_DENSITY_Q_NODES
export DEFAULT_PHASE_SHIFT_Q_MAX, DEFAULT_PHASE_SHIFT_Q_NODES
export DEFAULT_PHASE_SHIFT_OMEGA_MAX, DEFAULT_PHASE_SHIFT_OMEGA_NODES
export bose_distribution, meson_degeneracy
export stable_meson_number_density, stable_kpi_ratio, stable_kpi_scan
export strict_bw_meson_number_density, strict_bw_meson_density_summary
export strict_bw_qpole_meson_number_density, strict_bw_qpole_density_summary
export phase_shift_meson_number_density, phase_shift_meson_density_summary
export phase_shift_meson_number_density_derivative_reference, phase_shift_meson_density_derivative_reference_summary
export phase_shift_point_diagnostic

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
    elseif meson === :pi_plus || meson === :pi_minus || meson === :K_plus || meson === :K_minus
        return 1
    else
        throw(ArgumentError("Unsupported meson $(meson). Use :pi, :pi_plus, :pi_minus, :K, :K_plus, or :K_minus."))
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

@inline function _strict_bw_energy_kernel(ω::Float64, E::Float64, gamma::Float64)::Float64
    half_gamma = gamma / 2.0
    return half_gamma / ((ω - E)^2 + half_gamma^2)
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

@inline function _sorted_gauleg(a::Float64, b::Float64, n::Int)
    nodes, weights = gauleg(a, b, n)
    perm = sortperm(nodes)
    return nodes[perm], weights[perm]
end

@inline function _strict_bw_stage_symbol(stage::Symbol)::Symbol
    if stage === :stage1_reduced || stage === :reduced
        return :strict_bw_stage1_reduced
    elseif stage === :stage2_qpole || stage === :qpole
        return :strict_bw_stage2_qpole
    end
    throw(ArgumentError("strict BW stage must be :stage1_reduced/:reduced or :stage2_qpole/:qpole, got $(stage)"))
end

function _solve_qpole_at_q(
    meson::Symbol,
    q::Float64,
    quark_params,
    thermo_params;
    initial_mass::Float64,
    initial_gamma::Float64,
    solver_iterations::Int,
    residual_norm_max::Float64,
    require_converged::Bool,
)
    seeds = (
        (initial_mass, max(initial_gamma, 0.0)),
        (max(initial_mass, 1e-6), 0.0),
        (hypot(q, max(initial_mass, 1e-6)), max(initial_gamma, 0.0)),
    )

    best = nothing
    for (m0, g0) in seeds
        for method in (:trust_region, :newton)
            rr = try
                solve_meson_mass(
                    meson,
                    quark_params,
                    thermo_params;
                    k_norm=q,
                    initial_mass=Float64(m0),
                    initial_gamma=Float64(max(g0, 0.0)),
                    method=method,
                    iterations=solver_iterations,
                )
            catch
                nothing
            end
            rr === nothing && continue
            if best === nothing || Float64(rr.residual_norm) < Float64(best.residual_norm)
                best = rr
            end
            if isfinite(Float64(rr.residual_norm)) &&
               Float64(rr.residual_norm) <= residual_norm_max &&
               (!require_converged || Bool(rr.converged))
                return (
                    mass=Float64(rr.mass),
                    gamma=max(Float64(rr.gamma), 0.0),
                    residual_norm=Float64(rr.residual_norm),
                    converged=Bool(rr.converged),
                    accepted=true,
                    method=method,
                )
            end
        end
    end

    best === nothing && return (
        mass=NaN,
        gamma=NaN,
        residual_norm=Inf,
        converged=false,
        accepted=false,
        method=:failed,
    )

    return (
        mass=Float64(best.mass),
        gamma=max(Float64(best.gamma), 0.0),
        residual_norm=Float64(best.residual_norm),
        converged=Bool(best.converged),
        accepted=isfinite(Float64(best.residual_norm)) &&
                 Float64(best.residual_norm) <= residual_norm_max &&
                 (!require_converged || Bool(best.converged)),
        method=:best_effort,
    )
end

"""
    strict_bw_qpole_meson_number_density(meson, mass0, gamma0, quark_params, thermo_params; ...) -> NamedTuple

Stage-2 q 依赖复极点版本的 strict BW 单通道介子数密度。
"""
function strict_bw_qpole_meson_number_density(
    meson::Symbol,
    mass0::Float64,
    gamma0::Float64,
    quark_params,
    thermo_params;
    μ::Float64=0.0,
    degeneracy::Integer=1,
    qmax::Float64=DEFAULT_PHASE_SHIFT_Q_MAX,
    q_nodes::Int=DEFAULT_PHASE_SHIFT_Q_NODES,
    omega_max::Float64=DEFAULT_PHASE_SHIFT_OMEGA_MAX,
    omega_nodes::Int=DEFAULT_PHASE_SHIFT_OMEGA_NODES,
    gamma_zero_tol::Float64=1e-12,
    solver_iterations::Int=20,
    pole_residual_norm_max::Float64=1e-6,
    pole_require_converged::Bool=true,
)
    degeneracy > 0 || throw(ArgumentError("degeneracy must be positive, got $(degeneracy)"))
    _require_nonnegative("mass0", mass0)
    _require_nonnegative("gamma0", gamma0)
    _require_nonnegative("temperature T", Float64(thermo_params.T))
    qmax > 0.0 || throw(ArgumentError("qmax must be positive, got $(qmax)"))
    _require_positive_node_count("q_nodes", q_nodes)
    _require_positive_node_count("omega_nodes", omega_nodes)
    omega_max > 0.0 || throw(ArgumentError("omega_max must be positive, got $(omega_max)"))

    T = Float64(thermo_params.T)
    T == 0.0 && return (
        meson=meson,
        density=0.0,
        q_integral_estimate=0.0,
        omega_shell_at_qmax=0.0,
        q_values=Float64[],
        E_values=Float64[],
        gamma_values=Float64[],
        residual_norms=Float64[],
        converged_flags=Bool[],
        accepted_flags=Bool[],
        qmax=qmax,
        q_nodes=q_nodes,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        degeneracy=Int(degeneracy),
        mode=:strict_bw_stage2_qpole,
    )

    if gamma0 <= gamma_zero_tol
        density = stable_meson_number_density(
            mass0,
            T;
            μ=μ,
            degeneracy=degeneracy,
            qmax=qmax,
            num_q_nodes=q_nodes,
        )
        return (
            meson=meson,
            density=density,
            q_integral_estimate=density / Float64(degeneracy),
            omega_shell_at_qmax=0.0,
            q_values=Float64[],
            E_values=Float64[],
            gamma_values=Float64[],
            residual_norms=Float64[],
            converged_flags=Bool[],
            accepted_flags=Bool[],
            qmax=qmax,
            q_nodes=q_nodes,
            omega_max=omega_max,
            omega_nodes=omega_nodes,
            degeneracy=Int(degeneracy),
            mode=:stable_limit,
        )
    end

    qp = ensure_quark_params_has_A(quark_params, thermo_params)
    q_grid, q_w = _sorted_gauleg(0.0, qmax, q_nodes)
    ω_grid, ω_w = gauleg(0.0, omega_max, omega_nodes)

    q_values = Float64[]
    E_values = Float64[]
    gamma_values = Float64[]
    residual_norms = Float64[]
    converged_flags = Bool[]
    accepted_flags = Bool[]

    q_shell_weighted_sum = 0.0
    q_shell_at_qmax = NaN
    prev_mass = mass0
    prev_gamma = gamma0

    @inbounds for iq in eachindex(q_grid, q_w)
        q = Float64(q_grid[iq])
        pole = _solve_qpole_at_q(
            meson,
            q,
            qp,
            thermo_params;
            initial_mass=prev_mass,
            initial_gamma=prev_gamma,
            solver_iterations=solver_iterations,
            residual_norm_max=pole_residual_norm_max,
            require_converged=pole_require_converged,
        )

        E_q = pole.mass
        gamma_q = pole.gamma
        push!(q_values, q)
        push!(E_values, E_q)
        push!(gamma_values, gamma_q)
        push!(residual_norms, pole.residual_norm)
        push!(converged_flags, pole.converged)
        push!(accepted_flags, pole.accepted)

        if pole.accepted && isfinite(E_q) && isfinite(gamma_q) && E_q > μ
            omega_val = 0.0
            for iω in eachindex(ω_grid, ω_w)
                ω = ω_grid[iω]
                omega_val += ω_w[iω] * bose_distribution(ω, μ, T) * _strict_bw_energy_kernel(ω, E_q, gamma_q)
            end
            omega_val /= (2.0 * π)
            q_shell = (q^2 / (2.0 * π^2)) * omega_val
            q_shell_weighted_sum += q_w[iq] * q_shell
            if iq == length(q_grid)
                q_shell_at_qmax = q_shell
            end
            prev_mass = E_q
            prev_gamma = gamma_q
        else
            if iq == length(q_grid)
                q_shell_at_qmax = NaN
            end
        end
    end

    return (
        meson=meson,
        density=Float64(degeneracy) * q_shell_weighted_sum,
        q_integral_estimate=q_shell_weighted_sum,
        omega_shell_at_qmax=q_shell_at_qmax,
        q_values=q_values,
        E_values=E_values,
        gamma_values=gamma_values,
        residual_norms=residual_norms,
        converged_flags=converged_flags,
        accepted_flags=accepted_flags,
        qmax=qmax,
        q_nodes=q_nodes,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        degeneracy=Int(degeneracy),
        mode=:strict_bw_stage2_qpole,
    )
end

function strict_bw_qpole_density_summary(
    pi_mass::Float64,
    pi_gamma::Float64,
    k_mass::Float64,
    k_gamma::Float64,
    quark_params,
    thermo_params;
    pi_channel::Symbol=:pi,
    k_channel::Symbol=:K,
    μ_pi::Float64=0.0,
    μ_K::Float64=0.0,
    d_pi::Integer=meson_degeneracy(:pi),
    d_K::Integer=meson_degeneracy(:K),
    qmax::Float64=DEFAULT_PHASE_SHIFT_Q_MAX,
    q_nodes::Int=DEFAULT_PHASE_SHIFT_Q_NODES,
    omega_max::Float64=DEFAULT_PHASE_SHIFT_OMEGA_MAX,
    omega_nodes::Int=DEFAULT_PHASE_SHIFT_OMEGA_NODES,
    gamma_zero_tol::Float64=1e-12,
    solver_iterations::Int=20,
    pole_residual_norm_max::Float64=1e-6,
    pole_require_converged::Bool=true,
)
    pi_density = strict_bw_qpole_meson_number_density(
        pi_channel,
        pi_mass,
        pi_gamma,
        quark_params,
        thermo_params;
        μ=μ_pi,
        degeneracy=Int(d_pi),
        qmax=qmax,
        q_nodes=q_nodes,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        gamma_zero_tol=gamma_zero_tol,
        solver_iterations=solver_iterations,
        pole_residual_norm_max=pole_residual_norm_max,
        pole_require_converged=pole_require_converged,
    )
    k_density = strict_bw_qpole_meson_number_density(
        k_channel,
        k_mass,
        k_gamma,
        quark_params,
        thermo_params;
        μ=μ_K,
        degeneracy=Int(d_K),
        qmax=qmax,
        q_nodes=q_nodes,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        gamma_zero_tol=gamma_zero_tol,
        solver_iterations=solver_iterations,
        pole_residual_norm_max=pole_residual_norm_max,
        pole_require_converged=pole_require_converged,
    )

    return (
        n_pi=pi_density.density,
        n_K=k_density.density,
        kpi_ratio=iszero(pi_density.density) ? NaN : k_density.density / pi_density.density,
        pi_channel=pi_channel,
        k_channel=k_channel,
        pi_density=pi_density,
        k_density=k_density,
        qmax=qmax,
        q_nodes=q_nodes,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        gamma_zero_tol=gamma_zero_tol,
        solver_iterations=solver_iterations,
        pole_residual_norm_max=pole_residual_norm_max,
        pole_require_converged=pole_require_converged,
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

@inline function _phase_shift_omega_lower_bound(omega_min::Float64, μ::Float64)::Float64
    return omega_min > μ ? omega_min : nextfloat(μ)
end

@inline function _complex_phase(z::Complex{T})::T where {T<:Real}
    return atan(imag(z), real(z))
end

@inline function _phase_derivative_from_components(reD::T, imD::T, dre::T, dim::T) where {T<:Real}
    denom = reD * reD + imD * imD
    return (reD * dim - imD * dre) / denom
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

@inline function _phase_shift_scheme_symbol(scheme::Symbol)::Symbol
    if scheme === :current || scheme === :phase_e3 || scheme === :phase_shift_current
        return :phase_shift_current
    elseif scheme === :gbu || scheme === :gbu_reference || scheme === :generalized_bu || scheme === :phase_shift_gbu_reference
        return :phase_shift_gbu_reference
    end
    throw(ArgumentError("phase-shift scheme must be :current/:phase_e3 or :gbu/:gbu_reference/:generalized_bu, got $(scheme)"))
end

@inline function _gbu_phase_function(δ::T)::T where {T<:Real}
    return δ - 0.5 * sin(2.0 * δ)
end

@inline function _phase_shift_weighted_phase(δ::T, scheme::Symbol)::T where {T<:Real}
    scheme_sym = _phase_shift_scheme_symbol(scheme)
    if scheme_sym === :phase_shift_current
        return δ
    end
    return _gbu_phase_function(δ)
end

@inline function _phase_shift_weight_derivative_factor(δ::T, scheme::Symbol)::T where {T<:Real}
    scheme_sym = _phase_shift_scheme_symbol(scheme)
    if scheme_sym === :phase_shift_current
        return one(T)
    end
    return T(2) * sin(δ)^2
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
    elseif meson === :pi_plus
        return (
            channel=:P,
            m1=Float64(qp.m.u), m2=Float64(qp.m.d),
            μ1=Float64(qp.μ.u), μ2=Float64(qp.μ.d),
            A1=Float64(qp.A.u), A2=Float64(qp.A.d),
            num_s_quark=0,
        )
    elseif meson === :pi_minus
        return (
            channel=:P,
            m1=Float64(qp.m.d), m2=Float64(qp.m.u),
            μ1=Float64(qp.μ.d), μ2=Float64(qp.μ.u),
            A1=Float64(qp.A.d), A2=Float64(qp.A.u),
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
    elseif meson === :K_plus
        return (
            channel=:P,
            m1=Float64(qp.m.u), m2=Float64(qp.m.s),
            μ1=Float64(qp.μ.u), μ2=Float64(qp.μ.s),
            A1=Float64(qp.A.u), A2=Float64(qp.A.s),
            num_s_quark=1,
        )
    elseif meson === :K_minus
        return (
            channel=:P,
            m1=Float64(qp.m.s), m2=Float64(qp.m.u),
            μ1=Float64(qp.μ.s), μ2=Float64(qp.μ.u),
            A1=Float64(qp.A.s), A2=Float64(qp.A.u),
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

function _propagator_components(meson::Symbol, ω::T, q::Float64, qp, tp, K_coeffs; eta::Float64) where {T<:Real}
    pol = _simple_meson_pol_params(meson, qp)
    Π_re, Π_im = polarization_with_width(
        pol.channel, ω, 2.0 * eta, q,
        pol.m1, pol.m2,
        pol.μ1, pol.μ2,
        Float64(tp.T), Float64(tp.Φ), Float64(tp.Φbar), Float64(tp.ξ),
        pol.A1, pol.A2, pol.num_s_quark,
    )
    Π = Complex{T}(Π_re, Π_im)
    D = meson_propagator_simple(meson, K_coeffs, Π)
    return real(D), imag(D)
end

function _propagator_phase(meson::Symbol, ω::T, q::Float64, qp, tp, K_coeffs; eta::Float64) where {T<:Real}
    reD, imD = _propagator_components(meson, ω, q, qp, tp, K_coeffs; eta=eta)
    return atan(imD, reD)
end

function phase_shift_point_diagnostic(
    meson::Symbol,
    ω::Float64,
    q::Float64,
    quark_params,
    thermo_params;
    scheme::Symbol=:current,
    eta::Float64=DEFAULT_PHASE_SHIFT_ETA,
    fd_step::Float64=1e-5,
)
    fd_step > 0.0 || throw(ArgumentError("fd_step must be positive, got $(fd_step)"))
    qp = ensure_quark_params_has_A(quark_params, thermo_params)
    tp = thermo_params
    K_coeffs = _build_k_coeffs(qp)
    scheme_sym = _phase_shift_scheme_symbol(scheme)

    re0, im0 = _propagator_components(meson, ω, q, qp, tp, K_coeffs; eta=eta)
    phase0 = atan(im0, re0)
    weighted0 = _phase_shift_weighted_phase(phase0, scheme_sym)

    re_fun(x) = first(_propagator_components(meson, x, q, qp, tp, K_coeffs; eta=eta))
    im_fun(x) = last(_propagator_components(meson, x, q, qp, tp, K_coeffs; eta=eta))
    phase_fun(x) = _propagator_phase(meson, x, q, qp, tp, K_coeffs; eta=eta)
    weighted_fun(x) = _phase_shift_weighted_phase(phase_fun(x), scheme_sym)

    dre = ForwardDiff.derivative(re_fun, ω)
    dim = ForwardDiff.derivative(im_fun, ω)
    dphase_ad = ForwardDiff.derivative(phase_fun, ω)
    dphase_formula = _phase_derivative_from_components(re0, im0, dre, dim)
    dweighted_ad = ForwardDiff.derivative(weighted_fun, ω)

    dre_fd = (re_fun(ω + fd_step) - re_fun(ω - fd_step)) / (2fd_step)
    dim_fd = (im_fun(ω + fd_step) - im_fun(ω - fd_step)) / (2fd_step)
    dphase_fd = _phase_derivative_from_components(re0, im0, dre_fd, dim_fd)
    dphase_raw_fd = (phase_fun(ω + fd_step) - phase_fun(ω - fd_step)) / (2fd_step)
    dweighted_fd = (weighted_fun(ω + fd_step) - weighted_fun(ω - fd_step)) / (2fd_step)

    return (
        meson=meson,
        omega=ω,
        q=q,
        xi=Float64(tp.ξ),
        temperature=Float64(tp.T),
        eta=eta,
        fd_step=fd_step,
        scheme=scheme_sym,
        reD=Float64(re0),
        imD=Float64(im0),
        phase=Float64(phase0),
        weighted_phase=Float64(weighted0),
        dReD_domega=Float64(dre),
        dImD_domega=Float64(dim),
        dReD_fd=Float64(dre_fd),
        dImD_fd=Float64(dim_fd),
        dphase_ad=Float64(dphase_ad),
        dphase_formula=Float64(dphase_formula),
        dphase_fd=Float64(dphase_fd),
        dphase_raw_fd=Float64(dphase_raw_fd),
        dphase_abs_diff=abs(Float64(dphase_ad) - Float64(dphase_fd)),
        dphase_formula_abs_diff=abs(Float64(dphase_formula) - Float64(dphase_ad)),
        dphase_raw_fd_abs_diff=abs(Float64(dphase_ad) - Float64(dphase_raw_fd)),
        dweighted_ad=Float64(dweighted_ad),
        dweighted_fd=Float64(dweighted_fd),
        dweighted_abs_diff=abs(Float64(dweighted_ad) - Float64(dweighted_fd)),
    )
end

function phase_shift_meson_number_density_derivative_reference(
    meson::Symbol,
    quark_params,
    thermo_params;
    degeneracy::Integer=meson_degeneracy(meson),
    μ::Float64=0.0,
    scheme::Symbol=:current,
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
    omega_lower = _phase_shift_omega_lower_bound(omega_min, μ)
    omega_max > omega_lower || throw(ArgumentError("omega_max must exceed effective omega_min=$(omega_lower)"))
    scheme_sym = _phase_shift_scheme_symbol(scheme)

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
        scheme=scheme_sym,
        derivative_backend=:forwarddiff,
        max_formula_abs_diff=0.0,
    )

    qp = ensure_quark_params_has_A(quark_params, tp)
    K_coeffs = _build_k_coeffs(qp)
    q_grid, q_w = gauleg(0.0, qmax, q_nodes)
    omega_grid, omega_w = gauleg(omega_lower, omega_max, omega_nodes)

    function q_shell_from_q(q::Float64)
        omega_val = 0.0
        for iω in eachindex(omega_grid, omega_w)
            ω = Float64(omega_grid[iω])
            gω = bose_distribution(ω, μ, Float64(tp.T))
            re0, im0 = _propagator_components(meson, ω, q, qp, tp, K_coeffs; eta=eta)
            phase0 = atan(im0, re0)
            weight_factor = _phase_shift_weight_derivative_factor(phase0, scheme_sym)

            re_fun(x) = first(_propagator_components(meson, x, q, qp, tp, K_coeffs; eta=eta))
            im_fun(x) = last(_propagator_components(meson, x, q, qp, tp, K_coeffs; eta=eta))
            weighted_fun(x) = _phase_shift_weighted_phase(_propagator_phase(meson, x, q, qp, tp, K_coeffs; eta=eta), scheme_sym)

            dre = ForwardDiff.derivative(re_fun, ω)
            dim = ForwardDiff.derivative(im_fun, ω)
            dphase_formula = _phase_derivative_from_components(re0, im0, dre, dim)
            dweighted_formula = weight_factor * dphase_formula
            dweighted_ad = ForwardDiff.derivative(weighted_fun, ω)
            max_formula_abs_diff = max(max_formula_abs_diff, abs(Float64(dweighted_ad) - Float64(dweighted_formula)))

            omega_val += omega_w[iω] * gω * dweighted_ad
        end
        omega_val /= (2.0 * π)
        return (q^2 / (2.0 * π^2)) * omega_val
    end

    q_shell_weighted_sum = 0.0
    q_shell_at_qmax = 0.0
    max_formula_abs_diff = 0.0
    @inbounds for iq in eachindex(q_grid, q_w)
        q = q_grid[iq]
        q_shell = q_shell_from_q(q)
        q_shell_weighted_sum += q_w[iq] * q_shell
    end
    q_shell_at_qmax = q_shell_from_q(qmax)

    density = Float64(degeneracy) * q_shell_weighted_sum
    return (
        meson=meson,
        density=density,
        q_integral_estimate=q_shell_weighted_sum,
        omega_shell_at_qmax=q_shell_at_qmax,
        qmax=qmax,
        q_nodes=q_nodes,
        omega_min=omega_lower,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        degeneracy=Int(degeneracy),
        scheme=scheme_sym,
        derivative_backend=:forwarddiff,
        max_formula_abs_diff=max_formula_abs_diff,
    )
end

function phase_shift_meson_density_derivative_reference_summary(
    quark_params,
    thermo_params;
    pi_channel::Symbol=:pi,
    k_channel::Symbol=:K,
    μ_pi::Float64=0.0,
    μ_K::Float64=0.0,
    d_pi::Integer=meson_degeneracy(:pi),
    d_K::Integer=meson_degeneracy(:K),
    scheme::Symbol=:current,
    qmax::Float64=DEFAULT_PHASE_SHIFT_Q_MAX,
    q_nodes::Int=DEFAULT_PHASE_SHIFT_Q_NODES,
    omega_min::Float64=0.05,
    omega_max::Float64=DEFAULT_PHASE_SHIFT_OMEGA_MAX,
    omega_nodes::Int=DEFAULT_PHASE_SHIFT_OMEGA_NODES,
    eta::Float64=DEFAULT_PHASE_SHIFT_ETA,
)
    pi_density = phase_shift_meson_number_density_derivative_reference(
        pi_channel, quark_params, thermo_params;
        μ=μ_pi,
        degeneracy=Int(d_pi), scheme=scheme, qmax=qmax, q_nodes=q_nodes,
        omega_min=omega_min, omega_max=omega_max, omega_nodes=omega_nodes, eta=eta,
    )
    k_density = phase_shift_meson_number_density_derivative_reference(
        k_channel, quark_params, thermo_params;
        μ=μ_K,
        degeneracy=Int(d_K), scheme=scheme, qmax=qmax, q_nodes=q_nodes,
        omega_min=omega_min, omega_max=omega_max, omega_nodes=omega_nodes, eta=eta,
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
        scheme=_phase_shift_scheme_symbol(scheme),
        derivative_backend=:forwarddiff,
        max_formula_abs_diff=max(pi_density.max_formula_abs_diff, k_density.max_formula_abs_diff),
        pi_density=pi_density,
        k_density=k_density,
    )
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
    μ::Float64=0.0,
    scheme::Symbol=:current,
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
    omega_lower = _phase_shift_omega_lower_bound(omega_min, μ)
    omega_max > omega_lower || throw(ArgumentError("omega_max must exceed effective omega_min=$(omega_lower)"))
    scheme_sym = _phase_shift_scheme_symbol(scheme)

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
        scheme=scheme_sym,
    )

    qp = ensure_quark_params_has_A(quark_params, tp)
    K_coeffs = _build_k_coeffs(qp)
    q_grid, q_w = gauleg(0.0, qmax, q_nodes)
    omega_grid, omega_w = gauleg(omega_lower, omega_max, omega_nodes)

    q_shell_weighted_sum = 0.0
    q_shell_at_qmax = NaN
    @inbounds for iq in eachindex(q_grid, q_w)
        q = q_grid[iq]
        phases = [_propagator_phase(meson, ω, q, qp, tp, K_coeffs; eta=eta) for ω in omega_grid]
        phase_unwrapped = _unwrap_phases(phases)
        omega_val = 0.0
        for iω in eachindex(omega_grid, omega_w, phase_unwrapped)
            gω = bose_distribution(Float64(omega_grid[iω]), μ, Float64(tp.T))
            omega_val += omega_w[iω] * gω * (1.0 + gω) * _phase_shift_weighted_phase(phase_unwrapped[iω], scheme_sym)
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
        omega_min=omega_lower,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        degeneracy=Int(degeneracy),
        scheme=scheme_sym,
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
    μ_pi::Float64=0.0,
    μ_K::Float64=0.0,
    d_pi::Integer=meson_degeneracy(:pi),
    d_K::Integer=meson_degeneracy(:K),
    scheme::Symbol=:current,
    qmax::Float64=DEFAULT_PHASE_SHIFT_Q_MAX,
    q_nodes::Int=DEFAULT_PHASE_SHIFT_Q_NODES,
    omega_min::Float64=0.05,
    omega_max::Float64=DEFAULT_PHASE_SHIFT_OMEGA_MAX,
    omega_nodes::Int=DEFAULT_PHASE_SHIFT_OMEGA_NODES,
    eta::Float64=DEFAULT_PHASE_SHIFT_ETA,
)
    pi_density = phase_shift_meson_number_density(
        pi_channel, quark_params, thermo_params;
        μ=μ_pi,
        degeneracy=Int(d_pi),
        scheme=scheme,
        qmax=qmax,
        q_nodes=q_nodes,
        omega_min=omega_min,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        eta=eta,
    )
    k_density = phase_shift_meson_number_density(
        k_channel, quark_params, thermo_params;
        μ=μ_K,
        degeneracy=Int(d_K),
        scheme=scheme,
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
        scheme=_phase_shift_scheme_symbol(scheme),
        pi_density=pi_density,
        k_density=k_density,
    )
end

end # module MesonDensity
