raw"""
MesonThermodynamics

介子热力学最小实现入口。当前阶段聚焦 `π/K` 通道的 pressure 侧能力，
并与现有 meson density workflow 保持同口径并行演进：

- stable particle limit
- reduced strict BW
- current phase-shift
- generalized BU reference
"""
module MesonThermodynamics

using ..GaussLegendre: gauleg
using ..AFieldBuilder: ensure_quark_params_has_A
using ..MesonDensity: DEFAULT_PHASE_SHIFT_ETA,
                      DEFAULT_PHASE_SHIFT_OMEGA_MAX,
                      DEFAULT_PHASE_SHIFT_OMEGA_NODES,
                      DEFAULT_PHASE_SHIFT_Q_MAX,
                      DEFAULT_PHASE_SHIFT_Q_NODES,
                      _build_k_coeffs,
                      _phase_shift_omega_lower_bound,
                      _phase_shift_scheme_symbol,
                      _phase_shift_weighted_phase,
                      _propagator_phase,
                      _require_phase_shift_isotropic_xi,
                      _require_positive_node_count,
                      _strict_bw_kernel,
                      _unwrap_phases,
                      bose_distribution,
                      meson_degeneracy

export bosonic_log_pressure_factor
export stable_meson_pressure, stable_meson_pressure_summary
export strict_bw_meson_pressure, strict_bw_meson_pressure_summary
export phase_shift_meson_pressure, phase_shift_meson_pressure_summary

@inline function _require_nonnegative(name::AbstractString, value::Float64)
    value >= 0.0 && return
    throw(ArgumentError("$(name) must be nonnegative, got $(value)"))
end

@inline function _default_qmax(mass::Float64, T::Float64, μ::Float64)::Float64
    gap = max(mass - μ, 0.0)
    return max(8.0, 20.0 * T + 10.0 * gap, 8.0 * mass + 10.0 * T)
end

@inline function _normalize_thermo_namedtuple(tp)
    return (
        T=Float64(tp.T),
        Φ=Float64(tp.Φ),
        Φbar=Float64(tp.Φbar),
        ξ=Float64(tp.ξ),
    )
end

@inline function _normalize_quark_namedtuple(qp)
    return (
        m=(u=Float64(qp.m.u), d=Float64(qp.m.d), s=Float64(qp.m.s)),
        μ=(u=Float64(qp.μ.u), d=Float64(qp.μ.d), s=Float64(qp.μ.s)),
    )
end

function bosonic_log_pressure_factor(E::Float64, μ::Float64, T::Float64)::Float64
    _require_nonnegative("temperature T", T)
    T == 0.0 && return 0.0
    E > μ || throw(ArgumentError("Bosonic pressure kernel requires E > μ, got E=$(E), μ=$(μ)"))

    exponent = (E - μ) / T
    exponent > 700.0 && return 0.0
    return -T * log1p(-exp(-exponent))
end

function stable_meson_pressure(
    mass::Float64,
    T::Float64;
    μ::Float64=0.0,
    degeneracy::Integer=1,
    qmax::Union{Nothing,Float64}=nothing,
    num_q_nodes::Int=256,
)::Float64
    _require_nonnegative("mass", mass)
    _require_nonnegative("temperature T", T)
    degeneracy > 0 || throw(ArgumentError("degeneracy must be positive, got $(degeneracy)"))
    num_q_nodes > 1 || throw(ArgumentError("num_q_nodes must be > 1, got $(num_q_nodes)"))
    mass > μ || throw(ArgumentError("Stable meson pressure requires mass > μ, got mass=$(mass), μ=$(μ)"))
    T == 0.0 && return 0.0

    q_upper = qmax === nothing ? _default_qmax(mass, T, μ) : Float64(qmax)
    q_upper > 0.0 || throw(ArgumentError("qmax must be positive, got $(q_upper)"))

    nodes, weights = gauleg(0.0, q_upper, num_q_nodes)
    integral = 0.0
    @inbounds for i in eachindex(nodes, weights)
        q = nodes[i]
        E = hypot(q, mass)
        integral += weights[i] * q^2 * bosonic_log_pressure_factor(E, μ, T)
    end
    return Float64(degeneracy) * integral / (2.0 * π^2)
end

function stable_meson_pressure_summary(
    pi_mass::Float64,
    k_mass::Float64,
    T::Float64;
    μ_pi::Float64=0.0,
    μ_K::Float64=0.0,
    d_pi::Integer=meson_degeneracy(:pi),
    d_K::Integer=meson_degeneracy(:K),
    qmax_pi::Union{Nothing,Float64}=nothing,
    qmax_K::Union{Nothing,Float64}=nothing,
    num_q_nodes::Int=256,
)
    p_pi = stable_meson_pressure(
        pi_mass,
        T;
        μ=μ_pi,
        degeneracy=Int(d_pi),
        qmax=qmax_pi,
        num_q_nodes=num_q_nodes,
    )
    p_K = stable_meson_pressure(
        k_mass,
        T;
        μ=μ_K,
        degeneracy=Int(d_K),
        qmax=qmax_K,
        num_q_nodes=num_q_nodes,
    )
    total = p_pi + p_K
    return (
        P_pi=p_pi,
        P_K=p_K,
        P_meson=total,
        P_K_over_P_pi=iszero(p_pi) ? NaN : p_K / p_pi,
        μ_pi=μ_pi,
        μ_K=μ_K,
        d_pi=Int(d_pi),
        d_K=Int(d_K),
        qmax_pi=qmax_pi,
        qmax_K=qmax_K,
        num_q_nodes=num_q_nodes,
    )
end

function strict_bw_meson_pressure(
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
    mass > μ || throw(ArgumentError("strict BW pressure helper requires mass > μ, got mass=$(mass), μ=$(μ)"))

    T == 0.0 && return (
        pressure=0.0,
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
        pressure = stable_meson_pressure(
            mass,
            T;
            μ=μ,
            degeneracy=degeneracy,
            qmax=qmax,
            num_q_nodes=q_nodes,
        )
        return (
            pressure=pressure,
            q_integral_estimate=pressure / Float64(degeneracy),
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
            omega_val += dω_w[iω] * bosonic_log_pressure_factor(ω, μ, T) * _strict_bw_kernel(Δω, gamma)
        end
        omega_val /= (2.0 * π)
        q_shell = (q^2 / (2.0 * π^2)) * omega_val
        q_shell_weighted_sum += q_w[iq] * q_shell
        if iq == length(q_grid)
            q_shell_at_qmax = q_shell
        end
    end

    return (
        pressure=Float64(degeneracy) * q_shell_weighted_sum,
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

function strict_bw_meson_pressure_summary(
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
    pi_pressure = strict_bw_meson_pressure(
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
    k_pressure = strict_bw_meson_pressure(
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

    p_pi = Float64(pi_pressure.pressure)
    p_K = Float64(k_pressure.pressure)
    return (
        P_pi=p_pi,
        P_K=p_K,
        P_meson=p_pi + p_K,
        P_K_over_P_pi=iszero(p_pi) ? NaN : p_K / p_pi,
        pi_pressure=pi_pressure,
        k_pressure=k_pressure,
        qmax=qmax,
        q_nodes=q_nodes,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        gamma_zero_tol=gamma_zero_tol,
    )
end

function phase_shift_meson_pressure(
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

    tp = _normalize_thermo_namedtuple(thermo_params)
    _require_nonnegative("temperature T", Float64(tp.T))
    _require_phase_shift_isotropic_xi(Float64(tp.ξ))
    Float64(tp.T) > 0.0 || return (
        meson=meson,
        pressure=0.0,
        q_integral_estimate=0.0,
        omega_shell_at_qmax=0.0,
        qmax=qmax,
        q_nodes=q_nodes,
        omega_min=omega_lower,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        degeneracy=Int(degeneracy),
        scheme=scheme_sym,
    )

    qp = ensure_quark_params_has_A(_normalize_quark_namedtuple(quark_params), tp)
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
            omega_val += omega_w[iω] * gω * _phase_shift_weighted_phase(phase_unwrapped[iω], scheme_sym)
        end
        omega_val /= (2.0 * π)
        q_shell = (q^2 / (2.0 * π^2)) * omega_val
        q_shell_weighted_sum += q_w[iq] * q_shell
        if iq == length(q_grid)
            q_shell_at_qmax = q_shell
        end
    end

    pressure = Float64(degeneracy) * q_shell_weighted_sum
    return (
        meson=meson,
        pressure=pressure,
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

function phase_shift_meson_pressure_summary(
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
    pi_pressure = phase_shift_meson_pressure(
        pi_channel,
        quark_params,
        thermo_params;
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
    k_pressure = phase_shift_meson_pressure(
        k_channel,
        quark_params,
        thermo_params;
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

    p_pi = Float64(pi_pressure.pressure)
    p_K = Float64(k_pressure.pressure)
    return (
        T_fm=Float64(thermo_params.T),
        xi=Float64(thermo_params.ξ),
        pi_channel=pi_channel,
        k_channel=k_channel,
        d_pi=Int(d_pi),
        d_K=Int(d_K),
        P_pi=p_pi,
        P_K=p_K,
        P_meson=p_pi + p_K,
        P_K_over_P_pi=iszero(p_pi) ? NaN : p_K / p_pi,
        qmax=qmax,
        q_nodes=q_nodes,
        omega_min=omega_min,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        eta=eta,
        scheme=_phase_shift_scheme_symbol(scheme),
        pi_pressure=pi_pressure,
        k_pressure=k_pressure,
    )
end

end # module MesonThermodynamics
