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
                      DEFAULT_PHASE_SHIFT_PHASE_CONVENTION,
                      DEFAULT_PHASE_SHIFT_REAL_AXIS_MODE,
                      DEFAULT_PHASE_SHIFT_OMEGA_MAX,
                      DEFAULT_PHASE_SHIFT_OMEGA_NODES,
                      DEFAULT_PHASE_SHIFT_Q_MAX,
                      DEFAULT_PHASE_SHIFT_Q_NODES,
                      _build_k_coeffs,
                      _phase_convention_symbol,
                      _phase_shift_omega_lower_bound,
                      _phase_shift_scheme_symbol,
                      _phase_shift_weighted_phase,
                      _propagator_phase,
                      _real_axis_mode_symbol,
                      _require_phase_shift_isotropic_xi,
                      _require_positive_node_count,
                      _resolve_real_axis_config,
                      _strict_bw_kernel,
                      _unwrap_phases,
                      bose_distribution,
                      meson_degeneracy
using Main.Constants_PNJL: Λ_inv_fm

export bosonic_log_pressure_factor
export stable_meson_pressure, stable_meson_pressure_summary
export strict_bw_meson_pressure, strict_bw_meson_pressure_summary
export phase_shift_meson_pressure, phase_shift_meson_pressure_summary

const PHASE_UNWRAP_BRANCH_TOL = 1e-4

@inline function _require_nonnegative(name::AbstractString, value::Real)
    value >= 0.0 && return
    throw(ArgumentError("$(name) must be nonnegative, got $(value)"))
end

@inline function _default_qmax(mass::Real, T::Real, μ::Real)::Float64
    gap = max(Float64(mass - μ), 0.0)
    return max(8.0, 20.0 * Float64(T) + 10.0 * gap, 8.0 * Float64(mass) + 10.0 * Float64(T))
end

@inline function _normalize_thermo_namedtuple(tp)
    return (
        T=tp.T,
        Φ=tp.Φ,
        Φbar=tp.Φbar,
        ξ=tp.ξ,
    )
end

@inline function _normalize_quark_namedtuple(qp)
    base = (
        m=(u=Float64(qp.m.u), d=Float64(qp.m.d), s=Float64(qp.m.s)),
        μ=(u=Float64(qp.μ.u), d=Float64(qp.μ.d), s=Float64(qp.μ.s)),
    )
    if hasproperty(qp, :A)
        return merge(base, (A=(u=qp.A.u, d=qp.A.d, s=qp.A.s),))
    end
    return base
end

@inline function _resolve_channel_degeneracy(meson::Symbol, degeneracy::Union{Nothing,Integer})::Int
    return degeneracy === nothing ? meson_degeneracy(meson) : Int(degeneracy)
end

@inline function _normalize_ld_cutoff_mode(mode::Symbol)::Symbol
    if mode === :match_model_lambda || mode === :match_qmax || mode === :explicit
        return mode
    end
    throw(ArgumentError("ld_cutoff_mode must be :match_model_lambda, :match_qmax or :explicit, got $(mode)"))
end

@inline function _normalize_ld_threshold_mode(mode::Symbol)::Symbol
    if mode === :omega_lt_q || mode === :spacelike_strict
        return mode
    end
    throw(ArgumentError("ld_threshold_mode must be :omega_lt_q or :spacelike_strict, got $(mode)"))
end

@inline function _resolve_ld_cutoff(
    ld_cutoff::Union{Nothing,Real},
    qmax::Float64,
    mode::Symbol,
)::Float64
    mode_sym = _normalize_ld_cutoff_mode(mode)
    if mode_sym === :match_model_lambda
        return ld_cutoff === nothing ? min(Float64(Λ_inv_fm), qmax) : Float64(ld_cutoff)
    end
    if mode_sym === :match_qmax
        return ld_cutoff === nothing ? qmax : Float64(ld_cutoff)
    end
    ld_cutoff === nothing && throw(ArgumentError("ld_cutoff must be provided when ld_cutoff_mode = :explicit"))
    return Float64(ld_cutoff)
end

@inline function _is_ld_region(ω::Real, q::Real, threshold_mode::Symbol)::Bool
    mode_sym = _normalize_ld_threshold_mode(threshold_mode)
    if mode_sym === :omega_lt_q
        return ω < q
    end
    return ω * ω < q * q
end

function bosonic_log_pressure_factor(E::Real, μ::Real, T::Real)
    _require_nonnegative("temperature T", T)
    T == 0.0 && return 0.0
    E > μ || throw(ArgumentError("Bosonic pressure kernel requires E > μ, got E=$(E), μ=$(μ)"))

    exponent = (E - μ) / T
    exponent > 700.0 && return 0.0
    return -T * log1p(-exp(-exponent))
end

function stable_meson_pressure(
    mass::Float64,
    T::Real;
    μ::Real=0.0,
    degeneracy::Integer=1,
    qmax::Union{Nothing,Float64}=nothing,
    num_q_nodes::Int=256,
)
    _require_nonnegative("mass", mass)
    _require_nonnegative("temperature T", T)
    degeneracy > 0 || throw(ArgumentError("degeneracy must be positive, got $(degeneracy)"))
    num_q_nodes > 1 || throw(ArgumentError("num_q_nodes must be > 1, got $(num_q_nodes)"))
    mass > μ || throw(ArgumentError("Stable meson pressure requires mass > μ, got mass=$(mass), μ=$(μ)"))
    T == 0.0 && return 0.0

    q_upper = qmax === nothing ? _default_qmax(mass, Float64(T), Float64(μ)) : Float64(qmax)
    q_upper > 0.0 || throw(ArgumentError("qmax must be positive, got $(q_upper)"))

    nodes, weights = gauleg(0.0, q_upper, num_q_nodes)
    integral = 0.0
    @inbounds for i in eachindex(nodes, weights)
        q = nodes[i]
        E = hypot(q, mass)
        integral += weights[i] * q^2 * bosonic_log_pressure_factor(E, μ, T)
    end
    return degeneracy * integral / (2.0 * π^2)
end

function stable_meson_pressure_summary(
    pi_mass::Float64,
    k_mass::Float64,
    T::Real;
    μ_pi::Real=0.0,
    μ_K::Real=0.0,
    d_pi::Union{Nothing,Integer}=nothing,
    d_K::Union{Nothing,Integer}=nothing,
    pi_channel::Symbol=:pi,
    k_channel::Symbol=:K,
    qmax_pi::Union{Nothing,Float64}=nothing,
    qmax_K::Union{Nothing,Float64}=nothing,
    num_q_nodes::Int=256,
)
    d_pi_resolved = _resolve_channel_degeneracy(pi_channel, d_pi)
    d_K_resolved = _resolve_channel_degeneracy(k_channel, d_K)
    p_pi = stable_meson_pressure(
        pi_mass,
        T;
        μ=μ_pi,
        degeneracy=d_pi_resolved,
        qmax=qmax_pi,
        num_q_nodes=num_q_nodes,
    )
    p_K = stable_meson_pressure(
        k_mass,
        T;
        μ=μ_K,
        degeneracy=d_K_resolved,
        qmax=qmax_K,
        num_q_nodes=num_q_nodes,
    )
    total = p_pi + p_K
    return (
        P_pi=p_pi,
        P_K=p_K,
        P_meson=total,
        P_K_over_P_pi=iszero(p_pi) ? NaN : p_K / p_pi,
        pi_channel=pi_channel,
        k_channel=k_channel,
        μ_pi=μ_pi,
        μ_K=μ_K,
        d_pi=d_pi_resolved,
        d_K=d_K_resolved,
        qmax_pi=qmax_pi === nothing ? _default_qmax(pi_mass, Float64(T), Float64(μ_pi)) : qmax_pi,
        qmax_K=qmax_K === nothing ? _default_qmax(k_mass, Float64(T), Float64(μ_K)) : qmax_K,
        num_q_nodes=num_q_nodes,
    )
end

function strict_bw_meson_pressure(
    mass::Float64,
    gamma::Float64,
    T::Real;
    μ::Real=0.0,
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
        pressure=degeneracy * q_shell_weighted_sum,
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
    T::Real;
    μ_pi::Float64=0.0,
    μ_K::Float64=0.0,
    d_pi::Union{Nothing,Integer}=nothing,
    d_K::Union{Nothing,Integer}=nothing,
    pi_channel::Symbol=:pi,
    k_channel::Symbol=:K,
    qmax::Float64=DEFAULT_PHASE_SHIFT_Q_MAX,
    q_nodes::Int=DEFAULT_PHASE_SHIFT_Q_NODES,
    omega_max::Float64=DEFAULT_PHASE_SHIFT_OMEGA_MAX,
    omega_nodes::Int=DEFAULT_PHASE_SHIFT_OMEGA_NODES,
    gamma_zero_tol::Float64=1e-12,
)
    d_pi_resolved = _resolve_channel_degeneracy(pi_channel, d_pi)
    d_K_resolved = _resolve_channel_degeneracy(k_channel, d_K)
    pi_pressure = strict_bw_meson_pressure(
        pi_mass,
        pi_gamma,
        T;
        μ=μ_pi,
        degeneracy=d_pi_resolved,
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
        degeneracy=d_K_resolved,
        qmax=qmax,
        q_nodes=q_nodes,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        gamma_zero_tol=gamma_zero_tol,
    )

    p_pi = pi_pressure.pressure
    p_K = k_pressure.pressure
    return (
        P_pi=p_pi,
        P_K=p_K,
        P_meson=p_pi + p_K,
        P_K_over_P_pi=iszero(p_pi) ? NaN : p_K / p_pi,
        pi_channel=pi_channel,
        k_channel=k_channel,
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
    μ::Real=0.0,
    scheme::Symbol=:current,
    qmax::Float64=DEFAULT_PHASE_SHIFT_Q_MAX,
    q_nodes::Int=DEFAULT_PHASE_SHIFT_Q_NODES,
    omega_min::Float64=0.05,
    omega_max::Float64=DEFAULT_PHASE_SHIFT_OMEGA_MAX,
    omega_nodes::Int=DEFAULT_PHASE_SHIFT_OMEGA_NODES,
    eta::Float64=DEFAULT_PHASE_SHIFT_ETA,
    real_axis_mode::Symbol=DEFAULT_PHASE_SHIFT_REAL_AXIS_MODE,
    phase_convention::Symbol=DEFAULT_PHASE_SHIFT_PHASE_CONVENTION,
    ld_cutoff::Union{Nothing,Real}=nothing,
    ld_cutoff_mode::Symbol=:match_model_lambda,
    ld_threshold_mode::Symbol=:omega_lt_q,
    phase_unwrap_branch_tol::Real=PHASE_UNWRAP_BRANCH_TOL,
)
    degeneracy > 0 || throw(ArgumentError("degeneracy must be positive, got $(degeneracy)"))
    _require_positive_node_count("q_nodes", q_nodes)
    _require_positive_node_count("omega_nodes", omega_nodes)
    qmax > 0.0 || throw(ArgumentError("qmax must be positive, got $(qmax)"))
    phase_unwrap_branch_tol >= 0.0 || throw(ArgumentError("phase_unwrap_branch_tol must be nonnegative, got $(phase_unwrap_branch_tol)"))
    omega_lower = _phase_shift_omega_lower_bound(omega_min, μ)
    omega_max > omega_lower || throw(ArgumentError("omega_max must exceed effective omega_min=$(omega_lower)"))
    scheme_sym = _phase_shift_scheme_symbol(scheme)
    ld_cutoff_value = _resolve_ld_cutoff(ld_cutoff, qmax, ld_cutoff_mode)
    ld_cutoff_value > 0.0 || throw(ArgumentError("ld_cutoff must be positive, got $(ld_cutoff_value)"))
    ld_cutoff_value <= qmax || throw(ArgumentError("ld_cutoff must not exceed qmax=$(qmax), got $(ld_cutoff_value)"))
    ld_cutoff_mode_sym = _normalize_ld_cutoff_mode(ld_cutoff_mode)
    ld_threshold_mode_sym = _normalize_ld_threshold_mode(ld_threshold_mode)
    axis = _resolve_real_axis_config(_real_axis_mode_symbol(real_axis_mode), eta)
    phase_convention_sym = _phase_convention_symbol(phase_convention)

    tp = _normalize_thermo_namedtuple(thermo_params)
    _require_nonnegative("temperature T", tp.T)
    _require_phase_shift_isotropic_xi(Float64(tp.ξ))
    tp.T > 0.0 || return (
        meson=meson,
        pressure=0.0,
        pressure_qp=0.0,
        pressure_ld=0.0,
        q_integral_estimate=0.0,
        q_integral_estimate_qp=0.0,
        q_integral_estimate_ld=0.0,
        omega_shell_at_qmax=0.0,
        omega_shell_at_qmax_qp=0.0,
        omega_shell_at_qmax_ld=0.0,
        qmax=qmax,
        q_nodes=q_nodes,
        omega_min=omega_lower,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        degeneracy=Int(degeneracy),
        scheme=scheme_sym,
        eta=axis.eta,
        real_axis_mode=axis.mode,
        polarization_backend=axis.polarization_backend,
        phase_convention=phase_convention_sym,
        ld_cutoff=ld_cutoff_value,
        ld_cutoff_mode=ld_cutoff_mode_sym,
        ld_threshold_mode=ld_threshold_mode_sym,
    )

    qp = ensure_quark_params_has_A(_normalize_quark_namedtuple(quark_params), tp)
    K_coeffs = _build_k_coeffs(qp)
    q_grid, q_w = gauleg(0.0, qmax, q_nodes)
    omega_grid, omega_w = gauleg(omega_lower, omega_max, omega_nodes)

    q_shell_weighted_sum = zero(tp.T)
    q_shell_weighted_sum_qp = zero(tp.T)
    q_shell_weighted_sum_ld = zero(tp.T)
    q_shell_at_qmax = NaN
    q_shell_at_qmax_qp = NaN
    q_shell_at_qmax_ld = NaN
    @inbounds for iq in eachindex(q_grid, q_w)
        q = q_grid[iq]
        phases = [
            _propagator_phase(
                meson, ω, q, qp, tp, K_coeffs;
                eta=axis.eta,
                real_axis_mode=axis.mode,
                phase_convention=phase_convention_sym,
            ) for ω in omega_grid
        ]
        phase_unwrapped = _unwrap_phases(phases; branch_tol=phase_unwrap_branch_tol)
        omega_val = zero(tp.T)
        omega_val_qp = zero(tp.T)
        omega_val_ld = zero(tp.T)
        ld_active = q <= ld_cutoff_value
        for iω in eachindex(omega_grid, omega_w, phase_unwrapped)
            gω = bose_distribution(omega_grid[iω], μ, tp.T)
            weighted_phase = _phase_shift_weighted_phase(phase_unwrapped[iω], scheme_sym)
            contribution = omega_w[iω] * gω * weighted_phase
            is_ld_region = _is_ld_region(omega_grid[iω], q, ld_threshold_mode_sym)
            if ld_active && is_ld_region
                omega_val_ld += contribution
                omega_val += contribution
            elseif !is_ld_region
                omega_val_qp += contribution
                omega_val += contribution
            end
        end
        omega_val /= (2.0 * π)
        omega_val_qp /= (2.0 * π)
        omega_val_ld /= (2.0 * π)
        q_shell = (q^2 / (2.0 * π^2)) * omega_val
        q_shell_qp = (q^2 / (2.0 * π^2)) * omega_val_qp
        q_shell_ld = (q^2 / (2.0 * π^2)) * omega_val_ld
        q_shell_weighted_sum += q_w[iq] * q_shell
        q_shell_weighted_sum_qp += q_w[iq] * q_shell_qp
        q_shell_weighted_sum_ld += q_w[iq] * q_shell_ld
        if iq == length(q_grid)
            q_shell_at_qmax = q_shell
            q_shell_at_qmax_qp = q_shell_qp
            q_shell_at_qmax_ld = q_shell_ld
        end
    end

    pressure_qp = degeneracy * q_shell_weighted_sum_qp
    pressure_ld = degeneracy * q_shell_weighted_sum_ld
    pressure = pressure_qp + pressure_ld
    return (
        meson=meson,
        pressure=pressure,
        pressure_qp=pressure_qp,
        pressure_ld=pressure_ld,
        q_integral_estimate=q_shell_weighted_sum,
        q_integral_estimate_qp=q_shell_weighted_sum_qp,
        q_integral_estimate_ld=q_shell_weighted_sum_ld,
        omega_shell_at_qmax=q_shell_at_qmax,
        omega_shell_at_qmax_qp=q_shell_at_qmax_qp,
        omega_shell_at_qmax_ld=q_shell_at_qmax_ld,
        qmax=qmax,
        q_nodes=q_nodes,
        omega_min=omega_lower,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        degeneracy=Int(degeneracy),
        scheme=scheme_sym,
        eta=axis.eta,
        real_axis_mode=axis.mode,
        polarization_backend=axis.polarization_backend,
        phase_convention=phase_convention_sym,
        ld_cutoff=ld_cutoff_value,
        ld_cutoff_mode=ld_cutoff_mode_sym,
        ld_threshold_mode=ld_threshold_mode_sym,
    )
end

function phase_shift_meson_pressure_summary(
    quark_params,
    thermo_params;
    pi_channel::Symbol=:pi,
    k_channel::Symbol=:K,
    μ_pi::Float64=0.0,
    μ_K::Float64=0.0,
    d_pi::Union{Nothing,Integer}=nothing,
    d_K::Union{Nothing,Integer}=nothing,
    scheme::Symbol=:current,
    qmax::Float64=DEFAULT_PHASE_SHIFT_Q_MAX,
    q_nodes::Int=DEFAULT_PHASE_SHIFT_Q_NODES,
    omega_min::Float64=0.05,
    omega_max::Float64=DEFAULT_PHASE_SHIFT_OMEGA_MAX,
    omega_nodes::Int=DEFAULT_PHASE_SHIFT_OMEGA_NODES,
    eta::Float64=DEFAULT_PHASE_SHIFT_ETA,
    real_axis_mode::Symbol=DEFAULT_PHASE_SHIFT_REAL_AXIS_MODE,
    phase_convention::Symbol=DEFAULT_PHASE_SHIFT_PHASE_CONVENTION,
    ld_cutoff::Union{Nothing,Real}=nothing,
    ld_cutoff_mode::Symbol=:match_model_lambda,
    ld_threshold_mode::Symbol=:omega_lt_q,
    phase_unwrap_branch_tol::Real=PHASE_UNWRAP_BRANCH_TOL,
)
    d_pi_resolved = _resolve_channel_degeneracy(pi_channel, d_pi)
    d_K_resolved = _resolve_channel_degeneracy(k_channel, d_K)
    pi_pressure = phase_shift_meson_pressure(
        pi_channel,
        quark_params,
        thermo_params;
        μ=μ_pi,
        degeneracy=d_pi_resolved,
        scheme=scheme,
        qmax=qmax,
        q_nodes=q_nodes,
        omega_min=omega_min,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        eta=eta,
        real_axis_mode=real_axis_mode,
        phase_convention=phase_convention,
        ld_cutoff=ld_cutoff,
        ld_cutoff_mode=ld_cutoff_mode,
        ld_threshold_mode=ld_threshold_mode,
        phase_unwrap_branch_tol=phase_unwrap_branch_tol,
    )
    k_pressure = phase_shift_meson_pressure(
        k_channel,
        quark_params,
        thermo_params;
        μ=μ_K,
        degeneracy=d_K_resolved,
        scheme=scheme,
        qmax=qmax,
        q_nodes=q_nodes,
        omega_min=omega_min,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        eta=eta,
        real_axis_mode=real_axis_mode,
        phase_convention=phase_convention,
        ld_cutoff=ld_cutoff,
        ld_cutoff_mode=ld_cutoff_mode,
        ld_threshold_mode=ld_threshold_mode,
        phase_unwrap_branch_tol=phase_unwrap_branch_tol,
    )

    p_pi = pi_pressure.pressure
    p_K = k_pressure.pressure
    p_pi_qp = pi_pressure.pressure_qp
    p_pi_ld = pi_pressure.pressure_ld
    p_K_qp = k_pressure.pressure_qp
    p_K_ld = k_pressure.pressure_ld
    return (
        T_fm=thermo_params.T,
        xi=thermo_params.ξ,
        pi_channel=pi_channel,
        k_channel=k_channel,
        d_pi=d_pi_resolved,
        d_K=d_K_resolved,
        P_pi=p_pi,
        P_K=p_K,
        P_pi_qp=p_pi_qp,
        P_pi_ld=p_pi_ld,
        P_K_qp=p_K_qp,
        P_K_ld=p_K_ld,
        P_meson_qp=p_pi_qp + p_K_qp,
        P_meson_ld=p_pi_ld + p_K_ld,
        P_meson=p_pi + p_K,
        P_K_over_P_pi=iszero(p_pi) ? NaN : p_K / p_pi,
        qmax=qmax,
        q_nodes=q_nodes,
        omega_min=omega_min,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        eta=pi_pressure.eta,
        real_axis_mode=pi_pressure.real_axis_mode,
        polarization_backend=pi_pressure.polarization_backend,
        phase_convention=pi_pressure.phase_convention,
        scheme=_phase_shift_scheme_symbol(scheme),
        ld_cutoff=pi_pressure.ld_cutoff,
        ld_cutoff_mode=pi_pressure.ld_cutoff_mode,
        ld_threshold_mode=pi_pressure.ld_threshold_mode,
        pi_pressure=pi_pressure,
        k_pressure=k_pressure,
    )
end

end # module MesonThermodynamics
