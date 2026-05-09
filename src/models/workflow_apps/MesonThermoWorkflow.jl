module MesonThermoWorkflow

"""
meson thermo workflow 最小实现：

1. 复用 `solve_gap_and_meson_point` 作为唯一 meson 求解入口；
2. 在 workflow 层补齐 `P_meson / P_total / entropy / epsilon / trace_anomaly`；
3. 暂不抽象成新的公共总线，先服务 `π/K` canonical thermo 主线。
"""

const _RELAX_TIME_MODULE_PATH = normpath(joinpath(@__DIR__, "..", "..", "relaxtime", "RelaxTime.jl"))
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RELAX_TIME_MODULE_PATH)
end

using ..MesonMassWorkflow: solve_gap_and_meson_point
using ..WorkflowParamAdapters: normalize_quark_params, normalize_thermo_params
using Main.MesonDensity: meson_degeneracy,
                         phase_shift_meson_density_summary,
                         strict_bw_meson_density_summary,
                         stable_meson_number_density
using Main.MesonThermodynamics: phase_shift_meson_pressure_summary,
                                strict_bw_meson_pressure_summary,
                                stable_meson_pressure_summary
using Main.Constants_PNJL: ħc_MeV_fm
using ..Models: create_model, default_momentum_count, default_theta_count, model_thermo

export solve_meson_thermo_from_meson_point
export solve_gap_and_meson_thermo_point
export solve_strict_bw_meson_thermo_from_meson_point
export solve_gap_and_strict_bw_meson_thermo_point
export solve_phase_shift_meson_thermo_from_meson_point
export solve_gap_and_phase_shift_meson_thermo_point
export build_meson_thermo_contract_row

@inline function _require_result_field(result, field::Symbol)
    hasproperty(result, field) || throw(ArgumentError("meson workflow result missing required field: $(field)"))
    return getproperty(result, field)
end

@inline function _require_finite_mass(meson_results, meson::Symbol)::Float64
    haskey(meson_results, meson) || throw(ArgumentError("meson_results missing required channel $(meson)"))
    entry = meson_results[meson]
    hasproperty(entry, :mass) || throw(ArgumentError("meson_results[$(meson)] missing mass field"))
    mass = Float64(entry.mass)
    isfinite(mass) || throw(ArgumentError("meson_results[$(meson)].mass must be finite, got $(mass)"))
    mass >= 0.0 || throw(ArgumentError("meson_results[$(meson)].mass must be nonnegative, got $(mass)"))
    return mass
end

@inline function _require_finite_gamma(meson_results, meson::Symbol)::Float64
    haskey(meson_results, meson) || throw(ArgumentError("meson_results missing required channel $(meson)"))
    entry = meson_results[meson]
    hasproperty(entry, :gamma) || throw(ArgumentError("meson_results[$(meson)] missing gamma field"))
    gamma = Float64(entry.gamma)
    isfinite(gamma) || throw(ArgumentError("meson_results[$(meson)].gamma must be finite, got $(gamma)"))
    gamma >= 0.0 || throw(ArgumentError("meson_results[$(meson)].gamma must be nonnegative, got $(gamma)"))
    return gamma
end

@inline function _background_thermo_from_meson_point(
    meson_point;
    p_num::Int=default_momentum_count(),
    t_num::Int=default_theta_count(),
)
    model = create_model(:PNJL)
    equilibrium = _require_result_field(meson_point, :equilibrium)
    thermo_params = _require_result_field(meson_point, :thermo_params)
    pressure, rho_norm, entropy, energy = model_thermo(
        model,
        equilibrium.x_state,
        equilibrium.mu_vec,
        Float64(thermo_params.T);
        p_num=p_num,
        t_num=t_num,
        xi=Float64(thermo_params.ξ),
    )
    return (
        P_quark_meanfield=Float64(pressure),
        rho_norm=Float64(rho_norm),
        entropy_quark_meanfield=Float64(entropy),
        epsilon_quark_meanfield=Float64(energy),
    )
end

@inline function _common_contract_fields(
    meson_point,
    meson_pressure::Float64,
    density_summary,
    workflow::Symbol,
    phase_shift_variant,
    background;
    channel_set::AbstractString="pi,K",
    P_meson_qp::Float64=NaN,
    P_meson_ld::Float64=NaN,
    P_pi_qp::Float64=NaN,
    P_pi_ld::Float64=NaN,
    P_K_qp::Float64=NaN,
    P_K_ld::Float64=NaN,
    ld_cutoff=nothing,
    ld_cutoff_mode=nothing,
    ld_threshold_mode=nothing,
)
    thermo_params = _require_result_field(meson_point, :thermo_params)
    equilibrium = _require_result_field(meson_point, :equilibrium)
    P_bg = Float64(background.P_quark_meanfield)
    P_total = P_bg + meson_pressure
    return (
        T_fm=Float64(thermo_params.T),
        xi=Float64(thermo_params.ξ),
        workflow=workflow,
        channel_set=String(channel_set),
        P_meson=meson_pressure,
        P_meson_qp=P_meson_qp,
        P_meson_ld=P_meson_ld,
        P_quark_meanfield=P_bg,
        P_total=P_total,
        delta_P_vs_no_meson=meson_pressure,
        P_meson_over_P_total=iszero(P_total) ? NaN : meson_pressure / P_total,
        P_pi_qp=P_pi_qp,
        P_pi_ld=P_pi_ld,
        P_K_qp=P_K_qp,
        P_K_ld=P_K_ld,
        μ_pi=Float64(getproperty(density_summary, :μ_pi)),
        μ_K=Float64(getproperty(density_summary, :μ_K)),
        n_pi=Float64(getproperty(density_summary, :n_pi)),
        n_K=Float64(getproperty(density_summary, :n_K)),
        equilibrium_converged=Bool(equilibrium.converged),
        phase_structure=:unknown,
        phase_shift_variant=phase_shift_variant,
        ld_cutoff=ld_cutoff,
        ld_cutoff_mode=ld_cutoff_mode,
        ld_threshold_mode=ld_threshold_mode,
        entropy_quark_meanfield=Float64(background.entropy_quark_meanfield),
        epsilon_quark_meanfield=Float64(background.epsilon_quark_meanfield),
        rho_norm=Float64(background.rho_norm),
    )
end

function _attach_derived_thermo(
    base_result,
    pressure_at_temperature::Function;
    temperature_step_fm::Float64=(1.0 / ħc_MeV_fm),
)
    T_fm = Float64(base_result.T_fm)
    temperature_step_fm > 0.0 || throw(ArgumentError("temperature_step_fm must be positive, got $(temperature_step_fm)"))

    if T_fm == 0.0
        entropy_meson = 0.0
    elseif T_fm > temperature_step_fm
        P_plus = Float64(pressure_at_temperature(T_fm + temperature_step_fm))
        P_minus = Float64(pressure_at_temperature(T_fm - temperature_step_fm))
        entropy_meson = (P_plus - P_minus) / (2.0 * temperature_step_fm)
    else
        P_plus = Float64(pressure_at_temperature(T_fm + temperature_step_fm))
        entropy_meson = (P_plus - Float64(base_result.P_meson)) / temperature_step_fm
    end

    μ_term = Float64(base_result.μ_pi) * Float64(base_result.n_pi) +
             Float64(base_result.μ_K) * Float64(base_result.n_K)
    epsilon_meson = -Float64(base_result.P_meson) + T_fm * entropy_meson + μ_term
    entropy_total = Float64(base_result.entropy_quark_meanfield) + entropy_meson
    epsilon_total = Float64(base_result.epsilon_quark_meanfield) + epsilon_meson
    trace_anomaly = T_fm == 0.0 ? NaN : (epsilon_total - 3.0 * Float64(base_result.P_total)) / T_fm^4

    return merge(base_result, (
        entropy_meson=entropy_meson,
        entropy=entropy_total,
        epsilon_meson=epsilon_meson,
        epsilon=epsilon_total,
        trace_anomaly=trace_anomaly,
    ))
end

function solve_meson_thermo_from_meson_point(
    meson_point;
    pi_channel::Symbol=:pi,
    k_channel::Symbol=:K,
    μ_pi::Real=0.0,
    μ_K::Real=0.0,
    d_pi::Integer=meson_degeneracy(:pi),
    d_K::Integer=meson_degeneracy(:K),
    qmax_pi::Union{Nothing,Real}=nothing,
    qmax_K::Union{Nothing,Real}=nothing,
    num_q_nodes::Int=256,
    p_num::Int=default_momentum_count(),
    t_num::Int=default_theta_count(),
)
    meson_results = _require_result_field(meson_point, :meson_results)
    m_pi = _require_finite_mass(meson_results, pi_channel)
    m_K = _require_finite_mass(meson_results, k_channel)
    thermo_params = _require_result_field(meson_point, :thermo_params)

    density_summary = (
        μ_pi=Float64(μ_pi),
        μ_K=Float64(μ_K),
        n_pi=stable_meson_number_density(m_pi, Float64(thermo_params.T); μ=Float64(μ_pi), degeneracy=Int(d_pi), qmax=qmax_pi === nothing ? nothing : Float64(qmax_pi), num_q_nodes=num_q_nodes),
        n_K=stable_meson_number_density(m_K, Float64(thermo_params.T); μ=Float64(μ_K), degeneracy=Int(d_K), qmax=qmax_K === nothing ? nothing : Float64(qmax_K), num_q_nodes=num_q_nodes),
    )
    pressure_summary = stable_meson_pressure_summary(
        m_pi,
        m_K,
        Float64(thermo_params.T);
        μ_pi=Float64(μ_pi),
        μ_K=Float64(μ_K),
        d_pi=Int(d_pi),
        d_K=Int(d_K),
        qmax_pi=qmax_pi === nothing ? nothing : Float64(qmax_pi),
        qmax_K=qmax_K === nothing ? nothing : Float64(qmax_K),
        num_q_nodes=num_q_nodes,
    )
    background = _background_thermo_from_meson_point(meson_point; p_num=p_num, t_num=t_num)
    return merge(
        _common_contract_fields(
            meson_point,
            Float64(pressure_summary.P_meson),
            density_summary,
            :stable_meson_pressure,
            nothing,
            background,
        ),
        (
            pi_channel=pi_channel,
            k_channel=k_channel,
            m_pi=m_pi,
            m_K=m_K,
            d_pi=Int(d_pi),
            d_K=Int(d_K),
            qmax_pi=qmax_pi === nothing ? nothing : Float64(qmax_pi),
            qmax_K=qmax_K === nothing ? nothing : Float64(qmax_K),
            num_q_nodes=num_q_nodes,
            meson_pressure=pressure_summary,
            meson_density=density_summary,
        ),
    )
end

function solve_strict_bw_meson_thermo_from_meson_point(
    meson_point;
    pi_channel::Symbol=:pi,
    k_channel::Symbol=:K,
    μ_pi::Real=0.0,
    μ_K::Real=0.0,
    d_pi::Integer=meson_degeneracy(:pi),
    d_K::Integer=meson_degeneracy(:K),
    qmax::Float64=12.0,
    q_nodes::Int=48,
    omega_max::Float64=10.0,
    omega_nodes::Int=48,
    gamma_zero_tol::Float64=1e-12,
    p_num::Int=default_momentum_count(),
    t_num::Int=default_theta_count(),
)
    meson_results = _require_result_field(meson_point, :meson_results)
    m_pi = _require_finite_mass(meson_results, pi_channel)
    m_K = _require_finite_mass(meson_results, k_channel)
    gamma_pi = _require_finite_gamma(meson_results, pi_channel)
    gamma_K = _require_finite_gamma(meson_results, k_channel)

    pressure_summary = strict_bw_meson_pressure_summary(
        m_pi,
        gamma_pi,
        m_K,
        gamma_K,
        Float64(meson_point.thermo_params.T);
        μ_pi=Float64(μ_pi),
        μ_K=Float64(μ_K),
        d_pi=Int(d_pi),
        d_K=Int(d_K),
        qmax=qmax,
        q_nodes=q_nodes,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        gamma_zero_tol=gamma_zero_tol,
    )
    density_summary = merge(
        strict_bw_meson_density_summary(
            m_pi,
            gamma_pi,
            m_K,
            gamma_K,
            Float64(meson_point.thermo_params.T);
            μ_pi=Float64(μ_pi),
            μ_K=Float64(μ_K),
            d_pi=Int(d_pi),
            d_K=Int(d_K),
            qmax=qmax,
            q_nodes=q_nodes,
            omega_max=omega_max,
            omega_nodes=omega_nodes,
            gamma_zero_tol=gamma_zero_tol,
        ),
        (
            μ_pi=Float64(μ_pi),
            μ_K=Float64(μ_K),
        ),
    )
    background = _background_thermo_from_meson_point(meson_point; p_num=p_num, t_num=t_num)
    return merge(
        _common_contract_fields(
            meson_point,
            Float64(pressure_summary.P_meson),
            density_summary,
            :strict_bw_stage1_reduced_pressure,
            nothing,
            background,
        ),
        (
            pi_channel=pi_channel,
            k_channel=k_channel,
            m_pi=m_pi,
            m_K=m_K,
            gamma_pi=gamma_pi,
            gamma_K=gamma_K,
            d_pi=Int(d_pi),
            d_K=Int(d_K),
            qmax=qmax,
            q_nodes=q_nodes,
            omega_max=omega_max,
            omega_nodes=omega_nodes,
            gamma_zero_tol=gamma_zero_tol,
            meson_pressure=pressure_summary,
            meson_density=density_summary,
        ),
    )
end

function solve_phase_shift_meson_thermo_from_meson_point(
    meson_point;
    pi_channel::Symbol=:pi,
    k_channel::Symbol=:K,
    μ_pi::Float64=0.0,
    μ_K::Float64=0.0,
    d_pi::Integer=meson_degeneracy(:pi),
    d_K::Integer=meson_degeneracy(:K),
    scheme::Symbol=:current,
    qmax::Float64=12.0,
    q_nodes::Int=48,
    omega_min::Float64=0.05,
    omega_max::Float64=10.0,
    omega_nodes::Int=48,
    eta::Float64=1e-6,
    ld_cutoff::Union{Nothing,Float64}=nothing,
    ld_cutoff_mode::Symbol=:match_qmax,
    ld_threshold_mode::Symbol=:omega_lt_q,
    p_num::Int=default_momentum_count(),
    t_num::Int=default_theta_count(),
)
    quark_params = normalize_quark_params(_require_result_field(meson_point, :quark_params))
    thermo_params = normalize_thermo_params(_require_result_field(meson_point, :thermo_params))

    pressure_summary = phase_shift_meson_pressure_summary(
        quark_params,
        thermo_params;
        pi_channel=pi_channel,
        k_channel=k_channel,
        μ_pi=μ_pi,
        μ_K=μ_K,
        d_pi=Int(d_pi),
        d_K=Int(d_K),
        scheme=scheme,
        qmax=qmax,
        q_nodes=q_nodes,
        omega_min=omega_min,
        omega_max=omega_max,
        omega_nodes=omega_nodes,
        eta=eta,
        ld_cutoff=ld_cutoff,
        ld_cutoff_mode=ld_cutoff_mode,
        ld_threshold_mode=ld_threshold_mode,
    )
    density_summary = merge(
        phase_shift_meson_density_summary(
            quark_params,
            thermo_params;
            pi_channel=pi_channel,
            k_channel=k_channel,
            μ_pi=μ_pi,
            μ_K=μ_K,
            d_pi=Int(d_pi),
            d_K=Int(d_K),
            scheme=scheme,
            qmax=qmax,
            q_nodes=q_nodes,
            omega_min=omega_min,
            omega_max=omega_max,
            omega_nodes=omega_nodes,
            eta=eta,
        ),
        (
            μ_pi=Float64(μ_pi),
            μ_K=Float64(μ_K),
        ),
    )
    background = _background_thermo_from_meson_point(meson_point; p_num=p_num, t_num=t_num)
    return merge(
        _common_contract_fields(
            meson_point,
            Float64(pressure_summary.P_meson),
            density_summary,
            Symbol(pressure_summary.scheme),
            Symbol(pressure_summary.scheme),
            background;
            channel_set=string(pi_channel, ",", k_channel),
            P_meson_qp=Float64(pressure_summary.P_meson_qp),
            P_meson_ld=Float64(pressure_summary.P_meson_ld),
            P_pi_qp=Float64(pressure_summary.P_pi_qp),
            P_pi_ld=Float64(pressure_summary.P_pi_ld),
            P_K_qp=Float64(pressure_summary.P_K_qp),
            P_K_ld=Float64(pressure_summary.P_K_ld),
            ld_cutoff=Float64(pressure_summary.ld_cutoff),
            ld_cutoff_mode=pressure_summary.ld_cutoff_mode,
            ld_threshold_mode=pressure_summary.ld_threshold_mode,
        ),
        (
            pi_channel=pi_channel,
            k_channel=k_channel,
            m_pi=haskey(meson_point.meson_results, pi_channel) ? Float64(meson_point.meson_results[pi_channel].mass) : NaN,
            m_K=haskey(meson_point.meson_results, k_channel) ? Float64(meson_point.meson_results[k_channel].mass) : NaN,
            d_pi=Int(d_pi),
            d_K=Int(d_K),
            qmax=qmax,
            q_nodes=q_nodes,
            omega_min=omega_min,
            omega_max=omega_max,
            omega_nodes=omega_nodes,
            eta=eta,
            ld_cutoff=Float64(pressure_summary.ld_cutoff),
            ld_cutoff_mode=pressure_summary.ld_cutoff_mode,
            ld_threshold_mode=pressure_summary.ld_threshold_mode,
            meson_pressure=pressure_summary,
            meson_density=density_summary,
        ),
    )
end

function solve_gap_and_meson_thermo_point(
    T_fm::Real,
    mu_fm::Real;
    thermo_kwargs::NamedTuple=(;),
    derive_thermo::Bool=true,
    temperature_step_fm::Float64=(1.0 / ħc_MeV_fm),
    kwargs...,
)
    meson_point = solve_gap_and_meson_point(T_fm, mu_fm; kwargs...)
    result = solve_meson_thermo_from_meson_point(meson_point; thermo_kwargs...)
    if !derive_thermo
        return merge(meson_point, (meson_thermo=result,))
    end

    pressure_at_temperature = T_probe -> begin
        probe = solve_gap_and_meson_thermo_point(
            T_probe,
            mu_fm;
            thermo_kwargs=thermo_kwargs,
            derive_thermo=false,
            kwargs...,
        )
        return Float64(probe.meson_thermo.P_meson)
    end
    derived = _attach_derived_thermo(result, pressure_at_temperature; temperature_step_fm=temperature_step_fm)
    return merge(meson_point, (meson_thermo=derived,))
end

function solve_gap_and_strict_bw_meson_thermo_point(
    T_fm::Real,
    mu_fm::Real;
    thermo_kwargs::NamedTuple=(;),
    derive_thermo::Bool=true,
    temperature_step_fm::Float64=(1.0 / ħc_MeV_fm),
    kwargs...,
)
    meson_point = solve_gap_and_meson_point(T_fm, mu_fm; kwargs...)
    result = solve_strict_bw_meson_thermo_from_meson_point(meson_point; thermo_kwargs...)
    if !derive_thermo
        return merge(meson_point, (strict_bw_meson_thermo=result,))
    end

    pressure_at_temperature = T_probe -> begin
        probe = solve_gap_and_strict_bw_meson_thermo_point(
            T_probe,
            mu_fm;
            thermo_kwargs=thermo_kwargs,
            derive_thermo=false,
            kwargs...,
        )
        return Float64(probe.strict_bw_meson_thermo.P_meson)
    end
    derived = _attach_derived_thermo(result, pressure_at_temperature; temperature_step_fm=temperature_step_fm)
    return merge(meson_point, (strict_bw_meson_thermo=derived,))
end

function solve_gap_and_phase_shift_meson_thermo_point(
    T_fm::Real,
    mu_fm::Real;
    thermo_kwargs::NamedTuple=(;),
    derive_thermo::Bool=true,
    temperature_step_fm::Float64=(1.0 / ħc_MeV_fm),
    kwargs...,
)
    meson_point = solve_gap_and_meson_point(T_fm, mu_fm; kwargs...)
    result = solve_phase_shift_meson_thermo_from_meson_point(meson_point; thermo_kwargs...)
    if !derive_thermo
        return merge(meson_point, (phase_shift_meson_thermo=result,))
    end

    pressure_at_temperature = T_probe -> begin
        probe = solve_gap_and_phase_shift_meson_thermo_point(
            T_probe,
            mu_fm;
            thermo_kwargs=thermo_kwargs,
            derive_thermo=false,
            kwargs...,
        )
        return Float64(probe.phase_shift_meson_thermo.P_meson)
    end
    derived = _attach_derived_thermo(result, pressure_at_temperature; temperature_step_fm=temperature_step_fm)
    return merge(meson_point, (phase_shift_meson_thermo=derived,))
end

function build_meson_thermo_contract_row(point_result)
    result = hasproperty(point_result, :meson_thermo) ? point_result.meson_thermo :
             hasproperty(point_result, :strict_bw_meson_thermo) ? point_result.strict_bw_meson_thermo :
             hasproperty(point_result, :phase_shift_meson_thermo) ? point_result.phase_shift_meson_thermo :
             point_result
    equilibrium = _require_result_field(point_result, :equilibrium)
    T_fm = Float64(result.T_fm)
    return (
        T_MeV=T_fm * ħc_MeV_fm,
        muB_MeV=Float64(equilibrium.mu_vec[1]) * 3.0 * ħc_MeV_fm,
        workflow=String(Symbol(result.workflow)),
        channel_set=String(result.channel_set),
        P_meson=Float64(result.P_meson),
        P_meson_qp=hasproperty(result, :P_meson_qp) ? Float64(result.P_meson_qp) : NaN,
        P_meson_ld=hasproperty(result, :P_meson_ld) ? Float64(result.P_meson_ld) : NaN,
        P_total=Float64(result.P_total),
        P_quark_meanfield=Float64(result.P_quark_meanfield),
        epsilon=hasproperty(result, :epsilon) ? Float64(result.epsilon) : NaN,
        entropy=hasproperty(result, :entropy) ? Float64(result.entropy) : NaN,
        trace_anomaly=hasproperty(result, :trace_anomaly) ? Float64(result.trace_anomaly) : NaN,
        P_meson_over_P_total=Float64(result.P_meson_over_P_total),
        delta_P_vs_no_meson=Float64(result.delta_P_vs_no_meson),
        P_pi_qp=hasproperty(result, :P_pi_qp) ? Float64(result.P_pi_qp) : NaN,
        P_pi_ld=hasproperty(result, :P_pi_ld) ? Float64(result.P_pi_ld) : NaN,
        P_K_qp=hasproperty(result, :P_K_qp) ? Float64(result.P_K_qp) : NaN,
        P_K_ld=hasproperty(result, :P_K_ld) ? Float64(result.P_K_ld) : NaN,
        equilibrium_converged=Bool(result.equilibrium_converged),
        phase_structure=String(Symbol(result.phase_structure)),
        phase_shift_variant=result.phase_shift_variant === nothing ? "" : String(Symbol(result.phase_shift_variant)),
        ld_cutoff=hasproperty(result, :ld_cutoff) && result.ld_cutoff !== nothing ? Float64(result.ld_cutoff) : NaN,
        ld_cutoff_mode=hasproperty(result, :ld_cutoff_mode) && result.ld_cutoff_mode !== nothing ? String(Symbol(result.ld_cutoff_mode)) : "",
        ld_threshold_mode=hasproperty(result, :ld_threshold_mode) && result.ld_threshold_mode !== nothing ? String(Symbol(result.ld_threshold_mode)) : "",
    )
end

end # module
