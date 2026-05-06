"""
    CrossoverMesonDensityScan

沿手征 crossover line 组织介子数密度 workflow 扫描。

当前定位：
- 为 Friesen 2019 一类 `K^+/pi^+` / `K^-/pi^-` 曲线提供正式 workflow 入口；
- 继续复用 MesonMassWorkflow 的 continuation_state，避免脚本层并行拼装。
"""
module CrossoverMesonDensityScan

using Main.Constants_PNJL: ħc_MeV_fm
using Main.MesonDensity: DEFAULT_MESON_DENSITY_Q_NODES,
                         DEFAULT_PHASE_SHIFT_Q_MAX,
                         DEFAULT_PHASE_SHIFT_Q_NODES,
                         DEFAULT_PHASE_SHIFT_OMEGA_MAX,
                         DEFAULT_PHASE_SHIFT_OMEGA_NODES

using ..FlavorChemicalProfiles
using ..MesonChemicalProfiles
using ..MesonDensityWorkflow
using ..ScanCommon

export run_crossover_meson_density_scan, DEFAULT_CROSSOVER_MESON_DENSITY_OUTPUT_PATH

const DEFAULT_CROSSOVER_MESON_DENSITY_OUTPUT_PATH = normpath(joinpath(
    @__DIR__, "..", "..", "..",
    "data", "outputs", "results", "relaxtime", "meson_density", "crossover",
    "crossover_meson_density_scan.csv",
))

const HEADER = join((
    "muq_MeV",
    "muB_MeV",
    "T_MeV",
    "T_over_muB",
    "xi",
    "flavor_chemical_profile",
    "meson_chemical_profile",
    "regime",
    "crossover_method",
    "crossover_variable",
    "crossover_converged",
    "crossover_derivative",
    "rho",
    "mu_u_MeV",
    "mu_d_MeV",
    "mu_s_MeV",
    "pi_label",
    "k_label",
    "charge_resolved",
    "mu_pi_MeV",
    "mu_K_MeV",
    "d_pi",
    "d_K",
    "Phi",
    "Phibar",
    "m_u",
    "m_d",
    "m_s",
    "m_pi_MeV",
    "m_K_MeV",
    "gamma_pi_MeV",
    "gamma_K_MeV",
    "n_pi",
    "n_K",
    "kpi_ratio",
    "equilibrium_converged",
    "equilibrium_iterations",
    "equilibrium_residual_norm",
    "message",
), ",")

@inline function _regime_symbol(regime::Symbol)::Symbol
    if regime === :stable || regime === :stable_density
        return :stable
    elseif regime === :strict_bw_stage1 || regime === :strict_bw_reduced || regime === :stage1_reduced
        return :strict_bw_stage1
    elseif regime === :strict_bw_stage2 || regime === :strict_bw_qpole || regime === :stage2_qpole
        return :strict_bw_stage2
    elseif regime === :phase_shift_current || regime === :current || regime === :phase_e3
        return :phase_shift_current
    elseif regime === :phase_shift_gbu || regime === :gbu || regime === :generalized_bu
        return :phase_shift_gbu
    end
    throw(ArgumentError("unsupported crossover meson-density regime: $(regime)"))
end

@inline function _density_kwargs_for_profile(profile::MesonChemicalProfiles.MesonChemicalProfile, flavor_chemical)
    chemical = MesonChemicalProfiles.meson_chemical_profile_fm(profile; flavor_mev=flavor_chemical)
    return (
        pi_channel=chemical.pi_channel,
        k_channel=chemical.k_channel,
        μ_pi=chemical.mu_pi_fm,
        μ_K=chemical.mu_K_fm,
        d_pi=chemical.d_pi,
        d_K=chemical.d_K,
    )
end

@inline function _equilibrium_bool(result, field::Symbol)
    eq = getproperty(result, :equilibrium)
    return hasproperty(eq, field) ? Bool(getproperty(eq, field)) : false
end

@inline function _equilibrium_real(result, field::Symbol)
    eq = getproperty(result, :equilibrium)
    return hasproperty(eq, field) ? Float64(getproperty(eq, field)) : NaN
end

@inline function _density_payload(result, regime::Symbol)
    if regime === :stable
        return result.meson_density
    elseif regime === :strict_bw_stage1 || regime === :strict_bw_stage2
        return result.strict_bw_meson_density
    end
    return result.phase_shift_meson_density
end

function _solve_density_point(
    T_fm::Float64,
    muq_fm::Float64,
    xi::Float64,
    regime::Symbol,
    continuation_state,
    flavor_profile::FlavorChemicalProfiles.FlavorChemicalProfile,
    chemical_profile::MesonChemicalProfiles.MesonChemicalProfile;
    p_num::Int,
    t_num::Int,
    solver_kwargs::NamedTuple,
    mass_kwargs::NamedTuple,
    stable_num_q_nodes::Int,
    strict_bw_qmax::Float64,
    strict_bw_q_nodes::Int,
    strict_bw_omega_max::Float64,
    strict_bw_omega_nodes::Int,
    phase_shift_qmax::Float64,
    phase_shift_q_nodes::Int,
    phase_shift_omega_min::Float64,
    phase_shift_omega_max::Float64,
    phase_shift_omega_nodes::Int,
    phase_shift_eta::Float64,
)
    flavor_chemical = FlavorChemicalProfiles.flavor_mu_profile_fm(
        flavor_profile,
        muq_fm * ħc_MeV_fm,
    )
    common_density = _density_kwargs_for_profile(chemical_profile, flavor_chemical)
    flavor_override = flavor_profile.apply_to_equilibrium ? (
        flavor_chemical.mu_u_fm,
        flavor_chemical.mu_d_fm,
        flavor_chemical.mu_s_fm,
    ) : nothing
    common_kwargs = (
        xi=xi,
        mesons=(common_density.pi_channel, common_density.k_channel),
        continuation_state=continuation_state,
        mixed_branch_align=:strict_sign_binding,
        flavor_mu_override=flavor_override,
        p_num=p_num,
        t_num=t_num,
        solver_kwargs=solver_kwargs,
        mass_kwargs=mass_kwargs,
    )

    if regime === :stable
        result = MesonDensityWorkflow.solve_gap_and_meson_density_point(
            T_fm,
            muq_fm;
            common_kwargs...,
            density_kwargs=(; common_density..., num_q_nodes=stable_num_q_nodes),
        )
    elseif regime === :strict_bw_stage1
        result = MesonDensityWorkflow.solve_gap_and_strict_bw_meson_density_point(
            T_fm,
            muq_fm;
            common_kwargs...,
            density_kwargs=(;
                common_density...,
                stage=:stage1_reduced,
                qmax=strict_bw_qmax,
                q_nodes=strict_bw_q_nodes,
                omega_max=strict_bw_omega_max,
                omega_nodes=strict_bw_omega_nodes,
            ),
        )
    elseif regime === :strict_bw_stage2
        result = MesonDensityWorkflow.solve_gap_and_strict_bw_meson_density_point(
            T_fm,
            muq_fm;
            common_kwargs...,
            density_kwargs=(;
                common_density...,
                stage=:stage2_qpole,
                qmax=strict_bw_qmax,
                q_nodes=strict_bw_q_nodes,
                omega_max=strict_bw_omega_max,
                omega_nodes=strict_bw_omega_nodes,
            ),
        )
    elseif regime === :phase_shift_current
        result = MesonDensityWorkflow.solve_gap_and_phase_shift_meson_density_point(
            T_fm,
            muq_fm;
            common_kwargs...,
            density_kwargs=(;
                common_density...,
                scheme=:current,
                qmax=phase_shift_qmax,
                q_nodes=phase_shift_q_nodes,
                omega_min=phase_shift_omega_min,
                omega_max=phase_shift_omega_max,
                omega_nodes=phase_shift_omega_nodes,
                eta=phase_shift_eta,
            ),
        )
    else
        result = MesonDensityWorkflow.solve_gap_and_phase_shift_meson_density_point(
            T_fm,
            muq_fm;
            common_kwargs...,
            density_kwargs=(;
                common_density...,
                scheme=:generalized_bu,
                qmax=phase_shift_qmax,
                q_nodes=phase_shift_q_nodes,
                omega_min=phase_shift_omega_min,
                omega_max=phase_shift_omega_max,
                omega_nodes=phase_shift_omega_nodes,
                eta=phase_shift_eta,
            ),
        )
    end

    return result, flavor_chemical
end

function _write_failure_row(io, row, xi::Float64, flavor_profile, flavor_chemical, chem_profile, regime::Symbol, message::AbstractString)
    muq_MeV = Float64(row.mu_MeV)
    muB_MeV = 3.0 * muq_MeV
    T_MeV = Float64(row.T_crossover_MeV)
    x = muB_MeV > 0.0 ? T_MeV / muB_MeV : NaN
    values = (
        ScanCommon.fmt(muq_MeV),
        ScanCommon.fmt(muB_MeV),
        ScanCommon.fmt(T_MeV),
        ScanCommon.fmt(x),
        ScanCommon.fmt(xi),
        flavor_profile.profile_name,
        chem_profile.profile_name,
        String(regime),
        String(row.method),
        String(row.variable),
        string(Bool(row.converged)),
        ScanCommon.fmt(row.derivative),
        ScanCommon.fmt(row.rho),
        ScanCommon.fmt(flavor_chemical.mu_u_MeV),
        ScanCommon.fmt(flavor_chemical.mu_d_MeV),
        ScanCommon.fmt(flavor_chemical.mu_s_MeV),
        chem_profile.pi_label,
        chem_profile.k_label,
        string(chem_profile.charge_resolved),
        ScanCommon.fmt(chem_profile.mu_pi_MeV),
        ScanCommon.fmt(chem_profile.mu_K_MeV),
        string(chem_profile.d_pi),
        string(chem_profile.d_K),
        "NaN", "NaN",
        "NaN", "NaN", "NaN",
        "NaN", "NaN",
        "NaN", "NaN",
        "NaN", "NaN", "NaN",
        "false", "-1", "NaN",
        ScanCommon.quote_csv(ScanCommon.clean_message(message)),
    )
    println(io, join(values, ','))
end

function _write_success_row(io, row, xi::Float64, flavor_profile, flavor_chemical, chem_profile, regime::Symbol, result)
    density = _density_payload(result, regime)
    qp = result.quark_params
    tp = result.thermo_params
    meson_results = result.meson_results
    m_pi = Float64(density.m_pi) * ħc_MeV_fm
    m_K = Float64(density.m_K) * ħc_MeV_fm
    gamma_pi = haskey(meson_results, density.pi_channel) ? Float64(meson_results[density.pi_channel].gamma) * ħc_MeV_fm : NaN
    gamma_K = haskey(meson_results, density.k_channel) ? Float64(meson_results[density.k_channel].gamma) * ħc_MeV_fm : NaN
    muq_MeV = Float64(row.mu_MeV)
    muB_MeV = 3.0 * muq_MeV
    T_MeV = Float64(row.T_crossover_MeV)
    x = muB_MeV > 0.0 ? T_MeV / muB_MeV : NaN

    values = (
        ScanCommon.fmt(muq_MeV),
        ScanCommon.fmt(muB_MeV),
        ScanCommon.fmt(T_MeV),
        ScanCommon.fmt(x),
        ScanCommon.fmt(xi),
        flavor_profile.profile_name,
        chem_profile.profile_name,
        String(regime),
        String(row.method),
        String(row.variable),
        string(Bool(row.converged)),
        ScanCommon.fmt(row.derivative),
        ScanCommon.fmt(row.rho),
        ScanCommon.fmt(flavor_chemical.mu_u_MeV),
        ScanCommon.fmt(flavor_chemical.mu_d_MeV),
        ScanCommon.fmt(flavor_chemical.mu_s_MeV),
        chem_profile.pi_label,
        chem_profile.k_label,
        string(chem_profile.charge_resolved),
        ScanCommon.fmt(chem_profile.mu_pi_MeV),
        ScanCommon.fmt(chem_profile.mu_K_MeV),
        string(chem_profile.d_pi),
        string(chem_profile.d_K),
        ScanCommon.fmt(tp.Φ),
        ScanCommon.fmt(tp.Φbar),
        ScanCommon.fmt(qp.m.u),
        ScanCommon.fmt(qp.m.d),
        ScanCommon.fmt(qp.m.s),
        ScanCommon.fmt(m_pi),
        ScanCommon.fmt(m_K),
        ScanCommon.fmt(gamma_pi),
        ScanCommon.fmt(gamma_K),
        ScanCommon.fmt(density.n_pi),
        ScanCommon.fmt(density.n_K),
        ScanCommon.fmt(density.kpi_ratio),
        string(_equilibrium_bool(result, :converged)),
        hasproperty(result.equilibrium, :iterations) ? string(getproperty(result.equilibrium, :iterations)) : "-1",
        ScanCommon.fmt(_equilibrium_real(result, :residual_norm)),
        "",
    )
    println(io, join(values, ','))
end

function run_crossover_meson_density_scan(;
    mu_max_MeV::Float64,
    T_min_MeV::Float64,
    T_max_MeV::Float64,
    mu_min_MeV::Float64=10.0,
    xi::Float64=0.0,
    n_mu::Int=12,
    method::Symbol=:peak,
    variable::Symbol=:phi_u,
    model_kind::Symbol=:PNJL,
    solver_backend::Symbol=:models,
    flavor_chemical_profile_name::AbstractString="default",
    meson_chemical_profile_name::AbstractString="default",
    regime::Symbol=:stable,
    output_path::AbstractString=DEFAULT_CROSSOVER_MESON_DENSITY_OUTPUT_PATH,
    overwrite::Bool=true,
    p_num::Int=24,
    t_num::Int=8,
    max_iter::Int=40,
    stable_num_q_nodes::Int=DEFAULT_MESON_DENSITY_Q_NODES,
    strict_bw_qmax::Float64=DEFAULT_PHASE_SHIFT_Q_MAX,
    strict_bw_q_nodes::Int=DEFAULT_PHASE_SHIFT_Q_NODES,
    strict_bw_omega_max::Float64=DEFAULT_PHASE_SHIFT_OMEGA_MAX,
    strict_bw_omega_nodes::Int=DEFAULT_PHASE_SHIFT_OMEGA_NODES,
    phase_shift_qmax::Float64=DEFAULT_PHASE_SHIFT_Q_MAX,
    phase_shift_q_nodes::Int=DEFAULT_PHASE_SHIFT_Q_NODES,
    phase_shift_omega_min::Float64=0.05,
    phase_shift_omega_max::Float64=DEFAULT_PHASE_SHIFT_OMEGA_MAX,
    phase_shift_omega_nodes::Int=DEFAULT_PHASE_SHIFT_OMEGA_NODES,
    phase_shift_eta::Float64=1e-6,
)
    mu_max_MeV > 0.0 || throw(ArgumentError("mu_max_MeV must be positive, got $(mu_max_MeV)"))
    mu_min_MeV >= 0.0 || throw(ArgumentError("mu_min_MeV must be nonnegative, got $(mu_min_MeV)"))
    mu_max_MeV > mu_min_MeV || throw(ArgumentError("mu_max_MeV must exceed mu_min_MeV"))
    T_max_MeV > T_min_MeV || throw(ArgumentError("T_max_MeV must exceed T_min_MeV"))
    n_mu >= 3 || throw(ArgumentError("n_mu must be >= 3, got $(n_mu)"))
    regime_sym = _regime_symbol(regime)

    flavor_profile = FlavorChemicalProfiles.load_flavor_chemical_profile(profile=String(flavor_chemical_profile_name))
    chemical_profile = MesonChemicalProfiles.load_meson_chemical_profile(profile=String(meson_chemical_profile_name))
    line_rows = Main.Models.build_crossover_line(
        mu_max_MeV=mu_max_MeV,
        T_min_MeV=T_min_MeV,
        T_max_MeV=T_max_MeV,
        xi=xi,
        n_mu=n_mu,
        method=method,
        variable=variable,
        model_kind=model_kind,
        solver_backend=solver_backend,
    )

    selected = [
        row for row in line_rows
        if Bool(row.converged) && isfinite(Float64(row.mu_MeV)) && Float64(row.mu_MeV) >= mu_min_MeV &&
           isfinite(Float64(row.T_crossover_MeV)) && Float64(row.T_crossover_MeV) > 0.0
    ]
    isempty(selected) && throw(ArgumentError("no valid crossover points found in requested range"))

    mkpath(dirname(output_path))
    !overwrite && isfile(output_path) && throw(ArgumentError("output_path already exists: $(output_path)"))

    solver_kwargs = (iterations=max_iter,)
    mass_kwargs = (iterations=max_iter,)
    continuation_state = nothing
    open(output_path, "w") do io
        println(io, HEADER)
        for row in selected
            muq_fm = Float64(row.mu_MeV) / ħc_MeV_fm
            T_fm = Float64(row.T_crossover_MeV) / ħc_MeV_fm
            flavor_chemical = FlavorChemicalProfiles.flavor_mu_profile_fm(flavor_profile, Float64(row.mu_MeV))
            try
                result, flavor_chemical = _solve_density_point(
                    T_fm,
                    muq_fm,
                    xi,
                    regime_sym,
                    continuation_state,
                    flavor_profile,
                    chemical_profile;
                    p_num=p_num,
                    t_num=t_num,
                    solver_kwargs=solver_kwargs,
                    mass_kwargs=mass_kwargs,
                    stable_num_q_nodes=stable_num_q_nodes,
                    strict_bw_qmax=strict_bw_qmax,
                    strict_bw_q_nodes=strict_bw_q_nodes,
                    strict_bw_omega_max=strict_bw_omega_max,
                    strict_bw_omega_nodes=strict_bw_omega_nodes,
                    phase_shift_qmax=phase_shift_qmax,
                    phase_shift_q_nodes=phase_shift_q_nodes,
                    phase_shift_omega_min=phase_shift_omega_min,
                    phase_shift_omega_max=phase_shift_omega_max,
                    phase_shift_omega_nodes=phase_shift_omega_nodes,
                    phase_shift_eta=phase_shift_eta,
                )
                continuation_state = hasproperty(result, :continuation_state) ? result.continuation_state : continuation_state
                _write_success_row(io, row, xi, flavor_profile, flavor_chemical, chemical_profile, regime_sym, result)
            catch err
                _write_failure_row(io, row, xi, flavor_profile, flavor_chemical, chemical_profile, regime_sym, sprint(showerror, err))
            end
        end
    end

    return (
        output_path=output_path,
        points=length(selected),
        xi=xi,
        mu_min_MeV=mu_min_MeV,
        mu_max_MeV=mu_max_MeV,
        T_min_MeV=T_min_MeV,
        T_max_MeV=T_max_MeV,
        n_mu=n_mu,
        flavor_chemical_profile=flavor_profile.profile_name,
        meson_chemical_profile=chemical_profile.profile_name,
        regime=regime_sym,
        workflow_entry="Models.run_crossover_meson_density_scan",
    )
end

end # module CrossoverMesonDensityScan
