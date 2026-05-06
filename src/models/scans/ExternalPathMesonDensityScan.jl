"""
    ExternalPathMesonDensityScan

沿显式外部离散路径点列组织介子数密度 workflow 扫描。

当前定位：
- 为 Friesen 2019 phase-line 一类外部 `(T, μ_B)` 点列提供统一 workflow 入口；
- 保留 continuation_state 契约，避免脚本层自行拼装求解顺序；
- 与 freezeout / crossover 扫描并列，作为“外部路径驱动”的正式扫描入口。
"""
module ExternalPathMesonDensityScan

using Main.Constants_PNJL: ħc_MeV_fm
using ..FlavorChemicalProfiles
using ..MesonChemicalProfiles
using ..CrossoverMesonDensityScan
using ..ScanCommon

export run_external_path_meson_density_scan, DEFAULT_EXTERNAL_PATH_MESON_DENSITY_OUTPUT_PATH

const DEFAULT_EXTERNAL_PATH_MESON_DENSITY_OUTPUT_PATH = normpath(joinpath(
    @__DIR__, "..", "..", "..",
    "data", "outputs", "results", "relaxtime", "meson_density", "external_path",
    "external_path_meson_density_scan.csv",
))

const HEADER = join((
    "path_source",
    "path_case_id",
    "path_line_style",
    "path_point_index",
    "path_order_key",
    "muq_MeV",
    "muB_MeV",
    "T_MeV",
    "T_over_muB",
    "xi",
    "flavor_chemical_profile",
    "meson_chemical_profile",
    "regime",
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

@inline _hasprop(x, name::Symbol) = x isa NamedTuple ? haskey(x, name) : hasproperty(x, name)
@inline _prop(x, name::Symbol) = x isa NamedTuple ? getfield(x, name) : getproperty(x, name)
@inline _get(x, name::Symbol, default) = _hasprop(x, name) ? _prop(x, name) : default

@inline function _required_real(point, field::Symbol)
    _hasprop(point, field) || throw(ArgumentError("external path point missing required field: $(field)"))
    value = Float64(_prop(point, field))
    isfinite(value) || throw(ArgumentError("external path point field $(field) must be finite, got $(value)"))
    return value
end

@inline function _resolve_muB_MeV(point)
    if _hasprop(point, :muB_MeV)
        return _required_real(point, :muB_MeV)
    elseif _hasprop(point, :muB_GeV)
        return 1000.0 * _required_real(point, :muB_GeV)
    end
    throw(ArgumentError("external path point requires muB_MeV or muB_GeV"))
end

@inline function _resolve_T_MeV(point)
    if _hasprop(point, :T_MeV)
        return _required_real(point, :T_MeV)
    elseif _hasprop(point, :T_GeV)
        return 1000.0 * _required_real(point, :T_GeV)
    end
    throw(ArgumentError("external path point requires T_MeV or T_GeV"))
end

@inline function _normalized_point(point, idx::Int)
    muB_MeV = _resolve_muB_MeV(point)
    T_MeV = _resolve_T_MeV(point)
    return (
        path_source=String(_get(point, :source_fig, _get(point, :path_source, "external_path"))),
        path_case_id=String(_get(point, :case_id, _get(point, :path_case_id, "default_case"))),
        path_line_style=String(_get(point, :line_style, _get(point, :path_line_style, "solid"))),
        path_point_index=Int(round(Float64(_get(point, :point_index, idx)))),
        path_order_key=Float64(_get(point, :path_order_key, muB_MeV)),
        muB_MeV=muB_MeV,
        T_MeV=T_MeV,
    )
end

@inline _path_group_key(pt) = (pt.path_source, pt.path_case_id, pt.path_line_style)

@inline function _write_failure_row(io, pt, xi::Float64, flavor_profile, flavor_chemical, chem_profile, regime::Symbol, message::AbstractString)
    muq_MeV = pt.muB_MeV / 3.0
    x = pt.muB_MeV > 0.0 ? pt.T_MeV / pt.muB_MeV : NaN
    values = (
        pt.path_source,
        pt.path_case_id,
        pt.path_line_style,
        string(pt.path_point_index),
        ScanCommon.fmt(pt.path_order_key),
        ScanCommon.fmt(muq_MeV),
        ScanCommon.fmt(pt.muB_MeV),
        ScanCommon.fmt(pt.T_MeV),
        ScanCommon.fmt(x),
        ScanCommon.fmt(xi),
        flavor_profile.profile_name,
        chem_profile.profile_name,
        String(regime),
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

@inline function _write_success_row(io, pt, xi::Float64, flavor_profile, flavor_chemical, chem_profile, regime::Symbol, result)
    density = CrossoverMesonDensityScan._density_payload(result, regime)
    qp = result.quark_params
    tp = result.thermo_params
    meson_results = result.meson_results
    m_pi = Float64(density.m_pi) * ħc_MeV_fm
    m_K = Float64(density.m_K) * ħc_MeV_fm
    gamma_pi = haskey(meson_results, density.pi_channel) ? Float64(meson_results[density.pi_channel].gamma) * ħc_MeV_fm : NaN
    gamma_K = haskey(meson_results, density.k_channel) ? Float64(meson_results[density.k_channel].gamma) * ħc_MeV_fm : NaN
    muq_MeV = pt.muB_MeV / 3.0
    x = pt.muB_MeV > 0.0 ? pt.T_MeV / pt.muB_MeV : NaN
    values = (
        pt.path_source,
        pt.path_case_id,
        pt.path_line_style,
        string(pt.path_point_index),
        ScanCommon.fmt(pt.path_order_key),
        ScanCommon.fmt(muq_MeV),
        ScanCommon.fmt(pt.muB_MeV),
        ScanCommon.fmt(pt.T_MeV),
        ScanCommon.fmt(x),
        ScanCommon.fmt(xi),
        flavor_profile.profile_name,
        chem_profile.profile_name,
        String(regime),
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
        string(CrossoverMesonDensityScan._equilibrium_bool(result, :converged)),
        hasproperty(result.equilibrium, :iterations) ? string(getproperty(result.equilibrium, :iterations)) : "-1",
        ScanCommon.fmt(CrossoverMesonDensityScan._equilibrium_real(result, :residual_norm)),
        "",
    )
    println(io, join(values, ','))
end

function run_external_path_meson_density_scan(;
    points,
    xi::Float64=0.0,
    flavor_chemical_profile_name::AbstractString="default",
    meson_chemical_profile_name::AbstractString="default",
    regime::Symbol=:stable,
    output_path::AbstractString=DEFAULT_EXTERNAL_PATH_MESON_DENSITY_OUTPUT_PATH,
    overwrite::Bool=true,
    p_num::Int=24,
    t_num::Int=8,
    max_iter::Int=40,
    stable_num_q_nodes::Int=CrossoverMesonDensityScan.DEFAULT_MESON_DENSITY_Q_NODES,
    strict_bw_qmax::Float64=CrossoverMesonDensityScan.DEFAULT_PHASE_SHIFT_Q_MAX,
    strict_bw_q_nodes::Int=CrossoverMesonDensityScan.DEFAULT_PHASE_SHIFT_Q_NODES,
    strict_bw_omega_max::Float64=CrossoverMesonDensityScan.DEFAULT_PHASE_SHIFT_OMEGA_MAX,
    strict_bw_omega_nodes::Int=CrossoverMesonDensityScan.DEFAULT_PHASE_SHIFT_OMEGA_NODES,
    phase_shift_qmax::Float64=CrossoverMesonDensityScan.DEFAULT_PHASE_SHIFT_Q_MAX,
    phase_shift_q_nodes::Int=CrossoverMesonDensityScan.DEFAULT_PHASE_SHIFT_Q_NODES,
    phase_shift_omega_min::Float64=0.05,
    phase_shift_omega_max::Float64=CrossoverMesonDensityScan.DEFAULT_PHASE_SHIFT_OMEGA_MAX,
    phase_shift_omega_nodes::Int=CrossoverMesonDensityScan.DEFAULT_PHASE_SHIFT_OMEGA_NODES,
    phase_shift_eta::Float64=1e-6,
)
    raw_points = collect(points)
    isempty(raw_points) && throw(ArgumentError("points must not be empty"))
    normalized = [_normalized_point(pt, idx) for (idx, pt) in enumerate(raw_points)]
    sort!(normalized; by=pt -> (pt.path_source, pt.path_case_id, pt.path_line_style, pt.path_order_key, pt.path_point_index))
    regime_sym = CrossoverMesonDensityScan._regime_symbol(regime)

    flavor_profile = FlavorChemicalProfiles.load_flavor_chemical_profile(profile=String(flavor_chemical_profile_name))
    chemical_profile = MesonChemicalProfiles.load_meson_chemical_profile(profile=String(meson_chemical_profile_name))

    mkpath(dirname(output_path))
    !overwrite && isfile(output_path) && throw(ArgumentError("output_path already exists: $(output_path)"))

    solver_kwargs = (iterations=max_iter,)
    mass_kwargs = (iterations=max_iter,)
    continuation_state = nothing
    previous_group_key = nothing
    open(output_path, "w") do io
        println(io, HEADER)
        for pt in normalized
            group_key = _path_group_key(pt)
            if previous_group_key !== nothing && group_key != previous_group_key
                continuation_state = nothing
            end
            previous_group_key = group_key
            muq_MeV = pt.muB_MeV / 3.0
            muq_fm = muq_MeV / ħc_MeV_fm
            T_fm = pt.T_MeV / ħc_MeV_fm
            flavor_chemical = FlavorChemicalProfiles.flavor_mu_profile_fm(flavor_profile, muq_MeV)
            try
                result, flavor_chemical = CrossoverMesonDensityScan._solve_density_point(
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
                _write_success_row(io, pt, xi, flavor_profile, flavor_chemical, chemical_profile, regime_sym, result)
            catch err
                _write_failure_row(io, pt, xi, flavor_profile, flavor_chemical, chemical_profile, regime_sym, sprint(showerror, err))
            end
        end
    end

    return (
        output_path=output_path,
        points=length(normalized),
        xi=xi,
        flavor_chemical_profile=flavor_profile.profile_name,
        meson_chemical_profile=chemical_profile.profile_name,
        regime=regime_sym,
        workflow_entry="Models.run_external_path_meson_density_scan",
    )
end

end # module ExternalPathMesonDensityScan
