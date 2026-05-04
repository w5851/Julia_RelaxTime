"""
    FreezeoutMesonDensityScan

沿 freeze-out baseline 路径组织介子数密度 workflow 扫描，并复用
MesonMassWorkflow 的 continuation_state 作为唯一连续性契约。
"""
module FreezeoutMesonDensityScan

using Main.Constants_PNJL: ħc_MeV_fm
using Main.MesonDensity: DEFAULT_MESON_DENSITY_Q_NODES,
                         DEFAULT_PHASE_SHIFT_Q_MAX,
                         DEFAULT_PHASE_SHIFT_Q_NODES,
                         DEFAULT_PHASE_SHIFT_OMEGA_MAX,
                         DEFAULT_PHASE_SHIFT_OMEGA_NODES

using ..FreezeoutProfiles
using ..FreezeoutPathProfiles
using ..FlavorChemicalProfiles
using ..MesonChemicalProfiles
using ..MesonDensityWorkflow
using ..ScanCommon

export run_freezeout_meson_density_scan, DEFAULT_FREEZEOUT_MESON_DENSITY_OUTPUT_PATH

const DEFAULT_FREEZEOUT_MESON_DENSITY_OUTPUT_PATH = normpath(joinpath(
    @__DIR__, "..", "..", "..",
    "data", "outputs", "results", "relaxtime", "meson_density", "freezeout",
    "freezeout_meson_density_scan.csv",
))

const HEADER = join((
    "sqrt_s_NN_GeV",
    "muB_MeV",
    "xi",
    "freezeout_profile",
    "path_profile",
    "path_segment",
    "flavor_chemical_profile",
    "meson_chemical_profile",
    "regime",
    "T_MeV",
    "muq_MeV",
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

@inline function _validate_real_vector(name::Symbol, values)
    values isa AbstractVector{<:Real} || throw(ArgumentError("$(name) must be AbstractVector{<:Real}, got $(typeof(values))"))
    isempty(values) && throw(ArgumentError("$(name) must not be empty"))
    for (i, v) in pairs(values)
        isfinite(Float64(v)) || throw(ArgumentError("$(name)[$(i)] must be finite Real, got $(v)"))
    end
    return nothing
end

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
    throw(ArgumentError("unsupported freezeout meson-density regime: $(regime)"))
end

@inline function _schedule_key(pt, xi, freezeout_profile::String, path_profile::String, flavor_profile::String, chem_profile::String, regime::Symbol)
    return join((
        ScanCommon.fmt(pt.sqrt_s_NN_GeV),
        ScanCommon.fmt(pt.muB_MeV),
        ScanCommon.fmt(xi),
        freezeout_profile,
        path_profile,
        flavor_profile,
        chem_profile,
        String(regime),
    ), '|')
end

function _load_completed_keys(path::AbstractString)
    completed = Set{String}()
    open(path, "r") do io
        first_line = true
        for line in eachline(io)
            if first_line
                first_line = false
                continue
            end
            s = strip(line)
            isempty(s) && continue
            cols = split(s, ',')
            length(cols) < 8 && continue
            push!(completed, join((
                strip(cols[1]),
                strip(cols[2]),
                strip(cols[3]),
                strip(cols[4]),
                strip(cols[5]),
                strip(cols[6]),
                strip(cols[7]),
                strip(cols[8]),
            ), '|'))
        end
    end
    return completed
end

@inline function _density_kwargs_for_profile(profile::MesonChemicalProfiles.MesonChemicalProfile)
    chemical = MesonChemicalProfiles.meson_chemical_profile_fm(profile)
    return (
        μ_pi=chemical.mu_pi_fm,
        μ_K=chemical.mu_K_fm,
        d_pi=chemical.d_pi,
        d_K=chemical.d_K,
    )
end

function _solve_density_point(
    pt,
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
    flavor_chemical = FlavorChemicalProfiles.flavor_mu_profile_fm(flavor_profile, pt.muq_MeV)
    common_density = _density_kwargs_for_profile(chemical_profile)
    flavor_override = flavor_profile.apply_to_equilibrium ? (
        flavor_chemical.mu_u_fm,
        flavor_chemical.mu_d_fm,
        flavor_chemical.mu_s_fm,
    ) : nothing
    common_kwargs = (
        xi=xi,
        mesons=(:pi, :K),
        continuation_state=continuation_state,
        mixed_branch_align=:strict_sign_binding,
        flavor_mu_override=flavor_override,
        p_num=p_num,
        t_num=t_num,
        solver_kwargs=solver_kwargs,
        mass_kwargs=mass_kwargs,
    )

    if regime === :stable
        return MesonDensityWorkflow.solve_gap_and_meson_density_point(
            pt.T_fm,
            pt.muB_fm / 3.0;
            common_kwargs...,
            density_kwargs=(; common_density..., num_q_nodes=stable_num_q_nodes),
        )
    elseif regime === :strict_bw_stage1
        return MesonDensityWorkflow.solve_gap_and_strict_bw_meson_density_point(
            pt.T_fm,
            pt.muB_fm / 3.0;
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
        return MesonDensityWorkflow.solve_gap_and_strict_bw_meson_density_point(
            pt.T_fm,
            pt.muB_fm / 3.0;
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
        return MesonDensityWorkflow.solve_gap_and_phase_shift_meson_density_point(
            pt.T_fm,
            pt.muB_fm / 3.0;
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
    end

    return MesonDensityWorkflow.solve_gap_and_phase_shift_meson_density_point(
        pt.T_fm,
        pt.muB_fm / 3.0;
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

@inline function _density_payload(result, regime::Symbol)
    if regime === :stable
        return result.meson_density
    elseif regime === :strict_bw_stage1 || regime === :strict_bw_stage2
        return result.strict_bw_meson_density
    end
    return result.phase_shift_meson_density
end

@inline function _equilibrium_bool(result, field::Symbol)
    eq = getproperty(result, :equilibrium)
    return hasproperty(eq, field) ? Bool(getproperty(eq, field)) : false
end

@inline function _equilibrium_real(result, field::Symbol)
    eq = getproperty(result, :equilibrium)
    return hasproperty(eq, field) ? Float64(getproperty(eq, field)) : NaN
end

function _write_failure_row(io, pt, xi, freezeout_profile::String, flavor_profile, flavor_chemical, chem_profile, regime::Symbol, message::AbstractString)
    values = (
        ScanCommon.fmt(pt.sqrt_s_NN_GeV),
        ScanCommon.fmt(pt.muB_MeV),
        ScanCommon.fmt(xi),
        freezeout_profile,
        String(pt.path_profile),
        String(pt.path_segment),
        flavor_profile.profile_name,
        chem_profile.profile_name,
        String(regime),
        ScanCommon.fmt(pt.T_MeV),
        ScanCommon.fmt(pt.muq_MeV),
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

function _write_success_row(io, pt, xi, freezeout_profile::String, flavor_profile, flavor_chemical, chem_profile, regime::Symbol, result)
    density = _density_payload(result, regime)
    qp = result.quark_params
    tp = result.thermo_params
    meson_results = result.meson_results
    m_pi = Float64(density.m_pi) * ħc_MeV_fm
    m_K = Float64(density.m_K) * ħc_MeV_fm
    gamma_pi = haskey(meson_results, :pi) ? Float64(meson_results[:pi].gamma) * ħc_MeV_fm : NaN
    gamma_K = haskey(meson_results, :K) ? Float64(meson_results[:K].gamma) * ħc_MeV_fm : NaN

    values = (
        ScanCommon.fmt(pt.sqrt_s_NN_GeV),
        ScanCommon.fmt(pt.muB_MeV),
        ScanCommon.fmt(xi),
        freezeout_profile,
        String(pt.path_profile),
        String(pt.path_segment),
        flavor_profile.profile_name,
        chem_profile.profile_name,
        String(regime),
        ScanCommon.fmt(pt.T_MeV),
        ScanCommon.fmt(pt.muq_MeV),
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

function run_freezeout_meson_density_scan(;
    sqrt_s_NN_values,
    xi_values=[0.0],
    freezeout_profile_name::AbstractString="default",
    path_profile_name::AbstractString="baseline_freezeout",
    flavor_chemical_profile_name::AbstractString="default",
    meson_chemical_profile_name::AbstractString="default",
    regime::Symbol=:stable,
    output_path::AbstractString=DEFAULT_FREEZEOUT_MESON_DENSITY_OUTPUT_PATH,
    traversal::Symbol=:sqrts_ascending,
    overwrite::Bool=false,
    resume::Bool=true,
    p_num::Int=24,
    t_num::Int=8,
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
    progress_cb::Union{Nothing, Function}=nothing,
    solver_kwargs::NamedTuple=(; iterations=40),
    mass_kwargs::NamedTuple=(; iterations=40),
)
    _validate_real_vector(:sqrt_s_NN_values, sqrt_s_NN_values)
    _validate_real_vector(:xi_values, xi_values)
    regime_sym = _regime_symbol(regime)

    freezeout_profile = FreezeoutProfiles.load_freezeout_profile(profile=String(freezeout_profile_name))
    path_profile = FreezeoutPathProfiles.load_freezeout_path_profile(profile=String(path_profile_name))
    flavor_profile = FlavorChemicalProfiles.load_flavor_chemical_profile(profile=String(flavor_chemical_profile_name))
    chemical_profile = MesonChemicalProfiles.load_meson_chemical_profile(profile=String(meson_chemical_profile_name))
    points = FreezeoutPathProfiles.build_freezeout_path_points(
        sqrt_s_NN_values;
        freezeout_profile=freezeout_profile,
        path_profile=path_profile,
        traversal=traversal,
    )

    mkpath(dirname(output_path))
    completed = (resume && !overwrite && isfile(output_path)) ? _load_completed_keys(output_path) : Set{String}()
    io_mode = (overwrite || !isfile(output_path)) ? "w" : "a"

    stats = Dict(:total => 0, :success => 0, :failure => 0, :skipped => 0)

    open(output_path, io_mode) do io
        if io_mode == "w"
            println(io, HEADER)
        end

        for xi_raw in xi_values
            xi = Float64(xi_raw)
            continuation_state = nothing
            for pt in points
                stats[:total] += 1
                key = _schedule_key(pt, xi, freezeout_profile.profile_name, path_profile.profile_name, flavor_profile.profile_name, chemical_profile.profile_name, regime_sym)
                if key in completed
                    stats[:skipped] += 1
                    continue
                end

                try
                    flavor_chemical = FlavorChemicalProfiles.flavor_mu_profile_MeV(flavor_profile, pt.muq_MeV)
                    result = _solve_density_point(
                        pt,
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
                    continuation_state = result.continuation_state
                    _write_success_row(io, pt, xi, freezeout_profile.profile_name, flavor_profile, flavor_chemical, chemical_profile, regime_sym, result)
                    stats[:success] += 1
                    if progress_cb !== nothing
                        try
                            progress_cb((sqrt_s_NN_GeV=pt.sqrt_s_NN_GeV, T_MeV=pt.T_MeV, muB_MeV=pt.muB_MeV, xi=xi), result)
                        catch
                        end
                    end
                catch err
                    flavor_chemical = FlavorChemicalProfiles.flavor_mu_profile_MeV(flavor_profile, pt.muq_MeV)
                    _write_failure_row(io, pt, xi, freezeout_profile.profile_name, flavor_profile, flavor_chemical, chemical_profile, regime_sym, sprint(showerror, err))
                    stats[:failure] += 1
                end

                flush(io)
                push!(completed, key)
            end
        end
    end

    return (
        total=stats[:total],
        success=stats[:success],
        failure=stats[:failure],
        skipped=stats[:skipped],
        output=output_path,
        freezeout_profile=freezeout_profile.profile_name,
        path_profile=path_profile.profile_name,
        flavor_chemical_profile=flavor_profile.profile_name,
        meson_chemical_profile=chemical_profile.profile_name,
        regime=regime_sym,
        traversal=traversal,
    )
end

end # module FreezeoutMesonDensityScan
