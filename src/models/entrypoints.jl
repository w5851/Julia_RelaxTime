"""entrypoints.jl

Models 统一流程入口（阶段 C）：
- 扫描：run_tmu_scan / run_trho_scan
- 工作流：solve_gap_and_transport / solve_gap_and_meson_point

当前采用统一入口策略：
- 业务流程模块统一位于 `src/models/*` 域并由 `Models` 暴露入口；
- 调用方应通过 `Models` 模块显式访问工作流能力。
"""

export run_tmu_scan, run_trho_scan, build_default_rho_grid
export scan_workflow_migration_map, scan_workflow_migration_status
export default_scan_numeric_options, solve_pnjl_point
export solve_gap_and_transport, solve_transport_from_equilibrium
export solve_gap_and_meson_point
export solve_gas_liquid_point
export solve_rotation_point
export run_phase_pipeline, run_production_phase_pipeline, find_cep, build_phase_artifacts
export resolve_phase_output_target, promote_phase_artifacts
export normalize_pm_seed_pair, pm_next_seed_source, derive_pm_seed_pair, analyze_pm_branch_competition
export transport_workflow_module, meson_workflow_module
export rotation_workflow_module
export gas_liquid_workflow_module
export workflow_param_adapters_module
export pnjl_module
export magnetic_thermodynamics_module

@inline function _transport_workflow_module()
    workflow = TransportWorkflow
    if !(isdefined(workflow, :solve_gap_and_transport) && isdefined(workflow, :solve_transport_from_equilibrium))
        error("TransportWorkflow module loaded but required API is missing")
    end
    return workflow
end

@inline function _meson_workflow_module()
    workflow = MesonMassWorkflow
    if !isdefined(workflow, :solve_gap_and_meson_point)
        error("MesonMassWorkflow module loaded but required API (solve_gap_and_meson_point) is missing")
    end
    return workflow
end

@inline function scan_workflow_migration_map()
    return [
        (
            old_entry="scripts/pnjl/run_tmu_scan.jl",
            new_entry="Models.run_tmu_scan(...; model_kind=...)",
            status=:active,
            removal_wave=:D,
            removal_threshold="script callers fully migrate to model-driven entrypoint and parity regression remains stable",
        ),
        (
            old_entry="scripts/pnjl/run_dense_trho_scan.jl",
            new_entry="Models.run_trho_scan(...; model_kind=...)",
            status=:active,
            removal_wave=:D,
            removal_threshold="script callers fully migrate to model-driven entrypoint and parity regression remains stable",
        ),
        (
            old_entry="scripts/pnjl/run_adaptive_trho_scan.jl",
            new_entry="Models.run_trho_scan(...; model_kind=...)",
            status=:active,
            removal_wave=:D,
            removal_threshold="workflow routes no longer rely on script SOP forks and smoke/regression stay stable",
        ),
    ]
end

function scan_workflow_migration_status(old_entry::AbstractString)
    for entry in scan_workflow_migration_map()
        if entry.old_entry == old_entry
            return (
                old_entry=entry.old_entry,
                unified_entry=entry.new_entry,
                status=entry.status,
                removal_wave=entry.removal_wave,
                removal_threshold=entry.removal_threshold,
                route=:compat_adapter,
                deprecation_ready=true,
            )
        end
    end
    throw(ArgumentError("unknown scan/workflow compatibility entry: $(old_entry)"))
end

function run_tmu_scan(args...; kwargs...)
    return TmuScan.run_tmu_scan(args...; kwargs...)
end

function run_trho_scan(args...; kwargs...)
    return TrhoScan.run_trho_scan(args...; kwargs...)
end

function build_default_rho_grid(args...; kwargs...)
    return TrhoScan.build_default_rho_grid(args...; kwargs...)
end

@inline function default_scan_numeric_options()
    return (p_num=24, t_num=8)
end

@inline function _resolve_numeric_option(value, fallback::Int, name::Symbol)
    resolved = value === nothing ? fallback : Int(value)
    resolved > 0 || throw(ArgumentError("$(name) must be positive, got $(resolved)"))
    return resolved
end

@inline function _resolve_real_option(value, name::Symbol)
    val = Float64(value)
    isfinite(val) || throw(ArgumentError("$(name) must be finite, got $(value)"))
    return val
end

@inline function _resolve_temperature_mev(T_mev, t_mev)
    chosen = T_mev === nothing ? t_mev : T_mev
    chosen === nothing && throw(ArgumentError("Missing required parameter: T_mev"))
    return _resolve_real_option(chosen, :T_mev)
end

@inline function _resolve_mu_mev(mu_mev, mu, rho_target)
    chosen = mu_mev === nothing ? mu : mu_mev
    if rho_target === nothing
        chosen === nothing && throw(ArgumentError("Missing required parameter: mu_mev (for FixedMu mode)"))
    elseif chosen === nothing
        return 0.0
    end
    return _resolve_real_option(chosen, :mu_mev)
end

function solve_pnjl_point(;
    T_mev=nothing,
    t_mev=nothing,
    mu_mev=nothing,
    mu=nothing,
    rho_target=nothing,
    xi=0.0,
    p_num=nothing,
    t_num=nothing,
    allow_seed_fallback=true,
)
    numeric = default_scan_numeric_options()
    resolved_p_num = _resolve_numeric_option(p_num, numeric.p_num, :p_num)
    resolved_t_num = _resolve_numeric_option(t_num, numeric.t_num, :t_num)
    resolved_xi = _resolve_real_option(xi, :xi)
    resolved_T_mev = _resolve_temperature_mev(T_mev, t_mev)
    resolved_mu_mev = _resolve_mu_mev(mu_mev, mu, rho_target)
    T_fm = resolved_T_mev / Main.Constants_PNJL.ħc_MeV_fm

    model = create_model(:PNJL)

    if rho_target !== nothing
        resolved_rho_target = _resolve_real_option(rho_target, :rho_target)
        mu_seed_fm = resolved_mu_mev / Main.Constants_PNJL.ħc_MeV_fm
        seed_fallback_used = false
        seed_state = try
            st_seed = solve_gap(model, T_fm, mu_seed_fm; xi=resolved_xi, p_num=resolved_p_num, t_num=resolved_t_num)
            state_vector(st_seed)
        catch err
            if !Bool(allow_seed_fallback)
                rethrow(err)
            end
            seed_fallback_used = true
            @warn "PNJL fixed-rho seed solve failed; falling back to default seed" T_mev=resolved_T_mev mu_seed_mev=resolved_mu_mev rho_target=resolved_rho_target xi=resolved_xi exception=(err, catch_backtrace())
            [0.02, 0.02, 0.03, 0.5, 0.5]
        end

        solver_primary = NLsolveGapSolver(method=:newton, jacobian=:forward, xtol=1e-9, ftol=1e-9)
        solver_secondary = NLsolveGapSolver(method=:trust_region, jacobian=:forward, xtol=1e-9, ftol=1e-9)
        r = solve_fixedrho_constraint(
            model,
            T_fm,
            resolved_rho_target;
            seed_guess=seed_state,
            mu0=mu_seed_fm,
            solver_primary=solver_primary,
            solver_secondary=solver_secondary,
            xi=resolved_xi,
            p_num=resolved_p_num,
            t_num=resolved_t_num,
            residual_norm_max=1e-6,
        )
        return (
            converged=r.converged,
            omega=r.omega,
            pressure=r.pressure,
            rho_norm=r.rho_norm,
            entropy=r.entropy,
            energy=r.energy,
            iterations=r.iterations,
            residual_norm=r.residual_norm,
            xi=resolved_xi,
            seed_fallback_used=seed_fallback_used,
            x_state=collect(r.x_state),
            mu_vec=collect(r.mu_vec),
            masses=collect(r.masses),
        )
    end

    mu_fm = resolved_mu_mev / Main.Constants_PNJL.ħc_MeV_fm
    st = solve_gap(model, T_fm, mu_fm; xi=resolved_xi, p_num=resolved_p_num, t_num=resolved_t_num)
    x_state = state_vector(st)
    mu_vec = normalize_mu_vec(mu_fm)
    omega_val = omega(model, x_state, T_fm, mu_vec; p_num=resolved_p_num, t_num=resolved_t_num, xi=resolved_xi)
    pressure = -omega_val
    rho = number_densities(model, x_state, T_fm, mu_vec; p_num=resolved_p_num, t_num=resolved_t_num, xi=resolved_xi)
    rho_norm = sum(rho.quark) / 3.0
    masses = calculate_mass_vec(model, x_state)
    return (
        converged=true,
        omega=omega_val,
        pressure=pressure,
        rho_norm=rho_norm,
        entropy=NaN,
        energy=NaN,
        iterations=0,
        residual_norm=NaN,
        xi=resolved_xi,
        seed_fallback_used=false,
        x_state=collect(x_state),
        mu_vec=collect(mu_vec),
        masses=collect(masses),
    )
end

function solve_gap_and_transport(args...; kwargs...)
    return _transport_workflow_module().solve_gap_and_transport(args...; kwargs...)
end

function solve_transport_from_equilibrium(args...; kwargs...)
    return _transport_workflow_module().solve_transport_from_equilibrium(args...; kwargs...)
end

function solve_gap_and_meson_point(args...; kwargs...)
    return _meson_workflow_module().solve_gap_and_meson_point(args...; kwargs...)
end

function solve_gas_liquid_point(args...; kwargs...)
    return gas_liquid_workflow_module().solve_gas_liquid_point(args...; kwargs...)
end

function solve_rotation_point(args...; kwargs...)
    return rotation_workflow_module().solve_rotation_point(args...; kwargs...)
end

@inline transport_workflow_module() = _transport_workflow_module()
@inline meson_workflow_module() = _meson_workflow_module()
@inline gas_liquid_workflow_module() = GasLiquidWorkflow
@inline rotation_workflow_module() = RotationWorkflow
@inline pnjl_module() = @__MODULE__
@inline magnetic_thermodynamics_module() = MagneticThermodynamics

@inline function workflow_param_adapters_module()
    adapters = WorkflowParamAdapters
    if !isdefined(adapters, :normalize_quark_params)
        error("WorkflowParamAdapters module loaded but required API (normalize_quark_params) is missing")
    end
    return adapters
end
