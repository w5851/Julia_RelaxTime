"""entrypoints.jl

Models 统一流程入口（阶段 C）：
- 扫描：run_tmu_scan / run_trho_scan
- 工作流：solve_gap_and_transport / solve_gap_and_meson_point

当前采用统一入口策略：
- 业务流程模块统一位于 `src/models/*` 域并由 `Models` 暴露入口；
- 调用方应通过 `Models` 模块显式访问工作流能力。
"""

export run_tmu_scan, run_trho_scan, build_default_rho_grid
export solve_gap_and_transport, solve_transport_from_equilibrium
export solve_gap_and_meson_point
export run_phase_pipeline, find_cep, build_phase_artifacts
export resolve_phase_output_target, promote_phase_artifacts
export normalize_pm_seed_pair, pm_next_seed_source, derive_pm_seed_pair, analyze_pm_branch_competition
export _pm_branch_scan_fieldnames, _pm_branch_scan_record
export _pm_maxwell_reference_from_rows
export _pm_refine_transition_bracket
export _pm_interpolate_transition_mu, _pm_extract_endpoints, _pm_has_bistable_window
export _pm_compare_with_maxwell, _pm_pressure_crosscheck
export _pm_accept_branch_point, _pm_check_branch_continuity
export transport_workflow_module, meson_workflow_module
export workflow_param_adapters_module
export pnjl_module

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

function run_tmu_scan(args...; kwargs...)
    return TmuScan.run_tmu_scan(args...; kwargs...)
end

function run_trho_scan(args...; kwargs...)
    return TrhoScan.run_trho_scan(args...; kwargs...)
end

function build_default_rho_grid(args...; kwargs...)
    return TrhoScan.build_default_rho_grid(args...; kwargs...)
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

@inline transport_workflow_module() = _transport_workflow_module()
@inline meson_workflow_module() = _meson_workflow_module()
@inline pnjl_module() = @__MODULE__

@inline function workflow_param_adapters_module()
    adapters = WorkflowParamAdapters
    if !isdefined(adapters, :normalize_quark_params)
        error("WorkflowParamAdapters module loaded but required API (normalize_quark_params) is missing")
    end
    return adapters
end
