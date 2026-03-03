"""entrypoints.jl

Models 统一流程入口（阶段 C）：
- 扫描：run_tmu_scan / run_trho_scan
- 工作流：solve_gap_and_transport / solve_gap_and_meson_point

当前采用兼容转发策略：
- 业务流程模块仍位于 `src/pnjl/*`，但统一经 `Models` 暴露调用入口；
- 后续可以在不改调用方的前提下，把实现逐步迁入 `src/models`。
"""

export run_tmu_scan, run_trho_scan, build_default_rho_grid
export solve_gap_and_transport, solve_transport_from_equilibrium
export solve_gap_and_meson_point
export run_phase_pipeline, find_cep, build_phase_artifacts
export resolve_phase_output_target, promote_phase_artifacts
export transport_workflow_module, meson_workflow_module
export workflow_param_adapters_module
export pnjl_module

# Include-once helper
const _INCLUDE_ONCE_PATH = normpath(joinpath(@__DIR__, "..", "utils", "IncludeOnce.jl"))
if !isdefined(Main, :IncludeOnce)
    Base.include(Main, _INCLUDE_ONCE_PATH)
end
const IncludeOnce = Main.IncludeOnce

const _SCAN_ENTRY_MODULE_PATH = normpath(joinpath(@__DIR__, "scans", "ScanEntrypoints.jl"))
const _TRANSPORT_WORKFLOW_BRIDGE_PATH = normpath(joinpath(@__DIR__, "workflows", "TransportWorkflowBridge.jl"))
const _MESON_WORKFLOW_BRIDGE_PATH = normpath(joinpath(@__DIR__, "workflows", "MesonMassWorkflowBridge.jl"))
const _WORKFLOW_PARAM_ADAPTERS_BRIDGE_PATH = normpath(joinpath(@__DIR__, "workflows", "WorkflowParamAdaptersBridge.jl"))

@inline function _pnjl_module()
    scan = IncludeOnce.include_once!(Main, :ModelsScanEntrypoints, _SCAN_ENTRY_MODULE_PATH)
    pnjl = Base.invokelatest(scan.pnjl_module_ref)
    if !(isdefined(pnjl, :run_tmu_scan) && isdefined(pnjl, :run_trho_scan))
        error("PNJL module loaded but required API (run_tmu_scan/run_trho_scan) is missing")
    end
    return pnjl
end

@inline function _transport_workflow_module()
    bridge = IncludeOnce.include_once!(Main, :ModelsTransportWorkflowBridge, _TRANSPORT_WORKFLOW_BRIDGE_PATH)
    workflow = Base.invokelatest(bridge.module_ref)
    if !(isdefined(workflow, :solve_gap_and_transport) && isdefined(workflow, :solve_transport_from_equilibrium))
        error("TransportWorkflow module loaded but required API is missing")
    end
    return workflow
end

@inline function _meson_workflow_module()
    bridge = IncludeOnce.include_once!(Main, :ModelsMesonMassWorkflowBridge, _MESON_WORKFLOW_BRIDGE_PATH)
    workflow = Base.invokelatest(bridge.module_ref)
    if !isdefined(workflow, :solve_gap_and_meson_point)
        error("MesonMassWorkflow module loaded but required API (solve_gap_and_meson_point) is missing")
    end
    return workflow
end

@inline function _invoke_worldsafe(mod, symbol::Symbol, args...; kwargs...)
    fn = getproperty(mod, symbol)
    return Base.invokelatest(fn, args...; kwargs...)
end

function run_tmu_scan(args...; kwargs...)
    return _invoke_worldsafe(_pnjl_module(), :run_tmu_scan, args...; kwargs...)
end

function run_trho_scan(args...; kwargs...)
    return _invoke_worldsafe(_pnjl_module(), :run_trho_scan, args...; kwargs...)
end

function build_default_rho_grid(args...; kwargs...)
    return _invoke_worldsafe(_pnjl_module(), :build_default_rho_grid, args...; kwargs...)
end

function solve_gap_and_transport(args...; kwargs...)
    return _invoke_worldsafe(_transport_workflow_module(), :solve_gap_and_transport, args...; kwargs...)
end

function solve_transport_from_equilibrium(args...; kwargs...)
    return _invoke_worldsafe(_transport_workflow_module(), :solve_transport_from_equilibrium, args...; kwargs...)
end

function solve_gap_and_meson_point(args...; kwargs...)
    return _invoke_worldsafe(_meson_workflow_module(), :solve_gap_and_meson_point, args...; kwargs...)
end

@inline transport_workflow_module() = _transport_workflow_module()
@inline meson_workflow_module() = _meson_workflow_module()
@inline pnjl_module() = _pnjl_module()

@inline function workflow_param_adapters_module()
    bridge = IncludeOnce.include_once!(Main, :ModelsWorkflowParamAdaptersBridge, _WORKFLOW_PARAM_ADAPTERS_BRIDGE_PATH)
    adapters = Base.invokelatest(bridge.module_ref)
    if !isdefined(adapters, :normalize_quark_params)
        error("WorkflowParamAdapters module loaded but required API (normalize_quark_params) is missing")
    end
    return adapters
end
