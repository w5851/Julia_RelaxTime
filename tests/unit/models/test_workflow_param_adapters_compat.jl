using Test

const PROJECT_ROOT_WPA_COMPAT = normpath(joinpath(@__DIR__, "..", "..", ".."))

const _PT_PATH_WPA_COMPAT = normpath(joinpath(PROJECT_ROOT_WPA_COMPAT, "src", "ParameterTypes.jl"))
const _MODELS_PATH_WPA_COMPAT = normpath(joinpath(PROJECT_ROOT_WPA_COMPAT, "src", "models", "Models.jl"))
if !isdefined(Main, :ParameterTypes)
    Base.include(Main, _PT_PATH_WPA_COMPAT)
end
if !isdefined(Main, :Models)
    Base.include(Main, _MODELS_PATH_WPA_COMPAT)
end

const WorkflowParamAdaptersCompatModule = Main.Models.workflow_param_adapters_module()

using .WorkflowParamAdaptersCompatModule: as_relaxtime_inputs
using Main.ParameterTypes: QuarkParams, ThermoParams

@testset "WorkflowParamAdapters compat aliases" begin
    qp = QuarkParams((m=(u=0.3, d=0.31, s=0.5), μ=(u=0.0, d=0.0, s=0.0)))
    tp = ThermoParams((T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0))
    @test_deprecated WorkflowParamAdaptersCompatModule.as_legacy_inputs(qp, tp)
    @test as_relaxtime_inputs(qp, tp) == WorkflowParamAdaptersCompatModule.as_legacy_inputs(qp, tp)
end