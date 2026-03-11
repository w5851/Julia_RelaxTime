# WorkflowParamAdapters 模块单元测试

using Test

const PROJECT_ROOT_WPA = normpath(joinpath(@__DIR__, "..", "..", ".."))

const _PT_PATH_WPA = normpath(joinpath(PROJECT_ROOT_WPA, "src", "ParameterTypes.jl"))
const _MODELS_PATH_WPA = normpath(joinpath(PROJECT_ROOT_WPA, "src", "models", "Models.jl"))
if !isdefined(Main, :ParameterTypes)
    Base.include(Main, _PT_PATH_WPA)
end
if !isdefined(Main, :Models)
    Base.include(Main, _MODELS_PATH_WPA)
end

const WorkflowParamAdaptersModule = Main.Models.workflow_param_adapters_module()

using .WorkflowParamAdaptersModule: normalize_quark_params, normalize_thermo_params, as_relaxtime_inputs
using Main.ParameterTypes: QuarkParams, ThermoParams

@testset "WorkflowParamAdapters" begin
    @testset "normalize_quark_params QuarkParams 输入" begin
        qp = QuarkParams((m=(u=0.3, d=0.31, s=0.5), μ=(u=0.0, d=0.0, s=0.0)))
        result = normalize_quark_params(qp)
        @test result.m.u ≈ 0.3
        @test result.m.s ≈ 0.5
    end

    @testset "normalize_thermo_params ThermoParams 输入" begin
        tp = ThermoParams((T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0))
        result = normalize_thermo_params(tp)
        @test result.T ≈ 0.15
        @test result.ξ ≈ 0.0
    end

    @testset "as_relaxtime_inputs 返回正确结构" begin
        qp = QuarkParams((m=(u=0.3, d=0.31, s=0.5), μ=(u=0.0, d=0.0, s=0.0)))
        tp = ThermoParams((T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0))
        relaxtime_inputs = as_relaxtime_inputs(qp, tp)
        @test haskey(relaxtime_inputs, :quark_params)
        @test haskey(relaxtime_inputs, :thermo_params)
    end

    @testset "NamedTuple 输入已拒绝" begin
        @test_throws ArgumentError normalize_quark_params((m=(u=0.3, d=0.31, s=0.5), μ=(u=0.0, d=0.0, s=0.0)))
        @test_throws ArgumentError normalize_thermo_params((T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0))
    end
end
