# WorkflowParamAdapters 模块单元测试

using Test

const PROJECT_ROOT_WPA = normpath(joinpath(@__DIR__, "..", "..", ".."))

const _PT_PATH_WPA = normpath(joinpath(PROJECT_ROOT_WPA, "src", "ParameterTypes.jl"))
if !isdefined(Main, :ParameterTypes)
    Base.include(Main, _PT_PATH_WPA)
end

const _WPA_PATH = normpath(joinpath(PROJECT_ROOT_WPA, "src", "models", "workflows", "WorkflowParamAdapters.jl"))
if !isdefined(Main, :WorkflowParamAdapters)
    Base.include(Main, _WPA_PATH)
end

using Main.WorkflowParamAdapters: normalize_quark_params, normalize_thermo_params, as_legacy_inputs
using Main.ParameterTypes: QuarkParams, ThermoParams

@testset "WorkflowParamAdapters" begin
    @testset "normalize_quark_params NamedTuple 输入" begin
        qp_nt = (m=(u=0.3, d=0.31, s=0.5), μ=(u=0.0, d=0.0, s=0.0))
        result = normalize_quark_params(qp_nt)
        @test result.m.u ≈ 0.3
        @test result.m.s ≈ 0.5
    end

    @testset "normalize_thermo_params NamedTuple 输入" begin
        tp_nt = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
        result = normalize_thermo_params(tp_nt)
        @test result.T ≈ 0.15
        @test result.ξ ≈ 0.0
    end

    @testset "as_legacy_inputs 返回正确结构" begin
        qp = (m=(u=0.3, d=0.31, s=0.5), μ=(u=0.0, d=0.0, s=0.0))
        tp = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
        legacy = as_legacy_inputs(qp, tp)
        @test haskey(legacy, :quark_params)
        @test haskey(legacy, :thermo_params)
    end
end
