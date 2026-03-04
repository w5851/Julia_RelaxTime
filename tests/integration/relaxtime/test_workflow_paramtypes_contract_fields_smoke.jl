using Test

if !isdefined(Main, :Models)
    if !isdefined(Main, :Models)
    include("../../../src/models/Models.jl")
end
end
Main.Models.workflow_param_adapters_module()
using .WorkflowParamAdapters

@testset "Workflow ParamTypes contract fields smoke" begin
    @test_throws ArgumentError WorkflowParamAdapters.normalize_quark_params((
        μ=(u=0.0, d=0.0, s=0.0),
    ))

    @test_throws ArgumentError WorkflowParamAdapters.normalize_quark_params((
        m=(u=300.0, d=300.0),
        μ=(u=0.0, d=0.0, s=0.0),
    ))

    @test_throws ArgumentError WorkflowParamAdapters.normalize_quark_params((
        m=(u=300.0, d=300.0, s=500.0),
        μ=(u=0.0, d=0.0, s="bad"),
    ))

    @test_throws ArgumentError WorkflowParamAdapters.normalize_thermo_params((
        T=0.15,
        Φ=0.1,
        ξ=0.0,
    ))

    @test_throws ArgumentError WorkflowParamAdapters.normalize_thermo_params((
        T=0.15,
        Φ=0.1,
        Φbar=0.1,
        ξ=Inf,
    ))
end
