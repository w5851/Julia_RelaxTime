using Test

if !isdefined(Main, :Models)
    include("../../../src/models/Models.jl")
end
Main.Models.workflow_param_adapters_module()
using .WorkflowParamAdapters

@testset "Workflow ParamTypes depwarn smoke" begin
    q_nt = (
        m=(u=300.0, d=300.0, s=500.0),
        μ=(u=0.0, d=0.0, s=0.0),
    )
    t_nt = (
        T=0.15,
        Φ=0.1,
        Φbar=0.1,
        ξ=0.0,
    )

    @test_logs (:warn, r"quark_params is deprecated") WorkflowParamAdapters.normalize_quark_params(q_nt)
    @test_logs (:warn, r"thermo_params is deprecated") WorkflowParamAdapters.normalize_thermo_params(t_nt)
end
