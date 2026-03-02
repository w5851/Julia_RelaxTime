using Test

if !isdefined(Main, :Models)
    if !isdefined(Main, :Models)
    include("../../../src/models/Models.jl")
end
end
Main.Models.workflow_param_adapters_module()
using .WorkflowParamAdapters

@testset "Workflow ParamTypes depwarn smoke" begin
    WorkflowParamAdapters._QUARK_NAMEDTUPLE_DEPWARN_EMITTED[] = false
    WorkflowParamAdapters._THERMO_NAMEDTUPLE_DEPWARN_EMITTED[] = false

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

    @test_logs (:warn, r"Passing NamedTuple to quark_params is deprecated") WorkflowParamAdapters.normalize_quark_params(q_nt)
    @test_logs (:warn, r"Passing NamedTuple to thermo_params is deprecated") WorkflowParamAdapters.normalize_thermo_params(t_nt)
end
