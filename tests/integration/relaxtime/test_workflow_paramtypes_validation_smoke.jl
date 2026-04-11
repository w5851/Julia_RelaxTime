using Test

if !isdefined(Main, :Models)
    include("../../../src/models/Models.jl")
end

if !isdefined(Main, :TransportWorkflow)
    const TransportWorkflow = Main.Models.transport_workflow_module()
end
if !isdefined(Main, :MesonMassWorkflow)
    const MesonMassWorkflow = Main.Models.meson_workflow_module()
end

using .TransportWorkflow
using .MesonMassWorkflow

@testset "Workflow ParamTypes validation smoke" begin
    T = 0.15
    mu = 0.0
    xi = 0.0

    base = TransportWorkflow.EquilibriumFacade.solve_equilibrium_backend(
        T,
        mu;
        xi=xi,
        solver_backend=:models,
        p_num=8,
        t_num=4,
        seed_state=TransportWorkflow.PNJL.HADRON_SEED_5,
        models_residual_norm_max=1e-4,
    )

    params = MesonMassWorkflow.build_equilibrium_params(base, T, mu; xi=xi)

    @test_throws ArgumentError TransportWorkflow._A_from_equilibrium(
        T,
        123,
        params.thermo_params;
        a_builder_config=(p_nodes=8, p_max=8.0, cos_nodes=4, use_aniso=false),
    )

    @test_throws ArgumentError TransportWorkflow._A_from_equilibrium(
        T,
        params.quark_params,
        456;
        a_builder_config=(p_nodes=8, p_max=8.0, cos_nodes=4, use_aniso=false),
    )

    @test_throws ArgumentError MesonMassWorkflow._solve_meson_mass_with_retries(
        :pi,
        123,
        params.thermo_params;
        k_norm=0.0,
        mass_kwargs=(iterations=20,),
    )

    @test_throws ArgumentError MesonMassWorkflow._solve_meson_mass_with_retries(
        :pi,
        params.quark_params,
        456;
        k_norm=0.0,
        mass_kwargs=(iterations=20,),
    )
end
