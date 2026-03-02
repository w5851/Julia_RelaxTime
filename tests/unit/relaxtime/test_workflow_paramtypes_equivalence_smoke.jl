using Test

if !isdefined(Main, :Models)
    if !isdefined(Main, :Models)
    include("../../../src/models/Models.jl")
end
end
Main.Models.transport_workflow_module()
Main.Models.meson_workflow_module()

using .TransportWorkflow
using .MesonMassWorkflow
using Main.ParameterTypes: as_namedtuple

@testset "Workflow ParamTypes equivalence smoke" begin
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

    inputs = TransportWorkflow._transport_inputs_from_equilibrium(base, T, mu;
        xi=xi,
        p_num=8,
        t_num=4,
    )

    legacy_from_struct = TransportWorkflow.as_legacy_inputs(inputs.quark_params, inputs.thermo_params)
    legacy_from_nt = TransportWorkflow.as_legacy_inputs(
        as_namedtuple(inputs.quark_params),
        as_namedtuple(inputs.thermo_params),
    )

    @test legacy_from_struct.quark_params == legacy_from_nt.quark_params
    @test legacy_from_struct.thermo_params == legacy_from_nt.thermo_params

    A_struct = TransportWorkflow._A_from_equilibrium(
        T,
        inputs.quark_params,
        inputs.thermo_params;
        a_builder_config=(p_nodes=8, p_max=8.0, cos_nodes=4, use_aniso=false),
    )

    A_nt = TransportWorkflow._A_from_equilibrium(
        T,
        as_namedtuple(inputs.quark_params),
        as_namedtuple(inputs.thermo_params);
        a_builder_config=(p_nodes=8, p_max=8.0, cos_nodes=4, use_aniso=false),
    )

    @test all(isfinite, A_struct)
    @test all(isfinite, A_nt)
    @test all(k -> isapprox(getfield(A_struct, k), getfield(A_nt, k); rtol=1e-10, atol=1e-12), keys(A_struct))

    params = MesonMassWorkflow.build_equilibrium_params(base, T, mu; xi=xi)

    res_struct = MesonMassWorkflow._solve_meson_mass_with_retries(
        :pi,
        params.quark_params,
        params.thermo_params;
        k_norm=0.0,
        mass_kwargs=(iterations=20,),
    )

    res_nt = MesonMassWorkflow._solve_meson_mass_with_retries(
        :pi,
        as_namedtuple(params.quark_params),
        as_namedtuple(params.thermo_params);
        k_norm=0.0,
        mass_kwargs=(iterations=20,),
    )

    @test res_struct !== nothing
    @test res_nt !== nothing

    if res_struct !== nothing && res_nt !== nothing
        @test isfinite(res_struct.mass)
        @test isfinite(res_nt.mass)
        @test isapprox(res_struct.mass, res_nt.mass; rtol=1e-7, atol=1e-9)
    end
end
