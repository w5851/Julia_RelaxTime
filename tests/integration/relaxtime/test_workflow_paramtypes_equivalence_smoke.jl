using Test

if !isdefined(Main, :Models)
    include("../../../src/models/Models.jl")
end

const TransportWorkflow = Main.Models.transport_workflow_module()
const MesonMassWorkflow = Main.Models.meson_workflow_module()

using .TransportWorkflow
using .MesonMassWorkflow

@testset "Workflow ParamTypes typed-input smoke" begin
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

    legacy_inputs = TransportWorkflow.as_legacy_inputs(inputs.quark_params, inputs.thermo_params)
    @test haskey(legacy_inputs, :quark_params)
    @test haskey(legacy_inputs, :thermo_params)

    A_struct = TransportWorkflow._A_from_equilibrium(
        T,
        inputs.quark_params,
        inputs.thermo_params;
        a_builder_config=(p_nodes=8, p_max=8.0, cos_nodes=4, use_aniso=false),
    )

    @test all(isfinite, A_struct)
    @test all(k -> isfinite(getfield(A_struct, k)), keys(A_struct))

    params = MesonMassWorkflow.build_equilibrium_params(base, T, mu; xi=xi)

    res_struct = MesonMassWorkflow._solve_meson_mass_with_retries(
        :pi,
        params.quark_params,
        params.thermo_params;
        k_norm=0.0,
        mass_kwargs=(iterations=20,),
    )

    @test res_struct !== nothing

    if res_struct !== nothing
        @test isfinite(res_struct.mass)
    end
end
