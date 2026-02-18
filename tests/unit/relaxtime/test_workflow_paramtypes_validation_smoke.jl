using Test

include("../../../src/pnjl/workflows/TransportWorkflow.jl")
include("../../../src/pnjl/workflows/MesonMassWorkflow.jl")

using .TransportWorkflow
using .MesonMassWorkflow

@testset "Workflow ParamTypes validation smoke" begin
    T = 0.15
    mu = 0.0
    xi = 0.0

    base = TransportWorkflow.PNJL.solve(
        TransportWorkflow.PNJL.FixedMu(),
        T,
        mu;
        xi=xi,
        thermo_backend=:legacy,
        p_num=8,
        t_num=4,
        seed_strategy=TransportWorkflow.PNJL.DefaultSeed(phase_hint=:auto),
        iterations=30,
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
