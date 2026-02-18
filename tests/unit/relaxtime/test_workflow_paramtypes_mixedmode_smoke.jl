using Test

include("../../../src/pnjl/workflows/TransportWorkflow.jl")
include("../../../src/pnjl/workflows/MesonMassWorkflow.jl")

using .TransportWorkflow
using .MesonMassWorkflow
using Main.ParameterTypes: as_namedtuple

@testset "Workflow ParamTypes mixed-mode smoke" begin
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

    inputs = TransportWorkflow._transport_inputs_from_equilibrium(base, T, mu;
        xi=xi,
        thermo_backend=:legacy,
        p_num=8,
        t_num=4,
    )

    A_qs_tn = TransportWorkflow._A_from_equilibrium(
        T,
        inputs.quark_params,
        as_namedtuple(inputs.thermo_params);
        a_builder_config=(p_nodes=8, p_max=8.0, cos_nodes=4, use_aniso=false),
    )

    A_qn_ts = TransportWorkflow._A_from_equilibrium(
        T,
        as_namedtuple(inputs.quark_params),
        inputs.thermo_params;
        a_builder_config=(p_nodes=8, p_max=8.0, cos_nodes=4, use_aniso=false),
    )

    @test all(isfinite, A_qs_tn)
    @test all(isfinite, A_qn_ts)
    @test all(k -> isapprox(getfield(A_qs_tn, k), getfield(A_qn_ts, k); rtol=1e-10, atol=1e-12), keys(A_qs_tn))

    params = MesonMassWorkflow.build_equilibrium_params(base, T, mu; xi=xi)

    res_qs_tn = MesonMassWorkflow._solve_meson_mass_with_retries(
        :pi,
        params.quark_params,
        as_namedtuple(params.thermo_params);
        k_norm=0.0,
        mass_kwargs=(iterations=20,),
    )

    res_qn_ts = MesonMassWorkflow._solve_meson_mass_with_retries(
        :pi,
        as_namedtuple(params.quark_params),
        params.thermo_params;
        k_norm=0.0,
        mass_kwargs=(iterations=20,),
    )

    @test res_qs_tn !== nothing
    @test res_qn_ts !== nothing

    if res_qs_tn !== nothing && res_qn_ts !== nothing
        @test isfinite(res_qs_tn.mass)
        @test isfinite(res_qn_ts.mass)
        @test isapprox(res_qs_tn.mass, res_qn_ts.mass; rtol=1e-7, atol=1e-9)
    end
end
