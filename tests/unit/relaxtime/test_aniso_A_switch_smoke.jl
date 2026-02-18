using Test

if !isdefined(Main, :RelaxationTime)
    include("../../../src/relaxtime/RelaxationTime.jl")
end
if !isdefined(Main, :TransportWorkflow)
    include("../../../src/pnjl/workflows/TransportWorkflow.jl")
end

using .RelaxationTime
using .TransportWorkflow
using Main.ParameterTypes: QuarkParams, ThermoParams
using Main.OneLoopIntegrals: A
using Main.OneLoopIntegralsCorrection: A_aniso

@testset "Anisotropic A auto-switch smoke" begin
    quark_params = (
        m=(u=1.20, d=1.20, s=1.80),
        μ=(u=0.10, d=0.10, s=0.20),
    )

    thermo_iso = (T=0.15, Φ=0.45, Φbar=0.45, ξ=0.0)
    thermo_aniso = (T=0.15, Φ=0.45, Φbar=0.45, ξ=0.30)

    @testset "RelaxationTime.ensure_quark_params_has_A switches to A_aniso" begin
        q_iso = RelaxationTime.ensure_quark_params_has_A(quark_params, thermo_iso)
        q_aniso = RelaxationTime.ensure_quark_params_has_A(quark_params, thermo_aniso)

        @test hasproperty(q_iso, :A)
        @test hasproperty(q_aniso, :A)

        nodes_p, weights_p = RelaxationTime.AverageScatteringRate.gauleg(0.0, 20.0, 16)
        nodes_cos, weights_cos = RelaxationTime.AverageScatteringRate.gauleg(-1.0, 1.0, RelaxationTime.DEFAULT_ANGLE_NODES)

        A_u_iso_expected = A(quark_params.m.u, quark_params.μ.u, thermo_iso.T, thermo_iso.Φ, thermo_iso.Φbar, nodes_p, weights_p)
        A_u_aniso_expected = A_aniso(quark_params.m.u, quark_params.μ.u, thermo_aniso.T, thermo_aniso.Φ, thermo_aniso.Φbar,
                                     thermo_aniso.ξ, nodes_p, weights_p, nodes_cos, weights_cos)

        @test isapprox(q_iso.A.u, A_u_iso_expected; rtol=1e-11, atol=0.0)
        @test isapprox(q_aniso.A.u, A_u_aniso_expected; rtol=1e-11, atol=0.0)
        @test !isapprox(q_aniso.A.u, q_iso.A.u; rtol=1e-8, atol=0.0)
    end

    @testset "TransportWorkflow._A_from_equilibrium switches to A_aniso" begin
        q = QuarkParams((m=quark_params.m, μ=quark_params.μ))
        thermo_iso_struct = ThermoParams(thermo_iso)
        thermo_aniso_struct = ThermoParams(thermo_aniso)

        A_map_iso = TransportWorkflow._A_from_equilibrium(thermo_iso.T, q, thermo_iso_struct)
        A_map_aniso = TransportWorkflow._A_from_equilibrium(thermo_aniso.T, q, thermo_aniso_struct)

        nodes = TransportWorkflow.DEFAULT_MOMENTUM_NODES
        weights = TransportWorkflow.DEFAULT_MOMENTUM_WEIGHTS
        nodes_cos, weights_cos = RelaxationTime.AverageScatteringRate.gauleg(-1.0, 1.0, 4)

        A_u_iso_expected = A(quark_params.m.u, quark_params.μ.u, thermo_iso.T, thermo_iso.Φ, thermo_iso.Φbar, nodes, weights)
        A_u_aniso_expected = A_aniso(quark_params.m.u, quark_params.μ.u, thermo_aniso.T, thermo_aniso.Φ, thermo_aniso.Φbar,
                                     thermo_aniso.ξ, nodes, weights, nodes_cos, weights_cos)

        @test isapprox(A_map_iso.u, A_u_iso_expected; rtol=5e-7, atol=0.0)
        @test isapprox(A_map_aniso.u, A_u_aniso_expected; rtol=5e-7, atol=0.0)
        @test !isapprox(A_map_aniso.u, A_map_iso.u; rtol=1e-8, atol=0.0)
    end
end
