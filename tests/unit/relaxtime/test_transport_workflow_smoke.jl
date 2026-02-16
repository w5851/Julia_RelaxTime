using Test

include("../../../src/pnjl/workflows/TransportWorkflow.jl")
using .TransportWorkflow

@testset "TransportWorkflow smoke: solve_gap_and_transport (no tau/bulk)" begin
    T = 0.15
    mu = 0.0
    xi = 0.0

    # Provide tau so workflow doesn't spend time computing it.
    tau = (u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0, sbar=1.0)

    res = solve_gap_and_transport(
        T,
        mu;
        xi=xi,
        thermo_backend=:models,
        tau=tau,
        compute_tau=false,
        compute_bulk=false,
        p_num=8,
        t_num=4,
        solver_kwargs=(iterations=30,),
        transport_config=TransportIntegrationConfig(p_nodes=8, p_max=3.5),
    )

    @testset "provider injection smoke" begin
        prov = Main.Models.transport_provider(:models)

        kind = TransportWorkflow.EquilibriumFacade.pnjl_model_kind(:models)
        m = TransportWorkflow.ThermoFacade.get_models_model(kind)
        prov_model = Main.Models.transport_provider(m)

        qp_nt = Main.ParameterTypes.as_namedtuple(res.quark_params)
        tp_nt = Main.ParameterTypes.as_namedtuple(res.thermo_params)
        prov_prepared = Main.Models.prepare_transport_provider(
            prov,
            res.equilibrium;
            quark_params=qp_nt,
            thermo_params=tp_nt,
            masses=res.masses,
        )
        @test hasproperty(prov_prepared, :ctx)
        @test hasproperty(prov_prepared.ctx, :masses)
        qp_wrong = (m=(u=999.0, d=999.0, s=999.0), μ=qp_nt.μ)
        @test prov_prepared.mass_for_species(:u, qp_wrong, tp_nt) == prov_prepared.ctx.masses.u

        res_prov = solve_gap_and_transport(
            T,
            mu;
            xi=xi,
            thermo_backend=:models,
            tau=tau,
            compute_tau=false,
            compute_bulk=false,
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=30,),
            transport_config=TransportIntegrationConfig(p_nodes=8, p_max=3.5),
            provider=prov,
        )

        res_model = solve_gap_and_transport(
            T,
            mu;
            xi=xi,
            thermo_backend=:models,
            tau=tau,
            compute_tau=false,
            compute_bulk=false,
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=30,),
            transport_config=TransportIntegrationConfig(p_nodes=8, p_max=3.5),
            provider=prov_model,
        )

        res_prepared = solve_gap_and_transport(
            T,
            mu;
            xi=xi,
            thermo_backend=:models,
            tau=tau,
            compute_tau=false,
            compute_bulk=false,
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=30,),
            transport_config=TransportIntegrationConfig(p_nodes=8, p_max=3.5),
            provider=prov_prepared,
        )

        @test isapprox(res_prov.transport.eta, res.transport.eta; rtol=1e-12, atol=0.0)
        @test isapprox(res_prov.transport.sigma, res.transport.sigma; rtol=1e-12, atol=0.0)
        @test isapprox(res_model.transport.eta, res.transport.eta; rtol=1e-12, atol=0.0)
        @test isapprox(res_model.transport.sigma, res.transport.sigma; rtol=1e-12, atol=0.0)
        @test isapprox(res_model.transport.eta, res_prov.transport.eta; rtol=1e-12, atol=0.0)
        @test isapprox(res_model.transport.sigma, res_prov.transport.sigma; rtol=1e-12, atol=0.0)
        @test isapprox(res_prepared.transport.eta, res.transport.eta; rtol=1e-12, atol=0.0)
        @test isapprox(res_prepared.transport.sigma, res.transport.sigma; rtol=1e-12, atol=0.0)

        prov_tuple = TransportWorkflow.TransportCoefficients.default_transport_provider()
        res_tuple = solve_gap_and_transport(
            T,
            mu;
            xi=xi,
            thermo_backend=:models,
            tau=tau,
            compute_tau=false,
            compute_bulk=false,
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=30,),
            transport_config=TransportIntegrationConfig(p_nodes=8, p_max=3.5),
            provider=prov_tuple,
        )

        @test isapprox(res_tuple.transport.eta, res.transport.eta; rtol=1e-12, atol=0.0)
        @test isapprox(res_tuple.transport.sigma, res.transport.sigma; rtol=1e-12, atol=0.0)

        toy_provider = (
            energy_from_p=(p::Float64, m::Float64) -> sqrt(p * p + m * m),
            quark_distribution=(E::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64) -> 0.1,
            antiquark_distribution=(E::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64) -> 0.1,
            quark_distribution_aniso=(p::Float64, m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, c::Float64) -> 0.1,
            antiquark_distribution_aniso=(p::Float64, m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, c::Float64) -> 0.1,
        )

        res_toy = solve_gap_and_transport(
            T,
            mu;
            xi=xi,
            thermo_backend=:models,
            tau=tau,
            compute_tau=false,
            compute_bulk=false,
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=30,),
            transport_config=TransportIntegrationConfig(p_nodes=8, p_max=3.5),
            provider=toy_provider,
        )

        @test isfinite(res_toy.transport.eta)
        @test isfinite(res_toy.transport.sigma)
        @test !isapprox(res_toy.transport.eta, res.transport.eta; rtol=1e-6, atol=0.0)

        @testset "prefer_energy_aniso is plumbed from workflow" begin
            xi2 = 0.2

            toy_provider_aniso = (
                energy_from_p=(p::Float64, m::Float64) -> sqrt(p * p + m * m),
                energy_from_p_aniso=(p::Float64, m::Float64, ξ::Float64, c::Float64) -> sqrt(p * p + m * m + ξ * (p * c)^2),
                quark_distribution=(E::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64) -> 0.20,
                antiquark_distribution=(E::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64) -> 0.20,
                quark_distribution_aniso=(p::Float64, m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, c::Float64) -> 0.90,
                antiquark_distribution_aniso=(p::Float64, m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, c::Float64) -> 0.90,
            )

            res_pref_true = solve_gap_and_transport(
                T,
                mu;
                xi=xi2,
                thermo_backend=:models,
                tau=tau,
                compute_tau=false,
                compute_bulk=false,
                p_num=8,
                t_num=4,
                solver_kwargs=(iterations=30,),
                transport_config=TransportIntegrationConfig(p_nodes=8, p_max=3.5, cos_nodes=6),
                provider=toy_provider_aniso,
                transport_kwargs=(; prefer_energy_aniso=true),
            )

            res_pref_false = solve_gap_and_transport(
                T,
                mu;
                xi=xi2,
                thermo_backend=:models,
                tau=tau,
                compute_tau=false,
                compute_bulk=false,
                p_num=8,
                t_num=4,
                solver_kwargs=(iterations=30,),
                transport_config=TransportIntegrationConfig(p_nodes=8, p_max=3.5, cos_nodes=6),
                provider=toy_provider_aniso,
                transport_kwargs=(; prefer_energy_aniso=false),
            )

            @test isfinite(res_pref_true.transport.eta)
            @test isfinite(res_pref_false.transport.eta)
            @test !isapprox(res_pref_true.transport.eta, res_pref_false.transport.eta; rtol=1e-10, atol=0.0)
        end
    end

    @test haskey(res, :equilibrium)
    @test res.equilibrium.converged isa Bool

    @test length(res.masses) == 3
    @test all(isfinite, res.masses)

    @test isfinite(res.transport.eta)
    @test isfinite(res.transport.sigma)
    @test res.transport.eta >= 0
    @test res.transport.sigma >= 0

    @testset "post-equilibrium API: solve_transport_from_equilibrium" begin
        eq = TransportWorkflow.EquilibriumFacade.solve_equilibrium_backend(
            T,
            mu;
            xi=xi,
            thermo_backend=:models,
            solver_backend=:legacy,
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=30,),
        )

        res2 = solve_transport_from_equilibrium(
            eq,
            T,
            mu;
            xi=xi,
            thermo_backend=:models,
            tau=tau,
            compute_tau=false,
            compute_bulk=false,
            p_num=8,
            t_num=4,
            transport_config=TransportIntegrationConfig(p_nodes=8, p_max=3.5),
        )

        prov = Main.Models.transport_provider(:models)
        res2_prov = solve_transport_from_equilibrium(
            eq,
            T,
            mu;
            xi=xi,
            thermo_backend=:models,
            tau=tau,
            compute_tau=false,
            compute_bulk=false,
            p_num=8,
            t_num=4,
            transport_config=TransportIntegrationConfig(p_nodes=8, p_max=3.5),
            provider=prov,
        )

        @test res2.equilibrium.converged isa Bool
        @test all(isfinite, res2.masses)
        @test isfinite(res2.transport.eta)
        @test isfinite(res2.transport.sigma)
        @test res2.transport.eta >= 0
        @test res2.transport.sigma >= 0

        @test isapprox(res2_prov.transport.eta, res2.transport.eta; rtol=1e-12, atol=0.0)
        @test isapprox(res2_prov.transport.sigma, res2.transport.sigma; rtol=1e-12, atol=0.0)
    end
end
