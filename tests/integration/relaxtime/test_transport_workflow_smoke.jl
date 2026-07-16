using Test

if !isdefined(Main, :Models)
    include("../../../src/models/Models.jl")
end

if !isdefined(Main, :TransportWorkflow)
    const TransportWorkflow = Models.transport_workflow_module()
end
using .TransportWorkflow

@testset "TransportWorkflow smoke: solve_gap_and_transport (no tau/bulk)" begin
    T = 0.15
    mu = 0.0
    xi = 0.0

    # Provide tau so workflow doesn't spend time computing it.
    tau = (u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0, sbar=1.0)

    res = TransportWorkflow.solve_gap_and_transport(
        T,
        mu;
        xi=xi,
        tau=tau,
        compute_tau=false,
        compute_bulk=false,
        p_num=8,
        t_num=4,
        solver_kwargs=(iterations=30,),
        transport_config=TransportIntegrationConfig(p_nodes=8, p_max=3.5),
    )

    @testset "provider injection smoke" begin
        m = Models.create_model(:PNJL)
        prov = Models.transport_provider(m)
        prov_model = Models.transport_provider(m)

        qp_nt = Main.ParameterTypes.as_namedtuple(res.quark_params)
        tp_nt = Main.ParameterTypes.as_namedtuple(res.thermo_params)
        prov_prepared = Models.prepare_transport_provider(
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

        res_prov = TransportWorkflow.solve_gap_and_transport(
            T,
            mu;
            xi=xi,
            tau=tau,
            compute_tau=false,
            compute_bulk=false,
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=30,),
            transport_config=TransportIntegrationConfig(p_nodes=8, p_max=3.5),
            provider=prov,
        )

        res_model = TransportWorkflow.solve_gap_and_transport(
            T,
            mu;
            xi=xi,
            tau=tau,
            compute_tau=false,
            compute_bulk=false,
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=30,),
            transport_config=TransportIntegrationConfig(p_nodes=8, p_max=3.5),
            provider=prov_model,
        )

        res_prepared = TransportWorkflow.solve_gap_and_transport(
            T,
            mu;
            xi=xi,
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
        res_tuple = TransportWorkflow.solve_gap_and_transport(
            T,
            mu;
            xi=xi,
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

        res_toy = TransportWorkflow.solve_gap_and_transport(
            T,
            mu;
            xi=xi,
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

            res_pref_true = TransportWorkflow.solve_gap_and_transport(
                T,
                mu;
                xi=xi2,
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

            res_pref_false = TransportWorkflow.solve_gap_and_transport(
                T,
                mu;
                xi=xi2,
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
    @test isfinite(res.thermo_background.pressure)
    @test isfinite(res.thermo_background.entropy)
    @test isfinite(res.thermo_background.energy)
    @test isfinite(res.thermo_background.rho_mass)
    @test isnan(res.thermo_background.c_p)

    @test isfinite(res.transport.eta)
    @test isfinite(res.transport.sigma)
    @test res.transport.eta >= 0
    @test res.transport.sigma >= 0
    @test res.transport.lorenz_number isa Float64
    @test res.transport.lorentz_legacy isa Float64
    @test res.transport.viscous_conductive_coupling_ratio isa Float64
    @test isnan(res.transport.prandtl_number)

    @testset "post-equilibrium API: solve_transport_from_equilibrium" begin
        eq = TransportWorkflow.EquilibriumFacade.solve_equilibrium_backend(
            T,
            mu;
            xi=xi,
            solver_backend=:legacy,
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=30,),
        )

        res2 = TransportWorkflow.solve_transport_from_equilibrium(
            eq,
            T,
            mu;
            xi=xi,
            tau=tau,
            compute_tau=false,
            compute_bulk=false,
            p_num=8,
            t_num=4,
            transport_config=TransportIntegrationConfig(p_nodes=8, p_max=3.5),
        )

        prov = Models.transport_provider(Models.create_model(:PNJL))
        res2_prov = TransportWorkflow.solve_transport_from_equilibrium(
            eq,
            T,
            mu;
            xi=xi,
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
        @test isfinite(res2.thermo_background.pressure)
        @test isfinite(res2.thermo_background.entropy)
        @test isfinite(res2.thermo_background.energy)
        @test isfinite(res2.thermo_background.rho_mass)
        @test isnan(res2.thermo_background.c_p)
        @test isfinite(res2.transport.eta)
        @test isfinite(res2.transport.sigma)
        @test res2.transport.eta >= 0
        @test res2.transport.sigma >= 0
        @test res2.transport.viscous_conductive_coupling_ratio isa Float64
        @test isnan(res2.transport.prandtl_number)

        @test isapprox(res2_prov.transport.eta, res2.transport.eta; rtol=1e-12, atol=0.0)
        @test isapprox(res2_prov.transport.sigma, res2.transport.sigma; rtol=1e-12, atol=0.0)

        res2_bulk = TransportWorkflow.solve_transport_from_equilibrium(
            eq,
            T,
            mu;
            xi=xi,
            tau=tau,
            compute_tau=false,
            compute_bulk=true,
            p_num=8,
            t_num=4,
            transport_config=TransportIntegrationConfig(p_nodes=8, p_max=3.5),
        )
        @test res2_bulk.bulk_coeffs !== nothing
        @test res2_bulk.bulk_coeffs.base_state_source === :workflow_equilibrium
        @test res2_bulk.bulk_coeffs.primal_solve_count == 0
        @test res2_bulk.bulk_coeffs.jacobian_factorization_count == 1
        @test res2_bulk.bulk_coeffs.derivative_series_count == 3
        @test res2_bulk.bulk_coeffs.branch_locked
        @test all(isapprox.(res2_bulk.bulk_coeffs.masses, eq.masses; rtol=1e-4, atol=1e-8))
    end

    @testset "optional Pr background passthrough" begin
        res_pr = TransportWorkflow.solve_gap_and_transport(
            T,
            mu;
            xi=xi,
            tau=tau,
            compute_tau=false,
            compute_bulk=false,
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=30,),
            transport_config=TransportIntegrationConfig(p_nodes=8, p_max=3.5),
            c_p=1.8,
        )

        @test isfinite(res_pr.thermo_background.rho_mass)
        @test isapprox(res_pr.thermo_background.c_p, 1.8; rtol=1e-12, atol=0.0)
        @test res_pr.transport.prandtl_number isa Float64
    end

    @testset "anisotropic diffusion and derived transport chain" begin
        T_chain = 0.9
        mu_chain = 0.15
        xi_chain = 0.2
        cfg_chain = TransportIntegrationConfig(p_nodes=8, p_max=3.5, cos_nodes=6)
        res_chain = TransportWorkflow.solve_gap_and_transport(
            T_chain,
            mu_chain;
            xi=xi_chain,
            tau=tau,
            compute_tau=false,
            compute_bulk=false,
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=30,),
            transport_config=cfg_chain,
            c_p=1.8,
        )

        @test res_chain.tau == tau
        @test res_chain.tau_inv === nothing
        @test res_chain.rates === nothing
        @test isfinite(res_chain.transport.kappa_BB)
        @test isfinite(res_chain.transport.kappa_BQ)
        @test isfinite(res_chain.transport.kappa_BS)
        @test isfinite(res_chain.transport.kappa_QQ)
        @test isfinite(res_chain.transport.kappa_QS)
        @test isfinite(res_chain.transport.kappa_SS)
        @test isfinite(res_chain.transport.lambda)
        @test isfinite(res_chain.transport.lorenz_number)
        @test isfinite(res_chain.transport.prandtl_number)

        n_B = TransportWorkflow.TransportCoefficients.conserved_charge_densities(res_chain.densities).B
        enthalpy = res_chain.thermo_background.energy + res_chain.thermo_background.pressure
        expected_lambda = res_chain.transport.kappa_BB * (enthalpy / (n_B * T_chain))^2
        @test isapprox(res_chain.transport.lambda, expected_lambda; rtol=1e-12, atol=0.0)
        @test isapprox(
            res_chain.transport.lorenz_number,
            res_chain.transport.lambda / (res_chain.transport.sigma * T_chain);
            rtol=1e-12,
            atol=0.0,
        )
        @test isapprox(
            res_chain.transport.prandtl_number,
            res_chain.transport.eta * 1.8 / (res_chain.transport.lambda * res_chain.thermo_background.rho_mass);
            rtol=1e-12,
            atol=0.0,
        )

        provider_dist5 = (
            energy_from_p=(p::Float64, m::Float64) -> 2.0,
            energy_from_p_aniso=(p::Float64, m::Float64, ξ::Float64, c::Float64) -> 5.0,
            quark_distribution=(E::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64) -> 0.2,
            antiquark_distribution=(E::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64) -> 0.2,
            quark_distribution_aniso=(p::Float64, m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, c::Float64) -> 0.2,
            antiquark_distribution_aniso=(p::Float64, m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, c::Float64) -> 0.2,
            prefer_energy_aniso=true,
        )
        provider_dist7 = merge(provider_dist5, (
            energy_from_p_aniso=(p::Float64, m::Float64, ξ::Float64, c::Float64) -> 7.0,
        ))

        function solve_with_provider(provider)
            return TransportWorkflow.solve_transport_from_equilibrium(
                res_chain.equilibrium,
                T_chain,
                mu_chain;
                xi=xi_chain,
                tau=tau,
                compute_tau=false,
                compute_bulk=false,
                p_num=8,
                t_num=4,
                transport_config=cfg_chain,
                provider=provider,
                c_p=1.8,
            )
        end

        res_dist5 = solve_with_provider(provider_dist5)
        res_dist7 = solve_with_provider(provider_dist7)
        for field in (:eta, :sigma, :kappa_BB, :kappa_BQ, :kappa_BS, :kappa_QQ, :kappa_QS, :kappa_SS, :lambda, :lorenz_number, :prandtl_number)
            @test isapprox(
                getproperty(res_dist7.transport, field),
                getproperty(res_dist5.transport, field);
                rtol=1e-12,
                atol=0.0,
            )
        end
    end
end
