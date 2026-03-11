using Test

const _MODELS_PATH = abspath(joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl"))
if !isdefined(Main, :Models)
    Base.include(Main, _MODELS_PATH)
end
const TW = Main.Models.transport_workflow_module()

@testset "TransportWorkflow config smoke: prefer_energy_aniso via PHYSICS_PARAM_PROFILE" begin
    # This test validates toml-based defaults in the workflow layer.
    # It must NOT rely on module-load-time const caching, so it sets ENV at call time.

    withenv("PHYSICS_PARAM_PROFILE" => "unittest") do
        T = 0.15
        mu = 0.0
        xi = 0.2

        tau = (u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0, sbar=1.0)

        toy_provider_aniso = (
            energy_from_p=(p::Float64, m::Float64) -> sqrt(p * p + m * m),
            energy_from_p_aniso=(p::Float64, m::Float64, ξ::Float64, c::Float64) -> sqrt(p * p + m * m + ξ * (p * c)^2),
            quark_distribution=(E::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64) -> 0.20,
            antiquark_distribution=(E::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64) -> 0.20,
            quark_distribution_aniso=(p::Float64, m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, c::Float64) -> 0.90,
            antiquark_distribution_aniso=(p::Float64, m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, c::Float64) -> 0.90,
        )

        res_default = TW.solve_gap_and_transport(
            T,
            mu;
            xi=xi,
            tau=tau,
            compute_tau=false,
            compute_bulk=false,
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=30,),
            transport_config=TW.TransportIntegrationConfig(p_nodes=8, p_max=3.5, cos_nodes=6),
            provider=toy_provider_aniso,
        )

        res_true = TW.solve_gap_and_transport(
            T,
            mu;
            xi=xi,
            tau=tau,
            compute_tau=false,
            compute_bulk=false,
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=30,),
            transport_config=TW.TransportIntegrationConfig(p_nodes=8, p_max=3.5, cos_nodes=6),
            provider=toy_provider_aniso,
            prefer_energy_aniso=true,
        )

        res_false = TW.solve_gap_and_transport(
            T,
            mu;
            xi=xi,
            tau=tau,
            compute_tau=false,
            compute_bulk=false,
            p_num=8,
            t_num=4,
            solver_kwargs=(iterations=30,),
            transport_config=TW.TransportIntegrationConfig(p_nodes=8, p_max=3.5, cos_nodes=6),
            provider=toy_provider_aniso,
            prefer_energy_aniso=false,
        )

        @test isfinite(res_default.transport.eta)
        @test isfinite(res_true.transport.eta)
        @test isfinite(res_false.transport.eta)

        # unittest profile sets prefer_energy_aniso=false in config/physics/unittest.toml
        @test isapprox(res_default.transport.eta, res_false.transport.eta; rtol=1e-12, atol=0.0)
        @test !isapprox(res_default.transport.eta, res_true.transport.eta; rtol=1e-10, atol=0.0)

        @testset "explicit workflow_cache injection" begin
            cache = TW.WorkflowCache()
            res_custom_cache = TW.solve_gap_and_transport(
                T,
                mu;
                xi=xi,
                tau=tau,
                compute_tau=false,
                compute_bulk=false,
                p_num=8,
                t_num=4,
                solver_kwargs=(iterations=30,),
                transport_config=TW.TransportIntegrationConfig(p_nodes=8, p_max=3.5, cos_nodes=6),
                provider=toy_provider_aniso,
                workflow_cache=cache,
            )

            @test isfinite(res_custom_cache.transport.eta)
            @test haskey(cache.model_cache, :PNJL)
            @test haskey(cache.prefer_energy_aniso_cache, "unittest")

            _ = TW._A_from_equilibrium(
                res_custom_cache.thermo_params.T,
                res_custom_cache.quark_params,
                res_custom_cache.thermo_params;
                workflow_cache=cache,
            )
            @test haskey(cache.a_builder_config_cache, "unittest")

            TW.reset_transport_workflow_config_cache!(cache)
            @test haskey(cache.model_cache, :PNJL)
            @test isempty(cache.prefer_energy_aniso_cache)
            @test isempty(cache.a_builder_config_cache)
        end

        @testset "a_builder config from profile and explicit override" begin
            q = Main.ParameterTypes.QuarkParams((m=(u=1.2, d=1.2, s=1.8), μ=(u=0.1, d=0.1, s=0.2)))
            tp = Main.ParameterTypes.ThermoParams((T=0.15, Φ=0.45, Φbar=0.45, ξ=0.3))
            relaxtime_inputs = TW.as_relaxtime_inputs(q, tp)

            A_default = TW._A_from_equilibrium(tp.T, q, tp)
            A_expected_profile = Main.AFieldBuilder.build_A_triplet(
                relaxtime_inputs.quark_params,
                relaxtime_inputs.thermo_params;
                p_nodes=10,
                p_max=8.0,
                cos_nodes=6,
                use_aniso=true,
            )
            @test isapprox(A_default.u, A_expected_profile.u; rtol=1e-12, atol=0.0)
            @test isapprox(A_default.s, A_expected_profile.s; rtol=1e-12, atol=0.0)

            A_override = TW._A_from_equilibrium(
                tp.T,
                q,
                tp;
                a_builder_config=(p_nodes=12, p_max=7.0, cos_nodes=5, use_aniso=false),
            )
            A_expected_override = Main.AFieldBuilder.build_A_triplet(
                relaxtime_inputs.quark_params,
                relaxtime_inputs.thermo_params;
                p_nodes=12,
                p_max=7.0,
                cos_nodes=5,
                use_aniso=false,
            )
            @test isapprox(A_override.u, A_expected_override.u; rtol=1e-12, atol=0.0)
            @test isapprox(A_override.s, A_expected_override.s; rtol=1e-12, atol=0.0)
            @test !isapprox(A_override.u, A_default.u; rtol=1e-8, atol=0.0)
        end

        # Cache reset helper smoke:
        # - same Julia session
        # - switch PHYSICS_PARAM_PROFILE
        # - modify the toml on disk and verify reset forces re-read
        repo_root = abspath(joinpath(@__DIR__, "..", "..", ".."))
        unittest_toml = joinpath(repo_root, "config", "physics", "unittest.toml")

        @test get(ENV, "PHYSICS_PARAM_PROFILE", "<none>") == "unittest"

        original = read(unittest_toml, String)
        modified = replace(original, r"(?m)^prefer_energy_aniso\s*=\s*(true|false)\s*$" => "prefer_energy_aniso = true")
        @test modified != original

        try
            write(unittest_toml, modified)

            # Without resetting cache, default behavior should remain unchanged (still from cache).
            res_stale = TW.solve_gap_and_transport(
                T,
                mu;
                xi=xi,
                tau=tau,
                compute_tau=false,
                compute_bulk=false,
                p_num=8,
                t_num=4,
                solver_kwargs=(iterations=30,),
                transport_config=TW.TransportIntegrationConfig(p_nodes=8, p_max=3.5, cos_nodes=6),
                provider=toy_provider_aniso,
            )
            @test isapprox(res_stale.transport.eta, res_default.transport.eta; rtol=1e-12, atol=0.0)

            # After reset, it must re-read the file and reflect the new value.
            TW.reset_transport_workflow_config_cache!()
            res_refreshed = TW.solve_gap_and_transport(
                T,
                mu;
                xi=xi,
                tau=tau,
                compute_tau=false,
                compute_bulk=false,
                p_num=8,
                t_num=4,
                solver_kwargs=(iterations=30,),
                transport_config=TW.TransportIntegrationConfig(p_nodes=8, p_max=3.5, cos_nodes=6),
                provider=toy_provider_aniso,
            )
            @test isapprox(res_refreshed.transport.eta, res_true.transport.eta; rtol=1e-12, atol=0.0)
        finally
            write(unittest_toml, original)
            TW.reset_transport_workflow_config_cache!()
        end

        # Profile switch smoke (same Julia session): default profile sets prefer_energy_aniso=true.
        withenv("PHYSICS_PARAM_PROFILE" => "default") do
            TW.reset_transport_workflow_config_cache!()
            res_default_profile_default = TW.solve_gap_and_transport(
                T,
                mu;
                xi=xi,
                tau=tau,
                compute_tau=false,
                compute_bulk=false,
                p_num=8,
                t_num=4,
                solver_kwargs=(iterations=30,),
                transport_config=TW.TransportIntegrationConfig(p_nodes=8, p_max=3.5, cos_nodes=6),
                provider=toy_provider_aniso,
            )
            @test isapprox(res_default_profile_default.transport.eta, res_true.transport.eta; rtol=1e-12, atol=0.0)
        end
    end
end
