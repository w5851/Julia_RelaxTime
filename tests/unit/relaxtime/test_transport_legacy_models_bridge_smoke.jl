using Test
using Printf

if !isdefined(Main, :TransportWorkflow)
    include("../../../src/pnjl/workflows/TransportWorkflow.jl")
end
using .TransportWorkflow

@testset "TransportWorkflow smoke: legacy/models transport bridge (4 fixed points)" begin
    BASELINE_PATH = joinpath(@__DIR__, "..", "..", "baselines", "relaxtime", "baseline_transport_fixedpoints_v1.csv")

    function _point_key(T::Float64, mu::Float64, xi::Float64)
        return @sprintf("%.6f|%.6f|%.6f", T, mu, xi)
    end

    function _load_baseline(path::String)
        isfile(path) || error("baseline CSV not found: $path")
        table = Dict{String,NamedTuple{(:eta,:sigma,:zeta),NTuple{3,Float64}}}()
        lines = readlines(path)
        isempty(lines) && error("baseline CSV is empty: $path")

        for line in lines[2:end]
            s = strip(line)
            isempty(s) && continue
            startswith(s, "#") && continue
            cols = split(s, ',')
            length(cols) == 6 || error("invalid baseline row: $line")
            T = parse(Float64, strip(cols[1]))
            mu = parse(Float64, strip(cols[2]))
            xi = parse(Float64, strip(cols[3]))
            eta = parse(Float64, strip(cols[4]))
            sigma = parse(Float64, strip(cols[5]))
            zeta = parse(Float64, strip(cols[6]))
            table[_point_key(T, mu, xi)] = (eta=eta, sigma=sigma, zeta=zeta)
        end
        return table
    end

    points = (
        (T=0.75, mu=0.00),
        (T=0.90, mu=0.00),
        (T=1.05, mu=0.00),
        (T=0.90, mu=0.15),
    )

    baseline = _load_baseline(BASELINE_PATH)

    tau = (u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0, sbar=1.0)

    rtol_transport = 8e-2
    atol_transport = 1e-6

    function _run_bridge_case(T, mu, xi, tau, models_solver, rtol_transport, atol_transport)
        common_kwargs = (
            xi=xi,
            thermo_backend=:models,
            tau=tau,
            compute_tau=false,
            compute_bulk=true,
            p_num=8,
            t_num=4,
            transport_config=TransportIntegrationConfig(p_nodes=8, p_max=3.5),
        )

        res_legacy = TransportWorkflow.solve_gap_and_transport(
            T,
            mu;
            common_kwargs...,
            solver_backend=:legacy,
            solver_kwargs=(iterations=30,),
        )

        res_models = TransportWorkflow.solve_gap_and_transport(
            T,
            mu;
            common_kwargs...,
            solver_backend=:models,
            models_solver=models_solver,
            models_residual_norm_max=1e-4,
            seed_state=collect(res_legacy.equilibrium.x_state),
        )

        return res_legacy, res_models
    end

    models_solver = Main.Models.NLsolveGapSolver(
        method=:trust_region,
        jacobian=:finite,
        xtol=1e-10,
        ftol=1e-10,
    )

    function _check_bridge_points(points, xi, tau, models_solver, rtol_transport, atol_transport)
        for pt in points
            T = pt.T
            mu = pt.mu
            point_label = "T=$(T), mu=$(mu), xi=$(xi)"

            @testset "bridge point $(point_label)" begin
                res_legacy, res_models = _run_bridge_case(T, mu, xi, tau, models_solver, rtol_transport, atol_transport)

                @testset "convergence $(point_label)" begin
                    @test res_legacy.equilibrium.converged isa Bool
                    @test res_models.equilibrium.converged isa Bool
                end

                @testset "finite masses $(point_label)" begin
                    @test all(isfinite, res_legacy.masses)
                    @test all(isfinite, res_models.masses)
                end

                @testset "finite transport $(point_label)" begin
                    @test isfinite(res_legacy.transport.eta)
                    @test isfinite(res_models.transport.eta)
                    @test isfinite(res_legacy.transport.sigma)
                    @test isfinite(res_models.transport.sigma)
                    @test isfinite(res_legacy.transport.zeta)
                    @test isfinite(res_models.transport.zeta)
                end

                @testset "approx eta $(point_label)" begin
                    @test isapprox(res_models.transport.eta, res_legacy.transport.eta; rtol=rtol_transport, atol=atol_transport)
                end

                @testset "approx sigma $(point_label)" begin
                    @test isapprox(res_models.transport.sigma, res_legacy.transport.sigma; rtol=rtol_transport, atol=atol_transport)
                end

                @testset "approx zeta $(point_label)" begin
                    @test isapprox(res_models.transport.zeta, res_legacy.transport.zeta; rtol=rtol_transport, atol=atol_transport)
                end

                baseline_key = _point_key(Float64(T), Float64(mu), Float64(xi))
                @test haskey(baseline, baseline_key)
                expected = baseline[baseline_key]

                @testset "legacy vs baseline $(point_label)" begin
                    @test isapprox(res_legacy.transport.eta, expected.eta; rtol=rtol_transport, atol=atol_transport)
                    @test isapprox(res_legacy.transport.sigma, expected.sigma; rtol=rtol_transport, atol=atol_transport)
                    @test isapprox(res_legacy.transport.zeta, expected.zeta; rtol=rtol_transport, atol=atol_transport)
                end

                @testset "models vs baseline $(point_label)" begin
                    @test isapprox(res_models.transport.eta, expected.eta; rtol=rtol_transport, atol=atol_transport)
                    @test isapprox(res_models.transport.sigma, expected.sigma; rtol=rtol_transport, atol=atol_transport)
                    @test isapprox(res_models.transport.zeta, expected.zeta; rtol=rtol_transport, atol=atol_transport)
                end
            end
        end
    end

    _check_bridge_points(points, 0.0, tau, models_solver, rtol_transport, atol_transport)

    include_xi_case = get(ENV, "UNIT_BRIDGE_XI_CASE", "0") in ("1", "true", "TRUE", "yes", "YES")
    if include_xi_case
        xi_points = ((T=0.90, mu=0.00),)
        _check_bridge_points(xi_points, 0.2, tau, models_solver, rtol_transport, atol_transport)
    end
end
