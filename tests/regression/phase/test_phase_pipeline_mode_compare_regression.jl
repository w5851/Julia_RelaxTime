using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "Phase pipeline production vs research regression" begin
    common = (
        T_grid=[120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0],
        rho_grid=collect(0.1:0.1:3.0),
        xi=0.0,
        profile=:regression,
        solver_backend=:models,
        reverse_rho=true,
        seed_policy=:hybrid_continuity,
        p_num=8,
        t_num=4,
        iterations=40,
        compute_crossover=true,
        crossover_method=:inflection,
        crossover_variable=:phi_u,
        crossover_n_mu=8,
        # This regression isolates production/research orchestration parity.
        # The production-only rho geometry gate has dedicated coverage.
        rho_geometry_convergence=false,
        promote_reference=false,
    )

    research = Models.run_phase_pipeline(
        :PNJL;
        mode=:research,
        output_dir=mktempdir(),
        common...,
    )

    production = Models.run_phase_pipeline(
        :PNJL;
        mode=:production,
        output_dir=mktempdir(),
        common...,
    )

    @test production.cep.result_status == research.cep.result_status
    @test production.cep.result_status == :ambiguous
    @test isapprox(production.cep.T_last_first_order_MeV, research.cep.T_last_first_order_MeV; atol=0.3)

    @test length(production.first_order_boundary) == length(research.first_order_boundary)
    for (prod_row, ref_row) in zip(production.first_order_boundary, research.first_order_boundary)
        @test prod_row.converged == ref_row.converged
        @test isapprox(prod_row.T_MeV, ref_row.T_MeV; atol=1e-10)
        @test isapprox(prod_row.mu_transition_MeV, ref_row.mu_transition_MeV; atol=1e-10)
        @test isapprox(prod_row.rho_hadron, ref_row.rho_hadron; atol=1e-10)
        @test isapprox(prod_row.rho_quark, ref_row.rho_quark; atol=1e-10)
    end

    @test length(production.crossover_line) == length(research.crossover_line)
    prod_conv = count(r -> r.converged, production.crossover_line)
    ref_conv = count(r -> r.converged, research.crossover_line)
    @test abs(prod_conv - ref_conv) <= 1

    for (prod_row, ref_row) in zip(production.crossover_line, research.crossover_line)
        @test isapprox(prod_row.mu_MeV, ref_row.mu_MeV; atol=0.3)
        if prod_row.converged && ref_row.converged
            @test isapprox(prod_row.T_crossover_MeV, ref_row.T_crossover_MeV; atol=2.0, nans=true)
            # Production starts the crossover search above the first-order boundary,
            # while research scans the complete requested T grid. The resulting
            # crossover temperatures are nearby rather than coordinate-identical;
            # compare the density scale relatively because rho(T, mu) is steep near
            # the endpoint instead of imposing an absolute same-coordinate contract.
            @test isapprox(prod_row.rho, ref_row.rho; rtol=1e-2, atol=1e-10, nans=true)
        end
    end
end
