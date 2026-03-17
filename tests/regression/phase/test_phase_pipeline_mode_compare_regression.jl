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
        solver_backend=:legacy,
        reverse_rho=true,
        seed_policy=:hybrid_continuity,
        p_num=8,
        t_num=4,
        iterations=40,
        compute_crossover=true,
        crossover_method=:inflection,
        crossover_variable=:phi_u,
        crossover_n_mu=8,
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

    @test production.cep.found == research.cep.found
    @test isapprox(production.cep.T_cep_MeV, research.cep.T_cep_MeV; atol=0.3)
    @test isapprox(production.cep.mu_cep_MeV, research.cep.mu_cep_MeV; atol=0.3)

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
            @test isapprox(prod_row.rho, ref_row.rho; atol=0.01, nans=true)
        end
    end
end
