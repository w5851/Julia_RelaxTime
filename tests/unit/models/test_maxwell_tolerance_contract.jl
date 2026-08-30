using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "Maxwell tolerance contract" begin
    @test Models._derive_maxwell_solver_tol([1e-4, 5e-5]) == 5e-6
    @test Models._derive_maxwell_solver_tol([1e-4]; factor=0.25) == 2.5e-5
    @test_throws ArgumentError Models._derive_maxwell_solver_tol(Float64[])
    @test_throws ArgumentError Models._derive_maxwell_solver_tol([0.0])
    @test_throws ArgumentError Models._derive_maxwell_solver_tol([1e-4]; factor=0.0)

    @testset "production derives from active acceptance gates" begin
        cfg = Models.ProductionPipelineConfig(
            area_tol_good=1e-4,
            rho_geometry_convergence=true,
            rho_maxwell_area_tol=5e-5,
            adaptive_temperature=false,
            temperature_maxwell_area_tol=1e-6,
        )
        @test Models._production_maxwell_solver_tol(cfg) == 5e-6
        @test Models._production_maxwell_options(cfg).tol_area == 5e-6
        @test Models._production_maxwell_options(cfg).candidate_steps == 1024

        no_rho = Models.ProductionPipelineConfig(
            area_tol_good=1e-4,
            rho_geometry_convergence=false,
            rho_maxwell_area_tol=5e-5,
        )
        @test Models._production_maxwell_solver_tol(no_rho) == 1e-5

        with_temperature = Models.ProductionPipelineConfig(
            area_tol_good=1e-4,
            rho_geometry_convergence=false,
            adaptive_temperature=true,
            temperature_maxwell_area_tol=2e-5,
        )
        @test isapprox(Models._production_maxwell_solver_tol(with_temperature), 2e-6; rtol=1e-12)
    end

    @testset "solver result is checked against the derived tolerance" begin
        good = Models.MaxwellResult(true, 300.0, 0.2, 1.5, 4e-6, 8, Dict{Symbol, Any}())
        bad = Models.MaxwellResult(true, 300.0, 0.2, 1.5, 6e-6, 8, Dict{Symbol, Any}())
        @test Models._maxwell_result_satisfies_tol(good, 5e-6)
        @test !Models._maxwell_result_satisfies_tol(bad, 5e-6)
    end
end
