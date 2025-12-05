using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

using .PNJL.MaxwellRhoMu
using .PNJL.PhaseTransition

function synthetic_curve()
    mu_vals = [-40.0, -25.0, -15.0, -5.0, 0.0, 5.0, 10.0, 15.0, 25.0]
    rho_vals = [-0.6, -0.3, -0.1, 0.25, 0.55, 0.35, 0.15, 0.4, 0.65]
    return mu_vals, rho_vals
end

@testset "maxwell_rho_mu basic" begin
    mu_vals, rho_vals = synthetic_curve()
    result = maxwell_rho_mu(mu_vals, rho_vals; min_samples=5)
    @test result.converged
    @test result.mu_coex_MeV !== nothing
    @test result.rho_gas < result.rho_liquid
    @test result.area_residual ≈ 0.0 atol=1e-2
end

@testset "maxwell respects spinodal hint" begin
    mu_vals, rho_vals = synthetic_curve()
    hint = detect_s_shape(mu_vals, rho_vals; min_points=5)
    result = maxwell_rho_mu(mu_vals, rho_vals; spinodal_hint=hint, min_samples=5)
    @test result.converged
end

@testset "build_phase_boundary" begin
    mu_vals, rho_vals = synthetic_curve()
    curves = Dict(70.0 => (mu_vals, rho_vals), 90.0 => (mu_vals, rho_vals))
    results = build_phase_boundary(curves; min_samples=5)
    @test length(results) == 2
    @test all(res -> res.converged, values(results))
end

@testset "phase_boundary_from_rows filters xi" begin
    mu_vals, rho_vals = synthetic_curve()
    rows = [
        Dict("T_MeV" => "70", "xi" => "0.0", "rho" => string(rho), "mu_avg_MeV" => string(mu))
        for (mu, rho) in zip(mu_vals, rho_vals)
    ]
    rows_bad = [Dict("T_MeV" => "70", "xi" => "0.5", "rho" => "0.1", "mu_avg_MeV" => "5.0")]
    all_rows = vcat(rows, rows_bad)
    results = phase_boundary_from_rows(all_rows; xi=0.0, min_samples=5)
    @test haskey(results, 70.0)
    @test results[70.0].converged
end
