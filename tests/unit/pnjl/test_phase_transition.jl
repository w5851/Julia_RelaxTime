using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

using .PNJL.PhaseTransition

@testset "detect_s_shape" begin
    mu_vals = [-50.0, -20.0, -10.0, 0.0, 5.0, 15.0, 30.0]
    rho_vals = [-0.6, -0.25, -0.1, 0.2, 0.35, 0.30, 0.45]
    result = detect_s_shape(mu_vals, rho_vals)
    @test result.has_s_shape
    @test result.mu_spinodal_low !== nothing
    @test result.mu_spinodal_high !== nothing
end

@testset "monotonic curve" begin
    mu_vals = collect(-30.0:10.0:30.0)
    rho_vals = collect(-0.3:0.1:0.3)
    result = detect_s_shape(mu_vals, rho_vals)
    @test !result.has_s_shape
end

@testset "group_curves" begin
    rows = [
        Dict("T_MeV" => "80", "xi" => "0.0", "rho" => "0.1", "mu_avg_MeV" => "10"),
        Dict("T_MeV" => "80", "xi" => "0.0", "rho" => "0.2", "mu_avg_MeV" => "20"),
        Dict("T_MeV" => "90", "xi" => "0.1", "rho" => "0.3", "mu_avg_MeV" => "30"),
    ]
    grouped = group_curves_by_temperature(rows; xi=0.0)
    @test length(grouped) == 1
    @test haskey(grouped, 80.0)
    @test length(grouped[80.0]) == 2
end
