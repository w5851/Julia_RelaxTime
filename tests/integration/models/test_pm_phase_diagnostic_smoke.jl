using Test
using JSON3

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const CLI_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "pnjl", "diagnose_pm_phase.jl")

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "PM phase diagnostic CLI smoke" begin
    output_dir = mktempdir()
    cmd = `$(Base.julia_cmd()) --project=$(PROJECT_ROOT) $(CLI_SCRIPT) --T_values=130.9 --mu_start=290.9 --mu_stop=291.1 --mu_step=0.1 --output_dir=$(output_dir) --solver_backend=models --p_num=24 --t_num=8 --xi=0.0`
    run(cmd)

    @test isfile(joinpath(output_dir, "pm_branch_scan.csv"))
    @test isfile(joinpath(output_dir, "pm_phase_summary.json"))
    @test isfile(joinpath(output_dir, "pm_vs_maxwell.csv"))
end

@testset "PM phase diagnostic smoke" begin
    output_dir = mktempdir()
    result = Models.analyze_pm_branch_competition(
        T_values=[130.9],
        mu_grid=[290.9, 291.0, 291.1],
        xi=0.0,
        solver_backend=:models,
        p_num=24,
        t_num=8,
        output_dir=output_dir,
    )

    @test haskey(result, :branch_rows)
    @test haskey(result, :temperature_summaries)
    @test haskey(result, :comparison_rows)

    @test !isempty(result.branch_rows)
    @test any(row.branch == :hadron for row in result.branch_rows)
    @test any(row.branch == :quark for row in result.branch_rows)

    @test isfile(joinpath(output_dir, "pm_branch_scan.csv"))
    @test isfile(joinpath(output_dir, "pm_phase_summary.json"))
    @test isfile(joinpath(output_dir, "pm_vs_maxwell.csv"))

    summary = JSON3.read(read(joinpath(output_dir, "pm_phase_summary.json"), String))
    @test length(summary) == 1
    summary_row = summary[1]
    @test haskey(summary_row, "T_MeV")
    @test haskey(summary_row, "mu_transition_pm_MeV")
    @test haskey(summary_row, "hadron_endpoint_mu_MeV")
    @test haskey(summary_row, "quark_endpoint_mu_MeV")
    @test haskey(summary_row, "bistable_window_width_MeV")
    @test haskey(summary_row, "comparison_status")
    @test haskey(summary_row, "comparison_mu_tol_MeV")
    @test haskey(summary_row, "residual_accept_tol")
    @test haskey(summary_row, "continuity_x_tol")
    @test haskey(summary_row, "continuity_rho_tol")
    @test summary_row["mu_transition_pm_MeV"] === nothing || abs(summary_row["mu_transition_pm_MeV"] * 100 - round(summary_row["mu_transition_pm_MeV"] * 100)) <= 1e-8

    compare_lines = readlines(joinpath(output_dir, "pm_vs_maxwell.csv"))
    @test length(compare_lines) >= 2
    @test compare_lines[1] == "T_MeV,mu_transition_pm_MeV,mu_transition_maxwell_MeV,delta_mu_pm_minus_maxwell_MeV,comparison_status,branch_disappears_first,hadron_endpoint_mu_MeV,quark_endpoint_mu_MeV,bistable_window_width_MeV"
end
