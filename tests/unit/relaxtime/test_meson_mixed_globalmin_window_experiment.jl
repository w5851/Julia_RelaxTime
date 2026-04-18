using Test

const SCRIPT_PATH = joinpath(@__DIR__, "..", "..", "..", "scripts", "analysis", "relaxtime", "meson_mixed", "run_globalmin_window_experiment.jl")

if !isdefined(Main, :_parse_args)
    include(SCRIPT_PATH)
end

@testset "meson mixed globalmin window script" begin
    @test isfile(SCRIPT_PATH)

    opts = Main._parse_args([
        "--outdir", "tmp_out",
        "--libs", "nlopt,metaheuristics,evolutionary",
        "--tmin", "220",
        "--tmax", "224",
        "--tstep", "2",
        "--jump-threshold", "0.3",
    ])

    @test opts.outdir == "tmp_out"
    @test opts.libs == [:nlopt, :metaheuristics, :evolutionary]
    @test opts.T_min_MeV == 220.0
    @test opts.T_max_MeV == 224.0
    @test opts.T_step_MeV == 2.0
    @test opts.jump_threshold == 0.3
end

@testset "parser accepts next-batch libraries" begin
    opts = Main._parse_args([
        "--outdir", "tmp_out",
        "--libs", "cmaes,optim_samin",
    ])
    @test opts.libs == [:cmaes, :optim_samin]
end

@testset "jump summary counts threshold crossings" begin
    rows = [
        Dict{String,Any}("lib" => "nlopt", "T_MeV" => 220.0, "M_eta" => 2.0, "M_eta_prime" => 2.4, "Gamma_eta" => 0.1, "Gamma_eta_prime" => 0.2),
        Dict{String,Any}("lib" => "nlopt", "T_MeV" => 222.0, "M_eta" => 2.0, "M_eta_prime" => 2.8, "Gamma_eta" => 0.5, "Gamma_eta_prime" => 0.2),
        Dict{String,Any}("lib" => "nlopt", "T_MeV" => 224.0, "M_eta" => 2.2, "M_eta_prime" => 2.9, "Gamma_eta" => 0.55, "Gamma_eta_prime" => 0.65),
    ]

    summary = Main._build_jump_summary(rows, 0.3)
    @test length(summary) == 1
    s = summary[1]
    @test s["lib"] == "nlopt"
    @test s["jump_count_Gamma_eta"] == 1
    @test s["jump_count_Gamma_eta_prime"] == 1
    @test s["jump_count_M_eta_prime"] == 1
end
