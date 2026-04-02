using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SCRIPT_PATH = joinpath(PROJECT_ROOT, "scripts", "pnjl", "run_tmu_scan.jl")

@testset "run_tmu_scan script mode contract" begin
    if !isfile(SCRIPT_PATH)
        include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
        @test !isdefined(Main.Models, :scan_workflow_migration_status)
        return
    end

    include(SCRIPT_PATH)

    @testset "parse_args supports point mode" begin
        cfg = parse_args([
            "--mode=point",
            "--T_mev=150",
            "--mu_mev=0",
            "--xi=0.0",
            "--output=$(joinpath(mktempdir(), "tmu_point.csv"))",
            "--overwrite",
            "--no_phase_aware",
            "--p_num=12",
            "--t_num=4",
        ])

        @test cfg.mode == :point
        @test cfg.T_mev == 150.0
        @test cfg.mu_mev == 0.0
    end

    @testset "parse_args keeps scan mode" begin
        cfg = parse_args([
            "--mode=scan",
            "--T_min=140",
            "--T_max=140",
            "--T_step=10",
            "--mu_min=10",
            "--mu_max=10",
            "--mu_step=10",
        ])

        @test cfg.mode == :scan
        @test cfg.T_min == 140.0
        @test cfg.mu_min == 10.0
    end

    @testset "point mode executes single point" begin
        out = joinpath(mktempdir(), "tmu_point_exec.csv")
        stats = main([
            "--mode=point",
            "--T_mev=150",
            "--mu_mev=0",
            "--xi=0.0",
            "--output=$(out)",
            "--overwrite",
            "--no_phase_aware",
            "--p_num=12",
            "--t_num=4",
        ])

        @test isfile(out)
        @test stats.total == 1
        @test stats.success == 1

        lines = readlines(out)
        @test length(lines) == 2
    end
end
