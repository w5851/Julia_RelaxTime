using Test
using JSON3
using TOML

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const CLI_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "pnjl", "calculate_phase_structure.jl")

include(CLI_SCRIPT)

function _run_phase_cli_direct(; bracket_mode::String, mode::String, start::String="low")
    output_dir = mktempdir()
    cmd = `$(Base.julia_cmd()) --project=$(PROJECT_ROOT) $(CLI_SCRIPT) --model_kind=PNJL --mode=$(mode) --T_min=150 --T_max=160 --T_step=10 --rho_min=0.1 --rho_max=0.3 --rho_step=0.1 --output_dir=$(output_dir) --profile=smoke --solver_backend=legacy --cep_strategy=direct --cep_max_bisect_iter=1 --cep_max_refine_level=0 --cep_direct_bracket_mode=$(bracket_mode) --cep_direct_start=$(start) --cep_direct_expand_factor=2.0 --cep_direct_max_expand_steps=4 --cep_direct_fallback_scan=true`
    run(cmd)
    summary_path = joinpath(output_dir, "phase_summary.json")
    return output_dir, summary_path, JSON3.read(read(summary_path, String))
end

@testset "Phase CLI direct smoke" begin
    @test isfile(CLI_SCRIPT)

    output_dir, summary_path, summary = _run_phase_cli_direct(bracket_mode="directional", mode="research", start="low")

    @test isfile(summary_path)

    @test haskey(summary, "stats")
    @test haskey(summary, "stats_compare")

    stats_compare = summary["stats_compare"]
    @test haskey(stats_compare, "scan")
    @test haskey(stats_compare, "phase")
    @test haskey(stats_compare, "cep")

    cep_view = stats_compare["cep"]
    @test haskey(cep_view, "method")
    @test haskey(cep_view, "eval_count")
    @test haskey(cep_view, "unknown_count")
    @test haskey(cep_view, "unknown_rate")
end

@testset "Phase CLI direct bracket mode compare" begin
    @test isfile(CLI_SCRIPT)

    _, summary_path_directional, summary_directional = _run_phase_cli_direct(bracket_mode="directional", mode="research", start="low")
    _, summary_path_scan, summary_scan = _run_phase_cli_direct(bracket_mode="scan", mode="research", start="low")

    @test isfile(summary_path_directional)
    @test isfile(summary_path_scan)

    @test haskey(summary_directional, "cep")
    @test haskey(summary_scan, "cep")
    @test haskey(summary_directional, "stats_compare")
    @test haskey(summary_scan, "stats_compare")

    directional_cep = summary_directional["cep"]
    scan_cep = summary_scan["cep"]
    directional_cmp_cep = summary_directional["stats_compare"]["cep"]
    scan_cmp_cep = summary_scan["stats_compare"]["cep"]

    @test haskey(directional_cep, "method")
    @test haskey(scan_cep, "method")
    @test haskey(directional_cmp_cep, "eval_count")
    @test haskey(scan_cmp_cep, "eval_count")
    @test haskey(directional_cmp_cep, "unknown_rate")
    @test haskey(scan_cmp_cep, "unknown_rate")

    @test directional_cmp_cep["method"] == directional_cep["method"]
    @test scan_cmp_cep["method"] == scan_cep["method"]
end

@testset "Phase CLI rejects invalid mode" begin
    @test parse_args(String[]).mode == :production

    err = try
        parse_args(["--mode=invalid"])
        nothing
    catch exc
        exc
    end

    @test err isa ArgumentError

    output = sprint(showerror, err)
    @test occursin("invalid --mode=invalid", output)
    @test occursin("production", output)
    @test occursin("research", output)
end

@testset "Phase CLI auto-loads default template by model_kind" begin
    cfg_pnjl = parse_args(String[])
    @test cfg_pnjl.model_kind == :PNJL
    @test cfg_pnjl.config_path !== nothing
    @test occursin(joinpath("config", "models", "pnjl", "phase_pipeline_default.toml"), cfg_pnjl.config_path)
    @test cfg_pnjl.compute_crossover == true

    cfg_rpnjl = parse_args(["--model_kind=RPNJL"])
    @test cfg_rpnjl.model_kind == :RPNJL
    @test cfg_rpnjl.config_path !== nothing
    @test occursin(joinpath("config", "models", "rpnjl", "phase_pipeline_default.toml"), cfg_rpnjl.config_path)
    @test cfg_rpnjl.compute_crossover == true
end

@testset "Phase CLI run without --config uses default template" begin
    output_dir = mktempdir()
    cmd = `$(Base.julia_cmd()) --project=$(PROJECT_ROOT) $(CLI_SCRIPT) --model_kind=PNJL --mode=research --T_min=150 --T_max=150 --T_step=10 --rho_min=0.1 --rho_max=0.3 --rho_step=0.1 --output_dir=$(output_dir) --profile=smoke --solver_backend=legacy --iterations=10`
    run(cmd)

    manifest_path = joinpath(output_dir, "run_manifest.json")
    @test isfile(manifest_path)
    manifest = JSON3.read(read(manifest_path, String))
    @test haskey(manifest, "config_path")
    @test occursin(joinpath("config", "models", "pnjl", "phase_pipeline_default.toml"), String(manifest["config_path"]))
end

@testset "Phase CLI loads config file and allows CLI override" begin
    tmpdir = mktempdir()
    cfg_path = joinpath(tmpdir, "phase_cli.toml")
    open(cfg_path, "w") do io
        TOML.print(io, Dict(
            "phase_pipeline" => Dict(
                "mode" => "research",
                "model_kind" => "RPNJL",
                "T_min" => 140.0,
                "T_max" => 160.0,
                "T_step" => 5.0,
                "rho_min" => 0.2,
                "rho_max" => 0.8,
                "rho_step" => 0.2,
                "profile" => "smoke",
                "solver_backend" => "legacy",
                "seed_policy" => "hybrid_continuity",
                "reverse_rho" => true,
                "compute_crossover" => false,
                "promote_reference" => false,
            )
        ))
    end

    from_cfg = parse_args(["--config=$(cfg_path)"])
    @test from_cfg.mode == :research
    @test from_cfg.model_kind == :RPNJL
    @test from_cfg.T_min == 140.0
    @test from_cfg.rho_max == 0.8

    overridden = parse_args([
        "--config=$(cfg_path)",
        "--mode=production",
        "--model_kind=PNJL",
        "--T_min=150",
        "--rho_max=0.4",
    ])
    @test overridden.mode == :production
    @test overridden.model_kind == :PNJL
    @test overridden.T_min == 150.0
    @test overridden.rho_max == 0.4
end

@testset "Phase CLI writes run manifest" begin
    output_dir = mktempdir()
    cmd = `$(Base.julia_cmd()) --project=$(PROJECT_ROOT) $(CLI_SCRIPT) --model_kind=PNJL --mode=research --T_min=150 --T_max=150 --T_step=10 --rho_min=0.1 --rho_max=0.3 --rho_step=0.1 --output_dir=$(output_dir) --profile=smoke --solver_backend=legacy --iterations=10`
    run(cmd)

    manifest_path = joinpath(output_dir, "run_manifest.json")
    @test isfile(manifest_path)
    manifest = JSON3.read(read(manifest_path, String))
    @test haskey(manifest, "generated_at")
    @test haskey(manifest, "git_commit")
    @test haskey(manifest, "argv")
    @test haskey(manifest, "config_hash")
    @test haskey(manifest, "run_id")
    @test haskey(manifest, "artifact_paths")
end
