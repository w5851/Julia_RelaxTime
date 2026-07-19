using Test
using JSON3
using TOML

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const CLI_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "pnjl", "calculate_phase_structure.jl")

include(CLI_SCRIPT)

function _run_phase_cli_direct(; bracket_mode::String, mode::String, start::String="low")
    output_dir = mktempdir()
    cmd = `$(Base.julia_cmd()) --project=$(PROJECT_ROOT) $(CLI_SCRIPT) --model_kind=PNJL --mode=$(mode) --T_min=150 --T_max=160 --T_step=10 --rho_min=0.1 --rho_max=0.3 --rho_step=0.1 --output_dir=$(output_dir) --profile=smoke --solver_backend=models --cep_strategy=direct --cep_max_bisect_iter=1 --cep_max_refine_level=0 --cep_direct_bracket_mode=$(bracket_mode) --cep_direct_start=$(start) --cep_direct_expand_factor=2.0 --cep_direct_max_expand_steps=4 --cep_direct_fallback_scan=true`
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
    @test isnan(cfg_pnjl.crossover_T_max_MeV)
    @test cfg_pnjl.rho_geometry_convergence == true
    @test cfg_pnjl.adaptive_temperature == false

    cfg_rpnjl = parse_args(["--model_kind=RPNJL"])
    @test cfg_rpnjl.model_kind == :RPNJL
    @test cfg_rpnjl.config_path !== nothing
    @test occursin(joinpath("config", "models", "rpnjl", "phase_pipeline_default.toml"), cfg_rpnjl.config_path)
    @test cfg_rpnjl.compute_crossover == true
end

@testset "Phase CLI parses phase-grid convergence controls" begin
    cfg = parse_args([
        "--mode=production",
        "--crossover_T_max_MeV=260",
        "--unknown_budget=7",
        "--rho_geometry_convergence=true",
        "--rho_position_tol_MeV=0.025",
        "--rho_density_tol=0.0025",
        "--rho_maxwell_area_tol=5e-5",
        "--adaptive_temperature=true",
        "--temperature_max_refine_level=3",
        "--temperature_position_tol_MeV=0.05",
        "--temperature_density_tol=0.005",
        "--temperature_maxwell_area_tol=5e-5",
    ])

    @test cfg.crossover_T_max_MeV == 260.0
    @test cfg.unknown_budget == 7
    @test cfg.rho_geometry_convergence == true
    @test cfg.rho_position_tol_MeV == 0.025
    @test cfg.rho_density_tol == 0.0025
    @test cfg.rho_maxwell_area_tol == 5e-5
    @test cfg.adaptive_temperature == true
    @test cfg.temperature_max_refine_level == 3
    @test cfg.temperature_position_tol_MeV == 0.05
    @test cfg.temperature_density_tol == 0.005
    @test cfg.temperature_maxwell_area_tol == 5e-5

    @test_throws ArgumentError parse_args(["--mode=production", "--cep_max_refine_level=0"])
    @test_throws ArgumentError parse_args(["--rho_density_tol=0"])
    @test_throws ArgumentError parse_args(["--crossover_T_max_MeV=20"])
end

@testset "Phase CLI run without --config uses default template" begin
    output_dir = mktempdir()
    cmd = `$(Base.julia_cmd()) --project=$(PROJECT_ROOT) $(CLI_SCRIPT) --model_kind=PNJL --mode=research --T_min=150 --T_max=150 --T_step=10 --rho_min=0.1 --rho_max=0.3 --rho_step=0.1 --output_dir=$(output_dir) --profile=smoke --solver_backend=models --iterations=10`
    run(cmd)

    manifest_path = joinpath(output_dir, "run_manifest.json")
    @test isfile(manifest_path)
    manifest = JSON3.read(read(manifest_path, String))
    @test haskey(manifest, "config_path")
    @test occursin("config/models/pnjl/phase_pipeline_default.toml", String(manifest["config_path"]))
    @test !occursin('\\', String(manifest["config_path"]))
end

@testset "Phase CLI supports --preset=smoke" begin
    cfg = parse_args(["--preset=smoke"])
    @test cfg.mode == :research
    @test cfg.profile == :smoke
    @test cfg.solver_backend == :models
    @test cfg.iterations == 10
    @test cfg.p_num == 12
    @test cfg.t_num == 4
end

@testset "Phase CLI parses thermal quadrature controls" begin
    cfg = parse_args([
        "--thermo_quadrature_policy=rs_reduced_adaptive",
        "--thermo_quadrature_rtol=1e-7",
        "--thermo_quadrature_atol=1e-9",
        "--thermo_quadrature_maxevals=12345",
    ])
    @test cfg.thermo_quadrature_policy === :rs_reduced_adaptive
    @test cfg.thermo_quadrature_rtol == 1e-7
    @test cfg.thermo_quadrature_atol == 1e-9
    @test cfg.thermo_quadrature_maxevals == 12345
    @test_throws ArgumentError parse_args(["--thermo_quadrature_policy=unknown"])
    @test_throws ArgumentError parse_args(["--thermo_quadrature_rtol=Inf"])
end

@testset "Phase CLI keeps explicit args over --preset=smoke" begin
    cfg = parse_args([
        "--preset=smoke",
        "--iterations=77",
        "--mode=production",
    ])
    @test cfg.iterations == 77
    @test cfg.mode == :production
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
                "solver_backend" => "models",
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
    cmd = `$(Base.julia_cmd()) --project=$(PROJECT_ROOT) $(CLI_SCRIPT) --model_kind=PNJL --mode=research --T_min=150 --T_max=150 --T_step=10 --rho_min=0.1 --rho_max=0.3 --rho_step=0.1 --output_dir=$(output_dir) --profile=smoke --solver_backend=models --iterations=10`
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
    artifact_paths = manifest["artifact_paths"]
    for (_, p) in pairs(artifact_paths)
        p_str = String(p)
        @test !occursin('\\', p_str)
        @test !occursin(r"^[A-Za-z]:", p_str)
    end
    @test haskey(manifest, "effective_config")
    effective = manifest["effective_config"]
    @test String(effective["thermo_quadrature_policy"]) == "tensor_gauss"
    @test Float64(effective["thermo_quadrature_rtol"]) == 1e-8
    @test Float64(effective["thermo_quadrature_atol"]) == 1e-10
    @test Int(effective["thermo_quadrature_maxevals"]) == 10^7
    @test Int(effective["p_num"]) == 12
    @test Int(effective["t_num"]) == 4
    @test Int(effective["iterations"]) == 10
    @test Float64(effective["crossover_T_max_MeV"]) == 150.0
    @test Float64(effective["cep_tol"]) == 0.01
    @test Bool(effective["rho_geometry_convergence"]) == true
    @test Bool(effective["adaptive_temperature"]) == false
    @test Float64(effective["rho_position_tol_MeV"]) == 0.05
    @test Float64(effective["temperature_position_tol_MeV"]) == 0.10
end

@testset "Phase CLI manifest includes preset and effective config" begin
    output_dir = mktempdir()
    cmd = `$(Base.julia_cmd()) --project=$(PROJECT_ROOT) $(CLI_SCRIPT) --preset=smoke --iterations=77 --output_dir=$(output_dir)`
    run(cmd)

    manifest_path = joinpath(output_dir, "run_manifest.json")
    @test isfile(manifest_path)
    manifest = JSON3.read(read(manifest_path, String))

    @test haskey(manifest, "preset")
    @test String(manifest["preset"]) == "smoke"

    @test haskey(manifest, "effective_config")
    effective = manifest["effective_config"]
    @test haskey(effective, "iterations")
    @test Int(effective["iterations"]) == 77
    @test haskey(effective, "profile")
    @test String(effective["profile"]) == "smoke"
    @test String(effective["thermo_quadrature_policy"]) == "tensor_gauss"
end
