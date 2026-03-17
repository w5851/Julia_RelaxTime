using Test
using JSON3

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
