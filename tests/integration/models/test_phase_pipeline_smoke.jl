using Test
using JSON3

_models_entry = joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    include(_models_entry)
end
using .Models

@testset "Phase pipeline smoke" begin
    tmp = mktempdir()

    result = Models.run_phase_pipeline(
        :PNJL;
        mode=:research,
        T_grid=[150.0],
        rho_grid=[0.1, 0.2, 0.3],
        xi=0.0,
        output_dir=tmp,
        profile=:smoke,
        solver_backend=:models,
        p_num=12,
        t_num=4,
        thermo_quadrature_policy=:rs_reduced_adaptive,
        thermo_quadrature_rtol=1e-6,
        thermo_quadrature_atol=1e-8,
        iterations=80,
        promote_reference=false,
    )

    @test result.model_kind == :PNJL
    @test isfile(joinpath(tmp, "trho_scan.csv"))
    @test haskey(result.artifact_paths, "phase_summary")
    @test haskey(result.artifact_paths, "first_order_boundary")
    @test haskey(result.artifact_paths, "spinodal")
    @test haskey(result.artifact_paths, "crossover_line")
    @test haskey(result.artifact_paths, "phase_grid_convergence")
    @test haskey(result.artifact_paths, "phase_report")
    @test isfile(result.artifact_paths["phase_summary"])
    @test isfile(result.artifact_paths["first_order_boundary"])
    @test isfile(result.artifact_paths["spinodal"])
    @test isfile(result.artifact_paths["crossover_line"])
    @test isfile(result.artifact_paths["phase_grid_convergence"])
    @test isfile(result.artifact_paths["phase_report"])
    @test haskey(result.diagnostics, "scan_total")
    @test result.diagnostics["scan_total"] >= 1
    @test result.config_snapshot["thermo_quadrature_policy"] == "rs_reduced_adaptive"
    @test result.config_snapshot["thermo_quadrature_rtol"] == 1e-6
    @test result.config_snapshot["thermo_quadrature_atol"] == 1e-8
end

@testset "Phase pipeline production smoke" begin
    tmp = mktempdir()

    result = Models.run_phase_pipeline(
        :PNJL;
        mode=:production,
        T_grid=[150.0],
        rho_grid=[0.1, 0.2, 0.3],
        xi=0.0,
        output_dir=tmp,
        profile=:smoke,
        solver_backend=:models,
        p_num=12,
        t_num=4,
        iterations=10,
        promote_reference=false,
    )

    @test isfile(joinpath(tmp, "trho_scan.csv"))
    @test result.config_snapshot["mode"] == "production"
    @test haskey(result.artifact_paths, "phase_summary")
    @test haskey(result.artifact_paths, "phase_report")
    @test haskey(result.artifact_paths, "phase_grid_convergence")
    @test isfile(result.artifact_paths["phase_summary"])
    @test isfile(result.artifact_paths["phase_report"])
    @test isfile(result.artifact_paths["phase_grid_convergence"])
    @test result.config_snapshot["p_num"] == 12
    @test result.config_snapshot["t_num"] == 4
    @test result.config_snapshot["iterations"] == 10
    @test result.config_snapshot["rho_geometry_convergence"] == true

    scan_lines = readlines(joinpath(tmp, "trho_scan.csv"))
    scan_header = split(first(scan_lines), ',')
    key_indices = [findfirst(==(name), scan_header) for name in ("T_MeV", "rho", "xi")]
    @test all(index -> index !== nothing, key_indices)
    resolved_key_indices = Int[something(index) for index in key_indices]
    scan_keys = [join(split(line, ',')[resolved_key_indices], ',') for line in scan_lines[2:end]]
    @test length(scan_keys) == length(unique(scan_keys))

    report_text = read(result.artifact_paths["phase_report"], String)
    @test occursin("## Conclusion", report_text)
    @test occursin("- phase_structure:", report_text)
    @test occursin("- cep_result:", report_text)

    summary = JSON3.read(read(result.artifact_paths["phase_summary"], String))
    @test haskey(summary, "conclusion")
    conclusion = summary["conclusion"]
    @test haskey(conclusion, "phase_structure")
    @test haskey(conclusion, "cep_result")
end
