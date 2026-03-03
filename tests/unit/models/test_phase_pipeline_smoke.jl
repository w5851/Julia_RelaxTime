using Test

_models_entry = joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    include(_models_entry)
end
using .Models

@testset "Phase pipeline smoke" begin
    tmp = mktempdir()

    result = Models.run_phase_pipeline(
        :PNJL;
        T_grid=[150.0],
        rho_grid=[0.1, 0.2, 0.3],
        xi=0.0,
        output_dir=tmp,
        profile=:smoke,
        solver_backend=:legacy,
        p_num=12,
        t_num=4,
        iterations=80,
        promote_reference=false,
    )

    @test result.model_kind == :PNJL
    @test isfile(joinpath(tmp, "trho_scan.csv"))
    @test haskey(result.artifact_paths, "phase_summary")
    @test haskey(result.artifact_paths, "first_order_boundary")
    @test haskey(result.artifact_paths, "spinodal")
    @test haskey(result.artifact_paths, "crossover_line")
    @test haskey(result.artifact_paths, "phase_report")
    @test isfile(result.artifact_paths["phase_summary"])
    @test isfile(result.artifact_paths["first_order_boundary"])
    @test isfile(result.artifact_paths["spinodal"])
    @test isfile(result.artifact_paths["crossover_line"])
    @test isfile(result.artifact_paths["phase_report"])
    @test haskey(result.diagnostics, "scan_total")
    @test result.diagnostics["scan_total"] == 3
end
