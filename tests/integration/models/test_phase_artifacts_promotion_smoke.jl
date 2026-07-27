using Test

_models_entry = joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    include(_models_entry)
end
using .Models

@testset "Phase artifacts promotion smoke" begin
    temp_root = mktempdir()
    run_id = "phase_smoke_promotion"

    target = Models.resolve_phase_output_target(
        :PNJL;
        profile=:smoke,
        run_id=run_id,
        policy=:processed_first,
        project_root=temp_root,
    )

    @test endswith(normpath(target.run_dir), normpath(joinpath("data", "processed", "pnjl", "phase_diagram", run_id)))

    result = Models.run_phase_pipeline(
        :PNJL;
        mode=:research,
        T_grid=[150.0],
        rho_grid=[0.1, 0.2, 0.3],
        xi=0.0,
        output_dir=target.run_dir,
        profile=:smoke,
        run_id=run_id,
        solver_backend=:models,
        p_num=12,
        t_num=4,
        iterations=80,
        promote_reference=false,
    )

    @test isfile(joinpath(target.run_dir, "phase_summary.json"))
    @test isfile(joinpath(target.run_dir, "first_order_boundary.csv"))
    @test isfile(joinpath(target.run_dir, "spinodal.csv"))
    @test isfile(joinpath(target.run_dir, "crossover_line.csv"))
    @test isfile(joinpath(target.run_dir, "phase_report.md"))
    @test result.promotion_status.passed == false

    promoted = Models.promote_phase_artifacts(
        target.run_dir;
        reference_root=target.reference_root,
        gate_options=(; profile="smoke", expected_model_kind="PNJL", expected_schema_version="phase-v2"),
        write_reference=true,
    )

    @test promoted.passed
    @test promoted.baseline_id !== nothing
    @test promoted.reference_dir !== nothing
    @test isdir(promoted.reference_dir)
    @test isfile(joinpath(promoted.reference_dir, "phase_summary.json"))

    bad_dir = mktempdir()
    write(joinpath(bad_dir, "phase_summary.json"), "{}")
    for name in ("first_order_boundary.csv", "spinodal.csv", "crossover_line.csv", "phase_report.md")
        write(joinpath(bad_dir, name), "")
    end

    rejected = Models.promote_phase_artifacts(
        bad_dir;
        reference_root=target.reference_root,
        gate_options=(; expected_model_kind="PNJL", expected_schema_version="phase-v2"),
        write_reference=false,
    )

    @test !rejected.passed
    @test !isempty(rejected.failed_checks)
end
