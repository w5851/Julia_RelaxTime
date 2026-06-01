using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "multiseed candidate governance" begin
    @testset "priority is hard_constraint > pressure > residual > seed_index" begin
        candidates = [
            (pressure=10.0, residual_norm=1e-8, hard_constraint_ok=false, seed_index=1, converged=true, failed_constraints=Symbol[:residual_too_large]),
            (pressure=8.0, residual_norm=2e-7, hard_constraint_ok=true, seed_index=4, converged=true, failed_constraints=Symbol[]),
            (pressure=9.0, residual_norm=1e-7, hard_constraint_ok=true, seed_index=2, converged=true, failed_constraints=Symbol[]),
            (pressure=12.0, residual_norm=1e-7, hard_constraint_ok=true, seed_index=3, converged=true, failed_constraints=Symbol[]),
            (pressure=12.0, residual_norm=1e-7, hard_constraint_ok=true, seed_index=1, converged=true, failed_constraints=Symbol[]),
        ]

        selected = Models.select_pressure_max_candidate(candidates)
        @test selected.selection_reason == :pressure_max_under_constraints
        @test selected.selected_index == 5
        @test selected.selected_candidate.seed_index == 1
    end

    @testset "fallback pool keeps same tie-break ordering" begin
        candidates = [
            (pressure=7.0, residual_norm=1e-6, hard_constraint_ok=false, seed_index=5, converged=false, failed_constraints=Symbol[:a]),
            (pressure=9.0, residual_norm=1e-6, hard_constraint_ok=false, seed_index=7, converged=false, failed_constraints=Symbol[:b]),
            (pressure=9.0, residual_norm=1e-6, hard_constraint_ok=false, seed_index=3, converged=false, failed_constraints=Symbol[:c]),
            (pressure=9.0, residual_norm=2e-6, hard_constraint_ok=false, seed_index=1, converged=false, failed_constraints=Symbol[:d]),
        ]

        selected = Models.select_pressure_max_candidate(candidates)
        @test selected.selection_reason == :no_candidate_passed_constraints
        @test selected.selected_index == 3
        @test selected.selected_candidate.seed_index == 3
    end
end
