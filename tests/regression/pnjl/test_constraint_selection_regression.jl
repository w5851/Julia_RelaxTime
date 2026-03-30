using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "PNJL constraint selection regression" begin
    # 反例：高 pressure 但违反硬约束，不能被选中
    candidates = [
        (pressure=1.0, residual_norm=1e-8, hard_constraint_ok=true, failed_constraints=Symbol[], converged=true),
        (pressure=2.0, residual_norm=1e-10, hard_constraint_ok=false, failed_constraints=[:phi_out_of_range], converged=true),
    ]

    selected = Models.select_pressure_max_candidate(candidates)
    @test selected.selection_reason == :pressure_max_under_constraints
    @test selected.selected_index == 1
    @test selected.selected_candidate.hard_constraint_ok
end
