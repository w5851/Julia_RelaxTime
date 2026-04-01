using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "candidate governance stability" begin
    context = Models.build_candidate_context(Models.FixedMu();
        continuity_seed_available=true,
        phase_hint=:hadron,
        residual_norm_max=1e-6,
        prefer_continuity=true,
    )

    candidates = [
        (pressure=10.0, residual_norm=1e-8, hard_constraint_ok=true, failed_constraints=Symbol[], converged=true, seed_source=:continuity),
        (pressure=12.0, residual_norm=1e-8, hard_constraint_ok=false, failed_constraints=[:mass_nonpositive], converged=true, seed_source=:default),
        (pressure=9.5, residual_norm=1e-7, hard_constraint_ok=true, failed_constraints=Symbol[], converged=true, seed_source=:default),
    ]

    selected = Models.select_pressure_max_candidate(candidates, nothing, context)
    @test selected.selected_index == 1
    @test selected.selection_reason == :pressure_max_under_constraints
    @test selected.selected_candidate.seed_source == :continuity
end
