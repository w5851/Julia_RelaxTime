using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const P = Models.pnjl_module()

@testset "problem spec semantic modes" begin
    model = Models.create_model(:PNJL)
    mode = Models.FixedEntropy(0.5)
    spec = Models.build_problem_spec(mode)

    T_fm = 100.0 / 197.327
    seed = copy(P.HADRON_SEED_8)

    @testset "ground_state and constrained_manifold contracts" begin
        ground = spec.forward_solve(
            model,
            T_fm;
            seed_guess=seed,
            rho0=0.16,
            semantic_mode=:ground_state,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
            iterations=120,
        )

        manifold = spec.forward_solve(
            model,
            T_fm;
            seed_guess=seed,
            rho0=0.16,
            semantic_mode=:constrained_manifold,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
            iterations=120,
        )

        @test ground.selection_reason in (:pressure_max_under_constraints, :no_candidate_passed_constraints)
        @test manifold.selection_reason in (:residual_min_under_constraints, :no_candidate_passed_constraints)
        @test ground.candidate_count >= 1
        @test manifold.candidate_count >= 1
    end

    @testset "selector ordering is reproducible" begin
        r1 = spec.forward_solve(
            model,
            T_fm;
            seed_guess=seed,
            rho0=0.16,
            semantic_mode=:constrained_manifold,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
            iterations=120,
        )

        r2 = spec.forward_solve(
            model,
            T_fm;
            seed_guess=seed,
            rho0=0.16,
            semantic_mode=:constrained_manifold,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
            iterations=120,
        )

        @test r1.selected_index == r2.selected_index
        @test r1.selection_reason == r2.selection_reason
    end

    @testset "custom selector hook" begin
        custom_selector = candidates -> begin
            selected = candidates[end]
            return (
                selected_index=length(candidates),
                selected_candidate=selected,
                selection_reason=:custom_selector,
                eligible_indices=Int[length(candidates)],
            )
        end

        custom = spec.forward_solve(
            model,
            T_fm;
            seed_guess=seed,
            rho0=0.16,
            selector=custom_selector,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
            iterations=120,
        )

        @test custom.selection_reason == :custom_selector
    end
end
