using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Constants_PNJL)
    include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
end

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "problem spec contract" begin
    @test isdefined(Models, :ProblemSpec)
    @test isdefined(Models, :build_problem_spec)

    @testset "mode dimensions" begin
        modes = [
            Models.FixedMu(),
            Models.FixedRho(1.0),
            Models.FixedAsymmetricRho(1.0, 0.876, 0.0),
            Models.FixedEntropy(0.5),
            Models.FixedSigma(8.0),
        ]

        for mode in modes
            spec = Models.build_problem_spec(mode)
            @test spec isa Models.ProblemSpec
            @test spec.mode == mode
            @test spec.x_dim == Models.state_dim(mode)
            @test spec.theta_dim == Models.param_dim(mode)
            @test spec.residual! isa Function
            @test spec.forward_solve isa Function
            @test spec.conditions isa Function
            @test spec.unpack_solution isa Function
            @test spec.postprocess isa Function
            @test spec.hard_rules isa AbstractVector
            @test spec.selector isa Function
        end
    end

    @testset "fixedrho spec conditions and forward solve are wired" begin
        model = Models.create_model(:PNJL)
        mode = Models.FixedRho(0.2)
        spec = Models.build_problem_spec(mode)

        params = Models.GapParams(0.5, Models.cached_nodes(8, 4), 0.0)
        cond = spec.conditions(params)
        residual = cond([0.5], [-1.5, -1.5, -2.1, 0.2, 0.2, 1.5, 1.5, 1.5])
        @test length(residual) == 8

        solved = spec.forward_solve(model, 0.5; seed_guess=[-1.5, -1.5, -2.1, 0.2, 0.2, 1.5, 1.5, 1.5], p_num=8, t_num=4)
        @test solved isa NamedTuple
        @test haskey(solved, :converged)
        @test haskey(solved, :residual_norm)
        @test haskey(solved, :selection_reason)
        @test haskey(solved, :candidate_count)
    end

    @testset "solve_constraint fixedrho prefers problem_spec forward_solve" begin
        model = Models.create_model(:PNJL)
        mode = Models.FixedRho(0.2)
        spec = Models.ProblemSpec(
            mode;
            x_dim=8,
            theta_dim=1,
            forward_solve=(m, T_fm; fwd_kwargs...) -> (
                converged=true,
                residual_norm=0.0,
                pressure=0.0,
                used_problem_spec=true,
            ),
        )

        result = Models.solve_constraint(
            model,
            mode,
            0.5;
            problem_spec=spec,
        )

        @test result isa NamedTuple
        @test result.converged
        @test result.used_problem_spec
        @test result.residual_norm == 0.0
    end

    @testset "solve_constraint supports problem_spec override for other modes" begin
        model = Models.create_model(:PNJL)
        modes = (
            Models.FixedEntropy(0.5),
            Models.FixedSigma(8.0),
            Models.FixedAsymmetricRho(0.6, 1.0, 0.0),
        )

        for mode in modes
            spec = Models.ProblemSpec(
                mode;
                x_dim=Models.state_dim(mode),
                theta_dim=Models.param_dim(mode),
                forward_solve=(m, T_fm; fwd_kwargs...) -> (
                    converged=true,
                    residual_norm=0.0,
                    mode_tag=string(typeof(mode)),
                ),
            )

            result = Models.solve_constraint(model, mode, 0.5; problem_spec=spec)
            @test result.converged
            @test result.residual_norm == 0.0
            @test occursin("Fixed", result.mode_tag)
        end
    end

    @testset "solve_constraint validates problem_spec type" begin
        model = Models.create_model(:PNJL)
        mode = Models.FixedEntropy(0.5)
        @test_throws ArgumentError Models.solve_constraint(model, mode, 0.5; problem_spec=:invalid)
    end

    @testset "solve_constraint defaults to ProblemSpec chain and gates legacy fallback" begin
        model = Models.create_model(:PNJL)
        mode = Models.FixedEntropy(0.5)
        T_fm = 100.0 / 197.327
        seed = copy(Models.pnjl_module().HADRON_SEED_8)

        auto_spec = Models.solve_constraint(
            model,
            mode,
            T_fm;
            seed_guess=seed,
            rho0=0.16,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
            iterations=120,
        )
        @test haskey(auto_spec, :selection_reason)
        @test haskey(auto_spec, :candidate_count)

        @test_throws ArgumentError Models.solve_constraint(
            model,
            mode,
            T_fm;
            use_problem_spec=false,
            seed_guess=seed,
            rho0=0.16,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
            iterations=120,
        )

        legacy = Models.solve_constraint(
            model,
            mode,
            T_fm;
            use_problem_spec=false,
            allow_legacy_path=true,
            seed_guess=seed,
            rho0=0.16,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
            iterations=120,
        )
        @test !haskey(legacy, :selection_reason)
        @test !haskey(legacy, :candidate_count)

        @test_logs (:warn, r"use_problem_spec=false") Models.solve_constraint(
            model,
            mode,
            T_fm;
            use_problem_spec=false,
            allow_legacy_path=true,
            seed_guess=seed,
            rho0=0.16,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
            iterations=120,
        )

        @test_logs Models.solve_constraint(
            model,
            mode,
            T_fm;
            use_problem_spec=false,
            allow_legacy_path=true,
            warn_on_legacy_path=false,
            seed_guess=seed,
            rho0=0.16,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
            iterations=120,
        )

        @test_throws ArgumentError Models.solve_constraint(
            model,
            mode,
            T_fm;
            use_problem_spec=:yes,
            seed_guess=seed,
            rho0=0.16,
            p_num=8,
            t_num=4,
        )

        @test_throws ArgumentError Models.solve_constraint(
            model,
            mode,
            T_fm;
            use_problem_spec=false,
            allow_legacy_path=true,
            warn_on_legacy_path=:invalid,
            seed_guess=seed,
            rho0=0.16,
            p_num=8,
            t_num=4,
        )

        @test_throws ArgumentError Models.solve_constraint(
            model,
            mode,
            T_fm;
            use_problem_spec=false,
            allow_legacy_path=:invalid,
            seed_guess=seed,
            rho0=0.16,
            p_num=8,
            t_num=4,
        )
    end

    @testset "build_problem_spec non-rho forward_solve exposes governed metadata" begin
        model = Models.create_model(:PNJL)
        T_fm = 100.0 / 197.327
        seed = copy(Models.pnjl_module().HADRON_SEED_8)

        modes = (
            Models.FixedEntropy(0.5),
            Models.FixedSigma(10.0),
            Models.FixedAsymmetricRho(0.05, 1.0, 0.0),
        )

        for mode in modes
            spec = Models.build_problem_spec(mode)
            solved = spec.forward_solve(
                model,
                T_fm;
                seed_guess=seed,
                rho0=0.16,
                p_num=8,
                t_num=4,
                residual_norm_max=1e-6,
                iterations=120,
            )

            @test solved isa NamedTuple
            @test haskey(solved, :selection_reason)
            @test haskey(solved, :candidate_count)
            @test solved.selection_reason in (:pressure_max_under_constraints, :residual_min_under_constraints, :no_candidate_passed_constraints)
            @test solved.candidate_count >= 1
        end
    end

    @testset "problem spec non-rho semantic_mode selects candidate selector" begin
        model = Models.create_model(:PNJL)
        mode = Models.FixedEntropy(0.5)
        spec = Models.build_problem_spec(mode)
        T_fm = 100.0 / 197.327
        seed = copy(Models.pnjl_module().HADRON_SEED_8)

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
        @test ground.selection_reason in (:pressure_max_under_constraints, :no_candidate_passed_constraints)

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
        @test manifold.selection_reason in (:residual_min_under_constraints, :no_candidate_passed_constraints)

        @test_throws ArgumentError spec.forward_solve(
            model,
            T_fm;
            seed_guess=seed,
            rho0=0.16,
            semantic_mode=:bad_mode,
            p_num=8,
            t_num=4,
        )

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
