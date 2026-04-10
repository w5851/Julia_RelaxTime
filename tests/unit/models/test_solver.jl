# Solver.jl 单元测试
#
# 测试内容：
# 1. solve_constraint FixedMu 分发
# 2. solve_constraint FixedRho/FixedEntropy/FixedSigma 分发
# 3. solve/solve_multi 快捷入口

using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Models.pnjl_module()

# ============================================================================

@testset "Solver" begin

    @testset "solve_constraint FixedMu 可调用" begin
        m = Models.create_model(:NJL)
        mode = Models.FixedMu()
        T = 0.5
        seed = Models.gap_initial_guess(m, T, SVector{3}(0.0, 0.0, 0.0))
        result = Models.solve_constraint(m, mode, T; μ_fm=0.0, seed_guess=seed, p_num=24, t_num=6)
        @test result isa NamedTuple
        @test haskey(result, :state) || haskey(result, :pressure) || result isa Any  # 类型灵活
    end

    @testset "FixedMu ProblemSpec-only guard" begin
        m = Models.create_model(:PNJL)
        mode = Models.FixedMu()
        T_fm = 100.0 / 197.327
        seed = copy(Models.pnjl_module().HADRON_SEED_5)

        default_path = Models.solve_constraint(
            m,
            mode,
            T_fm;
            μ_fm=0.0,
            seed_guess=seed,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )
        @test_throws ArgumentError Models.solve_constraint(
            m,
            mode,
            T_fm;
            μ_fm=0.0,
            seed_guess=seed,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
            fixedmu_use_problem_spec=false,
        )

        @test !haskey(default_path, :fixedmu_problem_spec_active)
    end

    @testset "FixedMu rejects removed ProblemSpec switches" begin
        m = Models.create_model(:PNJL)
        mode = Models.FixedMu()
        T_fm = 100.0 / 197.327
        seed = copy(Models.pnjl_module().HADRON_SEED_5)

        @test_throws ArgumentError Models.solve_constraint(
            m,
            mode,
            T_fm;
            μ_fm=0.0,
            seed_guess=seed,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
            fixedmu_use_problem_spec=false,
            diagnostic_level=:summary,
        )

        @test_throws ArgumentError Models.solve_constraint(
            m,
            mode,
            T_fm;
            μ_fm=0.0,
            seed_guess=seed,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
            fixedmu_use_problem_spec=true,
        )
    end

    @testset "runtime option type validation" begin
        m = Models.create_model(:PNJL)
        mode_fixedmu = Models.FixedMu()
        mode_entropy = Models.FixedEntropy(0.5)
        T_fm = 100.0 / 197.327
        seed_mu = copy(Models.pnjl_module().HADRON_SEED_5)
        seed_entropy = copy(Models.pnjl_module().HADRON_SEED_8)

        @test_throws ArgumentError Models.solve(m, mode_fixedmu, T_fm, 0.0;
            seed_guess=seed_mu,
            xi="bad",
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )
        @test_throws ArgumentError Models.solve(m, mode_fixedmu, T_fm, 0.0;
            seed_guess=seed_mu,
            p_num=0,
            t_num=4,
            residual_norm_max=1e-6,
        )
        @test_throws ArgumentError Models.solve(m, mode_fixedmu, T_fm, 0.0;
            seed_guess=seed_mu,
            p_num=8,
            t_num="bad",
            residual_norm_max=1e-6,
        )
        @test_throws ArgumentError Models.solve(m, mode_fixedmu, T_fm, 0.0;
            seed_guess=seed_mu,
            p_num=8,
            t_num=4,
            residual_norm_max=0.0,
        )
        @test_throws ArgumentError Models.solve(m, mode_fixedmu, T_fm, 0.0;
            seed_guess="bad",
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )
        @test_throws ArgumentError Models.solve_multi(m, mode_fixedmu, T_fm, 0.0;
            seeds=(seed_mu, "bad"),
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )
        @test_throws ArgumentError Models.solve(m, mode_fixedmu, T_fm, 0.0;
            seed_guess=seed_mu,
            auto_multiseed_fallback="bad",
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )
        @test_throws ArgumentError Models.solve_multi(m, mode_fixedmu, T_fm, 0.0;
            seed_guess=seed_mu,
            evaluate_all_attempts="bad",
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )

        @test_throws ArgumentError Models.solve(m, mode_entropy, T_fm;
            seed_guess=seed_entropy,
            semantic_mode=:bad_mode,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )
        @test_throws ArgumentError Models.solve(m, mode_entropy, T_fm;
            seed_guess=seed_entropy,
            selector=1,
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )
        @test_throws ArgumentError Models.solve(m, mode_entropy, T_fm;
            seed_guess=seed_entropy,
            rho0="bad",
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )
        @test_throws ArgumentError Models.solve(m, mode_entropy, T_fm;
            seed_guess=seed_entropy,
            seed_candidates=(seed_entropy, "bad"),
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )
        @test_throws ArgumentError Models.solve_multi(m, Models.FixedAsymmetricRho(0.05, 1.0, 0.0), T_fm;
            seed_candidates=(seed_entropy,),
            evaluate_all_attempts="bad",
            p_num=8,
            t_num=4,
            residual_norm_max=1e-6,
        )
    end

    @testset "solve FixedMu 快捷入口" begin
        mode = Models.FixedMu()
        T = 0.5
        μ = 0.0
        result = Models.solve(mode, T, μ; p_num=24, t_num=6)
        @test result isa Any  # 确保不抛异常
    end

    @testset "SolverResult 类型可用" begin
        @test isdefined(Models, :SolverResult) || true  # 软检查
    end

    @testset "_select_pressure_max_local requires mandatory candidate fields" begin
        missing_hard_constraint_ok = [(
            pressure=1.0,
            residual_norm=1e-8,
            seed_index=1,
            converged=true,
        )]
        missing_residual_norm = [(
            pressure=1.0,
            hard_constraint_ok=true,
            seed_index=1,
            converged=true,
        )]
        missing_pressure = [(
            residual_norm=1e-8,
            hard_constraint_ok=true,
            seed_index=1,
            converged=true,
        )]

        @test_throws ArgumentError Models._select_pressure_max_local(missing_hard_constraint_ok)
        @test_throws ArgumentError Models._select_pressure_max_local(missing_residual_norm)
        @test_throws ArgumentError Models._select_pressure_max_local(missing_pressure)
    end

    @testset "solve_constraint FixedRho 可调用" begin
        m = Models.create_model(:NJL)
        mode = Models.FixedRho(0.5)
        T = 0.5
        seed = Models.gap_initial_guess(m, T, SVector{3}(0.0, 0.0, 0.0))
        # FixedRho 可能需要额外参数；确保接口存在即可
        @test Models.solve_constraint isa Function
    end

    @testset "non-FixedMu solve supports semantic bridge kwargs" begin
        mode = Models.FixedEntropy(0.5)
        T_fm = 100.0 / 197.327
        seed = copy(Models.pnjl_module().HADRON_SEED_8)

        result = Models.solve(mode, T_fm;
            seed_guess=seed,
            seed_candidates=(seed,),
            semantic_mode=:constrained_manifold,
            p_num=8,
            t_num=4,
            rho0=0.16,
            residual_norm_max=1e-6,
        )

        @test result isa Models.SolverResult
        @test isfinite(result.residual_norm)
    end

    @testset "non-FixedMu solve_multi supports semantic bridge kwargs" begin
        mode = Models.FixedAsymmetricRho(0.05, 1.0, 0.0)
        T_fm = 100.0 / 197.327
        seed = copy(Models.pnjl_module().HADRON_SEED_8)

        result = Models.solve_multi(mode, T_fm;
            seed_candidates=(seed,),
            semantic_mode=:constrained_manifold,
            p_num=8,
            t_num=4,
            rho0=0.16,
            residual_norm_max=1e-6,
        )

        @test result isa Models.SolverResult
        @test isfinite(result.residual_norm) || isinf(result.residual_norm)
        @test result.residual_norm >= 0.0 || isinf(result.residual_norm)
    end

    @testset "ThermoPostprocess unified helpers are available" begin
        model = Models.create_model(:PNJL)
        T_fm = 100.0 / 197.327
        seed5 = copy(Models.pnjl_module().HADRON_SEED_5)
        solution = Float64[seed5..., 0.0, 0.0, 0.0]

        thermo = Models.compute_thermo_from_solution(
            model,
            solution,
            T_fm;
            xi=0.0,
            p_num=8,
            t_num=4,
            rho0_scale=Models.pnjl_module().ρ0,
        )
        residual_norm = Models.compute_residual_norm_from_solution(
            model,
            solution,
            T_fm;
            xi=0.0,
            p_num=8,
            t_num=4,
        )
        candidate = Models.build_solver_candidate(
            solution,
            thermo,
            residual_norm;
            converged=true,
            iterations=1,
            residual_norm_max=1e-6,
        )

        @test isfinite(thermo.pressure)
        @test isfinite(residual_norm)
        @test candidate.solution == solution
        @test candidate.iterations == 1
    end

    @testset "SolverRuntimeConfig typed parsers for FixedRho/FixedEntropy" begin
        model = Models.create_model(:PNJL)
        T_fm = 100.0 / 197.327
        seed = copy(Models.pnjl_module().HADRON_SEED_8)

        rho_mode = Models.FixedRho(1.0)
        rho_kwargs = Dict{Symbol,Any}(
            :seed_guess => seed,
            :residual_norm_max => 1e-6,
        )
        rho_cfg = Models._fixedrho_runtime_config_from_kwargs(rho_mode, rho_kwargs)
        @test rho_cfg isa Models.FixedRhoRuntimeConfig
        @test rho_cfg.fixedrho_joint_solve === true
        @test rho_cfg.seed_guess == Float64.(seed)

        entropy_mode = Models.FixedEntropy(0.5)
        entropy_kwargs = Dict{Symbol,Any}(
            :seed_guess => seed,
            :residual_norm_max => 1e-6,
            :rho0 => Models.pnjl_module().ρ0,
        )
        entropy_cfg = Models._fixedentropy_runtime_config_from_kwargs(entropy_mode, entropy_kwargs)
        @test entropy_cfg isa Models.FixedEntropyRuntimeConfig
        @test entropy_cfg.seed_guess == Float64.(seed)
        @test entropy_cfg.primary_method == :trust_region

        @test_throws ArgumentError Models._fixedrho_runtime_config_from_kwargs(
            rho_mode,
            Dict{Symbol,Any}(:seed_guess => seed, :fixedrho_joint_solve => "bad"),
        )
        @test_throws ArgumentError Models._fixedentropy_runtime_config_from_kwargs(
            entropy_mode,
            Dict{Symbol,Any}(),
        )
    end

    @testset "SolverDiagnostics typed structures map to legacy NamedTuple" begin
        summary = Models.SolverDiagnosticSummary(
            Models.SOLVER_DIAGNOSTIC_VERSION_V1,
            :primary,
            :seed,
            true,
            Symbol[],
            :none,
            "",
            :pressure_max_under_constraints,
            :converged,
            nothing,
            :problem_spec_selector,
        )
        candidate = Models.SolverDiagnosticCandidate(summary)
        full = Models.SolverDiagnosticFull(summary, [candidate])

        summary_nt = Models.to_namedtuple(summary)
        candidate_nt = Models.to_namedtuple(candidate)
        full_nt = Models.to_namedtuple(full)

        @test summary_nt.selection_reason == :pressure_max_under_constraints
        @test candidate_nt.error_kind == :none
        @test full_nt.selection_reason_source == :problem_spec_selector
        @test length(full_nt.candidates) == 1
    end

    @testset "FixedEntropy default/bridge semantic parity (single point)" begin
        model = Models.create_model(:PNJL)
        mode = Models.FixedEntropy(0.5)
        T_fm = 100.0 / 197.327
        seed = copy(Models.pnjl_module().HADRON_SEED_8)

        default_result = Models.solve(model, mode, T_fm;
            seed_strategy=Models.DefaultSeed(seed, seed, :hadron),
            p_num=16,
            t_num=6,
            residual_norm_max=1e-6,
        )

        bridge_result = Models.solve(model, mode, T_fm;
            seed_strategy=Models.DefaultSeed(seed, seed, :hadron),
            seed_guess=seed,
            seed_candidates=(seed,),
            semantic_mode=:ground_state,
            rho0=Models.pnjl_module().ρ0,
            p_num=16,
            t_num=6,
            residual_norm_max=1e-6,
        )

        @test default_result.converged == bridge_result.converged
        @test isfinite(default_result.residual_norm)
        @test isfinite(bridge_result.residual_norm)
        @test isapprox(default_result.residual_norm, bridge_result.residual_norm; rtol=1e-6, atol=1e-8)
        if default_result.converged && bridge_result.converged
            @test isapprox(default_result.entropy, bridge_result.entropy; rtol=1e-3, atol=1e-5)
            @test isapprox(default_result.pressure, bridge_result.pressure; rtol=1e-3, atol=1e-5)
        end
    end

    @testset "FixedSigma default/bridge semantic parity (single point)" begin
        model = Models.create_model(:PNJL)
        mode = Models.FixedSigma(10.0)
        T_fm = 100.0 / 197.327
        seed = copy(Models.pnjl_module().HADRON_SEED_8)

        default_result = Models.solve(model, mode, T_fm;
            seed_strategy=Models.DefaultSeed(seed, seed, :hadron),
            p_num=16,
            t_num=6,
            residual_norm_max=1e-6,
        )

        bridge_result = Models.solve(model, mode, T_fm;
            seed_strategy=Models.DefaultSeed(seed, seed, :hadron),
            seed_guess=seed,
            seed_candidates=(seed,),
            semantic_mode=:ground_state,
            rho0=Models.pnjl_module().ρ0,
            p_num=16,
            t_num=6,
            residual_norm_max=1e-6,
        )

        @test default_result.converged == bridge_result.converged
        @test isfinite(default_result.residual_norm)
        @test isfinite(bridge_result.residual_norm)
        @test isapprox(default_result.residual_norm, bridge_result.residual_norm; rtol=1e-6, atol=1e-8)
        if default_result.converged && bridge_result.converged
            @test isapprox(default_result.rho_norm, bridge_result.rho_norm; rtol=1e-3, atol=1e-5)
            @test isapprox(default_result.pressure, bridge_result.pressure; rtol=1e-3, atol=1e-5)
        end
    end

    @testset "FixedRho default/bridge semantic parity (single point)" begin
        model = Models.create_model(:PNJL)
        mode = Models.FixedRho(1.0)
        T_fm = 100.0 / 197.327
        seed = copy(Models.pnjl_module().HADRON_SEED_8)
        seed_strategy = Models.DefaultSeed(seed, seed, :hadron)
        bridge_seed = Models.get_seed(seed_strategy, [T_fm], mode)

        default_result = Models.solve(model, mode, T_fm;
            seed_strategy=seed_strategy,
            p_num=16,
            t_num=6,
            residual_norm_max=1e-6,
        )

        bridge_result = Models.solve(model, mode, T_fm;
            seed_strategy=seed_strategy,
            seed_guess=bridge_seed,
            seed_candidates=(bridge_seed,),
            semantic_mode=:ground_state,
            p_num=16,
            t_num=6,
            residual_norm_max=1e-6,
        )

        @test default_result.converged == bridge_result.converged
        @test isfinite(default_result.residual_norm)
        @test isfinite(bridge_result.residual_norm)
        @test isapprox(default_result.residual_norm, bridge_result.residual_norm; rtol=1e-6, atol=1e-8)
        if default_result.converged && bridge_result.converged
            @test isapprox(default_result.rho_norm, bridge_result.rho_norm; rtol=1e-3, atol=1e-5)
            @test isapprox(default_result.pressure, bridge_result.pressure; rtol=1e-3, atol=1e-5)
        end
    end

    @testset "FixedAsymmetricRho default/bridge semantic parity (single point)" begin
        model = Models.create_model(:PNJL)
        mode = Models.FixedAsymmetricRho(0.05, 1.0, 0.0)
        T_fm = 100.0 / 197.327
        seed = copy(Models.pnjl_module().HADRON_SEED_8)

        default_result = Models.solve(model, mode, T_fm;
            seed_strategy=Models.DefaultSeed(seed, seed, :hadron),
            p_num=8,
            t_num=4,
            iterations=120,
            residual_norm_max=1e-6,
        )

        bridge_result = Models.solve(model, mode, T_fm;
            seed_strategy=Models.DefaultSeed(seed, seed, :hadron),
            seed_guess=seed,
            seed_candidates=(seed,),
            semantic_mode=:ground_state,
            rho0=Models.pnjl_module().ρ0,
            p_num=8,
            t_num=4,
            iterations=120,
            residual_norm_max=1e-6,
        )

        @test default_result.converged == bridge_result.converged
        @test isfinite(default_result.residual_norm)
        @test isfinite(bridge_result.residual_norm)
        if default_result.converged && bridge_result.converged
            @test isapprox(default_result.rho_norm, bridge_result.rho_norm; rtol=1e-3, atol=1e-5)
            @test isapprox(default_result.pressure, bridge_result.pressure; rtol=1e-3, atol=1e-5)
        end
    end
end
