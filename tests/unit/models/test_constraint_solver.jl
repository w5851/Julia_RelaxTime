# constraint_solver.jl 单元测试
#
# 测试内容：
# 1. default_mu0_from_seed / default_muvec0_from_seed
# 2. solve_fixedmu_constraint 基本调用
# 3. solve_fixedrho_constraint 基本调用
# 4. 返回值结构

using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Models.pnjl_module()

# ============================================================================

@testset "constraint_solver" begin

    # --- 辅助函数 ---
    @testset "default_mu0_from_seed" begin
        # 长度 < 6 → 默认 0.2
        @test Models.default_mu0_from_seed([1.0, 2.0, 3.0]) ≈ 0.2
        # 长度 >= 6 → 取 seed[6]
        seed6 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.5]
        @test Models.default_mu0_from_seed(seed6) ≈ 0.5
        # 长度 >= 8 → 取 mean(seed[6:8])
        seed8 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.6, 0.9]
        @test Models.default_mu0_from_seed(seed8) ≈ 0.6
    end

    @testset "default_muvec0_from_seed" begin
        seed8 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.6, 0.9]
        μv = Models.default_muvec0_from_seed(seed8)
        @test length(μv) == 3
        @test μv == [0.3, 0.6, 0.9]

        seed3 = [1.0, 2.0, 3.0]
        μv3 = Models.default_muvec0_from_seed(seed3)
        @test length(μv3) == 3
        @test μv3[1] ≈ 0.2
    end

    @testset "solution pack/unpack helpers" begin
        x_state = SVector{5}(-1.5, -1.5, -2.1, 0.2, 0.2)
        mu_vec = SVector{3}(1.5, 1.4, 1.3)

        solution = Models._pack_solution(x_state, mu_vec)
        @test length(solution) == 8

        x2, mu2 = Models._unpack_solution(solution; state_n=5, mu_n=3)
        @test x2 == x_state
        @test mu2 == mu_vec

        empty_candidate = Models._empty_candidate(state_n=5, mu_n=3, residual_norm_max=1e-6)
        @test length(empty_candidate.solution) == 8
        @test all(isnan, empty_candidate.solution)
        @test empty_candidate.failed_constraints == [:solver_failed]
        @test !empty_candidate.converged

        fixedmu_empty = Models._empty_candidate(state_n=5, mu_n=3, solution_n=5, residual_norm_max=1e-6)
        @test length(fixedmu_empty.solution) == 5
        @test all(isnan, fixedmu_empty.solution)
    end

    @testset "solve_constraint(FixedMu) NJL" begin
        m = Models.create_model(:NJL)
        T = 0.5
        seed = Float64.(Models.gap_initial_guess(m, T, SVector{3}(0.0, 0.0, 0.0)))
        result = Models.solve_constraint(m, Models.FixedMu(), T; μ_fm=0.0, seed_guess=seed, p_num=24, t_num=6)
        @test result isa NamedTuple
        @test haskey(result, :converged)
        @test haskey(result, :pressure)
        @test haskey(result, :rho_norm)
        @test haskey(result, :entropy)
        @test haskey(result, :energy)
        @test haskey(result, :masses)
        @test haskey(result, :residual_norm)
        @test isfinite(result.pressure)
        @test isfinite(result.omega)
        @test result.pressure ≈ -result.omega rtol=1e-10
    end

    @testset "solve_constraint(FixedMu) 默认多初值候选池" begin
        m = Models.create_model(:NJL)
        T = 0.5
        seed = Float64.(Models.gap_initial_guess(m, T, SVector{3}(0.0, 0.0, 0.0)))
        result = Models.solve_constraint(m, Models.FixedMu(), T; μ_fm=0.0, seed_guess=seed, p_num=24, t_num=6)
        @test haskey(result, :candidate_count)
        @test result.candidate_count >= 1
        @test haskey(result, :selection_reason)
        @test result.selection_reason in (:pressure_max_under_constraints, :no_candidate_passed_constraints)
    end

    @testset "solve_constraint(FixedRho) PNJL 归一化密度收敛" begin
        m = Models.create_model(:PNJL)
        mode = Models.FixedRho(0.2)
        T_fm = 150.0 / Main.Constants_PNJL.ħc_MeV_fm
        seed = [-1.5, -1.5, -2.1, 0.2, 0.2, 1.5, 1.5, 1.5]

        result = Models.solve_constraint(
            m,
            mode,
            T_fm;
            seed_guess=seed,
            p_num=12,
            t_num=4,
            residual_norm_max=1e-6,
        )

        @test result.converged
        @test isfinite(result.rho_norm)
        @test result.rho_norm ≈ 0.2 atol=1e-3
    end

    @testset "legacy fixed-* APIs removed" begin
        @test !isdefined(Models, :solve_fixedmu_constraint)
        @test !isdefined(Models, :solve_fixedrho_constraint)
        @test !isdefined(Models, :solve_fixedentropy_constraint)
        @test !isdefined(Models, :solve_fixedsigma_constraint)
        @test !isdefined(Models, :solve_fixedasymrho_constraint)
    end

    @testset "硬约束评估可解耦扩展" begin
        rules = [
            c -> (c.mass_min > 0.0, :mass_positive),
            c -> (c.user_tag != :ban, :user_custom),
        ]

        passed, failed = Models.evaluate_hard_constraints((mass_min=1.2, user_tag=:ok), rules)
        @test passed
        @test isempty(failed)

        passed2, failed2 = Models.evaluate_hard_constraints((mass_min=1.2, user_tag=:ban), rules)
        @test !passed2
        @test failed2 == [:user_custom]
    end

    @testset "pressure_max 仅在约束通过池内选优" begin
        c1 = (pressure=10.0, residual_norm=1e-8, hard_constraint_ok=true, failed_constraints=Symbol[], converged=true)
        c2 = (pressure=12.0, residual_norm=1e-8, hard_constraint_ok=false, failed_constraints=[:mass_nonpositive], converged=true)

        selected = Models.select_pressure_max_candidate([c1, c2])
        @test selected.selection_reason == :pressure_max_under_constraints
        @test selected.selected_index == 1
        @test selected.selected_candidate.pressure == 10.0
    end

    @testset "无约束通过候选时回退并标注原因" begin
        c1 = (pressure=8.0, residual_norm=1e-6, hard_constraint_ok=false, failed_constraints=[:mass_nonpositive], converged=true)
        c2 = (pressure=9.0, residual_norm=1e-6, hard_constraint_ok=false, failed_constraints=[:phi_out_of_range], converged=true)

        selected = Models.select_pressure_max_candidate([c1, c2])
        @test selected.selection_reason == :no_candidate_passed_constraints
        @test selected.selected_index == 2
        @test !selected.selected_candidate.hard_constraint_ok
    end

    @testset "residual spine guard" begin
        solver_path = joinpath(PROJECT_ROOT, "src", "models", "solver", "ConstraintSolver.jl")
        source = read(solver_path, String)
        @test occursin("_gap_norm_from_state", source)
        @test !occursin("gap_residual(model", source)
    end
end
