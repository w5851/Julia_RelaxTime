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

    # --- solve_fixedmu_constraint ---
    @testset "solve_fixedmu_constraint NJL" begin
        m = Models.create_model(:NJL)
        T = 0.5
        seed = Float64.(Models.gap_initial_guess(m, T, SVector{3}(0.0, 0.0, 0.0)))
        result = Models.solve_fixedmu_constraint(m, T, 0.0; seed_guess=seed, p_num=24, t_num=6)
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
    end

    @testset "solve_fixedmu_constraint 结果一致性 P=-Ω" begin
        m = Models.create_model(:NJL)
        T = 0.5
        seed = Float64.(Models.gap_initial_guess(m, T, SVector{3}(0.0, 0.0, 0.0)))
        result = Models.solve_fixedmu_constraint(m, T, 0.0; seed_guess=seed, p_num=24, t_num=6)
        @test result.pressure ≈ -result.omega rtol=1e-10
    end

    # --- solve_fixedrho_constraint 接口存在 ---
    @testset "solve_fixedrho_constraint 函数存在" begin
        @test isdefined(Models, :solve_fixedrho_constraint)
        @test Models.solve_fixedrho_constraint isa Function
    end

    @testset "solve_fixedentropy_constraint 函数存在" begin
        @test isdefined(Models, :solve_fixedentropy_constraint)
    end

    @testset "solve_fixedsigma_constraint 函数存在" begin
        @test isdefined(Models, :solve_fixedsigma_constraint)
    end

    @testset "solve_fixedasymrho_constraint 函数存在" begin
        @test isdefined(Models, :solve_fixedasymrho_constraint)
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
end
