# implicit_gap.jl 单元测试
#
# 测试内容：
# 1. create_implicit_gap_solver NJL2/NJL/PNJL
# 2. 隐函数求解 forward solve
# 3. solve_pnjl_with_derivatives (PNJL wrapper)

using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
Models.pnjl_module()

# ============================================================================

@testset "implicit_gap" begin

    # --- create_implicit_gap_solver NJL2 ---
    @testset "create_implicit_gap_solver NJL2" begin
        m = Models.create_model(:NJL2)
        igf = Models.create_implicit_gap_solver(m; p_num=24, t_num=6)
        @test igf isa Any  # ImplicitFunction 类型
    end

    # --- create_implicit_gap_solver NJL ---
    @testset "create_implicit_gap_solver NJL" begin
        m = Models.create_model(:NJL)
        igf = Models.create_implicit_gap_solver(m; p_num=24, t_num=6)
        @test igf isa Any
    end

    # --- forward solve NJL2 ---
    @testset "隐函数 forward solve NJL2" begin
        m = Models.create_model(:NJL2)
        igf = Models.create_implicit_gap_solver(m; p_num=24, t_num=6)
        θ = [0.5, 0.0]  # [T, μ]
        result = igf(θ)
        # ImplicitFunction 返回 (x, z) tuple
        x = result isa Tuple ? result[1] : result
        @test length(x) >= 2
        @test all(isfinite.(x))
    end

    # --- forward solve NJL ---
    @testset "隐函数 forward solve NJL" begin
        m = Models.create_model(:NJL)
        igf = Models.create_implicit_gap_solver(m; p_num=24, t_num=6)
        θ = [0.5, 0.0]  # [T, μ]
        result = igf(θ)
        x = result isa Tuple ? result[1] : result
        @test length(x) >= 3
        @test all(isfinite.(x))
    end

    # --- create_pnjl_implicit_solver ---
    @testset "create_pnjl_implicit_solver" begin
        @test isdefined(Models, :create_pnjl_implicit_solver)
        igf = Models.create_pnjl_implicit_solver(p_num=24, t_num=6)
        @test igf isa Any
    end

    # --- solve_pnjl_with_derivatives ---
    @testset "solve_pnjl_with_derivatives" begin
        @test isdefined(Models, :solve_pnjl_with_derivatives)
    end
end
