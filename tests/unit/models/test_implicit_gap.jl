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
    @testset "legacy factories are qualified compat-only" begin
        exported = names(Models)
        @test !(:create_implicit_gap_solver in exported)
        @test !(:create_pnjl_implicit_solver in exported)
    end

    # --- compat create_implicit_gap_solver NJL2 ---
    @testset "compat create_implicit_gap_solver NJL2" begin
        m = Models.create_model(:NJL2)
        igf = Models.create_implicit_gap_solver(m; p_num=24, t_num=6)
        @test igf isa Any  # ImplicitFunction 类型
    end

    # --- compat create_implicit_gap_solver NJL ---
    @testset "compat create_implicit_gap_solver NJL" begin
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

    @testset "derive_vec and derive_named parity" begin
        model = Models.create_model(:PNJL)
        theta_vec = [0.5, 0.1]
        theta_named = (T_fm=0.5, μ_fm=0.1)

        @test isdefined(Models, :derive_vec)
        @test isdefined(Models, :derive_named)

        old_result = Models.solve_pnjl_with_derivatives(theta_vec[1], theta_vec[2]; order=1, p_num=24, t_num=6)
        auto_result = Models.solve_pnjl_with_derivatives(theta_vec[1], theta_vec[2]; order=1, p_num=8, t_num=4, derivative_backend=:auto)
        td_result = Models.solve_pnjl_with_derivatives(theta_vec[1], theta_vec[2]; order=1, p_num=8, t_num=4, derivative_backend=:taylordiff)
        vec_result = Models.derive_vec(model, theta_vec; order=1, p_num=24, t_num=6)
        named_result = Models.derive_named(model, theta_named; order=1, p_num=24, t_num=6)

        @test all(isapprox.(vec_result.x, old_result.x; rtol=1e-7, atol=1e-9))
        @test all(isapprox.(vec_result.dx_dT, old_result.dx_dT; rtol=1e-6, atol=1e-8))
        @test all(isapprox.(vec_result.dx_dμ, old_result.dx_dμ; rtol=1e-6, atol=1e-8))
        @test all(isapprox.(td_result.x, auto_result.x; rtol=1e-12, atol=1e-12))
        @test all(isapprox.(td_result.dx_dT, auto_result.dx_dT; rtol=1e-12, atol=1e-12))
        @test all(isapprox.(td_result.dx_dμ, auto_result.dx_dμ; rtol=1e-12, atol=1e-12))
        @test_throws ArgumentError Models.solve_pnjl_with_derivatives(theta_vec[1], theta_vec[2]; order=1, p_num=8, t_num=4, derivative_backend=:forwarddiff)

        @test all(isapprox.(named_result.x, vec_result.x; rtol=1e-12, atol=1e-12))
        @test all(isapprox.(named_result.dx_dT, vec_result.dx_dT; rtol=1e-12, atol=1e-12))
        @test all(isapprox.(named_result.dx_dμ, vec_result.dx_dμ; rtol=1e-12, atol=1e-12))
    end
end
