# gap_solver.jl 单元测试
#
# 测试内容：
# 1. NLsolveGapSolver 构造
# 2. gap_state_dim 分发
# 3. gap_initial_guess 默认种子
# 4. gap_residual 残差计算
# 5. solve_gap 基本收敛

using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

# ============================================================================

@testset "gap_solver" begin

    # --- NLsolveGapSolver 构造 ---
    @testset "NLsolveGapSolver 默认" begin
        s = Models.NLsolveGapSolver()
        @test s.method === :newton
        @test s.xtol > 0
        @test s.ftol > 0
    end

    @testset "NLsolveGapSolver kwarg 重载" begin
        s = Models.NLsolveGapSolver(method=:trust_region, xtol=1e-6, ftol=1e-6, jacobian=:finite)
        @test s.method === :trust_region
        @test s.xtol ≈ 1e-6
        @test s.jacobian === :finite
    end

    # --- gap_state_dim ---
    @testset "gap_state_dim NJL=3" begin
        m = Models.create_model(:NJL)
        @test Models.gap_state_dim(m) == 3
    end

    @testset "gap_state_dim NJL2=2" begin
        m = Models.create_model(:NJL2)
        @test Models.gap_state_dim(m) == 2
    end

    @testset "gap_state_dim PNJL=5" begin
        Models.pnjl_module()
        m = Models.create_model(:PNJL)
        @test Models.gap_state_dim(m) == 5
    end

    # --- gap_initial_guess ---
    @testset "gap_initial_guess NJL 长度=3" begin
        m = Models.create_model(:NJL)
        g = Models.gap_initial_guess(m, 0.5, SVector{3}(0.0, 0.0, 0.0))
        @test length(g) == 3
        @test all(isfinite.(g))
    end

    @testset "gap_initial_guess PNJL 长度=5" begin
        m = Models.create_model(:PNJL)
        g = Models.gap_initial_guess(m, 0.5, SVector{3}(0.0, 0.0, 0.0))
        @test length(g) == 5
        @test all(isfinite.(g))
    end

    # --- gap_residual ---
    @testset "gap_residual NJL 返回 SVector{3}" begin
        m = Models.create_model(:NJL)
        x = SVector{3}(-1.84, -1.84, -2.23)
        T = 0.5
        μ = SVector{3}(0.0, 0.0, 0.0)
        r = Models.gap_residual(m, x, T, μ; p_num=24, t_num=6, xi=0.0)
        @test r isa SVector{3}
        @test all(isfinite.(r))
    end

    @testset "gap_residual PNJL 返回 SVector{5}" begin
        m = Models.create_model(:PNJL)
        x = SVector{5}(-1.84, -1.84, -2.23, 0.5, 0.5)
        T = 0.5
        μ = SVector{3}(0.0, 0.0, 0.0)
        r = Models.gap_residual(m, x, T, μ; p_num=24, t_num=6, xi=0.0)
        @test r isa SVector{5}
        @test all(isfinite.(r))
    end

    # --- solve_gap 基本收敛 ---
    @testset "solve_gap NJL 收敛" begin
        m = Models.create_model(:NJL)
        T = 0.5
        μ = SVector{3}(0.0, 0.0, 0.0)
        st = Models.solve_gap(m, T, μ; p_num=24, t_num=6, xi=0.0)
        @test st isa Models.MeanFieldState
        @test all(isfinite.(Models.state_vector(st)))
    end
end
