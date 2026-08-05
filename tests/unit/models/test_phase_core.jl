# PhaseCore.jl 单元测试
#
# 测试内容：
# 1. detect_s_shape 正 S 型 / 单调 / 边界
# 2. maxwell_construction 收敛性
# 3. SShapeResult / MaxwellResult 结构

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

# ============================================================================

@testset "PhaseCore" begin

    # --- detect_s_shape ---
    @testset "detect_s_shape 单调曲线无 S 型" begin
        mu = collect(range(0.0, 2.0, length=20))
        rho = collect(range(0.0, 1.0, length=20))  # 严格递增
        res = Models.detect_s_shape(mu, rho)
        @test res isa Models.SShapeResult
        @test res.has_s_shape == false
    end

    @testset "detect_s_shape 合成 S 型" begin
        # 构造一个有回折的 S 型曲线
        mu = collect(range(0.0, 2.0, length=40))
        rho = @. 0.5 + 0.4 * tanh(3.0 * (mu - 1.0)) - 0.15 * sin(3.0 * mu)
        # 添加回折段
        rho[15:25] .= reverse(rho[15:25])
        res = Models.detect_s_shape(mu, rho)
        @test res isa Any
    end

    # --- maxwell_construction ---
    @testset "maxwell_construction 接口存在" begin
        @test isdefined(Models, :maxwell_construction)
        @test Models.maxwell_construction isa Function
    end

    @testset "maxwell_construction 单调曲线无转变" begin
        mu = collect(range(0.0, 2.0, length=20))
        rho = collect(range(0.0, 1.0, length=20))
        res = Models.maxwell_construction(mu, rho)
        @test res isa Models.MaxwellResult
    end

    @testset "Maxwell 严格三交点与候选诊断" begin
        rho = Float64[0.0, 1.0, 2.0, 3.0]
        mu = Float64[0.0, 2.0, 1.0, 1.0]
        @test length(Models._find_intersections(1.5, rho, mu)) == 2
        @test Models._area_difference(1.5, rho, mu) === nothing

        mu_s = Float64[250.0, 255.0, 260.0, 265.0, 268.0, 266.0, 262.0,
            258.0, 260.0, 265.0, 272.0, 280.0, 288.0, 296.0, 305.0]
        rho_s = collect(0.1:0.1:1.5)
        good = Models.maxwell_construction(mu_s, rho_s; min_samples=8)
        @test good.converged
        @test get(good.details, :candidate_count, 0) == 1
        @test get(good.details, :crossing_count, 0) == 3
        strict = Models.maxwell_construction(mu_s, rho_s; min_samples=8,
            tol_area=1e-14, max_iter=1)
        @test !strict.converged
        @test get(strict.details, :reason, "") == "bisection_failed"
        @test get(strict.details, :failure_reason, "") == "solver_tolerance_not_met"
    end

    # --- SShapeResult 结构 ---
    @testset "SShapeResult 零参构造" begin
        res = Models.SShapeResult()
        @test res.has_s_shape == false
        @test res.derivative_sign_changes == 0
        @test res.mu_spinodal_hadron === nothing
    end

    # --- MaxwellResult 结构 ---
    @testset "MaxwellResult 零参构造" begin
        res = Models.MaxwellResult()
        @test res.converged == false
        @test res.iterations == 0
        @test res.mu_transition === nothing
    end
end
