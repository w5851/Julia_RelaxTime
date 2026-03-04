# AdaptiveRhoRefinement 模块单元测试

using Test

const PROJECT_ROOT_ARR = normpath(joinpath(@__DIR__, "..", "..", ".."))

const _ARR_PATH = normpath(joinpath(PROJECT_ROOT_ARR, "src", "models", "phase", "AdaptiveRhoRefinement.jl"))
if !isdefined(Main, :AdaptiveRhoRefinement)
    Base.include(Main, _ARR_PATH)
end

using Main.AdaptiveRhoRefinement: AdaptiveRhoConfig, suggest_refinement_points, merge_rho_values

@testset "AdaptiveRhoRefinement" begin
    @testset "AdaptiveRhoConfig 默认参数" begin
        cfg = AdaptiveRhoConfig()
        @test cfg.slope_tol > 0
        @test cfg.min_gap > 0
        @test cfg.max_points > 0
    end

    @testset "merge_rho_values 合并并排序" begin
        existing = [0.0, 0.5, 1.0]
        additions = [0.25, 0.75]
        merged = merge_rho_values(existing, additions)
        @test issorted(merged)
        @test length(merged) == 5
        @test 0.25 in merged
    end

    @testset "merge_rho_values 去重" begin
        existing = [0.0, 0.5, 1.0]
        additions = [0.5, 0.75]
        merged = merge_rho_values(existing, additions)
        @test count(==(0.5), merged) == 1
    end

    @testset "suggest_refinement_points 返回合理建议" begin
        rho_vals = collect(range(0.0, 1.0, length=10))
        # 创建一个有急剧变化的 mu 曲线
        mu_vals = [r < 0.5 ? 0.3 : 0.1 for r in rho_vals]
        points = suggest_refinement_points(rho_vals, mu_vals)
        @test points isa AbstractVector
        @test all(isfinite, points)
    end
end
