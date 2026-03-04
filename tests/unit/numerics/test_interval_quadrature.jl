# IntervalQuadratureStrategies 单元测试
#
# 测试内容：
# 1. IntegrationStrategy 枚举
# 2. clustered_gl_nodes 节点生成
# 3. power_left_nodes / power_right_nodes 变换
# 4. de_nodes (tanh-sinh) 变换
# 5. hybrid_nodes 分发
# 6. 积分精度验证

using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

# IntervalQuadratureStrategies 依赖 GaussLegendre 已被 include
include(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
include(joinpath(PROJECT_ROOT, "src", "integration", "IntervalQuadratureStrategies.jl"))

# ============================================================================
# 枚举与常量
# ============================================================================

@testset "IntervalQuadratureStrategies" begin

    @testset "IntegrationStrategy 枚举" begin
        @test STRATEGY_QUADGK isa IntegrationStrategy
        @test STRATEGY_INTERVAL_GL isa IntegrationStrategy
        @test STRATEGY_CLUSTER_GL isa IntegrationStrategy
        @test STRATEGY_HYBRID isa IntegrationStrategy
        @test DEFAULT_STRATEGY == STRATEGY_HYBRID
    end

    @testset "默认参数常量" begin
        @test DEFAULT_CLUSTER_BETA > 0
        @test DEFAULT_CLUSTER_N > 0
        @test DEFAULT_POWER_ALPHA > 0
        @test DEFAULT_DE_H > 0
    end

    # ============================================================================
    # clustered_gl_nodes
    # ============================================================================

    @testset "clustered_gl_nodes" begin
        a, b = 0.0, 1.0
        xs, wx = clustered_gl_nodes(a, b, 32)
        @test length(xs) == 32
        @test length(wx) == 32
        @test all(a .<= xs .<= b)
        @test all(isfinite.(wx))

        # 积分 ∫₀¹ 1 dx ≈ 1 (tanh 聚簇在高 beta 下有微小精度损失)
        @test sum(wx) ≈ 1.0 rtol=1e-3
    end

    # ============================================================================
    # power_left_nodes
    # ============================================================================

    @testset "power_left_nodes" begin
        a, b = 0.0, 1.0
        xs, wx = power_left_nodes(a, b, 32)
        @test length(xs) == 32
        @test all(a .<= xs .<= b)
        # 节点应聚簇于左端（更多点在 a 附近）
        mid = (a + b) / 2
        n_left = count(x -> x < mid, xs)
        @test n_left > 16  # 大多数节点在左半

        # 积分 ∫₀¹ x dx = 0.5
        @test sum(wx .* xs) ≈ 0.5 rtol=1e-4
    end

    # ============================================================================
    # power_right_nodes
    # ============================================================================

    @testset "power_right_nodes" begin
        a, b = 0.0, 1.0
        xs, wx = power_right_nodes(a, b, 32)
        @test length(xs) == 32
        @test all(a .<= xs .<= b)
        # 节点应聚簇于右端
        mid = (a + b) / 2
        n_right = count(x -> x > mid, xs)
        @test n_right > 16

        # 积分 ∫₀¹ x dx = 0.5
        @test sum(wx .* xs) ≈ 0.5 rtol=1e-4
    end

    # ============================================================================
    # de_nodes (tanh-sinh 变换)
    # ============================================================================

    @testset "de_nodes" begin
        a, b = 0.0, 1.0
        xs, wx = de_nodes(a, b, 32)
        @test length(xs) == length(wx)
        @test all(isfinite.(xs))
        @test all(isfinite.(wx))

        # 积分 ∫₀¹ 1 dx ≈ 1
        @test sum(wx) ≈ 1.0 rtol=1e-2
    end

    # ============================================================================
    # hybrid_nodes 分发
    # ============================================================================

    @testset "hybrid_nodes 分发" begin
        a, b = 0.0, 1.0
        n = 32

        @testset "SING_LEFT → power_left" begin
            xs_h, wx_h = hybrid_nodes(a, b, n, SING_LEFT)
            xs_p, wx_p = power_left_nodes(a, b, n)
            @test xs_h ≈ xs_p
            @test wx_h ≈ wx_p
        end

        @testset "SING_RIGHT → power_right" begin
            xs_h, wx_h = hybrid_nodes(a, b, n, SING_RIGHT)
            xs_p, wx_p = power_right_nodes(a, b, n)
            @test xs_h ≈ xs_p
            @test wx_h ≈ wx_p
        end

        @testset "SING_BOTH → de" begin
            xs_h, wx_h = hybrid_nodes(a, b, n, SING_BOTH)
            xs_d, wx_d = de_nodes(a, b, n)
            @test xs_h ≈ xs_d
            @test wx_h ≈ wx_d
        end
    end

    # ============================================================================
    # 精度验证：光滑函数
    # ============================================================================

    @testset "精度验证：∫₀¹ x² dx = 1/3" begin
        a, b = 0.0, 1.0
        f(x) = x^2

        xs, wx = clustered_gl_nodes(a, b, 32)
        @test sum(wx .* f.(xs)) ≈ 1/3 rtol=1e-2  # tanh 聚簇对光滑函数精度 ~0.1%

        xs2, wx2 = power_left_nodes(a, b, 32)
        @test sum(wx2 .* f.(xs2)) ≈ 1/3 rtol=1e-4
    end

    # ============================================================================
    # SingularityPosition 枚举
    # ============================================================================

    @testset "SingularityPosition" begin
        @test SING_NONE isa SingularityPosition
        @test SING_LEFT isa SingularityPosition
        @test SING_RIGHT isa SingularityPosition
        @test SING_BOTH isa SingularityPosition
    end
end
