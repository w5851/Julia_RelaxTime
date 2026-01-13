"""
测试极化函数缓存功能

测试内容：
1. 基本缓存功能：相同参数重复调用
2. 浮点数容差：微小差异参数的缓存命中
3. 不同参数：确保不同参数不会错误命中
4. 缓存统计：验证统计信息准确性
5. 内存管理：reset_cache!功能
"""

using Test

# 添加src路径
push!(LOAD_PATH, joinpath(@__DIR__, "../../src"))
push!(LOAD_PATH, joinpath(@__DIR__, "../../src/relaxtime"))

include("../../../src/relaxtime/PolarizationCache.jl")
using .PolarizationCache

@testset "极化函数缓存测试" begin
    
    # 测试参数
    channel = :P
    k0 = 100.0
    k_norm = 50.0
    m1 = 5.0
    m2 = 5.0
    μ1 = 300.0
    μ2 = 300.0
    T = 150.0
    Φ = 0.5
    Φbar = 0.5
    ξ = 0.0
    A1 = -50.0
    A2 = -50.0
    num_s = 0
    
    # 测试1：基本缓存功能
    @testset "基本缓存功能" begin
        reset_cache!()
        
        # 第一次调用
        result1 = polarization_aniso_cached(channel, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, num_s)
        stats1 = get_cache_stats()
        
        @test stats1.total_calls == 1
        @test stats1.cache_hits == 0
        @test stats1.cache_misses == 1
        @test stats1.hit_rate == 0.0
        @test stats1.cache_size == 1
        
        # 第二次调用相同参数
        result2 = polarization_aniso_cached(channel, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, num_s)
        stats2 = get_cache_stats()
        
        @test result1 == result2
        @test stats2.total_calls == 2
        @test stats2.cache_hits == 1
        @test stats2.cache_misses == 1
        @test stats2.hit_rate ≈ 0.5
        @test stats2.cache_size == 1  # 缓存条目数不变
    end
    
    # 测试2：浮点数容差
    @testset "浮点数容差" begin
        reset_cache!()
        
        # 第一次调用
        result1 = polarization_aniso_cached(channel, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, num_s)
        
        # 第二次调用：k0有微小差异（在容差范围内）
        k0_tiny_diff = k0 + 1e-13 * k0  # 相对差异 1e-13 < EPS_CACHE
        result2 = polarization_aniso_cached(channel, k0_tiny_diff, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, num_s)
        
        stats = get_cache_stats()
        @test stats.cache_hits == 1  # 应该命中缓存
        @test stats.cache_size == 1
    end
    
    # 测试3：不同参数不应命中缓存
    @testset "不同参数" begin
        reset_cache!()
        
        # 调用不同参数
        result1 = polarization_aniso_cached(:P, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, num_s)
        result2 = polarization_aniso_cached(:S, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, num_s)  # 不同通道
        result3 = polarization_aniso_cached(:P, k0 + 10.0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, num_s)  # 不同k0
        
        stats = get_cache_stats()
        @test stats.total_calls == 3
        @test stats.cache_hits == 0  # 全部未命中
        @test stats.cache_size == 3  # 三个不同的缓存条目
    end
    
    # 测试4：缓存统计准确性
    @testset "缓存统计" begin
        reset_cache!()
        
        # 重复调用相同参数10次
        for i in 1:10
            polarization_aniso_cached(channel, k0, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, num_s)
        end
        
        stats = get_cache_stats()
        @test stats.total_calls == 10
        @test stats.cache_hits == 9  # 第一次未命中，后9次命中
        @test stats.cache_misses == 1
        @test stats.hit_rate ≈ 0.9
        @test stats.cache_size == 1
    end
    
    # 测试5：reset_cache!功能
    @testset "重置缓存" begin
        # 先填充一些缓存
        for i in 1:5
            polarization_aniso_cached(channel, k0 + i, k_norm, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, num_s)
        end
        
        stats_before = get_cache_stats()
        @test stats_before.cache_size > 0
        @test stats_before.total_calls > 0
        
        # 重置缓存
        reset_cache!()
        
        stats_after = get_cache_stats()
        @test stats_after.cache_size == 0
        @test stats_after.total_calls == 0
        @test stats_after.cache_hits == 0
        @test stats_after.cache_misses == 0
        @test stats_after.hit_rate == 0.0
    end
    
    # 测试7：多参数组合场景（模拟实际使用）
    @testset "实际使用场景模拟" begin
        reset_cache!()
        
        # 模拟计算多个散射过程和散射道
        k0_values = [80.0, 100.0, 120.0, 100.0, 80.0]  # 有重复
        k_values = [40.0, 50.0, 60.0, 50.0, 40.0]      # 有重复
        
        for (k0_val, k_val) in zip(k0_values, k_values)
            # 计算P通道
            polarization_aniso_cached(:P, k0_val, k_val, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, num_s)
            # 计算S通道
            polarization_aniso_cached(:S, k0_val, k_val, m1, m2, μ1, μ2, T, Φ, Φbar, ξ, A1, A2, num_s)
        end
        
        stats = get_cache_stats()
        @test stats.total_calls == 10
        # 前3组(k0,k)是新的，后2组重复 → 6个新缓存，4次命中
        @test stats.cache_hits == 4
        @test stats.cache_size == 6  # 3种(k0,k) × 2种通道
        @test stats.hit_rate ≈ 0.4
    end
end

println("\n" * "="^70)
println("极化函数缓存测试完成！")
println("="^70)


