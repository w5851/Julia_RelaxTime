"""
测试 CauchyPV 模块

运行方式：
```julia
using Pkg
Pkg.activate(".")
include("test/test_cauchypv.jl")
```
"""

using Test

# 加载模块
include("../../../src/integration/CauchyPV.jl")
using .CauchyPV

@testset "CauchyPV 模块测试" begin
    
    @testset "基本功能测试" begin
        # 测试标准区间 [0, 2]，奇点在中间
        nodes, weights, n_actual = cpvi(0.0, 2.0, 1.0, 100)
        
        @test length(nodes) == n_actual
        @test length(weights) == n_actual
        @test all(0.0 .< nodes .< 2.0)     # 节点在区间内
        @test all(nodes .!= 1.0)            # 没有节点在奇点上
        @test all(weights .> 0.0)           # 权重为正
        @test n_actual > 0                  # 生成了节点
    end
    
    @testset "参数验证测试" begin
        # 测试无效节点数
        @test_throws ArgumentError cpvi(0.0, 1.0, 0.5, 1)
        @test_throws ArgumentError cpvi(0.0, 1.0, 0.5, 0)
        @test_throws ArgumentError cpvi(0.0, 1.0, 0.5, -5)
        
        # 测试无效区间
        @test_throws ArgumentError cpvi(1.0, 0.0, 0.5, 100)
        @test_throws ArgumentError cpvi(1.0, 1.0, 0.5, 100)
        
        # 测试无效奇点位置
        @test_throws ArgumentError cpvi(0.0, 1.0, 0.0, 100)   # 奇点在左端点
        @test_throws ArgumentError cpvi(0.0, 1.0, 1.0, 100)   # 奇点在右端点
        @test_throws ArgumentError cpvi(0.0, 1.0, -0.5, 100)  # 奇点在区间外
        @test_throws ArgumentError cpvi(0.0, 1.0, 1.5, 100)   # 奇点在区间外
    end
    
    @testset "权重和检验" begin
        # 对于常函数 f(x) = 1，积分应该接近区间长度
        # 但由于奇点附近的间隙，会略小于 (b-a)
        
        a, b, x0 = 0.0, 2.0, 1.0
        nodes, weights, n_actual = cpvi(a, b, x0, 200)
        
        weight_sum = sum(weights)
        @test isapprox(weight_sum, b - a, rtol=0.01)  # 误差在1%内
        @test weight_sum < (b - a)  # 由于间隙，应该略小
    end
    
    @testset "对称性测试 - 奇函数" begin
        # 对于对称区间和居中奇点，奇函数的主值积分应为零
        a, b, x0 = -1.0, 1.0, 0.0
        nodes, weights, n_actual = cpvi(a, b, x0, 500)
        
        # f(x) = x（奇函数）
        # P.V. ∫₋₁¹ x/x dx = P.V. ∫₋₁¹ 1 dx = 0（由于对称性）
        # 但这里 f(x)/x 在 x=0 处可导，实际上不是真正的奇点
        # 更好的测试：P.V. ∫₋₁¹ x³/(x) dx = P.V. ∫₋₁¹ x² dx
        
        f(x) = x^3
        result = sum(weights[i] * f(nodes[i]) / nodes[i] for i in 1:n_actual)
        @test isapprox(result, 2/3, rtol=1e-2)  # ∫₋₁¹ x² dx = 2/3，允许1%误差
    end
    
    @testset "对称性测试 - 严格奇函数" begin
        # 测试严格的奇函数：f(x) = x
        # P.V. ∫₋ₐᵃ x/(x-0) dx 应该非常接近 0
        a, b, x0 = -2.0, 2.0, 0.0
        nodes, weights, n_actual = cpvi(a, b, x0, 1000)
        
        # 由于对称性，左右贡献应该抵消
        left_nodes = nodes[nodes .< x0]
        right_nodes = nodes[nodes .> x0]
        
        @test length(left_nodes) > 0
        @test length(right_nodes) > 0
        
        # 检查节点的对称性（近似）
        # 最小负节点的绝对值应该接近最大正节点
        @test isapprox(abs(minimum(left_nodes)), maximum(right_nodes), rtol=0.1)
    end
    
    @testset "节点避让奇点测试" begin
        # 验证所有节点确实避开了奇点
        a, b, x0 = 0.0, 10.0, 5.0
        nodes, weights, n_actual = cpvi(a, b, x0, 200)
        
        h = (b - a) / 199.0  # 基本步长
        gap = 0.001 * h       # 期望的间隙
        
        # 所有节点应该距离奇点至少 gap
        min_distance = minimum(abs.(nodes .- x0))
        @test min_distance >= gap * 0.9  # 允许数值误差
    end
    
    @testset "节点单调性测试" begin
        # 节点应该是单调递增的
        a, b, x0 = 0.0, 5.0, 2.5
        nodes, weights, n_actual = cpvi(a, b, x0, 100)
        
        @test issorted(nodes)  # 检查节点是否单调递增
    end
    
    @testset "积分精度测试 - 可积函数" begin
        # 测试一个光滑函数除以 (x-x₀) 的积分
        # P.V. ∫₀² (x²)/(x-1) dx
        
        a, b, x0 = 0.0, 2.0, 1.0
        nodes, weights, n_actual = cpvi(a, b, x0, 500)
        
        # f(x) = x²
        # 可以通过多项式长除法：x²/(x-1) = x + 1 + 1/(x-1)
        # P.V. ∫₀² [x + 1 + 1/(x-1)] dx = [x²/2 + x]₀² + 0 = 4
        f(x) = x^2
        result = sum(weights[i] * f(nodes[i]) / (nodes[i] - x0) for i in 1:n_actual)
        @test isapprox(result, 4.0, rtol=1e-2)
    end
    
    @testset "节点数收敛性测试" begin
        # 测试随着节点数增加，积分精度提高
        a, b, x0 = 0.0, 2.0, 1.0
        f(x) = x^2
        
        results = Float64[]
        for n in [50, 100, 200, 500]
            nodes, weights, n_actual = cpvi(a, b, x0, n)
            result = sum(weights[i] * f(nodes[i]) / (nodes[i] - x0) for i in 1:n_actual)
            push!(results, result)
        end
        
        # 验证收敛到某个值（精度提高）
        errors = abs.(results .- 4.0)
        @test errors[end] < errors[1]  # 最后的误差小于第一个
        @test errors[end] < 0.1        # 最后的误差足够小
    end
    
    @testset "物理应用示例测试" begin
        # 模拟色散关系计算（自然单位制）
        T = 0.1      # 温度，fm⁻¹
        μ = 0.3      # 化学势，fm⁻¹
        k0 = 0.5     # 外动量，fm⁻¹
        Λf = 3.05    # 截断参数，fm⁻¹
        
        nodes, weights, n_actual = cpvi(0.0, Λf, k0, 200)
        
        # 简单的谱函数
        function spectral(k)
            ω = sqrt(k^2 + 0.3^2)
            return 1.0 / (exp((ω - μ) / T) + 1.0)
        end
        
        # 色散积分
        integral = sum(weights[i] * spectral(nodes[i]) / (k0 - nodes[i]) 
                      for i in 1:n_actual)
        
        # 验证结果在合理范围内
        @test isfinite(integral)      # 结果应该是有限的
        @test !isnan(integral)         # 不应该是 NaN
    end
    
    @testset "边界情况测试" begin
        # 测试奇点靠近左端点
        nodes1, weights1, n1 = cpvi(0.0, 10.0, 0.1, 100)
        @test n1 > 0
        @test all(nodes1 .>= 0.0)
        @test all(nodes1 .<= 10.0)
        
        # 测试奇点靠近右端点
        nodes2, weights2, n2 = cpvi(0.0, 10.0, 9.9, 100)
        @test n2 > 0
        @test all(nodes2 .>= 0.0)
        @test all(nodes2 .<= 10.0)
        
        # 测试很小的区间
        nodes3, weights3, n3 = cpvi(0.0, 0.1, 0.05, 50)
        @test n3 > 0
        @test isapprox(sum(weights3), 0.1, rtol=0.1)
    end
    
    @testset "不同区间测试" begin
        # 测试负区间
        nodes, weights, n_actual = cpvi(-5.0, -1.0, -3.0, 100)
        @test all(-5.0 .<= nodes .<= -1.0)
        @test all(nodes .!= -3.0)
        
        # 测试跨零区间
        nodes2, weights2, n2 = cpvi(-2.0, 3.0, 1.0, 100)
        @test all(-2.0 .<= nodes2 .<= 3.0)
        @test all(nodes2 .!= 1.0)
    end
    
    @testset "实际节点数测试" begin
        # 验证实际节点数与期望节点数的关系
        for n_expect in [50, 100, 200]
            nodes, weights, n_actual = cpvi(0.0, 2.0, 1.0, n_expect)
            # 实际节点数应该接近期望值（在合理范围内）
            @test 0.5 * n_expect < n_actual < 1.5 * n_expect
        end
    end
    
    @testset "简单解析解验证" begin
        # P.V. ∫₀² 1/(x-1) dx = [ln|x-1|]₀² = ln(1) - ln(1) = 0
        # 实际上：ln|2-1| - ln|0-1| = ln(1) - ln(1) = 0
        # 但主值定义：lim[ln(1-ε) - ln(-ε) + ln(1) - ln(ε)]
        # = lim[ln(1-ε) + ln(1) - 2ln(ε)] = ln(1) + ln(1) = 0（当我们取对数差的主值）
        
        a, b, x0 = 0.0, 2.0, 1.0
        nodes, weights, n_actual = cpvi(a, b, x0, 1000)
        
        f(x) = 1.0
        result = sum(weights[i] * f(nodes[i]) / (nodes[i] - x0) for i in 1:n_actual)
        
        # 由于数值方法和间隙，结果应该接近 0 但可能有小误差
        @test abs(result) < 0.1  # 允许较大的数值误差
    end
    
    @testset "权重结构测试" begin
        # 检查梯形法则的权重结构
        nodes, weights, n_actual = cpvi(0.0, 4.0, 2.0, 200)
        
        h = (4.0 - 0.0) / 199.0
        
        # 第一个和最后一个权重应该是 h/2
        @test isapprox(weights[1], h/2, rtol=1e-10)
        @test isapprox(weights[end], h/2, rtol=1e-10)
        
        # 大部分内部权重应该是 h（除了奇点附近）
        # 统计等于 h 的权重数量
        h_weights = count(w -> isapprox(w, h, rtol=1e-10), weights)
        half_h_weights = count(w -> isapprox(w, h/2, rtol=1e-10), weights)
        
        @test h_weights > 0
        @test half_h_weights >= 2  # 至少端点是 h/2
    end
    
end

# 运行所有测试并显示结果
println("\n" * "="^70)
println("CauchyPV 模块测试完成！")
println("="^70)


