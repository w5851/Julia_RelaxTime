"""
测试 GaussLegendre 模块

运行方式：
```julia
using Pkg
Pkg.activate(".")
Pkg.test("RelaxTime")
```

或在 REPL 中：
```julia
include("test/test_gausslegendre.jl")
```
"""

using Test

# 加载模块
include("../../../src/integration/GaussLegendre.jl")
using .GaussLegendre

@testset "GaussLegendre 模块测试" begin
    
    @testset "基本功能测试" begin
        # 测试标准区间 [0, 1]
        nodes, weights = gauleg(0.0, 1.0, 5)
        
        @test length(nodes) == 5
        @test length(weights) == 5
        @test all(0.0 .< nodes .< 1.0)  # 节点在区间内
        @test all(weights .> 0.0)       # 权重为正
    end
    
    @testset "参数验证测试" begin
        # 测试无效节点数
        @test_throws ArgumentError gauleg(0.0, 1.0, 0)
        @test_throws ArgumentError gauleg(0.0, 1.0, -5)
        
        # 测试无效区间
        @test_throws ArgumentError gauleg(1.0, 0.0, 5)
        @test_throws ArgumentError gauleg(1.0, 1.0, 5)
    end
    
    @testset "积分精度测试 - 多项式" begin
        # 高斯-勒让德积分可以精确积分 2n-1 次多项式
        n = 5  # 5个节点应该精确积分至多9次多项式
        
        # 测试 ∫₀¹ x² dx = 1/3
        nodes, weights = gauleg(0.0, 1.0, n)
        f1(x) = x^2
        result1 = sum(weights .* f1.(nodes))
        @test isapprox(result1, 1/3, rtol=1e-12)
        
        # 测试 ∫₀¹ x⁴ dx = 1/5
        f2(x) = x^4
        result2 = sum(weights .* f2.(nodes))
        @test isapprox(result2, 1/5, rtol=1e-12)
        
        # 测试 ∫₀¹ (3x² - 2x + 1) dx = 1
        f3(x) = 3x^2 - 2x + 1
        result3 = sum(weights .* f3.(nodes))
        @test isapprox(result3, 1.0, rtol=1e-12)
    end
    
    @testset "积分精度测试 - 超越函数" begin
        # 测试 ∫₀^π sin(x) dx = 2
        nodes, weights = gauleg(0.0, Float64(π), 20)
        f_sin(x) = sin(x)
        result_sin = sum(weights .* f_sin.(nodes))
        @test isapprox(result_sin, 2.0, rtol=1e-10)
        
        # 测试 ∫₀¹ exp(x) dx = e - 1
        nodes, weights = gauleg(0.0, 1.0, 15)
        f_exp(x) = exp(x)
        result_exp = sum(weights .* f_exp.(nodes))
        @test isapprox(result_exp, ℯ - 1, rtol=1e-10)
        
        # 测试 ∫₀^(π/2) cos(x) dx = 1
        nodes, weights = gauleg(0.0, Float64(π)/2, 20)
        f_cos(x) = cos(x)
        result_cos = sum(weights .* f_cos.(nodes))
        @test isapprox(result_cos, 1.0, rtol=1e-10)
    end
    
    @testset "不同区间测试" begin
        # 测试 ∫₋₁¹ x² dx = 2/3
        nodes, weights = gauleg(-1.0, 1.0, 10)
        f(x) = x^2
        result = sum(weights .* f.(nodes))
        @test isapprox(result, 2/3, rtol=1e-12)
        
        # 测试 ∫₂⁵ x dx = 21/2
        nodes, weights = gauleg(2.0, 5.0, 10)
        f2(x) = x
        result2 = sum(weights .* f2.(nodes))
        @test isapprox(result2, 21/2, rtol=1e-12)
    end
    
    @testset "权重和检验" begin
        # 对于常函数 f(x) = 1，积分应该等于区间长度
        # ∫ₐᵇ 1 dx = b - a
        
        # 区间 [0, 1]
        nodes, weights = gauleg(0.0, 1.0, 10)
        @test isapprox(sum(weights), 1.0, rtol=1e-12)
        
        # 区间 [-1, 1]
        nodes, weights = gauleg(-1.0, 1.0, 15)
        @test isapprox(sum(weights), 2.0, rtol=1e-12)
        
        # 区间 [2.5, 7.8]
        a, b = 2.5, 7.8
        nodes, weights = gauleg(a, b, 20)
        @test isapprox(sum(weights), b - a, rtol=1e-12)
    end
    
    @testset "节点数收敛性测试" begin
        # 测试随着节点数增加，积分精度提高
        # ∫₀¹ exp(x²) dx ≈ 1.4626517459
        exact = 1.4626517459
        f(x) = exp(x^2)
        
        errors = Float64[]
        for n in [5, 10, 20, 50]
            nodes, weights = gauleg(0.0, 1.0, n)
            result = sum(weights .* f.(nodes))
            push!(errors, abs(result - exact))
        end
        
        # 验证误差递减（前几次明显递减）
        @test errors[2] < errors[1]
        @test errors[3] < errors[2]
        # 当误差达到机器精度附近时，可能不再严格递减
        @test errors[4] < 1e-10  # 50个节点应该达到很高精度
    end
    
    @testset "物理应用示例测试" begin
        # 模拟物理计算场景（自然单位制）
        T = 0.1   # 温度，fm⁻¹
        μ = 0.3   # 化学势，fm⁻¹
        Λf = 3.05 # 截断参数，fm⁻¹
        
        # 费米-狄拉克分布积分
        fermi_dirac(k) = 1.0 / (exp((k - μ) / T) + 1.0)
        
        nodes, weights = gauleg(0.0, Λf, 50)
        integral = sum(weights .* fermi_dirac.(nodes))
        
        # 验证结果在合理范围内
        @test 0.0 < integral < Λf  # 积分值应该在 0 到 Λf 之间
        @test isfinite(integral)    # 结果应该是有限的
    end
    
    @testset "边界情况测试" begin
        # 测试很小的区间
        nodes, weights = gauleg(0.0, 1e-10, 5)
        @test isapprox(sum(weights), 1e-10, rtol=1e-12)
        
        # 测试很大的区间
        nodes, weights = gauleg(0.0, 1e10, 5)
        @test isapprox(sum(weights), 1e10, rtol=1e-12)
        
        # 测试单个节点
        nodes, weights = gauleg(0.0, 1.0, 1)
        @test length(nodes) == 1
        @test isapprox(sum(weights), 1.0, rtol=1e-12)
    end
    
    @testset "对称性测试" begin
        # 对于对称区间和对称函数，节点应该关于中心对称
        n = 10
        nodes, weights = gauleg(-1.0, 1.0, n)
        
        # 检查节点的对称性
        center = 0.0
        for i in 1:div(n, 2)
            @test isapprox(nodes[i] + nodes[end-i+1], 2*center, atol=1e-12)
            @test isapprox(weights[i], weights[end-i+1], rtol=1e-12)
        end
    end
    
end

# 运行所有测试并显示结果
println("\n" * "="^70)
println("GaussLegendre 模块测试完成！")
println("="^70)


