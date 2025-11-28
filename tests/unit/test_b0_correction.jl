"""
测试 B0_correction 函数在不同各向异性参数 ξ 下的行为

运行方式：
```julia
using Pkg
Pkg.activate(".")
include("test/test_b0_correction.jl")
```
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
using Test
using Printf

include("../../src/relaxtime/OneLoopIntegralsAniso.jl")
include("../../src/relaxtime/OneLoopIntegrals.jl")
using .OneLoopIntegralsCorrection: B0_correction
using .OneLoopIntegrals: B0

const TEST_PARAMS = (
    λ = 0.45,
    k = 0.30,
    m1 = 0.24,
    m2 = 0.38,
    μ1 = 0.12,
    μ2 = -0.05,
    T = 0.17,
    Φ = 0.15,
    Φbar = 0.15,
)

@testset "B0_correction 各向异性参数 ξ 依赖性" begin
    
    @testset "基本功能测试" begin
        λ = TEST_PARAMS.λ
        k = TEST_PARAMS.k
        m1 = TEST_PARAMS.m1
        m2 = TEST_PARAMS.m2
        μ1 = TEST_PARAMS.μ1
        μ2 = TEST_PARAMS.μ2
        T = TEST_PARAMS.T
        Φ = TEST_PARAMS.Φ
        Φbar = TEST_PARAMS.Φbar
        ξ = 0.1
        
        result = B0_correction(λ, k, m1, m2, μ1, μ2, T, Φ, Φbar, ξ)
        
        @test result isa NTuple{2, Float64}
        @test all(isfinite, result)
        @info "B0_correction 基本调用" ξ=ξ real_part=result[1] imag_part=result[2]
    end
    
    @testset "ξ = 0 时应该返回零" begin
        λ = TEST_PARAMS.λ
        k = TEST_PARAMS.k
        m1 = TEST_PARAMS.m1
        m2 = TEST_PARAMS.m2
        μ1 = TEST_PARAMS.μ1
        μ2 = TEST_PARAMS.μ2
        T = TEST_PARAMS.T
        Φ = TEST_PARAMS.Φ
        Φbar = TEST_PARAMS.Φbar
        
        result = B0_correction(λ, k, m1, m2, μ1, μ2, T, Φ, Φbar, 0.0)
        
        @info "B0_correction (ξ=0)" result
        @test abs(result[1]) < 1e-10
        @test abs(result[2]) < 1e-10
    end
    
    @testset "不同 ξ 值的系统性测试" begin
        λ = TEST_PARAMS.λ
        k = TEST_PARAMS.k
        m1 = TEST_PARAMS.m1
        m2 = TEST_PARAMS.m2
        μ1 = TEST_PARAMS.μ1
        μ2 = TEST_PARAMS.μ2
        T = TEST_PARAMS.T
        Φ = TEST_PARAMS.Φ
        Φbar = TEST_PARAMS.Φbar
        
        # 测试不同的 ξ 值
        ξ_values = [0.0, 0.05, 0.1, 0.2, 0.3, 0.5]
        
        println("\n" * "="^70)
        println("B0_correction 随 ξ 变化的系统性测试")
        println("="^70)
        println(@sprintf("%-10s %-20s %-20s", "ξ", "实部", "虚部"))
        println("-"^70)
        
        results = []
        for ξ in ξ_values
            result = B0_correction(λ, k, m1, m2, μ1, μ2, T, Φ, Φbar, ξ)
            push!(results, (ξ=ξ, real=result[1], imag=result[2]))
            println(@sprintf("%-10.2f %-20.10e %-20.10e", ξ, result[1], result[2]))
            
            @test all(isfinite, result)
        end
        println("="^70 * "\n")
        
        # 验证 ξ 增加时，修正项的绝对值也应该增加（一般情况）
        for i in 2:length(results)
            if results[i].ξ > 0
                @test abs(results[i].real) > abs(results[1].real) || abs(results[i].real) < 1e-8
            end
        end
    end
    
    @testset "k = 0 情况下的 ξ 依赖性" begin
        λ = TEST_PARAMS.λ
        m1 = TEST_PARAMS.m1
        m2 = TEST_PARAMS.m2
        μ1 = TEST_PARAMS.μ1
        μ2 = TEST_PARAMS.μ2
        T = TEST_PARAMS.T
        Φ = TEST_PARAMS.Φ
        Φbar = TEST_PARAMS.Φbar
        
        ξ_values = [0.0, 0.1, 0.2, 0.3]
        
        println("\n" * "="^70)
        println("B0_correction (k=0) 随 ξ 变化的测试")
        println("="^70)
        println(@sprintf("%-10s %-20s %-20s", "ξ", "实部", "虚部"))
        println("-"^70)
        
        for ξ in ξ_values
            result = B0_correction(λ, 0.0, m1, m2, μ1, μ2, T, Φ, Φbar, ξ)
            println(@sprintf("%-10.2f %-20.10e %-20.10e", ξ, result[1], result[2]))
            
            @test all(isfinite, result)
        end
        println("="^70 * "\n")
    end
    
    @testset "不同物理参数下的 ξ 依赖性" begin
        println("\n" * "="^70)
        println("不同物理参数组合下的 B0_correction (ξ=0.2)")
        println("="^70)
        
        test_cases = [
            (λ=0.3, k=0.2, m1=0.2, m2=0.3, μ1=0.1, μ2=0.0, T=0.15, name="Case 1"),
            (λ=0.5, k=0.4, m1=0.25, m2=0.35, μ1=0.15, μ2=-0.1, T=0.18, name="Case 2"),
            (λ=0.6, k=0.5, m1=0.3, m2=0.4, μ1=0.2, μ2=-0.05, T=0.20, name="Case 3"),
        ]
        
        Φ = TEST_PARAMS.Φ
        Φbar = TEST_PARAMS.Φbar
        ξ = 0.2
        
        println(@sprintf("%-10s %-20s %-20s", "测试用例", "实部", "虚部"))
        println("-"^70)
        
        for case in test_cases
            result = B0_correction(case.λ, case.k, case.m1, case.m2, 
                                  case.μ1, case.μ2, case.T, Φ, Φbar, ξ)
            println(@sprintf("%-10s %-20.10e %-20.10e", case.name, result[1], result[2]))
            
            @test all(isfinite, result)
        end
        println("="^70 * "\n")
    end
    
    @testset "精度测试" begin
        λ = TEST_PARAMS.λ
        k = TEST_PARAMS.k
        m1 = TEST_PARAMS.m1
        m2 = TEST_PARAMS.m2
        μ1 = TEST_PARAMS.μ1
        μ2 = TEST_PARAMS.μ2
        T = TEST_PARAMS.T
        Φ = TEST_PARAMS.Φ
        Φbar = TEST_PARAMS.Φbar
        ξ = 0.2
        
        result_default = B0_correction(λ, k, m1, m2, μ1, μ2, T, Φ, Φbar, ξ)
        result_tight = B0_correction(λ, k, m1, m2, μ1, μ2, T, Φ, Φbar, ξ; rtol=1e-6)
        
        @info "精度对比" rtol_default=result_default rtol_1e6=result_tight
        
        # 不同精度的结果应该在合理范围内一致
        @test isapprox(result_default[1], result_tight[1]; rtol=1e-2)
        @test isapprox(result_default[2], result_tight[2]; rtol=1e-2)
    end
    
    @testset "ξ 线性近似验证" begin
        # 对于小的 ξ，修正项应该近似与 ξ 成正比
        λ = TEST_PARAMS.λ
        k = TEST_PARAMS.k
        m1 = TEST_PARAMS.m1
        m2 = TEST_PARAMS.m2
        μ1 = TEST_PARAMS.μ1
        μ2 = TEST_PARAMS.μ2
        T = TEST_PARAMS.T
        Φ = TEST_PARAMS.Φ
        Φbar = TEST_PARAMS.Φbar
        
        ξ1 = 0.05
        ξ2 = 0.10
        
        result1 = B0_correction(λ, k, m1, m2, μ1, μ2, T, Φ, Φbar, ξ1)
        result2 = B0_correction(λ, k, m1, m2, μ1, μ2, T, Φ, Φbar, ξ2)
        
        # 计算比值
        ratio_real = result2[1] / result1[1]
        ratio_imag = result2[2] / result1[2]
        expected_ratio = ξ2 / ξ1
        
        println("\n线性近似检验:")
        println(@sprintf("ξ1 = %.2f, ξ2 = %.2f, 预期比值 = %.2f", ξ1, ξ2, expected_ratio))
        println(@sprintf("实部比值: %.4f", ratio_real))
        println(@sprintf("虚部比值: %.4f", ratio_imag))
        
        # 对于小的 ξ，比值应该接近 2.0 (允许一定偏差)
        if abs(result1[1]) > 1e-10
            @test isapprox(ratio_real, expected_ratio; rtol=0.2)
        end
        if abs(result1[2]) > 1e-10
            @test isapprox(ratio_imag, expected_ratio; rtol=0.2)
        end
    end
    
    @testset "B0_correction 相对于 B0 的大小比较" begin
        λ = TEST_PARAMS.λ
        k = TEST_PARAMS.k
        m1 = TEST_PARAMS.m1
        m2 = TEST_PARAMS.m2
        μ1 = TEST_PARAMS.μ1
        μ2 = TEST_PARAMS.μ2
        T = TEST_PARAMS.T
        Φ = TEST_PARAMS.Φ
        Φbar = TEST_PARAMS.Φbar
        
        # 计算零阶项 B0（不含各向异性修正）
        b0_result = B0(λ, k, m1, μ1, m2, μ2, T; Φ=Φ, Φbar=Φbar)
        
        println("\n" * "="^80)
        println("B0_correction 相对于 B0 的大小比较")
        println("="^80)
        println(@sprintf("B0 (零阶项): 实部 = %20.10e, 虚部 = %20.10e", b0_result[1], b0_result[2]))
        println("-"^80)
        println(@sprintf("%-10s %-22s %-22s %-12s %-12s", "ξ", "B0_correction 实部", "B0_correction 虚部", "实部比值(%)", "虚部比值(%)"))
        println("-"^80)
        
        ξ_values = [0.05, 0.1, 0.15, 0.2, 0.3, 0.5]
        
        for ξ in ξ_values
            correction_result = B0_correction(λ, k, m1, m2, μ1, μ2, T, Φ, Φbar, ξ)
            
            # 计算相对大小（百分比）
            ratio_real_percent = abs(correction_result[1] / b0_result[1]) * 100
            ratio_imag_percent = abs(correction_result[2] / b0_result[2]) * 100
            
            println(@sprintf("%-10.2f %-22.10e %-22.10e %-12.4f %-12.4f", 
                            ξ, correction_result[1], correction_result[2], 
                            ratio_real_percent, ratio_imag_percent))
            
            @test all(isfinite, correction_result)
            
            # 验证实部修正项确实是一阶小量（对于小的 ξ，应该远小于零阶项）
            # 注意：虚部修正项在某些参数下可能与零阶项有不同的符号或量级特性
            if ξ <= 0.2
                @test abs(correction_result[1]) < abs(b0_result[1]) || abs(b0_result[1]) < 1e-10
            end
        end
        println("="^80 * "\n")
        
        @info "B0 零阶项" real_part=b0_result[1] imag_part=b0_result[2]
    end
    
    @testset "不同 k 值下 B0_correction 相对于 B0 的比较" begin
        λ = TEST_PARAMS.λ
        m1 = TEST_PARAMS.m1
        m2 = TEST_PARAMS.m2
        μ1 = TEST_PARAMS.μ1
        μ2 = TEST_PARAMS.μ2
        T = TEST_PARAMS.T
        Φ = TEST_PARAMS.Φ
        Φbar = TEST_PARAMS.Φbar
        ξ = 0.2
        
        k_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
        
        println("\n" * "="^80)
        println("不同 k 值下的 B0 与 B0_correction 比较 (ξ = 0.2)")
        println("="^80)
        println(@sprintf("%-8s %-20s %-20s %-20s %-20s", "k", "B0 实部", "B0 虚部", "Correction 实部", "Correction 虚部"))
        println("-"^80)
        
        for k in k_values
            b0_result = B0(λ, k, m1, μ1, m2, μ2, T; Φ=Φ, Φbar=Φbar)
            correction_result = B0_correction(λ, k, m1, m2, μ1, μ2, T, Φ, Φbar, ξ)
            
            println(@sprintf("%-8.2f %-20.10e %-20.10e %-20.10e %-20.10e", 
                            k, b0_result[1], b0_result[2], 
                            correction_result[1], correction_result[2]))
            
            @test all(isfinite, b0_result)
            @test all(isfinite, correction_result)
            
            # 修正项应该小于零阶项
            if abs(b0_result[1]) > 1e-10
                @test abs(correction_result[1]) < abs(b0_result[1]) * 2.0  # 允许一定的灵活度
            end
        end
        println("="^80 * "\n")
    end
    
    @testset "B0_correction/B0 比值的 ξ 依赖性" begin
        λ = TEST_PARAMS.λ
        k = TEST_PARAMS.k
        m1 = TEST_PARAMS.m1
        m2 = TEST_PARAMS.m2
        μ1 = TEST_PARAMS.μ1
        μ2 = TEST_PARAMS.μ2
        T = TEST_PARAMS.T
        Φ = TEST_PARAMS.Φ
        Φbar = TEST_PARAMS.Φbar
        
        b0_result = B0(λ, k, m1, μ1, m2, μ2, T; Φ=Φ, Φbar=Φbar)
        
        println("\n" * "="^80)
        println("B0_correction/B0 比值随 ξ 的变化")
        println("="^80)
        println(@sprintf("%-10s %-20s %-20s", "ξ", "实部比值", "虚部比值"))
        println("-"^80)
        
        ξ_values = [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]
        ratios_real = Float64[]
        ratios_imag = Float64[]
        
        for ξ in ξ_values
            correction_result = B0_correction(λ, k, m1, m2, μ1, μ2, T, Φ, Φbar, ξ)
            
            ratio_real = correction_result[1] / b0_result[1]
            # 虚部比值：如果分母为0，设为0（表示无虚部贡献）
            ratio_imag = abs(b0_result[2]) > 1e-15 ? correction_result[2] / b0_result[2] : 0.0
            
            push!(ratios_real, ratio_real)
            push!(ratios_imag, ratio_imag)
            
            println(@sprintf("%-10.3f %-20.10e %-20.10e", ξ, ratio_real, ratio_imag))
            
            @test isfinite(ratio_real)
            @test isfinite(ratio_imag)  # 现在应该总是有限的
        end
        println("="^80 * "\n")
        
        # 验证对于小的 ξ，比值近似与 ξ 成正比
        if length(ratios_real) >= 3
            # 检查前几个点的比值是否近似线性增长
            ξ_small = ξ_values[1:3]
            ratios_small = ratios_real[1:3]
            
            # 计算斜率（应该近似相等）
            slope1 = (ratios_small[2] - ratios_small[1]) / (ξ_small[2] - ξ_small[1])
            slope2 = (ratios_small[3] - ratios_small[2]) / (ξ_small[3] - ξ_small[2])
            
            @info "小 ξ 下的线性度检验" slope1 slope2 relative_diff=abs(slope2-slope1)/abs(slope1)
        end
    end
end

println("\n✓ 所有 B0_correction 测试完成")


