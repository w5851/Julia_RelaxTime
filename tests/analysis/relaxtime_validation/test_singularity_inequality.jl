"""
测试 singularity_k_positive 函数返回的能量区间是否满足原始不等式

验证区间内的 E 满足：
|λ² + 2λE - k² + m² - m'²| < 2k√(E² - m²)

运行方式：
```julia
using Pkg
Pkg.activate(".")
include("test/test_singularity_inequality.jl")
```
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
using Test
using Printf

include("../../../src/relaxtime/OneLoopIntegrals.jl")
using .OneLoopIntegrals: singularity_k_positive

"""
检查能量 E 是否满足原始不等式
|λ² + 2λE - k² + m² - m'²| < 2k√(E² - m²)
允许边界点有数值误差（< 1e-10）
"""
function check_inequality(E::Float64, λ::Float64, k::Float64, m::Float64, m_prime::Float64; 
                         tol::Float64=1e-10)
    # 左边: |λ² + 2λE - k² + m² - m'²|
    left_side = abs(λ^2 + 2*λ*E - k^2 + m^2 - m_prime^2)
    
    # 右边: 2k√(E² - m²)
    # 确保 E >= m (物理约束)
    if E < m
        return false, left_side, 0.0, "E < m"
    end
    
    p = sqrt(E^2 - m^2)
    right_side = 2*k*p
    
    # 检查不等式，允许边界点的数值误差
    diff = left_side - right_side
    satisfied = diff < tol  # 允许边界点 (left ≈ right)
    
    return satisfied, left_side, right_side, satisfied ? "✓" : "✗"
end

# 积分区间的默认设置 (需要与实际物理场景一致)
const DEFAULT_EMIN = 0.24  # 通常是 m
const DEFAULT_EMAX = 10.0  # 足够大的上限

"""
测试一个区间内的多个点是否都满足不等式
"""
function test_interval(interval::Tuple{Float64,Float64}, λ::Float64, k::Float64, 
                       m::Float64, m_prime::Float64; n_points::Int=20)
    E_min, E_max = interval
    
    println("\n  测试区间 [$(E_min), $(E_max)]:")
    println("  " * "-"^70)
    
    # 在区间内均匀采样
    E_values = range(E_min, E_max, length=n_points)
    
    all_satisfied = true
    violations = 0
    
    for E in E_values
        satisfied, left, right, status = check_inequality(E, λ, k, m, m_prime)
        
        if !satisfied
            all_satisfied = false
            violations += 1
            println(@sprintf("  E = %.6f: %s (左边=%.6e, 右边=%.6e, 差值=%.6e)", 
                           E, status, left, right, left - right))
        end
    end
    
    if all_satisfied
        println("  ✓ 区间内所有 $(n_points) 个测试点都满足不等式")
    else
        println("  ✗ 区间内有 $(violations)/$(n_points) 个点违反不等式")
    end
    
    return all_satisfied
end

"""
测试区间边界是否接近等号 (边界条件)
"""
function test_boundary(E::Float64, λ::Float64, k::Float64, m::Float64, m_prime::Float64)
    left_side = abs(λ^2 + 2*λ*E - k^2 + m^2 - m_prime^2)
    p = sqrt(E^2 - m^2)
    right_side = 2*k*p
    
    # 边界处应该接近相等
    diff = abs(left_side - right_side)
    is_boundary = diff < 1e-6  # 允许数值误差
    
    return is_boundary, left_side, right_side, diff
end

@testset "singularity_k_positive 不等式验证" begin
    
    @testset "情况1: a > 0 (λ² > k²)" begin
        # 使用标准测试参数
        λ = 0.45
        k = 0.30
        m = 0.24
        m_prime = 0.38
        Emin = m
        Emax = 10.0
        
        a = λ^2 - k^2
        
        println("\n" * "="^80)
        println("情况1: a > 0 (λ² > k²)")
        println("="^80)
        println(@sprintf("参数: λ=%.3f, k=%.3f, m=%.3f, m'=%.3f", λ, k, m, m_prime))
        println(@sprintf("a = λ² - k² = %.6f > 0", a))
        println(@sprintf("积分区间: [%.3f, %.3f]", Emin, Emax))
        
        intervals, sign_type = singularity_k_positive(λ, k, m, m_prime, Emin, Emax)
        
        println(@sprintf("返回: sign_type = %s, 区间数 = %d", sign_type, length(intervals)))
        
        if sign_type == :none
            println("✓ 判别式 < 0，无奇点区间（预期行为）")
            @test isempty(intervals)
        elseif sign_type == :between
            println("期望：奇点在两根之间 (E₁ < E < E₂)")
            @test length(intervals) == 1
            
            if !isempty(intervals)
                interval = intervals[1]
                E1, E2 = interval
                
                # 测试边界
                println("\n  边界测试:")
                is_bound1, left1, right1, diff1 = test_boundary(E1, λ, k, m, m_prime)
                is_bound2, left2, right2, diff2 = test_boundary(E2, λ, k, m, m_prime)
                
                println(@sprintf("  E₁ = %.6f: 边界条件 %s (|左-右| = %.3e)", 
                               E1, is_bound1 ? "✓" : "✗", diff1))
                println(@sprintf("  E₂ = %.6f: 边界条件 %s (|左-右| = %.3e)", 
                               E2, is_bound2 ? "✓" : "✗", diff2))
                
                # 测试区间内的点
                all_satisfied = test_interval(interval, λ, k, m, m_prime)
                @test all_satisfied
            end
        end
        println("="^80)
    end
    
    @testset "情况2: a < 0 (λ² < k²)" begin
        # k > λ 的情况
        λ = 0.45
        k = 0.50  # k > λ
        m = 0.24
        m_prime = 0.38
        Emin = m
        Emax = 10.0
        
        a = λ^2 - k^2
        
        println("\n" * "="^80)
        println("情况2: a < 0 (λ² < k²)")
        println("="^80)
        println(@sprintf("参数: λ=%.3f, k=%.3f, m=%.3f, m'=%.3f", λ, k, m, m_prime))
        println(@sprintf("a = λ² - k² = %.6f < 0", a))
        println(@sprintf("积分区间: [%.3f, %.3f]", Emin, Emax))
        
        intervals, sign_type = singularity_k_positive(λ, k, m, m_prime, Emin, Emax)
        
        println(@sprintf("返回: sign_type = %s, 区间数 = %d", sign_type, length(intervals)))
        
        if sign_type == :none
            println("✓ 判别式 < 0，无奇点区间")
            @test isempty(intervals)
        elseif sign_type == :outside
            println("期望：奇点在两根之外 (E < E₁ 或 E > E₂)")
            @test length(intervals) >= 1  # 可能有1个或2个区间
            
            for (i, interval) in enumerate(intervals)
                println("\n区间 $i:")
                all_satisfied = test_interval(interval, λ, k, m, m_prime)
                @test all_satisfied
            end
        end
        println("="^80)
    end
    
    @testset "情况3: a = 0 (λ² = k²)" begin
        # λ = k 的情况
        λ = 0.45
        k = 0.45  # k = λ
        m = 0.24
        m_prime = 0.38  # m' > m
        Emin = m
        Emax = 10.0
        
        a = λ^2 - k^2
        
        println("\n" * "="^80)
        println("情况3: a = 0 (λ² = k²)")
        println("="^80)
        println(@sprintf("参数: λ=%.3f, k=%.3f, m=%.3f, m'=%.3f", λ, k, m, m_prime))
        println(@sprintf("a = λ² - k² = %.6f ≈ 0", a))
        println(@sprintf("条件检查: λ > 0? %s, m' > m? %s", λ > 0, m_prime > m))
        println(@sprintf("积分区间: [%.3f, %.3f]", Emin, Emax))
        
        intervals, sign_type = singularity_k_positive(λ, k, m, m_prime, Emin, Emax)
        
        println(@sprintf("返回: sign_type = %s, 区间数 = %d", sign_type, length(intervals)))
        
        if λ > 0 && m_prime > m
            # a=0 时，如果满足条件，应该返回区间（sign_type 可能是 :none 但有区间）
            @test length(intervals) == 1
            
            if !isempty(intervals)
                interval = intervals[1]
                E_min, E_max = interval
                
                # 计算临界能量
                E_crit = (m_prime^2 - m^2)/(4*λ) + λ*m^2/(m_prime^2 - m^2)
                println(@sprintf("  临界能量 E_crit = %.6f", E_crit))
                println(@sprintf("  返回区间: [%.6f, %.6f]", E_min, E_max))
                
                # E_min 应该接近 E_crit
                @test abs(E_min - E_crit) < 1e-6
                
                all_satisfied = test_interval(interval, λ, k, m, m_prime)
                @test all_satisfied
            end
        else
            println("  不满足条件 (λ > 0 且 m' > m)，无解")
            @test sign_type == :none
            @test isempty(intervals)
        end
        println("="^80)
    end
    
    @testset "边缘情况测试" begin
        println("\n" * "="^80)
        println("边缘情况测试")
        println("="^80)
        
        # 测试 1: 判别式接近零
        println("\n子测试 1: 判别式接近零")
        λ = 0.5
        k = 0.3
        m = 0.3
        m_prime = 0.31  # 使判别式接近零
        Emin = m
        Emax = 10.0
        
        intervals, sign_type = singularity_k_positive(λ, k, m, m_prime, Emin, Emax)
        println(@sprintf("  λ=%.2f, k=%.2f, m=%.2f, m'=%.2f", λ, k, m, m_prime))
        println(@sprintf("  返回: sign_type=%s, 区间数=%d", sign_type, length(intervals)))
        
        if sign_type != :none
            for interval in intervals
                all_satisfied = test_interval(interval, λ, k, m, m_prime, n_points=10)
                @test all_satisfied
            end
        end
        
        # 测试 2: λ 很小但 > 0
        println("\n子测试 2: λ 很小但 > 0")
        λ = 0.01
        k = 0.01
        m = 0.24
        m_prime = 0.38
        Emin = m
        Emax = 10.0
        
        intervals, sign_type = singularity_k_positive(λ, k, m, m_prime, Emin, Emax)
        println(@sprintf("  λ=%.2f, k=%.2f, m=%.2f, m'=%.2f", λ, k, m, m_prime))
        println(@sprintf("  返回: sign_type=%s, 区间数=%d", sign_type, length(intervals)))
        
        if sign_type != :none
            for interval in intervals
                all_satisfied = test_interval(interval, λ, k, m, m_prime, n_points=10)
                @test all_satisfied
            end
        end
        
        # 测试 3: m' 接近 m
        println("\n子测试 3: m' 接近 m")
        λ = 0.45
        k = 0.30
        m = 0.24
        m_prime = 0.2401  # m' 非常接近 m
        Emin = m
        Emax = 10.0
        
        intervals, sign_type = singularity_k_positive(λ, k, m, m_prime, Emin, Emax)
        println(@sprintf("  λ=%.2f, k=%.2f, m=%.4f, m'=%.4f", λ, k, m, m_prime))
        println(@sprintf("  返回: sign_type=%s, 区间数=%d", sign_type, length(intervals)))
        
        if sign_type != :none
            for interval in intervals
                all_satisfied = test_interval(interval, λ, k, m, m_prime, n_points=10)
                @test all_satisfied
            end
        end
        
        println("="^80)
    end
    
    @testset "多组参数的系统性测试" begin
        println("\n" * "="^80)
        println("多组参数的系统性测试")
        println("="^80)
        
        test_cases = [
            # (λ, k, m, m_prime, description)
            (0.45, 0.30, 0.24, 0.38, "标准参数 (a>0)"),
            (0.45, 0.50, 0.24, 0.38, "k>λ (a<0)"),
            (0.45, 0.45, 0.24, 0.38, "k=λ (a=0)"),
            (0.60, 0.40, 0.20, 0.35, "大λ (a>0)"),
            (0.30, 0.50, 0.20, 0.35, "小λ (a<0)"),
            (0.35, 0.35, 0.15, 0.30, "λ=k, 小质量"),
            (0.50, 0.20, 0.30, 0.40, "小k (a>0)"),
            (0.20, 0.50, 0.30, 0.40, "大k (a<0)"),
        ]
        
        all_passed = true
        
        for (i, (λ, k, m, m_prime, desc)) in enumerate(test_cases)
            println("\n测试用例 $i: $desc")
            println(@sprintf("  参数: λ=%.2f, k=%.2f, m=%.2f, m'=%.2f", λ, k, m, m_prime))
            
            a = λ^2 - k^2
            println(@sprintf("  a = λ² - k² = %.6f", a))
            
            Emin = m
            Emax = 10.0
            
            intervals, sign_type = singularity_k_positive(λ, k, m, m_prime, Emin, Emax)
            println(@sprintf("  结果: sign_type=%s, 区间数=%d", sign_type, length(intervals)))
            
            case_passed = true
            
            if sign_type != :none && !isempty(intervals)
                for (j, interval) in enumerate(intervals)
                    E_min, E_max = interval
                    println(@sprintf("    区间 %d: [%.6f, %.6f]", j, E_min, E_max))
                    
                    # 只测试几个点以节省时间
                    satisfied = test_interval(interval, λ, k, m, m_prime, n_points=5)
                    if !satisfied
                        case_passed = false
                        all_passed = false
                    end
                end
            else
                println("    无奇点区间")
            end
            
            @test case_passed
            println("  " * (case_passed ? "✓ 通过" : "✗ 失败"))
        end
        
        println("\n" * "="^80)
        if all_passed
            println("✓ 所有测试用例都通过了不等式验证")
        else
            println("✗ 部分测试用例未通过不等式验证")
        end
        println("="^80)
    end
end

println("\n✓ singularity_k_positive 不等式验证测试完成")



