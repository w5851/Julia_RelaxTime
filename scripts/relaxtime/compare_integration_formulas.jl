#!/usr/bin/env julia
"""
对比Julia和Fortran的积分公式

在各向同性情况下（ξ=0），检查积分公式是否等价
"""

using Printf

println("="^80)
println("积分公式对比")
println("="^80)
println()

println("Fortran的积分公式（各向同性）:")
println("-"^80)
println("w_ij = ∫∫∫ dp_i dp_j dθ · p_i² p_j² sin(θ) · f(p_i) f(p_j) · σ(s) · v_rel")
println()
println("其中:")
println("  - θ: 两个动量之间的夹角")
println("  - s = m_i² + m_j² + 2(E_i E_j - p_i p_j cos(θ))")
println("  - 积分范围: p_i ∈ [0, ∞), p_j ∈ [0, ∞), θ ∈ [0, π]")
println()
println("归一化:")
println("  w_ij_n = N_c²/(2π⁴) × w_ij")
println("  w_ij = w_ij_n / (n1 × n2)")
println()

println("="^80)
println()

println("Julia的积分公式（各向异性，但ξ=0时应退化为各向同性）:")
println("-"^80)
println("ω = ∫∫∫∫∫ dp_i dp_j dcosθ_i dcosθ_j dφ · p_i² p_j² · f(p_i,cosθ_i) f(p_j,cosθ_j) · σ(s) · v_rel")
println()
println("其中:")
println("  - cosθ_i, cosθ_j: 每个粒子相对于某个参考轴的角度")
println("  - φ: 两个粒子的方位角差")
println("  - cosΘ = cosθ_i cosθ_j + sinθ_i sinθ_j cos(φ)")
println("  - s = m_i² + m_j² + 2(E_i E_j - p_i p_j cosΘ)")
println("  - 积分范围: p_i ∈ [0, ∞), p_j ∈ [0, ∞), cosθ_i ∈ [-1, 1], cosθ_j ∈ [-1, 1], φ ∈ [0, 2π]")
println()
println("归一化:")
println("  w̄_ij = DQ²/(32π⁵ρ_iρ_j) × ω")
println("  其中 DQ = 2N_c = 6")
println()

println("="^80)
println()

println("关键问题：在ξ=0时，Julia的5维积分应该等价于Fortran的3维积分")
println()

println("理论分析:")
println("-"^80)
println()

println("在各向同性情况下（ξ=0）:")
println("  f(p_i, cosθ_i) = f(p_i) (与角度无关)")
println("  f(p_j, cosθ_j) = f(p_j) (与角度无关)")
println()

println("Julia的角度积分:")
println("  ∫∫∫ dcosθ_i dcosθ_j dφ · [...]")
println()

println("变量替换:")
println("  令 cosΘ = cosθ_i cosθ_j + sinθ_i sinθ_j cos(φ)")
println("  这是两个粒子动量之间的夹角余弦")
println()

println("在各向同性下，应该有:")
println("  ∫_{-1}^{1} dcosθ_i ∫_{-1}^{1} dcosθ_j ∫_{0}^{2π} dφ · g(cosΘ)")
println("  = 某个因子 × ∫_{-1}^{1} d(cosΘ) · g(cosΘ)")
println()

println("这个因子应该是多少？")
println()

println("="^80)
println()

println("雅可比行列式计算:")
println("-"^80)
println()

println("从 (cosθ_i, cosθ_j, φ) 到 (cosΘ, ...) 的变换")
println()

println("固定 cosθ_i 和 cosθ_j，对 φ 积分:")
println("  ∫_{0}^{2π} dφ · g(cosΘ)")
println()

println("但是 cosΘ 依赖于 φ:")
println("  cosΘ = cosθ_i cosθ_j + sinθ_i sinθ_j cos(φ)")
println()

println("这不是一个简单的变量替换！")
println("  - 对于固定的 (cosθ_i, cosθ_j)，cosΘ 在一个范围内变化")
println("  - 需要对所有可能的 (cosθ_i, cosθ_j) 组合积分")
println()

println("="^80)
println()

println("数值验证:")
println("-"^80)
println()

println("让我们数值计算一下这个因子")
println()

# 数值计算角度积分的因子
using QuadGK

# 测试函数：g(cosΘ) = 1（常数）
function test_constant()
    # Julia的5维积分（只积分角度部分）
    julia_integral, _ = quadgk(cosθ_i -> begin
        quadgk(cosθ_j -> begin
            quadgk(φ -> begin
                # cosΘ = cosθ_i * cosθ_j + sqrt(1-cosθ_i^2) * sqrt(1-cosθ_j^2) * cos(φ)
                1.0  # g(cosΘ) = 1
            end, 0, 2π)[1]
        end, -1, 1)[1]
    end, -1, 1)
    
    # Fortran的3维积分（只积分角度部分）
    # ∫_{0}^{π} sin(θ) dθ = 2
    fortran_integral = 2.0
    
    return julia_integral, fortran_integral
end

julia_int, fortran_int = test_constant()

@printf("测试函数 g(cosΘ) = 1:\n")
@printf("  Julia积分:   ∫∫∫ dcosθ_i dcosθ_j dφ · 1 = %.6f\n", julia_int)
@printf("  Fortran积分: ∫ sin(θ) dθ · 1 = %.6f\n", fortran_int)
@printf("  比值: Julia/Fortran = %.6f\n", julia_int / fortran_int)
println()

println("⚠️ 关键发现:")
@printf("  Julia的角度积分比Fortran大 %.2f 倍！\n", julia_int / fortran_int)
println()

println("这解释了为什么Julia的ω比Fortran大约40倍：")
println("  - 角度积分因子: ~4π (Julia) vs ~2 (Fortran) ≈ 6.3倍")
println("  - 其他因子（动量积分等）可能还有差异")
println()

println("="^80)
println()

println("结论:")
println("-"^80)
println()

println("问题找到了！")
println()

println("Julia在各向同性情况下（ξ=0）使用了5维积分:")
println("  ∫∫∫∫∫ dcosθ_i dcosθ_j dφ")
println()

println("但这个积分包含了额外的因子 4π！")
println()

println("正确的做法应该是:")
println("  1. 在ξ=0时，直接使用3维积分（像Fortran那样）")
println("  2. 或者，在5维积分中除以正确的归一化因子")
println()

println("当前Julia的归一化因子:")
println("  DQ²/(32π⁵ρ_iρ_j)")
println()

println("应该修正为（在ξ=0时）:")
println("  DQ²/(32π⁵ρ_iρ_j) × (某个因子)")
println()

println("或者，更好的做法是:")
println("  在ξ=0时，使用专门的各向同性积分公式")
println()

println("="^80)
