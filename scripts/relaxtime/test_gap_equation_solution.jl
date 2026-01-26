#!/usr/bin/env julia
"""
测试 Julia 的能隙方程求解,并与 Fortran 结果对比

目标: 确认在 T=300 MeV, μ=2 MeV 下,Julia 和 Fortran 的序参量是否一致
"""

using Printf

# 包含必要的模块
include("../../src/pnjl/PNJL.jl")
using .PNJL

include("../../src/Constants_PNJL.jl")
using .Constants_PNJL

# 测试参数
const ħc = 0.1973269804  # GeV·fm
T_MeV = 300.0
μ_B_MeV = 6.0  # 重子化学势 = 3 × 夸克化学势
μ_MeV = μ_B_MeV / 3.0  # 夸克化学势
T = T_MeV / (ħc * 1000.0)
μ = μ_MeV / (ħc * 1000.0)

println("=" ^ 80)
println("Julia 能隙方程求解测试")
println("=" ^ 80)
println()
println("测试条件:")
@printf("  T = %.6f fm⁻¹ (%.1f MeV)\n", T, T_MeV)
@printf("  μ_B = %.6f fm⁻¹ (%.1f MeV)\n", 3*μ, μ_B_MeV)
@printf("  μ_q = %.6f fm⁻¹ (%.1f MeV)\n", μ, μ_MeV)
println()

println("=" ^ 80)
println("求解能隙方程...")
println("=" ^ 80)
println()

# 求解能隙方程
try
    result = PNJL.solve(
        PNJL.FixedMu(),
        T,
        μ;  # 夸克化学势
        xi=0.0,
        p_num=64,
        t_num=16
    )
    
    println("求解成功!")
    @printf("  收敛: %s\n", result.converged)
    @printf("  迭代次数: %d\n", result.iterations)
    @printf("  残差: %.6e\n", result.residual_norm)
    println()
    
    # 提取解
    x_state = result.x_state
    φ_u = x_state[1]
    φ_d = x_state[2]
    φ_s = x_state[3]
    Φ = x_state[4]
    Φbar = x_state[5]
    
    # 计算质量
    masses = result.masses
    m_u = masses[1]
    m_d = masses[2]
    m_s = masses[3]
    
    println("=" ^ 80)
    println("Julia 的解")
    println("=" ^ 80)
    println()
    
    @printf("序参量:\n")
    @printf("  φ_u = %.6e fm⁻³\n", φ_u)
    @printf("  φ_d = %.6e fm⁻³\n", φ_d)
    @printf("  φ_s = %.6e fm⁻³\n", φ_s)
    println()
    
    @printf("Polyakov 环:\n")
    @printf("  Φ = %.8f\n", Φ)
    @printf("  Φbar = %.8f\n", Φbar)
    println()
    
    @printf("有效质量:\n")
    @printf("  m_u = %.5f fm⁻¹ (%.2f MeV)\n", m_u, m_u * ħc * 1000.0)
    @printf("  m_d = %.5f fm⁻¹ (%.2f MeV)\n", m_d, m_d * ħc * 1000.0)
    @printf("  m_s = %.5f fm⁻¹ (%.2f MeV)\n", m_s, m_s * ħc * 1000.0)
    println()
    
    println("=" ^ 80)
    println("Fortran 的解 (从调试输出)")
    println("=" ^ 80)
    println()
    
    m_u_fortran = 0.04051  # fm⁻¹
    m_s_fortran = 1.0323   # fm⁻¹
    Phi_fortran = 0.99999994
    Phibar_fortran = 0.99999994
    
    @printf("有效质量:\n")
    @printf("  m_u = %.5f fm⁻¹ (%.2f MeV)\n", m_u_fortran, m_u_fortran * ħc * 1000.0)
    @printf("  m_s = %.5f fm⁻¹ (%.2f MeV)\n", m_s_fortran, m_s_fortran * ħc * 1000.0)
    println()
    
    @printf("Polyakov 环:\n")
    @printf("  Φ = %.8f\n", Phi_fortran)
    @printf("  Φbar = %.8f\n", Phibar_fortran)
    println()
    
    println("=" ^ 80)
    println("对比")
    println("=" ^ 80)
    println()
    
    @printf("m_u 差异:\n")
    @printf("  Julia: %.5f fm⁻¹\n", m_u)
    @printf("  Fortran: %.5f fm⁻¹\n", m_u_fortran)
    @printf("  绝对差异: %.5f fm⁻¹\n", abs(m_u - m_u_fortran))
    @printf("  相对差异: %.2f%%\n", 100.0 * abs(m_u - m_u_fortran) / m_u_fortran)
    println()
    
    @printf("m_s 差异:\n")
    @printf("  Julia: %.5f fm⁻¹\n", m_s)
    @printf("  Fortran: %.5f fm⁻¹\n", m_s_fortran)
    @printf("  绝对差异: %.5f fm⁻¹\n", abs(m_s - m_s_fortran))
    @printf("  相对差异: %.2f%%\n", 100.0 * abs(m_s - m_s_fortran) / m_s_fortran)
    println()
    
    @printf("Φ 差异:\n")
    @printf("  Julia: %.8f\n", Φ)
    @printf("  Fortran: %.8f\n", Phi_fortran)
    @printf("  绝对差异: %.8e\n", abs(Φ - Phi_fortran))
    println()
    
    println("=" ^ 80)
    println("结论")
    println("=" ^ 80)
    println()
    
    if abs(m_u - m_u_fortran) / m_u_fortran < 0.01 && abs(m_s - m_s_fortran) / m_s_fortran < 0.01
        println("✅ 能隙方程的解一致 (差异 < 1%)")
        println()
        println("这意味着:")
        println("  1. Julia 和 Fortran 的能隙方程形式一致")
        println("  2. 序参量和 Polyakov 环的值一致")
        println("  3. A 函数的 31% 差异不是由能隙方程引起的")
        println()
        println("下一步: 需要检查 A 函数本身的计算")
    else
        println("❌ 能隙方程的解不一致!")
        println()
        println("这意味着:")
        println("  1. Julia 和 Fortran 的能隙方程可能有差异")
        println("  2. 或者初始猜测/求解器参数不同")
        println("  3. 或者 A 函数的计算不同,导致能隙方程的解不同")
        println()
        println("下一步: 需要检查能隙方程的形式和 A 函数的计算")
    end
    println()
    
catch e
    println("❌ 求解失败!")
    println("错误信息: ", e)
    println()
    println("可能的原因:")
    println("  1. 初始猜测不合适")
    println("  2. 求解器参数需要调整")
    println("  3. 能隙方程在这个参数点没有解")
    println()
end
